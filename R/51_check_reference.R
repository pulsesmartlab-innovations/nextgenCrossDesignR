# Check lines as REFERENCES, not filters. A check is a benchmark genotype (a released variety,
# a commercial check) that is never crossed: check x check would be a self, which
# ng_make_pairs() does not enumerate, so a check has no variance and can never be a row in the
# cross table. It contributes a scalar per trait, on the same mean_source the cross means use.
# See docs/design/2026-09-03-check-reference-lines-design.md.

# Build a validated per-trait check spec (one check per trait). `direction` NA is filled from
# `trait_direction`: increase -> the breeder wants the mid-parent ABOVE the check, so a cross is
# on the wrong side when it is below. There is no `basis` argument: the check always follows the
# run's per-trait mean_source, which is what keeps the reference on the plotted scale.
ng_trait_check_spec <- function(trait, check, direction = NA, trait_direction = NULL) {
  trait <- as.character(trait); check <- as.character(check)
  n <- length(trait)
  if (!n) ng_stop("ng_trait_check_spec needs at least one trait")
  if (length(check) != n) ng_stop("check must be one per trait")
  direction <- tolower(as.character(rep_len(direction, n)))
  need <- is.na(direction) | !nzchar(direction) | direction == "na"
  if (any(need)) {
    if (is.null(trait_direction))
      ng_stop("direction is NA and no trait_direction supplied to resolve it for: ",
              paste(trait[need], collapse = ", "))
    td <- tolower(as.character(trait_direction[trait[need]]))
    if (anyNA(td)) ng_stop("trait_direction has no entry for: ",
                           paste(trait[need][is.na(td)], collapse = ", "))
    direction[need] <- ifelse(td == "increase", "below",
                       ifelse(td == "decrease", "above", NA_character_))
  }
  if (any(!direction %in% c("above", "below")))
    ng_stop("direction must resolve to 'above' or 'below'")
  dup <- unique(trait[duplicated(trait)])
  if (length(dup))
    ng_stop("ng_trait_check_spec: one check per trait; duplicate trait(s): ",
            paste(dup, collapse = ", "))
  data.frame(trait = trait, check = check, reject_if = direction, stringsAsFactors = FALSE)
}

# Align a separate check genotype matrix to the marker set already fixed by parent QC. A check
# matrix that silently disagrees on marker set or column order produces a wrong-but-plausible
# reference value, so every disagreement is an error rather than an intersection.
ng_align_check_geno <- function(check_geno, marker_names, ploidy = 2L) {
  marker_names <- as.character(marker_names)
  check_geno <- as.matrix(check_geno)
  if (is.null(colnames(check_geno))) {
    ng_stop("check_geno must have marker column names to verify alignment with the parent ",
            "marker set; an unnamed matrix cannot be checked for allele-order agreement")
  }
  check_geno <- ng_as_numeric_matrix(check_geno, "check_geno")
  if (is.null(rownames(check_geno))) ng_stop("check_geno must have check ids as rownames")
  miss <- setdiff(marker_names, colnames(check_geno))
  if (length(miss)) {
    ng_stop(sprintf(
      "check_geno is missing %d of %d markers used by the fitted effects (e.g. %s). The check ",
      length(miss), length(marker_names),
      paste(utils::head(miss, 5L), collapse = ", ")),
      "file must carry the same marker set, coding, and reference allele as the parent genotypes.")
  }
  out <- check_geno[, marker_names, drop = FALSE]
  ploidy <- as.integer(ploidy[[1L]])
  rng <- suppressWarnings(range(out, na.rm = TRUE))
  if (any(is.finite(rng)) && (min(rng) < 0 || max(rng) > ploidy)) {
    ng_stop(sprintf(
      "check_geno dosage values fall outside [0, %d] (observed %g..%g); the check file must use ",
      ploidy, rng[[1L]], rng[[2L]]),
      "the same allele coding as the parent genotypes")
  }
  out
}

# Resolve each check's value onto the SAME source the cross means used. ng_choose_mean_source()
# returns one of "GEBV", "GEBV_low_reliability", "GEBV_uncalibrated", "BLUP", "BLUE",
# "adjusted_pheno" -- per trait. A GEBV* source predicts from the check's markers with the
# trait's fitted effects; a phenotypic source reads the check's own record. A check with no
# record on a phenotypic source is NA (not evaluable) and must NEVER silently fall back to a
# GEBV, which would put the reference on a different scale from the axis it is drawn on.
ng_check_reference_value <- function(source, check_geno_aligned, effects, check_records = NULL) {
  source <- as.character(source)[[1L]]
  ids <- rownames(check_geno_aligned)
  if (startsWith(source, "GEBV")) {
    return(stats::setNames(ng_predict_gebv(check_geno_aligned, effects), ids))
  }
  rec <- if (is.null(check_records)) NULL else check_records[[source]]
  if (is.null(rec)) return(stats::setNames(rep(NA_real_, length(ids)), ids))
  v <- suppressWarnings(as.numeric(rec[ids]))
  stats::setNames(v, ids)
}

# Attach per-trait check reference columns to a scored cross table. This function NEVER changes
# nrow(scores) and NEVER reorders it: a check informs the breeder, it does not decide for them.
# (The 0.14.0 predecessor, ng_apply_trait_checks(), dropped rows -- that is the behaviour this
# replaces.)
ng_attach_check_reference <- function(scores, spec, trait_values = NULL, check_values,
                                      k_progeny,
                                      mean_suffix = "_mean", var_suffix = "_pmv_used") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  n <- nrow(scores)
  k_progeny <- as.integer(k_progeny[[1L]])
  ok_mat <- matrix(NA, nrow = n, ncol = nrow(spec), dimnames = list(NULL, spec$trait))
  n_not_evaluable <- list()
  n_wrong <- list()
  # The cross table names its per-trait columns with the SANITISED trait name
  # (ng_run_cp_clean_trait_name: make.names + dots to underscores), while the spec carries the
  # raw name the breeder typed. Look columns up by the key, report by the raw name.
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  for (k in seq_len(nrow(spec))) {
    tr <- spec$trait[[k]]; ck <- spec$check[[k]]; kk <- key[[k]]
    # reject_if == "below" means the breeder wants the mid-parent ABOVE the check (increase
    # trait); sgn flips the comparison so that a positive margin always means "better".
    sgn <- if (identical(spec$reject_if[[k]], "below")) 1 else -1
    mcol <- paste0(kk, mean_suffix); scol <- paste0(kk, var_suffix)
    if (!(mcol %in% names(scores))) ng_stop("scores missing mean column for check trait: ", mcol)
    mu <- suppressWarnings(as.numeric(scores[[mcol]]))
    v <- if (scol %in% names(scores)) suppressWarnings(as.numeric(scores[[scol]])) else rep(NA_real_, n)
    cv <- check_values[[tr]]
    tau <- suppressWarnings(as.numeric(if (!is.null(cv) && ck %in% names(cv)) cv[[ck]] else NA_real_))
    if (!length(tau)) tau <- NA_real_
    margin <- sgn * (mu - tau)
    ok <- margin >= 0                       # ties are not violations
    scores[[paste0(kk, "_check_id")]] <- ck
    scores[[paste0(kk, "_check_value")]] <- tau
    scores[[paste0(kk, "_vs_check")]] <- margin
    scores[[paste0(kk, "_check_ok")]] <- ok
    p <- rep(NA_real_, n)
    usable <- is.finite(mu) & is.finite(v) & v >= 0 & is.finite(tau)
    if (any(usable)) {
      # sgn folds the decrease case into the same closed form: negating both mu and tau turns
      # P(at least one of k progeny >= tau) into P(at least one <= tau).
      p[usable] <- ng_p_superior_progeny(sgn * mu[usable], sqrt(v[usable]),
                                         sgn * tau, k_progeny)
    }
    scores[[paste0(kk, "_p_beat_check")]] <- p
    ok_mat[, k] <- ok
    # Per TRAIT (mirroring n_wrong below), not summed across traits: sum(is.na(ok)) counts the
    # CROSSES for which THIS trait's check was not evaluable (either the check's tau itself is
    # NA, in which case every cross is NA, or an individual cross's own mean/variance is NA).
    # Aggregating across traits into one scalar (the previous behaviour) mixed unrelated traits'
    # not-evaluable crosses together, which cannot honestly be rendered as "N of M crosses" for
    # any single trait.
    n_not_evaluable[[tr]] <- sum(is.na(ok))
    n_wrong[[tr]] <- sum(ok %in% FALSE)
  }
  # Three-valued logic via base R's own Kleene `all()`: FALSE wins outright (a real failure),
  # NA propagates when nothing resolves the row (all-NA, or a mix of NA and TRUE with no FALSE --
  # in the latter case we genuinely do not know whether the row passes), and TRUE only when every
  # evaluated check passed and none were NA. This is what actually implements "NA never counts as
  # a failure" without also turning an all-NA (never evaluated) row into a false affirmative pass
  # -- the bug the previous `!any(r %in% FALSE)` had (NA %in% FALSE is FALSE, so a row of all-NA
  # silently reported checks_all_ok = TRUE).
  scores$checks_all_ok <- apply(ok_mat, 1L, all)
  attr(scores, "check_reference_diagnostics") <- list(
    active = spec, n_wrong_side = n_wrong,
    n_not_evaluable = n_not_evaluable, n_candidates = n)
  rownames(scores) <- NULL
  scores
}

# Translate a check spec into the (tau_lower, tau_upper) pair that
# ng_p_superior_progeny_multitrait() consumes: an increase trait bounds the progeny from below
# at its check, a decrease trait bounds it from above. A check with no evaluable value leaves
# that trait unbounded rather than dropping it, so the joint probability stays defined.
ng_check_tau_bounds <- function(spec, check_values) {
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  tr <- spec$trait
  lo <- stats::setNames(rep(-Inf, length(tr)), tr)
  hi <- stats::setNames(rep(Inf, length(tr)), tr)
  for (k in seq_along(tr)) {
    cv <- check_values[[tr[[k]]]]
    ck <- spec$check[[k]]
    tau <- suppressWarnings(as.numeric(if (!is.null(cv) && ck %in% names(cv)) cv[[ck]] else NA_real_))
    if (!length(tau) || !is.finite(tau)) next
    if (identical(spec$reject_if[[k]], "below")) lo[[tr[[k]]]] <- tau else hi[[tr[[k]]]] <- tau
  }
  list(tau_lower = lo, tau_upper = hi)
}

# P(a progeny beats EVERY check at once). Thin wrapper over the existing multi-trait threshold
# machinery: the check values ARE the tau bounds, so no new probability model is introduced.
ng_attach_joint_check_probability <- function(scores, spec, check_values, k_progeny,
                                              mean_suffix = "_mean", var_suffix = "_pmv_used",
                                              cross_trait_cov = NULL) {
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  b <- ng_check_tau_bounds(spec, check_values)
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  trait_specs <- data.frame(
    trait = spec$trait,
    mean_col = paste0(key, mean_suffix),
    var_col = paste0(key, var_suffix),
    stringsAsFactors = FALSE)
  ng_add_p_superior_progeny_multitrait(
    scores, trait_specs,
    tau_lower = as.numeric(b$tau_lower[spec$trait]),
    tau_upper = as.numeric(b$tau_upper[spec$trait]),
    k_progeny = k_progeny, cross_trait_cov = cross_trait_cov,
    out_col = "p_beat_all_checks")
}

# Convert a check phenotype table -- the same shape as the phenotype file a breeder already
# supplies -- into the check_records structure ng_check_reference_value() consumes. The SOURCE
# KEY is the run's own resolved mean source, not anything named in the input: the check is
# apples-to-apples with the parents' mean by construction. Returns NULL when the run resolved to
# a GEBV source, because the value is then predicted from the check's markers instead.
ng_check_records_from_pheno <- function(check_pheno, id_col, trait_columns, source) {
  source <- as.character(source)[[1L]]
  if (startsWith(source, "GEBV")) return(NULL)
  check_pheno <- as.data.frame(check_pheno, stringsAsFactors = FALSE, check.names = FALSE)
  if (!(id_col %in% names(check_pheno))) {
    ng_stop("check_pheno is missing its id column: ", id_col)
  }
  ids <- as.character(check_pheno[[id_col]])
  out <- list()
  for (tr in names(trait_columns)) {
    col <- trait_columns[[tr]]
    v <- if (col %in% names(check_pheno)) {
      suppressWarnings(as.numeric(check_pheno[[col]]))
    } else {
      rep(NA_real_, length(ids))
    }
    out[[tr]] <- stats::setNames(list(stats::setNames(v, ids)), source)
  }
  out
}
