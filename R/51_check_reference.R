# Check lines as REFERENCES, not filters. A check is a benchmark genotype (a released variety,
# a commercial check) that is never crossed: check x check would be a self, which
# ng_make_pairs() does not enumerate, so a check has no variance and can never be a row in the
# cross table. It contributes a scalar per trait, on the same mean_source the cross means use.
# See docs/design/2026-09-03-check-reference-lines-design.md.

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
                                      mean_suffix = "_mean", sd_suffix = "_pmv_used") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  n <- nrow(scores)
  k_progeny <- as.integer(k_progeny[[1L]])
  ok_mat <- matrix(NA, nrow = n, ncol = nrow(spec), dimnames = list(NULL, spec$trait))
  n_not_evaluable <- 0L
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
    mcol <- paste0(kk, mean_suffix); scol <- paste0(kk, sd_suffix)
    if (!(mcol %in% names(scores))) ng_stop("scores missing mean column for check trait: ", mcol)
    mu <- suppressWarnings(as.numeric(scores[[mcol]]))
    v <- if (scol %in% names(scores)) suppressWarnings(as.numeric(scores[[scol]])) else rep(NA_real_, n)
    tau <- suppressWarnings(as.numeric(check_values[[tr]][[ck]]))
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
    n_not_evaluable <- n_not_evaluable + sum(is.na(ok))
    n_wrong[[tr]] <- sum(ok %in% FALSE)
  }
  # NA never counts as a failure: an unevaluable check is reported, not held against a cross.
  scores$checks_all_ok <- !apply(ok_mat, 1L, function(r) any(r %in% FALSE))
  attr(scores, "check_reference_diagnostics") <- list(
    active = spec, n_wrong_side = n_wrong,
    n_not_evaluable = as.integer(n_not_evaluable), n_candidates = n)
  rownames(scores) <- NULL
  scores
}
