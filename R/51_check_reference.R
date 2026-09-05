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
#
# WHAT THIS VALIDATES, PRECISELY (design doc section 6, "Marker alignment"):
#   - marker set: every marker used by the fitted effects is present in check_geno (hard error,
#     naming the missing count), and the check matrix is subset + reordered to that exact set/
#     order -- so a caller can never silently score against a shuffled or partial marker vector.
#   - dosage RANGE: every value falls in [0, ploidy]. This catches a wrong PLOIDY assumption and
#     most wrong-coding files (e.g. {-1,0,1} or {0,1} dosages), because those ranges collide with
#     [0, ploidy] only by coincidence.
# WHAT THIS DOES NOT, AND CANNOT, VALIDATE:
#   - reference-allele agreement. A check_geno file coded against the OPPOSITE reference allele
#     at a subset of markers (dosage = ploidy - true_dosage for just those columns) produces
#     values that are STILL inside [0, ploidy] at every marker -- the range check cannot
#     distinguish "correctly coded" from "flipped at some markers" because both produce valid-
#     looking dosages. This is exactly the "wrong but entirely plausible" failure mode the design
#     doc calls out: a check GEBV in the right ballpark, in the wrong place, with no warning.
#     Detecting it would require an independent ground truth for the reference allele per marker
#     (e.g. cross-referencing against the parent genotype matrix's own allele-frequency direction,
#     which this function is not given) -- there is no rule over check_geno's dosages alone that
#     can tell a flipped subset of markers from a genuinely different (but valid) check genotype.
#     Relaxing check_geno to optional, or adding a frequency-based heuristic, is a separate,
#     larger decision and is deliberately out of scope here.
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
#
# D1 fix: `var_suffix` (default "_pmv_used") is PMV = VPM + posterior marker-effect uncertainty
# (R/03_metrics.R), and that uncertainty component is SHARED across a family's k progeny (same
# beta-hat scores every progeny), not an independent per-progeny draw. Raising PMV to the k-th
# power -- the previous behaviour here -- is invalid and anti-conservative (see
# ng_p_superior_progeny_pev(), R/30_posterior_prediction.R, for the derivation and the
# closed-form-vs-truth numbers). `vpm_suffix` (default "_vpm") recovers the independent-across-
# progeny variance; PEV per row is v - vpm. When `<key>_vpm` is entirely absent from `scores`,
# vpm is aliased to v for every row of this trait, which makes PEV exactly 0 and reproduces the
# historical closed form exactly (not an approximation) -- this keeps callers that never attached
# a `_vpm` column (e.g. a hand-built `scores` fixture) working unchanged. Where `_vpm` IS present
# but a given row's value is non-finite or negative, that row alone is left not-evaluable (its p
# is NA, like any other not-evaluable row) rather than guessing a variance; counted in
# `n_pev_unavailable`.
ng_attach_check_reference <- function(scores, spec, trait_values = NULL, check_values,
                                      k_progeny,
                                      mean_suffix = "_mean", var_suffix = "_pmv_used",
                                      vpm_suffix = "_vpm") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  n <- nrow(scores)
  k_progeny <- as.integer(k_progeny[[1L]])
  ok_mat <- matrix(NA, nrow = n, ncol = nrow(spec), dimnames = list(NULL, spec$trait))
  n_not_evaluable <- list()
  n_wrong <- list()
  n_pev_unavailable <- list()
  # The cross table names its per-trait columns with the SANITISED trait name
  # (ng_run_cp_clean_trait_name: make.names + dots to underscores), while the spec carries the
  # raw name the breeder typed. Look columns up by the key, report by the raw name.
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  for (k in seq_len(nrow(spec))) {
    tr <- spec$trait[[k]]; ck <- spec$check[[k]]; kk <- key[[k]]
    # reject_if == "below" means the breeder wants the mid-parent ABOVE the check (increase
    # trait); sgn flips the comparison so that a positive margin always means "better".
    sgn <- if (identical(spec$reject_if[[k]], "below")) 1 else -1
    mcol <- paste0(kk, mean_suffix); scol <- paste0(kk, var_suffix); vcol <- paste0(kk, vpm_suffix)
    if (!(mcol %in% names(scores))) ng_stop("scores missing mean column for check trait: ", mcol)
    mu <- suppressWarnings(as.numeric(scores[[mcol]]))
    v <- if (scol %in% names(scores)) suppressWarnings(as.numeric(scores[[scol]])) else rep(NA_real_, n)
    vpm_present <- vcol %in% names(scores)
    vpm <- if (vpm_present) suppressWarnings(as.numeric(scores[[vcol]])) else v
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
    # mu/tau evaluability is independent of the variance source; vpm_ok is the D1 requirement
    # that the independent-across-progeny variance be recoverable at all (v itself, when the
    # _vpm column is absent, per the fallback above).
    core_ok <- is.finite(mu) & is.finite(tau)
    vpm_ok <- is.finite(vpm) & vpm >= 0
    usable <- core_ok & vpm_ok
    pev_ok <- vpm_ok & is.finite(v) & (v >= vpm - 1e-8)
    pev <- rep(0, n)
    pev[usable] <- ifelse(pev_ok[usable], pmax(v[usable] - vpm[usable], 0), 0)
    if (any(usable)) {
      # sgn folds the decrease case into the same closed form: negating both mu and tau turns
      # P(at least one of k progeny >= tau) into P(at least one <= tau); delta ~ N(0, PEV) is
      # symmetric about 0, so integrating it over the sgn-flipped axis answers the same mirrored
      # question unmodified (see ng_p_superior_progeny_pev()'s docstring).
      p[usable] <- ng_p_superior_progeny_pev(sgn * mu[usable], vpm[usable], pev[usable],
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
    # Rows where mu/tau were otherwise fine (core_ok) but VPM could not be recovered as a valid
    # non-negative number -- these rows are left NA in _p_beat_check rather than guessing a
    # variance; distinct from n_not_evaluable, which is driven by the check/ok flag (mean/tau),
    # not by the variance source.
    n_pev_unavailable[[tr]] <- sum(core_ok & !vpm_ok)
  }
  # Three-valued logic via base R's own Kleene `all()`: FALSE wins outright (a real failure),
  # NA propagates when nothing resolves the row (all-NA, or a mix of NA and TRUE with no FALSE --
  # in the latter case we genuinely do not know whether the row passes), and TRUE only when every
  # evaluated check passed and none were NA. This is what actually implements "NA never counts as
  # a failure" without also turning an all-NA (never evaluated) row into a false affirmative pass
  # -- the bug the previous `!any(r %in% FALSE)` had (NA %in% FALSE is FALSE, so a row of all-NA
  # silently reported checks_all_ok = TRUE).
  scores$checks_all_ok <- apply(ok_mat, 1L, all)
  # INTEGRATION: a dimensionless per-cross violation COUNT (never a trait-unit quantity, which is
  # exactly why it may combine across traits) -- how many active checks this cross falls on the
  # wrong side of. NA-aware: an unevaluable check (NA in ok_mat) counts as neither pass nor fail,
  # i.e. contributes 0 to the sum; a row where NOTHING was evaluable is NA (unknown), never a
  # false-affirmative 0. Feeds ng_rank_cross_priority()'s optional check_weight component.
  evaluated <- !is.na(ok_mat)
  viol <- matrix(0, nrow(ok_mat), ncol(ok_mat))
  viol[evaluated] <- as.numeric(ok_mat[evaluated] %in% FALSE)
  any_evaluated <- rowSums(evaluated) > 0
  check_violation <- rep(NA_real_, n)
  if (any(any_evaluated)) {
    check_violation[any_evaluated] <- rowSums(viol, na.rm = TRUE)[any_evaluated]
  }
  scores$check_violation <- as.integer(check_violation)
  attr(scores, "check_reference_diagnostics") <- list(
    active = spec, n_wrong_side = n_wrong,
    n_not_evaluable = n_not_evaluable, n_pev_unavailable = n_pev_unavailable, n_candidates = n)
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
#
# D2 fix: `var_suffix` now defaults to "_vpm", not "_pmv_used". PMV mixes in the SHARED posterior
# marker-effect uncertainty that ng_p_superior_progeny_pev() (R/30_posterior_prediction.R, D1)
# integrates out of the per-trait marginals; using it as the joint's diagonal made the joint
# incommensurable with the (corrected) marginals -- empirically able to EXCEED them, which is
# impossible for a true joint-vs-marginal pair. VPM is the SAME diagonal
# ng_cross_trait_within_family_cov() (the exact within-family cross-trait covariance, a'Ra)
# already puts on the diagonal when `cross_trait_cov` is supplied, so both the exact-covariance
# branch and this fallback (population G_hat / independence) branch of
# ng_add_p_superior_progeny_multitrait() now agree with each other AND with the marginals'
# variance. What the joint does NOT do (integrating the multivariate case over posterior effect
# uncertainty is not attempted here) is documented, not silently different: the joint is
# CONDITIONAL on the point-estimated marker effects (delta = 0), while the D1-corrected marginals
# integrate delta ~ N(0, PEV). This is recorded in the returned attribute so a caller inspecting
# the result is told, not left to infer it.
#
# D3 fix: ng_check_tau_bounds() leaves an unevaluable trait's tau_lower/tau_upper at (-Inf, Inf)
# -- unbounded, not excluded. With NOTHING evaluable, pmvnorm over the whole space returns 1, the
# same false affirmative checks_all_ok was fixed to avoid. p_beat_all_checks promises the progeny
# beats EVERY active check; when even one active check cannot be evaluated that promise cannot be
# honestly reported under this name, so the column is NA_real_ rather than silently covering only
# the evaluable subset.
#
# D7 fix: a real-data run (147 parents x 3189 markers, two correlated `decrease` traits, most
# crosses on the wrong side of the check) produced p_beat_all_checks > a single-trait
# <trait>_p_beat_check -- impossible for a true joint-vs-marginal pair, since beating every check
# is a subset of beating any one of them. The cause: the D1 fix made the MARGINAL integrate shared
# posterior effect uncertainty (delta ~ N(0, PEV)) while this joint stayed conditional on
# delta = 0 (VPM diagonal only) -- two different probability models describing the same crosses.
# `pmv_suffix` (default "_pmv_used", matching ng_attach_check_reference's own default for the SAME
# quantity) recovers PEV = PMV - VPM per trait per cross, exactly the way the marginal does; when
# supplied and genuinely positive somewhere, ng_add_p_superior_progeny_multitrait() integrates the
# SAME shared delta (now a per-trait vector) via Monte Carlo
# (ng_p_superior_progeny_multitrait_pev(), R/33). Sigma_c (the VPM/exact-cov diagonal) is
# UNCHANGED -- only the mean vector shifts by delta per draw, mirroring the marginal exactly.
# Wherever PMV is absent or PEV collapses to 0 everywhere, this reproduces the pre-D7 value
# EXACTLY (see tests/check_reference_joint_pev.R), so a caller that never attached a `_pmv_used`
# column (e.g. a hand-built `scores` fixture, or Task 4's fixture in tests/check_reference.R)
# keeps working unchanged.
#
# PROVABLE CEILING (not a new approximation): for ANY fixed delta, the rectangle event "one
# progeny clears every check at once" is a SUBSET of "one progeny clears check t alone", for every
# t, so the family-level probability 1-(1-p_one)^k is monotone in that subset relationship too --
# p_beat_all_checks(delta) <= p_beat_check_t(delta) for every delta, hence (integrating the SAME
# delta_t marginal on both sides) E_delta[joint] <= E_{delta_t}[marginal_t] EXACTLY, in the true
# (infinite-draw) population quantities. The Monte Carlo estimate above can occasionally overshoot
# that provable ceiling by a small finite-sample amount (see the fix report's MC-error table); this
# block computes each trait's marginal via the SAME closed form ng_attach_check_reference() uses
# (ng_p_superior_progeny_pev(), sgn-mirrored per spec$reject_if) and clips the joint to the
# row-wise minimum. This can only move the reported estimate CLOSER to the unknown true value,
# never further from it, and it is what makes the invariant hold EXACTLY (not just "usually,
# within Monte Carlo noise") for every row where the marginals themselves are evaluable.
ng_attach_joint_check_probability <- function(scores, spec, check_values, k_progeny,
                                              mean_suffix = "_mean", var_suffix = "_vpm",
                                              pmv_suffix = "_pmv_used",
                                              cross_trait_cov = NULL) {
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  b <- ng_check_tau_bounds(spec, check_values)
  evaluable_trait <- is.finite(as.numeric(b$tau_lower[spec$trait])) |
    is.finite(as.numeric(b$tau_upper[spec$trait]))
  if (!length(evaluable_trait) || !all(evaluable_trait)) {
    scores <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
    scores$p_beat_all_checks <- NA_real_
    return(scores)
  }
  scores_df <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  trait_specs <- data.frame(
    trait = spec$trait,
    mean_col = paste0(key, mean_suffix),
    var_col = paste0(key, var_suffix),
    stringsAsFactors = FALSE)
  # Per-cross, per-trait PEV = PMV - VPM, mirroring ng_attach_check_reference(): 0 wherever the
  # PMV column is absent, non-finite, or (numerically) below VPM, rather than NA -- an invalid or
  # unavailable PEV for one trait/row falls back to "no correction for that trait/row", it never
  # blocks the row's joint probability the way an unevaluable CHECK (tau) does.
  n <- nrow(scores_df)
  t_n <- nrow(trait_specs)
  pev_mat <- matrix(0, nrow = n, ncol = t_n)
  for (kk in seq_len(t_n)) {
    vcol <- trait_specs$var_col[[kk]]
    pcol <- paste0(key[[kk]], pmv_suffix)
    if (pcol %in% names(scores_df) && vcol %in% names(scores_df)) {
      vpmv <- suppressWarnings(as.numeric(scores_df[[vcol]]))
      pmvv <- suppressWarnings(as.numeric(scores_df[[pcol]]))
      ok <- is.finite(vpmv) & is.finite(pmvv) & (pmvv >= vpmv - 1e-8)
      pev_mat[ok, kk] <- pmax(pmvv[ok] - vpmv[ok], 0)
    }
  }
  pev_mat_mc <- if (any(is.finite(pev_mat) & pev_mat > 0)) pev_mat else NULL
  out <- ng_add_p_superior_progeny_multitrait(
    scores_df, trait_specs,
    tau_lower = as.numeric(b$tau_lower[spec$trait]),
    tau_upper = as.numeric(b$tau_upper[spec$trait]),
    k_progeny = k_progeny, cross_trait_cov = cross_trait_cov,
    pev_mat = pev_mat_mc,
    out_col = "p_beat_all_checks")
  # The provable ceiling: each trait's own PEV-integrated marginal, computed the same way
  # ng_attach_check_reference() computes <trait>_p_beat_check (same closed form, same sgn
  # convention), independent of whether the caller separately attached those columns.
  marg_mat <- matrix(NA_real_, nrow = n, ncol = t_n)
  for (kk in seq_len(t_n)) {
    tr <- spec$trait[[kk]]; ck <- spec$check[[kk]]
    sgn <- if (identical(spec$reject_if[[kk]], "below")) 1 else -1
    mu_k <- suppressWarnings(as.numeric(scores_df[[trait_specs$mean_col[[kk]]]]))
    vpm_k <- suppressWarnings(as.numeric(scores_df[[trait_specs$var_col[[kk]]]]))
    cv <- check_values[[tr]]
    tau_k <- suppressWarnings(as.numeric(if (!is.null(cv) && ck %in% names(cv)) cv[[ck]] else NA_real_))
    ok <- is.finite(mu_k) & is.finite(vpm_k) & vpm_k >= 0 & is.finite(tau_k)
    if (any(ok)) {
      marg_mat[ok, kk] <- ng_p_superior_progeny_pev(sgn * mu_k[ok], vpm_k[ok], pev_mat[ok, kk],
                                                     sgn * tau_k, k_progeny)
    }
  }
  marg_min <- apply(marg_mat, 1L, function(r) if (anyNA(r)) NA_real_ else min(r))
  out$p_beat_all_checks <- pmin(out$p_beat_all_checks, marg_min)
  attr(out, "p_beat_all_checks_note") <- if (!is.null(pev_mat_mc)) paste(
    "p_beat_all_checks integrates shared posterior marker-effect uncertainty (PEV, a diagonal",
    "per-trait approximation -- no cross-trait PEV covariance is estimated) via a fixed",
    sprintf("%d-draw seeded Monte Carlo (D7), the same shared-delta model", NG_JOINT_PEV_MC_DRAWS),
    "<trait>_p_beat_check (D1) integrates via Gauss-Hermite, then is clipped to the row-wise",
    "minimum of the (independently recomputed) per-trait marginals -- a provable ceiling, not a",
    "new approximation -- so the invariant p_beat_all_checks <= min(marginals) holds exactly.") else paste(
    "p_beat_all_checks is conditional on the point-estimated marker effects (VPM diagonal);",
    "no _pmv_used column was available to recover PEV, so this does not integrate posterior",
    "marker-effect uncertainty and is not comparable to a PEV-integrated <trait>_p_beat_check.")
  out
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
