# Polyploid (any-ploidy) mate design, built by COMPOSITION of the ploidy-general stack.
# The package separates cross scoring from mate allocation (FRAMEWORK.md); the allocator
# (ng_optimize_mating_plan) and all its controls -- strategy dial, target_coancestry, committed
# matings, group permission/quotas, cost/budget/logistics, min-use, the evolution optimizer --
# operate purely on a candidate-cross table + a kinship matrix and are ploidy-agnostic. So the
# only ploidy-specific piece needed is a scorer that emits that standard contract. These two
# functions add exactly that plus a one-call entry point, and reuse everything downstream.
#
# Scope: MEAN-based (mid-parent GEBV). Evidence (VALIDATED_STATE.md): the mean criterion is at
# least as good as usefulness for realized gain at practical horizons, so a mean scorer is a
# correct, useful default. The autopolyploid within-family VARIANCE kernel (double reduction +
# segregation) is deliberately NOT included here -- it is a separate, unvalidated derivation.
# Kinship uses the correct allele-frequency-based polyploid GRM (ng_polyploid_grm, R/47).

# Analytic, no-simulation, any-ploidy ADDITIVE cross scorer. Emits parent1/parent2 with the
# mid-parent breeding value (poly_mean), the within-family additive segregation variance
# (poly_var, from the progeny-moment table), the usefulness (poly_usefulness = mean + i*SD),
# and pair_kinship, plus a parent_kinship (polyploid GRM) attribute. This is the standard candidate-cross
# contract the allocator consumes; select on poly_mean (gain) or poly_usefulness (variance-aware).
ng_polyploid_score_crosses <- function(dosage, effects, ploidy = 2L, pairs = NULL,
                                  grm_method = c("vanraden", "yang"),
                                  selection_prop = 0.10, double_reduction = 0) {
  grm_method <- match.arg(grm_method)
  geno <- ng_polyploid_as_dosage_matrix(dosage, ploidy = ploidy, name = "dosage")
  ids <- rownames(geno)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) ng_stop("dosage must have parent row names")
  if (anyDuplicated(ids)) ng_stop("dosage parent row names must be unique")
  eff <- suppressWarnings(as.numeric(effects))
  if (length(eff) != ncol(geno) || any(!is.finite(eff))) {
    ng_stop("effects must provide one finite value per marker column (", ncol(geno), ")")
  }
  gebv <- as.numeric(geno %*% eff); names(gebv) <- ids
  # Correct allele-frequency-based polyploid GRM (VanRaden/Yang generalized to ploidy), not the
  # ploidy-midpoint shortcut. See ng_polyploid_grm (R/47).
  parent_kinship <- ng_polyploid_grm(geno, ploidy = ploidy, method = grm_method)

  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = FALSE)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) ng_stop("pairs must contain parent1 and parent2")
  p1 <- as.character(pairs$parent1); p2 <- as.character(pairs$parent2)
  unknown <- setdiff(unique(c(p1, p2)), ids)
  if (length(unknown)) ng_stop("pairs contains unknown parent IDs: ", paste(unknown, collapse = ", "))
  if (any(p1 == p2)) ng_stop("pairs must not contain self-crosses")

  # Within-family additive segregation variance from the progeny-moment table:
  # Var = sum_k a_k^2 Var(X_k). The between-locus term is zero (R = I) because autopolyploid
  # parental phase is not identifiable from dosage; averaged over the phase configurations
  # consistent with the dosages that covariance is EXACTLY zero, so this is unbiased rather than
  # approximate -- but it cannot separate two crosses that differ only in linkage phase.
  # See R/49 header for the full statement.
  mt <- ng_polyploid_progeny_moment_table(ploidy, double_reduction = double_reduction)
  Mi <- geno; storage.mode(Mi) <- "integer"
  poly_var <- vapply(seq_along(p1), function(i) {
    ij <- cbind(Mi[p1[i], ] + 1L, Mi[p2[i], ] + 1L)
    sum(eff^2 * mt$varX[ij])
  }, numeric(1))
  intensity <- ng_selection_intensity(selection_prop)

  out <- data.frame(
    parent1 = p1, parent2 = p2,
    poly_mean = (gebv[p1] + gebv[p2]) / 2,
    poly_var = poly_var,
    poly_usefulness = (gebv[p1] + gebv[p2]) / 2 + intensity * sqrt(pmax(poly_var, 0)),
    poly_parent1_gebv = gebv[p1], poly_parent2_gebv = gebv[p2],
    variance_model = "unlinked_phase_marginalized",
    pair_kinship = ng_poly4x_pair_coancestry(parent_kinship, pairs),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  attr(out, "parent_kinship") <- parent_kinship
  attr(out, "parent_gebv") <- gebv
  attr(out, "ploidy") <- as.integer(ploidy)
  out
}

# One-call any-ploidy mate design: score (mean + kinship) -> native allocator with the full
# control suite (passed via ...) and ploidy-aware QC. Supply per-marker `effects`, or per-parent
# Supply per-marker `effects`, or per-parent `phenotype` to estimate additive effects (ridge).
ng_polyploid_design_crosses <- function(dosage,
                                   n_crosses,
                                   ploidy = 2L,
                                   effects = NULL,
                                   phenotype = NULL,
                                   pairs = NULL,
                                   max_crosses_per_parent = 4L,
                                   ridge_seed = 1L,
                                   run_qc = TRUE,
                                   qc = list(),
                                   dominance = FALSE,
                                   gain = c("mean", "usefulness"),
                                   selection_prop = 0.10,
                                   double_reduction = 0,
                                   grm_method = c("vanraden", "yang"),
                                   ...) {
  gain <- match.arg(gain); grm_method <- match.arg(grm_method)
  # Ploidy-aware QC first (dosage range 0..ploidy, missingness, MAF, monomorphic, duplicates);
  # cleans the dosage matrix before scoring. Set run_qc = FALSE to skip.
  if (isTRUE(run_qc)) {
    qc_res <- do.call(ng_polyploid_qc, c(list(dosage = dosage, ploidy = ploidy), qc))
    dosage <- qc_res$clean
    if (!is.null(effects) && length(effects) != ncol(dosage)) {
      # keep effects aligned to the markers QC retained
      kept <- !qc_res$marker_report$dropped
      effects <- as.numeric(effects)[kept]
    }
  }
  geno <- ng_polyploid_as_dosage_matrix(dosage, ploidy = ploidy, name = "dosage")

  if (isTRUE(dominance)) {
    # OPTIONAL additive + dominance path (clonal / heterosis crops): estimate both marker-effect
    # components and score crosses on genotypic value = mid-parent breeding value + heterosis, with
    # a within-family additive+dominance variance. Requires a phenotype to estimate dominance.
    if (is.null(phenotype)) ng_stop("dominance = TRUE needs a phenotype to estimate dominance effects")
    y <- suppressWarnings(as.numeric(phenotype))
    names(y) <- if (!is.null(names(phenotype))) names(phenotype) else rownames(geno)
    fit <- ng_polyploid_fit_effects(geno, y[rownames(geno)], ploidy = ploidy,
                                    model = "additive_dominance", seed = ridge_seed)
    scores <- ng_polyploid_score_crosses_dominance(fit, geno, pairs = pairs, selection_prop = selection_prop,
                                              double_reduction = double_reduction, grm_method = grm_method)
    gain_col <- if (identical(gain, "usefulness")) "cross_usefulness" else "cross_mean"
  } else {
    if (is.null(effects)) {
      if (is.null(phenotype)) ng_stop("supply either effects (per marker) or phenotype (per parent)")
      y <- suppressWarnings(as.numeric(phenotype))
      names(y) <- if (!is.null(names(phenotype))) names(phenotype) else rownames(geno)
      y <- y[rownames(geno)]                     # align phenotype to (possibly QC-filtered) samples
      effects <- ng_fit_ridge_effects(geno, y, seed = ridge_seed)$beta
    }
    scores <- ng_polyploid_score_crosses(geno, effects, ploidy = ploidy, pairs = pairs, grm_method = grm_method,
                                    selection_prop = selection_prop, double_reduction = double_reduction)
    gain_col <- if (identical(gain, "usefulness")) "poly_usefulness" else "poly_mean"
  }
  parent_kinship <- attr(scores, "parent_kinship")

  # Native allocation with all forwarded controls (strategy, target_coancestry, committed_crosses,
  # group_permission/quota, cost/budget/logistics, method = "evolution", ...).
  plan <- ng_optimize_mating_plan(scores, n_crosses = n_crosses, gain_col = gain_col,
                                  parent_kinship = parent_kinship, max_crosses_per_parent = max_crosses_per_parent, ...)
  s <- attr(plan, "summary")
  s$ploidy <- as.integer(ploidy); s$poly_gain_col <- gain_col; s$dominance <- isTRUE(dominance)
  attr(plan, "summary") <- s
  attr(plan, "ploidy") <- as.integer(ploidy)
  if (isTRUE(run_qc)) attr(plan, "qc") <- qc_res$summary
  plan
}
