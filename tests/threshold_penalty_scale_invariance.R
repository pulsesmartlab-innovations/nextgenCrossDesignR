helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# DEFECT 3 (multi-trait audit, R/19_multi_trait_selection.R): the soft threshold penalty
# measured its deficit on `threshold_column` (trait units, e.g. <trait>_mean) but divided it by
# the robust spread of `column` -- the RANKING metric. Those are different quantities in
# different units, so a knob with nothing to do with thresholds (trait_value_metric) rescaled
# the penalty: the audit measured a mean violation of 0.479 when the ranking column was
# `protein_usefulness` versus 18.0 when it was `protein_pmv` on identical data with identical
# thresholds -- a ~38x swing in the deficit and a 52x swing in the applied penalty.
#
# The invariant asserted here: with the threshold column and the thresholds held fixed, the
# per-trait violation must not change when the RANKING column changes.

set.seed(11)
n <- 40L
protein_mean <- seq(60, 80, length.out = n)           # trait units; thresholds live here
# Ranking columns on wildly different scales. `usefulness` is in trait units (spread ~ O(1));
# `pmv` is a VARIANCE (spread ~ O(0.02)). Under the defect the pmv run divides the same
# trait-unit deficit by a variance IQR and blows the penalty up by ~40x.
protein_usefulness <- protein_mean + 1.755 * sqrt(seq(0.5, 2.0, length.out = n))
protein_pmv <- seq(0.02, 0.09, length.out = n)

base <- data.frame(
  parent1 = paste0("P", seq_len(n)),
  parent2 = paste0("Q", seq_len(n)),
  protein_mean = protein_mean,
  stringsAsFactors = FALSE
)

traits_for <- function(rank_col) {
  data.frame(
    trait = "protein",
    column = rank_col,
    threshold_column = "protein_mean",
    direction = "maximize",
    weight = 1,
    min_value = 72,          # bites on the lower ~60% of the pool
    max_value = NA_real_,
    stringsAsFactors = FALSE
  )
}

s_use <- base; s_use$protein_usefulness <- protein_usefulness
s_pmv <- base; s_pmv$protein_pmv <- protein_pmv

scored_use <- ng_add_multitrait_score(s_use, traits_for("protein_usefulness"), method = "weighted")
scored_pmv <- ng_add_multitrait_score(s_pmv, traits_for("protein_pmv"), method = "weighted")

v_use <- scored_use$multi_trait_threshold_violation
v_pmv <- scored_pmv$multi_trait_threshold_violation

cat(sprintf("mean violation, ranking column = protein_usefulness : %.6f\n", mean(v_use)))
cat(sprintf("mean violation, ranking column = protein_pmv        : %.6f\n", mean(v_pmv)))

# The thresholds bite: this test would be vacuous if nothing violated.
stopifnot(sum(v_use > 0) > 5L)

# THE INVARIANT. Same threshold column, same thresholds, different ranking metric -> identical
# violation. Pre-fix this ratio was ~38.
stopifnot(isTRUE(all.equal(v_use, v_pmv, tolerance = 1e-12)))

# The deficit is scaled by the dispersion of the column it is measured on, so it is exactly the
# trait-unit shortfall divided by the robust spread of `protein_mean`.
expected_scale <- IQR(protein_mean) / 1.349
expected <- pmax(72 - protein_mean, 0) / expected_scale
stopifnot(isTRUE(all.equal(as.numeric(v_use), as.numeric(expected), tolerance = 1e-12)))

# NON-REGRESSION: when threshold_column IS column (the default for every direct caller of
# ng_add_multitrait_score / ng_breeder_selection_objective), the scale is unchanged and the
# result must be bit-identical to 0.25.0's behaviour -- deficit / IQR(column)/1.349.
same_col <- data.frame(
  parent1 = base$parent1, parent2 = base$parent2,
  protein_value = protein_mean, stringsAsFactors = FALSE
)
tr_same <- data.frame(trait = "protein", column = "protein_value", direction = "maximize",
                      weight = 1, min_value = 72, max_value = NA_real_,
                      stringsAsFactors = FALSE)
scored_same <- ng_add_multitrait_score(same_col, tr_same, method = "weighted")
stopifnot(isTRUE(all.equal(as.numeric(scored_same$multi_trait_threshold_violation),
                           as.numeric(expected), tolerance = 1e-12)))

# A degenerate threshold column (zero IQR, zero MAD, zero SD, zero range) must not divide by
# zero: ng_multitrait_value_scale() falls back to 1 and the deficit stays finite.
flat <- data.frame(parent1 = base$parent1, parent2 = base$parent2,
                   protein_mean = rep(50, n), protein_usefulness = protein_usefulness,
                   stringsAsFactors = FALSE)
scored_flat <- ng_add_multitrait_score(flat, traits_for("protein_usefulness"), method = "weighted")
stopifnot(all(is.finite(scored_flat$multi_trait_threshold_violation)))
stopifnot(isTRUE(all.equal(unique(round(scored_flat$multi_trait_threshold_violation, 10)), 22)))

cat("threshold_penalty_scale_invariance: OK\n")
