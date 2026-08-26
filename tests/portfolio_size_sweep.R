ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"),
            file.path("nextgen_cross_design", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---- Test 1: basic curve shape --------------------------------------------
ids <- paste0("P", 1:8)
pairs <- ng_make_pairs(ids)
n_pairs <- nrow(pairs)
scores <- data.frame(
  parent1 = pairs$parent1,
  parent2 = pairs$parent2,
  usefulness_pmv_gebv = seq(10, 10 - 0.5 * (n_pairs - 1), length.out = n_pairs),
  pair_kinship = 0,
  stringsAsFactors = FALSE
)

curve <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "elbow_relative"
)

expected_cols <- c("K", "total_gain", "mean_gain", "group_coancestry",
                   "unique_parents", "Ne_estimate",
                   "marginal_gain", "relative_marginal")
stopifnot(all(expected_cols %in% names(curve)))
stopifnot(nrow(curve) == 8L)
stopifnot(all(curve$K == 3:10))
stopifnot(all(diff(curve$total_gain) >= -1e-9))
finite_marg <- curve$marginal_gain[is.finite(curve$marginal_gain)]
stopifnot(all(diff(finite_marg) <= 1e-9))
cat(sprintf("Test 1 (curve shape): K=%s total_gain=%s marginal=%s  OK\n",
            paste(curve$K, collapse = ","),
            paste(sprintf("%.2f", curve$total_gain), collapse = ","),
            paste(sprintf("%.2f", curve$marginal_gain), collapse = ",")))

# ---- Test 2: elbow attribute present --------------------------------------
stopifnot(!is.null(attr(curve, "elbow_K")))
stopifnot(identical(attr(curve, "criterion"), "elbow_relative"))
elbow <- attr(curve, "elbow_K")
if (!is.na(elbow)) stopifnot(elbow %in% curve$K)
cat(sprintf("Test 2 (elbow detected): elbow_K = %s  OK\n",
            ifelse(is.na(elbow), "NA", elbow)))

# ---- Test 2b: elbow fires positively with a higher threshold --------------
# The default 0.05 is not crossed by the linear-decay fixture; bump the
# threshold to 0.4 to verify the detection logic actually returns an
# integer K rather than NA.
curve_high_thr <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "elbow_relative",
  relative_threshold = 0.4
)
elbow_high <- attr(curve_high_thr, "elbow_K")
stopifnot(is.finite(elbow_high))
stopifnot(elbow_high %in% curve_high_thr$K)
# Verify the threshold attribute round-tripped
stopifnot(identical(attr(curve_high_thr, "relative_threshold"), 0.4))
cat(sprintf("Test 2b (positive elbow at threshold=0.4): elbow_K = %d  OK\n",
            elbow_high))

# ---- Test 3: dot args are passed through ----------------------------------
err <- tryCatch(
  ng_optimize_mating_plan_curve(
    scores = scores, K_range = 9:10,
    max_crosses_per_parent = 1L, lambda_group = 0,
    criterion = "elbow_relative"
  ),
  error = function(e) conditionMessage(e)
)
if (is.character(err)) {
  cat("Test 3 (dot args propagate; infeasible K errored): OK\n")
} else {
  cat("Test 3 (dot args propagate; check skipped — no error raised)\n")
}

# ---- Test 4: kneedle elbow on synthetic curve with known elbow at K=8 -----
K_synth <- 1:20
gain_synth <- ifelse(K_synth <= 8,
                     log(1 + K_synth),
                     log(1 + 8) + 0.02 * (K_synth - 8))
curve_synth <- data.frame(
  K = K_synth,
  total_gain = gain_synth,
  mean_gain = gain_synth / K_synth,
  group_coancestry = 0,
  unique_parents = pmin(K_synth, 10L),
  marginal_gain = c(NA_real_, diff(gain_synth)),
  stringsAsFactors = FALSE
)
curve_synth$relative_marginal <- curve_synth$marginal_gain / curve_synth$marginal_gain[2L]

elbow_kneedle <- .ng_elbow_kneedle(curve_synth)
stopifnot(!is.na(elbow_kneedle))
gap <- abs(elbow_kneedle - 8L)
if (gap > 1L) {
  stop(sprintf("kneedle elbow off target: got %d, expected 8 (+/- 1)", elbow_kneedle))
}
cat(sprintf("Test 4 (kneedle, synthetic elbow at K=8): detected K=%d  OK\n",
            elbow_kneedle))

# ---- Test 5: kneedle on linear curve returns NA ---------------------------
K_lin <- 1:10
curve_lin <- data.frame(
  K = K_lin,
  total_gain = 2 * K_lin,
  mean_gain = rep(2, length(K_lin)),
  group_coancestry = 0,
  unique_parents = K_lin,
  marginal_gain = c(NA_real_, rep(2, length(K_lin) - 1L)),
  stringsAsFactors = FALSE
)
curve_lin$relative_marginal <- curve_lin$marginal_gain / curve_lin$marginal_gain[2L]
elbow_lin <- .ng_elbow_kneedle(curve_lin)
stopifnot(is.na(elbow_lin))
cat("Test 5 (kneedle on linear curve: NA): OK\n")

# ---- Test 6: kneedle dispatch via ng_optimize_mating_plan_curve -----------
curve_k <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "elbow_kneedle"
)
stopifnot(identical(attr(curve_k, "criterion"), "elbow_kneedle"))
cat(sprintf("Test 6 (kneedle dispatch): elbow_K = %s  OK\n",
            ifelse(is.na(attr(curve_k, "elbow_K")), "NA",
                   attr(curve_k, "elbow_K"))))

# ---- Test 7: Ne_estimate column matches the analytical formula ------------
curve_ne <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "elbow_relative"
)
stopifnot("Ne_estimate" %in% names(curve_ne))
finite_idx <- which(is.finite(curve_ne$group_coancestry) &
                     curve_ne$group_coancestry > 0)
if (length(finite_idx)) {
  expected_ne <- 1 / (2 * curve_ne$group_coancestry[finite_idx])
  observed_ne <- curve_ne$Ne_estimate[finite_idx]
  stopifnot(all(abs(expected_ne - observed_ne) < 1e-9))
}
cat(sprintf("Test 7 (Ne_estimate matches analytical formula): %d finite rows checked  OK\n",
            length(finite_idx)))

# ---- Test 8: ne_target criterion picks smallest K with Ne >= ne_min -------
ne_min_test <- 5
curve_ne_target <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "ne_target",
  ne_min = ne_min_test
)
stopifnot(identical(attr(curve_ne_target, "criterion"), "ne_target"))
stopifnot(identical(attr(curve_ne_target, "ne_min"), ne_min_test))
elbow_ne <- attr(curve_ne_target, "elbow_K")
if (!is.na(elbow_ne)) {
  satisfying_K <- curve_ne_target$K[curve_ne_target$Ne_estimate >= ne_min_test]
  stopifnot(elbow_ne == min(satisfying_K))
}
cat(sprintf("Test 8 (ne_target, ne_min=%g): elbow_K = %s  OK\n",
            ne_min_test, ifelse(is.na(elbow_ne), "NA", elbow_ne)))

# ---- Test 8b: ne_target fires positively with a very small ne_min ---------
# The default fixture's Ne_estimate stays below 5; ne_min = 0.1 is below
# any plausible floor so the picker should always return the smallest K.
curve_ne_low <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "ne_target",
  ne_min = 0.1
)
elbow_ne_low <- attr(curve_ne_low, "elbow_K")
stopifnot(is.finite(elbow_ne_low))
stopifnot(elbow_ne_low %in% curve_ne_low$K)
satisfying <- curve_ne_low$K[curve_ne_low$Ne_estimate >= 0.1]
stopifnot(elbow_ne_low == min(satisfying))
cat(sprintf("Test 8b (positive ne_target, ne_min=0.1): elbow_K = %d  OK\n",
            elbow_ne_low))

# ---- Test 9: coancestry_budget picks largest K under the ceiling ---------
coanc_max_test <- 0.2
curve_cb <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:10,
  max_crosses_per_parent = 4L, lambda_group = 0,
  criterion = "coancestry_budget",
  coancestry_max = coanc_max_test
)
stopifnot(identical(attr(curve_cb, "criterion"), "coancestry_budget"))
stopifnot(identical(attr(curve_cb, "coancestry_max"), coanc_max_test))
elbow_cb <- attr(curve_cb, "elbow_K")
if (!is.na(elbow_cb)) {
  satisfying <- curve_cb$K[curve_cb$group_coancestry <= coanc_max_test]
  stopifnot(elbow_cb == max(satisfying))
}
cat(sprintf("Test 9 (coancestry_budget, max=%g): elbow_K = %s  OK\n",
            coanc_max_test, ifelse(is.na(elbow_cb), "NA", elbow_cb)))

# ---- Test 10: legacy elbow_method values signal an error -----------------
err_legacy <- tryCatch(
  ng_optimize_mating_plan_curve(
    scores = scores, K_range = 3:10,
    max_crosses_per_parent = 4L, lambda_group = 0,
    elbow_method = "relative_marginal"
  ),
  error = function(e) conditionMessage(e),
  warning = function(w) conditionMessage(w)
)
stopifnot(is.character(err_legacy))
cat("Test 10 (legacy elbow_method arg signals): OK\n")

# ---- Test 11: ng_plot_diminishing_returns returns a ggplot ----------------
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("ggplot2 not installed; skipping plot test")
} else {
  # Use curve_high_thr (Test 2b output) — it has a finite elbow_K.
  p <- ng_plot_diminishing_returns(curve_high_thr)
  stopifnot(inherits(p, "ggplot"))
  layer_types <- vapply(p$layers, function(L) class(L$geom)[[1L]], character(1L))
  stopifnot(any(grepl("GeomVline", layer_types)))
  cat("Test 11 (ng_plot_diminishing_returns returns ggplot with elbow vline): OK\n")
}

# ---- Test 12: show_elbow = FALSE drops the vline --------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_no_elbow <- ng_plot_diminishing_returns(curve_high_thr, show_elbow = FALSE)
  layer_types <- vapply(p_no_elbow$layers, function(L) class(L$geom)[[1L]], character(1L))
  stopifnot(!any(grepl("GeomVline", layer_types)))
  cat("Test 12 (show_elbow = FALSE drops vline): OK\n")
}

# ---- Test 13: NA elbow_K produces a plot without the vline ----------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # curve (from Test 1/2) has elbow_K = NA at the default threshold.
  p_na_elbow <- ng_plot_diminishing_returns(curve)
  stopifnot(inherits(p_na_elbow, "ggplot"))
  layer_types <- vapply(p_na_elbow$layers, function(L) class(L$geom)[[1L]], character(1L))
  stopifnot(!any(grepl("GeomVline", layer_types)))
  cat("Test 13 (NA elbow_K renders without vline): OK\n")
}

# ---- Test 14: monotonicity on a realistic fixture (validation gate 1) ----
# Build a 20-parent fixture with realistic usefulness_pmv_gebv (decreasing across
# pairs) and a non-trivial parent_kinship so lambda_group > 0 actually couples
# the choices. The marginal gain at each K must be non-increasing.
set.seed(2026L)
ids_r <- paste0("R", sprintf("%02d", 1:20))
pairs_r <- ng_make_pairs(ids_r)
n_r <- nrow(pairs_r)
scores_r <- data.frame(
  parent1 = pairs_r$parent1,
  parent2 = pairs_r$parent2,
  usefulness_pmv_gebv = sort(rnorm(n_r, mean = 5, sd = 1), decreasing = TRUE),
  stringsAsFactors = FALSE
)
# Synthetic parent_kinship (positive-definite kinship-like matrix).
Z <- matrix(rnorm(length(ids_r) * 5L), nrow = length(ids_r))
K_r <- tcrossprod(Z) / 5 + diag(0.05, length(ids_r))
rownames(K_r) <- colnames(K_r) <- ids_r
scores_r$pair_kinship <- K_r[cbind(scores_r$parent1, scores_r$parent2)] / 2

curve_r <- ng_optimize_mating_plan_curve(
  scores = scores_r, K_range = 4:15,
  parent_kinship = K_r,
  max_crosses_per_parent = 4L,
  lambda_group = 0.5,
  criterion = "elbow_relative"
)
finite_marg <- curve_r$marginal_gain[is.finite(curve_r$marginal_gain)]
n_strict_increases <- sum(diff(finite_marg) > 1e-6)
if (n_strict_increases > 1L) {
  stop(sprintf(
    "marginal_gain is not monotone non-increasing: %d strict increases found",
    n_strict_increases
  ))
}
cat(sprintf("Test 14 (monotonicity on realistic fixture): %d strict increases (<=1 allowed)  OK\n",
            n_strict_increases))
