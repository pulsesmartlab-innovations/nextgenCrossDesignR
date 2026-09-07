helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# Robust posterior allocation at an arbitrary quantile, and direction correctness.
#
# 1. An exact non-default robustness quantile is served from the SAME posterior draws as the
#    credible interval, with no allow_normal_approximation opt-in.
# 2. Asking for a robustness quantile does NOT move the reported credible interval: ci_level
#    keeps owning _post_lower / _post_upper (the collateral-damage guard).
# 3. Direction correctness: for a minimize-oriented ranked value the conservative tail is the
#    UPPER one (1 - q), and the optimizer must prefer LOWER robust values.
# 4. The existing refusal still fires for a quantile that genuinely cannot be served exactly.

set.seed(20260906L)
n <- 24L
m <- 40L
geno <- matrix(2L * rbinom(n * m, 1L, 0.45), nrow = n, ncol = m)
ids <- paste0("L", seq_len(n))
markers <- paste0("M", seq_len(m))
rownames(geno) <- ids; colnames(geno) <- markers
beta_true <- rnorm(m, sd = 0.10)
y <- as.numeric(geno %*% beta_true + rnorm(n, sd = 0.5))
names(y) <- ids
mk <- data.frame(
  marker = markers,
  chr = rep(1:2, length.out = m),
  pos_cm = rep(seq(0, 90, length.out = m / 2), 2)
)
adj <- setNames(as.numeric(y), ids)

post_fit <- ng_fit_ridge_effects_posterior(
  geno = geno, y = y, ids = ids, lambda = 6, kfold = 5L,
  n_draws = 300L, method = "closed_form", seed = 11L
)

predict_args <- list(
  geno = geno, posterior_effects = post_fit, marker_map = mk, ids = ids,
  adjusted_pheno = adj, selection_prop = 0.10, target = "DH",
  recomb_model = "haldane", use_cpp = FALSE, top_n_targets = c(6L, 10L)
)
gain_col <- "usefulness_pmv_gebv"

# ---- (1) baseline vs robustness-quantile run: the CI must be identical ------
base_scores <- do.call(ng_posterior_cross_predict, predict_args)
rq <- 0.25
rq_scores <- do.call(ng_posterior_cross_predict, c(predict_args, list(robustness_quantile = rq)))

q_lo_col <- ng_posterior_quantile_col(gain_col, rq)
q_hi_col <- ng_posterior_quantile_col(gain_col, 1 - rq)
stopifnot(identical(q_lo_col, "usefulness_pmv_gebv_post_q0_25"))
stopifnot(identical(q_hi_col, "usefulness_pmv_gebv_post_q0_75"))
stopifnot(all(c(q_lo_col, q_hi_col) %in% names(rq_scores)))
# Both tails are cached from one request, so the same run serves either direction.
stopifnot(!any(c(q_lo_col, q_hi_col) %in% names(base_scores)))

lo <- paste0(gain_col, "_post_lower")
hi <- paste0(gain_col, "_post_upper")
mn <- paste0(gain_col, "_post_mean")
# COLLATERAL-DAMAGE GUARD: requesting a robustness quantile leaves the reported 95% credible
# interval bit-for-bit unchanged. This is the whole reason ci_level was not repurposed.
for (col in c(lo, hi, mn, "pmv_post_lower", "pmv_post_upper", "ranked_value_post_sd")) {
  if (!isTRUE(all.equal(base_scores[[col]], rq_scores[[col]], tolerance = 0))) {
    stop("requesting robustness_quantile changed the reported column: ", col)
  }
}
stopifnot(identical(attr(rq_scores, "posterior")$ci_level, 0.95))
stopifnot(isTRUE(all.equal(attr(rq_scores, "posterior")$robustness_quantile, rq)))
pq <- attr(rq_scores, "posterior")$posterior_quantiles
stopifnot(is.data.frame(pq), nrow(pq) == 2L)
stopifnot(isTRUE(all.equal(sort(pq$prob), c(0.25, 0.75))))

# The cached quantile is a genuine interior quantile: strictly inside the 95% interval and
# ordered against the median-ish mean, and NOT a copy of a CI tail.
stopifnot(all(rq_scores[[lo]] <= rq_scores[[q_lo_col]] + 1e-12))
stopifnot(all(rq_scores[[q_lo_col]] <= rq_scores[[q_hi_col]] + 1e-12))
stopifnot(all(rq_scores[[q_hi_col]] <= rq_scores[[hi]] + 1e-12))
stopifnot(max(abs(rq_scores[[q_lo_col]] - rq_scores[[lo]])) > 1e-8)

# ---- (2) the exact quantile is now reachable without a normal approximation --
parent_kinship <- ng_parent_kinship(geno)
plan_args <- list(
  posterior_scores = rq_scores, n_crosses = 8L, parent_kinship = parent_kinship,
  gain_col = gain_col, max_crosses_per_parent = Inf, lambda_group = 0,
  lambda_parent_use = 0, method = "greedy_local"
)
plan_max <- do.call(ng_optimize_robust_mating_plan,
                    c(plan_args, list(robustness_quantile = rq)))
s_max <- attr(plan_max, "summary")
stopifnot(nrow(plan_max) == 8L)
stopifnot(isFALSE(s_max$robustness_quantile_is_normal_approximation))
stopifnot(identical(s_max$robust_quantile_source, q_lo_col))
stopifnot(identical(s_max$robust_direction, "maximize"))
stopifnot(isTRUE(all.equal(s_max$robust_tail_probability, 0.25)))
# Requirement 3: the aggregation approximation is recorded, not implied away.
stopifnot(identical(s_max$robustness_quantile_aggregation, "sum_of_per_cross_quantiles"))
stopifnot(grepl("not additive", s_max$robustness_quantile_aggregation_note, fixed = TRUE))

# Unconstrained, the maximize plan is exactly the n crosses with the LARGEST q0.25.
key <- function(df) sort(paste(df$parent1, df$parent2, sep = "__"))
expect_max <- rq_scores[order(rq_scores[[q_lo_col]], decreasing = TRUE)[1:8], ]
stopifnot(identical(key(plan_max), key(expect_max)))

# ---- (3) DIRECTION CORRECTNESS ---------------------------------------------
# Build a genuinely minimize-oriented ranked value: mean - i * SD, exactly what
# ng_run_cp_trait_value() returns for a decrease trait (disease, lodging). Lower is better,
# so the pessimistic case is the UPPER posterior tail.
i_int <- ng_selection_intensity(0.10)
min_scores <- do.call(ng_posterior_cross_predict, c(predict_args, list(
  robustness_quantile = rq,
  direction = "minimize",
  value_fun = function(sc) sc$cross_mean_gebv - i_int * sqrt(pmax(sc$pmv, 0))
)))
stopifnot(identical(attr(min_scores, "posterior")$direction, "minimize"))

min_args <- list(
  posterior_scores = min_scores, n_crosses = 8L, parent_kinship = parent_kinship,
  gain_col = gain_col, max_crosses_per_parent = Inf, lambda_group = 0,
  lambda_parent_use = 0, method = "greedy_local", robustness_quantile = rq
)
plan_min <- do.call(ng_optimize_robust_mating_plan, c(min_args, list(direction = "minimize")))
s_min <- attr(plan_min, "summary")
stopifnot(identical(s_min$robust_direction, "minimize"))
# Conservative tail for a minimize trait is 1 - q, i.e. the UPPER tail.
stopifnot(isTRUE(all.equal(s_min$robust_tail_probability, 0.75)))
stopifnot(identical(s_min$robust_quantile_source, q_hi_col))
stopifnot(isFALSE(s_min$robustness_quantile_is_normal_approximation))
stopifnot(isTRUE(s_min$robust_objective_is_negated))

# The minimize plan is exactly the n crosses with the SMALLEST upper-tail value: good at the
# pessimistic tail under a lower-is-better metric.
expect_min <- min_scores[order(min_scores[[q_hi_col]], decreasing = FALSE)[1:8], ]
stopifnot(identical(key(plan_min), key(expect_min)))

# Same table, maximize direction: a different plan, and demonstrably the anti-robust one for a
# minimize trait -- it is worse (higher) on the pessimistic upper tail than the minimize plan.
plan_min_wrong <- do.call(ng_optimize_robust_mating_plan, c(min_args, list(direction = "maximize")))
stopifnot(!identical(key(plan_min), key(plan_min_wrong)))
worst_ok   <- mean(plan_min[[q_hi_col]])
worst_bad  <- mean(plan_min_wrong[[q_hi_col]])
if (!(worst_ok < worst_bad)) {
  stop(sprintf("minimize-direction plan is not better at the pessimistic tail (%.6f vs %.6f)",
               worst_ok, worst_bad))
}
# The un-negated conservative value is reported on its native scale.
stopifnot(isTRUE(all.equal(s_min$robust_mean_value, worst_ok)))
stopifnot(isTRUE(all.equal(s_min$robust_total_value, sum(plan_min[[q_hi_col]]))))

# Top-N orientation follows direction too: under minimize the stable top-N is the N SMALLEST.
stopifnot(abs(sum(min_scores$posterior_topn_prob_6, na.rm = TRUE) - 6) < 0.05)
top6 <- order(min_scores$posterior_topn_prob_6, decreasing = TRUE)[1:6]
if (mean(min_scores$ranked_value_post_mean[top6]) >= mean(min_scores$ranked_value_post_mean)) {
  stop("minimize-direction posterior_topn_prob is not selecting the smallest ranked values")
}
# The top-N column bakes in an orientation, so a mismatched plan direction is refused.
topn_mismatch <- tryCatch({
  ng_optimize_robust_mating_plan(
    posterior_scores = min_scores, n_crosses = 6L, parent_kinship = parent_kinship,
    objective = "posterior_topn_prob", top_n_target = 6L, direction = "maximize",
    max_crosses_per_parent = 3L, method = "greedy_local")
  NA_character_
}, error = function(e) conditionMessage(e))
stopifnot(grepl("same direction", topn_mismatch, fixed = TRUE))

# ---- (4) the refusal still fires for a quantile that cannot be served exactly -
refusal <- tryCatch({
  do.call(ng_optimize_robust_mating_plan, c(plan_args, list(robustness_quantile = 0.40)))
  NA_character_
}, error = function(e) conditionMessage(e))
stopifnot(grepl("allow_normal_approximation", refusal, fixed = TRUE))
stopifnot(grepl("robustness_quantile = 0.4", refusal, fixed = TRUE))
# ... and the opt-in approximation is still available, still flagged.
plan_approx <- withCallingHandlers(
  do.call(ng_optimize_robust_mating_plan,
          c(plan_args, list(robustness_quantile = 0.40, allow_normal_approximation = TRUE))),
  warning = function(w) invokeRestart("muffleWarning")
)
s_approx <- attr(plan_approx, "summary")
stopifnot(isTRUE(s_approx$robustness_quantile_is_normal_approximation))
stopifnot(identical(s_approx$robust_quantile_source, "normal_approximation"))

# The default (NULL) still resolves to the exact cached lower CI tail, unchanged behaviour.
plan_default <- do.call(ng_optimize_robust_mating_plan, plan_args)
s_default <- attr(plan_default, "summary")
stopifnot(isFALSE(s_default$robustness_quantile_is_normal_approximation))
stopifnot(identical(s_default$robust_quantile_source, lo))
stopifnot(isTRUE(all.equal(s_default$robustness_quantile, (1 - 0.95) / 2)))

# ---- (5) orientation helper: pure-variance metrics stay higher-is-better ----
stopifnot(identical(ng_run_cp_value_orientation("decrease", "usefulness"), "minimize"))
stopifnot(identical(ng_run_cp_value_orientation("decrease", "mean"), "minimize"))
stopifnot(identical(ng_run_cp_value_orientation("increase", "usefulness"), "maximize"))
# A minimize trait ranked on family variance is still higher-is-better: more variance is more
# opportunity whichever way the trait points (ng_run_cp_trait_value applies no sign there).
stopifnot(identical(ng_run_cp_value_orientation("decrease", "pmv"), "maximize"))
stopifnot(identical(ng_run_cp_value_orientation("decrease", "vpm"), "maximize"))
stopifnot(identical(ng_run_cp_value_orientation("decrease", "parent_distance"), "maximize"))

cat("robust_quantile_direction: all checks passed\n")
cat(sprintf("  exact q=%.2f served from draws (%s); reported 95%% CI unchanged\n",
            rq, s_max$robust_quantile_source))
cat(sprintf("  minimize plan uses tail %.2f (%s); pessimistic-tail mean %.4f vs %.4f for the ",
            s_min$robust_tail_probability, s_min$robust_quantile_source, worst_ok, worst_bad))
cat("direction-blind plan\n")
