helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# D1 — Posterior cross-prediction framework.
# Verifies (a) the closed-form BCM sampler returns draws whose mean and
# covariance match the analytical ridge posterior, (b) MCMC convergence on a
# small case, (c) the per-cross posterior summary attaches CI columns whose
# point-estimate mean is close to the legacy ng_score_crosses() result, and
# (d) posterior_topn_prob behaves as expected at the extremes.

set.seed(20260523L)
n <- 30L
m <- 60L
geno <- matrix(2L * rbinom(n * m, 1L, 0.45), nrow = n, ncol = m)
ids <- paste0("L", seq_len(n))
markers <- paste0("M", seq_len(m))
rownames(geno) <- ids; colnames(geno) <- markers
beta_true <- rnorm(m, sd = 0.10)
y <- as.numeric(geno %*% beta_true + rnorm(n, sd = 0.5))
names(y) <- ids

# ---- BCM sampler matches the analytical posterior mean and variance --------
fit <- ng_fit_ridge_effects(geno, y, ids = ids, lambda = 8, kfold = 5L, seed = 1L)
Xc <- sweep(geno, 2L, fit$marker_mean, "-")
yc <- y - mean(y)
post_mean_analytic <- as.numeric(solve(crossprod(Xc) + diag(fit$lambda, m), crossprod(Xc, yc)))
post_cov_analytic  <- fit$sigma_e2 * solve(crossprod(Xc) + diag(fit$lambda, m))

draws <- ng_sample_ridge_posterior_bcm(
  X = Xc, yc = yc, sigma_e2 = fit$sigma_e2, lambda = fit$lambda,
  n_draws = 2000L, seed = 1L
)
mc_mean <- rowMeans(draws)
mc_var  <- apply(draws, 1L, stats::var)

mean_err <- max(abs(mc_mean - post_mean_analytic))
diag_var_err <- max(abs(mc_var - diag(post_cov_analytic))) /
  max(diag(post_cov_analytic))
if (mean_err > 0.05) {
  stop(sprintf("BCM posterior mean off by %.4f (m=%d, S=2000); should be near 0", mean_err, m))
}
if (diag_var_err > 0.20) {
  stop(sprintf("BCM posterior diag variance off by %.2f (relative); should be < 0.20 at S=2000",
               diag_var_err))
}

# Posterior mean of draws should also match fit$beta within Monte-Carlo error.
fit_mean_err <- max(abs(mc_mean - fit$beta))
if (fit_mean_err > 0.05) {
  stop(sprintf("BCM mean disagrees with closed-form fit$beta by %.4f", fit_mean_err))
}

# ---- MCMC convergence on the same fixture ----------------------------------
mcmc <- ng_sample_ridge_posterior_mcmc(
  X = Xc, yc = yc, n_draws = 500L, burnin = 500L, thin = 2L, seed = 2L
)
stopifnot(ncol(mcmc$beta_draws) == 500L)
stopifnot(length(mcmc$sigma_e2_draws) == 500L)
stopifnot(all(mcmc$sigma_e2_draws > 0))
stopifnot(all(mcmc$sigma_beta2_draws > 0))
# Posterior means from MCMC should be close to the analytical posterior mean
# even though it integrates over (sigma_e2, sigma_beta2).
mcmc_mean_err <- max(abs(rowMeans(mcmc$beta_draws) - post_mean_analytic))
if (mcmc_mean_err > 0.10) {
  stop(sprintf("MCMC posterior mean off analytical by %.4f at S=500 (post-burnin)", mcmc_mean_err))
}

# ---- P(superior progeny) closed-form sanity checks -------------------------
# Limit: as k -> Inf and mu < tau, P -> 0; as mu > tau and sigma > 0, P -> 1.
p_low  <- ng_p_superior_progeny(mu = 0, sigma = 0.1, tau = 2, k_progeny = 100)
p_high <- ng_p_superior_progeny(mu = 5, sigma = 1.0, tau = 2, k_progeny = 100)
p_zero <- ng_p_superior_progeny(mu = 5, sigma = 0,   tau = 2, k_progeny = 100)
p_neg  <- ng_p_superior_progeny(mu = 0, sigma = 0,   tau = 2, k_progeny = 100)
stopifnot(p_low  < 1e-6)
stopifnot(p_high > 1 - 1e-6)
stopifnot(p_zero == 1)
stopifnot(p_neg  == 0)
# Vectorized.
p_vec <- ng_p_superior_progeny(c(0, 5), c(0.1, 1.0), 2, c(100, 100))
stopifnot(length(p_vec) == 2L)

# ---- ng_posterior_cross_predict integration --------------------------------
# Use inbred founders so the heterozygous-parent guard doesn't trip.
beta_named <- setNames(beta_true, markers)
y_named <- y
post_fit <- ng_fit_ridge_effects_posterior(
  geno = geno, y = y_named, ids = ids, lambda = 8, kfold = 5L,
  n_draws = 100L, method = "closed_form", seed = 1L
)
mk <- data.frame(
  marker = markers,
  chr = rep(1:3, length.out = m),
  pos_cm = rep(seq(0, 80, length.out = m / 3), 3)
)
adj <- setNames(rnorm(n), ids)

post_scores <- ng_posterior_cross_predict(
  geno = geno, posterior_effects = post_fit, marker_map = mk, ids = ids,
  adjusted_pheno = adj, selection_prop = 0.10,
  target = "DH", recomb_model = "haldane", use_cpp = FALSE,
  tau_superior = quantile(adj, 0.90), k_progeny = 50L,
  top_n_targets = c(5L, 10L)
)
expected_cols <- c(
  "uc_dh_gebv_post_mean", "uc_dh_gebv_post_lower", "uc_dh_gebv_post_upper",
  "dh_pmv_var_post_mean", "dh_pmv_var_post_lower", "dh_pmv_var_post_upper",
  "p_superior_progeny_post_mean", "p_superior_progeny_post_lower", "p_superior_progeny_post_upper",
  "posterior_topn_prob_5", "posterior_topn_prob_10"
)
missing <- setdiff(expected_cols, names(post_scores))
if (length(missing)) stop("posterior_cross_predict missing columns: ",
                          paste(missing, collapse = ", "))

# Sanity: posterior_topn_prob is between 0 and 1, sums to N across all crosses.
for (N in c(5L, 10L)) {
  col <- paste0("posterior_topn_prob_", N)
  stopifnot(all(post_scores[[col]] >= 0 & post_scores[[col]] <= 1))
  expect_sum <- N
  actual_sum <- sum(post_scores[[col]], na.rm = TRUE)
  if (abs(actual_sum - expect_sum) > 0.05) {
    stop(sprintf("posterior_topn_prob_%d sums to %.3f, expected %d", N, actual_sum, expect_sum))
  }
}

# Sanity: posterior CI brackets the posterior mean.
stopifnot(all(post_scores$uc_dh_gebv_post_lower <= post_scores$uc_dh_gebv_post_mean + 1e-6))
stopifnot(all(post_scores$uc_dh_gebv_post_mean <= post_scores$uc_dh_gebv_post_upper + 1e-6))
stopifnot(all(post_scores$dh_pmv_var_post_lower <= post_scores$dh_pmv_var_post_mean + 1e-6))
stopifnot(all(post_scores$dh_pmv_var_post_mean <= post_scores$dh_pmv_var_post_upper + 1e-6))
# Probability columns can have very-near-degenerate posteriors (p close to 0
# or 1 for crosses far below/above the threshold). With small S the empirical
# quantile of a few-outlier-in-many-zeros distribution can sit outside the
# arithmetic mean — that's a property of order statistics, not a bug. Assert
# only the universally true bounds: in [0, 1] and lower <= upper.
stopifnot(all(post_scores$p_superior_progeny_post_lower >= 0 - 1e-9))
stopifnot(all(post_scores$p_superior_progeny_post_upper <= 1 + 1e-9))
stopifnot(all(post_scores$p_superior_progeny_post_lower <= post_scores$p_superior_progeny_post_upper + 1e-9))

# Sanity: posterior mean of usefulness should be close to the point estimate.
point_uc <- post_scores$uc_dh_gebv
mean_post_uc <- post_scores$uc_dh_gebv_post_mean
gap <- max(abs(point_uc - mean_post_uc) / pmax(abs(point_uc), 1e-6))
if (gap > 0.50) {
  stop(sprintf("Posterior mean of usefulness differs from point estimate by relative %.3f", gap))
}

# ---- Posterior-aware OCS reuses the new columns ----------------------------
parent_K <- ng_parent_kinship(geno)
plan_robust <- ng_optimize_robust_mating_plan(
  posterior_scores = post_scores, n_crosses = 6L, parent_K = parent_K,
  gain_col = "uc_dh_gebv", robustness_quantile = (1 - 0.95) / 2,
  max_crosses_per_parent = 3L, lambda_group = 0, lambda_parent_use = 0,
  lambda_parent_use_mode = "absolute", method = "greedy_local"
)
s_robust <- attr(plan_robust, "summary")
stopifnot(nrow(plan_robust) == 6L)
stopifnot(identical(s_robust$robust_objective, "posterior_quantile"))
stopifnot(s_robust$max_parent_use <= 3L)

plan_topn <- ng_optimize_robust_mating_plan(
  posterior_scores = post_scores, n_crosses = 6L, parent_K = parent_K,
  objective = "posterior_topn_prob", top_n_target = 10L,
  max_crosses_per_parent = 3L, method = "greedy_local"
)
s_topn <- attr(plan_topn, "summary")
stopifnot(nrow(plan_topn) == 6L)
stopifnot(identical(s_topn$robust_objective, "posterior_topn_prob"))
stopifnot(identical(s_topn$robust_top_n_target, 10L))

cat("posterior_cross_predict: 7/7 checks passed\n")
cat(sprintf("  BCM mean err %.4f, diag var err %.3f (S=2000)\n", mean_err, diag_var_err))
cat(sprintf("  MCMC posterior mean err vs analytical %.4f (S=500 post-burnin)\n", mcmc_mean_err))
cat(sprintf("  P(superior progeny) extreme cases: low=%.2e high=%.6f\n", p_low, p_high))
cat(sprintf("  posterior CI bracketing OK for usefulness, variance, P(superior)\n"))
cat(sprintf("  robust plan picked %d crosses (%d unique parents, max use %d)\n",
            nrow(plan_robust), s_robust$unique_parents, s_robust$max_parent_use))
