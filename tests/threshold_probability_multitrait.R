ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"),
            file.path("nextgen_cross_design", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  message("mvtnorm not installed; skipping multivariate threshold tests")
  quit(status = 0)
}

# ---- Test 1: reduction to univariate ---------------------------------------
# A 2-trait case where trait 2's interval is (-Inf, +Inf) must reduce exactly
# to the univariate result for trait 1.
mu_2  <- c(1.5, 0.0)
S_2   <- diag(c(0.25, 0.16))
lo_2  <- c(1.0, -Inf)
hi_2  <- c(3.0, +Inf)
k     <- 50L
p_mt  <- ng_p_superior_progeny_multitrait(mu_2, S_2, lo_2, hi_2, k)
p_one_univariate <- pnorm(3.0, 1.5, 0.5) - pnorm(1.0, 1.5, 0.5)
p_uni <- 1 - (1 - p_one_univariate)^k
gap <- abs(p_mt - p_uni)
if (gap > 1e-8) {
  stop(sprintf("multitrait did not reduce to univariate: mt=%.10f uni=%.10f gap=%g",
               p_mt, p_uni, gap))
}
cat(sprintf("Test 1 (reduction): mt=%.6f uni=%.6f gap=%g  OK\n",
            p_mt, p_uni, gap))

# ---- Test 2: upper = +Inf reduces to Pr(G >= tau) -------------------------
mu_1  <- 2.5
sd_1  <- 0.8
tau_1 <- 2.0
p_mt2 <- ng_p_superior_progeny_multitrait(
  mu = mu_1, Sigma_c = matrix(sd_1^2, 1, 1),
  tau_lower = tau_1, tau_upper = +Inf, k_progeny = k
)
p_st  <- ng_p_superior_progeny(mu = mu_1, sigma = sd_1, tau = tau_1, k_progeny = k)
gap2 <- abs(p_mt2 - p_st)
if (gap2 > 1e-8) {
  stop(sprintf("upper=+Inf did not reduce to ng_p_superior_progeny: mt=%.10f st=%.10f gap=%g",
               p_mt2, p_st, gap2))
}
cat(sprintf("Test 2 (upper=+Inf reduction): mt=%.6f st=%.6f gap=%g  OK\n",
            p_mt2, p_st, gap2))

# ---- Test 3: empty rectangle -> 0 -----------------------------------------
p_empty <- ng_p_superior_progeny_multitrait(
  mu = c(0, 0), Sigma_c = diag(2), tau_lower = c(1, 1),
  tau_upper = c(0.5, 0.5), k_progeny = 10
)
stopifnot(p_empty == 0)
cat("Test 3 (empty rectangle): OK\n")

# ---- Test 4: degenerate Sigma_c (zero variance on one trait) --------------
p_deg <- ng_p_superior_progeny_multitrait(
  mu = c(2.0, 1.5),
  Sigma_c = diag(c(0.25, 0)),
  tau_lower = c(1.5, 1.0),
  tau_upper = c(3.0, 2.0),
  k_progeny = 50
)
p_one_t1 <- pnorm(3.0, 2.0, 0.5) - pnorm(1.5, 2.0, 0.5)
p_t1_max <- 1 - (1 - p_one_t1)^50
stopifnot(abs(p_deg - p_t1_max) < 1e-8)
cat("Test 4 (degenerate Sigma diagonal): OK\n")

# ---- Test 5: Sigma_c constructor diagonal mode ----------------------------
v <- c(0.4, 0.9, 0.16)
S_diag <- ng_build_cross_trait_covariance(per_trait_var = v, G_hat = NULL)
stopifnot(is.matrix(S_diag), nrow(S_diag) == 3, ncol(S_diag) == 3)
stopifnot(all(abs(diag(S_diag) - v) < 1e-12))
stopifnot(all(abs(S_diag[upper.tri(S_diag)]) < 1e-12))
cat("Test 5 (Sigma_c diagonal mode): OK\n")

# ---- Test 6: Sigma_c constructor with G_hat -------------------------------
# Take a known correlation matrix R, set per-trait var v; expect Sigma_c
# off-diagonals = sqrt(v_i v_j) * R[i, j].
R <- matrix(c(1, 0.5, -0.3,
              0.5, 1, 0.2,
              -0.3, 0.2, 1), 3, 3)
S_full <- ng_build_cross_trait_covariance(per_trait_var = v, G_hat = R)
expected_off <- sqrt(v[1] * v[2]) * R[1, 2]
stopifnot(abs(S_full[1, 2] - expected_off) < 1e-12)
stopifnot(all(abs(diag(S_full) - v) < 1e-12))
ev <- eigen(S_full, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev) >= -1e-10)
cat("Test 6 (Sigma_c with G_hat): OK\n")

# ---- Test 7: non-PSD G_hat is rejected ------------------------------------
R_bad <- matrix(c(1, 0.9, 0.9,
                  0.9, 1, -0.9,
                  0.9, -0.9, 1), 3, 3)
bad_g_error <- tryCatch(
  ng_build_cross_trait_covariance(per_trait_var = v, G_hat = R_bad),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("not positive semidefinite", bad_g_error, fixed = TRUE))
cat("Test 7 (non-PSD G rejected): OK\n")

# ---- Test 5b: negative per_trait_var is rejected --------------------------
err <- tryCatch(
  ng_build_cross_trait_covariance(per_trait_var = c(0.4, -0.1, 0.16)),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("negative", err))
cat("Test 5b (negative variance rejected): OK\n")

# ---- Test 8: score-table augmentation -------------------------------------
scores_t <- data.frame(
  parent1 = c("P1", "P1", "P2"),
  parent2 = c("P2", "P3", "P3"),
  yield_mean   = c(2.5, 1.8, 3.0),
  yield_var    = c(0.4, 0.5, 0.3),
  disease_mean = c(1.2, 0.5, 0.3),
  disease_var  = c(0.16, 0.20, 0.12),
  stringsAsFactors = FALSE
)
trait_specs_t <- data.frame(
  trait    = c("yield", "disease"),
  mean_col = c("yield_mean", "disease_mean"),
  var_col  = c("yield_var",  "disease_var"),
  stringsAsFactors = FALSE
)
tau_lo_t <- c(2.0, -Inf)
tau_hi_t <- c(+Inf, 1.0)
out <- ng_add_p_superior_progeny_multitrait(
  scores = scores_t, trait_specs = trait_specs_t,
  tau_lower = tau_lo_t, tau_upper = tau_hi_t,
  k_progeny = 50L, G_hat = NULL
)
stopifnot("p_superior_progeny_mt" %in% names(out))
stopifnot(nrow(out) == nrow(scores_t))
stopifnot(all(out$p_superior_progeny_mt >= 0 - 1e-9))
stopifnot(all(out$p_superior_progeny_mt <= 1 + 1e-9))
# Cross 3 has highest yield mean (3.0) AND lowest disease mean (0.3): both
# directions favored (yield >= 2.0, disease <= 1.0), so it should have the
# highest probability of all three.
stopifnot(out$p_superior_progeny_mt[3] >= out$p_superior_progeny_mt[2])
# Cross 3 should also dominate Cross 1 (Cross 1's yield mean is 2.5, just
# above the 2.0 threshold; Cross 3 has yield mean 3.0 and a much lower
# disease mean). A future regression that flips trait direction would
# pass the Cross 3 vs Cross 2 check (Cross 2 has yield 1.8, below the
# threshold) but fail this stronger assertion.
stopifnot(out$p_superior_progeny_mt[3] >= out$p_superior_progeny_mt[1])
cat(sprintf("Test 8 (augmentation): p_mt = %s  OK\n",
            paste(sprintf("%.4f", out$p_superior_progeny_mt), collapse = ", ")))

# ---- Test 8b: row-level result matches direct call -------------------------
S3 <- ng_build_cross_trait_covariance(per_trait_var = c(0.3, 0.12), G_hat = NULL)
p3_direct <- ng_p_superior_progeny_multitrait(
  mu = c(3.0, 0.3), Sigma_c = S3,
  tau_lower = tau_lo_t, tau_upper = tau_hi_t, k_progeny = 50
)
stopifnot(abs(out$p_superior_progeny_mt[3] - p3_direct) < 1e-12)
cat("Test 8b (per-row matches direct call): OK\n")

# ---- Test 8c: missing column in scores -------------------------------------
err <- tryCatch(
  ng_add_p_superior_progeny_multitrait(
    scores = scores_t[, c("parent1", "parent2", "yield_mean", "yield_var")],
    trait_specs = trait_specs_t,
    tau_lower = tau_lo_t, tau_upper = tau_hi_t, k_progeny = 50L
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("missing columns", err))
cat("Test 8c (missing column rejected): OK\n")

# ---- Test 8d: length mismatch ---------------------------------------------
err2 <- tryCatch(
  ng_add_p_superior_progeny_multitrait(
    scores = scores_t, trait_specs = trait_specs_t,
    tau_lower = c(2.0), tau_upper = c(+Inf, 1.0), k_progeny = 50L
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("one entry per", err2))
cat("Test 8d (length mismatch rejected): OK\n")

# ---- Test 8e: negative variance in scores propagates -----------------------
# The augmentation should not silently floor negative variances; T3's
# constructor guard must fire from this path.
scores_bad <- scores_t
scores_bad$yield_var[2] <- -0.1
err3 <- tryCatch(
  ng_add_p_superior_progeny_multitrait(
    scores = scores_bad, trait_specs = trait_specs_t,
    tau_lower = tau_lo_t, tau_upper = tau_hi_t, k_progeny = 50L
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("negative", err3))
cat("Test 8e (negative variance in scores rejected): OK\n")

# ---- Test 9: posterior multi-trait threshold CI ---------------------------
# Verifies ng_posterior_multitrait_cross_predict gains the four tau args and
# emits p_superior_progeny_mt_post_mean/_lower/_upper.
set.seed(20260525L)
n9 <- 24L; m9 <- 300L
ids9 <- paste0("L", seq_len(n9))
markers9 <- paste0("M", seq_len(m9))
geno9 <- 2L * matrix(rbinom(n9 * m9, 1L, 0.4), n9, m9,
                     dimnames = list(ids9, markers9))
marker_map9 <- data.frame(
  marker = markers9,
  chr = rep(1:3, length.out = m9),
  pos_cm = rep(seq(0, 100, length.out = m9 / 3L), 3L),
  stringsAsFactors = FALSE
)
Y9 <- cbind(
  yield   = as.numeric(geno9 %*% rnorm(m9, sd = 0.05) + rnorm(n9, sd = 0.4)),
  disease = as.numeric(geno9 %*% rnorm(m9, sd = 0.04) + rnorm(n9, sd = 0.4))
)
rownames(Y9) <- ids9

traits9 <- ng_multitrait_spec(
  trait = c("yield", "disease"),
  direction = c("maximize", "minimize"),
  economic_weight = c(1, 1)
)

# Use only a small subset of crosses so the test stays fast.
pair_ids9 <- ids9[seq_len(10L)]
pairs9 <- ng_make_pairs(pair_ids9)
geno9_use <- geno9[pair_ids9, , drop = FALSE]
Y9_use <- Y9[pair_ids9, , drop = FALSE]

# 0.29.0: this fixture used genetic_covariance = diag(c(yield = 1, disease = 1)) -- unit
# genetic variance against phenotypes whose variance is ~0.92 and ~0.97, i.e. implied
# h2 = 1.083 and 1.026. ng_posterior_multitrait_cross_predict() pairs the caller's G with a
# P estimated from Y when no P is given, so the pair is now refused by the P - G guard, and
# rightly: it was never a possible P = G + R. The threshold property under test does not
# depend on G's scale, so take G as half the phenotypic covariance of the same rows the run
# uses -- a uniform h2 of 0.5, with R = P - G = 0.5 P positive definite by construction.
G9 <- 0.5 * as.matrix(ng_estimate_phenotypic_covariance(Y9_use, shrinkage = "auto"))
dimnames(G9) <- list(c("yield", "disease"), c("yield", "disease"))

post9 <- ng_posterior_multitrait_cross_predict(
  geno = geno9_use, Y = Y9_use, traits = traits9,
  ids = pair_ids9, pairs = pairs9,
  marker_map = marker_map9,
  n_draws = 30L, posterior_method = "closed_form",
  genetic_covariance_method = "beta_posterior",
  genetic_covariance = G9,
  index_method = "economic_index", value_mode = "mean",
  use_cpp = FALSE, seed = 7L,
  tau_lower_vec = c(0, -Inf),
  tau_upper_vec = c(+Inf, 0),
  threshold_k_progeny = 50L
)

required_cols <- c("p_superior_progeny_mt_post_mean",
                   "p_superior_progeny_mt_post_lower",
                   "p_superior_progeny_mt_post_upper")
missing_cols <- setdiff(required_cols, names(post9))
if (length(missing_cols)) {
  stop("posterior_multitrait_cross_predict missing columns: ",
       paste(missing_cols, collapse = ", "))
}
stopifnot(all(post9$p_superior_progeny_mt_post_lower >= 0 - 1e-9))
stopifnot(all(post9$p_superior_progeny_mt_post_upper <= 1 + 1e-9))
stopifnot(all(post9$p_superior_progeny_mt_post_lower <=
              post9$p_superior_progeny_mt_post_upper + 1e-9))
stopifnot(all(post9$p_superior_progeny_mt_post_lower <=
              post9$p_superior_progeny_mt_post_mean + 1e-9))
stopifnot(all(post9$p_superior_progeny_mt_post_mean <=
              post9$p_superior_progeny_mt_post_upper + 1e-9))

# Regression guard: posterior CI should have non-trivial width somewhere
# (catches a future saturation collapse where all draws return identical
# values).
ci_widths <- post9$p_superior_progeny_mt_post_upper - post9$p_superior_progeny_mt_post_lower
stopifnot(max(ci_widths, na.rm = TRUE) > 1e-3)

# Verify the posterior metadata records the threshold settings.
meta9 <- attr(post9, "posterior_multitrait")
stopifnot(!is.null(meta9$tau_lower_vec))
stopifnot(!is.null(meta9$tau_upper_vec))
stopifnot(identical(meta9$threshold_k_progeny, 50L))

cat(sprintf("Test 9 (posterior CI): %d crosses, mean range [%.4f, %.4f]  OK\n",
            nrow(post9),
            min(post9$p_superior_progeny_mt_post_mean, na.rm = TRUE),
            max(post9$p_superior_progeny_mt_post_mean, na.rm = TRUE)))

# ---- Test 10: empirical agreement vs simulated progeny (k=1, n=5000) ------
# For one synthetic cross, draw 5000 progeny from N(mu, Sigma_c) and compare
# the empirical fraction satisfying A to the predicted Pr.
set.seed(42L)
mu_e    <- c(2.5, 1.0)
S_e     <- matrix(c(0.40, 0.10,
                    0.10, 0.16), 2, 2)
tau_lo_e <- c(2.0, -Inf)
tau_hi_e <- c(+Inf, 1.2)
k_e      <- 1L
n_sim    <- 5000L
chol_S   <- chol(S_e)
z        <- matrix(rnorm(2L * n_sim), nrow = 2L)
sims     <- t(chol_S) %*% z + mu_e
in_A     <- colSums(sims >= tau_lo_e & sims <= tau_hi_e) == 2L
empirical <- mean(in_A)
predicted <- ng_p_superior_progeny_multitrait(
  mu = mu_e, Sigma_c = S_e,
  tau_lower = tau_lo_e, tau_upper = tau_hi_e,
  k_progeny = k_e
)
gap <- abs(empirical - predicted)
if (gap > 0.02) {
  stop(sprintf("empirical/predicted gap too large: empirical=%.4f predicted=%.4f gap=%.4f",
               empirical, predicted, gap))
}
cat(sprintf("Test 10 (empirical agreement, n=%d): empirical=%.4f predicted=%.4f gap=%.4f  OK\n",
            n_sim, empirical, predicted, gap))

# ---- Test 11: empirical agreement at k=20 (order-statistic regime) --------
n_families <- 1000L
hits <- logical(n_families)
for (f in seq_len(n_families)) {
  z <- matrix(rnorm(2L * 20L), nrow = 2L)
  sims <- t(chol_S) %*% z + mu_e
  hits[f] <- any(colSums(sims >= tau_lo_e & sims <= tau_hi_e) == 2L)
}
empirical_k20 <- mean(hits)
predicted_k20 <- ng_p_superior_progeny_multitrait(
  mu = mu_e, Sigma_c = S_e,
  tau_lower = tau_lo_e, tau_upper = tau_hi_e,
  k_progeny = 20L
)
gap_k20 <- abs(empirical_k20 - predicted_k20)
if (gap_k20 > 0.05) {
  stop(sprintf("k=20 empirical gap too large (families=%d): empirical=%.4f predicted=%.4f gap=%.4f",
               n_families, empirical_k20, predicted_k20, gap_k20))
}
cat(sprintf("Test 11 (k=20, families=%d): empirical=%.4f predicted=%.4f gap=%.4f  OK\n",
            n_families, empirical_k20, predicted_k20, gap_k20))

# ---- Test 12: invalid G_hat is rejected (Fix 1) ---------------------------
G_bad_cov <- matrix(c(0.04, 0.05,
                      0.05, 0.04), 2, 2)  # |G[1,2]| > sqrt(0.04*0.04) = 0.04
err_g <- tryCatch(
  ng_build_cross_trait_covariance(per_trait_var = c(1, 1), G_hat = G_bad_cov),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("not a valid covariance", err_g))
cat("Test 12 (invalid G_hat covariance rejected): OK\n")

# ---- Test 13: fractional k_progeny < 1 is rejected (Fix 2) ---------------
err_k <- tryCatch(
  ng_p_superior_progeny_multitrait(
    mu = c(0), Sigma_c = matrix(1, 1, 1),
    tau_lower = -1, tau_upper = 1, k_progeny = 0.5
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("k_progeny", err_k))
cat("Test 13 (fractional k_progeny rejected): OK\n")

# ---- Test 14: empty mu vector is rejected (Fix 3) ------------------------
err_mu <- tryCatch(
  ng_p_superior_progeny_multitrait(
    mu = numeric(0), Sigma_c = matrix(0, 0, 0),
    tau_lower = numeric(0), tau_upper = numeric(0), k_progeny = 10
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("length >= 1", err_mu))
cat("Test 14 (empty mu rejected): OK\n")
