helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# T1.3: ridge posterior marker-effect variance must equal
#   sigma_e2 * diag((X'X + lambda I)^{-1}), NOT sigma_e2 / (diag(X'X) + lambda).
# Under LD between markers the two differ; the latter underestimates posterior variance.

# ---- Test 1: correlated markers, compared against explicit primal inverse ---------------
set.seed(2026L)
n <- 30L
m <- 8L
Z <- matrix(rnorm(n * 3L), nrow = n)
loadings <- matrix(rnorm(3L * m, sd = 0.6), nrow = 3L)
noise <- matrix(rnorm(n * m, sd = 0.2), nrow = n)
geno <- Z %*% loadings + noise          # markers correlated through 3 latent factors
ids <- paste0("L", seq_len(n))
markers <- paste0("M", seq_len(m))
rownames(geno) <- ids; colnames(geno) <- markers
y <- as.numeric(geno %*% rnorm(m, sd = 0.4) + rnorm(n, sd = 0.5))
names(y) <- ids

lambda <- 7
fit <- ng_fit_ridge_effects(geno, y, ids = ids, lambda = lambda, kfold = 5L, seed = 2026L)

Xc <- sweep(geno, 2L, colMeans(geno), "-")
primal_diag <- diag(solve(crossprod(Xc) + diag(lambda, m)))
expected <- fit$sigma_e2 * primal_diag

stopifnot(length(fit$beta_var) == m)
err <- max(abs(fit$beta_var - expected))
if (err > 1e-8) stop(sprintf("ridge beta_var deviates from primal diag(inv): max abs err = %g", err))

# ---- Test 2: new formula must be >= old (broken) formula, with strict > under LD --------
old_beta_var <- fit$sigma_e2 / (colSums(Xc * Xc) + lambda)
stopifnot(all(fit$beta_var + 1e-12 >= old_beta_var))
ratio <- mean(fit$beta_var / old_beta_var)
if (ratio < 1.05) {
  stop(sprintf("expected new beta_var to noticeably exceed old under LD; mean ratio = %.3f", ratio))
}

# ---- Test 3: under heavy regularization both formulas approach sigma_e2 / lambda ---------
# As lambda -> Inf, A^{-1} -> 0, so x_k' A^{-1} x_k -> 0 and the dual identity yields
# beta_var_k -> sigma_e2 / lambda. The broken formula also -> sigma_e2 / lambda, so the
# gap closes asymptotically; we use this as an upper-bound sanity check.
big_lambda <- 1e6
fit3 <- ng_fit_ridge_effects(geno, y, ids = ids, lambda = big_lambda, kfold = 5L, seed = 99L)
asymptotic <- fit3$sigma_e2 / big_lambda
if (max(abs(fit3$beta_var - asymptotic)) / asymptotic > 1e-3) {
  stop("under heavy regularization, beta_var should approach sigma_e2 / lambda")
}

# ---- Test 4: df_eff must equal trace((X'X+lambda I)^{-1} X'X), in (0, m] -----------------
# Sanity-check the identity sigma_e2 * (n - df_eff) = sum(resid^2) using the public
# predict path so we exercise the same intercept handling as the fitter.
fitted_pkg <- ng_predict_gebv(geno, fit)
resid_pkg <- y - fitted_pkg
rss <- sum(resid_pkg * resid_pkg)
inferred_df <- length(y) - rss / fit$sigma_e2
df_eff_primal <- sum(diag(solve(crossprod(Xc) + diag(lambda, m), crossprod(Xc))))
stopifnot(inferred_df > 0)
stopifnot(inferred_df <= m + 1e-6)
if (abs(inferred_df - df_eff_primal) > 1e-6) {
  stop(sprintf("df_eff mismatch: inferred=%.6f primal=%.6f", inferred_df, df_eff_primal))
}

cat("ridge_beta_var: 4/4 checks passed\n")
cat(sprintf("  mean ratio new/old beta_var under LD: %.3f\n", ratio))
cat(sprintf("  effective df: %.3f (n=%d, m=%d)\n", inferred_df, n, m))
cat(sprintf("  asymptotic beta_var at lambda=1e6: %.3e\n", asymptotic))
