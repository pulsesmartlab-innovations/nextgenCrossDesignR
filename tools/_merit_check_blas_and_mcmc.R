# Merit checks for the remaining two ports:
#   (a) BLAS-backed full-posterior kernel: would calling dgemv from C++
#       meaningfully beat the current plain-nested-loop C++ version?
#   (b) MCMC Gibbs sampler: how slow is it at realistic scale, and does
#       porting it matter?
ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"),
            file.path("nextgen_cross_design", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(20260524L)

# ---- (a) BLAS-backed full-posterior merit ---------------------------------
# Time a single matrix-vector multiply (R's BLAS dgemv) at m=1500. The
# current C++ kernel does 2 mat-vec per pair as plain nested loops with
# the column-skip optimization. If R's BLAS dgemv is much faster than the
# current C++ loop on a SINGLE call, the BLAS-backed port would help.
m <- 1500L
R_mat <- matrix(rnorm(m * m), m, m); R_mat <- (R_mat + t(R_mat)) / 2; diag(R_mat) <- 1
a <- rnorm(m)

# Time R's %*% (= dgemv) for 1000 multiplies.
t_blas <- system.time(replicate(1000L, R_mat %*% a))["elapsed"]
cat(sprintf("R's BLAS dgemv at m=%d: %.4fs / 1000 calls = %.2f us/call\n",
            m, t_blas, t_blas * 1e6 / 1000))

# Compare to current full-posterior C++ kernel: it does 2 mat-vec PER PAIR.
# So for n_pairs candidate pairs, current C++ does 2 * n_pairs mat-vec ops.
# If BLAS dgemv per call is x microseconds, the BLAS ceiling for the full-
# posterior at n_pairs pairs is 2 * n_pairs * x microseconds.
n_pairs_ref <- 1770L  # matches the speedup-bench fixture (n=60)
blas_ceiling_seconds <- 2 * n_pairs_ref * t_blas * 1e-3   # since t_blas is for 1000 calls
cat(sprintf("BLAS-ceiling estimate for full-posterior at n_pairs=%d, m=%d: ~%.2fs\n",
            n_pairs_ref, m, blas_ceiling_seconds))

# The current C++ kernel takes ~17s on this fixture; if BLAS-ceiling is e.g.
# 4s, BLAS-backed port would give ~4x improvement = worthwhile. If BLAS-
# ceiling is ~10s, the gain is marginal.
cat("\n## Merit verdict for BLAS-backed full-posterior\n")
current_cpp_estimate <- 17.0  # observed in prior microbench
gain_factor <- current_cpp_estimate / max(blas_ceiling_seconds, 1e-3)
if (gain_factor < 2) {
  cat(sprintf("VERDICT: BLAS ceiling = %.2fs, current C++ = %.2fs, gain = %.1fx — SKIP\n",
              blas_ceiling_seconds, current_cpp_estimate, gain_factor))
  cat("R loop overhead, not the mat-vec, was the bottleneck. C++ already captured it.\n")
} else if (gain_factor < 5) {
  cat(sprintf("VERDICT: BLAS ceiling = %.2fs, current C++ = %.2fs, gain = %.1fx — moderate\n",
              blas_ceiling_seconds, current_cpp_estimate, gain_factor))
  cat("Port if posterior cross prediction is a common workflow; defer otherwise.\n")
} else {
  cat(sprintf("VERDICT: BLAS ceiling = %.2fs, current C++ = %.2fs, gain = %.1fx — PORT\n",
              blas_ceiling_seconds, current_cpp_estimate, gain_factor))
}

# ---- (b) MCMC sampler merit ------------------------------------------------
# Time the R reference for a small but realistic MCMC run. If it's seconds,
# port has low merit. If it's minutes/hours at scale, port has high merit.
n <- 60L
geno <- 2L * matrix(rbinom(n * m, 1L, 0.45), n, m)
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n, sd = 0.5))
ids <- paste0("P", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))
names(y) <- ids
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 5L, seed = 1L)
Xc <- sweep(geno, 2L, fit$marker_mean, "-")
storage.mode(Xc) <- "double"
yc <- y - mean(y)

# Time MCMC: 100 draws + 100 burnin (small run)
cat(sprintf("\nMCMC merit fixture: n=%d, m=%d. Running 100 burnin + 100 draws (R reference).\n", n, m))
t_mcmc <- system.time(
  ng_sample_ridge_posterior_mcmc(X = Xc, yc = yc, n_draws = 100L,
                                 burnin = 100L, thin = 1L,
                                 lambda_init = fit$lambda,
                                 sigma_e2_init = fit$sigma_e2, seed = 1L)
)["elapsed"]
cat(sprintf("MCMC R reference: %.2fs for 200 total iters\n", t_mcmc))
cat(sprintf("Per iteration: %.3fs\n", t_mcmc / 200))

# Project to realistic run: 500 draws + 500 burnin = 1000 iter
t_full <- t_mcmc * 1000 / 200
cat(sprintf("Projected R time for 500 burnin + 500 draws (typical): %.1fs\n", t_full))

# C++ port expected speedup: each iter is dominated by:
#   - A = XX' + lambda I    (BLAS dsyrk; doesn't benefit from C++ port)
#   - chol(A)               (LAPACK dpotrf; doesn't benefit)
#   - BCM beta draw         (ALREADY in C++ via ng_bcm_posterior_sampler_cpp)
#   - sigma_e2 / sigma_beta2 draws (trivial)
# The BCM draw was 1.9x in our microbench. Most of MCMC's per-iter cost is
# BLAS dsyrk + chol which are already in LAPACK. Expected MCMC speedup: ~2-3x.
cat("\n## Merit verdict for MCMC port\n")
if (t_full < 10) {
  cat(sprintf("VERDICT: full MCMC run only %.1fs — SKIP. R version is already fast enough.\n",
              t_full))
} else if (t_full < 120) {
  cat(sprintf("VERDICT: full MCMC run ~%.0fs — moderate merit. Port if MCMC is common; defer if rare.\n",
              t_full))
} else {
  cat(sprintf("VERDICT: full MCMC run ~%.0fs (~%.1f min) — PORT. Multi-minute waits hurt usability.\n",
              t_full, t_full / 60))
}
