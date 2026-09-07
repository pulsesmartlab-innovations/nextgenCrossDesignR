helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 0.29.0 -- R/31_genetic_covariance.R must not depend on the COLUMN ORDER of Y.
#
# This is the same defect tests/trait_order_invariance.R (0.28.0) asserts for the
# runner and the multi-trait posterior, applied to the four remaining
# position-derived seed sites in R/31:
#
#   ng_genetic_cov_two_stage_ridge()          seed + j  -> lambda CV fold split
#   ng_posterior_genetic_covariance()
#     beta_posterior      per-trait draws     seed + j  -> posterior innovations
#     parametric_bootstrap per-trait fit      seed + j  -> lambda CV fold split
#     parametric_bootstrap residual draws     one shared stream consumed in
#                                             COLUMN ORDER
#
# The assertion, exactly as in trait_order_invariance.R: permute the columns of
# Y, permute the answer back by NAME, and require the same numbers. Everything
# the seeds control is required to match at tolerance = 0.

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }

set.seed(31029L)
n <- 60L
m <- 150L
geno <- matrix(stats::rbinom(n * m, 2L, 0.3), n, m,
               dimnames = list(paste0("L", seq_len(n)), paste0("M", seq_len(m))))
b1 <- stats::rnorm(m, sd = 0.20)
b2 <- 0.5 * b1 + stats::rnorm(m, sd = 0.20)
b3 <- stats::rnorm(m, sd = 0.30)
g <- scale(geno %*% cbind(b1, b2, b3))
Y <- g + matrix(stats::rnorm(n * 3L), n, 3L)
# Deliberately NOT in alphabetical order, so a by-name reorder is a real reorder.
colnames(Y) <- c("yield", "protein", "disease")
rownames(Y) <- rownames(geno)
tn <- colnames(Y)
perm <- c(3L, 1L, 2L)              # disease, yield, protein
Y_perm <- Y[, perm, drop = FALSE]

back <- function(M) M[tn, tn, drop = FALSE]

# ============================================================================
# 1. ng_genetic_cov_two_stage_ridge(): everything the seeds control is IDENTICAL
# ============================================================================
run_ts <- function(Ym) suppressWarnings(ng_genetic_cov_two_stage_ridge(
  geno = geno, Y = Ym, kfold = 5L, seed = 11L, return_diagnostics = TRUE))
A <- run_ts(Y)
B <- run_ts(Y_perm)

stopifnot(identical(as.numeric(A$lambdas[tn]), as.numeric(B$lambdas[tn])))
stopifnot(identical(as.numeric(A$cv_predictive_r2[tn]), as.numeric(B$cv_predictive_r2[tn])))
stopifnot(identical(as.numeric(A$sigma_e2[tn]), as.numeric(B$sigma_e2[tn])))
stopifnot(identical(as.numeric(A$sigma_g2[tn]), as.numeric(B$sigma_g2[tn])))
ok("two_stage_ridge: per-trait lambda, CV r2, sigma_e2 and sigma_g2 identical at tolerance = 0")

stopifnot(identical(as.numeric(A$beta_hat[, tn]), as.numeric(B$beta_hat[, tn])))
stopifnot(identical(as.numeric(back(A$R_beta)), as.numeric(back(B$R_beta))))
stopifnot(identical(as.numeric(back(A$G_raw)), as.numeric(back(B$G_raw))))
ok("two_stage_ridge: marker effects, beta correlation matrix and the raw G identical at tolerance = 0")

# The only quantity that is NOT bit-identical is the PSD projection of G_raw. It
# is a LAPACK eigendecomposition of a permuted matrix, so the difference is
# floating-point reassociation inside eigen()/nearPD, not seed dependence -- and
# it is bounded far below anything a breeder could see. Documented, not hidden.
gap <- max(abs(back(A$G_hat) - back(B$G_hat)))
rel <- gap / max(1, max(abs(A$G_hat)))
cat(sprintf("  G_hat (after the PSD projection) differs by %.3e absolute, %.3e relative\n",
            gap, rel))
stopifnot(rel < 1e-12)
ok("two_stage_ridge: G_hat agrees to LAPACK reassociation noise (< 1e-12 relative) after the PSD projection")

# ...and through the exported entry point.
run_pub <- function(Ym) suppressWarnings(ng_estimate_genetic_covariance(
  geno = geno, Y = Ym, method = "two_stage_ridge", kfold = 5L, seed = 11L))
GA <- run_pub(Y)
GB <- run_pub(Y_perm)
stopifnot(identical(dimnames(GA), list(tn, tn)))
stopifnot(max(abs(back(as.matrix(GA)) - back(as.matrix(GB)))) /
            max(1, max(abs(GA))) < 1e-12)
stopifnot(identical(as.numeric(attr(GA, "genetic_variance")[tn]),
                    as.numeric(attr(GB, "genetic_variance")[tn])))
stopifnot(identical(as.numeric(attr(GA, "residual_variance")[tn]),
                    as.numeric(attr(GB, "residual_variance")[tn])))
ok("ng_estimate_genetic_covariance(two_stage_ridge): same verdict through the exported entry point")

# ============================================================================
# 2. ng_posterior_genetic_covariance(): both methods
# ============================================================================
run_post <- function(Ym, method) suppressWarnings(ng_posterior_genetic_covariance(
  geno = geno, Y = Ym, n_draws = 6L, method = method, kfold = 5L, seed = 7L,
  allow_heuristic = TRUE))

for (meth in c("beta_posterior", "parametric_bootstrap")) {
  DA <- run_post(Y, meth)
  DB <- run_post(Y_perm, meth)
  # Reorder every draw of the permuted run back into the reference trait order.
  reord <- DB[tn, tn, , drop = FALSE]
  gap <- max(abs(DA - reord))
  rel <- gap / max(1, max(abs(DA)))
  cat(sprintf("  %s: max |draw difference| = %.3e absolute, %.3e relative\n", meth, gap, rel))
  stopifnot(rel < 1e-12)
  mean_gap <- max(abs(back(attr(DA, "G_mean")) - back(attr(DB, "G_mean"))))
  stopifnot(mean_gap / max(1, max(abs(attr(DA, "G_mean")))) < 1e-12)
  ok(sprintf("ng_posterior_genetic_covariance(%s): every draw is order-invariant", meth))
}

# ============================================================================
# 3. The bug is real: the distinct-stream property the fix must NOT break
# ============================================================================
# Traits must still get DIFFERENT streams. If the fix had simply given every
# trait one shared seed, the per-trait posterior draws would become identical
# innovations and manufacture cross-trait correlation. Assert the streams differ.
DA <- run_post(Y, "beta_posterior")
off <- apply(DA, 3L, function(M) M[1L, 2L] / sqrt(M[1L, 1L] * M[2L, 2L]))
stopifnot(stats::sd(off) > 0)
stopifnot(max(abs(off)) < 1 + 1e-8)
cat(sprintf("  beta_posterior implied r(yield, protein) across draws: sd = %.4f, range = [%.3f, %.3f]\n",
            stats::sd(off), min(off), max(off)))
ok("per-trait posterior streams remain DISTINCT (draw-to-draw correlation still varies)")

# ============================================================================
# 4. Identity, not position: renaming a trait DOES move its posterior stream
# ============================================================================
# The complement of the invariance above. The streams are keyed on the trait
# NAME, so a rename is expected to change the draws -- that is what makes the
# key an identity and not a constant.
Y_renamed <- Y
colnames(Y_renamed)[[1L]] <- "yield_kg"
DR <- run_post(Y_renamed, "beta_posterior")
stopifnot(!identical(as.numeric(DR[1L, 1L, ]), as.numeric(DA[1L, 1L, ])))
ok("renaming a trait moves its stream -- the key is the trait's identity, not a constant")

cat(sprintf("genetic_covariance_trait_order_invariance.R: PASS (%d checks)\n", checks))
