helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# D2: full off-diagonal posterior PMV. With Sigma_beta = diag(beta_var) the new
# pmv_full_posterior must exactly reproduce the legacy diagonal column.
# With a non-trivial off-diagonal Sigma_beta the full posterior PMV must add
#   d' (R o Sigma_beta) d   to the VPM, matching the genomicMateSelectR
# formulation.

# ---- Test 1: diagonal Sigma_beta reproduces the legacy diagonal column -------------
set.seed(2026L)
n <- 8L
m <- 12L
geno <- matrix(2L * rbinom(n * m, 1, 0.5), nrow = n)
ids <- paste0("P", seq_len(n))
markers <- paste0("M", seq_len(m))
rownames(geno) <- ids; colnames(geno) <- markers
beta <- setNames(rnorm(m, sd = 0.1), markers)
beta_var <- setNames(runif(m, min = 0.001, max = 0.01), markers)
Sigma_diag <- diag(beta_var)
rownames(Sigma_diag) <- markers; colnames(Sigma_diag) <- markers

mm <- ng_prepare_marker_map(
  data.frame(marker = markers,
             chr = rep(1:2, length.out = m),
             pos_cm = rep(seq(0, 100, length.out = m / 2), 2)),
  marker_ids = markers, model = "haldane"
)
sorted <- ng_sort_by_map(geno, beta, beta_var, mm)
pairs <- ng_make_pairs(ids)

diag_dh <- ng_dh_recomb_variance_pairs_dense(
  geno = sorted$geno, beta = sorted$effects, beta_var = sorted$beta_var,
  marker_map = sorted$marker_map, ids = ids, pairs = pairs,
  recomb_model = "haldane", target = "DH"
)
full_diag <- ng_dh_recomb_variance_pairs_full_posterior(
  geno = sorted$geno, beta = sorted$effects, beta_cov_full = Sigma_diag,
  marker_map = sorted$marker_map, ids = ids, pairs = pairs,
  target = "DH", recomb_model = "haldane"
)
stopifnot("pmv_full_posterior" %in% names(full_diag))
err_vpm <- max(abs(full_diag$vpm - diag_dh$vpm))
err_pmv <- max(abs(full_diag$pmv - diag_dh$pmv))
err_full <- max(abs(full_diag$pmv_full_posterior - diag_dh$pmv))
if (err_vpm > 1e-9 || err_pmv > 1e-9 || err_full > 1e-9) {
  stop(sprintf("diagonal Sigma_beta: VPM err=%g, PMV err=%g, full err=%g",
               err_vpm, err_pmv, err_full))
}

# ---- Test 2: 4-parent x 8-marker fixture with off-diagonal Sigma_beta ----------------
# Construct Sigma_beta = sigma_e2 * (X'X + lambda I)^{-1} on a real ridge fit
# (return_beta_cov_full = TRUE). Verify (a) diag(Sigma_beta) matches the
# scalar beta_var column, (b) Sigma_beta is symmetric, (c) the full-posterior
# PMV is >= the diagonal PMV for crosses where the off-diagonal structure
# aligns with the d_k d_l (1 - 2r) terms.
set.seed(77L)
n2 <- 24L
m2 <- 8L
Z <- matrix(rnorm(n2 * 3L), nrow = n2)
loadings <- matrix(rnorm(3L * m2, sd = 0.6), nrow = 3L)
geno2 <- Z %*% loadings + matrix(rnorm(n2 * m2, sd = 0.2), nrow = n2)
ids2 <- paste0("L", seq_len(n2))
markers2 <- paste0("M", seq_len(m2))
rownames(geno2) <- ids2; colnames(geno2) <- markers2
y2 <- as.numeric(geno2 %*% rnorm(m2, sd = 0.4) + rnorm(n2, sd = 0.5))
names(y2) <- ids2

fit_full <- ng_fit_ridge_effects(geno2, y2, ids = ids2, lambda = 3,
                                 kfold = 4L, seed = 77L,
                                 return_beta_cov_full = TRUE)
stopifnot(!is.null(fit_full$beta_cov_full))
stopifnot(isTRUE(all.equal(dim(fit_full$beta_cov_full), c(m2, m2))))
stopifnot(max(abs(fit_full$beta_cov_full - t(fit_full$beta_cov_full))) < 1e-10)
# diag(Sigma_beta) must equal beta_var to ~1e-8
diag_err <- max(abs(diag(fit_full$beta_cov_full) - fit_full$beta_var))
if (diag_err > 1e-8) {
  stop(sprintf("diag(beta_cov_full) != beta_var, max abs err = %g", diag_err))
}
# Sigma_beta should match the primal (sigma_e2 * (X'X + lambda I)^{-1}).
Xc <- sweep(geno2, 2L, colMeans(geno2), "-")
primal <- fit_full$sigma_e2 * solve(crossprod(Xc) + diag(fit_full$lambda, m2))
prim_err <- max(abs(fit_full$beta_cov_full - primal))
if (prim_err > 1e-7) {
  stop(sprintf("beta_cov_full does not match primal inverse, max abs err = %g", prim_err))
}

# Use a small parent sample of 4 for the per-cross PMV check.
parents <- ids2[1:4]
geno_p <- geno2[parents, , drop = FALSE]
# Round to {0, 2} to satisfy the inbred-parent contract used by ng_score_crosses.
geno_p[] <- ifelse(geno_p > median(geno_p), 2L, 0L)
mm2 <- ng_prepare_marker_map(
  data.frame(marker = markers2, chr = c(1, 1, 1, 1, 2, 2, 2, 2),
             pos_cm = c(0, 10, 25, 50, 0, 15, 35, 70)),
  marker_ids = markers2, model = "haldane"
)
pairs_p <- ng_make_pairs(parents)
sorted_p <- ng_sort_by_map(geno_p, fit_full$beta, fit_full$beta_var, mm2)

full_pp <- ng_dh_recomb_variance_pairs_full_posterior(
  geno = sorted_p$geno, beta = sorted_p$effects,
  beta_cov_full = fit_full$beta_cov_full,
  marker_map = sorted_p$marker_map, ids = parents, pairs = pairs_p,
  target = "DH", recomb_model = "haldane"
)
# Build the diagonal-only PMV for the same fit by zeroing the off-diagonals.
Sigma_only_diag <- diag(diag(fit_full$beta_cov_full))
rownames(Sigma_only_diag) <- markers2; colnames(Sigma_only_diag) <- markers2
diag_pp <- ng_dh_recomb_variance_pairs_full_posterior(
  geno = sorted_p$geno, beta = sorted_p$effects,
  beta_cov_full = Sigma_only_diag,
  marker_map = sorted_p$marker_map, ids = parents, pairs = pairs_p,
  target = "DH", recomb_model = "haldane"
)
# VPM (a'Ra) and diagonal pmv must be invariant to the off-diagonal Sigma.
stopifnot(max(abs(full_pp$vpm - diag_pp$vpm)) < 1e-10)
stopifnot(max(abs(full_pp$pmv - diag_pp$pmv)) < 1e-10)
# At least one pair must move (off-diagonals are non-zero and d_k d_l are
# non-trivial); the column must exist and not collapse to NA.
stopifnot(all(is.finite(full_pp$pmv_full_posterior)))
delta_full <- full_pp$pmv_full_posterior - diag_pp$pmv_full_posterior
if (max(abs(delta_full)) < 1e-8) {
  stop("off-diagonal Sigma_beta produced no change in pmv_full_posterior")
}

# ---- Test 3: 2-locus closed-form sanity check --------------------------------------
# Two pure inbreds, 1 chromosome, markers at 0 and 30 cM, equal positive effects,
# non-zero off-diagonal posterior covariance, zero diagonal posterior variance
# so the analytic check isolates the off-diagonal contribution.
geno3 <- rbind(P1 = c(2, 2), P2 = c(0, 0))
colnames(geno3) <- c("M1", "M2")
beta3 <- setNames(c(0.4, 0.6), c("M1", "M2"))
d_cm <- 30
r3 <- ng_meiosis_r(d_cm, model = "haldane")
decay3 <- 1 - 2 * r3  # this is the value placed at R[1, 2]
cov12 <- 0.05         # off-diagonal posterior covariance
Sigma3 <- matrix(c(0, cov12, cov12, 0), nrow = 2)
rownames(Sigma3) <- c("M1", "M2"); colnames(Sigma3) <- c("M1", "M2")
mm3 <- ng_prepare_marker_map(
  data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, d_cm)),
  marker_ids = c("M1", "M2"), model = "haldane"
)
pairs3 <- data.frame(parent1 = "P1", parent2 = "P2", stringsAsFactors = FALSE)
full3 <- ng_dh_recomb_variance_pairs_full_posterior(
  geno = geno3, beta = beta3, beta_cov_full = Sigma3,
  marker_map = mm3, ids = rownames(geno3), pairs = pairs3,
  target = "DH", recomb_model = "haldane"
)
# d_k = 0.5 * (2 - 0) = 1 for both markers, so the closed form reduces to
#   PMV_full = beta_1^2 + beta_2^2 + 2 (beta_1 beta_2 + cov12) (1 - 2 r)
# because diag(Sigma) = 0 (no per-marker variance).
expect_vpm <- beta3[1]^2 + beta3[2]^2 + 2 * beta3[1] * beta3[2] * decay3
expect_pmv_full <- beta3[1]^2 + beta3[2]^2 +
  2 * (beta3[1] * beta3[2] + cov12) * decay3
err_vpm3 <- abs(full3$vpm - expect_vpm)
err_full3 <- abs(full3$pmv_full_posterior - expect_pmv_full)
if (err_vpm3 > 1e-10) {
  stop(sprintf("2-locus VPM disagrees with closed form: err=%g", err_vpm3))
}
if (err_full3 > 1e-10) {
  stop(sprintf("2-locus full-posterior PMV disagrees with closed form: err=%g", err_full3))
}
# diag(Sigma)=0 here, so the diagonal-only PMV must equal the VPM.
stopifnot(abs(full3$pmv - full3$vpm) < 1e-10)
# With cov12 > 0 and a' R a using a positive (1-2r), the full-posterior PMV
# must exceed the diagonal-only PMV.
stopifnot(full3$pmv_full_posterior > full3$pmv)

# ---- Test 4: ng_score_crosses wiring (posterior_cov_full = ...) --------------------
# Default (NULL) -> pmv_full_posterior column is NA. With the matrix
# supplied, the column is finite and matches the kernel.
effects_p <- list(beta = fit_full$beta, beta_var = fit_full$beta_var,
                  beta_cov_full = fit_full$beta_cov_full,
                  reliability = 0.5, intercept = 0,
                  marker_mean = colMeans(geno_p))
adj_p <- setNames(rnorm(length(parents)), parents)
sc_legacy <- ng_score_crosses(
  geno = geno_p, effects = effects_p,
  marker_map = mm2, ids = parents, adjusted_pheno = adj_p,
  selection_prop = 0.10, recomb_model = "haldane", use_cpp = FALSE
)
stopifnot("pmv_full_posterior" %in% names(sc_legacy))
stopifnot(all(is.na(sc_legacy$pmv_full_posterior)))

sc_full <- ng_score_crosses(
  geno = geno_p, effects = effects_p,
  marker_map = mm2, ids = parents, adjusted_pheno = adj_p,
  selection_prop = 0.10, recomb_model = "haldane", use_cpp = FALSE,
  posterior_cov_full = fit_full$beta_cov_full
)
stopifnot(all(is.finite(sc_full$pmv_full_posterior)))
# The wired path must agree with a direct kernel call (after aligning by pair).
key_l <- paste(sc_legacy$parent1, sc_legacy$parent2, sep = "x")
key_f <- paste(sc_full$parent1, sc_full$parent2, sep = "x")
stopifnot(identical(key_l, key_f))
# Diagonal columns are unchanged across the two calls.
stopifnot(max(abs(sc_legacy$vpm - sc_full$vpm)) < 1e-10)
stopifnot(max(abs(sc_legacy$pmv - sc_full$pmv)) < 1e-10)
# Full posterior PMV must be >= VPM (extra variance is non-negative when Sigma
# is positive semi-definite, which the ridge dual identity guarantees).
stopifnot(all(sc_full$pmv_full_posterior + 1e-10 >= sc_full$vpm))

cat("posterior_pmv_full: 4/4 checks passed\n")
cat(sprintf("  max |diag-only - legacy diagonal PMV| = %.3e\n", err_full))
cat(sprintf("  max |off-diag Sigma effect on pmv_full_posterior| = %.3e\n",
            max(abs(delta_full))))
cat(sprintf("  2-locus closed-form full PMV err = %.3e\n", err_full3))
