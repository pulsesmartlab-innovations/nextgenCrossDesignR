helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# D3: built-in genetic-covariance estimator G_hat and the companion phenotypic
# covariance estimator P_hat. Validates:
#   1. two_stage_ridge recovers G_true to a loose Frobenius tolerance,
#   2. G_hat is symmetric PSD with positive diagonals,
#   3. ng_estimate_phenotypic_covariance with shrinkage = "auto" returns a
#      symmetric PSD P_hat whose off-diagonals are shrunk relative to the
#      unshrunken sample covariance, with diagonals approximately equal to
#      the sample variances,
#   4. G_hat and P_hat plug into ng_add_multitrait_score(method =
#      "economic_index") and trigger the Smith-Hazel solve path.

set.seed(2026L)

n <- 100L
m <- 200L
t_traits <- 3L

# --- Simulate synthetic 3-trait scenario ------------------------------------
maf <- runif(m, 0.1, 0.45)
geno <- matrix(0, nrow = n, ncol = m)
for (k in seq_len(m)) geno[, k] <- rbinom(n, size = 2L, prob = maf[k])
rownames(geno) <- paste0("L", seq_len(n))
colnames(geno) <- paste0("M", seq_len(m))

# Target G_true: diag(c(1, 2, 0.5)) with off-diagonal correlations 0.3, -0.4, 0.1.
sd_g <- c(1, sqrt(2), sqrt(0.5))
R_true <- matrix(c(1.0,  0.3, -0.4,
                   0.3,  1.0,  0.1,
                  -0.4,  0.1,  1.0), nrow = 3L, byrow = TRUE)
G_true <- diag(sd_g) %*% R_true %*% diag(sd_g)
dimnames(G_true) <- list(c("trait_a", "trait_b", "trait_c"),
                         c("trait_a", "trait_b", "trait_c"))

# Per-marker scaling so that the realized additive variance from marker effects
# matches diag(G_true). Effects are independent across markers but correlated
# across traits via G_true.
p_hat <- colMeans(geno) / 2
w_k <- 2 * p_hat * (1 - p_hat)
denom <- sum(w_k)
# Draw effects: rows = markers, cols = traits, with cross-trait covariance
# G_true / denom (so that sum_k w_k * Cov(beta_k) ~ G_true).
eig <- eigen(G_true, symmetric = TRUE)
sqrt_G <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
beta_true <- matrix(rnorm(m * t_traits), nrow = m) %*% (sqrt_G / sqrt(denom))
# Center geno (per-marker mean 2*p_hat) so we match the package's centering.
Xc <- sweep(geno, 2L, 2 * p_hat, "-")
G_breeding <- Xc %*% beta_true
# Residual error with mild correlation; h2 ~ 0.5-0.6 per trait.
E_sd <- sqrt(diag(G_true))
E <- matrix(rnorm(n * t_traits), nrow = n) %*% diag(E_sd)
Y <- G_breeding + E
colnames(Y) <- c("trait_a", "trait_b", "trait_c")
rownames(Y) <- rownames(geno)

# --- Estimator: two-stage ridge ---------------------------------------------
G_hat <- ng_estimate_genetic_covariance(
  geno = geno, Y = Y,
  method = "two_stage_ridge",
  kfold = 5L, seed = 2026L,
  return_diagnostics = TRUE
)

# Symmetry, positive diagonals, PSD-ness.
stopifnot(identical(dim(G_hat), c(t_traits, t_traits)))
stopifnot(max(abs(G_hat - t(G_hat))) < 1e-8)
stopifnot(all(diag(G_hat) > 0))
ev_hat <- eigen(G_hat, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev_hat) > -1e-8)

# Frobenius-norm recovery.
froberr <- sqrt(sum((G_hat - G_true)^2))
fronorm <- sqrt(sum(G_true^2))
ratio <- froberr / fronorm
if (ratio > 0.5) {
  stop(sprintf("two_stage_ridge G_hat too far from G_true: Frobenius ratio = %.3f", ratio))
}

# Diagnostics attributes.
R_hat <- attr(G_hat, "genetic_correlation")
stopifnot(!is.null(R_hat))
stopifnot(max(abs(diag(R_hat) - 1)) < 1e-8)
stopifnot(identical(attr(G_hat, "n_used"), n))
stopifnot(identical(attr(G_hat, "method"), "two_stage_ridge"))

# --- Estimator: phenotypic covariance with Ledoit-Wolf shrinkage ------------
P_none <- ng_estimate_phenotypic_covariance(Y, shrinkage = "none")
P_auto <- ng_estimate_phenotypic_covariance(Y, shrinkage = "auto")
stopifnot(max(abs(P_auto - t(P_auto))) < 1e-8)
ev_p <- eigen(P_auto, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev_p) > -1e-8)
sample_var <- diag(stats::var(Y))
stopifnot(max(abs(diag(P_auto) - sample_var) / sample_var) < 0.10)
off_none <- P_none[upper.tri(P_none)]
off_auto <- P_auto[upper.tri(P_auto)]
if (sum(abs(off_auto)) >= sum(abs(off_none))) {
  stop("Ledoit-Wolf shrinkage should reduce off-diagonal magnitudes")
}
intensity <- attr(P_auto, "intensity")
stopifnot(intensity >= 0 && intensity <= 1)

# --- Integration: Smith-Hazel route through ng_add_multitrait_score ---------
n_cross <- 30L
parents <- paste0("Par", seq_len(20L))
pair_idx <- t(combn(seq_along(parents), 2L))[seq_len(n_cross), , drop = FALSE]
scores <- data.frame(
  parent1 = parents[pair_idx[, 1L]],
  parent2 = parents[pair_idx[, 2L]],
  trait_a = rnorm(n_cross, mean = 5, sd = 1),
  trait_b = rnorm(n_cross, mean = 10, sd = sqrt(2)),
  trait_c = rnorm(n_cross, mean = 2, sd = sqrt(0.5)),
  stringsAsFactors = FALSE
)
traits <- ng_multitrait_spec(
  trait = c("trait_a", "trait_b", "trait_c"),
  direction = c("maximize", "maximize", "maximize"),
  economic_weight = c(2, 1, 1)
)
scored <- ng_add_multitrait_score(
  scores = scores, traits = traits, method = "economic_index",
  phenotypic_covariance = P_auto,
  genetic_covariance = G_hat
)
meta <- attr(scored, "multi_trait")
stopifnot(identical(meta$economic_index_cov_source, "smith_hazel"))
stopifnot(grepl("P\\^\\{-1\\} G a", meta$economic_index_cov_solve_form))
stopifnot(all(is.finite(meta$economic_index_coefficients)))

cat("genetic_covariance_estimator: 4/4 checks passed\n")
cat(sprintf("  two_stage_ridge Frobenius error: %.4f (||G_true||_F = %.4f, ratio = %.3f)\n",
            froberr, fronorm, ratio))
cat(sprintf("  Ledoit-Wolf shrinkage intensity: %.3f\n", intensity))
cat(sprintf("  smith_hazel cov_source confirmed via ng_add_multitrait_score\n"))
