helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# D3: built-in genetic-covariance estimator G_hat and the companion phenotypic
# covariance estimator P_hat.
#
# ============================================================================
# 0.29.0 -- WHAT THIS FILE ASSERTS, AND WHAT IT ONLY CHARACTERISES
# ============================================================================
#
# Through 0.28.0 this file asserted that two_stage_ridge recovers G_true to a
# relative Frobenius error <= 0.5 on the COVARIANCE matrix, at one hard-coded
# simulation seed. That criterion was measured, over 200 independently simulated
# datasets from this exact generating model, to be met by 12% of them
# (mean relative error 1.88, sd 2.06, median 1.26, max 19.3). Seed 2026 passed by
# chance. A 12%-pass criterion certifies nothing: it fails on any reseeding, and
# it would keep passing on an estimator that had become much worse.
#
# Diagnosing it separates the estimator into two halves that behave completely
# differently, exactly as R/31's own comments predict:
#
#   OFF-DIAGONALS are the Pearson correlation of the per-trait ridge beta
#   vectors. Pearson correlation is invariant to per-trait multiplicative
#   shrinkage, so it survives the ridge attenuation. Measured over the same 200
#   datasets, the relative Frobenius error of the genetic CORRELATION matrix is
#   mean 0.226, sd 0.067, median 0.221, max 0.479 -- 200/200 below 0.6.
#
#   DIAGONALS are the GBLUP lambda inversion sigma_g2 = sigma_e2 * denom / lambda.
#   That is a heuristic, and it is where all the error lives: the worst per-trait
#   genetic variance is out by a factor of exp(1.82) ~ 6.2 on average and by up to
#   exp(6.94) ~ 1000. It is also biased UPWARD often enough that G_hat exceeds the
#   companion P_hat on 87% of datasets, i.e. it implies h2 > 1.
#
# This file therefore now:
#   * ASSERTS the correlation-structure property, which genuinely holds, with a
#     bound justified by that 200-dataset measurement, over SEVERAL independent
#     datasets rather than one;
#   * CHARACTERISES the covariance-scale error -- printed, and bounded only
#     loosely enough to document the known bias, never tightened to a lucky seed;
#   * ASSERTS that the pair (G_hat, P_hat) is REFUSED by the 0.29.0 P - G guard on
#     this simulation, because it implies h2 > 1. That refusal is the guard doing
#     its job on the package's own estimator, and it is what the Smith-Hazel
#     integration below now has to route around.
#
# The estimator's structural guarantees (symmetry, PSD, positive diagonals,
# provenance attributes) are asserted unconditionally, as before.
# ============================================================================

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }

t_traits <- 3L
trait_names <- c("trait_a", "trait_b", "trait_c")
sd_g <- c(1, sqrt(2), sqrt(0.5))
R_true <- matrix(c(1.0,  0.3, -0.4,
                   0.3,  1.0,  0.1,
                  -0.4,  0.1,  1.0), nrow = 3L, byrow = TRUE)
G_true <- diag(sd_g) %*% R_true %*% diag(sd_g)
dimnames(G_true) <- dimnames(R_true) <- list(trait_names, trait_names)
# The generating model's own phenotypic covariance: residual variance equals the
# genetic variance per trait, so every true h2 is exactly 0.5 and P_true - G_true
# is PD by construction.
E_var <- diag(G_true)
P_true <- G_true + diag(E_var)
dimnames(P_true) <- list(trait_names, trait_names)

# One draw from the generating model used throughout the file.
simulate_panel <- function(seed) {
  set.seed(seed)
  n <- 100L
  m <- 200L
  maf <- stats::runif(m, 0.1, 0.45)
  geno <- matrix(0, nrow = n, ncol = m)
  for (k in seq_len(m)) geno[, k] <- stats::rbinom(n, size = 2L, prob = maf[k])
  rownames(geno) <- paste0("L", seq_len(n))
  colnames(geno) <- paste0("M", seq_len(m))
  p_hat <- colMeans(geno) / 2
  denom <- sum(2 * p_hat * (1 - p_hat))
  # Effects independent across markers, correlated across traits via G_true, and
  # scaled so sum_k w_k Cov(beta_k) ~ G_true.
  eig <- eigen(G_true, symmetric = TRUE)
  sqrt_G <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  beta_true <- matrix(stats::rnorm(m * t_traits), nrow = m) %*% (sqrt_G / sqrt(denom))
  Xc <- sweep(geno, 2L, 2 * p_hat, "-")   # the package's centering
  E <- matrix(stats::rnorm(n * t_traits), nrow = n) %*% diag(sqrt(E_var))
  Y <- Xc %*% beta_true + E
  colnames(Y) <- trait_names
  rownames(Y) <- rownames(geno)
  list(geno = geno, Y = Y, n = n, m = m)
}

fro_ratio <- function(A, B) sqrt(sum((A - B)^2)) / sqrt(sum(B^2))

# ---- Reference panel (the historical seed) ---------------------------------
panel <- simulate_panel(2026L)
geno <- panel$geno
Y <- panel$Y
n <- panel$n

G_hat <- ng_estimate_genetic_covariance(
  geno = geno, Y = Y,
  method = "two_stage_ridge",
  kfold = 5L, seed = 2026L,
  return_diagnostics = TRUE
)

# ============================================================================
# ASSERTED: structural guarantees -- unconditional, on every input
# ============================================================================
stopifnot(identical(dim(G_hat), c(t_traits, t_traits)))
stopifnot(max(abs(G_hat - t(G_hat))) < 1e-8)
stopifnot(all(diag(G_hat) > 0))
ev_hat <- eigen(G_hat, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev_hat) > -1e-8)
R_hat <- attr(G_hat, "genetic_correlation")
stopifnot(!is.null(R_hat))
stopifnot(max(abs(diag(R_hat) - 1)) < 1e-8)
stopifnot(max(abs(R_hat[upper.tri(R_hat)])) <= 1 + 1e-8)
stopifnot(identical(attr(G_hat, "n_used"), n))
stopifnot(identical(attr(G_hat, "method"), "two_stage_ridge"))
ok("G_hat is symmetric PSD with positive diagonals, unit correlation diagonal and correct provenance")

# ============================================================================
# ASSERTED: the correlation-structure property, over SEVERAL datasets
# ============================================================================
# The bound is 0.6 on the relative Frobenius error of the genetic CORRELATION
# matrix. Justification, from 200 independently simulated datasets of this exact
# model: mean 0.226, sd 0.067, max 0.479 -- so 0.6 sits about 5.6 standard
# deviations above the mean and 25% above the worst case observed. It is a bound
# on a property that holds, not a threshold fitted to a seed. The assertion runs
# over 8 independent datasets so a single lucky draw cannot carry it.
cor_seeds <- as.integer(seq(3001L, 3008L))
cor_ratios <- vapply(cor_seeds, function(s) {
  p <- simulate_panel(s)
  Gh <- ng_estimate_genetic_covariance(p$geno, p$Y, method = "two_stage_ridge",
                                       kfold = 5L, seed = s)
  fro_ratio(ng_genetic_cov_to_correlation(matrix(as.numeric(Gh), t_traits, t_traits)),
            R_true)
}, numeric(1L))
cat(sprintf("  genetic-correlation recovery over %d independent datasets: %s\n",
            length(cor_seeds), paste(sprintf("%.3f", cor_ratios), collapse = ", ")))
cat(sprintf("    mean = %.3f, max = %.3f (asserted bound 0.6; 200-dataset reference max 0.479)\n",
            mean(cor_ratios), max(cor_ratios)))
if (max(cor_ratios) > 0.6) {
  stop(sprintf(paste0("two_stage_ridge no longer recovers the genetic CORRELATION structure: ",
                      "worst relative Frobenius error %.3f over %d datasets exceeds 0.6 ",
                      "(reference distribution over 200 datasets: mean 0.226, sd 0.067, max 0.479)"),
               max(cor_ratios), length(cor_seeds)))
}
ok("two_stage_ridge recovers the genetic CORRELATION matrix to <= 0.6 relative Frobenius error on every dataset")

# ============================================================================
# CHARACTERISING (not a guard): the covariance-SCALE error
# ============================================================================
# Printed so a reader can see what the estimator does, and bounded only loosely
# enough to catch a total collapse. It is deliberately NOT tightened to a seed.
# The pre-0.29.0 criterion for reference: cov_ratio <= 0.5, met by 12% of
# datasets. This is characterisation, not certification.
cov_ratio <- fro_ratio(matrix(as.numeric(G_hat), t_traits, t_traits), G_true)
diag_logratio <- max(abs(log(diag(G_hat) / diag(G_true))))
cat(sprintf("  CHARACTERISATION at seed 2026: covariance relative Frobenius error = %.3f\n", cov_ratio))
cat(sprintf("    (200-dataset reference: mean 1.88, sd 2.06, median 1.26, max 19.3; P(<= 0.5) = 0.12)\n"))
cat(sprintf("    worst per-trait genetic variance off by a factor of %.2f (GBLUP lambda inversion)\n",
            exp(diag_logratio)))
cat(sprintf("    diag(G_hat) = [%s] vs diag(G_true) = [%s]\n",
            paste(sprintf("%.3f", diag(G_hat)), collapse = ", "),
            paste(sprintf("%.3f", diag(G_true)), collapse = ", ")))
# A collapse ceiling only: three orders of magnitude on the variance scale.
stopifnot(cov_ratio < 1e3, exp(diag_logratio) < 1e3)
ok("covariance-scale error is CHARACTERISED and printed, not certified (the 0.5 criterion is withdrawn)")

# ============================================================================
# Estimator: phenotypic covariance with Ledoit-Wolf shrinkage -- unchanged
# ============================================================================
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
ok("P_hat is symmetric PSD, shrunk off-diagonals, diagonals within 10% of the sample variances")

# ============================================================================
# ASSERTED (0.29.0): the P - G guard REFUSES the package's own (G_hat, P_hat)
# ============================================================================
# This is the same finding as the withdrawn accuracy criterion, seen from the
# other side. The lambda-inversion diagonal overshoots, so diag(G_hat) exceeds
# diag(P_hat) and the pair implies h2 > 1 -- on 87% of datasets from this model.
# Before 0.29.0 that pair was accepted and solved, and the resulting Smith-Hazel
# weights b = P^{-1} G a were computed from a P = G + R decomposition that does
# not exist. The guard now refuses it, by name and by number.
h2_implied <- diag(G_hat) / diag(P_auto)
cat(sprintf("  implied per-trait h2 from (G_hat, P_hat): [%s]\n",
            paste(sprintf("%.3f", h2_implied), collapse = ", ")))
stopifnot(any(h2_implied > 1))
msg <- tryCatch({
  ng_multitrait_validate_cov_pair(
    ng_multitrait_validate_cov(P_auto, trait_names, "phenotypic_covariance"),
    ng_multitrait_validate_cov(matrix(as.numeric(G_hat), t_traits, t_traits,
                                      dimnames = list(trait_names, trait_names)),
                               trait_names, "genetic_covariance"))
  NA_character_
}, error = function(e) conditionMessage(e))
stopifnot(!is.na(msg))
stopifnot(grepl("heritability above 1", msg, fixed = TRUE))
cat("  guard message:", substr(msg, 1L, 200L), "...\n")
ok("the 0.29.0 P - G guard refuses the package's own (G_hat, P_hat) pair on this simulation")

# ============================================================================
# Integration: Smith-Hazel route through ng_add_multitrait_score
# ============================================================================
# The integration is exercised with a VALID pair -- the generating model's own
# G_true and P_true, whose h2 is 0.5 for every trait -- because a pair that fails
# the guard is, correctly, no longer accepted. What is being tested here is that
# a valid P/G pair reaches the Smith-Hazel solve, not the estimator's accuracy.
set.seed(2026L)
n_cross <- 30L
parents <- paste0("Par", seq_len(20L))
pair_idx <- t(combn(seq_along(parents), 2L))[seq_len(n_cross), , drop = FALSE]
scores <- data.frame(
  parent1 = parents[pair_idx[, 1L]],
  parent2 = parents[pair_idx[, 2L]],
  trait_a = stats::rnorm(n_cross, mean = 5, sd = 1),
  trait_b = stats::rnorm(n_cross, mean = 10, sd = sqrt(2)),
  trait_c = stats::rnorm(n_cross, mean = 2, sd = sqrt(0.5)),
  stringsAsFactors = FALSE
)
traits <- ng_multitrait_spec(
  trait = trait_names,
  direction = c("maximize", "maximize", "maximize"),
  economic_weight = c(2, 1, 1)
)
scored <- ng_add_multitrait_score(
  scores = scores, traits = traits, method = "economic_index",
  phenotypic_covariance = P_true,
  genetic_covariance = G_true
)
meta <- attr(scored, "multi_trait")
stopifnot(identical(meta$economic_index_cov_source, "smith_hazel"))
stopifnot(grepl("P\\^\\{-1\\} G a", meta$economic_index_cov_solve_form))
stopifnot(all(is.finite(meta$economic_index_coefficients)))
ok("a guard-valid (G, P) pair reaches the Smith-Hazel solve through ng_add_multitrait_score")

cat(sprintf("genetic_covariance_estimator: PASS (%d checks)\n", checks))
cat(sprintf("  ASSERTED: genetic-correlation recovery, worst of %d datasets = %.3f (bound 0.6)\n",
            length(cor_seeds), max(cor_ratios)))
cat(sprintf("  CHARACTERISED: covariance-scale relative Frobenius error = %.3f (no bound asserted)\n",
            cov_ratio))
cat(sprintf("  Ledoit-Wolf shrinkage intensity: %.3f\n", intensity))
