helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# T1.6 + T2.7: the multi-trait index solve must
#  1. use a single ridge (not double-applied),
#  2. accept user-supplied genetic_covariance (Pesek-Baker) and
#     phenotypic_covariance (Smith-Hazel) and route them to the correct
#     formula via cov_source diagnostics,
#  3. autoscale the threshold penalty so a 1-SD violation costs a
#     predictable fraction of the index IQR.

set.seed(2026)

# ---- Synthetic candidate-cross trait table (yield, disease, lodging) ---------------------
n <- 60L
ids1 <- paste0("P", seq_len(20L))
pairs <- t(combn(ids1, 2L))[seq_len(n), , drop = FALSE]
yield   <- rnorm(n, mean = 10, sd = 2.0)
disease <- rnorm(n, mean = 30, sd = 5.0) - 0.3 * yield     # mild negative correlation
lodging <- rnorm(n, mean = 20, sd = 4.0) + 0.4 * yield     # mild positive correlation
scores <- data.frame(
  parent1 = pairs[, 1L], parent2 = pairs[, 2L],
  yield = yield, disease = disease, lodging = lodging,
  stringsAsFactors = FALSE
)
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "lodging"),
  direction = c("maximize", "minimize", "minimize"),
  economic_weight = c(2, 1, 1)
)

# ---- Missing quantitative-genetic matrices are rejected ----------------------------------
missing_both <- tryCatch(
  ng_add_multitrait_score(scores, traits, method = "economic_index"),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("requires both phenotypic_covariance", missing_both, fixed = TRUE))

# ---- G alone is not a Smith-Hazel index ---------------------------------------------------
G <- diag(c(4, 25, 16))                   # arbitrary diagonal G
dimnames(G) <- list(traits$trait, traits$trait)
missing_P <- tryCatch(
  ng_add_multitrait_score(scores, traits, method = "economic_index", genetic_covariance = G),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("requires both phenotypic_covariance", missing_P, fixed = TRUE))

# ---- User-supplied P AND G (Smith-Hazel) -------------------------------------------------
P <- diag(c(8, 60, 35))                   # arbitrary diagonal P
dimnames(P) <- list(traits$trait, traits$trait)
scored_sh <- ng_add_multitrait_score(
  scores, traits, method = "economic_index",
  phenotypic_covariance = P, genetic_covariance = G
)
meta_sh <- attr(scored_sh, "multi_trait")
stopifnot(identical(meta_sh$economic_index_cov_source, "smith_hazel"))
stopifnot(grepl("P\\^\\{-1\\} G a", meta_sh$economic_index_cov_solve_form))

# ---- Smith-Hazel predicted response uses G, not P ----------------------------------------
# With diagonal P and G, b = P^{-1} G a / |.|, predicted = G b (genetic-units response).
# The index is applied to value_z = sign*(raw - center)/scale, so raw-unit external
# covariances enter as M_z = L M L with L = diag(sign/scale). The predicted response is
# therefore the value_z-mapped G_z %*% b divided by the index SD, not the raw-unit
# G %*% b. L is rebuilt here from
# first principles (scale = IQR/1.349) rather than via ng_multitrait_cov_to_value_z, so this
# stays an independent check of the mapping instead of restating it.
sh_scale <- vapply(traits$trait, function(tt) stats::IQR(scores[[tt]]) / 1.349, numeric(1))
sh_sign <- ifelse(traits$direction == "maximize", 1, -1)
L_sh <- sh_sign / sh_scale
G_z <- G * outer(L_sh, L_sh)
P_z <- P * outer(L_sh, L_sh)

b_sh <- meta_sh$economic_index_coefficients
pred_sh <- meta_sh$economic_index_predicted_response
sigma_i <- sqrt(as.numeric(crossprod(b_sh, P_z %*% b_sh)))
expected_response <- as.numeric(G_z %*% b_sh) / sigma_i
if (max(abs(pred_sh - expected_response)) > 1e-8) {
  stop("Smith-Hazel predicted response should equal G_z %*% b / sigma_I")
}
# The defect this guards against: reporting the P-projection instead of the G-projection.
if (max(abs(pred_sh - as.numeric(P_z %*% b_sh))) < 1e-8) {
  stop("Smith-Hazel predicted response must not equal P %*% b")
}

# ---- Single ridge: coefficients equal the one-penalty Smith-Hazel solve ------------------
fit_small <- ng_multitrait_economic_index_fit(
  value_z = matrix(rnorm(n * 3, sd = 0.5), nrow = n),
  traits = traits, ridge = 1e-6,
  phenotypic_covariance = P, genetic_covariance = G
)
fit_huge <- ng_multitrait_economic_index_fit(
  value_z = matrix(rnorm(n * 3, sd = 0.5), nrow = n),
  traits = traits, ridge = 1e6,
  phenotypic_covariance = P, genetic_covariance = G
)
stopifnot(all(is.finite(fit_small$coefficients)))
stopifnot(all(is.finite(fit_huge$coefficients)))
one_ridge_reference <- function(ridge, target) {
  rhs <- as.numeric(G %*% target)
  b <- as.numeric(solve(P + diag(ridge * mean(diag(P)), nrow(P)), rhs))
  b / sum(abs(b))
}
ref_small <- one_ridge_reference(1e-6, fit_small$target)
ref_huge <- one_ridge_reference(1e6, fit_huge$target)
stopifnot(max(abs(fit_small$coefficients - ref_small)) < 1e-10)
stopifnot(max(abs(fit_huge$coefficients - ref_huge)) < 1e-10)

# ---- Threshold autoscale: a violation of one trait-SD should bind --------------------
traits_thr <- ng_multitrait_spec(
  trait = c("yield", "disease"),
  direction = c("maximize", "minimize"),
  weight = c(1, 1),
  min_value = c(8, NA),
  max_value = c(NA, 25)
)
scores_thr <- scores[, c("parent1", "parent2", "yield", "disease")]
scored_thr_no_auto <- ng_add_multitrait_score(
  scores_thr, traits_thr, method = "weighted",
  threshold_penalty_weight = 1.0, threshold_penalty_autoscale = FALSE
)
scored_thr_auto <- ng_add_multitrait_score(
  scores_thr, traits_thr, method = "weighted",
  threshold_penalty_weight = 1.0, threshold_penalty_autoscale = TRUE
)
meta_auto <- attr(scored_thr_auto, "multi_trait")
stopifnot(isTRUE(meta_auto$threshold_penalty_autoscale))
stopifnot(is.finite(meta_auto$effective_threshold_penalty))
# With autoscale, the effective penalty should be on the same order of
# magnitude as the IQR of the unpenalized weighted score (within a factor of 5).
score_no_penalty <- scored_thr_no_auto$multi_trait_score +
  scored_thr_no_auto$multi_trait_threshold_violation  # add back unpenalized
score_iqr <- stats::IQR(score_no_penalty[is.finite(score_no_penalty)], na.rm = TRUE)
ratio <- meta_auto$effective_threshold_penalty / max(score_iqr, 1e-8)
if (ratio < 0.5 || ratio > 5) {
  stop(sprintf("autoscaled penalty (%.3f) is not on the same order as score IQR (%.3f)",
               meta_auto$effective_threshold_penalty, score_iqr))
}

cat("multitrait_index_math: 5/5 checks passed\n")
cat(sprintf("  covariance contract covered: missing P/G rejected; Smith-Hazel P + G solved\n"))
cat("  ridge solve check: package coefficients equal the single-penalty closed form\n")
cat(sprintf("  autoscaled penalty/IQR ratio: %.3f\n", ratio))
