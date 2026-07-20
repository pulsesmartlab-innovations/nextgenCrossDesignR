helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(5)
n <- 24; m <- 200
ids <- paste0("P", seq_len(n))
geno <- matrix(2L * rbinom(n * m, 1, 0.45), n, m, dimnames = list(ids, paste0("M", seq_len(m))))
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n)); names(y) <- ids
mm <- data.frame(marker = paste0("M", seq_len(m)), chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 120, length.out = ceiling(m / 5)), 5)[seq_len(m)])
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 0L, seed = 5)
sc <- ng_score_crosses(geno, fit, marker_map = mm, ids = ids, adjusted_pheno = y,
                       selection_prop = 0.10, use_cpp = FALSE)
K <- ng_parent_kinship(geno)

if (requireNamespace("lpSolve", quietly = TRUE)) {
  # Normal MIP: solves, no fallback.
  p1 <- ng_optimize_mating_plan(sc, n_crosses = 10, parent_kinship = K, max_crosses_per_parent = 4,
                                lambda_group = 1, lambda_parent_use = 1, method = "mip_contribution")
  s1 <- attr(p1, "summary")
  stopifnot(nrow(p1) == 10, isFALSE(s1$mip_fallback))

  # Oversized size guard -> never hangs, falls back to greedy, records mip_fallback,
  # still returns a valid n_crosses feasible plan.
  p2 <- suppressWarnings(ng_optimize_mating_plan(
    sc, n_crosses = 10, parent_kinship = K, max_crosses_per_parent = 4,
    lambda_group = 1, lambda_parent_use = 1, method = "mip_contribution",
    mip_max_binary_vars = 50, mip_time_limit = 2))
  s2 <- attr(p2, "summary")
  stopifnot(nrow(p2) == 10, isTRUE(s2$mip_fallback), s2$max_parent_use <= 4)

  # method="auto" also degrades gracefully under the size guard.
  p3 <- suppressWarnings(ng_optimize_mating_plan(
    sc, n_crosses = 10, parent_kinship = K, max_crosses_per_parent = 4,
    lambda_group = 1, lambda_parent_use = 1, method = "auto", mip_max_binary_vars = 50))
  stopifnot(nrow(p3) == 10)
  cat("optimizer_mip_guard: normal solve + oversized fallback + auto degrade OK\n")
} else {
  cat("lpSolve unavailable; skipping MIP guard test\n")
}
