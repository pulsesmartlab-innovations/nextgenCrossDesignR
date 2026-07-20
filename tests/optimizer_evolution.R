ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(9)
n <- 40; m <- 300
ids <- paste0("P", seq_len(n))
geno <- matrix(2L * rbinom(n * m, 1, 0.45), n, m, dimnames = list(ids, paste0("M", seq_len(m))))
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n)); names(y) <- ids
mm <- data.frame(marker = paste0("M", seq_len(m)), chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 120, length.out = ceiling(m / 5)), 5)[seq_len(m)])
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 0L, seed = 9)
sc <- ng_score_crosses(geno, fit, marker_map = mm, ids = ids, adjusted_pheno = y,
                       selection_prop = 0.10, use_cpp = FALSE)
K <- ng_parent_kinship(geno)

obj <- function(p, lg, lpu) {
  s <- attr(p, "summary")
  s$mean_gain * nrow(p) - lg * s$group_coancestry - lpu * s$parent_use_sq
}

grid <- list(c(lg = 0, lpu = 0), c(lg = 1, lpu = 0), c(lg = 1, lpu = 2), c(lg = 2, lpu = 1))
for (cfg in grid) {
  lg <- cfg[["lg"]]; lpu <- cfg[["lpu"]]
  common <- list(scores = sc, n_crosses = 12, parent_kinship = K, max_crosses_per_parent = 5,
                 lambda_group = lg, lambda_parent_use = lpu, lambda_parent_use_mode = "absolute")
  g  <- do.call(ng_optimize_mating_plan, c(common, method = "greedy_local"))
  r  <- do.call(ng_optimize_mating_plan, c(common, method = "repair_local"))
  e1 <- do.call(ng_optimize_mating_plan, c(common, method = "evolution", evol_seed = 123L,
                                           evol_iterations = 80L, evol_solutions = 60L))
  e2 <- do.call(ng_optimize_mating_plan, c(common, method = "evolution", evol_seed = 123L,
                                           evol_iterations = 80L, evol_solutions = 60L))
  # feasible
  stopifnot(nrow(e1) == 12, attr(e1, "summary")$max_parent_use <= 5)
  # reproducible under a fixed seed
  stopifnot(isTRUE(all.equal(obj(e1, lg, lpu), obj(e2, lg, lpu))))
  # warm-start + elitism => never worse than greedy or repair (small tolerance)
  stopifnot(obj(e1, lg, lpu) >= obj(g, lg, lpu) - 1e-6)
  stopifnot(obj(e1, lg, lpu) >= obj(r, lg, lpu) - 1e-6)
  # within a small relative gap of the exact-ish MIP where available
  if (requireNamespace("lpSolve", quietly = TRUE)) {
    mp <- do.call(ng_optimize_mating_plan, c(common, method = "mip_contribution"))
    if (!isTRUE(attr(mp, "summary")$mip_fallback)) {
      om <- obj(mp, lg, lpu); oe <- obj(e1, lg, lpu)
      stopifnot(oe >= om - 0.05 * abs(om) - 1e-6)
    }
  }
}
cat("optimizer_evolution: reproducible, feasible, >= greedy/repair, competitive with MIP\n")
