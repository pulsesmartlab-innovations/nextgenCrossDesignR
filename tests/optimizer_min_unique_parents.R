ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11)
n <- 20; m <- 200
ids <- paste0("P", seq_len(n))
geno <- matrix(2L * rbinom(n * m, 1, 0.45), n, m, dimnames = list(ids, paste0("M", seq_len(m))))
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n)); names(y) <- ids
mm <- data.frame(marker = paste0("M", seq_len(m)), chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 120, length.out = ceiling(m / 5)), 5)[seq_len(m)])
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 0L, seed = 11)
sc <- ng_score_crosses(geno, fit, marker_map = mm, ids = ids, adjusted_pheno = y,
                       selection_prop = 0.10, use_cpp = FALSE)
K <- ng_parent_kinship(geno)

# With a generous per-parent cap the naive top-gain plan concentrates on few parents.
# greedy_local and repair_local must now actively ENFORCE min_unique_parents (they
# previously only warned), reaching the target while staying feasible.
min_unique <- 15L
for (meth in c("greedy_local", "repair_local")) {
  p_free <- ng_optimize_mating_plan(sc, n_crosses = 10, parent_K = K,
                                    max_crosses_per_parent = 10, lambda_group = 0, method = meth)
  p_con <- ng_optimize_mating_plan(sc, n_crosses = 10, parent_K = K,
                                   max_crosses_per_parent = 10, lambda_group = 0,
                                   min_unique_parents = min_unique, method = meth)
  s_free <- attr(p_free, "summary"); s_con <- attr(p_con, "summary")
  stopifnot(nrow(p_con) == 10)
  stopifnot(s_con$max_parent_use <= 10)
  stopifnot(s_con$unique_parents >= min_unique)          # constraint actually enforced
  stopifnot(s_con$unique_parents > s_free$unique_parents) # and it did something
}
cat("optimizer_min_unique_parents: greedy_local + repair_local enforce the constraint\n")
