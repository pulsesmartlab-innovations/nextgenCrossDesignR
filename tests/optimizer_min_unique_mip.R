# min_unique_parents is a HARD lower bound on the number of distinct parents used. It must be
# honored by EVERY optimizer path. Regression: the plain linear MIP (ng_mip_linear) has no such
# constraint, so `method = "mip_linear"` (and `method = "auto"` with lambda_group = 0, which
# routes there) silently produced plans with fewer unique parents than requested.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(1234)
n <- 24L; m <- 120L; ids <- sprintf("P%02d", seq_len(n)); mk <- sprintf("S%03d", seq_len(m))
G <- matrix(2L * rbinom(n * m, 1, 0.5), n, m, dimnames = list(ids, mk)); storage.mode(G) <- "double"
y <- setNames(as.numeric(G %*% rnorm(m, 0, 0.1)) + rnorm(n), ids)
map <- data.frame(marker = mk, chr = rep(1:6, each = m / 6),
                  pos_cm = rep(seq(0, 50, length.out = m / 6), 6))
eff <- ng_fit_ridge_effects(G, y, ids)
pairs <- as.data.frame(t(utils::combn(ids, 2)), stringsAsFactors = FALSE)
names(pairs) <- c("parent1", "parent2")
sc <- ng_score_crosses(G, eff, map, ids, pairs, target = "DH", assume_inbred = TRUE)
sc$gain <- sc$usefulness_pmv_gebv
pk <- ng_parent_kinship(G)

min_unique <- 18L
uniq <- function(plan) length(unique(c(as.character(plan$parent1), as.character(plan$parent2))))

# Every method must deliver at least `min_unique` distinct parents (the plan is feasible:
# 20 crosses x 2 parents = 40 slots, max 8/parent, so >= 18 distinct is achievable).
for (meth in c("auto", "mip_linear", "mip_contribution", "greedy_local", "evolution")) {
  plan <- ng_optimize_mating_plan(
    sc, n_crosses = 20L, gain_col = "gain", parent_kinship = pk,
    max_crosses_per_parent = 8L, min_unique_parents = min_unique,
    lambda_group = 0, lambda_mating = 0, method = meth
  )
  u <- uniq(plan)
  cat(sprintf("  %-16s unique parents = %d (>= %d: %s)\n", meth, u, min_unique, u >= min_unique))
  stopifnot(u >= min_unique)
}

cat("optimizer_min_unique_mip.R: PASS\n")
