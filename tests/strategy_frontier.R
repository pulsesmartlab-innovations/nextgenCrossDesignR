# Gain-diversity balancing: strategy / diversity-emphasis frontier control.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(101)
n_parents <- 24L
parents <- sprintf("P%02d", seq_len(n_parents))

# Synthetic candidate crosses with a gain column and a pair_kinship column, plus a
# parent relationship matrix (VanRaden-style, PSD by construction).
pairs <- t(utils::combn(n_parents, 2L))
scores <- data.frame(
  parent1 = parents[pairs[, 1]],
  parent2 = parents[pairs[, 2]],
  stringsAsFactors = FALSE
)
scores$usefulness_pmv_gebv <- rnorm(nrow(scores), 10, 2)
L <- matrix(rnorm(n_parents * n_parents, 0, 0.3), n_parents, n_parents)
K <- crossprod(L) / n_parents
diag(K) <- diag(K) + 1
K <- K / mean(diag(K))
dimnames(K) <- list(parents, parents)
scores$pair_kinship <- K[cbind(match(scores$parent1, parents), match(scores$parent2, parents))]

n_crosses <- 20L

hi <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K, strategy = "high_gain")
ba <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K, strategy = "balanced")
di <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K, strategy = "diversity")

sh <- attr(hi, "summary")
sb <- attr(ba, "summary")
sd <- attr(di, "summary")

# Each returns a valid plan of the requested size.
stopifnot(nrow(hi) == n_crosses, nrow(ba) == n_crosses, nrow(di) == n_crosses)

# Strategy metadata is recorded.
stopifnot(identical(sh$strategy, "high_gain"),
          isTRUE(all.equal(sh$diversity_emphasis, 15)),
          isTRUE(all.equal(sb$diversity_emphasis, 45)),
          isTRUE(all.equal(sd$diversity_emphasis, 75)))

# Monotonicity across strategies: high_gain should not have LESS gain than diversity,
# and diversity should not have MORE coancestry than high_gain (the whole point of the
# emphasis control). Use >= / <= to tolerate frontier flats.
stopifnot(sh$mean_gain >= sd$mean_gain - 1e-8)
stopifnot(sd$group_coancestry <= sh$group_coancestry + 1e-8)
stopifnot(sb$mean_gain >= sd$mean_gain - 1e-8, sh$mean_gain >= sb$mean_gain - 1e-8)

# Numeric diversity_emphasis path works and achieved emphasis tracks the target.
t70 <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K, diversity_emphasis = 70)
s70 <- attr(t70, "summary")
stopifnot(is.finite(s70$achieved_emphasis), abs(s70$emphasis_gap) <= 25)

# ng_choose_frontier_point extremes: max_gain picks the highest-gain frontier row,
# min_coancestry the lowest-coancestry row.
sweep <- ng_pareto_mate_allocation(scores, n_crosses, parent_kinship = K)
fr <- sweep$frontier
stopifnot(ng_choose_frontier_point(fr, mode = "max_gain") == which.max(fr$mean_gain))
stopifnot(ng_choose_frontier_point(fr, mode = "min_coancestry") ==
            order(fr$group_coancestry, -fr$mean_gain)[[1L]])

# The refactored proxy chooser agrees with the shared mapper for a mid emphasis.
idx_shared <- ng_choose_frontier_point(fr, emphasis = 45, mode = "target")
idx_proxy <- ng_alphamate_style_choose_frontier_plan(fr, mode = "ModeOptTarget1", target_degree = 45)
stopifnot(identical(idx_shared, idx_proxy))

cat("strategy frontier control test passed\n")
