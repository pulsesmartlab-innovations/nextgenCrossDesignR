# Merit check for porting ng_local_swap to C++. Measures:
#   1. Wall time of ng_local_swap in a realistic greedy_local OCS run.
#   2. ng_local_swap as a fraction of ng_optimize_mating_plan total time.
#   3. Estimated speedup ceiling (R loop overhead estimate vs intrinsic compute).
# Decides: port if ng_local_swap >20% of total OCS wall time on realistic n.
ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"),
            file.path("nextgen_cross_design", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(20260524L)
n <- 80L; m <- 1500L
ids <- paste0("P", seq_len(n))
markers <- paste0("M", seq_len(m))
geno <- 2L * matrix(rbinom(n * m, 1L, 0.45), n, m, dimnames = list(ids, markers))
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n, sd = 0.5))
names(y) <- ids
mk <- data.frame(marker = markers, chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 100, length.out = m / 5L), 5L))
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 5L, seed = 1L)
parent_kinship <- ng_parent_kinship(geno)
adj <- setNames(rnorm(n), ids)
scores <- ng_score_crosses(geno = geno, effects = fit, marker_map = mk, ids = ids,
                           adjusted_pheno = adj, selection_prop = 0.10,
                           recomb_model = "haldane", target = "DH", use_cpp = TRUE)
scores$.linear_gain <- scores$usefulness_pmv_gebv
n_pairs <- nrow(scores)
cat(sprintf("Fixture: n=%d parents, m=%d markers, %d candidate pairs.\n", n, m, n_pairs))

# ---- Cost A: ng_local_swap alone on a seeded greedy plan ------------------
ord <- order(scores$.linear_gain, decreasing = TRUE)
selected_init <- integer(0)
counts <- setNames(integer(n), ids)
for (idx in ord) {
  p <- c(scores$parent1[idx], scores$parent2[idx])
  if (any(counts[p] >= 4L)) next
  selected_init <- c(selected_init, idx)
  counts[p] <- counts[p] + 1L
  if (length(selected_init) == 20L) break
}
t_swap <- system.time(
  ng_local_swap(scores = scores, selected = selected_init,
                parents = ids, parent_kinship = parent_kinship,
                max_crosses_per_parent = 4L, lambda_group = 1,
                local_iter = 2000L)
)["elapsed"]
cat(sprintf("ng_local_swap (n_crosses=20, local_iter=2000, lambda_group=1): %.2fs\n", t_swap))

# ---- Cost B: full ng_optimize_mating_plan(method='greedy_local') ----------
t_total <- system.time(
  ng_optimize_mating_plan(scores = scores, n_crosses = 20L, parent_kinship = parent_kinship,
                          gain_col = "usefulness_pmv_gebv", method = "greedy_local",
                          max_crosses_per_parent = 4L, lambda_group = 1,
                          local_iter = 2000L)
)["elapsed"]
cat(sprintf("ng_optimize_mating_plan(method=greedy_local) total: %.2fs\n", t_total))
cat(sprintf("ng_local_swap share: %.0f%% of greedy_local total\n",
            100 * t_swap / max(t_total, 1e-6)))

# ---- Cost C: how often is greedy_local actually used? ---------------------
# In the v0.2.0 frontier-policy benchmark grid, the OCS methods (ocs_lps,
# meta_router, frontier_policy, etc.) use method='mip_contribution' by
# default. greedy_local is invoked for var_simple_topn fallback and the
# ng_allocator (no longer in the frontier method set). So this matters most
# for users who explicitly request method='greedy_local' or who hit the
# lpSolve-unavailable fallback.
cat("\n## Per-iteration objective cost breakdown\n")
t_obj <- system.time(replicate(10000L,
  ng_plan_objective(scores, selected_init, parent_kinship, lambda_group = 1)
))["elapsed"]
cat(sprintf("ng_plan_objective x 10000 calls (n_crosses=20, n_parents=%d): %.2fs\n", n, t_obj))
cat(sprintf("Per call: %.2f microseconds\n", t_obj * 1e6 / 10000))

# ---- Decision ------------------------------------------------------------
share_pct <- 100 * t_swap / max(t_total, 1e-6)
cat("\n## Merit verdict\n")
if (share_pct < 30) {
  cat(sprintf("VERDICT: ng_local_swap is %.0f%% of greedy_local total.\n", share_pct))
  cat("Most time is in setup (kinship, etc.), not the swap loop. Marginal win.\n")
} else if (share_pct < 60) {
  cat(sprintf("VERDICT: ng_local_swap is %.0f%% of greedy_local total — moderate merit.\n", share_pct))
  cat("Port worthwhile if greedy_local is a common code path; skip otherwise.\n")
} else {
  cat(sprintf("VERDICT: ng_local_swap dominates greedy_local at %.0f%% — port it.\n", share_pct))
}
