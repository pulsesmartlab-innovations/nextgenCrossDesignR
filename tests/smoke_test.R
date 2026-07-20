helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(10)
n <- 40
m <- 300
ids <- paste0("P", seq_len(n))
markers <- paste0("M", seq_len(m))
geno <- matrix(2L * rbinom(n * m, 1, 0.45), nrow = n, ncol = m,
               dimnames = list(ids, markers))
true_beta <- rnorm(m, sd = 0.08)
y <- as.numeric(geno %*% true_beta + rnorm(n, sd = 1.0))
names(y) <- ids
marker_map <- data.frame(
  marker = markers,
  chr = rep(seq_len(5), length.out = m),
  pos_cm = rep(seq(0, 120, length.out = ceiling(m / 5)), 5)[seq_len(m)]
)

fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 3, seed = 10)
scores <- ng_score_crosses(geno, fit, marker_map = marker_map, ids = ids,
                           adjusted_pheno = y, selection_prop = 0.10,
                           use_cpp = FALSE)
stopifnot(nrow(scores) == n * (n - 1) / 2)
stopifnot(all(is.finite(scores$usefulness_pmv_gebv)))
stopifnot(all(scores$pmv >= -1e-8))
stopifnot(all(scores$pmv + 1e-12 >= scores$vpm))  # PMV >= VPM
scores$parent_distance_cal <- scores$parent_distance
scores$vpm_cal <- scores$vpm
scores$pmv_cal <- scores$pmv
scores$etk_parent_distance_cal <- scores$cross_mean + ng_selection_intensity(0.10) * sqrt(pmax(scores$parent_distance_cal, 0))
scores$etk_vpm_gebv_cal <- scores$cross_mean_gebv + ng_selection_intensity(0.10) * sqrt(pmax(scores$vpm_cal, 0))
scores$etk_vpm_blend_cal <- scores$cross_mean_blend + ng_selection_intensity(0.10) * sqrt(pmax(scores$vpm_cal, 0))
scores$etk_pmv_blend_cal <- scores$cross_mean_blend + ng_selection_intensity(0.10) * sqrt(pmax(scores$pmv_cal, 0))
scores <- ng_add_adaptive_stack_scores(scores, n_parents = n, n_crosses = 12)
stopifnot(all(is.finite(scores$ng_adaptive_score)))
stopifnot(all(is.finite(scores$ng_adaptive_score_raw)))
stopifnot(all(is.finite(scores$ng_adaptive_var)))
stopifnot(all(is.finite(scores$ng_adaptive_fallback_weight)))
stopifnot(is.finite(scores$ng_adaptive_parent_use_input[1]))
scores <- ng_add_meta_portfolio_scores(scores, n_parents = n, n_crosses = 12)
stopifnot(all(is.finite(scores$ng_meta_score)))
stopifnot(all(is.finite(scores$ng_meta_var)))
stopifnot(all(is.finite(scores$ng_meta_method_history_weight)))
stopifnot(all(is.finite(scores$ng_meta_leader_score)))
stopifnot(all(is.finite(scores$ng_meta_leader_var)))
stopifnot(all(nzchar(scores$ng_meta_leader_col)))
stopifnot(is.finite(scores$ng_meta_parent_use_input[1]))

K <- ng_parent_kinship(geno)
plan <- ng_optimize_mating_plan(scores, n_crosses = 12, parent_kinship = K,
                                max_crosses_per_parent = 3,
                                lambda_group = 1,
                                method = "greedy_local")
s <- attr(plan, "summary")
stopifnot(nrow(plan) == 12)
stopifnot(s$max_parent_use <= 3)
print(s[c("mean_gain", "group_coancestry", "unique_parents", "max_parent_use")])

ocs <- ng_optimize_mating_plan(scores, n_crosses = 12, parent_kinship = K,
                               gain_col = "usefulness_pmv_gebv",
                               max_crosses_per_parent = 5,
                               lambda_group = 2,
                               lambda_parent_use = 2,
                               lambda_parent_use_mode = "adaptive",
                               method = "mip_contribution",
                               ocs_iter = 3)
ocs_s <- attr(ocs, "summary")
stopifnot(nrow(ocs) == 12)
stopifnot(ocs_s$max_parent_use <= 5)
stopifnot(is.finite(ocs_s$parent_use_sq))
stopifnot(ocs_s$lambda_parent_use_mode == "adaptive")
stopifnot(is.finite(ocs_s$score_scale))
ocs_fam <- ng_allocate_family_sizes(ocs, total_progeny = 120,
                                    min_progeny = 4,
                                    max_progeny = 20,
                                    value_col = "usefulness_pmv_gebv")
stopifnot(sum(ocs_fam$n_progeny) == 120)
stopifnot(min(ocs_fam$n_progeny) >= 4)
stopifnot(max(ocs_fam$n_progeny) <= 20)
marg_fam <- ng_allocate_family_sizes(ocs, total_progeny = 120,
                                     min_progeny = 4,
                                     max_progeny = 20,
                                     method = "marginal_topk",
                                     mean_col = "cross_mean",
                                     var_col = "pmv",
                                     selected_top_n = 20)
stopifnot(sum(marg_fam$n_progeny) == 120)
stopifnot(min(marg_fam$n_progeny) >= 4)
stopifnot(max(marg_fam$n_progeny) <= 20)
print(ocs_s[c("mean_gain", "group_coancestry", "parent_use_sq", "unique_parents", "max_parent_use")])

adaptive <- ng_optimize_mating_plan(scores, n_crosses = 12, parent_kinship = K,
                                    gain_col = "ng_adaptive_score",
                                    max_crosses_per_parent = 5,
                                    lambda_group = 1,
                                    lambda_parent_use = scores$ng_adaptive_parent_use_input[1],
                                    lambda_parent_use_mode = "adaptive",
                                    method = "mip_contribution",
                                    ocs_iter = 3)
adaptive_s <- attr(adaptive, "summary")
stopifnot(nrow(adaptive) == 12)
stopifnot(adaptive_s$max_parent_use <= 5)
stopifnot(is.finite(adaptive_s$score_scale))
print(adaptive_s[c("mean_gain", "group_coancestry", "parent_use_sq", "unique_parents", "max_parent_use")])

meta <- ng_optimize_mating_plan(scores, n_crosses = 12, parent_kinship = K,
                                gain_col = "ng_meta_score",
                                max_crosses_per_parent = 5,
                                lambda_group = 1,
                                lambda_parent_use = scores$ng_meta_parent_use_input[1],
                                lambda_parent_use_mode = "adaptive",
                                method = "mip_contribution",
                                ocs_iter = 3)
meta_s <- attr(meta, "summary")
stopifnot(nrow(meta) == 12)
stopifnot(meta_s$max_parent_use <= 5)
stopifnot(is.finite(meta_s$score_scale))
print(meta_s[c("mean_gain", "group_coancestry", "parent_use_sq", "unique_parents", "max_parent_use")])

balanced <- ng_optimize_balanced_usefulness(scores, n_crosses = 12,
                                            gain_col = "usefulness_vpm_gebv",
                                            diversity_col = "parent_distance_cal",
                                            parent_kinship = K,
                                            max_crosses_per_parent = 5,
                                            lambda_group = 1,
                                            lambda_parent_use = 2,
                                            lambda_parent_use_mode = "adaptive",
                                            method = "mip_contribution",
                                            ocs_iter = 3)
balanced_s <- attr(balanced, "summary")
stopifnot(nrow(balanced) == 12)
stopifnot(balanced_s$max_parent_use <= 5)
stopifnot(balanced_s$balanced_gain_col == "usefulness_vpm_gebv")
stopifnot(balanced_s$balanced_diversity_col == "parent_distance_cal")
stopifnot(is.na(balanced_s$balanced_min_unique_used))
stopifnot(all(is.finite(balanced$.balanced_usefulness_gain)))
print(balanced_s[c("mean_gain", "group_coancestry", "parent_use_sq", "unique_parents", "max_parent_use")])
