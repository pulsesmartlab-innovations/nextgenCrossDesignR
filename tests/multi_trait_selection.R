helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

scores <- data.frame(
  parent1 = c("A", "A", "B", "D", "D", "E"),
  parent2 = c("B", "C", "C", "E", "F", "F"),
  yield = c(110, 104, 95, 116, 102, 90),
  lodging = c(80, 2, 10, 70, 12, 8),
  disease = c(80, 2, 8, 70, 10, 6),
  pair_kinship = c(0.20, 0.04, 0.06, 0.22, 0.05, 0.08),
  stringsAsFactors = FALSE
)

weighted_traits <- ng_multitrait_spec(
  trait = c("yield", "lodging", "disease"),
  direction = c("maximize", "minimize", "minimize"),
  weight = c(2, 1, 1)
)

weighted <- ng_add_multitrait_score(
  scores = scores,
  traits = weighted_traits,
  method = "weighted"
)
weighted_top <- ng_multitrait_select_topn(
  scores = scores,
  traits = weighted_traits,
  n_crosses = 1L,
  method = "weighted"
)
stopifnot("multi_trait_score" %in% names(weighted))
stopifnot(all(c("multi_trait_yield_z", "multi_trait_lodging_z", "multi_trait_disease_z") %in% names(weighted)))
stopifnot(identical(weighted_top$parent1[[1]], "A"))
stopifnot(identical(weighted_top$parent2[[1]], "C"))

auto_traits <- ng_multitrait_spec(
  trait = c("yield", "lodging", "disease"),
  direction = c("maximize", "minimize", "minimize")
)
auto_scored <- ng_add_multitrait_score(
  scores = scores,
  traits = auto_traits,
  method = "auto"
)
auto_meta <- attr(auto_scored, "multi_trait")
stopifnot(identical(auto_meta$method, "auto"))
stopifnot(length(auto_meta$weights) == 3L)
stopifnot(max(abs(auto_meta$weights - rep(1 / 3, 3))) < 1e-12)

threshold_scores <- data.frame(
  parent1 = c("A", "A", "B", "C"),
  parent2 = c("B", "C", "D", "D"),
  yield = c(125, 112, 104, 100),
  lodging = c(70, 10, 15, 6),
  pair_kinship = c(0.20, 0.04, 0.06, 0.08),
  stringsAsFactors = FALSE
)
threshold_traits <- ng_multitrait_spec(
  trait = c("yield", "lodging"),
  direction = c("maximize", "minimize"),
  weight = c(4, 1),
  max_value = c(NA, 20),
  threshold_weight = c(1, 4)
)
soft_top <- ng_multitrait_select_topn(
  scores = threshold_scores,
  traits = threshold_traits,
  n_crosses = 1L,
  method = "threshold",
  threshold_penalty_weight = 2
)
stopifnot(identical(soft_top$parent1[[1]], "A"))
stopifnot(identical(soft_top$parent2[[1]], "C"))
stopifnot(soft_top$multi_trait_threshold_violation[[1]] == 0)

strict_top <- ng_multitrait_select_topn(
  scores = threshold_scores,
  traits = threshold_traits,
  n_crosses = 3L,
  method = "threshold",
  strict_thresholds = TRUE
)
stopifnot(nrow(strict_top) == 3L)
stopifnot(all(strict_top$lodging <= 20))
stopifnot(all(is.finite(strict_top$multi_trait_score)))

desired_scores <- data.frame(
  parent1 = c("A", "A", "B", "C", "D", "E"),
  parent2 = c("B", "C", "D", "D", "F", "F"),
  yield = c(130, 125, 118, 112, 106, 100),
  disease = c(90, 80, 30, 20, 10, 5),
  dry_matter = c(41, 42, 45, 47, 48, 49),
  pair_kinship = c(0.20, 0.18, 0.06, 0.04, 0.05, 0.07),
  stringsAsFactors = FALSE
)
desired_traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "dry_matter"),
  direction = c("maximize", "minimize", "maximize"),
  desired_change = c(5, 35, 3),
  economic_weight = c(1, 4, 2)
)
desired_cov_names <- desired_traits$trait
P_desired <- diag(c(100, 400, 16))
G_desired <- diag(c(50, 200, 8))
dimnames(P_desired) <- dimnames(G_desired) <- list(desired_cov_names, desired_cov_names)
desired_scored <- ng_add_multitrait_score(
  scores = desired_scores,
  traits = desired_traits,
  method = "desired_gain",
  phenotypic_covariance = P_desired,
  genetic_covariance = G_desired
)
desired_meta <- attr(desired_scored, "multi_trait")
stopifnot(identical(desired_meta$method, "desired_gain"))
stopifnot(all(c("desired_gain_coefficients", "desired_gain_target",
                "desired_gain_predicted_response", "economic_weights") %in% names(desired_meta)))
stopifnot(all(is.finite(desired_meta$desired_gain_coefficients)))
stopifnot(all(names(desired_meta$desired_gain_coefficients) == desired_traits$trait))
stopifnot(desired_meta$economic_weights[["disease"]] > desired_meta$economic_weights[["yield"]])
stopifnot(desired_meta$desired_gain_target[["disease"]] > desired_meta$desired_gain_target[["yield"]])
stopifnot(sum(abs(desired_meta$desired_gain_coefficients - desired_meta$weights)) > 0.05)
target <- desired_meta$desired_gain_target
predicted <- desired_meta$desired_gain_predicted_response
stopifnot(sum(target * predicted) / sqrt(sum(target * target) * sum(predicted * predicted)) > 0.95)

desired_top <- ng_multitrait_select_topn(
  scores = desired_scores,
  traits = desired_traits,
  n_crosses = 1L,
  method = "desired_gain",
  phenotypic_covariance = P_desired,
  genetic_covariance = G_desired
)
weighted_compare_traits <- desired_traits
weighted_compare_traits$weight <- c(1, 4, 2)
weighted_like_top <- ng_multitrait_select_topn(
  scores = desired_scores,
  traits = weighted_compare_traits,
  n_crosses = 1L,
  method = "weighted"
)
stopifnot(desired_top$disease[[1]] < 50)
stopifnot(desired_top$dry_matter[[1]] >= 45)

economic_scored <- ng_add_multitrait_score(
  scores = desired_scores,
  traits = desired_traits,
  method = "economic_index",
  phenotypic_covariance = P_desired,
  genetic_covariance = G_desired
)
economic_meta <- attr(economic_scored, "multi_trait")
stopifnot(identical(economic_meta$method, "economic_index"))
stopifnot(all(c("economic_index_coefficients", "economic_index_target",
                "economic_index_predicted_response", "economic_weights") %in% names(economic_meta)))
stopifnot(all(is.finite(economic_meta$economic_index_coefficients)))
stopifnot(all(names(economic_meta$economic_index_coefficients) == desired_traits$trait))
economic_scales <- vapply(
  desired_traits$column,
  function(nm) ng_multitrait_value_scale(desired_scores[[nm]]), numeric(1)
)
expected_economic_target <- economic_meta$economic_weights * economic_scales
expected_economic_target <- expected_economic_target / sum(abs(expected_economic_target))
stopifnot(max(abs(economic_meta$economic_index_target - expected_economic_target)) < 1e-12)
stopifnot(sum(abs(economic_meta$economic_index_coefficients - economic_meta$weights)) > 0.05)
economic_weighted_scored <- ng_add_multitrait_score(
  scores = desired_scores,
  traits = weighted_compare_traits,
  method = "weighted"
)
stopifnot(max(abs(economic_scored$multi_trait_score - economic_weighted_scored$multi_trait_score)) > 0.01)
stopifnot(sum(economic_meta$economic_index_target * economic_meta$economic_index_predicted_response) > 0)

desired_parents <- sort(unique(c(desired_scores$parent1, desired_scores$parent2)))
desired_parent_K <- diag(length(desired_parents))
rownames(desired_parent_K) <- colnames(desired_parent_K) <- desired_parents
desired_scores_plan <- desired_scores
desired_scores_plan$pair_kinship <- 0
desired_plan <- ng_optimize_multitrait_mating_plan(
  scores = desired_scores_plan,
  traits = desired_traits,
  n_crosses = 2L,
  parent_kinship = desired_parent_K,
  multitrait_method = "desired_gain",
  optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L,
  phenotypic_covariance = P_desired,
  genetic_covariance = G_desired
)
desired_summary <- attr(desired_plan, "summary")
stopifnot(identical(desired_summary$multitrait_method, "desired_gain"))
stopifnot(length(desired_summary$multitrait_desired_gain_coefficients) == nrow(desired_traits))
stopifnot(length(desired_summary$multitrait_desired_gain_target) == nrow(desired_traits))
stopifnot(length(desired_summary$multitrait_desired_gain_predicted_response) == nrow(desired_traits))
stopifnot(all(is.finite(desired_summary$multitrait_desired_gain_coefficients)))
stopifnot(all(names(desired_summary$multitrait_desired_gain_coefficients) == desired_traits$trait))
stopifnot(sum(desired_summary$multitrait_desired_gain_target * desired_summary$multitrait_desired_gain_predicted_response) > 0)

economic_plan <- ng_optimize_multitrait_mating_plan(
  scores = desired_scores_plan,
  traits = desired_traits,
  n_crosses = 2L,
  parent_kinship = desired_parent_K,
  multitrait_method = "economic_index",
  optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L,
  phenotypic_covariance = P_desired,
  genetic_covariance = G_desired
)
economic_summary <- attr(economic_plan, "summary")
stopifnot(identical(economic_summary$multitrait_method, "economic_index"))
stopifnot(length(economic_summary$multitrait_economic_index_coefficients) == nrow(desired_traits))
stopifnot(all(is.finite(economic_summary$multitrait_economic_index_coefficients)))

policy <- ng_multitrait_policy_select()
stopifnot(identical(policy$method, "auto"))
stopifnot(identical(policy$primary_method, "auto"))
stopifnot(identical(policy$family, "rank_threshold"))
stopifnot(!isTRUE(policy$is_fallback))
threshold_policy <- ng_multitrait_policy_select(available_methods = "threshold")
stopifnot(identical(threshold_policy$method, "threshold"))
stopifnot(isTRUE(threshold_policy$is_fallback))
economic_policy <- ng_multitrait_policy_select(method = "economic_index")
stopifnot(identical(economic_policy$method, "economic_index"))
stopifnot(identical(economic_policy$family, "economic_index"))

parents <- sort(unique(c(scores$parent1, scores$parent2)))
parent_kinship <- diag(length(parents))
rownames(parent_kinship) <- colnames(parent_kinship) <- parents
scores_plan <- scores
scores_plan$pair_kinship <- 0
plan <- ng_optimize_multitrait_mating_plan(
  scores = scores_plan,
  traits = weighted_traits,
  n_crosses = 2L,
  parent_kinship = parent_kinship,
  multitrait_method = "weighted",
  optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L
)
summary <- attr(plan, "summary")
stopifnot(nrow(plan) == 2L)
stopifnot(identical(summary$gain_col, "multi_trait_score"))
stopifnot(identical(summary$multitrait_method, "weighted"))
stopifnot(identical(summary$multitrait_trait_count, 3L))

cat("multi-trait selection tests passed\n")
