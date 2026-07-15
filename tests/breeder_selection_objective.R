helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

scores <- data.frame(
  parent1 = c("A", "A", "B", "C", "D", "E"),
  parent2 = c("B", "C", "D", "D", "F", "F"),
  yield = c(110, 112, 118, 104, 121, 95),
  disease = c(85, 8, 42, 12, 90, 5),
  lodging = c(70, 9, 25, 11, 80, 4),
  quality = c(40, 45, 44, 42, 46, 39),
  pair_kinship = c(0.20, 0.05, 0.08, 0.06, 0.21, 0.07),
  stringsAsFactors = FALSE
)

stopifnot(is.function(get("ng_breeder_selection_objective", mode = "function")))
stopifnot(is.function(get("ng_score_breeder_objective", mode = "function")))
stopifnot(is.function(get("ng_optimize_breeder_selection_plan", mode = "function")))

missing_direction <- tryCatch(
  ng_breeder_selection_objective(
    trait = c("yield", "disease"),
    column = c("yield", "disease")
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("direction", missing_direction))

obj_auto <- ng_breeder_selection_objective(
  trait = c("yield", "disease", "lodging"),
  direction = c("maximize", "minimize", "minimize"),
  method = "auto"
)
stopifnot(inherits(obj_auto, "ng_breeder_selection_objective"))
stopifnot(identical(obj_auto$method, "auto"))
stopifnot(identical(obj_auto$diagnostics$method_reason, "equal_weight_rank_default"))
stopifnot(max(abs(obj_auto$diagnostics$resolved_weights - rep(1 / 3, 3))) < 1e-12)

scored <- ng_score_breeder_objective(scores, obj_auto)
stopifnot("multi_trait_score" %in% names(scored))
stopifnot("multi_trait_threshold_policy" %in% names(scored))
stopifnot(identical(attr(scored, "breeder_objective")$method, "auto"))

top <- scored[order(-scored$multi_trait_score), ][1L, ]
stopifnot(top$disease < 20)
stopifnot(top$lodging < 20)

obj_yield_wrong <- ng_breeder_selection_objective(
  trait = c("yield", "disease"),
  direction = c("minimize", "minimize"),
  method = "auto"
)
wrong <- ng_score_breeder_objective(scores, obj_yield_wrong)
stopifnot(!identical(
  paste(scored$parent1[order(-scored$multi_trait_score)][1:3],
        scored$parent2[order(-scored$multi_trait_score)][1:3]),
  paste(wrong$parent1[order(-wrong$multi_trait_score)][1:3],
        wrong$parent2[order(-wrong$multi_trait_score)][1:3])
))

obj_threshold <- ng_breeder_selection_objective(
  trait = c("yield", "disease"),
  direction = c("maximize", "minimize"),
  weight = c(2, 1),
  min_value = c(100, NA),
  max_value = c(NA, 20),
  threshold_weight = c(1, 4),
  threshold_policy = "soft",
  method = "weighted"
)
soft <- ng_score_breeder_objective(scores, obj_threshold, threshold_penalty_weight = 2)
stopifnot(any(soft$multi_trait_threshold_violation > 0))
stopifnot(all(is.finite(soft$multi_trait_score)))
stopifnot(identical(unique(soft$multi_trait_threshold_policy), "soft"))

obj_strict <- ng_breeder_selection_objective(
  trait = c("yield", "disease"),
  direction = c("maximize", "minimize"),
  min_value = c(120, NA),
  max_value = c(NA, 4),
  threshold_policy = "strict",
  method = "auto"
)
strict_err <- tryCatch(
  ng_optimize_breeder_selection_plan(
    scores = scores,
    objective = obj_strict,
    n_crosses = 2L,
    parent_K = diag(6)
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("feasible", strict_err))

obj_econ <- ng_breeder_selection_objective(
  trait = c("yield", "disease", "quality"),
  direction = c("maximize", "minimize", "maximize"),
  economic_weight = c(1, 4, 2),
  method = "auto"
)
stopifnot(identical(obj_econ$method, "economic_index"))
stopifnot(identical(obj_econ$method_requested, "auto"))
econ_scored <- ng_score_breeder_objective(scores, obj_econ)
econ_meta <- attr(econ_scored, "multi_trait")
stopifnot(all(c("economic_index_coefficients", "economic_index_target") %in% names(econ_meta)))
stopifnot(all(names(econ_meta$economic_index_coefficients) == obj_econ$traits$trait))

obj_desired <- ng_breeder_selection_objective(
  trait = c("yield", "disease", "quality"),
  direction = c("maximize", "minimize", "maximize"),
  desired_change = c(5, 40, 2),
  economic_weight = c(1, 4, 2),
  method = "auto"
)
stopifnot(identical(obj_desired$method, "desired_gain"))
desired_scored <- ng_score_breeder_objective(scores, obj_desired)
desired_meta <- attr(desired_scored, "multi_trait")
stopifnot(all(c("desired_gain_coefficients", "desired_gain_target") %in% names(desired_meta)))
stopifnot(all(names(desired_meta$desired_gain_coefficients) == obj_desired$traits$trait))

parents <- sort(unique(c(scores$parent1, scores$parent2)))
parent_K <- diag(length(parents))
rownames(parent_K) <- colnames(parent_K) <- parents
plan <- ng_optimize_breeder_selection_plan(
  scores = scores,
  objective = obj_threshold,
  n_crosses = 2L,
  parent_K = parent_K,
  optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L
)
summary <- attr(plan, "summary")
stopifnot(nrow(plan) == 2L)
stopifnot(identical(summary$breeder_objective_method, "weighted"))
stopifnot(identical(summary$breeder_objective_threshold_policy, "soft"))
stopifnot(is.finite(summary$breeder_objective_zero_violation_selected))

cat("breeder selection objective tests passed\n")
