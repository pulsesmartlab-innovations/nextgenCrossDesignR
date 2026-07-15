helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

scenario <- ng_multitrait_validation_scenario(n_parents = 12L, seed = 17L)
stopifnot(is.data.frame(scenario$scores))
stopifnot(nrow(scenario$scores) == 66L)
stopifnot(all(c(
  "parent1", "parent2", "pred_yield", "pred_disease", "pred_quality",
  "realized_yield", "realized_disease", "realized_quality", "pair_kinship"
) %in% names(scenario$scores)))
stopifnot(is.data.frame(scenario$traits))
stopifnot(identical(scenario$traits$trait, c("yield", "disease", "quality")))

result <- ng_run_multitrait_validation(
  n_parents = 12L,
  n_crosses = 3L,
  seed = 17L,
  methods = c("auto", "weighted", "economic_index", "desired_gain", "threshold"),
  allocator = "topn"
)
summary <- result$summary
selections <- result$selections
stopifnot(is.data.frame(summary))
stopifnot(is.data.frame(selections))
stopifnot(all(c("auto", "weighted", "economic_index", "desired_gain", "threshold") %in% summary$method))
stopifnot(all(summary$selected_crosses == 3L))
stopifnot(all(is.finite(summary$mean_realized_index)))
stopifnot(diff(range(summary$mean_realized_index)) > 0.01)
stopifnot(all(c(
  "mean_realized_yield", "mean_realized_disease", "mean_realized_quality",
  "mean_predicted_score", "unique_parents", "max_parent_use"
) %in% names(summary)))
stopifnot(nrow(selections) == 15L)
stopifnot(all(c("method", "selection_rank", "realized_index") %in% names(selections)))

weighted_pairs <- selections[selections$method == "weighted", c("parent1", "parent2")]
economic_pairs <- selections[selections$method == "economic_index", c("parent1", "parent2")]
stopifnot(!identical(
  paste(weighted_pairs$parent1, weighted_pairs$parent2),
  paste(economic_pairs$parent1, economic_pairs$parent2)
))
economic_row <- summary[summary$method == "economic_index", , drop = FALSE]
desired_row <- summary[summary$method == "desired_gain", , drop = FALSE]
stopifnot(is.finite(economic_row$coef_disease[[1]]))
stopifnot(is.finite(desired_row$coef_disease[[1]]))

tmp <- tempfile("ng_multitrait_validation_")
dir.create(tmp, recursive = TRUE)
written <- ng_run_multitrait_validation(
  n_parents = 10L,
  n_crosses = 3L,
  seed = 23L,
  methods = c("weighted", "economic_index"),
  allocator = "topn",
  output_dir = tmp,
  prefix = "contract"
)
stopifnot(file.exists(file.path(tmp, "contract_summary.csv")))
stopifnot(file.exists(file.path(tmp, "contract_selections.csv")))
stopifnot(file.exists(file.path(tmp, "contract_scores.csv")))
stopifnot(nrow(written$summary) == 2L)

cat("multi-trait validation smoke tests passed\n")
