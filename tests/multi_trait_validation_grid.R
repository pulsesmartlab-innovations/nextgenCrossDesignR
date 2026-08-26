source(file.path("tests", "helper_load.R"))

methods <- c("weighted", "economic_index", "desired_gain")
grid <- ng_run_multitrait_validation_grid(
  parent_sizes = c(8L, 10L),
  reps = 2L,
  n_crosses = 3L,
  seed = 101L,
  methods = methods,
  allocator = "topn"
)

summary <- grid$summary
selections <- grid$selections
winners <- grid$winner_summary
scores <- grid$scores
config <- grid$config

stopifnot(is.data.frame(summary))
stopifnot(is.data.frame(selections))
stopifnot(is.data.frame(winners))
stopifnot(is.data.frame(scores))
stopifnot(is.data.frame(config))
stopifnot(nrow(summary) == length(c(8L, 10L)) * 2L * length(methods))
stopifnot(nrow(selections) == length(c(8L, 10L)) * 2L * length(methods) * 3L)
stopifnot(nrow(scores) == (28L + 45L) * 2L)
stopifnot(nrow(config) == length(c(8L, 10L)) * 2L)
stopifnot(all(c("n_parents", "rep", "seed", "allocator") %in% names(summary)))
stopifnot(all(c("n_parents", "rep", "seed", "allocator") %in% names(selections)))
stopifnot(all(c("n_parents", "rep", "seed", "allocator") %in% names(scores)))
stopifnot(all(c("n_parents", "rep", "seed", "allocator", "methods", "n_crosses") %in% names(config)))
stopifnot(all(summary$method %in% methods))
stopifnot(all(selections$method %in% methods))
stopifnot(all(summary$selected_crosses == 3L))
stopifnot(all(is.finite(summary$mean_realized_index)))
stopifnot(identical(sort(unique(summary$n_parents)), c(8L, 10L)))
stopifnot(identical(sort(unique(summary$rep)), c(1L, 2L)))

expected_metrics <- c(
  "mean_realized_index",
  "mean_realized_yield",
  "mean_realized_disease",
  "mean_realized_quality",
  "unique_parents",
  "max_parent_use"
)
stopifnot(all(expected_metrics %in% winners$metric))
stopifnot(all(winners$method %in% methods))
stopifnot(all(is.finite(winners$value)))
stopifnot(identical(sort(unique(winners$n_parents)), c(8L, 10L)))
stopifnot(all(winners$reps == 2L))

tie_summary <- data.frame(
  method = c("auto", "weighted", "economic_index"),
  method_family = c("weighted", "weighted", "economic_index"),
  n_parents = c(8L, 8L, 8L),
  rep = c(1L, 1L, 1L),
  mean_realized_index = c(1, 1, 0.5),
  mean_realized_yield = c(2, 2, 3),
  mean_realized_disease = c(10, 10, 12),
  mean_realized_quality = c(5, 5, 4),
  unique_parents = c(4, 4, 3),
  max_parent_use = c(2, 2, 3),
  stringsAsFactors = FALSE
)
tie_winners <- ng_multitrait_validation_winner_summary(
  tie_summary,
  metrics = c("mean_realized_index", "mean_realized_yield")
)
index_tie <- tie_winners[tie_winners$metric == "mean_realized_index", , drop = FALSE]
yield_tie <- tie_winners[tie_winners$metric == "mean_realized_yield", , drop = FALSE]
stopifnot(nrow(index_tie) == 1L)
stopifnot(identical(index_tie$method[[1]], "auto"))
stopifnot(identical(index_tie$tied_methods[[1]], "auto,weighted"))
stopifnot(index_tie$tied_method_count[[1]] == 2L)
stopifnot(identical(yield_tie$tied_methods[[1]], "economic_index"))
stopifnot(yield_tie$tied_method_count[[1]] == 1L)

ocs_grid <- ng_run_multitrait_validation_grid(
  parent_sizes = 8L,
  reps = 1L,
  n_crosses = 2L,
  seed = 303L,
  methods = "weighted",
  allocator = "ocs",
  ocs_lambda_group = 0.20
)
stopifnot("lambda_group" %in% names(ocs_grid$summary))
stopifnot(all(ocs_grid$summary$lambda_group > 0))
stopifnot(all(is.finite(ocs_grid$summary$group_coancestry)))

tmp <- tempfile("ng_multitrait_validation_grid_")
dir.create(tmp, recursive = TRUE)
written <- ng_run_multitrait_validation_grid(
  parent_sizes = c(8L, 10L),
  reps = 2L,
  n_crosses = 3L,
  seed = 211L,
  methods = methods,
  allocator = "topn",
  output_dir = tmp,
  prefix = "grid_contract"
)
stopifnot(file.exists(file.path(tmp, "grid_contract_summary.csv")))
stopifnot(file.exists(file.path(tmp, "grid_contract_selections.csv")))
stopifnot(file.exists(file.path(tmp, "grid_contract_scores.csv")))
stopifnot(file.exists(file.path(tmp, "grid_contract_winner_summary.csv")))
stopifnot(file.exists(file.path(tmp, "grid_contract_config.csv")))
stopifnot(nrow(written$winner_summary) == length(c(8L, 10L)) * length(expected_metrics))
stopifnot(is.data.frame(written$scores))
stopifnot(is.data.frame(written$config))

cat("multi-trait validation grid tests passed\n")
