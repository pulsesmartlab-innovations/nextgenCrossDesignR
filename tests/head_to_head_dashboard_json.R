root_candidates <- c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("jsonlite unavailable; skipping dashboard JSON contract test")
  quit(status = 0)
}

summary <- data.frame(
  scenario = rep("cassava_diploid", 3L),
  crop = rep("cassava", 3L),
  harness_model = rep("diploid", 3L),
  n_parents = rep(8L, 3L),
  rep = rep(1L, 3L),
  method = c("nextgen_auto_ocs", "nextgen_economic_index_ocs", "popvar_style_weighted_topn"),
  benchmark_role = c("candidate", "candidate", "baseline"),
  external_tool = c("NextGen", "NextGen", "PopVar"),
  implementation = c("native_multitrait", "native_multitrait", "style_proxy"),
  exact_external_status = c("not_applicable", "not_applicable", "not_run_multitrait_ci_proxy"),
  fallback_reason = c("", "", "PopVar exact external multi-trait run was not available."),
  selected_crosses = rep(2L, 3L),
  mean_realized_index = c(1.18, 1.32, 1.05),
  mean_realized_yield = c(114, 116, 112),
  mean_realized_disease = c(31, 29, 35),
  mean_realized_quality = c(42, 43, 40),
  unique_parents = c(4L, 4L, 3L),
  max_parent_use = c(1L, 1L, 2L),
  mean_pair_kinship = c(0.07, 0.08, 0.12),
  group_coancestry = c(0.06, 0.07, 0.11),
  stringsAsFactors = FALSE
)

selections <- data.frame(
  scenario = rep("cassava_diploid", 6L),
  n_parents = rep(8L, 6L),
  rep = rep(1L, 6L),
  method = rep(summary$method, each = 2L),
  benchmark_role = rep(summary$benchmark_role, each = 2L),
  external_tool = rep(summary$external_tool, each = 2L),
  implementation = rep(summary$implementation, each = 2L),
  selection_rank = rep(1:2, 3L),
  parent1 = c("P01", "P02", "P01", "P03", "P01", "P01"),
  parent2 = c("P04", "P05", "P04", "P06", "P03", "P04"),
  multi_trait_score = c(2.1, 1.9, 2.4, 2.0, 1.8, 1.7),
  realized_index = c(1.2, 1.16, 1.35, 1.29, 1.08, 1.02),
  realized_yield = c(115, 113, 117, 115, 112, 112),
  realized_disease = c(30, 32, 29, 29, 34, 36),
  realized_quality = c(42, 42, 43, 43, 40, 40),
  pair_kinship = c(0.05, 0.09, 0.07, 0.09, 0.12, 0.12),
  stringsAsFactors = FALSE
)

comparisons <- data.frame(
  scenario = rep("cassava_diploid", 3L),
  crop = rep("cassava", 3L),
  harness_model = rep("diploid", 3L),
  n_parents = rep(8L, 3L),
  rep = rep(1L, 3L),
  metric = c("mean_realized_index", "mean_realized_yield", "mean_realized_disease"),
  direction = c("maximize", "maximize", "minimize"),
  method = rep("nextgen_economic_index_ocs", 3L),
  method_role = rep("candidate", 3L),
  method_external_tool = rep("NextGen", 3L),
  method_implementation = rep("native_multitrait", 3L),
  method_exact_external_status = rep("not_applicable", 3L),
  method_fallback_reason = rep("", 3L),
  baseline_method = rep("popvar_style_weighted_topn", 3L),
  baseline_external_tool = rep("PopVar", 3L),
  baseline_implementation = rep("style_proxy", 3L),
  baseline_exact_external_status = rep("not_run_multitrait_ci_proxy", 3L),
  baseline_fallback_reason = rep("Fixture baseline is a style proxy.", 3L),
  value = c(1.32, 116, 29),
  baseline_value = c(1.05, 112, 35),
  delta = c(0.27, 4, 6),
  better_than_baseline = rep(TRUE, 3L),
  tied_with_baseline = rep(FALSE, 3L),
  stringsAsFactors = FALSE
)

winner_summary <- data.frame(
  n_parents = c(8L, 8L),
  metric = c("mean_realized_index", "mean_realized_disease"),
  direction = c("maximize", "minimize"),
  method = c("nextgen_economic_index_ocs", "nextgen_economic_index_ocs"),
  tied_methods = c("nextgen_economic_index_ocs", "nextgen_economic_index_ocs"),
  tied_method_count = c(1L, 1L),
  value = c(1.32, 29),
  reps = c(1L, 1L),
  stringsAsFactors = FALSE
)

config <- data.frame(
  scenario = "cassava_diploid",
  crop = "cassava",
  harness_model = "diploid",
  n_parents = 8L,
  n_crosses = 2L,
  methods = paste(summary$method, collapse = ","),
  baseline_methods = "popvar_style_weighted_topn",
  stringsAsFactors = FALSE
)

tmp <- tempfile("ng_dashboard_json_")
dir.create(tmp, recursive = TRUE)
out_path <- file.path(tmp, "dashboard.json")
returned <- ng_write_head_to_head_dashboard_json(
  output_path = out_path,
  summary = summary,
  selections = selections,
  comparisons = comparisons,
  winner_summary = winner_summary,
  method_registry = ng_head_to_head_default_methods(),
  config = config,
  prefix = "fixture",
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)

stopifnot(identical(normalizePath(returned, winslash = "/", mustWork = TRUE),
                    normalizePath(out_path, winslash = "/", mustWork = TRUE)))
payload <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
stopifnot(identical(payload$schema_version, "ng_head_to_head_dashboard.v1"))
stopifnot(identical(payload$source$prefix, "fixture"))
stopifnot(identical(payload$metric_directions$mean_realized_disease, "minimize"))
stopifnot(identical(payload$context$scenarios[[1]], "cassava_diploid"))
stopifnot(identical(payload$context$crops[[1]], "cassava"))
stopifnot(identical(payload$recommended_method, "nextgen_economic_index_ocs"))
stopifnot(length(payload$method_summary) == 3L)
stopifnot(length(payload$summary) == nrow(summary))
stopifnot(length(payload$selections) == nrow(selections))
stopifnot(length(payload$comparisons) == nrow(comparisons))
stopifnot(length(payload$winner_summary) == nrow(winner_summary))
stopifnot(length(payload$method_registry) >= 3L)
stopifnot(length(payload$config) == 1L)
stopifnot(identical(payload$method_summary[[1]]$method, "nextgen_economic_index_ocs"))
stopifnot(isTRUE(all.equal(as.numeric(payload$comparisons[[3]]$delta), 6)))

prefix <- "fixture"
write.csv(summary, file.path(tmp, paste0(prefix, "_summary.csv")), row.names = FALSE)
write.csv(selections, file.path(tmp, paste0(prefix, "_selections.csv")), row.names = FALSE)
write.csv(comparisons, file.path(tmp, paste0(prefix, "_comparisons.csv")), row.names = FALSE)
write.csv(winner_summary, file.path(tmp, paste0(prefix, "_winner_summary.csv")), row.names = FALSE)
write.csv(ng_head_to_head_default_methods(), file.path(tmp, paste0(prefix, "_method_registry.csv")), row.names = FALSE)
write.csv(config, file.path(tmp, paste0(prefix, "_config.csv")), row.names = FALSE)

cli_out <- file.path(tmp, "fixture_cli_dashboard.json")
cli_env <- c(
  NG_HEAD_TO_HEAD_PREFIX = prefix,
  NG_HEAD_TO_HEAD_OUTPUT_DIR = tmp,
  NG_HEAD_TO_HEAD_DASHBOARD_JSON_OUT = cli_out
)
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
old_cli_env <- Sys.getenv(names(cli_env), unset = NA_character_)
restore_cli_env <- function() {
  for (name in names(old_cli_env)) {
    if (is.na(old_cli_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old_cli_env[[name]], name)))
    }
  }
}
on.exit(restore_cli_env(), add = TRUE)
do.call(Sys.setenv, as.list(cli_env))
cli_output <- system2(
  rscript,
  normalizePath(file.path(root, "tools", "export_head_to_head_dashboard_json.R"), winslash = "/", mustWork = FALSE),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cli_output, "status")
if (is.null(status)) status <- 0L
if (!identical(as.integer(status), 0L)) print(cli_output)
stopifnot(length(status) == 1L, is.finite(status), status == 0L)
stopifnot(file.exists(cli_out))
cli_payload <- jsonlite::fromJSON(cli_out, simplifyVector = FALSE)
stopifnot(identical(cli_payload$schema_version, "ng_head_to_head_dashboard.v1"))
stopifnot(identical(cli_payload$recommended_method, "nextgen_economic_index_ocs"))

cat("head-to-head dashboard JSON tests passed\n")
