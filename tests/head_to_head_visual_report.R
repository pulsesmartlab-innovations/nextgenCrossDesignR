# Repository root FIRST. Checking `nextgen_cross_design/` ahead of getwd() made
# these tests load a STALE 0.19.0 copy of the package that sits in the working
# tree under exactly that name -- so they validated a package eleven versions old
# while appearing to cover the current one. The ones that failed were the lucky
# case; the ones that passed gave false assurance. ng_load() already resolves in
# this order; only these hand-rolled preambles inverted it.
root_candidates <- c(
  getwd(),
  file.path(".."),
  file.path(getwd(), "nextgen_cross_design"),
  file.path("..", "nextgen_cross_design")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

summary <- data.frame(
  scenario = rep("cassava_diploid", 4L),
  crop = rep("cassava", 4L),
  harness_model = rep("diploid", 4L),
  n_parents = rep(8L, 4L),
  rep = rep(1L, 4L),
  method = c(
    "nextgen_auto_ocs",
    "nextgen_economic_index_ocs",
    "popvar_style_weighted_topn",
    "alphamate_style_weighted_ocs"
  ),
  benchmark_role = c("candidate", "candidate", "baseline", "baseline"),
  external_tool = c("NextGen", "NextGen", "PopVar", "AlphaMate"),
  implementation = c("native_multitrait", "native_multitrait", "style_proxy", "style_proxy"),
  exact_external_status = c(
    "not_applicable",
    "not_applicable",
    "not_run_multitrait_ci_proxy",
    "not_run_multitrait_ci_proxy"
  ),
  fallback_reason = c(
    "",
    "",
    "PopVar exact external multi-trait run was not available for this fixture.",
    "AlphaMate exact external multi-trait run was not available for this fixture."
  ),
  selected_crosses = rep(2L, 4L),
  mean_realized_index = c(1.18, 1.32, 1.05, 1.20),
  mean_realized_yield = c(114, 116, 112, 115),
  mean_realized_disease = c(31, 29, 35, 34),
  mean_realized_quality = c(42, 43, 40, 41),
  unique_parents = c(4L, 4L, 3L, 4L),
  max_parent_use = c(1L, 1L, 2L, 1L),
  mean_pair_kinship = c(0.07, 0.08, 0.12, 0.09),
  group_coancestry = c(0.06, 0.07, 0.11, 0.08),
  stringsAsFactors = FALSE
)

selections <- data.frame(
  scenario = rep("cassava_diploid", 8L),
  n_parents = rep(8L, 8L),
  rep = rep(1L, 8L),
  method = rep(summary$method, each = 2L),
  benchmark_role = rep(summary$benchmark_role, each = 2L),
  external_tool = rep(summary$external_tool, each = 2L),
  implementation = rep(summary$implementation, each = 2L),
  selection_rank = rep(1:2, 4L),
  parent1 = c("P01", "P02", "P01", "P03", "P01", "P01", "P02", "P04"),
  parent2 = c("P04", "P05", "P04", "P06", "P03", "P04", "P05", "P07"),
  multi_trait_score = c(2.1, 1.9, 2.4, 2.0, 1.8, 1.7, 2.2, 1.8),
  realized_index = c(1.2, 1.16, 1.35, 1.29, 1.08, 1.02, 1.23, 1.17),
  realized_yield = c(115, 113, 117, 115, 112, 112, 116, 114),
  realized_disease = c(30, 32, 29, 29, 34, 36, 33, 35),
  realized_quality = c(42, 42, 43, 43, 40, 40, 41, 41),
  pair_kinship = c(0.05, 0.09, 0.07, 0.09, 0.12, 0.12, 0.08, 0.10),
  stringsAsFactors = FALSE
)

comparisons <- data.frame(
  scenario = rep("cassava_diploid", 6L),
  crop = rep("cassava", 6L),
  harness_model = rep("diploid", 6L),
  n_parents = rep(8L, 6L),
  rep = rep(1L, 6L),
  metric = c(
    "mean_realized_index",
    "mean_realized_yield",
    "mean_realized_disease",
    "mean_realized_index",
    "mean_realized_yield",
    "mean_realized_disease"
  ),
  direction = c("maximize", "maximize", "minimize", "maximize", "maximize", "minimize"),
  method = rep("nextgen_economic_index_ocs", 6L),
  method_role = rep("candidate", 6L),
  method_external_tool = rep("NextGen", 6L),
  method_implementation = rep("native_multitrait", 6L),
  method_exact_external_status = rep("not_applicable", 6L),
  method_fallback_reason = rep("", 6L),
  baseline_method = rep(c("popvar_style_weighted_topn", "alphamate_style_weighted_ocs"), each = 3L),
  baseline_external_tool = rep(c("PopVar", "AlphaMate"), each = 3L),
  baseline_implementation = rep("style_proxy", 6L),
  baseline_exact_external_status = rep("not_run_multitrait_ci_proxy", 6L),
  baseline_fallback_reason = rep("Fixture baseline is a style proxy.", 6L),
  value = c(1.32, 116, 29, 1.32, 116, 29),
  baseline_value = c(1.05, 112, 35, 1.20, 115, 34),
  delta = c(0.27, 4, 6, 0.12, 1, 5),
  better_than_baseline = rep(TRUE, 6L),
  tied_with_baseline = rep(FALSE, 6L),
  stringsAsFactors = FALSE
)

winner_summary <- data.frame(
  n_parents = c(8L, 8L, 8L, 8L),
  metric = c(
    "mean_realized_index",
    "mean_realized_yield",
    "mean_realized_disease",
    "group_coancestry"
  ),
  direction = c("maximize", "maximize", "minimize", "minimize"),
  method = c(
    "nextgen_economic_index_ocs",
    "nextgen_economic_index_ocs",
    "nextgen_economic_index_ocs",
    "nextgen_auto_ocs"
  ),
  tied_methods = c(
    "nextgen_economic_index_ocs",
    "nextgen_economic_index_ocs",
    "nextgen_economic_index_ocs",
    "nextgen_auto_ocs"
  ),
  tied_method_count = rep(1L, 4L),
  value = c(1.32, 116, 29, 0.06),
  reps = rep(1L, 4L),
  stringsAsFactors = FALSE
)

tmp <- tempfile("ng_visual_report_")
dir.create(tmp, recursive = TRUE)
out_path <- file.path(tmp, "fixture_visual_report.html")
returned <- ng_write_head_to_head_visual_report(
  output_path = out_path,
  summary = summary,
  selections = selections,
  comparisons = comparisons,
  winner_summary = winner_summary,
  method_registry = ng_head_to_head_default_methods(),
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)

stopifnot(identical(normalizePath(returned, winslash = "/", mustWork = TRUE),
                    normalizePath(out_path, winslash = "/", mustWork = TRUE)))
html <- paste(readLines(out_path, warn = FALSE), collapse = "\n")
stopifnot(grepl("Head-to-Head Multi-Trait Decision Report", html, fixed = TRUE))
stopifnot(grepl("Selection Index Response", html, fixed = TRUE))
stopifnot(grepl("Trait Direction Response", html, fixed = TRUE))
stopifnot(grepl("Gain-Diversity Frontier", html, fixed = TRUE))
stopifnot(grepl("Parent Contribution", html, fixed = TRUE))
stopifnot(grepl("Mate Allocation Matrix", html, fixed = TRUE))
stopifnot(grepl("Benchmark Evidence Caveats", html, fixed = TRUE))
stopifnot(grepl("data-panel=\"gain-diversity-frontier\"", html, fixed = TRUE))
stopifnot(grepl("Increase yield", html, fixed = TRUE))
stopifnot(grepl("Decrease disease", html, fixed = TRUE))
stopifnot(grepl("style proxy", html, fixed = TRUE))
stopifnot(grepl("P01 x P04", html, fixed = TRUE))
stopifnot(grepl("<svg", html, fixed = TRUE))

prefix <- "fixture"
write.csv(summary, file.path(tmp, paste0(prefix, "_summary.csv")), row.names = FALSE)
write.csv(selections, file.path(tmp, paste0(prefix, "_selections.csv")), row.names = FALSE)
write.csv(comparisons, file.path(tmp, paste0(prefix, "_comparisons.csv")), row.names = FALSE)
write.csv(winner_summary, file.path(tmp, paste0(prefix, "_winner_summary.csv")), row.names = FALSE)
write.csv(ng_head_to_head_default_methods(), file.path(tmp, paste0(prefix, "_method_registry.csv")), row.names = FALSE)

cli_out <- file.path(tmp, "fixture_cli_visual_report.html")
cli_env <- c(
  NG_HEAD_TO_HEAD_PREFIX = prefix,
  NG_HEAD_TO_HEAD_OUTPUT_DIR = tmp,
  NG_HEAD_TO_HEAD_VISUAL_REPORT_OUT = cli_out
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
  normalizePath(file.path(root, "tools", "render_head_to_head_visual_report.R"), winslash = "/", mustWork = FALSE),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cli_output, "status")
if (is.null(status)) status <- 0L
if (!identical(as.integer(status), 0L)) print(cli_output)
stopifnot(length(status) == 1L, is.finite(status), status == 0L)
stopifnot(file.exists(cli_out))
cli_html <- paste(readLines(cli_out, warn = FALSE), collapse = "\n")
stopifnot(grepl("cassava_diploid", cli_html, fixed = TRUE))
stopifnot(grepl("nextgen_economic_index_ocs", cli_html, fixed = TRUE))

cat("head-to-head visual report tests passed\n")
