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

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  message("AlphaSimR unavailable; skipping head-to-head benchmark test")
  quit(status = 0)
}

methods <- ng_head_to_head_default_methods()
stopifnot(is.data.frame(methods))
stopifnot(all(c(
  "label", "method", "allocator", "role", "external_tool", "implementation", "fallback_reason"
) %in% names(methods)))
stopifnot(all(c(
  "nextgen_auto_ocs",
  "nextgen_economic_index_ocs",
  "nextgen_desired_gain_ocs",
  "popvar_style_weighted_topn",
  "simplemate_style_threshold_topn",
  "alphamate_style_weighted_ocs"
) %in% methods$label))
stopifnot(all(methods$role[methods$external_tool != "NextGen"] == "baseline"))
stopifnot(all(methods$implementation[methods$external_tool != "NextGen"] == "style_proxy"))
stopifnot(all(nzchar(methods$fallback_reason[methods$implementation == "style_proxy"])))
stopifnot("SimpleMating" %in% methods$external_tool)

tmp <- tempfile("ng_head_to_head_")
dir.create(tmp, recursive = TRUE)
bench <- ng_run_head_to_head_benchmark(
  scenarios = c("compact_selfing", "potato_tetraploid_stress"),
  parent_sizes = 6L,
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 91L,
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L,
  output_dir = tmp,
  prefix = "head_to_head_contract"
)

summary <- bench$summary
selections <- bench$selections
scores <- bench$scores
comparisons <- bench$comparisons
winners <- bench$winner_summary
config <- bench$config
registry <- bench$method_registry

stopifnot(is.data.frame(summary))
stopifnot(is.data.frame(selections))
stopifnot(is.data.frame(scores))
stopifnot(is.data.frame(comparisons))
stopifnot(is.data.frame(winners))
stopifnot(is.data.frame(config))
stopifnot(is.data.frame(registry))
stopifnot(nrow(summary) == 2L * nrow(methods))
stopifnot(nrow(selections) == nrow(summary) * 2L)
stopifnot(nrow(scores) == 2L * choose(6L, 2L))
stopifnot(nrow(config) == 2L)
stopifnot(nrow(registry) == nrow(methods))
stopifnot(all(summary$method %in% methods$label))
stopifnot(all(summary$benchmark_role %in% c("candidate", "baseline")))
stopifnot(all(c("external_tool", "implementation", "allocator", "source_method") %in% names(summary)))
stopifnot(all(c("fallback_reason") %in% names(summary)))
stopifnot(all(nzchar(summary$fallback_reason[summary$implementation == "style_proxy"])))
stopifnot(all(summary$selected_crosses == 2L))
stopifnot(all(is.finite(summary$mean_realized_index)))
stopifnot(all(c(
  "scenario", "n_parents", "rep", "method", "baseline_method", "metric", "delta",
  "method_implementation", "method_exact_external_status", "method_fallback_reason",
  "baseline_implementation", "baseline_exact_external_status", "baseline_fallback_reason"
) %in% names(comparisons)))
stopifnot(all(comparisons$method_role == "candidate"))
stopifnot(all(comparisons$method_implementation == "native_multitrait"))
stopifnot(all(comparisons$baseline_implementation == "style_proxy"))
stopifnot(all(nzchar(comparisons$baseline_fallback_reason)))
stopifnot(any(comparisons$baseline_method == "popvar_style_weighted_topn"))
stopifnot(any(comparisons$baseline_method == "simplemate_style_threshold_topn"))
stopifnot(any(comparisons$baseline_method == "alphamate_style_weighted_ocs"))
stopifnot(all(is.finite(comparisons$delta)))
stopifnot(all(comparisons$direction %in% c("maximize", "minimize")))
stopifnot(all(c("tied_methods", "tied_method_count") %in% names(winners)))
stopifnot(all(c("methods", "baseline_methods", "n_crosses") %in% names(config)))

expected_files <- paste0(
  "head_to_head_contract_",
  c("summary", "selections", "scores", "comparisons", "winner_summary", "config", "method_registry"),
  ".csv"
)
stopifnot(all(file.exists(file.path(tmp, expected_files))))

cli_tmp <- tempfile("ng_head_to_head_cli_")
dir.create(cli_tmp, recursive = TRUE)
cli_prefix <- "head_to_head_cli_contract"
cli_env <- c(
  NG_HEAD_TO_HEAD_PREFIX = cli_prefix,
  NG_HEAD_TO_HEAD_OUTPUT_DIR = cli_tmp,
  NG_HEAD_TO_HEAD_SCENARIOS = "compact_selfing",
  NG_HEAD_TO_HEAD_PARENT_SIZES = "6",
  NG_HEAD_TO_HEAD_REPS = "1",
  NG_HEAD_TO_HEAD_CROSSES = "2",
  NG_HEAD_TO_HEAD_REALIZED_PROGENY = "3",
  NG_HEAD_TO_HEAD_SEED = "97",
  NG_HEAD_TO_HEAD_METHODS = "nextgen_auto_ocs,popvar_style_weighted_topn,alphamate_style_weighted_ocs",
  NG_HEAD_TO_HEAD_N_FOUNDERS = "10",
  NG_HEAD_TO_HEAD_N_CHR = "2",
  NG_HEAD_TO_HEAD_SEG_SITES = "30",
  NG_HEAD_TO_HEAD_SNP_PER_CHR = "8",
  NG_HEAD_TO_HEAD_QTL_PER_CHR = "2"
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
  normalizePath(file.path(root, "tools", "run_head_to_head_benchmark.R"), winslash = "/", mustWork = FALSE),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cli_output, "status")
if (is.null(status)) status <- 0L
if (!identical(as.integer(status), 0L)) print(cli_output)
stopifnot(length(status) == 1L, is.finite(status), status == 0L)
stopifnot(file.exists(file.path(cli_tmp, paste0(cli_prefix, "_comparisons.csv"))))
cli_summary <- read.csv(file.path(cli_tmp, paste0(cli_prefix, "_summary.csv")), stringsAsFactors = FALSE)
cli_registry <- read.csv(file.path(cli_tmp, paste0(cli_prefix, "_method_registry.csv")), stringsAsFactors = FALSE)
stopifnot(nrow(cli_summary) == 3L)
stopifnot(identical(cli_registry$label, c(
  "nextgen_auto_ocs",
  "popvar_style_weighted_topn",
  "alphamate_style_weighted_ocs"
)))

cat("head-to-head benchmark tests passed\n")
