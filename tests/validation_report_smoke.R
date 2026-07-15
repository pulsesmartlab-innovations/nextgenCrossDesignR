helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

tmp <- tempfile("ng_validation_report_")
dir.create(tmp, recursive = TRUE)

write_fixture <- function(name, data) {
  write.csv(data, file.path(tmp, name), row.names = FALSE)
}

write_fixture("diagnostic_family_5k_status.csv", data.frame(
  rep = c(1L, 1L),
  n_parents = c(20L, 30L),
  n_scores = c(80L, 80L),
  popvar_status = c("ok", "ok"),
  simple_status = c("ok", "ok"),
  gms_status = c("ok", "ok"),
  popvar_error = c(NA, NA),
  simple_mpv_error = c(NA, NA),
  simple_usefa_error = c(NA, NA),
  gms_error = c(NA, NA),
  effect_reliability = c(0.42, 0.48),
  in_sample_reliability = c(0.81, 0.83)
))

write_fixture("diagnostic_parent_grid_5k_winner_summary.csv", data.frame(
  n_parents = c(20L, 20L, 20L, 30L, 30L, 30L),
  metric = rep(c("mean_gv", "top10_gv", "max_gv"), 2L),
  method = c(
    "ng_recomb_gebv_ocs10_lps2",
    "ng_recomb_gebv_ocs10_lps2",
    "ng_recomb_gebv_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2"
  ),
  value = c(5.10, 5.90, 6.01, 5.20, 6.05, 6.10),
  top_crosses = c(5L, 5L, 5L, 5L, 5L, 5L),
  effect_training_n = c(400L, 400L, 400L, 400L, 400L, 400L)
))

write_fixture("diagnostic_allocator_5k_winner_summary.csv", data.frame(
  n_parents = c(20L, 30L),
  metric = c("top10_gv", "top10_gv"),
  method = c("ng_ocs_mip10_lps1", "popvar_musp_ocs10_lps1"),
  value = c(4.17, 4.12),
  top_crosses = c(5L, 5L),
  effect_training_n = c(400L, 400L)
))

write_fixture("diagnostic_family_5k_metric_summary_avg.csv", data.frame(
  n_parents = c(20L, 20L, 30L),
  target = c("realized_top10", "realized_var", "realized_top10"),
  score = c("uc_dh", "dh_pmv_var", "uc_dh"),
  spearman = c(0.91, 0.55, 0.89),
  top_overlap = c(0.82, 0.35, 0.80)
))

out_path <- file.path(tmp, "validation_report.md")
script <- normalizePath(file.path(root, "tools", "summarize_validation_report.R"), mustWork = TRUE)
old_env <- Sys.getenv(c("NG_VALIDATION_REPORT_RESULTS_DIR", "NG_VALIDATION_REPORT_OUT"), unset = NA_character_)
on.exit({
  for (name in names(old_env)) {
    if (is.na(old_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, stats::setNames(as.list(old_env[[name]]), name))
    }
  }
}, add = TRUE)
Sys.setenv(
  NG_VALIDATION_REPORT_RESULTS_DIR = normalizePath(tmp, winslash = "/", mustWork = TRUE),
  NG_VALIDATION_REPORT_OUT = normalizePath(out_path, winslash = "/", mustWork = FALSE)
)

source(script, local = new.env(parent = globalenv()))

stopifnot(file.exists(out_path))
report <- readLines(out_path, warn = FALSE)
stopifnot(any(grepl("# Framework Validation Report", report, fixed = TRUE)))
stopifnot(any(grepl("Family status: all 2 calibration runs completed with PopVar, SimpleMating, and GMS status ok.", report, fixed = TRUE)))
stopifnot(any(grepl("Parent-size top10 frontier", report, fixed = TRUE)))
stopifnot(any(grepl("ng_recomb_gebv_ocs10_lps2", report, fixed = TRUE)))
stopifnot(any(grepl("ng_meta_selector_ocs10_lps2", report, fixed = TRUE)))
stopifnot(any(grepl("Do not promote a universal winner from mixed parent-size evidence.", report, fixed = TRUE)))

cat("validation report smoke test passed\n")
