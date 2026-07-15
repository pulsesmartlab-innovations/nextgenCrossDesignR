find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(start, "nextgen_cross_design"),
    file.path(dirname(start), "nextgen_cross_design"),
    start,
    dirname(start),
    file.path(dirname(dirname(start)), "nextgen_cross_design"),
    dirname(dirname(start))
  ), winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "load.R"))) return(candidate)
  }
  stop("Could not locate nextgen_cross_design root containing R/load.R", call. = FALSE)
}

root <- find_project_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  message("AlphaSimR unavailable; skipping poly4x runner smoke test")
  quit(save = "no", status = 0)
}

prefix <- paste0("poly4x_runner_smoke_", Sys.getpid())
script <- file.path(root, "tools", "run_poly4x_benchmark.R")
env_values <- c(
  NG_POLY4X_PREFIX = prefix,
  NG_POLY4X_SCENARIO = "potato_autotetraploid_4x",
  NG_POLY4X_REPS = "1",
  NG_POLY4X_CYCLES = "1",
  NG_POLY4X_N_PARENTS = "8",
  NG_POLY4X_N_FOUNDERS = "14",
  NG_POLY4X_N_CHR = "2",
  NG_POLY4X_SEG_SITES = "24",
  NG_POLY4X_SNP_PER_CHR = "6",
  NG_POLY4X_QTL_PER_CHR = "4",
  NG_POLY4X_TOP_CROSSES = "3",
  NG_POLY4X_SCORE_PROGENY_PER_CROSS = "4",
  NG_POLY4X_REALIZED_PROGENY_PER_CROSS = "5",
  NG_POLY4X_METHODS = "ng_poly4x_var_topn,ng_poly4x_usefulness_topn,ng_poly4x_ocs",
  NG_POLY4X_POLICY_MODES = "gain"
)
old_env <- Sys.getenv(names(env_values), unset = NA_character_)
restore_env <- function() {
  for (name in names(old_env)) {
    if (is.na(old_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old_env[[name]], name)))
    }
  }
}
on.exit(restore_env(), add = TRUE)
on.exit(unlink(file.path("results", paste0(prefix, "_*.csv"))), add = TRUE)

do.call(Sys.setenv, as.list(env_values))
status <- system2("Rscript", c("--vanilla", normalizePath(script, winslash = "/", mustWork = TRUE)))
stopifnot(identical(status, 0L))

families_path <- file.path("results", paste0(prefix, "_families.csv"))
metrics_path <- file.path("results", paste0(prefix, "_metrics.csv"))
selections_path <- file.path("results", paste0(prefix, "_selections.csv"))
summary_path <- file.path("results", paste0(prefix, "_selection_summary.csv"))
stopifnot(file.exists(families_path))
stopifnot(file.exists(metrics_path))
stopifnot(file.exists(selections_path))
stopifnot(file.exists(summary_path))

families <- read.csv(families_path, stringsAsFactors = FALSE)
metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
selections <- read.csv(selections_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)

stopifnot(nrow(families) == 12L)
stopifnot(nrow(metrics) >= 6L)
stopifnot(nrow(selections) == 12L)
stopifnot(nrow(summary) == 4L)
stopifnot(all(families$ploidy == 4L))
stopifnot(all(c("parent1", "parent2", "rep", "cycle", "method", "used_dh") %in% names(selections)))
stopifnot(all(selections$used_dh == 0L))
stopifnot(all(summary$used_dh == 0L))
stopifnot(all(c("ng_poly4x_var_topn", "ng_poly4x_usefulness_topn", "ng_poly4x_ocs") %in% summary$method))
stopifnot("poly4x_policy_gain" %in% unique(selections$method))
stopifnot(all(c("poly4x_policy_mode", "poly4x_policy_source") %in% names(selections)))
policy_rows <- selections[selections$method == "poly4x_policy_gain", ]
stopifnot(nrow(policy_rows) == 3L)
stopifnot(all(policy_rows$poly4x_policy_mode == "gain"))
stopifnot(all(c("pred_poly4x_policy_mode", "pred_poly4x_policy_source") %in% names(families)))
policy_family_rows <- families[families$method == "poly4x_policy_gain", ]
stopifnot(nrow(policy_family_rows) == 3L)
stopifnot(all(policy_family_rows$pred_poly4x_policy_mode == "gain"))
stopifnot(all(is.finite(families$pred_poly4x_usefulness)))
stopifnot(all(is.finite(families$realized_top10_gv)))

cat("poly4x runner smoke tests passed\n")
