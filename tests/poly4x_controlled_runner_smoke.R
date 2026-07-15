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
  message("AlphaSimR unavailable; skipping poly4x controlled runner smoke test")
  quit(save = "no", status = 0)
}

prefix <- paste0("poly4x_controlled_runner_smoke_", Sys.getpid())
script <- file.path(root, "tools", "run_poly4x_controlled_ocs.R")
env_values <- c(
  NG_POLY4X_CONTROLLED_PREFIX = prefix,
  NG_POLY4X_CONTROLLED_SCENARIO = "potato_autotetraploid_4x",
  NG_POLY4X_CONTROLLED_REPS = "1",
  NG_POLY4X_CONTROLLED_CYCLES = "2",
  NG_POLY4X_CONTROLLED_N_PARENTS = "8",
  NG_POLY4X_CONTROLLED_N_FOUNDERS = "14",
  NG_POLY4X_CONTROLLED_N_CHR = "2",
  NG_POLY4X_CONTROLLED_SEG_SITES = "24",
  NG_POLY4X_CONTROLLED_SNP_PER_CHR = "6",
  NG_POLY4X_CONTROLLED_QTL_PER_CHR = "4",
  NG_POLY4X_CONTROLLED_TOP_CROSSES = "3",
  NG_POLY4X_CONTROLLED_SCORE_PROGENY_PER_CROSS = "4",
  NG_POLY4X_CONTROLLED_REALIZED_PROGENY_PER_CROSS = "5",
  NG_POLY4X_CONTROLLED_CONFIGS = "usefulness,default,strongdiv,gain,diversity,ocs",
  NG_POLY4X_CONTROLLED_ADVANCE_CONFIG = "gain"
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

stopifnot(file.exists(script))

do.call(Sys.setenv, as.list(env_values))
status <- system2("Rscript", c("--vanilla", normalizePath(script, winslash = "/", mustWork = TRUE)))
stopifnot(identical(status, 0L))

metrics_path <- file.path("results", paste0(prefix, "_metrics.csv"))
selections_path <- file.path("results", paste0(prefix, "_selections.csv"))
summary_path <- file.path("results", paste0(prefix, "_selection_summary.csv"))
families_path <- file.path("results", paste0(prefix, "_families.csv"))
stopifnot(file.exists(metrics_path))
stopifnot(file.exists(selections_path))
stopifnot(file.exists(summary_path))
stopifnot(file.exists(families_path))

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
selections <- read.csv(selections_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)
families <- read.csv(families_path, stringsAsFactors = FALSE)

stopifnot(nrow(metrics) >= 9L)
stopifnot(nrow(selections) == 36L)
stopifnot(nrow(summary) == 12L)
stopifnot(nrow(families) == 36L)
stopifnot(all(metrics$ploidy == 4L))
stopifnot(all(families$ploidy == 4L))
stopifnot(all(selections$used_dh == 0L))
stopifnot(all(summary$used_dh == 0L))
stopifnot(all(c("ng_poly4x_usefulness_topn", "ng_poly4x_ocs") %in% unique(selections$method)))
stopifnot(all(c("usefulness", "default", "strongdiv", "gain", "diversity", "ocs") %in% unique(selections$config)))
stopifnot(all(c("advance_config") %in% names(metrics)))
stopifnot(all(c("advance_config") %in% names(summary)))
stopifnot(all(c("advance_config") %in% names(families)))
stopifnot(all(metrics$advance_config == "gain"))
stopifnot(all(summary$advance_config == "gain"))
stopifnot(all(families$advance_config == "gain"))
stopifnot(all(c("poly4x_policy_mode", "poly4x_policy_source") %in% names(selections)))
policy_selections <- selections[selections$config %in% c("gain", "diversity", "ocs"), ]
stopifnot(nrow(policy_selections) == 18L)
stopifnot(all(policy_selections$poly4x_policy_mode %in% c("gain", "diversity", "ocs")))
stopifnot(all(c("poly4x_policy_mode", "poly4x_policy_source") %in% names(summary)))
policy_summary <- summary[summary$config %in% c("gain", "diversity", "ocs"), ]
stopifnot(nrow(policy_summary) == 6L)
stopifnot(all(policy_summary$poly4x_policy_mode %in% c("gain", "diversity", "ocs")))
stopifnot(all(c("pred_poly4x_policy_mode", "pred_poly4x_policy_source") %in% names(families)))
policy_families <- families[families$config %in% c("gain", "diversity", "ocs"), ]
stopifnot(nrow(policy_families) == 18L)
stopifnot(all(policy_families$pred_poly4x_policy_mode %in% c("gain", "diversity", "ocs")))
stopifnot(all(c("cache_key", "cache_reused") %in% names(families)))
stopifnot(all(families$cycle %in% c(1L, 2L)))

dupes <- unique(families$cache_key[duplicated(families$cache_key)])
for (key in dupes) {
  rows <- families[families$cache_key == key, ]
  stopifnot(any(rows$cache_reused))
  stopifnot(length(unique(rows$realized_mean_gv)) == 1L)
  stopifnot(length(unique(rows$realized_top10_gv)) == 1L)
  stopifnot(length(unique(rows$realized_max_gv)) == 1L)
}

cat("poly4x controlled runner smoke tests passed\n")
