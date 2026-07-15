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
  message("AlphaSimR unavailable; skipping poly4x grid runner contract test")
  quit(save = "no", status = 0)
}
if (!nzchar(Sys.which("powershell"))) {
  message("powershell unavailable; skipping poly4x grid runner contract test")
  quit(save = "no", status = 0)
}

prefix <- paste0("poly4x_grid_contract_", Sys.getpid())
script <- file.path(root, "tools", "run_poly4x_grid.ps1")
scenario <- "potato_autotetraploid_4x"
parent_size <- "8"
run_prefix <- paste(prefix, scenario, paste0("p", parent_size), sep = "_")

env_values <- c(
  NG_POLY4X_GRID_PREFIX = prefix,
  NG_POLY4X_GRID_SCENARIOS = scenario,
  NG_POLY4X_GRID_PARENT_SIZES = parent_size,
  NG_POLY4X_REPS = "1",
  NG_POLY4X_CYCLES = "1",
  NG_POLY4X_N_FOUNDERS = "14",
  NG_POLY4X_N_CHR = "2",
  NG_POLY4X_SEG_SITES = "24",
  NG_POLY4X_SNP_PER_CHR = "6",
  NG_POLY4X_QTL_PER_CHR = "4",
  NG_POLY4X_TOP_CROSSES = "3",
  NG_POLY4X_SCORE_PROGENY_PER_CROSS = "4",
  NG_POLY4X_REALIZED_PROGENY_PER_CROSS = "5"
)
managed_env <- c(names(env_values), "NG_POLY4X_METHODS", "NG_POLY4X_POLICY_MODES")
old_env <- Sys.getenv(managed_env, unset = NA_character_)
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
on.exit(unlink(file.path("results", paste0(run_prefix, "_*.csv"))), add = TRUE)

Sys.unsetenv(c("NG_POLY4X_METHODS", "NG_POLY4X_POLICY_MODES"))
do.call(Sys.setenv, as.list(env_values))
status <- system2("powershell", c(
  "-ExecutionPolicy", "Bypass",
  "-File", normalizePath(script, winslash = "\\", mustWork = TRUE)
))
stopifnot(identical(status, 0L))

summary_path <- file.path("results", paste0(run_prefix, "_selection_summary.csv"))
selections_path <- file.path("results", paste0(run_prefix, "_selections.csv"))
families_path <- file.path("results", paste0(run_prefix, "_families.csv"))
stopifnot(file.exists(summary_path))
stopifnot(file.exists(selections_path))
stopifnot(file.exists(families_path))

summary <- read.csv(summary_path, stringsAsFactors = FALSE)
selections <- read.csv(selections_path, stringsAsFactors = FALSE)
families <- read.csv(families_path, stringsAsFactors = FALSE)
expected_policy_methods <- paste0("poly4x_policy_", c("gain", "diversity", "ocs"))

stopifnot(all(expected_policy_methods %in% summary$method))
stopifnot(all(expected_policy_methods %in% selections$method))
stopifnot(all(expected_policy_methods %in% families$method))
stopifnot(all(c("poly4x_policy_mode", "poly4x_policy_source") %in% names(selections)))
stopifnot(all(c("pred_poly4x_policy_mode", "pred_poly4x_policy_source") %in% names(families)))

cat("poly4x grid runner contract tests passed\n")
