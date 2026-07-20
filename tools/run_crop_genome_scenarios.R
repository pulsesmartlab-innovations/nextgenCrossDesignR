local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

root <- normalizePath(file.path(getwd(), "nextgen_cross_design"), mustWork = FALSE)
if (!dir.exists(root)) root <- normalizePath(file.path(".."), mustWork = TRUE)
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

bind_rows_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    for (m in missing) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  do.call(rbind, x)
}

run_grid <- function(env) {
  script <- file.path(root, "tools", "run_parent_size_grid.ps1")
  status <- ng_with_env_vars(env, {
    system2(
      "powershell",
      c("-ExecutionPolicy", "Bypass", "-File", normalizePath(script, winslash = "\\", mustWork = TRUE))
    )
  })
  if (!identical(status, 0L)) {
    ng_stop("crop genome grid failed with exit code ", status)
  }
  invisible(TRUE)
}

prefix <- env_chr("NG_CROP_GRID_PREFIX", "crop_genome_frontier_grid")
scenario_names <- env_chr("NG_CROP_GRID_SCENARIOS", "compact_selfing,maize_like,large_diploid_approx")
parent_sizes <- env_chr("NG_CROP_GRID_PARENT_SIZES", "20,60,80")
reps <- env_chr("NG_CROP_GRID_REPS", "1")
cycles <- env_chr("NG_CROP_GRID_CYCLES", "2")
use_cpp <- env_chr("NG_CROP_GRID_USE_CPP", "1")
shortlist_multiplier <- env_chr("NG_CROP_GRID_EXTERNAL_SHORTLIST_MULTIPLIER", "10")
methods <- env_chr(
  "NG_CROP_GRID_METHODS",
  paste(c(
    "var_simple_topn",
    "var_simple_ocs10_lps2",
    "popvar_uc_ocs10_lps1",
    "simple_usefa_ocs10_lps1",
    "ng_recomb_gebv_ocs10_lps2",
    "ng_pmv_blend_balanced_ocs10_lps2",
    "ng_meta_portfolio_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2",
    "ng_meta_router_ocs10_lps2",
    "ng_frontier_policy_ocs10_lps2",
    "ng_crop_aware_policy_ocs10_lps2"
  ), collapse = ",")
)

scenarios <- ng_crop_genome_select(scenario_names)
dir.create("results", showWarnings = FALSE)

scenario_prefixes <- character()
for (i in seq_len(nrow(scenarios))) {
  scenario <- scenarios[i, , drop = FALSE]
  scenario_prefix <- paste(prefix, scenario$scenario[[1]], sep = "_")
  scenario_prefixes <- c(scenario_prefixes, scenario_prefix)
  env <- c(
    ng_crop_genome_env(scenario),
    NG_GRID_PREFIX = scenario_prefix,
    NG_GRID_PARENT_SIZES = parent_sizes,
    NG_GRID_REPS = reps,
    NG_GRID_CYCLES = cycles,
    NG_GRID_USE_CPP = use_cpp,
    NG_GRID_EXTERNAL_SHORTLIST_MULTIPLIER = shortlist_multiplier,
    NG_GRID_METHODS = methods,
    NG_ALPHASIMR_THREADS = env_chr("NG_CROP_GRID_ALPHASIMR_THREADS", "1"),
    NG_SHARED_SCORING = env_chr("NG_CROP_GRID_SHARED_SCORING", "1")
  )
  message("Running crop genome scenario: ", scenario$scenario[[1]])
  message("  prefix=", scenario_prefix)
  message("  parent_sizes=", parent_sizes, " reps=", reps, " cycles=", cycles)
  run_grid(env)
}

combine_suffix <- function(suffix) {
  rows <- lapply(seq_along(scenario_prefixes), function(i) {
    path <- file.path("results", paste0(scenario_prefixes[[i]], "_", suffix, ".csv"))
    if (!file.exists(path)) return(NULL)
    d <- read.csv(path, stringsAsFactors = FALSE)
    if (!nrow(d)) return(NULL)
    d$crop_scenario <- scenarios$scenario[[i]]
    d$crop <- scenarios$crop[[i]]
    d$crop_description <- scenarios$description[[i]]
    d$crop_ploidy_label <- scenarios$ploidy_label[[i]]
    d$crop_harness_model <- scenarios$harness_model[[i]]
    d$crop_validation_scope <- scenarios$validation_scope[[i]]
    d$crop_n_chr <- scenarios$n_chr[[i]]
    d$crop_genome_length_m <- scenarios$genome_length_m[[i]]
    d$crop_snp_per_chr <- scenarios$snp_per_chr[[i]]
    d$crop_qtl_per_chr <- scenarios$qtl_per_chr[[i]]
    d$crop_phenotype_h2 <- scenarios$phenotype_h2[[i]]
    d
  })
  out <- bind_rows_fill(rows)
  write.csv(out, file.path("results", paste0(prefix, "_", suffix, ".csv")), row.names = FALSE)
  out
}

overall_avg <- combine_suffix("overall_avg")
comparison_avg <- combine_suffix("comparison_avg")
winner_summary <- combine_suffix("winner_summary")
selection_avg <- combine_suffix("selection_avg")

scenario_config <- scenarios
scenario_config$grid_prefix <- scenario_prefixes
write.csv(scenario_config, file.path("results", paste0(prefix, "_scenario_config.csv")), row.names = FALSE)

message("Crop genome scenario summary:")
if (nrow(winner_summary)) {
  top10 <- winner_summary[winner_summary$metric == "top10_gv", , drop = FALSE]
  print(top10[order(top10$crop_scenario, top10$n_parents), ], row.names = FALSE)
} else if (nrow(overall_avg)) {
  print(overall_avg[order(overall_avg$crop_scenario, overall_avg$n_parents, overall_avg$method), ], row.names = FALSE)
}
message("Wrote combined crop scenario outputs with prefix: results/", prefix, "_*.csv")
