local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

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

read_cfg <- function(prefix) {
  path <- file.path("results", paste0(prefix, "_config.csv"))
  if (!file.exists(path)) return(data.frame(name = character(), value = character()))
  read.csv(path, stringsAsFactors = FALSE)
}

cfg_value <- function(cfg, name, default = NA_character_) {
  hit <- cfg$value[cfg$name == name]
  if (length(hit)) hit[[1]] else default
}

grid_prefix <- env_chr("NG_GRID_PREFIX", "nextgen_parent_grid")
parent_sizes <- as.integer(strsplit(env_chr("NG_GRID_PARENT_SIZES", "20,30,40,50,60,70,80"), ",", fixed = TRUE)[[1]])
parent_sizes <- parent_sizes[is.finite(parent_sizes)]

comparisons <- list()
selection_summaries <- list()
overalls <- list()
configs <- list()

for (n_parents in parent_sizes) {
  prefix <- paste0(grid_prefix, "_", n_parents, "p")
  cfg <- read_cfg(prefix)
  comp_path <- file.path("results", paste0(prefix, "_comparison_vs_var_simple.csv"))
  sel_path <- file.path("results", paste0(prefix, "_selection_summary.csv"))
  overall_path <- file.path("results", paste0(prefix, "_overall.csv"))
  if (file.exists(comp_path)) {
    comp <- read.csv(comp_path, stringsAsFactors = FALSE)
    if (nrow(comp)) {
      comp$n_parents <- n_parents
      comp$top_crosses <- as.integer(cfg_value(cfg, "top_crosses"))
      comp$effect_training_n <- as.integer(cfg_value(cfg, "effect_training_n"))
      comparisons[[prefix]] <- comp
    }
  }
  if (file.exists(sel_path)) {
    sel <- read.csv(sel_path, stringsAsFactors = FALSE)
    if (nrow(sel)) {
      sel$n_parents <- n_parents
      sel$top_crosses <- as.integer(cfg_value(cfg, "top_crosses"))
      sel$effect_training_n <- as.integer(cfg_value(cfg, "effect_training_n"))
      selection_summaries[[prefix]] <- sel
    }
  }
  if (file.exists(overall_path)) {
    overall <- read.csv(overall_path, stringsAsFactors = FALSE)
    if (nrow(overall)) {
      overall$n_parents <- n_parents
      overall$top_crosses <- as.integer(cfg_value(cfg, "top_crosses"))
      overall$effect_training_n <- as.integer(cfg_value(cfg, "effect_training_n"))
      overalls[[prefix]] <- overall
    }
  }
  if (nrow(cfg)) {
    cfg$scenario <- prefix
    cfg$n_parents <- n_parents
    configs[[prefix]] <- cfg
  }
}

comparison <- bind_rows_fill(comparisons)
selection_summary <- bind_rows_fill(selection_summaries)
overall <- bind_rows_fill(overalls)
config <- bind_rows_fill(configs)

avg_comparison <- if (nrow(comparison)) {
  d <- comparison[comparison$cycle > 0, , drop = FALSE]
  stats::aggregate(
    cbind(delta_mean_gv, delta_top10_gv, delta_max_gv, delta_transgressive_rate) ~
      n_parents + top_crosses + effect_training_n + method,
    d,
    mean
  )
} else data.frame()

final_comparison <- if (nrow(comparison)) {
  do.call(rbind, by(comparison, comparison$n_parents, function(d) {
    d[d$cycle == max(d$cycle), , drop = FALSE]
  }))
} else data.frame()

avg_overall <- if (nrow(overall)) {
  d <- overall[overall$cycle > 0, , drop = FALSE]
  stats::aggregate(
    cbind(mean_gv, top10_gv, max_gv, transgressive_rate) ~
      n_parents + top_crosses + effect_training_n + method,
    d,
    mean
  )
} else data.frame()

winner_rows <- list()
if (nrow(avg_overall)) {
  metrics <- c("mean_gv", "top10_gv", "max_gv", "transgressive_rate")
  for (n in sort(unique(avg_overall$n_parents))) {
    d <- avg_overall[avg_overall$n_parents == n, , drop = FALSE]
    for (metric in metrics) {
      i <- which.max(d[[metric]])
      winner_rows[[paste(n, metric, sep = "_")]] <- data.frame(
        n_parents = n,
        metric = metric,
        method = d$method[i],
        value = d[[metric]][i],
        top_crosses = d$top_crosses[i],
        effect_training_n = d$effect_training_n[i],
        stringsAsFactors = FALSE
      )
    }
  }
}
winner_summary <- bind_rows_fill(winner_rows)

avg_selection <- if (nrow(selection_summary)) {
  selection_summary$lambda_parent_use_mode[is.na(selection_summary$lambda_parent_use_mode)] <- "none"
  selector_cols <- c(
    "meta_selector_confidence",
    "meta_selector_leader_margin",
    "meta_selector_portfolio_score_delta",
    "meta_router_score",
    "meta_router_history_score",
    "meta_router_prior_score",
    "meta_router_plan_score",
    "meta_router_gain_score",
    "meta_router_diversity_score",
    "meta_router_reliability_score",
    "meta_router_history_n",
    "meta_router_history_cycles",
    "meta_router_reliability",
    "meta_router_guard_override",
    "meta_router_guard_score",
    "meta_router_original_guard_score",
    "meta_router_guard_advantage",
    "meta_router_guard_router_penalty",
    "frontier_policy_fallback"
  )
  for (col in selector_cols) {
    if (!(col %in% names(selection_summary))) selection_summary[[col]] <- NA_real_
  }
  stats::aggregate(
    cbind(unique_parents, max_parent_use, parent_use_sq, group_coancestry,
          mean_pair_kinship, lambda_parent_use, lambda_parent_use_input, score_scale,
          meta_selector_confidence, meta_selector_leader_margin,
          meta_selector_portfolio_score_delta,
          meta_router_score, meta_router_history_score, meta_router_prior_score, meta_router_plan_score,
          meta_router_gain_score, meta_router_diversity_score,
          meta_router_reliability_score, meta_router_history_n,
          meta_router_history_cycles, meta_router_reliability,
          meta_router_guard_override, meta_router_guard_score,
          meta_router_original_guard_score, meta_router_guard_advantage,
          meta_router_guard_router_penalty, frontier_policy_fallback) ~
      n_parents + top_crosses + effect_training_n + method + lambda_parent_use_mode,
    selection_summary,
    mean,
    na.rm = TRUE,
    na.action = stats::na.pass
  )
} else data.frame()

dir.create("results", showWarnings = FALSE)
write.csv(comparison, file.path("results", paste0(grid_prefix, "_comparison_all.csv")), row.names = FALSE)
write.csv(avg_comparison, file.path("results", paste0(grid_prefix, "_comparison_avg.csv")), row.names = FALSE)
write.csv(final_comparison, file.path("results", paste0(grid_prefix, "_comparison_final.csv")), row.names = FALSE)
write.csv(overall, file.path("results", paste0(grid_prefix, "_overall_all.csv")), row.names = FALSE)
write.csv(avg_overall, file.path("results", paste0(grid_prefix, "_overall_avg.csv")), row.names = FALSE)
write.csv(winner_summary, file.path("results", paste0(grid_prefix, "_winner_summary.csv")), row.names = FALSE)
write.csv(selection_summary, file.path("results", paste0(grid_prefix, "_selection_all.csv")), row.names = FALSE)
write.csv(avg_selection, file.path("results", paste0(grid_prefix, "_selection_avg.csv")), row.names = FALSE)
write.csv(config, file.path("results", paste0(grid_prefix, "_config_all.csv")), row.names = FALSE)

message("Average comparison by parent size:")
if (nrow(avg_comparison)) print(avg_comparison[order(avg_comparison$n_parents, -avg_comparison$delta_mean_gv), ], row.names = FALSE)
message("Average absolute metrics by parent size:")
if (nrow(avg_overall)) print(avg_overall[order(avg_overall$n_parents, -avg_overall$mean_gv), ], row.names = FALSE)
message("Winner summary:")
if (nrow(winner_summary)) print(winner_summary[order(winner_summary$n_parents, winner_summary$metric), ], row.names = FALSE)
message("Average selection summary:")
if (nrow(avg_selection)) print(avg_selection[order(avg_selection$n_parents, avg_selection$parent_use_sq), ], row.names = FALSE)
