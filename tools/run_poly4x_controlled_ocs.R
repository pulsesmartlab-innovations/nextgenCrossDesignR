# Standalone controlled autotetraploid 4x OCS comparison runner.

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
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)
ng_poly4x_require_alphasimr()
suppressPackageStartupMessages(library(AlphaSimR))

env_value <- function(name) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value)) return(NA_character_)
  trimws(value)
}

env_int <- function(name, default) {
  value <- env_value(name)
  if (is.na(value) || !nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (
    length(parsed) != 1L ||
      is.na(parsed) ||
      !is.finite(parsed) ||
      parsed != round(parsed) ||
      parsed < -.Machine$integer.max ||
      parsed > .Machine$integer.max
  ) {
    ng_stop("Environment variable ", name, " must be a single finite integer-like value when set")
  }
  as.integer(parsed)
}

env_num <- function(name, default) {
  value <- env_value(name)
  if (is.na(value) || !nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || is.na(parsed) || !is.finite(parsed)) {
    ng_stop("Environment variable ", name, " must be a single finite numeric value when set")
  }
  parsed
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

env_bool <- function(name, default) {
  value <- tolower(env_value(name))
  if (is.na(value) || !nzchar(value)) return(default)
  if (value %in% c("1", "true", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "no", "n")) return(FALSE)
  ng_stop("Environment variable ", name, " must be one of 1/true/yes/y/0/false/no/n when set")
}

env_csv <- function(name, default) {
  value <- env_chr(name, default)
  out <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
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
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

merge_pop_list <- function(pop_list) {
  pop_list <- Filter(Negate(is.null), pop_list)
  if (!length(pop_list)) ng_stop("No poly4x progeny populations to merge")
  Reduce(c, pop_list)
}

poly4x_controlled_population_metrics <- function(pop, base_best, rep_id, cycle, method, config, advance_config) {
  ng_poly4x_assert_pop4x(pop)
  g <- as.numeric(AlphaSimR::gv(pop)[, 1])
  top_n <- max(1L, ceiling(0.10 * length(g)))
  data.frame(
    rep = rep_id,
    cycle = cycle,
    method = method,
    config = config,
    advance_config = advance_config,
    n_ind = length(g),
    ploidy = as.integer(unique(pop@ploidy)),
    used_dh = 0L,
    mean_gv = mean(g),
    max_gv = max(g),
    top10_gv = mean(utils::head(sort(g, decreasing = TRUE), top_n)),
    var_gv = if (length(g) > 1L) stats::var(g) else 0,
    transgressive_rate = mean(g > base_best),
    stringsAsFactors = FALSE
  )
}

poly4x_controlled_select <- function(config, scores, cfg, parent_K) {
  if (config %in% ng_poly4x_policy_modes()) {
    selected <- ng_poly4x_policy(
      scores = scores,
      n_crosses = cfg$top_crosses,
      mode = config,
      parent_K = parent_K
    )
    selected$config <- config
    return(selected)
  }
  if (identical(config, "usefulness")) {
    selected <- ng_poly4x_usefulness_topn(scores, n_crosses = cfg$top_crosses)
    selected$config <- "usefulness"
    return(selected)
  }
  if (identical(config, "default")) {
    selected <- ng_poly4x_ocs(
      scores,
      cfg$top_crosses,
      parent_K,
      max_crosses_per_parent = 4L,
      lambda_group = 0.5,
      lambda_parent_use = 1.0
    )
    selected$config <- "default"
    return(selected)
  }
  if (identical(config, "strongdiv")) {
    selected <- ng_poly4x_ocs(
      scores,
      cfg$top_crosses,
      parent_K,
      max_crosses_per_parent = 3L,
      lambda_group = 1.0,
      lambda_parent_use = 2.0
    )
    selected$config <- "strongdiv"
    return(selected)
  }
  ng_stop("Unknown controlled config: ", config)
}

poly4x_controlled_realize_selected <- function(parent_pop, selected, cfg, sim_param,
                                               rep_id, cycle, method, config, cache) {
  families <- vector("list", nrow(selected))
  progeny <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    parent1 <- selected$parent1[[i]]
    parent2 <- selected$parent2[[i]]
    cache_key <- paste(
      rep_id,
      cycle,
      ng_poly4x_pair_key(parent1, parent2),
      cfg$realized_progeny_per_cross,
      sep = "|"
    )
    cache_reused <- exists(cache_key, envir = cache, inherits = FALSE)
    if (cache_reused) {
      cached <- get(cache_key, envir = cache, inherits = FALSE)
      family <- cached$family
      stats <- cached$stats
    } else {
      family <- ng_poly4x_make_family(
        parent_pop = parent_pop,
        parent1 = parent1,
        parent2 = parent2,
        n_progeny = cfg$realized_progeny_per_cross,
        sim_param = sim_param,
        seed = ng_poly4x_pair_seed(cfg$seed + rep_id * 100000L + cycle * 1000L, parent1, parent2)
      )
      stats <- ng_poly4x_family_stats(family)
      assign(cache_key, list(family = family, stats = stats), envir = cache)
    }
    family_for_branch <- unserialize(serialize(family, NULL))
    family_for_branch@id <- paste0(
      "R", rep_id, "_C", cycle, "_", config, "_F", i, "_P",
      seq_len(AlphaSimR::nInd(family_for_branch))
    )
    families[[i]] <- data.frame(
      rep = rep_id,
      cycle = cycle,
      method = method,
      config = config,
      advance_config = cfg$advance_config,
      family_index = i,
      parent1 = parent1,
      parent2 = parent2,
      ploidy = as.integer(unique(family@ploidy)),
      used_dh = 0L,
      cache_key = cache_key,
      cache_reused = cache_reused,
      pred_poly4x_mean = selected$poly4x_mean[[i]],
      pred_poly4x_var = selected$poly4x_var[[i]],
      pred_poly4x_top10 = selected$poly4x_top10[[i]],
      pred_poly4x_max = selected$poly4x_max[[i]],
      pred_poly4x_usefulness = selected$poly4x_usefulness[[i]],
      pred_poly4x_pair_coancestry = selected$poly4x_pair_coancestry[[i]],
      pred_pair_kinship = selected$pair_kinship[[i]],
      pred_poly4x_score_n = selected$poly4x_score_n[[i]],
      pred_poly4x_policy_mode = if ("poly4x_policy_mode" %in% names(selected)) selected$poly4x_policy_mode[[i]] else NA_character_,
      pred_poly4x_policy_reason = if ("poly4x_policy_reason" %in% names(selected)) selected$poly4x_policy_reason[[i]] else NA_character_,
      pred_poly4x_policy_scope = if ("poly4x_policy_scope" %in% names(selected)) selected$poly4x_policy_scope[[i]] else NA_character_,
      pred_poly4x_policy_source = if ("poly4x_policy_source" %in% names(selected)) selected$poly4x_policy_source[[i]] else NA_character_,
      realized_mean_gv = stats$mean_gv,
      realized_var_gv = stats$var_gv,
      realized_top10_gv = stats$top10_gv,
      realized_max_gv = stats$max_gv,
      stringsAsFactors = FALSE
    )
    progeny[[i]] <- family_for_branch
  }
  list(pop = merge_pop_list(progeny), families = bind_rows_fill(families))
}

cfg <- list(
  prefix = env_chr("NG_POLY4X_CONTROLLED_PREFIX", "poly4x_controlled_ocs"),
  scenario = env_chr("NG_POLY4X_CONTROLLED_SCENARIO", "potato_autotetraploid_4x"),
  configs = env_csv("NG_POLY4X_CONTROLLED_CONFIGS", "usefulness,default,strongdiv"),
  reps = env_int("NG_POLY4X_CONTROLLED_REPS", 3L),
  cycles = env_int("NG_POLY4X_CONTROLLED_CYCLES", 2L),
  n_parents = env_int("NG_POLY4X_CONTROLLED_N_PARENTS", 40L),
  top_crosses = env_int("NG_POLY4X_CONTROLLED_TOP_CROSSES", 6L),
  score_progeny_per_cross = env_int("NG_POLY4X_CONTROLLED_SCORE_PROGENY_PER_CROSS", 20L),
  realized_progeny_per_cross = env_int("NG_POLY4X_CONTROLLED_REALIZED_PROGENY_PER_CROSS", 20L),
  selection_prop = env_num("NG_POLY4X_CONTROLLED_SELECTION_PROP", 0.10),
  advance_config = env_chr("NG_POLY4X_CONTROLLED_ADVANCE_CONFIG", "usefulness"),
  seed = env_int("NG_POLY4X_CONTROLLED_SEED", env_int("NG_SEED", 6605L)),
  include_digenic = env_bool("NG_POLY4X_CONTROLLED_INCLUDE_DIGENIC", TRUE)
)

scenario <- ng_poly4x_select_scenario(cfg$scenario)
override_fields <- c(
  n_founders = "NG_POLY4X_CONTROLLED_N_FOUNDERS",
  n_chr = "NG_POLY4X_CONTROLLED_N_CHR",
  seg_sites = "NG_POLY4X_CONTROLLED_SEG_SITES",
  snp_per_chr = "NG_POLY4X_CONTROLLED_SNP_PER_CHR",
  qtl_per_chr = "NG_POLY4X_CONTROLLED_QTL_PER_CHR"
)
for (field in names(override_fields)) {
  scenario[[field]] <- env_int(override_fields[[field]], scenario[[field]])
}

valid_configs <- c("usefulness", "default", "strongdiv", ng_poly4x_policy_modes())
unknown_configs <- setdiff(cfg$configs, valid_configs)
if (length(unknown_configs)) {
  ng_stop("Unknown controlled config(s): ", paste(unknown_configs, collapse = ", "))
}
duplicate_configs <- unique(cfg$configs[duplicated(cfg$configs)])
if (length(duplicate_configs)) {
  ng_stop("NG_POLY4X_CONTROLLED_CONFIGS must not contain duplicates: ", paste(duplicate_configs, collapse = ", "))
}
if (!cfg$advance_config %in% cfg$configs) {
  ng_stop("NG_POLY4X_CONTROLLED_ADVANCE_CONFIG must be one of NG_POLY4X_CONTROLLED_CONFIGS")
}
if (cfg$top_crosses * cfg$realized_progeny_per_cross < cfg$n_parents) {
  ng_stop("top_crosses * realized_progeny_per_cross must be at least n_parents")
}

dir.create("results", showWarnings = FALSE)

rep_results <- vector("list", cfg$reps)
for (rep_id in seq_len(cfg$reps)) {
  message("poly4x controlled replicate ", rep_id, " of ", cfg$reps)
  setup <- ng_poly4x_setup_simparam(
    scenario = scenario,
    seed = cfg$seed + rep_id,
    include_digenic = cfg$include_digenic
  )
  base_parents <- ng_poly4x_make_parent_pop(setup, n_parents = cfg$n_parents)
  base_parents@id <- paste0("R", rep_id, "_C0_P", seq_len(AlphaSimR::nInd(base_parents)))
  base_gv <- as.numeric(AlphaSimR::gv(base_parents)[, 1])
  base_best <- max(base_gv)
  current_parents <- unserialize(serialize(base_parents, NULL))

  metrics <- list()
  selected_out <- list()
  selection_summary_out <- list()
  family_out <- list()
  family_cache <- new.env(parent = emptyenv())

  for (config in cfg$configs) {
    method <- if (identical(config, "usefulness")) "ng_poly4x_usefulness_topn" else "ng_poly4x_ocs"
    metrics[[paste(config, 0L, sep = "_")]] <- poly4x_controlled_population_metrics(
      current_parents, base_best, rep_id, 0L, method, config, cfg$advance_config
    )
  }

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    current_parents@id <- paste0("R", rep_id, "_C", cycle - 1L, "_P", seq_len(AlphaSimR::nInd(current_parents)))
    scores <- ng_poly4x_score_crosses(
      parent_pop = current_parents,
      sim_param = setup$sim_param,
      n_score_progeny = cfg$score_progeny_per_cross,
      selection_prop = cfg$selection_prop,
      seed = cfg$seed + rep_id * 10000L + cycle * 100L
    )
    parent_K <- attr(scores, "parent_K")
    next_pops <- setNames(vector("list", length(cfg$configs)), cfg$configs)
    for (config in cfg$configs) {
      selected <- poly4x_controlled_select(config, scores, cfg, parent_K)
      method <- unique(selected$poly4x_method)
      if (length(method) != 1L) ng_stop("controlled selection produced ambiguous method labels")
      selected$rep <- rep_id
      selected$cycle <- cycle
      selected$method <- method
      selected$used_dh <- 0L
      selected$advance_config <- cfg$advance_config
      selected_out[[paste(rep_id, config, cycle, sep = "_")]] <- selected

      selected_counts <- ng_parent_counts(selected, rownames(parent_K))
      selected_contrib <- selected_counts / sum(selected_counts)
      plan_summary <- attr(selected, "summary")
      selection_summary_out[[paste(rep_id, config, cycle, sep = "_")]] <- data.frame(
        rep = rep_id,
        cycle = cycle,
        method = method,
        config = config,
        advance_config = cfg$advance_config,
        n_crosses = nrow(selected),
        unique_parents = sum(selected_counts > 0),
        max_parent_use = max(selected_counts),
        parent_use_sq = sum(selected_contrib * selected_contrib),
        group_coancestry = ng_group_coancestry(selected_counts, parent_K),
        mean_pair_kinship = mean(selected$pair_kinship, na.rm = TRUE),
        mean_pred_poly4x_usefulness = mean(selected$poly4x_usefulness, na.rm = TRUE),
        used_dh = 0L,
        lambda_group = if (!is.null(plan_summary) && "lambda_group" %in% names(plan_summary)) plan_summary$lambda_group else NA_real_,
        lambda_parent_use = if (!is.null(plan_summary) && "lambda_parent_use" %in% names(plan_summary)) plan_summary$lambda_parent_use else NA_real_,
        poly4x_policy_mode = if (!is.null(plan_summary) && "poly4x_policy_mode" %in% names(plan_summary)) plan_summary$poly4x_policy_mode else if ("poly4x_policy_mode" %in% names(selected)) selected$poly4x_policy_mode[[1]] else NA_character_,
        poly4x_policy_reason = if (!is.null(plan_summary) && "poly4x_policy_reason" %in% names(plan_summary)) plan_summary$poly4x_policy_reason else if ("poly4x_policy_reason" %in% names(selected)) selected$poly4x_policy_reason[[1]] else NA_character_,
        poly4x_policy_scope = if (!is.null(plan_summary) && "poly4x_policy_scope" %in% names(plan_summary)) plan_summary$poly4x_policy_scope else if ("poly4x_policy_scope" %in% names(selected)) selected$poly4x_policy_scope[[1]] else NA_character_,
        poly4x_policy_source = if (!is.null(plan_summary) && "poly4x_policy_source" %in% names(plan_summary)) plan_summary$poly4x_policy_source else if ("poly4x_policy_source" %in% names(selected)) selected$poly4x_policy_source[[1]] else NA_character_,
        stringsAsFactors = FALSE
      )

      realized <- poly4x_controlled_realize_selected(
        parent_pop = current_parents,
        selected = selected,
        cfg = cfg,
        sim_param = setup$sim_param,
        rep_id = rep_id,
        cycle = cycle,
        method = method,
        config = config,
        cache = family_cache
      )
      family_out[[paste(rep_id, config, cycle, sep = "_")]] <- realized$families

      progeny <- realized$pop
      g <- as.numeric(AlphaSimR::gv(progeny)[, 1])
      next_idx <- order(g, decreasing = TRUE)[seq_len(cfg$n_parents)]
      next_pop <- progeny[next_idx]
      next_pop@id <- paste0("R", rep_id, "_C", cycle, "_P", seq_len(AlphaSimR::nInd(next_pop)))
      next_pops[[config]] <- next_pop
      metrics[[paste(config, cycle, sep = "_")]] <- poly4x_controlled_population_metrics(
        next_pop, base_best, rep_id, cycle, method, config, cfg$advance_config
      )
    }
    current_parents <- next_pops[[cfg$advance_config]]
  }

  rep_results[[rep_id]] <- list(
    metrics = bind_rows_fill(metrics),
    selections = bind_rows_fill(selected_out),
    selection_summary = bind_rows_fill(selection_summary_out),
    families = bind_rows_fill(family_out)
  )
}

metrics <- bind_rows_fill(lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
selection_summary <- bind_rows_fill(lapply(rep_results, `[[`, "selection_summary"))
families <- bind_rows_fill(lapply(rep_results, `[[`, "families"))

prefix <- file.path("results", cfg$prefix)
write.csv(metrics, paste0(prefix, "_metrics.csv"), row.names = FALSE)
write.csv(selections, paste0(prefix, "_selections.csv"), row.names = FALSE)
write.csv(selection_summary, paste0(prefix, "_selection_summary.csv"), row.names = FALSE)
write.csv(families, paste0(prefix, "_families.csv"), row.names = FALSE)

message("poly4x controlled runner complete: wrote ", normalizePath("results", winslash = "/", mustWork = FALSE))
