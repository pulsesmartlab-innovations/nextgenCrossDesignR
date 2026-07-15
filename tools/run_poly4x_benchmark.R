# Standalone autotetraploid 4x benchmark runner.

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

poly4x_population_metrics <- function(pop, base_best, rep_id, cycle, method) {
  ng_poly4x_assert_pop4x(pop)
  g <- as.numeric(AlphaSimR::gv(pop)[, 1])
  top_n <- max(1L, ceiling(0.10 * length(g)))
  data.frame(
    rep = rep_id,
    cycle = cycle,
    method = method,
    n_ind = length(g),
    ploidy = as.integer(unique(pop@ploidy)),
    mean_gv = mean(g),
    max_gv = max(g),
    top10_gv = mean(utils::head(sort(g, decreasing = TRUE), top_n)),
    var_gv = if (length(g) > 1L) stats::var(g) else 0,
    transgressive_rate = mean(g > base_best),
    stringsAsFactors = FALSE
  )
}

poly4x_select_method <- function(method, scores, cfg, parent_K) {
  if (grepl("^poly4x_policy_", method)) {
    mode <- sub("^poly4x_policy_", "", method)
    selected <- ng_poly4x_policy(
      scores = scores,
      n_crosses = cfg$top_crosses,
      mode = mode,
      parent_K = parent_K
    )
    selected$poly4x_method <- method
    return(selected)
  }
  if (identical(method, "ng_poly4x_var_topn")) {
    return(ng_poly4x_var_topn(scores, n_crosses = cfg$top_crosses))
  }
  if (identical(method, "ng_poly4x_usefulness_topn")) {
    return(ng_poly4x_usefulness_topn(scores, n_crosses = cfg$top_crosses))
  }
  if (identical(method, "ng_poly4x_ocs")) {
    return(ng_poly4x_ocs(
      scores = scores,
      n_crosses = cfg$top_crosses,
      parent_K = parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      lambda_group = cfg$lambda_group,
      lambda_parent_use = cfg$lambda_parent_use
    ))
  }
  if (
    identical(method, "alphamate_opt") ||
      grepl("^alphamate_opt[0-9]+(?:p[0-9]+)?$", method, perl = TRUE)
  ) {
    target_degree <- cfg$alphamate_target_degree
    if (grepl("^alphamate_opt[0-9]+(?:p[0-9]+)?$", method, perl = TRUE)) {
      value <- sub("^alphamate_opt", "", method)
      value <- gsub("p", ".", value, fixed = TRUE)
      target_degree <- as.numeric(value)
    }
    selected <- ng_select_alphamate(
      scores = scores,
      criterion_col = "poly4x_usefulness",
      n_crosses = cfg$top_crosses,
      parent_K = parent_K,
      executable = cfg$alphamate_executable,
      runtime_path = cfg$alphamate_runtime_path,
      target_degree = target_degree,
      max_contributions = cfg$alphamate_max_contributions,
      number_of_parents = cfg$alphamate_number_of_parents,
      evol_solutions = cfg$alphamate_evol_solutions,
      evol_iterations = cfg$alphamate_evol_iterations,
      evol_stop = cfg$alphamate_evol_stop,
      n_threads = cfg$alphamate_threads,
      selfing_weight = cfg$alphamate_selfing_weight,
      workdir = file.path("results", paste0(cfg$prefix, "_", method, "_alphamate")),
      keep_files = cfg$alphamate_keep_files,
      mode = cfg$alphamate_mode
    )
    selected$poly4x_method <- method
    return(selected)
  }
  ng_stop("Unknown poly4x selection method: ", method)
}

poly4x_realize_selected <- function(parent_pop, selected, cfg, sim_param, rep_id, cycle, method) {
  families <- vector("list", nrow(selected))
  progeny <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    family <- ng_poly4x_make_family(
      parent_pop = parent_pop,
      parent1 = selected$parent1[[i]],
      parent2 = selected$parent2[[i]],
      n_progeny = cfg$realized_progeny_per_cross,
      sim_param = sim_param,
      seed = cfg$seed + rep_id * 100000L + cycle * 1000L + i
    )
    family@id <- paste0(method, "_R", rep_id, "_C", cycle, "_F", i, "_P", seq_len(AlphaSimR::nInd(family)))
    stats <- ng_poly4x_family_stats(family)
    families[[i]] <- data.frame(
      rep = rep_id,
      cycle = cycle,
      method = method,
      family_index = i,
      parent1 = selected$parent1[[i]],
      parent2 = selected$parent2[[i]],
      ploidy = as.integer(unique(family@ploidy)),
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
    progeny[[i]] <- family
  }
  list(pop = merge_pop_list(progeny), families = bind_rows_fill(families))
}

cfg <- list(
  prefix = env_chr("NG_POLY4X_PREFIX", "poly4x_benchmark"),
  scenario = env_chr("NG_POLY4X_SCENARIO", "potato_autotetraploid_4x"),
  methods = env_csv("NG_POLY4X_METHODS", "ng_poly4x_var_topn,ng_poly4x_usefulness_topn,ng_poly4x_ocs"),
  policy_modes = env_csv("NG_POLY4X_POLICY_MODES", ""),
  reps = env_int("NG_POLY4X_REPS", 3L),
  cycles = env_int("NG_POLY4X_CYCLES", 5L),
  n_parents = env_int("NG_POLY4X_N_PARENTS", 80L),
  top_crosses = env_int("NG_POLY4X_TOP_CROSSES", 10L),
  score_progeny_per_cross = env_int("NG_POLY4X_SCORE_PROGENY_PER_CROSS", 25L),
  realized_progeny_per_cross = env_int("NG_POLY4X_REALIZED_PROGENY_PER_CROSS", 40L),
  selection_prop = env_num("NG_POLY4X_SELECTION_PROP", 0.10),
  seed = env_int("NG_POLY4X_SEED", env_int("NG_SEED", 4404L)),
  include_digenic = env_bool("NG_POLY4X_INCLUDE_DIGENIC", TRUE),
  max_crosses_per_parent = env_int("NG_POLY4X_MAX_CROSSES_PER_PARENT", 4L),
  lambda_group = env_num("NG_POLY4X_LAMBDA_GROUP", 0.5),
  lambda_parent_use = env_num("NG_POLY4X_LAMBDA_PARENT_USE", 1.0),
  alphamate_executable = env_chr("NG_POLY4X_ALPHAMATE_EXECUTABLE", ng_alphamate_default_executable()),
  alphamate_runtime_path = env_chr("NG_POLY4X_ALPHAMATE_RUNTIME_PATH", Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")),
  alphamate_target_degree = env_num("NG_POLY4X_ALPHAMATE_TARGET_DEGREE", 45),
  alphamate_max_contributions = env_int("NG_POLY4X_ALPHAMATE_MAX_CONTRIBUTIONS", NA_integer_),
  alphamate_number_of_parents = env_int("NG_POLY4X_ALPHAMATE_NUMBER_OF_PARENTS", NA_integer_),
  alphamate_evol_solutions = env_int("NG_POLY4X_ALPHAMATE_EVOL_SOLUTIONS", 100L),
  alphamate_evol_iterations = env_int("NG_POLY4X_ALPHAMATE_EVOL_ITERATIONS", 1000L),
  alphamate_evol_stop = env_int("NG_POLY4X_ALPHAMATE_EVOL_STOP", 200L),
  alphamate_threads = env_int("NG_POLY4X_ALPHAMATE_THREADS", 1L),
  alphamate_selfing_weight = env_num("NG_POLY4X_ALPHAMATE_SELFING_WEIGHT", -1e6),
  alphamate_keep_files = env_bool("NG_POLY4X_ALPHAMATE_KEEP_FILES", FALSE),
  alphamate_mode = env_chr("NG_POLY4X_ALPHAMATE_MODE", "ModeOptTarget1")
)

if (!is.finite(cfg$alphamate_max_contributions)) cfg$alphamate_max_contributions <- NULL
if (!is.finite(cfg$alphamate_number_of_parents)) cfg$alphamate_number_of_parents <- NULL
if (length(cfg$policy_modes)) {
  cfg$policy_modes <- vapply(cfg$policy_modes, ng_poly4x_normalize_policy_mode, character(1))
  cfg$methods <- unique(c(cfg$methods, paste0("poly4x_policy_", cfg$policy_modes)))
}

scenario <- ng_poly4x_select_scenario(cfg$scenario)
override_fields <- c(
  n_founders = "NG_POLY4X_N_FOUNDERS",
  n_chr = "NG_POLY4X_N_CHR",
  seg_sites = "NG_POLY4X_SEG_SITES",
  snp_per_chr = "NG_POLY4X_SNP_PER_CHR",
  qtl_per_chr = "NG_POLY4X_QTL_PER_CHR"
)
for (field in names(override_fields)) {
  scenario[[field]] <- env_int(override_fields[[field]], scenario[[field]])
}

if (cfg$top_crosses * cfg$realized_progeny_per_cross < cfg$n_parents) {
  ng_stop("top_crosses * realized_progeny_per_cross must be at least n_parents")
}

dir.create("results", showWarnings = FALSE)

rep_results <- vector("list", cfg$reps)
for (rep_id in seq_len(cfg$reps)) {
  message("poly4x replicate ", rep_id, " of ", cfg$reps)
  setup <- ng_poly4x_setup_simparam(
    scenario = scenario,
    seed = cfg$seed + rep_id,
    include_digenic = cfg$include_digenic
  )
  base_parents <- ng_poly4x_make_parent_pop(setup, n_parents = cfg$n_parents)
  base_parents@id <- paste0("R", rep_id, "_C0_P", seq_len(AlphaSimR::nInd(base_parents)))
  base_gv <- as.numeric(AlphaSimR::gv(base_parents)[, 1])
  base_best <- max(base_gv)

  branches <- setNames(vector("list", length(cfg$methods)), cfg$methods)
  for (method in cfg$methods) branches[[method]] <- unserialize(serialize(base_parents, NULL))

  metrics <- list()
  selected_out <- list()
  selection_summary_out <- list()
  family_out <- list()

  for (method in cfg$methods) {
    metrics[[paste(method, 0L, sep = "_")]] <- poly4x_population_metrics(
      branches[[method]], base_best, rep_id, 0L, method
    )
  }

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    for (method in cfg$methods) {
      parent_pop <- branches[[method]]
      parent_pop@id <- paste0(method, "_R", rep_id, "_C", cycle - 1L, "_P", seq_len(AlphaSimR::nInd(parent_pop)))
      scores <- ng_poly4x_score_crosses(
        parent_pop = parent_pop,
        sim_param = setup$sim_param,
        n_score_progeny = cfg$score_progeny_per_cross,
        selection_prop = cfg$selection_prop,
        seed = cfg$seed + rep_id * 10000L + cycle * 100L
      )
      parent_K <- attr(scores, "parent_K")
      selected <- poly4x_select_method(method, scores, cfg, parent_K)
      selected$rep <- rep_id
      selected$cycle <- cycle
      selected$method <- method
      selected$used_dh <- 0L
      selected_out[[paste(rep_id, method, cycle, sep = "_")]] <- selected

      selected_counts <- ng_parent_counts(selected, rownames(parent_K))
      selected_contrib <- selected_counts / sum(selected_counts)
      plan_summary <- attr(selected, "summary")
      selection_summary_out[[paste(rep_id, method, cycle, sep = "_")]] <- data.frame(
        rep = rep_id,
        cycle = cycle,
        method = method,
        n_crosses = nrow(selected),
        unique_parents = sum(selected_counts > 0),
        max_parent_use = max(selected_counts),
        parent_use_sq = sum(selected_contrib * selected_contrib),
        group_coancestry = ng_group_coancestry(selected_counts, parent_K),
        mean_pair_kinship = mean(selected$pair_kinship, na.rm = TRUE),
        mean_pred_poly4x_usefulness = mean(selected$poly4x_usefulness, na.rm = TRUE),
        used_dh = 0L,
        lambda_group = if (!is.null(plan_summary) && "lambda_group" %in% names(plan_summary)) plan_summary$lambda_group else cfg$lambda_group,
        lambda_parent_use = if (!is.null(plan_summary) && "lambda_parent_use" %in% names(plan_summary)) plan_summary$lambda_parent_use else cfg$lambda_parent_use,
        poly4x_policy_mode = if (!is.null(plan_summary) && "poly4x_policy_mode" %in% names(plan_summary)) plan_summary$poly4x_policy_mode else if ("poly4x_policy_mode" %in% names(selected)) selected$poly4x_policy_mode[[1]] else NA_character_,
        poly4x_policy_reason = if (!is.null(plan_summary) && "poly4x_policy_reason" %in% names(plan_summary)) plan_summary$poly4x_policy_reason else if ("poly4x_policy_reason" %in% names(selected)) selected$poly4x_policy_reason[[1]] else NA_character_,
        poly4x_policy_scope = if (!is.null(plan_summary) && "poly4x_policy_scope" %in% names(plan_summary)) plan_summary$poly4x_policy_scope else if ("poly4x_policy_scope" %in% names(selected)) selected$poly4x_policy_scope[[1]] else NA_character_,
        poly4x_policy_source = if (!is.null(plan_summary) && "poly4x_policy_source" %in% names(plan_summary)) plan_summary$poly4x_policy_source else if ("poly4x_policy_source" %in% names(selected)) selected$poly4x_policy_source[[1]] else NA_character_,
        stringsAsFactors = FALSE
      )

      realized <- poly4x_realize_selected(
        parent_pop = parent_pop,
        selected = selected,
        cfg = cfg,
        sim_param = setup$sim_param,
        rep_id = rep_id,
        cycle = cycle,
        method = method
      )
      family_out[[paste(rep_id, method, cycle, sep = "_")]] <- realized$families

      progeny <- realized$pop
      g <- as.numeric(AlphaSimR::gv(progeny)[, 1])
      next_idx <- order(g, decreasing = TRUE)[seq_len(cfg$n_parents)]
      next_pop <- progeny[next_idx]
      next_pop@id <- paste0(method, "_R", rep_id, "_C", cycle, "_P", seq_len(AlphaSimR::nInd(next_pop)))
      branches[[method]] <- next_pop
      metrics[[paste(method, cycle, sep = "_")]] <- poly4x_population_metrics(
        next_pop, base_best, rep_id, cycle, method
      )
    }
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

message("poly4x benchmark complete: wrote ", normalizePath("results", winslash = "/", mustWork = FALSE))
