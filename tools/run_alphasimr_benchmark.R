# Standalone AlphaSimR benchmark for the fresh nextgen_cross_design track.
#
# Default target: 80 DH parents, 5K markers, 3 phenotype reps, 5 cycles.
# Override with NG_* environment variables, for example:
#   NG_REPS=1 NG_CYCLES=1 Rscript tools/run_alphasimr_benchmark.R

local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

root <- normalizePath(file.path(getwd(), "nextgen_cross_design"), mustWork = FALSE)
if (!dir.exists(root)) root <- normalizePath(file.path(".."), mustWork = TRUE)
source(file.path(root, "R", "load.R"))
ng_use_cpp <- tolower(trimws(Sys.getenv("NG_USE_CPP", unset = "0"))) %in% c("1", "true", "yes", "y")
ng_load(root, use_cpp = ng_use_cpp)

suppressPackageStartupMessages(library(AlphaSimR))

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.numeric(value)
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

env_bool <- function(name, default) {
  value <- tolower(trimws(Sys.getenv(name, unset = NA_character_)))
  if (is.na(value) || !nzchar(value)) return(default)
  value %in% c("1", "true", "yes", "y")
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
  do.call(rbind, x)
}

merge_pop_list <- function(pop_list) {
  pop_list <- Filter(Negate(is.null), pop_list)
  if (!length(pop_list)) stop("No populations to merge.", call. = FALSE)
  Reduce(c, pop_list)
}

make_dh_base_parents <- function(founder_pop, n_parents, sim_param) {
  n_dh_per_founder <- ceiling(n_parents / nInd(founder_pop))
  dh <- makeDH(founder_pop, nDH = n_dh_per_founder, keepParents = FALSE, simParam = sim_param)
  dh[seq_len(n_parents)]
}

make_training_population <- function(parent_pop, target_n, sim_param) {
  if (target_n <= nInd(parent_pop)) return(parent_pop)
  extra_n <- target_n - nInd(parent_pop)
  p1 <- sample.int(nInd(parent_pop), extra_n, replace = TRUE)
  p2 <- sample.int(nInd(parent_pop), extra_n, replace = TRUE)
  same <- p1 == p2
  while (any(same)) {
    p2[same] <- sample.int(nInd(parent_pop), sum(same), replace = TRUE)
    same <- p1 == p2
  }
  f1 <- makeCross(parent_pop, cbind(p1, p2), nProgeny = 1, simParam = sim_param)
  aux <- makeDH(f1, nDH = 1, keepParents = FALSE, simParam = sim_param)
  c(parent_pop, aux)
}

with_rng_seed <- function(seed, expr) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  if (!is.null(seed)) set.seed(seed)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

phenotype <- function(pop, h2 = 0.5, reps = 3, seed = NULL) {
  g <- as.numeric(gv(pop)[, 1])
  if (h2 >= 0.999) return(g)
  vg <- stats::var(g)
  ve <- if (is.finite(vg) && vg > 0) vg * (1 - h2) / h2 / max(1, reps) else 1
  with_rng_seed(seed, g + stats::rnorm(length(g), sd = sqrt(ve)))
}

numeric_fingerprint <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  idx <- seq_along(x)
  c(
    length = length(x),
    sum = sum(x),
    weighted = sum(x * ((idx %% 104729L) + 1L)),
    weighted2 = sum((x + 1) * ((idx %% 8191L) + 3L)^2)
  )
}

population_fingerprint <- function(geno, pop) {
  g <- as.numeric(gv(pop)[, 1])
  fp <- c(dim = dim(geno), numeric_fingerprint(geno), gv = numeric_fingerprint(g))
  paste(format(fp, digits = 17, scientific = TRUE), collapse = "|")
}

stable_benchmark_seed <- function(..., seed = 1L) {
  txt <- paste(..., collapse = "|")
  bytes <- utf8ToInt(txt)
  h <- as.numeric(seed %% 2147483647L)
  for (b in bytes) {
    h <- (h * 131 + b) %% 2147483647
  }
  as.integer(max(1, h))
}

population_realization_key <- function(pop) {
  fp <- numeric_fingerprint(as.numeric(gv(pop)[, 1]))
  paste(format(fp, digits = 17, scientific = TRUE), collapse = "|")
}

history_cache_key <- function(method, cfg, calibration_history, global_calibration_history) {
  if (!isTRUE(cfg$shared_scoring)) return(paste0("method:", method))
  if (cfg$calibration_pool == "global") return("global")
  h <- calibration_history[[method]]
  if (is.null(h) || !nrow(h)) return("branch_empty")
  paste0("branch:", method)
}

method_parent_ids <- function(method, rep_id, cycle, n) {
  paste0(method, "_R", rep_id, "_C", cycle - 1L, "_P", seq_len(n))
}

shared_parent_ids <- function(rep_id, cycle, group_id, n) {
  paste0("shared_R", rep_id, "_C", cycle - 1L, "_G", group_id, "_P", seq_len(n))
}

shared_training_ids <- function(rep_id, cycle, group_id, n) {
  paste0("shared_R", rep_id, "_C", cycle - 1L, "_G", group_id, "_T", seq_len(n))
}

remap_score_parent_ids <- function(scores, from_ids, to_ids) {
  id_map <- setNames(as.character(to_ids), as.character(from_ids))
  out <- scores
  out$parent1 <- unname(id_map[as.character(out$parent1)])
  out$parent2 <- unname(id_map[as.character(out$parent2)])
  if (anyNA(out$parent1) || anyNA(out$parent2)) {
    stop("Shared score parent IDs could not be remapped to method IDs.", call. = FALSE)
  }
  out
}

rename_kinship_ids <- function(K, ids) {
  dimnames(K) <- list(as.character(ids), as.character(ids))
  K
}

population_metrics <- function(pop, base_best, cycle, rep_id, method) {
  g <- as.numeric(gv(pop)[, 1])
  data.frame(
    rep = rep_id,
    method = method,
    cycle = cycle,
    n_ind = length(g),
    mean_gv = mean(g),
    max_gv = max(g),
    top10_gv = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
    var_gv = stats::var(g),
    transgressive_rate = mean(g > base_best),
    stringsAsFactors = FALSE
  )
}

method_cap <- function(method, prefix, default) {
  m <- regexec(paste0("^", prefix, "([0-9]+)?"), method)
  hit <- regmatches(method, m)[[1]]
  if (length(hit) >= 2L && nzchar(hit[2])) return(as.integer(hit[2]))
  default
}

method_parent_penalty <- function(method, default_value, default_mode) {
  parse_num <- function(x) as.numeric(gsub("p", ".", x, fixed = TRUE))
  m_abs <- regexec("_lpu([0-9]+(?:p[0-9]+)?)$", method, perl = TRUE)
  hit_abs <- regmatches(method, m_abs)[[1]]
  if (length(hit_abs) >= 2L && nzchar(hit_abs[2])) {
    return(list(value = parse_num(hit_abs[2]), mode = "absolute"))
  }
  m_scale <- regexec("_lps([0-9]+(?:p[0-9]+)?)$", method, perl = TRUE)
  hit_scale <- regmatches(method, m_scale)[[1]]
  if (length(hit_scale) >= 2L && nzchar(hit_scale[2])) {
    return(list(value = parse_num(hit_scale[2]), mode = "adaptive"))
  }
  list(value = default_value, mode = default_mode)
}

apply_family_allocation <- function(selected, cfg, value_col, var_col) {
  method <- match.arg(cfg$family_allocation_method, c("marginal_topk", "score_weighted"))
  ng_allocate_family_sizes(
    selected,
    total_progeny = cfg$total_progeny,
    min_progeny = cfg$min_progeny_per_cross,
    max_progeny = cfg$max_progeny_per_cross,
    value_col = value_col,
    power = cfg$family_size_power,
    method = method,
    mean_col = "cross_mean",
    var_col = var_col,
    selected_top_n = cfg$family_selected_top_n,
    selected_top_prop = cfg$family_selected_top_prop
  )
}

select_top_by_score <- function(scores, score_col, n_crosses) {
  if (!(score_col %in% names(scores))) stop("scores missing ", score_col, call. = FALSE)
  ok <- is.finite(scores[[score_col]])
  if (sum(ok) < n_crosses) {
    stop("Not enough finite scores in ", score_col, " for requested crosses.", call. = FALSE)
  }
  o <- order(scores[[score_col]][ok], decreasing = TRUE)
  scores[which(ok)[o[seq_len(n_crosses)]], , drop = FALSE]
}

first_finite_score_col <- function(scores, cols) {
  cols <- unique(trimws(as.character(cols)))
  cols <- cols[nzchar(cols)]
  for (col in cols) {
    if (col %in% names(scores) && any(is.finite(suppressWarnings(as.numeric(scores[[col]]))))) {
      return(col)
    }
  }
  NA_character_
}

select_by_named_score <- function(method, score_cols, scores, cfg) {
  if (!(method %in% names(score_cols))) return(NULL)
  select_top_by_score(scores, unname(score_cols[[method]]), cfg$top_crosses)
}

select_ocs_by_prefix <- function(method, prefix_score_cols, scores, cfg, parent_K) {
  for (prefix in names(prefix_score_cols)) {
    if (method == prefix || grepl(paste0("^", prefix, "[0-9]*(_lp[us][0-9p]+)?$"), method, perl = TRUE)) {
      return(select_with_ocs(method, prefix, unname(prefix_score_cols[[prefix]]), scores, cfg, parent_K))
    }
  }
  NULL
}

select_with_simplemating <- function(method, prefix, score_col, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$max_crosses_per_parent)
  ng_select_simplemating(
    scores = scores,
    score_col = score_col,
    n_crosses = cfg$top_crosses,
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    min_crosses_per_parent = cfg$simplemating_min_cross,
    max_crosses_to_search = cfg$simplemating_max_search,
    culling_pairwise_k = if (is.na(cfg$simplemating_culling_k)) NULL else cfg$simplemating_culling_k
  )
}

select_with_alphamate <- function(method, scores, cfg, parent_K) {
  criterion_col <- first_finite_score_col(scores, cfg$alphamate_criterion_cols)
  if (!is.character(criterion_col) || !nzchar(criterion_col) || is.na(criterion_col)) {
    stop("No finite AlphaMate criterion column found", call. = FALSE)
  }
  target_degree <- cfg$alphamate_target_degree
  mode <- "ModeOptTarget1"
  if (grepl("^alphamate_opt[0-9]+(?:p[0-9]+)?$", method, perl = TRUE)) {
    value <- sub("^alphamate_opt", "", method)
    value <- gsub("p", ".", value, fixed = TRUE)
    target_degree <- as.numeric(value)
  } else if (identical(method, "alphamate_maxcriterion")) {
    mode <- "ModeMaxCriterion"
  } else if (identical(method, "alphamate_mincoancestry")) {
    mode <- "ModeMinCoancestry"
  }
  number_of_parents <- cfg$alphamate_number_of_parents
  if (!is.finite(number_of_parents)) number_of_parents <- NULL
  max_contributions <- cfg$alphamate_max_contributions
  if (!is.finite(max_contributions)) max_contributions <- cfg$ocs_max_crosses_per_parent
  workdir <- tempfile(
    paste0(gsub("[^A-Za-z0-9_]+", "_", method), "_"),
    tmpdir = normalizePath("results", winslash = "/", mustWork = FALSE)
  )
  ng_select_alphamate(
    scores = scores,
    criterion_col = criterion_col,
    n_crosses = cfg$top_crosses,
    parent_K = parent_K,
    executable = cfg$alphamate_executable,
    runtime_path = cfg$alphamate_runtime_path,
    target_degree = target_degree,
    max_contributions = max_contributions,
    number_of_parents = number_of_parents,
    evol_solutions = cfg$alphamate_evol_solutions,
    evol_iterations = cfg$alphamate_evol_iterations,
    evol_stop = cfg$alphamate_evol_stop,
    n_threads = cfg$alphamate_threads,
    workdir = workdir,
    keep_files = cfg$alphamate_keep_files,
    mode = mode
  )
}

select_with_ocs <- function(method, prefix, gain_col, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$ocs_max_crosses_per_parent)
  pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
  ng_optimize_mating_plan(
    scores = scores,
    n_crosses = cfg$top_crosses,
    gain_col = gain_col,
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    max_pair_kinship = cfg$max_pair_kinship,
    lambda_group = cfg$lambda_group,
    lambda_mating = cfg$lambda_mating,
    lambda_parent_use = pp$value,
    lambda_parent_use_mode = pp$mode,
    method = "mip_contribution",
    local_iter = cfg$local_iter,
    ocs_iter = cfg$ocs_iter
  )
}

select_with_balanced_usefulness <- function(method, prefix, gain_col, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$ocs_max_crosses_per_parent)
  pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
  ng_optimize_balanced_usefulness(
    scores = scores,
    n_crosses = cfg$top_crosses,
    gain_col = gain_col,
    diversity_col = cfg$balanced_diversity_col,
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    min_unique_parents = if (isTRUE(cfg$balanced_auto_min_unique)) "auto" else cfg$balanced_min_unique_parents,
    max_pair_kinship = cfg$max_pair_kinship,
    lambda_group = cfg$lambda_group,
    lambda_mating = cfg$lambda_mating,
    lambda_parent_use = pp$value,
    lambda_parent_use_mode = pp$mode,
    gain_weight = cfg$balanced_gain_weight,
    diversity_weight = cfg$balanced_diversity_weight,
    pair_kinship_weight = cfg$balanced_pair_kinship_weight,
    method = "mip_contribution",
    local_iter = cfg$local_iter,
    ocs_iter = cfg$ocs_iter
  )
}

select_with_adaptive_stack <- function(method, prefix, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$ocs_max_crosses_per_parent)
  pp <- method_parent_penalty(method, NA_real_, cfg$lambda_parent_use_mode)
  stack <- attr(scores, "adaptive_stack")
  if (!is.finite(pp$value)) {
    pp$value <- if (!is.null(stack) && is.finite(stack$parent_use_input)) {
      stack$parent_use_input
    } else {
      cfg$lambda_parent_use
    }
    pp$mode <- "adaptive"
  }
  selected <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = cfg$top_crosses,
    gain_col = "ng_adaptive_score",
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    max_pair_kinship = cfg$max_pair_kinship,
    lambda_group = cfg$lambda_group,
    lambda_mating = cfg$lambda_mating,
    lambda_parent_use = pp$value,
    lambda_parent_use_mode = pp$mode,
    method = "mip_contribution",
    local_iter = cfg$local_iter,
    ocs_iter = cfg$ocs_iter
  )
  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  if (!is.null(stack)) {
    s$adaptive_score_cols <- paste(stack$score_cols, collapse = ",")
    s$adaptive_weights <- paste(names(stack$weights), signif(stack$weights, 4), sep = ":", collapse = ",")
    s$adaptive_reliability <- stack$reliability
    s$adaptive_history_weight <- stack$history_weight
    s$adaptive_history_n <- stack$history_n
    s$adaptive_champion_weight <- stack$champion_weight
    s$adaptive_fallback_col <- stack$fallback_col
    s$adaptive_fallback_weight <- stack$fallback_weight
    s$adaptive_parent_use_input <- stack$parent_use_input
  }
  attr(selected, "summary") <- s
  selected
}

select_with_meta_portfolio <- function(method, prefix, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$ocs_max_crosses_per_parent)
  pp <- method_parent_penalty(method, NA_real_, cfg$lambda_parent_use_mode)
  meta <- attr(scores, "meta_portfolio")
  if (!is.finite(pp$value)) {
    pp$value <- if (!is.null(meta) && is.finite(meta$parent_use_input)) {
      meta$parent_use_input
    } else {
      cfg$lambda_parent_use
    }
    pp$mode <- "adaptive"
  }
  selected <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = cfg$top_crosses,
    gain_col = "ng_meta_score",
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    max_pair_kinship = cfg$max_pair_kinship,
    lambda_group = cfg$lambda_group,
    lambda_mating = cfg$lambda_mating,
    lambda_parent_use = pp$value,
    lambda_parent_use_mode = pp$mode,
    method = "mip_contribution",
    local_iter = cfg$local_iter,
    ocs_iter = cfg$ocs_iter
  )
  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  if (!is.null(meta)) {
    s$meta_score_cols <- paste(meta$score_cols, collapse = ",")
    s$meta_weights <- paste(names(meta$weights), signif(meta$weights, 4), sep = ":", collapse = ",")
    s$meta_reliability <- meta$reliability
    s$meta_history_weight <- meta$history_weight
    s$meta_history_n <- meta$history_n
    s$meta_method_history_weight <- meta$method_history_weight
    s$meta_method_history_n <- meta$method_history_n
    s$meta_method_weights <- paste(names(meta$method_weights), signif(meta$method_weights, 4), sep = ":", collapse = ",")
    s$meta_leader_col <- meta$leader_col
    s$meta_leader_weight <- meta$leader_weight
    s$meta_champion_weight <- meta$champion_weight
    s$meta_parent_use_input <- meta$parent_use_input
  }
  attr(selected, "summary") <- s
  selected
}

meta_selector_confidence <- function(meta, cfg) {
  if (is.null(meta)) {
    return(list(confidence = 0, margin = NA_real_, second_weight = NA_real_))
  }
  leader_weight <- as.numeric(meta$leader_weight)
  if (!length(leader_weight) || !is.finite(leader_weight)) leader_weight <- 0
  method_w <- as.numeric(meta$method_weights)
  method_w <- method_w[is.finite(method_w) & method_w > 0]
  n_scores <- max(2L, length(meta$score_cols))
  second_weight <- if (length(method_w) >= 2L) {
    sort(method_w, decreasing = TRUE)[[2]]
  } else {
    1 / n_scores
  }
  margin <- leader_weight - second_weight
  method_history_weight <- as.numeric(meta$method_history_weight)
  if (!length(method_history_weight) || !is.finite(method_history_weight)) method_history_weight <- 0

  equal_weight <- 1 / n_scores
  weight_component <- max(0, min(1, (leader_weight - equal_weight) / max(1e-6, 1 - equal_weight)))
  margin_component <- max(0, min(1, margin / max(0.05, cfg$meta_selector_margin_scale)))
  history_component <- max(0, min(1, method_history_weight / 0.65))
  confidence <- 0.50 * weight_component + 0.30 * margin_component + 0.20 * history_component
  list(confidence = confidence, margin = margin, second_weight = second_weight)
}

add_meta_selector_columns <- function(selected, decision, gain_col, allocator) {
  selected$ng_meta_selector_decision <- decision
  selected$ng_meta_selector_gain_col <- gain_col
  selected$ng_meta_selector_allocator <- allocator
  selected$ng_meta_selector_score <- if (identical(decision, "portfolio") &&
                                        "ng_meta_score" %in% names(selected)) {
    selected$ng_meta_score
  } else if ("ng_meta_leader_score" %in% names(selected)) {
    selected$ng_meta_leader_score
  } else if ("ng_meta_score" %in% names(selected)) {
    selected$ng_meta_score
  } else {
    NA_real_
  }
  selected$ng_meta_selector_var <- if (identical(decision, "portfolio") &&
                                      "ng_meta_var" %in% names(selected)) {
    selected$ng_meta_var
  } else if ("ng_meta_leader_var" %in% names(selected)) {
    selected$ng_meta_leader_var
  } else if ("ng_meta_var" %in% names(selected)) {
    selected$ng_meta_var
  } else {
    NA_real_
  }
  selected
}

select_with_meta_selector <- function(method, prefix, scores, cfg, parent_K) {
  cap <- method_cap(method, prefix, cfg$ocs_max_crosses_per_parent)
  pp <- method_parent_penalty(method, NA_real_, cfg$lambda_parent_use_mode)
  meta <- attr(scores, "meta_portfolio")
  leader_gain_col <- if (!is.null(meta) && !is.null(meta$leader_col) &&
                         meta$leader_col %in% names(scores)) {
    meta$leader_col
  } else {
    "ng_meta_score"
  }
  if (!is.finite(pp$value)) {
    pp$value <- if (!is.null(meta) && is.finite(meta$parent_use_input)) {
      meta$parent_use_input
    } else {
      cfg$lambda_parent_use
    }
    pp$mode <- "adaptive"
  }
  balanced_leader <- leader_gain_col %in% c("etk_dh_pmv_var_blend_cal")
  leader_selected <- if (isTRUE(balanced_leader)) {
    select_with_balanced_usefulness(
      method = method,
      prefix = prefix,
      gain_col = leader_gain_col,
      scores = scores,
      cfg = cfg,
      parent_K = parent_K
    )
  } else {
    ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = leader_gain_col,
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = pp$value,
      lambda_parent_use_mode = pp$mode,
      method = "mip_contribution",
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    )
  }
  portfolio_selected <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = cfg$top_crosses,
    gain_col = "ng_meta_score",
    parent_K = parent_K,
    max_crosses_per_parent = cap,
    max_pair_kinship = cfg$max_pair_kinship,
    lambda_group = cfg$lambda_group,
    lambda_mating = cfg$lambda_mating,
    lambda_parent_use = pp$value,
    lambda_parent_use_mode = pp$mode,
    method = "mip_contribution",
    local_iter = cfg$local_iter,
    ocs_iter = cfg$ocs_iter
  )
  gate <- meta_selector_confidence(meta, cfg)
  portfolio_score <- if ("ng_meta_score" %in% names(portfolio_selected)) {
    mean(portfolio_selected$ng_meta_score, na.rm = TRUE)
  } else {
    NA_real_
  }
  leader_portfolio_score <- if ("ng_meta_score" %in% names(leader_selected)) {
    mean(leader_selected$ng_meta_score, na.rm = TRUE)
  } else {
    NA_real_
  }
  portfolio_score_delta <- portfolio_score - leader_portfolio_score
  leader_weight <- if (!is.null(meta) && !is.null(meta$leader_weight) &&
                       is.finite(as.numeric(meta$leader_weight))) {
    as.numeric(meta$leader_weight)
  } else {
    0
  }
  leader_strong <- is.finite(gate$confidence) &&
    gate$confidence >= cfg$meta_selector_min_confidence &&
    leader_weight >= cfg$meta_selector_min_leader_weight
  use_portfolio <- is.finite(portfolio_score_delta) &&
    portfolio_score_delta > cfg$meta_selector_portfolio_delta &&
    !isTRUE(leader_strong)
  selected <- if (isTRUE(use_portfolio)) {
    add_meta_selector_columns(portfolio_selected, "portfolio", "ng_meta_score", "mip_contribution")
  } else {
    add_meta_selector_columns(
      leader_selected,
      "leader",
      leader_gain_col,
      if (isTRUE(balanced_leader)) "balanced_usefulness" else "mip_contribution"
    )
  }
  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  s$meta_selector_gain_col <- selected$ng_meta_selector_gain_col[[1]]
  s$meta_selector_allocator <- selected$ng_meta_selector_allocator[[1]]
  s$meta_selector_decision <- selected$ng_meta_selector_decision[[1]]
  s$meta_selector_confidence <- gate$confidence
  s$meta_selector_leader_margin <- gate$margin
  s$meta_selector_portfolio_score_delta <- portfolio_score_delta
  if (!is.null(meta)) {
    s$meta_score_cols <- paste(meta$score_cols, collapse = ",")
    s$meta_weights <- paste(names(meta$weights), signif(meta$weights, 4), sep = ":", collapse = ",")
    s$meta_reliability <- meta$reliability
    s$meta_history_weight <- meta$history_weight
    s$meta_history_n <- meta$history_n
    s$meta_method_history_weight <- meta$method_history_weight
    s$meta_method_history_n <- meta$method_history_n
    s$meta_method_weights <- paste(names(meta$method_weights), signif(meta$method_weights, 4), sep = ":", collapse = ",")
    s$meta_leader_col <- meta$leader_col
    s$meta_leader_weight <- meta$leader_weight
    s$meta_champion_weight <- meta$champion_weight
    s$meta_parent_use_input <- meta$parent_use_input
  }
  attr(selected, "summary") <- s
  selected
}

meta_router_family_patterns <- function() {
  c(
    var_simple = "^var_simple_ocs|^var_simple_select|^var_simple_allocator|^var_simple_mip|^var_simple_repair",
    recomb_gebv = "^ng_recomb_gebv_ocs",
    pmv_balanced = "^ng_pmv_blend_balanced_ocs",
    portfolio = "^ng_meta_portfolio",
    adaptive_stack = "^ng_adaptive_stack",
    popvar_uc = "^popvar_uc_ocs|^popvar_uc_select",
    simple_usefa = "^simple_usefa_ocs|^simple_usefa_select",
    external_uc = "^popvar_uc_ocs|^popvar_uc_select|^simple_usefa_ocs|^simple_usefa_select"
  )
}

meta_router_history_stats <- function(history,
                                      families,
                                      target_col = "realized_top10",
                                      min_n = 10L,
                                      recent_cycles = 3L) {
  families <- unique(trimws(as.character(families)))
  families <- families[nzchar(families)]
  out <- data.frame(
    family = families,
    history_z = NA_real_,
    history_n = 0L,
    history_cycles = 0L,
    stringsAsFactors = FALSE
  )
  if (!length(families) || is.null(history) || !nrow(history) ||
      !(target_col %in% names(history)) || !("method" %in% names(history))) {
    return(out)
  }

  patterns <- meta_router_family_patterns()
  families <- intersect(families, names(patterns))
  if (!length(families)) return(out)

  h <- as.data.frame(history, stringsAsFactors = FALSE)
  cyc <- if ("cycle" %in% names(h)) {
    suppressWarnings(as.integer(h$cycle))
  } else {
    rep(1L, nrow(h))
  }
  if (any(is.finite(cyc))) {
    keep_cycles <- sort(unique(cyc[is.finite(cyc)]), decreasing = TRUE)
    keep_cycles <- keep_cycles[seq_len(min(length(keep_cycles), max(1L, as.integer(recent_cycles))))]
    h <- h[cyc %in% keep_cycles, , drop = FALSE]
    cyc <- cyc[cyc %in% keep_cycles]
  }

  y <- suppressWarnings(as.numeric(h[[target_col]]))
  method_chr <- as.character(h$method)
  routed_family <- if ("pred_ng_meta_router_family" %in% names(h)) {
    as.character(h$pred_ng_meta_router_family)
  } else {
    rep(NA_character_, nrow(h))
  }
  frontier_family <- if ("pred_ng_frontier_policy_family" %in% names(h)) {
    as.character(h$pred_ng_frontier_policy_family)
  } else {
    rep(NA_character_, nrow(h))
  }
  use_frontier_family <- (is.na(routed_family) | !nzchar(routed_family)) &
    !is.na(frontier_family) & nzchar(frontier_family)
  routed_family[use_frontier_family] <- frontier_family[use_frontier_family]
  hit_map <- lapply(families, function(fam) {
    routed_hit <- routed_family == fam
    if (identical(fam, "external_uc")) {
      routed_hit <- routed_family %in% c("popvar_uc", "simple_usefa", "external_uc")
    }
    routed_hit[is.na(routed_hit)] <- FALSE
    hit <- grepl(patterns[[fam]], method_chr, perl = TRUE) | routed_hit
    hit[is.na(hit)] <- FALSE
    hit
  })
  names(hit_map) <- families
  n_by_family <- vapply(hit_map, function(hit) sum(hit & is.finite(y)), integer(1))
  cycle_id <- ifelse(is.finite(cyc), cyc, 1L)
  cycle_levels <- sort(unique(cycle_id[is.finite(cycle_id)]))
  if (!length(cycle_levels)) cycle_levels <- 1L

  cycle_scores <- matrix(NA_real_, nrow = length(cycle_levels), ncol = length(families))
  colnames(cycle_scores) <- families
  rownames(cycle_scores) <- as.character(cycle_levels)
  for (i in seq_along(cycle_levels)) {
    in_cycle <- cycle_id == cycle_levels[[i]]
    for (fam in families) {
      if (n_by_family[[fam]] < min_n) next
      ok <- in_cycle & hit_map[[fam]] & is.finite(y)
      if (any(ok)) cycle_scores[i, fam] <- mean(y[ok], na.rm = TRUE)
    }
    ok_cols <- is.finite(cycle_scores[i, ])
    if (sum(ok_cols) >= 2L && stats::sd(cycle_scores[i, ok_cols]) > 0) {
      cycle_scores[i, ok_cols] <- ng_standardize(cycle_scores[i, ok_cols])
    } else {
      cycle_scores[i, ok_cols] <- NA_real_
    }
  }

  recency_rank <- rank(cycle_levels, ties.method = "first")
  recency_w <- 0.65^(max(recency_rank) - recency_rank)
  recency_w <- recency_w / sum(recency_w)
  for (fam in families) {
    row <- match(fam, out$family)
    if (is.na(row)) next
    ok <- is.finite(cycle_scores[, fam])
    out$history_n[row] <- n_by_family[[fam]]
    out$history_cycles[row] <- sum(ok)
    if (any(ok)) {
      out$history_z[row] <- stats::weighted.mean(cycle_scores[ok, fam], recency_w[ok], na.rm = TRUE)
    }
  }
  out
}

meta_router_first_score_col <- function(scores, cols) {
  cols <- unique(cols[nzchar(cols)])
  for (col in cols) {
    if (col %in% names(scores) && any(is.finite(suppressWarnings(as.numeric(scores[[col]]))))) {
      return(col)
    }
  }
  NA_character_
}

meta_router_variance_col <- function(family, gain_col) {
  if (identical(family, "portfolio")) {
    "ng_meta_var"
  } else if (identical(family, "adaptive_stack")) {
    "ng_adaptive_var"
  } else if (identical(gain_col, "popvar_uc")) {
    "popvar_varG"
  } else if (identical(gain_col, "simple_usefa")) {
    "simple_usefa_var"
  } else {
    ng_variance_col_for_score(gain_col)
  }
}

meta_router_candidate_def <- function(family, scores) {
  family <- tolower(trimws(as.character(family[[1]])))
  gain_col <- switch(
    family,
    var_simple = meta_router_first_score_col(scores, c("etk_var_simple_cal", "uc_var_simple", "var_simple")),
    recomb_gebv = meta_router_first_score_col(scores, c("etk_dh_recomb_var_gebv_cal", "uc_recomb_gebv")),
    pmv_balanced = meta_router_first_score_col(scores, c("etk_dh_pmv_var_blend_cal", "etk_dh_pmv_var_adj_cal")),
    portfolio = meta_router_first_score_col(scores, c("ng_meta_score")),
    adaptive_stack = meta_router_first_score_col(scores, c("ng_adaptive_score")),
    popvar_uc = meta_router_first_score_col(scores, c("popvar_uc")),
    simple_usefa = meta_router_first_score_col(scores, c("simple_usefa")),
    external_uc = meta_router_first_score_col(scores, c("popvar_uc", "simple_usefa")),
    NA_character_
  )
  if (!is.character(gain_col) || !nzchar(gain_col) || is.na(gain_col)) return(NULL)
  allocator <- if (identical(family, "pmv_balanced")) {
    "balanced_usefulness"
  } else {
    "mip_contribution"
  }
  list(
    family = family,
    gain_col = gain_col,
    allocator = allocator,
    var_col = meta_router_variance_col(family, gain_col)
  )
}

meta_router_selected_index <- function(scores, selected) {
  score_key <- paste(as.character(scores$parent1), as.character(scores$parent2), sep = "\r")
  selected_key <- paste(as.character(selected$parent1), as.character(selected$parent2), sep = "\r")
  idx <- match(selected_key, score_key)
  missing <- is.na(idx)
  if (any(missing)) {
    score_key_rev <- paste(as.character(scores$parent2), as.character(scores$parent1), sep = "\r")
    idx[missing] <- match(selected_key[missing], score_key_rev)
  }
  idx
}

meta_router_plan_mean_rank <- function(scores, selected, score_col) {
  if (!(score_col %in% names(scores)) || !(score_col %in% names(selected))) return(NA_real_)
  idx <- meta_router_selected_index(scores, selected)
  z <- ng_rank_normal_score(scores[[score_col]])
  ok <- !is.na(idx) & is.finite(z[idx])
  if (any(ok)) return(mean(z[idx][ok], na.rm = TRUE))
  z_sel <- ng_rank_normal_score(selected[[score_col]])
  if (any(is.finite(z_sel))) mean(z_sel, na.rm = TRUE) else NA_real_
}

meta_router_component_z <- function(x) {
  x <- as.numeric(x)
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (sum(ok) >= 2L && stats::sd(x[ok]) > 0) {
    out[ok] <- ng_standardize(x[ok])
  }
  out
}

meta_router_plan_diversity <- function(summary, n_crosses) {
  if (is.null(summary)) return(NA_real_)
  value <- function(name, default = NA_real_) {
    x <- summary[[name]]
    if (is.null(x) || !length(x) || !is.finite(as.numeric(x[[1]]))) default else as.numeric(x[[1]])
  }
  unique_parents <- value("unique_parents", 0)
  parent_use_sq <- value("parent_use_sq", 1)
  group_coancestry <- value("group_coancestry", 0)
  mean_pair_kinship <- value("mean_pair_kinship", 0)
  unique_scale <- unique_parents / max(1, 2 * as.numeric(n_crosses))
  0.60 * unique_scale - 0.25 * parent_use_sq - 0.10 * group_coancestry - 0.05 * mean_pair_kinship
}

meta_router_reliability_component <- function(family, reliability) {
  rel <- as.numeric(reliability[[1]])
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))
  centered <- rel - 0.50
  if (family %in% c("recomb_gebv", "pmv_balanced", "adaptive_stack")) {
    centered
  } else if (family %in% c("portfolio")) {
    0.50 * centered
  } else if (family %in% c("var_simple", "popvar_uc", "simple_usefa", "external_uc")) {
    -0.35 * centered
  } else {
    0
  }
}

meta_router_prior_component <- function(family, n_parents, n_crosses, reliability) {
  rel <- as.numeric(reliability[[1]])
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))
  n_parents <- as.numeric(n_parents[[1]])
  n_crosses <- as.numeric(n_crosses[[1]])
  if (!is.finite(n_crosses) || n_crosses <= 0) n_crosses <- 10
  if (!is.finite(n_parents) || n_parents <= 0) n_parents <- 2 * n_crosses
  pressure <- max(0, min(1, n_crosses / max(1, n_parents)))

  small_parent <- 1 / (1 + exp((n_parents - 32) / 6))
  mid_parent <- exp(-((n_parents - 45) / 18)^2)
  large_parent <- 1 / (1 + exp(-(n_parents - 62) / 8))
  low_rel <- 1 - rel

  if (identical(family, "var_simple")) {
    0.70 * small_parent + 0.20 * low_rel + 0.20 * pressure - 0.20 * large_parent
  } else if (identical(family, "pmv_balanced")) {
    0.52 * mid_parent + 0.18 * rel + 0.10 * pressure - 0.10 * small_parent
  } else if (identical(family, "recomb_gebv")) {
    0.22 * large_parent + 0.22 * rel + 0.10 * mid_parent - 0.08 * small_parent
  } else if (identical(family, "portfolio")) {
    0.48 * large_parent + 0.16 * mid_parent + 0.12 * rel - 0.05 * small_parent
  } else if (family %in% c("popvar_uc", "simple_usefa", "external_uc")) {
    0.16 * mid_parent + 0.08 * low_rel + 0.06 * pressure
  } else {
    0
  }
}

meta_router_build_candidate <- function(family, method, prefix, scores, cfg, parent_K) {
  def <- meta_router_candidate_def(family, scores)
  if (is.null(def)) return(NULL)
  plan <- tryCatch({
    if (identical(def$family, "pmv_balanced")) {
      select_with_balanced_usefulness(
        method = method,
        prefix = prefix,
        gain_col = def$gain_col,
        scores = scores,
        cfg = cfg,
        parent_K = parent_K
      )
    } else if (identical(def$family, "portfolio")) {
      select_with_meta_portfolio(method, prefix, scores, cfg, parent_K)
    } else if (identical(def$family, "adaptive_stack")) {
      select_with_adaptive_stack(method, prefix, scores, cfg, parent_K)
    } else {
      select_with_ocs(method, prefix, def$gain_col, scores, cfg, parent_K)
    }
  }, error = function(e) {
    attr(e, "router_family") <- family
    e
  })
  if (inherits(plan, "error")) return(NULL)
  list(def = def, plan = plan)
}

meta_router_plan_var <- function(selected, var_col) {
  if (is.character(var_col) && nzchar(var_col) && var_col %in% names(selected)) {
    return(pmax(as.numeric(selected[[var_col]]), 0))
  }
  if ("ng_meta_var" %in% names(selected)) return(pmax(as.numeric(selected$ng_meta_var), 0))
  if ("var_simple_cal" %in% names(selected)) return(pmax(as.numeric(selected$var_simple_cal), 0))
  rep(NA_real_, nrow(selected))
}

add_meta_router_columns <- function(selected, best, best_diag, guard = NULL) {
  if (is.null(guard)) {
    guard <- list(
      override = FALSE,
      selected_family = best$def$family,
      original_family = best$def$family,
      guard_score = NA_real_,
      original_guard_score = NA_real_,
      advantage = NA_real_,
      router_penalty = NA_real_,
      reason = "not_evaluated"
    )
  }
  selected$ng_meta_router_family <- best$def$family
  selected$ng_meta_router_gain_col <- best$def$gain_col
  selected$ng_meta_router_allocator <- best$def$allocator
  selected$ng_meta_router_score <- best_diag$router_score
  selected$ng_meta_router_history_score <- best_diag$history_component
  selected$ng_meta_router_prior_score <- best_diag$prior_z
  selected$ng_meta_router_plan_score <- best_diag$plan_z
  selected$ng_meta_router_gain_score <- best_diag$gain_z
  selected$ng_meta_router_diversity_score <- best_diag$diversity_z
  selected$ng_meta_router_reliability_score <- best_diag$reliability_component
  selected$ng_meta_router_var <- meta_router_plan_var(selected, best$def$var_col)
  selected$ng_meta_router_guard_override <- as.integer(isTRUE(guard$override))
  selected$ng_meta_router_original_family <- guard$original_family
  selected$ng_meta_router_guard_score <- guard$guard_score
  selected$ng_meta_router_original_guard_score <- guard$original_guard_score
  selected$ng_meta_router_guard_advantage <- guard$advantage
  selected$ng_meta_router_guard_router_penalty <- guard$router_penalty
  selected$ng_meta_router_guard_reason <- guard$reason
  selected
}

select_with_meta_router <- function(method, prefix, scores, cfg, parent_K) {
  families <- unique(trimws(as.character(cfg$meta_router_families)))
  families <- families[nzchar(families)]
  if (!length(families)) families <- c("recomb_gebv", "pmv_balanced", "portfolio", "popvar_uc", "simple_usefa", "var_simple")
  history <- attr(scores, "meta_method_history")
  router_min_history_n <- min(
    as.integer(cfg$meta_router_min_history_n),
    max(3L, as.integer(cfg$top_crosses))
  )
  if (!is.finite(router_min_history_n) || router_min_history_n < 1L) {
    router_min_history_n <- 10L
  }
  history_stats <- meta_router_history_stats(
    history = history,
    families = families,
    target_col = cfg$meta_router_target,
    min_n = router_min_history_n,
    recent_cycles = cfg$meta_router_recent_cycles
  )
  candidates <- Filter(Negate(is.null), lapply(families, meta_router_build_candidate,
                                                method = method, prefix = prefix,
                                                scores = scores, cfg = cfg, parent_K = parent_K))
  if (!length(candidates)) {
    return(select_with_meta_selector(method, prefix, scores, cfg, parent_K))
  }

  meta <- attr(scores, "meta_portfolio")
  reliability <- if (!is.null(meta) && is.finite(as.numeric(meta$reliability))) {
    as.numeric(meta$reliability)
  } else if ("effect_reliability" %in% names(scores)) {
    suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  } else {
    NA_real_
  }

  diag <- do.call(rbind, lapply(seq_along(candidates), function(i) {
    cand <- candidates[[i]]
    fam <- cand$def$family
    s <- attr(cand$plan, "summary")
    hist_row <- history_stats[match(fam, history_stats$family), , drop = FALSE]
    hist_z <- if (nrow(hist_row)) hist_row$history_z[[1]] else NA_real_
    hist_n <- if (nrow(hist_row)) hist_row$history_n[[1]] else 0L
    hist_cycles <- if (nrow(hist_row)) hist_row$history_cycles[[1]] else 0L
    data.frame(
      candidate_i = i,
      family = fam,
      gain_col = cand$def$gain_col,
      allocator = cand$def$allocator,
      plan_meta_score = meta_router_plan_mean_rank(scores, cand$plan, "ng_meta_score"),
      plan_gain_score = meta_router_plan_mean_rank(scores, cand$plan, cand$def$gain_col),
      diversity_raw = meta_router_plan_diversity(s, cfg$top_crosses),
      history_z = hist_z,
      history_n = hist_n,
      history_cycles = hist_cycles,
      prior_raw = meta_router_prior_component(
        fam,
        n_parents = if (!is.null(parent_K)) nrow(parent_K) else length(unique(c(scores$parent1, scores$parent2))),
        n_crosses = cfg$top_crosses,
        reliability = reliability
      ),
      reliability_component = meta_router_reliability_component(fam, reliability),
      stringsAsFactors = FALSE
    )
  }))

  diag$prior_z <- meta_router_component_z(diag$prior_raw)
  diag$plan_z <- meta_router_component_z(diag$plan_meta_score)
  diag$gain_z <- meta_router_component_z(diag$plan_gain_score)
  diag$diversity_z <- meta_router_component_z(diag$diversity_raw)
  history_confidence <- diag$history_n / (diag$history_n + max(1, router_min_history_n))
  history_confidence[!is.finite(history_confidence)] <- 0
  diag$history_component <- ifelse(is.finite(diag$history_z), diag$history_z * history_confidence, 0)

  diag$router_score <-
    cfg$meta_router_history_weight * diag$history_component +
    cfg$meta_router_prior_weight * diag$prior_z +
    cfg$meta_router_plan_weight * diag$plan_z +
    cfg$meta_router_gain_weight * diag$gain_z +
    cfg$meta_router_diversity_weight * diag$diversity_z +
    cfg$meta_router_reliability_weight * diag$reliability_component
  diag$router_score[!is.finite(diag$router_score)] <- -Inf
  raw_best_idx <- which.max(diag$router_score)
  guard <- ng_meta_router_regret_guard(
    diag,
    selected_index = raw_best_idx,
    enabled = cfg$meta_router_regret_guard,
    plan_weight = cfg$meta_router_regret_plan_weight,
    gain_weight = cfg$meta_router_regret_gain_weight,
    min_advantage = cfg$meta_router_regret_min_advantage,
    max_router_penalty = cfg$meta_router_regret_max_router_penalty
  )
  best_idx <- guard$selected_index
  best <- candidates[[diag$candidate_i[[best_idx]]]]
  selected <- add_meta_router_columns(best$plan, best, diag[best_idx, , drop = FALSE], guard = guard)

  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  s$meta_router_family <- diag$family[[best_idx]]
  s$meta_router_original_family <- guard$original_family
  s$meta_router_gain_col <- diag$gain_col[[best_idx]]
  s$meta_router_allocator <- diag$allocator[[best_idx]]
  s$meta_router_score <- diag$router_score[[best_idx]]
  s$meta_router_history_score <- diag$history_component[[best_idx]]
  s$meta_router_prior_score <- diag$prior_z[[best_idx]]
  s$meta_router_plan_score <- diag$plan_z[[best_idx]]
  s$meta_router_gain_score <- diag$gain_z[[best_idx]]
  s$meta_router_diversity_score <- diag$diversity_z[[best_idx]]
  s$meta_router_reliability_score <- diag$reliability_component[[best_idx]]
  s$meta_router_history_n <- diag$history_n[[best_idx]]
  s$meta_router_history_min_n <- router_min_history_n
  s$meta_router_history_cycles <- diag$history_cycles[[best_idx]]
  s$meta_router_reliability <- reliability
  s$meta_router_guard_override <- as.integer(isTRUE(guard$override))
  s$meta_router_guard_score <- guard$guard_score
  s$meta_router_original_guard_score <- guard$original_guard_score
  s$meta_router_guard_advantage <- guard$advantage
  s$meta_router_guard_router_penalty <- guard$router_penalty
  s$meta_router_guard_reason <- guard$reason
  s$meta_router_candidates <- paste(diag$family, signif(diag$router_score, 4), sep = ":", collapse = ",")
  s$meta_router_candidate_details <- paste(
    diag$family,
    paste0(
      "hist=", signif(diag$history_component, 3),
      ";prior=", signif(diag$prior_z, 3),
      ";plan=", signif(diag$plan_z, 3),
      ";gain=", signif(diag$gain_z, 3),
      ";div=", signif(diag$diversity_z, 3),
      ";rel=", signif(diag$reliability_component, 3)
    ),
    sep = "[",
    collapse = "],"
  )
  if (nzchar(s$meta_router_candidate_details)) s$meta_router_candidate_details <- paste0(s$meta_router_candidate_details, "]")
  attr(selected, "summary") <- s
  selected
}

frontier_policy_object <- function(cfg) {
  source <- trimws(as.character(cfg$frontier_policy_source))
  if (!nzchar(source)) source <- ng_frontier_default_source()
  spec <- trimws(as.character(cfg$frontier_policy_spec))
  if (nzchar(spec)) {
    ng_frontier_policy_from_spec(spec, source = source)
  } else {
    ng_frontier_default_policy(source = source)
  }
}

add_frontier_policy_columns <- function(selected, policy_choice, selected_method, errors = character()) {
  selected_method <- as.character(selected_method)[[1]]
  selected_family <- ng_frontier_method_family(selected_method)
  fallback <- !identical(selected_method, policy_choice$primary_method)
  error_log <- paste(errors, collapse = " | ")
  reason <- policy_choice$reason
  if (length(errors)) reason <- paste(reason, "candidate_fallback", sep = "+")

  selected$ng_frontier_policy_method <- selected_method
  selected$ng_frontier_policy_family <- selected_family
  selected$ng_frontier_policy_primary_method <- policy_choice$primary_method
  selected$ng_frontier_policy_primary_family <- policy_choice$primary_family
  selected$ng_frontier_policy_source <- policy_choice$source
  selected$ng_frontier_policy_reason <- reason
  selected$ng_frontier_policy_band <- policy_choice$band_label
  selected$ng_frontier_policy_fallback <- as.integer(fallback)
  selected$ng_frontier_policy_candidates <- paste(policy_choice$candidate_methods, collapse = ",")
  selected$ng_frontier_policy_error_log <- error_log

  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  s$frontier_policy_method <- selected_method
  s$frontier_policy_family <- selected_family
  s$frontier_policy_primary_method <- policy_choice$primary_method
  s$frontier_policy_primary_family <- policy_choice$primary_family
  s$frontier_policy_source <- policy_choice$source
  s$frontier_policy_reason <- reason
  s$frontier_policy_band <- policy_choice$band_label
  s$frontier_policy_fallback <- as.integer(fallback)
  s$frontier_policy_candidates <- paste(policy_choice$candidate_methods, collapse = ",")
  s$frontier_policy_error_log <- error_log
  attr(selected, "summary") <- s
  selected
}

select_with_frontier_policy <- function(method, scores, cfg, parent_K) {
  n_parents <- if (!is.null(parent_K)) {
    nrow(parent_K)
  } else {
    length(unique(c(as.character(scores$parent1), as.character(scores$parent2))))
  }
  policy <- frontier_policy_object(cfg)
  fallback_method <- trimws(as.character(cfg$frontier_policy_fallback_method))
  policy_choice <- ng_frontier_policy_select(
    n_parents = n_parents,
    policy = policy,
    fallback_method = fallback_method,
    source = cfg$frontier_policy_source
  )
  candidates <- policy_choice$candidate_methods
  candidates <- candidates[!grepl("^ng_frontier_policy", candidates, perl = TRUE)]
  candidates <- unique(candidates[nzchar(candidates)])
  if (!length(candidates)) ng_stop("frontier policy has no dispatchable candidate methods")

  errors <- character()
  for (candidate in candidates) {
    selected <- tryCatch({
      select_branch(candidate, scores, cfg, parent_K)
    }, error = function(e) {
      e
    })
    if (!inherits(selected, "error")) {
      return(add_frontier_policy_columns(
        selected = selected,
        policy_choice = policy_choice,
        selected_method = candidate,
        errors = errors
      ))
    }
    errors <- c(errors, paste0(candidate, ": ", conditionMessage(selected)))
  }
  ng_stop("frontier policy could not select a candidate: ", paste(errors, collapse = " | "))
}

add_crop_policy_columns <- function(selected, policy_choice, selected_method, errors = character()) {
  selected_method <- as.character(selected_method)[[1]]
  selected_family <- ng_frontier_method_family(selected_method)
  fallback <- !identical(selected_method, policy_choice$primary_method)
  error_log <- paste(errors, collapse = " | ")
  reason <- policy_choice$reason
  if (length(errors)) reason <- paste(reason, "candidate_fallback", sep = "+")

  selected$ng_crop_policy_method <- selected_method
  selected$ng_crop_policy_family <- selected_family
  selected$ng_crop_policy_primary_method <- policy_choice$primary_method
  selected$ng_crop_policy_primary_family <- policy_choice$primary_family
  selected$ng_crop_policy_source <- policy_choice$source
  selected$ng_crop_policy_reason <- reason
  selected$ng_crop_policy_band <- policy_choice$band_label
  selected$ng_crop_policy_fallback <- as.integer(fallback)
  selected$ng_crop_policy_candidates <- paste(policy_choice$candidate_methods, collapse = ",")
  selected$ng_crop_policy_error_log <- error_log
  selected$ng_crop_policy_crop_scenario <- policy_choice$crop_scenario
  selected$ng_crop_policy_crop <- policy_choice$crop
  selected$ng_crop_policy_harness_model <- policy_choice$harness_model
  selected$ng_crop_policy_validation_scope <- policy_choice$validation_scope

  s <- attr(selected, "summary")
  if (is.null(s)) s <- list()
  s$crop_policy_method <- selected_method
  s$crop_policy_family <- selected_family
  s$crop_policy_primary_method <- policy_choice$primary_method
  s$crop_policy_primary_family <- policy_choice$primary_family
  s$crop_policy_source <- policy_choice$source
  s$crop_policy_reason <- reason
  s$crop_policy_band <- policy_choice$band_label
  s$crop_policy_fallback <- as.integer(fallback)
  s$crop_policy_candidates <- paste(policy_choice$candidate_methods, collapse = ",")
  s$crop_policy_error_log <- error_log
  s$crop_policy_crop_scenario <- policy_choice$crop_scenario
  s$crop_policy_crop <- policy_choice$crop
  s$crop_policy_harness_model <- policy_choice$harness_model
  s$crop_policy_validation_scope <- policy_choice$validation_scope
  attr(selected, "summary") <- s
  selected
}

select_with_crop_aware_policy <- function(method, scores, cfg, parent_K) {
  n_parents <- if (!is.null(parent_K)) {
    nrow(parent_K)
  } else {
    length(unique(c(as.character(scores$parent1), as.character(scores$parent2))))
  }
  policy_choice <- ng_crop_aware_policy_select(
    n_parents = n_parents,
    crop_scenario = cfg$crop_scenario,
    crop = cfg$crop,
    harness_model = cfg$crop_harness_model,
    validation_scope = cfg$crop_validation_scope,
    source = cfg$crop_policy_source
  )
  candidates <- policy_choice$candidate_methods
  candidates <- candidates[!grepl("^ng_crop_aware_policy", candidates, perl = TRUE)]
  candidates <- unique(candidates[nzchar(candidates)])
  if (!length(candidates)) ng_stop("crop-aware policy has no dispatchable candidate methods")

  errors <- character()
  for (candidate in candidates) {
    selected <- tryCatch({
      select_branch(candidate, scores, cfg, parent_K)
    }, error = function(e) {
      e
    })
    if (!inherits(selected, "error")) {
      return(add_crop_policy_columns(
        selected = selected,
        policy_choice = policy_choice,
        selected_method = candidate,
        errors = errors
      ))
    }
    errors <- c(errors, paste0(candidate, ": ", conditionMessage(selected)))
  }
  ng_stop("crop-aware policy could not select a candidate: ", paste(errors, collapse = " | "))
}

select_branch <- function(method, scores, cfg, parent_K) {
  i <- ng_selection_intensity(cfg$selection_prop)
  var_mean <- if ("cross_mean_adjusted_pheno" %in% names(scores) &&
                  any(is.finite(scores$cross_mean_adjusted_pheno))) {
    scores$cross_mean_adjusted_pheno
  } else {
    scores$cross_mean
  }
  scores$uc_var_simple <- var_mean + i * sqrt(pmax(scores$var_simple, 0))
  if (method == "var_simple_topn") {
    selected <- scores[order(scores$uc_var_simple, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "var_simple_etk_topn") {
    selected <- scores[order(scores$etk_var_simple_cal, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "popvar_mu_topn") {
    return(select_top_by_score(scores, "popvar_mu", cfg$top_crosses))
  }
  if (method == "popvar_var_topn") {
    return(select_top_by_score(scores, "popvar_varG", cfg$top_crosses))
  }
  if (method == "popvar_uc_topn") {
    return(select_top_by_score(scores, "popvar_uc", cfg$top_crosses))
  }
  if (method == "popvar_musp_topn") {
    return(select_top_by_score(scores, "popvar_musp_high", cfg$top_crosses))
  }
  if (method == "simple_mpv_topn") {
    return(select_top_by_score(scores, "simple_mpv", cfg$top_crosses))
  }
  if (method == "simple_usefa_topn") {
    return(select_top_by_score(scores, "simple_usefa", cfg$top_crosses))
  }
  if (method == "simple_mpv_select" || grepl("^simple_mpv_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "simple_mpv_select", "simple_mpv", scores, cfg, parent_K))
  }
  if (method == "simple_usefa_select" || grepl("^simple_usefa_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "simple_usefa_select", "simple_usefa", scores, cfg, parent_K))
  }
  if (method == "var_simple_select" || grepl("^var_simple_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "var_simple_select", "etk_var_simple_cal", scores, cfg, parent_K))
  }
  if (method == "ng_hybrid_select" || grepl("^ng_hybrid_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "ng_hybrid_select", "etk_dh_pmv_var_blend_cal", scores, cfg, parent_K))
  }
  if (method == "popvar_musp_select" || grepl("^popvar_musp_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "popvar_musp_select", "popvar_musp_high", scores, cfg, parent_K))
  }
  if (method == "popvar_uc_select" || grepl("^popvar_uc_select[0-9]*$", method)) {
    return(select_with_simplemating(method, "popvar_uc_select", "popvar_uc", scores, cfg, parent_K))
  }
  if (method == "alphamate_opt" ||
      grepl("^alphamate_opt[0-9]+(?:p[0-9]+)?$", method, perl = TRUE) ||
      method %in% c("alphamate_maxcriterion", "alphamate_mincoancestry")) {
    return(select_with_alphamate(method, scores, cfg, parent_K))
  }
  if (method == "var_simple_allocator") {
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "uc_var_simple",
      parent_K = parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      method = "greedy_local",
      local_iter = cfg$local_iter
    ))
  }
  if (grepl("^var_simple_repair[0-9]+$", method)) {
    cap <- as.integer(sub("^var_simple_repair", "", method))
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "uc_var_simple",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = 0,
      lambda_mating = cfg$lambda_mating,
      method = "repair_local",
      local_iter = cfg$local_iter
    ))
  }
  if (grepl("^var_simple_mip[0-9]+$", method)) {
    cap <- as.integer(sub("^var_simple_mip", "", method))
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "uc_var_simple",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = 0,
      lambda_mating = cfg$lambda_mating,
      method = "mip_linear",
      local_iter = cfg$local_iter
    ))
  }
  if (method == "var_simple_ocs" || grepl("^var_simple_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    cap <- method_cap(method, "var_simple_ocs", cfg$ocs_max_crosses_per_parent)
    pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_var_simple_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = pp$value,
      lambda_parent_use_mode = pp$mode,
      method = "mip_contribution",
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    ))
  }
  if (method == "var_simple_ocs_fam" || grepl("^var_simple_ocs_fam[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    cap <- method_cap(method, "var_simple_ocs_fam", cfg$ocs_max_crosses_per_parent)
    pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
    selected <- ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_var_simple_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = pp$value,
      lambda_parent_use_mode = pp$mode,
      method = "mip_contribution",
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    )
    return(apply_family_allocation(selected, cfg, "etk_var_simple_cal", "var_simple_cal"))
  }
  if (method == "popvar_musp_ocs" || grepl("^popvar_musp_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "popvar_musp_ocs", "popvar_musp_high", scores, cfg, parent_K))
  }
  if (method == "popvar_uc_ocs" || grepl("^popvar_uc_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "popvar_uc_ocs", "popvar_uc", scores, cfg, parent_K))
  }
  if (method == "popvar_mu_ocs" || grepl("^popvar_mu_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "popvar_mu_ocs", "popvar_mu", scores, cfg, parent_K))
  }
  if (method == "simple_mpv_ocs" || grepl("^simple_mpv_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "simple_mpv_ocs", "simple_mpv", scores, cfg, parent_K))
  }
  if (method == "simple_usefa_ocs" || grepl("^simple_usefa_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "simple_usefa_ocs", "simple_usefa", scores, cfg, parent_K))
  }
  if (method == "ng_portfolio_topn") {
    return(select_top_by_score(scores, "ng_portfolio_score", cfg$top_crosses))
  }
  if (method == "ng_portfolio_ocs" || grepl("^ng_portfolio_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_ocs(method, "ng_portfolio_ocs", "ng_portfolio_score", scores, cfg, parent_K))
  }
  if (method == "ng_adaptive_stack_topn") {
    return(select_top_by_score(scores, "ng_adaptive_score", cfg$top_crosses))
  }
  if (method == "ng_adaptive_stack_ocs" || grepl("^ng_adaptive_stack_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_adaptive_stack(method, "ng_adaptive_stack_ocs", scores, cfg, parent_K))
  }
  if (method == "ng_meta_portfolio_topn") {
    return(select_top_by_score(scores, "ng_meta_score", cfg$top_crosses))
  }
  if (method == "ng_meta_portfolio_ocs" || grepl("^ng_meta_portfolio_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_meta_portfolio(method, "ng_meta_portfolio_ocs", scores, cfg, parent_K))
  }
  if (method == "ng_meta_selector_topn") {
    return(select_top_by_score(scores, "ng_meta_leader_score", cfg$top_crosses))
  }
  if (method == "ng_meta_selector_ocs" || grepl("^ng_meta_selector_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_meta_selector(method, "ng_meta_selector_ocs", scores, cfg, parent_K))
  }
  if (method == "ng_meta_router_ocs" || grepl("^ng_meta_router_ocs[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    return(select_with_meta_router(method, "ng_meta_router_ocs", scores, cfg, parent_K))
  }
  if (method == "ng_frontier_policy" || grepl("^ng_frontier_policy(_ocs[0-9]*(_lp[us][0-9p]+)?)?$", method, perl = TRUE)) {
    return(select_with_frontier_policy(method, scores, cfg, parent_K))
  }
  if (method == "ng_crop_aware_policy" || grepl("^ng_crop_aware_policy(_ocs[0-9]*(_lp[us][0-9p]+)?)?$", method, perl = TRUE)) {
    return(select_with_crop_aware_policy(method, scores, cfg, parent_K))
  }
  balanced_gain_cols <- c(
    ng_useful_balanced_ocs = cfg$balanced_gain_col,
    ng_recomb_gebv_balanced_ocs = "uc_recomb_gebv",
    ng_recomb_blend_balanced_ocs = "etk_dh_recomb_var_blend_cal",
    ng_pmv_blend_balanced_ocs = "etk_dh_pmv_var_blend_cal",
    ng_portfolio_balanced_ocs = "ng_portfolio_score"
  )
  for (prefix in names(balanced_gain_cols)) {
    if (method == prefix || grepl(paste0("^", prefix, "[0-9]*(_lp[us][0-9p]+)?$"), method, perl = TRUE)) {
      return(select_with_balanced_usefulness(
        method = method,
        prefix = prefix,
        gain_col = unname(balanced_gain_cols[[prefix]]),
        scores = scores,
        cfg = cfg,
        parent_K = parent_K
      ))
    }
  }
  mean_source_ocs <- c(
    ng_recomb_ocs = "etk_dh_recomb_var_cal",
    ng_recomb_gebv_ocs = "etk_dh_recomb_var_gebv_cal",
    ng_recomb_adj_ocs = "etk_dh_recomb_var_adj_cal",
    ng_recomb_blend_ocs = "etk_dh_recomb_var_blend_cal",
    ng_pmv_ocs = "etk_dh_pmv_var_cal",
    ng_pmv_gebv_ocs = "etk_dh_pmv_var_gebv_cal",
    ng_pmv_adj_ocs = "etk_dh_pmv_var_adj_cal",
    ng_pmv_blend_ocs = "etk_dh_pmv_var_blend_cal"
  )
  selected <- select_ocs_by_prefix(method, mean_source_ocs, scores, cfg, parent_K)
  if (!is.null(selected)) return(selected)
  mean_source_topn <- c(
    ng_recomb_topn = "uc_recomb",
    ng_recomb_gebv_topn = "uc_recomb_gebv",
    ng_recomb_adj_topn = "uc_recomb_adj",
    ng_recomb_blend_topn = "uc_recomb_blend",
    ng_pmv_topn = "uc_dh",
    ng_pmv_gebv_topn = "uc_dh_gebv",
    ng_pmv_adj_topn = "uc_dh_adj",
    ng_pmv_blend_topn = "uc_dh_blend",
    ng_recomb_cal_topn = "etk_dh_recomb_var_cal",
    ng_recomb_gebv_cal_topn = "etk_dh_recomb_var_gebv_cal",
    ng_recomb_adj_cal_topn = "etk_dh_recomb_var_adj_cal",
    ng_recomb_blend_cal_topn = "etk_dh_recomb_var_blend_cal",
    ng_pmv_cal_topn = "etk_dh_pmv_var_cal",
    ng_pmv_gebv_cal_topn = "etk_dh_pmv_var_gebv_cal",
    ng_pmv_adj_cal_topn = "etk_dh_pmv_var_adj_cal",
    ng_pmv_blend_cal_topn = "etk_dh_pmv_var_blend_cal"
  )
  selected <- select_by_named_score(method, mean_source_topn, scores, cfg)
  if (!is.null(selected)) return(selected)
  if (method == "ng_uc_topn") {
    selected <- scores[order(scores$uc_dh_gebv, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_recomb_topn") {
    selected <- scores[order(scores$uc_recomb, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_recomb_gebv_topn") {
    selected <- scores[order(scores$uc_recomb_gebv, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_pmv_gebv_topn") {
    selected <- scores[order(scores$uc_dh_gebv, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_pmv_topn") {
    selected <- scores[order(scores$uc_dh, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_pmv_cal_topn") {
    selected <- scores[order(scores$etk_dh_pmv_var_cal, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_recomb_cal_topn") {
    selected <- scores[order(scores$etk_dh_recomb_var_cal, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_recomb_gebv_cal_topn") {
    selected <- scores[order(scores$etk_dh_recomb_var_gebv_cal, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_pmv_gebv_cal_topn") {
    selected <- scores[order(scores$etk_dh_pmv_var_gebv_cal, decreasing = TRUE), , drop = FALSE]
    return(utils::head(selected, cfg$top_crosses))
  }
  if (method == "ng_allocator") {
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "uc_dh_gebv",
      parent_K = parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      method = "greedy_local",
      local_iter = cfg$local_iter
    ))
  }
  if (method == "ng_allocator_hybrid") {
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "uc_dh_blend",
      parent_K = parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      method = "greedy_local",
      local_iter = cfg$local_iter
    ))
  }
  if (method == "ng_cal_allocator") {
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_dh_pmv_var_blend_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      method = "greedy_local",
      local_iter = cfg$local_iter
    ))
  }
  if (grepl("^ng_cal_repair[0-9]+$", method)) {
    cap <- as.integer(sub("^ng_cal_repair", "", method))
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_dh_pmv_var_blend_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = 0,
      lambda_mating = cfg$lambda_mating,
      method = "repair_local",
      local_iter = cfg$local_iter
    ))
  }
  if (grepl("^ng_cal_mip[0-9]+$", method)) {
    cap <- as.integer(sub("^ng_cal_mip", "", method))
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_dh_pmv_var_blend_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = 0,
      lambda_mating = cfg$lambda_mating,
      method = "mip_linear",
      local_iter = cfg$local_iter
    ))
  }
  if (method == "ng_ocs_mip" || grepl("^ng_ocs_mip[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    cap <- method_cap(method, "ng_ocs_mip", cfg$ocs_max_crosses_per_parent)
    pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
    return(ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_dh_pmv_var_blend_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = pp$value,
      lambda_parent_use_mode = pp$mode,
      method = "mip_contribution",
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    ))
  }
  if (method == "ng_ocs_fam" || grepl("^ng_ocs_fam[0-9]*(_lp[us][0-9p]+)?$", method, perl = TRUE)) {
    cap <- method_cap(method, "ng_ocs_fam", cfg$ocs_max_crosses_per_parent)
    pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)
    selected <- ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$top_crosses,
      gain_col = "etk_dh_pmv_var_blend_cal",
      parent_K = parent_K,
      max_crosses_per_parent = cap,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = pp$value,
      lambda_parent_use_mode = pp$mode,
      method = "mip_contribution",
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    )
    return(apply_family_allocation(selected, cfg, "etk_dh_pmv_var_blend_cal", "dh_pmv_var_cal"))
  }
  stop("Unknown method: ", method, call. = FALSE)
}

make_selected_progeny <- function(parent_pop,
                                  selected,
                                  progeny_per_cross,
                                  sim_param,
                                  method,
                                  rep_id = 1L,
                                  cycle = 1L,
                                  seed = 1L,
                                  realization_cache = NULL) {
  parent_ids <- parent_pop@id
  pop_key <- population_realization_key(parent_pop)
  progeny <- vector("list", nrow(selected))
  families <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- match(selected$parent1[i], parent_ids)
    p2 <- match(selected$parent2[i], parent_ids)
    if (is.na(p1) || is.na(p2)) stop("Selected parent ID missing from population.", call. = FALSE)
    n_progeny <- if ("n_progeny" %in% names(selected) && is.finite(selected$n_progeny[i])) {
      as.integer(selected$n_progeny[i])
    } else {
      as.integer(progeny_per_cross)
    }
    if (!is.finite(n_progeny) || n_progeny < 1L) stop("Invalid n_progeny for selected cross.", call. = FALSE)
    family_seed <- stable_benchmark_seed(
      seed, rep_id, cycle, pop_key, min(p1, p2), max(p1, p2), n_progeny
    )
    family_key <- paste(pop_key, min(p1, p2), max(p1, p2), n_progeny, family_seed, sep = "|")
    if (!is.null(realization_cache) && exists(family_key, envir = realization_cache, inherits = FALSE)) {
      dh <- get(family_key, envir = realization_cache, inherits = FALSE)
    } else {
      dh <- with_rng_seed(family_seed, {
        f1 <- makeCross(parent_pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
        makeDH(f1, nDH = n_progeny, keepParents = FALSE, simParam = sim_param)
      })
      if (!is.null(realization_cache)) {
        assign(family_key, dh, envir = realization_cache)
      }
    }
    g <- as.numeric(gv(dh)[, 1])
    balanced_summary <- attr(selected, "summary")
    balanced_gain_col <- if (!is.null(balanced_summary) && !is.null(balanced_summary$balanced_gain_col)) {
      balanced_summary$balanced_gain_col
    } else {
      NA_character_
    }
    method_for_pred <- if (grepl("^ng_crop_aware_policy", method, perl = TRUE) &&
                           "ng_crop_policy_method" %in% names(selected)) {
      crop_method <- as.character(selected$ng_crop_policy_method[i])
      if (grepl("^ng_frontier_policy", crop_method, perl = TRUE) &&
          "ng_frontier_policy_method" %in% names(selected)) {
        as.character(selected$ng_frontier_policy_method[i])
      } else {
        crop_method
      }
    } else if (grepl("^ng_frontier_policy", method, perl = TRUE) &&
               "ng_frontier_policy_method" %in% names(selected)) {
      as.character(selected$ng_frontier_policy_method[i])
    } else {
      method
    }
    pred_var <- if (method_for_pred %in% c("var_simple_etk_topn") ||
                    grepl("^var_simple_ocs", method_for_pred) ||
                    grepl("^var_simple_select[0-9]*$", method_for_pred)) {
      selected$var_simple_cal[i]
    } else if (grepl("^popvar_", method_for_pred)) {
      selected$popvar_varG[i]
    } else if (grepl("^simple_usefa", method_for_pred)) {
      selected$simple_usefa_var[i]
    } else if (grepl("^simple_mpv", method_for_pred)) {
      selected$var_simple[i]
    } else if (grepl("^alphamate", method_for_pred)) {
      selected$var_simple[i]
    } else if (grepl("balanced_ocs", method_for_pred)) {
      if (identical(balanced_gain_col, "ng_portfolio_score") && "ng_portfolio_var" %in% names(selected)) {
        selected$ng_portfolio_var[i]
      } else if (identical(balanced_gain_col, "popvar_uc") && "popvar_varG" %in% names(selected)) {
        selected$popvar_varG[i]
      } else if (identical(balanced_gain_col, "simple_usefa") && "simple_usefa_var" %in% names(selected)) {
        selected$simple_usefa_var[i]
      } else if (isTRUE(grepl("pmv", balanced_gain_col, fixed = TRUE))) {
        selected$dh_pmv_var_cal[i]
      } else {
        selected$dh_recomb_var_cal[i]
      }
    } else if (grepl("^ng_portfolio", method_for_pred)) {
      selected$ng_portfolio_var[i]
    } else if (grepl("^ng_adaptive_stack", method_for_pred)) {
      selected$ng_adaptive_var[i]
    } else if (grepl("^ng_meta_portfolio", method_for_pred)) {
      selected$ng_meta_var[i]
    } else if (grepl("^ng_meta_router", method_for_pred)) {
      if ("ng_meta_router_var" %in% names(selected)) {
        selected$ng_meta_router_var[i]
      } else {
        selected$ng_meta_var[i]
      }
    } else if (grepl("^ng_meta_selector", method_for_pred)) {
      if ("ng_meta_selector_var" %in% names(selected)) {
        selected$ng_meta_selector_var[i]
      } else {
        selected$ng_meta_leader_var[i]
      }
    } else if (method_for_pred %in% c("var_simple_topn", "var_simple_allocator") ||
               grepl("^var_simple_repair[0-9]+$", method_for_pred) ||
               grepl("^var_simple_mip[0-9]+$", method_for_pred)) {
      selected$var_simple[i]
    } else if (grepl("^ng_recomb", method_for_pred)) {
      if (grepl("_cal_topn$|_ocs", method_for_pred)) selected$dh_recomb_var_cal[i] else selected$dh_recomb_var[i]
    } else if (grepl("^ng_pmv", method_for_pred)) {
      if (grepl("_cal_topn$|_ocs", method_for_pred)) selected$dh_pmv_var_cal[i] else selected$dh_pmv_var[i]
    } else if (method_for_pred %in% c("ng_cal_allocator") ||
               grepl("^ng_hybrid_select[0-9]*$", method_for_pred) ||
               grepl("^ng_cal_repair[0-9]+$", method_for_pred) ||
               grepl("^ng_cal_mip[0-9]+$", method_for_pred)) {
      selected$dh_pmv_var_cal[i]
    } else {
      selected$dh_pmv_var[i]
    }
    families[[i]] <- data.frame(
      cross_rank = i,
      parent1 = as.character(selected$parent1[i]),
      parent2 = as.character(selected$parent2[i]),
      n_progeny = n_progeny,
      realization_seed = family_seed,
      pred_mean = selected$cross_mean[i],
      pred_mean_gebv = selected$cross_mean_gebv[i],
      pred_mean_adjusted_pheno = selected$cross_mean_adjusted_pheno[i],
      pred_mean_blend = selected$cross_mean_blend[i],
      pred_var = pred_var,
      pred_var_simple = selected$var_simple[i],
      pred_dh_recomb_var = selected$dh_recomb_var[i],
      pred_dh_pmv_var = selected$dh_pmv_var[i],
      pred_var_simple_cal = selected$var_simple_cal[i],
      pred_dh_recomb_var_cal = selected$dh_recomb_var_cal[i],
      pred_dh_pmv_var_cal = selected$dh_pmv_var_cal[i],
      pred_uc_recomb = selected$uc_recomb[i],
      pred_uc_recomb_gebv = selected$uc_recomb_gebv[i],
      pred_uc_recomb_adj = selected$uc_recomb_adj[i],
      pred_uc_recomb_blend = selected$uc_recomb_blend[i],
      pred_uc_dh_gebv = selected$uc_dh_gebv[i],
      pred_uc_dh_adj = selected$uc_dh_adj[i],
      pred_uc_dh_blend = selected$uc_dh_blend[i],
      pred_etk_var_simple_cal = selected$etk_var_simple_cal[i],
      pred_etk_dh_recomb_var_cal = selected$etk_dh_recomb_var_cal[i],
      pred_etk_dh_pmv_var_cal = selected$etk_dh_pmv_var_cal[i],
      pred_etk_dh_pmv_var_blend_cal = selected$etk_dh_pmv_var_blend_cal[i],
      pred_etk_var_simple_gebv_cal = selected$etk_var_simple_gebv_cal[i],
      pred_etk_dh_recomb_var_gebv_cal = selected$etk_dh_recomb_var_gebv_cal[i],
      pred_etk_dh_pmv_var_gebv_cal = selected$etk_dh_pmv_var_gebv_cal[i],
      pred_etk_dh_recomb_var_adj_cal = selected$etk_dh_recomb_var_adj_cal[i],
      pred_etk_dh_pmv_var_adj_cal = selected$etk_dh_pmv_var_adj_cal[i],
      pred_etk_dh_recomb_var_blend_cal = selected$etk_dh_recomb_var_blend_cal[i],
      pred_ng_portfolio_score = if ("ng_portfolio_score" %in% names(selected)) selected$ng_portfolio_score[i] else NA_real_,
      pred_ng_portfolio_var = if ("ng_portfolio_var" %in% names(selected)) selected$ng_portfolio_var[i] else NA_real_,
      pred_ng_portfolio_reliability = if ("ng_portfolio_reliability" %in% names(selected)) selected$ng_portfolio_reliability[i] else NA_real_,
      pred_ng_portfolio_history_weight = if ("ng_portfolio_history_weight" %in% names(selected)) selected$ng_portfolio_history_weight[i] else NA_real_,
      pred_ng_adaptive_score = if ("ng_adaptive_score" %in% names(selected)) selected$ng_adaptive_score[i] else NA_real_,
      pred_ng_adaptive_score_raw = if ("ng_adaptive_score_raw" %in% names(selected)) selected$ng_adaptive_score_raw[i] else NA_real_,
      pred_ng_adaptive_var = if ("ng_adaptive_var" %in% names(selected)) selected$ng_adaptive_var[i] else NA_real_,
      pred_ng_adaptive_reliability = if ("ng_adaptive_reliability" %in% names(selected)) selected$ng_adaptive_reliability[i] else NA_real_,
      pred_ng_adaptive_history_weight = if ("ng_adaptive_history_weight" %in% names(selected)) selected$ng_adaptive_history_weight[i] else NA_real_,
      pred_ng_adaptive_fallback_weight = if ("ng_adaptive_fallback_weight" %in% names(selected)) selected$ng_adaptive_fallback_weight[i] else NA_real_,
      pred_ng_adaptive_fallback_col = if ("ng_adaptive_fallback_col" %in% names(selected)) selected$ng_adaptive_fallback_col[i] else NA_character_,
      pred_ng_adaptive_parent_use_input = if ("ng_adaptive_parent_use_input" %in% names(selected)) selected$ng_adaptive_parent_use_input[i] else NA_real_,
      pred_ng_meta_score = if ("ng_meta_score" %in% names(selected)) selected$ng_meta_score[i] else NA_real_,
      pred_ng_meta_var = if ("ng_meta_var" %in% names(selected)) selected$ng_meta_var[i] else NA_real_,
      pred_ng_meta_reliability = if ("ng_meta_reliability" %in% names(selected)) selected$ng_meta_reliability[i] else NA_real_,
      pred_ng_meta_history_weight = if ("ng_meta_history_weight" %in% names(selected)) selected$ng_meta_history_weight[i] else NA_real_,
      pred_ng_meta_method_history_weight = if ("ng_meta_method_history_weight" %in% names(selected)) selected$ng_meta_method_history_weight[i] else NA_real_,
      pred_ng_meta_leader_score = if ("ng_meta_leader_score" %in% names(selected)) selected$ng_meta_leader_score[i] else NA_real_,
      pred_ng_meta_leader_var = if ("ng_meta_leader_var" %in% names(selected)) selected$ng_meta_leader_var[i] else NA_real_,
      pred_ng_meta_leader_col = if ("ng_meta_leader_col" %in% names(selected)) selected$ng_meta_leader_col[i] else NA_character_,
      pred_ng_meta_selector_score = if ("ng_meta_selector_score" %in% names(selected)) selected$ng_meta_selector_score[i] else NA_real_,
      pred_ng_meta_selector_var = if ("ng_meta_selector_var" %in% names(selected)) selected$ng_meta_selector_var[i] else NA_real_,
      pred_ng_meta_selector_decision = if ("ng_meta_selector_decision" %in% names(selected)) selected$ng_meta_selector_decision[i] else NA_character_,
      pred_ng_meta_selector_gain_col = if ("ng_meta_selector_gain_col" %in% names(selected)) selected$ng_meta_selector_gain_col[i] else NA_character_,
      pred_ng_meta_router_score = if ("ng_meta_router_score" %in% names(selected)) selected$ng_meta_router_score[i] else NA_real_,
      pred_ng_meta_router_var = if ("ng_meta_router_var" %in% names(selected)) selected$ng_meta_router_var[i] else NA_real_,
      pred_ng_meta_router_family = if ("ng_meta_router_family" %in% names(selected)) selected$ng_meta_router_family[i] else NA_character_,
      pred_ng_meta_router_gain_col = if ("ng_meta_router_gain_col" %in% names(selected)) selected$ng_meta_router_gain_col[i] else NA_character_,
      pred_ng_meta_router_allocator = if ("ng_meta_router_allocator" %in% names(selected)) selected$ng_meta_router_allocator[i] else NA_character_,
      pred_ng_meta_router_history_score = if ("ng_meta_router_history_score" %in% names(selected)) selected$ng_meta_router_history_score[i] else NA_real_,
      pred_ng_meta_router_prior_score = if ("ng_meta_router_prior_score" %in% names(selected)) selected$ng_meta_router_prior_score[i] else NA_real_,
      pred_ng_meta_router_plan_score = if ("ng_meta_router_plan_score" %in% names(selected)) selected$ng_meta_router_plan_score[i] else NA_real_,
      pred_ng_meta_router_gain_score = if ("ng_meta_router_gain_score" %in% names(selected)) selected$ng_meta_router_gain_score[i] else NA_real_,
      pred_ng_meta_router_diversity_score = if ("ng_meta_router_diversity_score" %in% names(selected)) selected$ng_meta_router_diversity_score[i] else NA_real_,
      pred_ng_meta_router_reliability_score = if ("ng_meta_router_reliability_score" %in% names(selected)) selected$ng_meta_router_reliability_score[i] else NA_real_,
      pred_ng_meta_router_guard_override = if ("ng_meta_router_guard_override" %in% names(selected)) selected$ng_meta_router_guard_override[i] else NA_integer_,
      pred_ng_meta_router_guard_score = if ("ng_meta_router_guard_score" %in% names(selected)) selected$ng_meta_router_guard_score[i] else NA_real_,
      pred_ng_meta_router_original_guard_score = if ("ng_meta_router_original_guard_score" %in% names(selected)) selected$ng_meta_router_original_guard_score[i] else NA_real_,
      pred_ng_meta_router_guard_advantage = if ("ng_meta_router_guard_advantage" %in% names(selected)) selected$ng_meta_router_guard_advantage[i] else NA_real_,
      pred_ng_meta_router_guard_router_penalty = if ("ng_meta_router_guard_router_penalty" %in% names(selected)) selected$ng_meta_router_guard_router_penalty[i] else NA_real_,
      pred_ng_meta_router_original_family = if ("ng_meta_router_original_family" %in% names(selected)) selected$ng_meta_router_original_family[i] else NA_character_,
      pred_ng_meta_router_guard_reason = if ("ng_meta_router_guard_reason" %in% names(selected)) selected$ng_meta_router_guard_reason[i] else NA_character_,
      pred_ng_frontier_policy_method = if ("ng_frontier_policy_method" %in% names(selected)) selected$ng_frontier_policy_method[i] else NA_character_,
      pred_ng_frontier_policy_family = if ("ng_frontier_policy_family" %in% names(selected)) selected$ng_frontier_policy_family[i] else NA_character_,
      pred_ng_frontier_policy_primary_method = if ("ng_frontier_policy_primary_method" %in% names(selected)) selected$ng_frontier_policy_primary_method[i] else NA_character_,
      pred_ng_frontier_policy_primary_family = if ("ng_frontier_policy_primary_family" %in% names(selected)) selected$ng_frontier_policy_primary_family[i] else NA_character_,
      pred_ng_frontier_policy_source = if ("ng_frontier_policy_source" %in% names(selected)) selected$ng_frontier_policy_source[i] else NA_character_,
      pred_ng_frontier_policy_reason = if ("ng_frontier_policy_reason" %in% names(selected)) selected$ng_frontier_policy_reason[i] else NA_character_,
      pred_ng_frontier_policy_band = if ("ng_frontier_policy_band" %in% names(selected)) selected$ng_frontier_policy_band[i] else NA_character_,
      pred_ng_frontier_policy_fallback = if ("ng_frontier_policy_fallback" %in% names(selected)) selected$ng_frontier_policy_fallback[i] else NA_integer_,
      pred_ng_frontier_policy_candidates = if ("ng_frontier_policy_candidates" %in% names(selected)) selected$ng_frontier_policy_candidates[i] else NA_character_,
      pred_ng_frontier_policy_error_log = if ("ng_frontier_policy_error_log" %in% names(selected)) selected$ng_frontier_policy_error_log[i] else NA_character_,
      pred_ng_crop_policy_method = if ("ng_crop_policy_method" %in% names(selected)) selected$ng_crop_policy_method[i] else NA_character_,
      pred_ng_crop_policy_family = if ("ng_crop_policy_family" %in% names(selected)) selected$ng_crop_policy_family[i] else NA_character_,
      pred_ng_crop_policy_primary_method = if ("ng_crop_policy_primary_method" %in% names(selected)) selected$ng_crop_policy_primary_method[i] else NA_character_,
      pred_ng_crop_policy_primary_family = if ("ng_crop_policy_primary_family" %in% names(selected)) selected$ng_crop_policy_primary_family[i] else NA_character_,
      pred_ng_crop_policy_source = if ("ng_crop_policy_source" %in% names(selected)) selected$ng_crop_policy_source[i] else NA_character_,
      pred_ng_crop_policy_reason = if ("ng_crop_policy_reason" %in% names(selected)) selected$ng_crop_policy_reason[i] else NA_character_,
      pred_ng_crop_policy_band = if ("ng_crop_policy_band" %in% names(selected)) selected$ng_crop_policy_band[i] else NA_character_,
      pred_ng_crop_policy_fallback = if ("ng_crop_policy_fallback" %in% names(selected)) selected$ng_crop_policy_fallback[i] else NA_integer_,
      pred_ng_crop_policy_candidates = if ("ng_crop_policy_candidates" %in% names(selected)) selected$ng_crop_policy_candidates[i] else NA_character_,
      pred_ng_crop_policy_error_log = if ("ng_crop_policy_error_log" %in% names(selected)) selected$ng_crop_policy_error_log[i] else NA_character_,
      pred_ng_crop_policy_crop_scenario = if ("ng_crop_policy_crop_scenario" %in% names(selected)) selected$ng_crop_policy_crop_scenario[i] else NA_character_,
      pred_ng_crop_policy_crop = if ("ng_crop_policy_crop" %in% names(selected)) selected$ng_crop_policy_crop[i] else NA_character_,
      pred_ng_crop_policy_harness_model = if ("ng_crop_policy_harness_model" %in% names(selected)) selected$ng_crop_policy_harness_model[i] else NA_character_,
      pred_ng_crop_policy_validation_scope = if ("ng_crop_policy_validation_scope" %in% names(selected)) selected$ng_crop_policy_validation_scope[i] else NA_character_,
      pred_ng_meta_champion_weight = if ("ng_meta_champion_weight" %in% names(selected)) selected$ng_meta_champion_weight[i] else NA_real_,
      pred_ng_meta_parent_use_input = if ("ng_meta_parent_use_input" %in% names(selected)) selected$ng_meta_parent_use_input[i] else NA_real_,
      pred_balanced_usefulness_score = if (".balanced_usefulness_gain" %in% names(selected)) selected$.balanced_usefulness_gain[i] else NA_real_,
      pred_balanced_gain_col = balanced_gain_col,
      pred_balanced_diversity_col = if (!is.null(balanced_summary) && !is.null(balanced_summary$balanced_diversity_col)) balanced_summary$balanced_diversity_col else NA_character_,
      pred_popvar_varG = if ("popvar_varG" %in% names(selected)) selected$popvar_varG[i] else NA_real_,
      pred_popvar_mu = if ("popvar_mu" %in% names(selected)) selected$popvar_mu[i] else NA_real_,
      pred_popvar_uc = if ("popvar_uc" %in% names(selected)) selected$popvar_uc[i] else NA_real_,
      pred_popvar_musp_high = if ("popvar_musp_high" %in% names(selected)) selected$popvar_musp_high[i] else NA_real_,
      pred_simple_usefa_var = if ("simple_usefa_var" %in% names(selected)) selected$simple_usefa_var[i] else NA_real_,
      pred_simple_usefa = if ("simple_usefa" %in% names(selected)) selected$simple_usefa[i] else NA_real_,
      pair_kinship = selected$pair_kinship[i],
      realized_mean = mean(g),
      realized_var = stats::var(g),
      realized_top10 = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
      realized_max = max(g),
      stringsAsFactors = FALSE
    )
    progeny[[i]] <- dh
  }
  list(pop = merge_pop_list(progeny), families = bind_rows_fill(families))
}

cfg <- list(
  seed = env_int("NG_SEED", 9001L),
  reps = env_int("NG_REPS", 3L),
  cycles = env_int("NG_CYCLES", 5L),
  n_founders = env_int("NG_N_FOUNDERS", 120L),
  n_parents = env_int("NG_N_PARENTS", 80L),
  n_chr = env_int("NG_N_CHR", 5L),
  genome_length_m = env_num("NG_GENOME_LENGTH_M", 1),
  seg_sites = env_int("NG_SEG_SITES", 1200L),
  snp_per_chr = env_int("NG_SNP_PER_CHR", 1000L),
  qtl_per_chr = env_int("NG_QTL_PER_CHR", 40L),
  top_crosses = env_int("NG_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("NG_PROGENY_PER_CROSS", 40L),
  phenotype_h2 = env_num("NG_PHENO_H2", 0.5),
  phenotype_reps = env_int("NG_PHENO_REPS", 3L),
  effect_training_n = env_int("NG_EFFECT_TRAINING_N", 80L),
  effect_h2_prior = env_num("NG_EFFECT_H2_PRIOR", 0.5),
  effect_kfold = env_int("NG_EFFECT_KFOLD", 0L),
  selection_prop = env_num("NG_SELECTION_PROP", 0.10),
  min_effect_reliability = env_num("NG_MIN_EFFECT_RELIABILITY", 0.35),
  max_crosses_per_parent = env_int("NG_MAX_CROSSES_PER_PARENT", 4L),
  max_pair_kinship = env_num("NG_MAX_PAIR_KINSHIP", Inf),
  lambda_group = env_num("NG_LAMBDA_GROUP", 1.0),
  lambda_mating = env_num("NG_LAMBDA_MATING", 0.0),
  lambda_parent_use = env_num("NG_LAMBDA_PARENT_USE", 2.0),
  lambda_parent_use_mode = env_chr("NG_LAMBDA_PARENT_USE_MODE", "adaptive"),
  ocs_iter = env_int("NG_OCS_ITER", 5L),
  ocs_max_crosses_per_parent = env_int("NG_OCS_MAX_CROSSES_PER_PARENT", 10L),
  balanced_gain_col = env_chr("NG_BALANCED_GAIN_COL", "uc_recomb_gebv"),
  balanced_diversity_col = env_chr("NG_BALANCED_DIVERSITY_COL", "var_simple_cal"),
  balanced_gain_weight = env_num("NG_BALANCED_GAIN_WEIGHT", 1.0),
  balanced_diversity_weight = env_num("NG_BALANCED_DIVERSITY_WEIGHT", 0.30),
  balanced_pair_kinship_weight = env_num("NG_BALANCED_PAIR_KINSHIP_WEIGHT", 0.10),
  balanced_min_unique_parents = env_int("NG_BALANCED_MIN_UNIQUE_PARENTS", NA_integer_),
  balanced_auto_min_unique = env_bool("NG_BALANCED_AUTO_MIN_UNIQUE", FALSE),
  total_progeny = env_int("NG_TOTAL_PROGENY", NA_integer_),
  min_progeny_per_cross = env_int("NG_MIN_PROGENY_PER_CROSS", 10L),
  max_progeny_per_cross = env_int("NG_MAX_PROGENY_PER_CROSS", 80L),
  family_size_power = env_num("NG_FAMILY_SIZE_POWER", 1.0),
  family_allocation_method = env_chr("NG_FAMILY_ALLOCATION_METHOD", "marginal_topk"),
  family_selected_top_n = env_int("NG_FAMILY_SELECTED_TOP_N", NA_integer_),
  family_selected_top_prop = env_num("NG_FAMILY_SELECTED_TOP_PROP", 0.10),
  local_iter = env_int("NG_LOCAL_ITER", 1000L),
  use_cpp = ng_use_cpp,
  alphasimr_threads = env_int("NG_ALPHASIMR_THREADS", 1L),
  calibration_min_n = env_int("NG_CALIBRATION_MIN_N", 20L),
  calibration_pool = env_chr("NG_CALIBRATION_POOL", "branch"),
  portfolio_target = env_chr("NG_PORTFOLIO_TARGET", "realized_top10"),
  portfolio_min_history_n = env_int("NG_PORTFOLIO_MIN_HISTORY_N", 20L),
  adaptive_score_cols = env_csv(
    "NG_ADAPTIVE_SCORE_COLS",
    "uc_recomb_gebv,etk_dh_recomb_var_gebv_cal,etk_dh_recomb_var_blend_cal,etk_dh_pmv_var_blend_cal,etk_var_simple_cal"
  ),
  adaptive_target = env_chr("NG_ADAPTIVE_TARGET", "realized_top10"),
  adaptive_min_history_n = env_int("NG_ADAPTIVE_MIN_HISTORY_N", 20L),
  adaptive_history_weight = env_num("NG_ADAPTIVE_HISTORY_WEIGHT", NA_real_),
  adaptive_fallback_col = env_chr("NG_ADAPTIVE_FALLBACK_COL", "etk_var_simple_cal"),
  adaptive_fallback_max_weight = env_num("NG_ADAPTIVE_FALLBACK_MAX_WEIGHT", 0.35),
  meta_score_cols = env_csv(
    "NG_META_SCORE_COLS",
    "ng_adaptive_score,uc_recomb_gebv,etk_dh_pmv_var_blend_cal,etk_var_simple_cal,popvar_uc,simple_usefa"
  ),
  meta_target = env_chr("NG_META_TARGET", "realized_top10"),
  meta_min_history_n = env_int("NG_META_MIN_HISTORY_N", 20L),
  meta_method_min_history_n = env_int("NG_META_METHOD_MIN_HISTORY_N", 10L),
  meta_method_recent_cycles = env_int("NG_META_METHOD_RECENT_CYCLES", 3L),
  meta_history_weight = env_num("NG_META_HISTORY_WEIGHT", NA_real_),
  meta_selector_min_confidence = env_num("NG_META_SELECTOR_MIN_CONFIDENCE", 0.45),
  meta_selector_min_leader_weight = env_num("NG_META_SELECTOR_MIN_LEADER_WEIGHT", 0.45),
  meta_selector_margin_scale = env_num("NG_META_SELECTOR_MARGIN_SCALE", 0.25),
  meta_selector_portfolio_delta = env_num("NG_META_SELECTOR_PORTFOLIO_DELTA", 0.05),
  meta_router_families = env_csv(
    "NG_META_ROUTER_FAMILIES",
    "recomb_gebv,pmv_balanced,portfolio,popvar_uc,simple_usefa,var_simple"
  ),
  meta_router_target = env_chr("NG_META_ROUTER_TARGET", "realized_top10"),
  meta_router_min_history_n = env_int("NG_META_ROUTER_MIN_HISTORY_N", 10L),
  meta_router_recent_cycles = env_int("NG_META_ROUTER_RECENT_CYCLES", 3L),
  meta_router_history_weight = env_num("NG_META_ROUTER_HISTORY_WEIGHT", 0.50),
  meta_router_prior_weight = env_num("NG_META_ROUTER_PRIOR_WEIGHT", 0.25),
  meta_router_plan_weight = env_num("NG_META_ROUTER_PLAN_WEIGHT", 0.05),
  meta_router_gain_weight = env_num("NG_META_ROUTER_GAIN_WEIGHT", 0.10),
  meta_router_diversity_weight = env_num("NG_META_ROUTER_DIVERSITY_WEIGHT", 0.12),
  meta_router_reliability_weight = env_num("NG_META_ROUTER_RELIABILITY_WEIGHT", 0.03),
  meta_router_regret_guard = env_bool("NG_META_ROUTER_REGRET_GUARD", TRUE),
  meta_router_regret_plan_weight = env_num("NG_META_ROUTER_REGRET_PLAN_WEIGHT", 0.60),
  meta_router_regret_gain_weight = env_num("NG_META_ROUTER_REGRET_GAIN_WEIGHT", 0.40),
  meta_router_regret_min_advantage = env_num("NG_META_ROUTER_REGRET_MIN_ADVANTAGE", 0.75),
  meta_router_regret_max_router_penalty = env_num("NG_META_ROUTER_REGRET_MAX_ROUTER_PENALTY", 0.40),
  frontier_policy_spec = env_chr("NG_FRONTIER_POLICY_SPEC", ""),
  frontier_policy_source = env_chr("NG_FRONTIER_POLICY_SOURCE", ng_frontier_default_source()),
  frontier_policy_fallback_method = env_chr("NG_FRONTIER_POLICY_FALLBACK_METHOD", "ng_meta_router_ocs10_lps2"),
  crop_scenario = env_chr("NG_CROP_SCENARIO", ""),
  crop = env_chr("NG_CROP", ""),
  crop_harness_model = env_chr("NG_CROP_HARNESS_MODEL", ""),
  crop_validation_scope = env_chr("NG_CROP_VALIDATION_SCOPE", ""),
  crop_policy_source = env_chr("NG_CROP_POLICY_SOURCE", ng_crop_aware_policy_default_source()),
  topk = env_int("NG_TOPK", NA_integer_),
  topk_prop = env_num("NG_TOPK_PROP", 0.10),
  shared_scoring = env_bool("NG_SHARED_SCORING", TRUE),
  external_tail_prop = env_num("NG_EXTERNAL_TAIL_PROP", 0.10),
  external_shortlist_n = env_int("NG_EXTERNAL_SHORTLIST_N", NA_integer_),
  external_shortlist_multiplier = env_num("NG_EXTERNAL_SHORTLIST_MULTIPLIER", NA_real_),
  external_shortlist_score_col = env_chr("NG_EXTERNAL_SHORTLIST_SCORE_COL", "etk_dh_pmv_var_blend_cal,uc_dh_blend,var_simple,mpv"),
  alphamate_executable = env_chr("NG_ALPHAMATE_EXE", ng_alphamate_default_executable()),
  alphamate_runtime_path = env_chr("NG_ALPHAMATE_RUNTIME_PATH", ""),
  alphamate_criterion_cols = env_csv("NG_ALPHAMATE_CRITERION_COLS", "cross_mean_adjusted_pheno,cross_mean_blend,cross_mean_gebv,cross_mean"),
  alphamate_target_degree = env_num("NG_ALPHAMATE_TARGET_DEGREE", 45),
  alphamate_max_contributions = env_int("NG_ALPHAMATE_MAX_CONTRIBUTIONS", NA_integer_),
  alphamate_number_of_parents = env_int("NG_ALPHAMATE_NUMBER_OF_PARENTS", NA_integer_),
  alphamate_evol_solutions = env_int("NG_ALPHAMATE_EVOL_SOLUTIONS", 100L),
  alphamate_evol_iterations = env_int("NG_ALPHAMATE_EVOL_ITERATIONS", 1000L),
  alphamate_evol_stop = env_int("NG_ALPHAMATE_EVOL_STOP", 200L),
  alphamate_threads = env_int("NG_ALPHAMATE_THREADS", 1L),
  alphamate_keep_files = env_bool("NG_ALPHAMATE_KEEP_FILES", FALSE),
  simplemating_min_cross = env_int("NG_SIMPLEMATING_MIN_CROSS", 1L),
  simplemating_max_search = env_int("NG_SIMPLEMATING_MAX_SEARCH", 100000L),
  simplemating_culling_k = env_num("NG_SIMPLEMATING_CULLING_K", NA_real_),
  simplemating_threads = env_int("NG_SIMPLEMATING_THREADS", 1L),
  methods = strsplit(env_chr(
    "NG_METHODS",
    "var_simple_topn,var_simple_etk_topn,ng_pmv_blend_cal_topn,ng_ocs_mip10_lps1,ng_ocs_mip10_lps2,ng_ocs_mip10_lps4"
  ), ",", fixed = TRUE)[[1]],
  output_prefix = env_chr("NG_OUTPUT_PREFIX", "nextgen_cross_design")
)
cfg$methods <- unique(trimws(cfg$methods[nzchar(cfg$methods)]))
cfg$calibration_pool <- match.arg(tolower(trimws(cfg$calibration_pool)), c("branch", "global"))
cfg$lambda_parent_use_mode <- match.arg(tolower(trimws(cfg$lambda_parent_use_mode)),
                                        c("absolute", "adaptive"))
cfg$family_allocation_method <- match.arg(tolower(trimws(cfg$family_allocation_method)),
                                          c("marginal_topk", "score_weighted"))
if (!is.finite(cfg$simplemating_culling_k)) cfg$simplemating_culling_k <- NA_real_
if (!is.finite(cfg$total_progeny)) cfg$total_progeny <- cfg$top_crosses * cfg$progeny_per_cross
if (!is.finite(cfg$family_selected_top_n)) cfg$family_selected_top_n <- cfg$n_parents
if (!is.finite(cfg$alphasimr_threads) || cfg$alphasimr_threads < 1L) cfg$alphasimr_threads <- 1L
if (!is.finite(cfg$genome_length_m) || cfg$genome_length_m <= 0) cfg$genome_length_m <- 1

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(name = names(cfg), value = vapply(cfg, function(x) paste(x, collapse = ","), character(1))),
          file.path("results", paste0(cfg$output_prefix, "_config.csv")), row.names = FALSE)

message("nextgen_cross_design AlphaSimR config:")
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(rep_id) {
  message("Replicate ", rep_id, " of ", cfg$reps)
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = cfg$genome_length_m)
  SP <- SimParam$new(founder)
  SP$nThreads <- max(1L, as.integer(cfg$alphasimr_threads))
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)

  base_pop <- newPop(founder, simParam = SP)
  base_parents <- make_dh_base_parents(base_pop, cfg$n_parents, SP)
  base_parents@id <- paste0("P0_", seq_len(nInd(base_parents)))
  base_best <- max(as.numeric(gv(base_parents)[, 1]))

  snp_map <- getSnpMap(simParam = SP)
  marker_names <- paste0("Chr", snp_map$chr, "_", snp_map$id)
  marker_map <- data.frame(marker = marker_names, chr = snp_map$chr, pos_cm = snp_map$pos * 100,
                           stringsAsFactors = FALSE)

  branches <- setNames(vector("list", length(cfg$methods)), cfg$methods)
  for (method in cfg$methods) branches[[method]] <- base_parents

  metrics <- list()
  selected_out <- list()
  selection_summary_out <- list()
  family_out <- list()
  timing_out <- list()
  effect_out <- list()
  calibration_out <- list()
  calibration_history <- setNames(vector("list", length(cfg$methods)), cfg$methods)
  global_calibration_history <- NULL
  realization_cache <- new.env(parent = emptyenv())

  for (method in cfg$methods) {
    metrics[[paste(method, 0, sep = "_")]] <- population_metrics(branches[[method]], base_best, 0L, rep_id, method)
  }

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    cycle_global_updates <- list()
    cycle_inputs <- setNames(vector("list", length(cfg$methods)), cfg$methods)
    for (method in cfg$methods) {
      parent_pop <- branches[[method]]
      ids <- method_parent_ids(method, rep_id, cycle, nInd(parent_pop))
      parent_pop@id <- ids
      geno <- pullSnpGeno(parent_pop, simParam = SP)
      colnames(geno) <- marker_names
      rownames(geno) <- ids
      parent_y <- phenotype(
        parent_pop,
        h2 = cfg$phenotype_h2,
        reps = cfg$phenotype_reps,
        seed = cfg$seed + rep_id * 100000L + cycle * 1000L + 11L
      )
      names(parent_y) <- ids
      parent_key <- population_fingerprint(geno, parent_pop)
      hist_key <- history_cache_key(method, cfg, calibration_history, global_calibration_history)
      cache_key <- paste(parent_key, hist_key, sep = "\n")
      if (!isTRUE(cfg$shared_scoring)) cache_key <- paste(cache_key, method, sep = "\n")
      cycle_inputs[[method]] <- list(
        method = method,
        parent_pop = parent_pop,
        ids = ids,
        geno = geno,
        parent_y = parent_y,
        cache_key = cache_key
      )
    }

    cache_keys <- vapply(cycle_inputs, function(x) x$cache_key, character(1))
    unique_cache_keys <- unique(cache_keys[cfg$methods])
    score_groups <- setNames(lapply(unique_cache_keys, function(key) {
      cfg$methods[cache_keys[cfg$methods] == key]
    }), unique_cache_keys)
    score_cache <- vector("list", length(score_groups))
    names(score_cache) <- names(score_groups)

    for (group_id in seq_along(score_groups)) {
      group_methods <- score_groups[[group_id]]
      representative <- cycle_inputs[[group_methods[[1]]]]
      score_seed <- stable_benchmark_seed(
        "score_group",
        cfg$seed,
        rep_id,
        cycle,
        representative$cache_key,
        seed = cfg$seed
      )
      shared_ids <- shared_parent_ids(rep_id, cycle, group_id, nInd(representative$parent_pop))
      scoring_pop <- representative$parent_pop
      scoring_pop@id <- shared_ids
      geno <- representative$geno
      rownames(geno) <- shared_ids
      parent_y <- representative$parent_y
      names(parent_y) <- shared_ids

      t_eff <- system.time({
        training_pop <- with_rng_seed(
          stable_benchmark_seed("training_pop", score_seed, seed = cfg$seed),
          make_training_population(scoring_pop, cfg$effect_training_n, SP)
        )
        training_ids <- shared_training_ids(rep_id, cycle, group_id, nInd(training_pop))
        training_pop@id <- training_ids
        training_geno <- pullSnpGeno(training_pop, simParam = SP)
        colnames(training_geno) <- marker_names
        rownames(training_geno) <- training_ids
        training_y <- phenotype(
          training_pop,
          h2 = cfg$phenotype_h2,
          reps = cfg$phenotype_reps,
          seed = stable_benchmark_seed("training_y", score_seed, seed = cfg$seed)
        )
        names(training_y) <- training_ids
        effects <- ng_fit_ridge_effects(training_geno, training_y, ids = training_ids,
                                        h2_prior = cfg$effect_h2_prior,
                                        kfold = cfg$effect_kfold,
                                        seed = stable_benchmark_seed("ridge", score_seed, seed = cfg$seed))
      })[["elapsed"]]

      t_score <- system.time({
        scores <- ng_score_crosses(
          geno = geno,
          effects = effects,
          marker_map = marker_map,
          ids = shared_ids,
          adjusted_pheno = parent_y,
          selection_prop = cfg$selection_prop,
          min_effect_reliability = cfg$min_effect_reliability,
          use_cpp = cfg$use_cpp
        )
        history_for_calibration <- if (cfg$calibration_pool == "global") {
          global_calibration_history
        } else {
          calibration_history[[group_methods[[1]]]]
        }
        calibrators <- ng_fit_family_variance_calibrators(
          history_for_calibration,
          min_n = cfg$calibration_min_n
        )
        scores <- ng_apply_family_variance_calibrators(
          scores,
          calibrators = calibrators,
          n_progeny = cfg$progeny_per_cross,
          top_k = if (is.na(cfg$topk)) NULL else cfg$topk,
          top_prop = cfg$topk_prop
        )
        scores <- ng_add_local_portfolio_scores(
          scores,
          history = history_for_calibration,
          target_col = cfg$portfolio_target,
          min_history_n = cfg$portfolio_min_history_n
        )
        scores <- ng_add_external_baseline_scores(
          scores = scores,
          geno = geno,
          effects = effects,
          marker_map = marker_map,
          adjusted_pheno = parent_y,
          methods = group_methods,
          selection_prop = cfg$selection_prop,
          popvar_tail_p = cfg$external_tail_prop,
          simplemating_threads = cfg$simplemating_threads,
          n_crosses = cfg$top_crosses,
          shortlist_n = if (is.na(cfg$external_shortlist_n)) NULL else cfg$external_shortlist_n,
          shortlist_multiplier = if (is.na(cfg$external_shortlist_multiplier)) NULL else cfg$external_shortlist_multiplier,
          shortlist_score_col = cfg$external_shortlist_score_col
        )
        scores <- ng_add_adaptive_stack_scores(
          scores,
          history = history_for_calibration,
          target_col = cfg$adaptive_target,
          min_history_n = cfg$adaptive_min_history_n,
          score_cols = cfg$adaptive_score_cols,
          n_parents = length(shared_ids),
          n_crosses = cfg$top_crosses,
          history_weight = if (is.na(cfg$adaptive_history_weight)) NULL else cfg$adaptive_history_weight,
          fallback_col = cfg$adaptive_fallback_col,
          fallback_max_weight = cfg$adaptive_fallback_max_weight
        )
        scores <- ng_add_meta_portfolio_scores(
          scores,
          history = history_for_calibration,
          method_history = global_calibration_history,
          target_col = cfg$meta_target,
          min_history_n = cfg$meta_min_history_n,
          method_min_history_n = cfg$meta_method_min_history_n,
          method_recent_cycles = cfg$meta_method_recent_cycles,
          score_cols = cfg$meta_score_cols,
          n_parents = length(shared_ids),
          n_crosses = cfg$top_crosses,
          history_weight = if (is.na(cfg$meta_history_weight)) NULL else cfg$meta_history_weight
        )
        attr(scores, "meta_method_history") <- global_calibration_history
      })[["elapsed"]]

      score_cache[[group_id]] <- list(
        scores = scores,
        effects = effects,
        parent_K = ng_parent_kinship(geno),
        shared_ids = shared_ids,
        elapsed_effects = t_eff,
        elapsed_score = t_score,
        group_methods = group_methods,
        group_id = group_id
      )
    }

    for (method in cfg$methods) {
      input <- cycle_inputs[[method]]
      bundle <- score_cache[[input$cache_key]]
      group_size <- length(bundle$group_methods)
      effects <- bundle$effects
      scores <- remap_score_parent_ids(bundle$scores, bundle$shared_ids, input$ids)
      parent_K <- rename_kinship_ids(bundle$parent_K, input$ids)
      t_eff <- bundle$elapsed_effects / group_size
      t_score <- bundle$elapsed_score / group_size

      t_select <- system.time({
        selected <- select_branch(method, scores, cfg, parent_K)
      })[["elapsed"]]
      plan_summary <- attr(selected, "summary")
      pp <- method_parent_penalty(method, cfg$lambda_parent_use, cfg$lambda_parent_use_mode)

      selected$rep <- rep_id
      selected$cycle <- cycle
      selected$method <- method
      selected_out[[paste(rep_id, method, cycle, sep = "_")]] <- selected
      selected_counts <- ng_parent_counts(selected, rownames(parent_K))
      selected_contrib <- selected_counts / sum(selected_counts)
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
        lambda_parent_use = if (!is.null(plan_summary)) plan_summary$lambda_parent_use else if (grepl("ocs", method)) pp$value else NA_real_,
        lambda_parent_use_input = if (!is.null(plan_summary)) plan_summary$lambda_parent_use_input else if (grepl("ocs", method)) pp$value else NA_real_,
        lambda_parent_use_mode = if (!is.null(plan_summary)) plan_summary$lambda_parent_use_mode else if (grepl("ocs", method)) pp$mode else NA_character_,
        score_scale = if (!is.null(plan_summary)) plan_summary$score_scale else NA_real_,
        balanced_gain_col = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_gain_col)) plan_summary$balanced_gain_col else NA_character_,
        balanced_diversity_col = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_diversity_col)) plan_summary$balanced_diversity_col else NA_character_,
        balanced_gain_weight = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_gain_weight)) plan_summary$balanced_gain_weight else NA_real_,
        balanced_diversity_weight = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_diversity_weight)) plan_summary$balanced_diversity_weight else NA_real_,
        balanced_pair_kinship_weight = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_pair_kinship_weight)) plan_summary$balanced_pair_kinship_weight else NA_real_,
        balanced_min_unique_used = if (!is.null(plan_summary) && !is.null(plan_summary$balanced_min_unique_used)) plan_summary$balanced_min_unique_used else NA_integer_,
        adaptive_score_cols = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_score_cols)) plan_summary$adaptive_score_cols else NA_character_,
        adaptive_weights = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_weights)) plan_summary$adaptive_weights else NA_character_,
        adaptive_reliability = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_reliability)) plan_summary$adaptive_reliability else NA_real_,
        adaptive_history_weight = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_history_weight)) plan_summary$adaptive_history_weight else NA_real_,
        adaptive_history_n = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_history_n)) plan_summary$adaptive_history_n else NA_integer_,
        adaptive_champion_weight = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_champion_weight)) plan_summary$adaptive_champion_weight else NA_real_,
        adaptive_fallback_col = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_fallback_col)) plan_summary$adaptive_fallback_col else NA_character_,
        adaptive_fallback_weight = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_fallback_weight)) plan_summary$adaptive_fallback_weight else NA_real_,
        adaptive_parent_use_input = if (!is.null(plan_summary) && !is.null(plan_summary$adaptive_parent_use_input)) plan_summary$adaptive_parent_use_input else NA_real_,
        meta_score_cols = if (!is.null(plan_summary) && !is.null(plan_summary$meta_score_cols)) plan_summary$meta_score_cols else NA_character_,
        meta_weights = if (!is.null(plan_summary) && !is.null(plan_summary$meta_weights)) plan_summary$meta_weights else NA_character_,
        meta_reliability = if (!is.null(plan_summary) && !is.null(plan_summary$meta_reliability)) plan_summary$meta_reliability else NA_real_,
        meta_history_weight = if (!is.null(plan_summary) && !is.null(plan_summary$meta_history_weight)) plan_summary$meta_history_weight else NA_real_,
        meta_history_n = if (!is.null(plan_summary) && !is.null(plan_summary$meta_history_n)) plan_summary$meta_history_n else NA_integer_,
        meta_method_history_weight = if (!is.null(plan_summary) && !is.null(plan_summary$meta_method_history_weight)) plan_summary$meta_method_history_weight else NA_real_,
        meta_method_history_n = if (!is.null(plan_summary) && !is.null(plan_summary$meta_method_history_n)) plan_summary$meta_method_history_n else NA_integer_,
        meta_method_weights = if (!is.null(plan_summary) && !is.null(plan_summary$meta_method_weights)) plan_summary$meta_method_weights else NA_character_,
        meta_leader_col = if (!is.null(plan_summary) && !is.null(plan_summary$meta_leader_col)) plan_summary$meta_leader_col else NA_character_,
        meta_leader_weight = if (!is.null(plan_summary) && !is.null(plan_summary$meta_leader_weight)) plan_summary$meta_leader_weight else NA_real_,
        meta_selector_gain_col = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_gain_col)) plan_summary$meta_selector_gain_col else NA_character_,
        meta_selector_allocator = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_allocator)) plan_summary$meta_selector_allocator else NA_character_,
        meta_selector_decision = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_decision)) plan_summary$meta_selector_decision else NA_character_,
        meta_selector_confidence = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_confidence)) plan_summary$meta_selector_confidence else NA_real_,
        meta_selector_leader_margin = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_leader_margin)) plan_summary$meta_selector_leader_margin else NA_real_,
        meta_selector_portfolio_score_delta = if (!is.null(plan_summary) && !is.null(plan_summary$meta_selector_portfolio_score_delta)) plan_summary$meta_selector_portfolio_score_delta else NA_real_,
        meta_router_family = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_family)) plan_summary$meta_router_family else NA_character_,
        meta_router_original_family = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_original_family)) plan_summary$meta_router_original_family else NA_character_,
        meta_router_gain_col = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_gain_col)) plan_summary$meta_router_gain_col else NA_character_,
        meta_router_allocator = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_allocator)) plan_summary$meta_router_allocator else NA_character_,
      meta_router_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_score)) plan_summary$meta_router_score else NA_real_,
      meta_router_history_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_history_score)) plan_summary$meta_router_history_score else NA_real_,
      meta_router_prior_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_prior_score)) plan_summary$meta_router_prior_score else NA_real_,
      meta_router_plan_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_plan_score)) plan_summary$meta_router_plan_score else NA_real_,
      meta_router_gain_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_gain_score)) plan_summary$meta_router_gain_score else NA_real_,
        meta_router_diversity_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_diversity_score)) plan_summary$meta_router_diversity_score else NA_real_,
        meta_router_reliability_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_reliability_score)) plan_summary$meta_router_reliability_score else NA_real_,
        meta_router_history_n = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_history_n)) plan_summary$meta_router_history_n else NA_integer_,
        meta_router_history_cycles = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_history_cycles)) plan_summary$meta_router_history_cycles else NA_integer_,
        meta_router_reliability = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_reliability)) plan_summary$meta_router_reliability else NA_real_,
        meta_router_guard_override = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_guard_override)) plan_summary$meta_router_guard_override else NA_integer_,
        meta_router_guard_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_guard_score)) plan_summary$meta_router_guard_score else NA_real_,
        meta_router_original_guard_score = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_original_guard_score)) plan_summary$meta_router_original_guard_score else NA_real_,
        meta_router_guard_advantage = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_guard_advantage)) plan_summary$meta_router_guard_advantage else NA_real_,
        meta_router_guard_router_penalty = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_guard_router_penalty)) plan_summary$meta_router_guard_router_penalty else NA_real_,
        meta_router_guard_reason = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_guard_reason)) plan_summary$meta_router_guard_reason else NA_character_,
        meta_router_candidates = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_candidates)) plan_summary$meta_router_candidates else NA_character_,
        meta_router_candidate_details = if (!is.null(plan_summary) && !is.null(plan_summary$meta_router_candidate_details)) plan_summary$meta_router_candidate_details else NA_character_,
        frontier_policy_method = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_method)) plan_summary$frontier_policy_method else NA_character_,
        frontier_policy_family = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_family)) plan_summary$frontier_policy_family else NA_character_,
        frontier_policy_primary_method = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_primary_method)) plan_summary$frontier_policy_primary_method else NA_character_,
        frontier_policy_primary_family = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_primary_family)) plan_summary$frontier_policy_primary_family else NA_character_,
        frontier_policy_source = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_source)) plan_summary$frontier_policy_source else NA_character_,
        frontier_policy_reason = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_reason)) plan_summary$frontier_policy_reason else NA_character_,
        frontier_policy_band = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_band)) plan_summary$frontier_policy_band else NA_character_,
        frontier_policy_fallback = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_fallback)) plan_summary$frontier_policy_fallback else NA_integer_,
        frontier_policy_candidates = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_candidates)) plan_summary$frontier_policy_candidates else NA_character_,
        frontier_policy_error_log = if (!is.null(plan_summary) && !is.null(plan_summary$frontier_policy_error_log)) plan_summary$frontier_policy_error_log else NA_character_,
        crop_policy_method = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_method)) plan_summary$crop_policy_method else NA_character_,
        crop_policy_family = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_family)) plan_summary$crop_policy_family else NA_character_,
        crop_policy_primary_method = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_primary_method)) plan_summary$crop_policy_primary_method else NA_character_,
        crop_policy_primary_family = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_primary_family)) plan_summary$crop_policy_primary_family else NA_character_,
        crop_policy_source = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_source)) plan_summary$crop_policy_source else NA_character_,
        crop_policy_reason = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_reason)) plan_summary$crop_policy_reason else NA_character_,
        crop_policy_band = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_band)) plan_summary$crop_policy_band else NA_character_,
        crop_policy_fallback = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_fallback)) plan_summary$crop_policy_fallback else NA_integer_,
        crop_policy_candidates = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_candidates)) plan_summary$crop_policy_candidates else NA_character_,
        crop_policy_error_log = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_error_log)) plan_summary$crop_policy_error_log else NA_character_,
        crop_policy_crop_scenario = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_crop_scenario)) plan_summary$crop_policy_crop_scenario else NA_character_,
        crop_policy_crop = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_crop)) plan_summary$crop_policy_crop else NA_character_,
        crop_policy_harness_model = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_harness_model)) plan_summary$crop_policy_harness_model else NA_character_,
        crop_policy_validation_scope = if (!is.null(plan_summary) && !is.null(plan_summary$crop_policy_validation_scope)) plan_summary$crop_policy_validation_scope else NA_character_,
        alphamate_mode = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_mode)) plan_summary$alphamate_mode else NA_character_,
        alphamate_target_degree = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_target_degree)) plan_summary$alphamate_target_degree else NA_real_,
        alphamate_criterion_col = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_criterion_col)) plan_summary$alphamate_criterion_col else NA_character_,
        alphamate_exit_code = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_exit_code)) plan_summary$alphamate_exit_code else NA_integer_,
        alphamate_executable = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_executable)) plan_summary$alphamate_executable else NA_character_,
        alphamate_runtime_path = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_runtime_path)) plan_summary$alphamate_runtime_path else NA_character_,
        alphamate_workdir = if (!is.null(plan_summary) && !is.null(plan_summary$alphamate_workdir)) plan_summary$alphamate_workdir else NA_character_,
        meta_champion_weight = if (!is.null(plan_summary) && !is.null(plan_summary$meta_champion_weight)) plan_summary$meta_champion_weight else NA_real_,
        meta_parent_use_input = if (!is.null(plan_summary) && !is.null(plan_summary$meta_parent_use_input)) plan_summary$meta_parent_use_input else NA_real_,
        total_progeny = if ("n_progeny" %in% names(selected)) {
          sum(selected$n_progeny, na.rm = TRUE)
        } else {
          cfg$top_crosses * cfg$progeny_per_cross
        },
        stringsAsFactors = FALSE
      )

      timing_out[[paste(rep_id, method, cycle, "effects", sep = "_")]] <- data.frame(
        rep = rep_id, cycle = cycle, method = method, stage = "effects", elapsed_sec = t_eff,
        shared_group_id = bundle$group_id,
        shared_group_size = group_size,
        shared_group_methods = paste(bundle$group_methods, collapse = ","),
        stringsAsFactors = FALSE
      )
      timing_out[[paste(rep_id, method, cycle, "score", sep = "_")]] <- data.frame(
        rep = rep_id, cycle = cycle, method = method, stage = "score", elapsed_sec = t_score,
        shared_group_id = bundle$group_id,
        shared_group_size = group_size,
        shared_group_methods = paste(bundle$group_methods, collapse = ","),
        stringsAsFactors = FALSE
      )
      timing_out[[paste(rep_id, method, cycle, "select", sep = "_")]] <- data.frame(
        rep = rep_id, cycle = cycle, method = method, stage = "select", elapsed_sec = t_select,
        shared_group_id = bundle$group_id,
        shared_group_size = group_size,
        shared_group_methods = paste(bundle$group_methods, collapse = ","),
        stringsAsFactors = FALSE
      )
      effect_out[[paste(rep_id, method, cycle, sep = "_")]] <- data.frame(
        rep = rep_id, cycle = cycle, method = method,
        effect_reliability = effects$reliability,
        in_sample_reliability = effects$in_sample_reliability,
        ridge_lambda = effects$lambda,
        mean_source = unique(scores$mean_source)[1],
        stringsAsFactors = FALSE
      )
      cal_sum <- ng_calibrator_summary(attr(scores, "variance_calibrators"))
      cal_sum$rep <- rep_id
      cal_sum$cycle <- cycle
      cal_sum$method <- method
      cal_sum$topk_intensity <- attr(scores, "topk_intensity")
      calibration_out[[paste(rep_id, method, cycle, sep = "_")]] <- cal_sum

      progeny_result <- make_selected_progeny(
        input$parent_pop,
        selected,
        cfg$progeny_per_cross,
        SP,
        method,
        rep_id = rep_id,
        cycle = cycle,
        seed = cfg$seed,
        realization_cache = realization_cache
      )
      families <- progeny_result$families
      families$rep <- rep_id
      families$cycle <- cycle
      families$method <- method
      family_out[[paste(rep_id, method, cycle, sep = "_")]] <- families
      calibration_history[[method]] <- bind_rows_fill(list(calibration_history[[method]], families))
      cycle_global_updates[[method]] <- families

      progeny <- progeny_result$pop
      g <- as.numeric(gv(progeny)[, 1])
      keep <- order(g, decreasing = TRUE)[seq_len(min(cfg$n_parents, length(g)))]
      next_pop <- progeny[keep]
      next_pop@id <- paste0(method, "_R", rep_id, "_C", cycle, "_P", seq_len(nInd(next_pop)))
      branches[[method]] <- next_pop
      metrics[[paste(method, cycle, sep = "_")]] <- population_metrics(next_pop, base_best, cycle, rep_id, method)
    }
    global_calibration_history <- bind_rows_fill(list(
      global_calibration_history,
      bind_rows_fill(cycle_global_updates)
    ))
  }

  list(metrics = bind_rows_fill(metrics),
       selections = bind_rows_fill(selected_out),
       selection_summary = bind_rows_fill(selection_summary_out),
       families = bind_rows_fill(family_out),
       timings = bind_rows_fill(timing_out),
       effects = bind_rows_fill(effect_out),
       calibrations = bind_rows_fill(calibration_out))
})

metrics <- bind_rows_fill(lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
selection_summary <- bind_rows_fill(lapply(rep_results, `[[`, "selection_summary"))
families <- bind_rows_fill(lapply(rep_results, `[[`, "families"))
timings <- bind_rows_fill(lapply(rep_results, `[[`, "timings"))
effects <- bind_rows_fill(lapply(rep_results, `[[`, "effects"))
calibration_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "calibrations"))

overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ method + cycle,
  metrics,
  mean
)

baseline <- overall[overall$method == "var_simple_topn",
                    c("cycle", "mean_gv", "max_gv", "top10_gv", "transgressive_rate"),
                    drop = FALSE]
comparison <- merge(overall, baseline, by = "cycle", suffixes = c("", "_var_simple"))
comparison$delta_mean_gv <- comparison$mean_gv - comparison$mean_gv_var_simple
comparison$delta_max_gv <- comparison$max_gv - comparison$max_gv_var_simple
comparison$delta_top10_gv <- comparison$top10_gv - comparison$top10_gv_var_simple
comparison$delta_transgressive_rate <- comparison$transgressive_rate - comparison$transgressive_rate_var_simple

calibration <- if (nrow(families)) {
  do.call(rbind, by(families, list(families$method, families$cycle), function(d) {
    out <- ng_validate_metric_calibration(
      scores = data.frame(parent1 = d$parent1, parent2 = d$parent2,
                          var_simple = d$pred_var_simple,
                          dh_recomb_var = d$pred_dh_recomb_var,
                          dh_pmv_var = d$pred_dh_pmv_var,
                          var_simple_cal = d$pred_var_simple_cal,
                          dh_recomb_var_cal = d$pred_dh_recomb_var_cal,
                          dh_pmv_var_cal = d$pred_dh_pmv_var_cal),
      realized = data.frame(parent1 = d$parent1, parent2 = d$parent2,
                            realized_var = d$realized_var,
                            realized_mean = d$realized_mean),
      pred_cols = c(
        "var_simple", "dh_recomb_var", "dh_pmv_var",
        "var_simple_cal", "dh_recomb_var_cal", "dh_pmv_var_cal"
      )
    )
    out$method <- d$method[1]
    out$cycle <- d$cycle[1]
    out
  }))
} else data.frame()

timing_summary <- if (nrow(timings)) {
  stats::aggregate(elapsed_sec ~ method + stage, timings, mean)
} else data.frame()

family_summary <- if (nrow(families)) {
  summary_cols <- c(
    "pred_var", "pred_var_simple", "pred_dh_recomb_var", "pred_dh_pmv_var",
    "pred_var_simple_cal", "pred_dh_recomb_var_cal", "pred_dh_pmv_var_cal",
    "pred_balanced_usefulness_score", "pred_ng_adaptive_score", "pred_ng_adaptive_score_raw",
    "pred_ng_adaptive_var", "pred_ng_adaptive_fallback_weight",
    "pred_ng_meta_score", "pred_ng_meta_var", "pred_ng_meta_method_history_weight",
    "pred_ng_meta_leader_score", "pred_ng_meta_leader_var",
    "pred_ng_meta_selector_score", "pred_ng_meta_selector_var",
    "pred_ng_meta_router_score", "pred_ng_meta_router_var",
    "pred_ng_meta_router_history_score", "pred_ng_meta_router_prior_score", "pred_ng_meta_router_plan_score",
    "pred_ng_meta_router_gain_score", "pred_ng_meta_router_diversity_score",
    "pred_ng_meta_router_reliability_score", "pred_ng_meta_router_guard_override",
    "pred_ng_meta_router_guard_score", "pred_ng_meta_router_original_guard_score",
    "pred_ng_meta_router_guard_advantage", "pred_ng_meta_router_guard_router_penalty",
    "pred_ng_frontier_policy_fallback", "pred_ng_crop_policy_fallback",
    "realized_var", "realized_mean", "realized_top10", "realized_max", "pair_kinship", "n_progeny"
  )
  summary_cols <- intersect(summary_cols, names(families))
  summary_cols <- summary_cols[vapply(summary_cols, function(col) {
    any(is.finite(as.numeric(families[[col]])))
  }, logical(1))]
  stats::aggregate(
    families[, summary_cols, drop = FALSE],
    by = families[c("method", "cycle")],
    FUN = mean,
    na.rm = TRUE
  )
} else data.frame()

prefix <- file.path("results", cfg$output_prefix)
write.csv(metrics, paste0(prefix, "_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, paste0(prefix, "_overall.csv"), row.names = FALSE)
write.csv(comparison, paste0(prefix, "_comparison_vs_var_simple.csv"), row.names = FALSE)
write.csv(selections, paste0(prefix, "_selected_crosses.csv"), row.names = FALSE)
write.csv(selection_summary, paste0(prefix, "_selection_summary.csv"), row.names = FALSE)
write.csv(families, paste0(prefix, "_families.csv"), row.names = FALSE)
write.csv(family_summary, paste0(prefix, "_family_summary.csv"), row.names = FALSE)
write.csv(calibration, paste0(prefix, "_calibration.csv"), row.names = FALSE)
write.csv(timings, paste0(prefix, "_timings.csv"), row.names = FALSE)
write.csv(timing_summary, paste0(prefix, "_timing_summary.csv"), row.names = FALSE)
write.csv(effects, paste0(prefix, "_effect_diagnostics.csv"), row.names = FALSE)
write.csv(calibration_diagnostics, paste0(prefix, "_variance_calibrators.csv"), row.names = FALSE)

message("Overall:")
print(overall)
message("Comparison versus var_simple_topn:")
print(comparison[order(comparison$cycle, comparison$method), , drop = FALSE])
message("Calibration:")
print(calibration)
message("Timing summary:")
print(timing_summary)
message("Selection summary:")
print(selection_summary)
