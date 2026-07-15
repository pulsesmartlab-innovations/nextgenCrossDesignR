ng_head_to_head_default_methods <- function() {
  data.frame(
    label = c(
      "nextgen_auto_ocs",
      "nextgen_economic_index_ocs",
      "nextgen_desired_gain_ocs",
      "popvar_style_weighted_topn",
      "simplemate_style_threshold_topn",
      "alphamate_style_weighted_ocs"
    ),
    method = c(
      "auto",
      "economic_index",
      "desired_gain",
      "weighted",
      "threshold",
      "weighted"
    ),
    allocator = c(
      "ocs",
      "ocs",
      "ocs",
      "topn",
      "topn",
      "ocs"
    ),
    role = c(
      "candidate",
      "candidate",
      "candidate",
      "baseline",
      "baseline",
      "baseline"
    ),
    external_tool = c(
      "NextGen",
      "NextGen",
      "NextGen",
      "PopVar",
      "SimpleMating",
      "AlphaMate"
    ),
    implementation = c(
      "native_multitrait",
      "native_multitrait",
      "native_multitrait",
      "style_proxy",
      "style_proxy",
      "style_proxy"
    ),
    exact_external_status = c(
      "not_applicable",
      "not_applicable",
      "not_applicable",
      "not_run_multitrait_ci_proxy",
      "not_run_multitrait_ci_proxy",
      "not_run_multitrait_ci_proxy"
    ),
    fallback_reason = c(
      "",
      "",
      "",
      "PopVar package scoring is optional and single-trait oriented in this repo; CI uses an economic-weight scalar proxy for multi-trait head-to-head contracts.",
      "SimpleMating package scoring is optional and single-trait oriented in this repo; CI uses threshold/usefulness-style multi-trait ranking as the proxy.",
      "AlphaMate executable is optional and optimizes a scalar criterion; CI uses the same scalar weighted index through the native relationship-aware OCS allocator."
    ),
    description = c(
      "NextGen default rank-normalized threshold-aware OCS policy for unknown weights.",
      "NextGen covariance-aware economic-index OCS policy.",
      "NextGen covariance-aware desired-gain OCS policy.",
      "PopVar-style scalar cross ranking using fixed economic weights and no allocation constraints.",
      "SimpleMating-style threshold/usefulness ranking using multi-trait threshold penalties.",
      "AlphaMate-style scalar index optimized with relationship-aware OCS allocation."
    ),
    stringsAsFactors = FALSE
  )
}

ng_head_to_head_normalize_methods <- function(methods = NULL) {
  required <- c(
    "label", "method", "allocator", "role", "external_tool",
    "implementation", "exact_external_status", "fallback_reason", "description"
  )
  registry <- ng_head_to_head_default_methods()
  if (is.null(methods)) {
    out <- registry
  } else if (is.data.frame(methods)) {
    out <- as.data.frame(methods, stringsAsFactors = FALSE)
    missing <- setdiff(required, names(out))
    if (length(missing)) ng_stop("head-to-head method registry missing columns: ", paste(missing, collapse = ", "))
    out <- out[required]
  } else {
    labels <- unique(trimws(as.character(methods)))
    labels <- labels[nzchar(labels) & !is.na(labels)]
    missing <- setdiff(labels, registry$label)
    if (length(missing)) ng_stop("unknown head-to-head method label(s): ", paste(missing, collapse = ", "))
    out <- registry[match(labels, registry$label), , drop = FALSE]
  }

  out$label <- trimws(as.character(out$label))
  out$method <- trimws(tolower(as.character(out$method)))
  out$allocator <- trimws(tolower(as.character(out$allocator)))
  out$role <- trimws(tolower(as.character(out$role)))
  out$external_tool <- trimws(as.character(out$external_tool))
  out$implementation <- trimws(as.character(out$implementation))
  out$exact_external_status <- trimws(as.character(out$exact_external_status))
  out$fallback_reason <- trimws(as.character(out$fallback_reason))
  out$description <- trimws(as.character(out$description))
  if (any(!nzchar(out$label) | is.na(out$label))) ng_stop("head-to-head method labels must be non-empty")
  if (any(duplicated(out$label))) ng_stop("head-to-head method labels must be unique")
  valid_methods <- c("auto", "weighted", "economic_index", "desired_gain", "threshold")
  invalid_method <- setdiff(out$method, valid_methods)
  if (length(invalid_method)) ng_stop("unsupported head-to-head multi-trait method(s): ", paste(invalid_method, collapse = ", "))
  invalid_allocator <- setdiff(out$allocator, c("topn", "ocs"))
  if (length(invalid_allocator)) ng_stop("unsupported head-to-head allocator(s): ", paste(invalid_allocator, collapse = ", "))
  invalid_role <- setdiff(out$role, c("candidate", "baseline"))
  if (length(invalid_role)) ng_stop("unsupported head-to-head method role(s): ", paste(invalid_role, collapse = ", "))
  rownames(out) <- NULL
  out
}

ng_head_to_head_metric_directions <- function() {
  c(
    mean_realized_index = "maximize",
    mean_realized_yield = "maximize",
    mean_realized_disease = "minimize",
    mean_realized_quality = "maximize",
    unique_parents = "maximize",
    max_parent_use = "minimize",
    mean_pair_kinship = "minimize",
    group_coancestry = "minimize"
  )
}

ng_head_to_head_compare_summary <- function(summary,
                                            baseline_methods = NULL,
                                            metrics = NULL,
                                            tolerance = 1e-12) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  required <- c("scenario", "crop", "harness_model", "n_parents", "rep", "method")
  missing <- setdiff(required, names(summary))
  if (length(missing)) ng_stop("head-to-head summary missing columns: ", paste(missing, collapse = ", "))
  if (is.null(baseline_methods)) {
    if ("benchmark_role" %in% names(summary)) {
      baseline_methods <- unique(as.character(summary$method[summary$benchmark_role == "baseline"]))
    } else {
      baseline_methods <- character()
    }
  }
  baseline_methods <- unique(trimws(as.character(baseline_methods)))
  baseline_methods <- baseline_methods[nzchar(baseline_methods) & !is.na(baseline_methods)]
  if (!length(baseline_methods)) ng_stop("at least one baseline method is required for head-to-head comparisons")

  directions <- ng_head_to_head_metric_directions()
  if (is.null(metrics)) metrics <- names(directions)
  metrics <- unique(trimws(as.character(metrics)))
  metrics <- metrics[nzchar(metrics) & !is.na(metrics)]
  missing_metrics <- setdiff(metrics, names(summary))
  if (length(missing_metrics)) {
    metrics <- setdiff(metrics, missing_metrics)
  }
  metrics <- metrics[metrics %in% names(directions)]
  if (!length(metrics)) ng_stop("no supported head-to-head comparison metrics are present in summary")

  keys <- unique(summary[, required[required != "method"], drop = FALSE])
  rows <- list()
  k <- 0L
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, , drop = FALSE]
    in_key <- summary$scenario == key$scenario[[1]] &
      summary$n_parents == key$n_parents[[1]] &
      summary$rep == key$rep[[1]]
    if ("crop" %in% names(summary)) in_key <- in_key & summary$crop == key$crop[[1]]
    if ("harness_model" %in% names(summary)) in_key <- in_key & summary$harness_model == key$harness_model[[1]]
    part <- summary[in_key, , drop = FALSE]
    for (baseline in baseline_methods) {
      baseline_row <- part[part$method == baseline, , drop = FALSE]
      if (nrow(baseline_row) != 1L) next
      for (metric in metrics) {
        baseline_value <- suppressWarnings(as.numeric(baseline_row[[metric]][[1]]))
        if (!is.finite(baseline_value)) next
        for (j in seq_len(nrow(part))) {
          method <- as.character(part$method[[j]])
          if (identical(method, baseline)) next
          method_role <- if ("benchmark_role" %in% names(part)) part$benchmark_role[[j]] else NA_character_
          if (!identical(method_role, "candidate")) next
          value <- suppressWarnings(as.numeric(part[[metric]][[j]]))
          if (!is.finite(value)) next
          direction <- directions[[metric]]
          delta <- if (identical(direction, "minimize")) baseline_value - value else value - baseline_value
          k <- k + 1L
          rows[[k]] <- data.frame(
            scenario = key$scenario[[1]],
            crop = key$crop[[1]],
            harness_model = key$harness_model[[1]],
            n_parents = key$n_parents[[1]],
            rep = key$rep[[1]],
            seed = if ("seed" %in% names(part)) part$seed[[j]] else NA_integer_,
            metric = metric,
            direction = direction,
            method = method,
            method_role = method_role,
            method_external_tool = if ("external_tool" %in% names(part)) part$external_tool[[j]] else NA_character_,
            method_implementation = if ("implementation" %in% names(part)) part$implementation[[j]] else NA_character_,
            method_exact_external_status = if ("exact_external_status" %in% names(part)) part$exact_external_status[[j]] else NA_character_,
            method_fallback_reason = if ("fallback_reason" %in% names(part)) part$fallback_reason[[j]] else NA_character_,
            baseline_method = baseline,
            baseline_external_tool = if ("external_tool" %in% names(baseline_row)) baseline_row$external_tool[[1]] else NA_character_,
            baseline_implementation = if ("implementation" %in% names(baseline_row)) baseline_row$implementation[[1]] else NA_character_,
            baseline_exact_external_status = if ("exact_external_status" %in% names(baseline_row)) baseline_row$exact_external_status[[1]] else NA_character_,
            baseline_fallback_reason = if ("fallback_reason" %in% names(baseline_row)) baseline_row$fallback_reason[[1]] else NA_character_,
            value = value,
            baseline_value = baseline_value,
            delta = delta,
            better_than_baseline = delta > tolerance,
            tied_with_baseline = abs(delta) <= tolerance,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      scenario = character(),
      crop = character(),
      harness_model = character(),
      n_parents = integer(),
      rep = integer(),
      seed = integer(),
      metric = character(),
      direction = character(),
      method = character(),
      method_role = character(),
      method_external_tool = character(),
      method_implementation = character(),
      method_exact_external_status = character(),
      method_fallback_reason = character(),
      baseline_method = character(),
      baseline_external_tool = character(),
      baseline_implementation = character(),
      baseline_exact_external_status = character(),
      baseline_fallback_reason = character(),
      value = numeric(),
      baseline_value = numeric(),
      delta = numeric(),
      better_than_baseline = logical(),
      tied_with_baseline = logical(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

ng_head_to_head_run_case <- function(scenario_name,
                                     n_parents,
                                     rep,
                                     run_seed,
                                     methods,
                                     n_crosses,
                                     realized_progeny,
                                     n_founders,
                                     n_chr,
                                     seg_sites,
                                     snp_per_chr,
                                     qtl_per_chr,
                                     genome_length_m,
                                     prediction_noise,
                                     alphasimr_threads,
                                     ocs_lambda_group,
                                     ocs_lambda_mating) {
  scenario_obj <- ng_multitrait_crop_validation_scenario(
    scenario = scenario_name,
    n_parents = n_parents,
    n_founders = n_founders,
    n_chr = n_chr,
    seg_sites = seg_sites,
    snp_per_chr = snp_per_chr,
    qtl_per_chr = qtl_per_chr,
    genome_length_m = genome_length_m,
    realized_progeny = realized_progeny,
    seed = run_seed,
    prediction_noise = prediction_noise,
    alphasimr_threads = alphasimr_threads
  )

  summaries <- list()
  selections <- list()
  for (i in seq_len(nrow(methods))) {
    spec <- methods[i, , drop = FALSE]
    run <- ng_run_multitrait_validation(
      scores = scenario_obj$scores,
      traits = scenario_obj$traits,
      realized_cols = scenario_obj$realized_cols,
      n_crosses = n_crosses,
      seed = run_seed,
      methods = spec$method[[1]],
      allocator = spec$allocator[[1]],
      parent_K = scenario_obj$parent_K,
      ocs_lambda_group = ocs_lambda_group,
      ocs_lambda_mating = ocs_lambda_mating
    )
    s <- run$summary
    s$method <- spec$label[[1]]
    s$method_family <- ng_multitrait_method_family(spec$method[[1]])
    s$source_method <- spec$method[[1]]
    s$allocator <- spec$allocator[[1]]
    s$benchmark_role <- spec$role[[1]]
    s$external_tool <- spec$external_tool[[1]]
    s$implementation <- spec$implementation[[1]]
    s$exact_external_status <- spec$exact_external_status[[1]]
    s$fallback_reason <- spec$fallback_reason[[1]]
    s$scenario <- as.character(scenario_obj$config$scenario[[1]])
    s$crop <- as.character(scenario_obj$config$crop[[1]])
    s$harness_model <- as.character(scenario_obj$config$harness_model[[1]])
    s$n_parents <- n_parents
    s$rep <- rep
    s$seed <- run_seed
    s$realized_progeny <- realized_progeny
    summaries[[i]] <- s

    selected <- run$selections
    selected$method <- spec$label[[1]]
    selected$source_method <- spec$method[[1]]
    selected$allocator <- spec$allocator[[1]]
    selected$benchmark_role <- spec$role[[1]]
    selected$external_tool <- spec$external_tool[[1]]
    selected$implementation <- spec$implementation[[1]]
    selected$exact_external_status <- spec$exact_external_status[[1]]
    selected$fallback_reason <- spec$fallback_reason[[1]]
    selected$scenario <- as.character(scenario_obj$config$scenario[[1]])
    selected$crop <- as.character(scenario_obj$config$crop[[1]])
    selected$harness_model <- as.character(scenario_obj$config$harness_model[[1]])
    selected$n_parents <- n_parents
    selected$rep <- rep
    selected$seed <- run_seed
    selected$realized_progeny <- realized_progeny
    selections[[i]] <- selected
  }

  scores <- scenario_obj$scores
  scores$scenario <- as.character(scenario_obj$config$scenario[[1]])
  scores$crop <- as.character(scenario_obj$config$crop[[1]])
  scores$harness_model <- as.character(scenario_obj$config$harness_model[[1]])
  scores$n_parents <- n_parents
  scores$rep <- rep
  scores$seed <- run_seed
  scores$realized_progeny <- realized_progeny

  config <- scenario_obj$config
  config$rep <- rep
  config$n_crosses <- n_crosses
  config$methods <- paste(methods$label, collapse = ",")
  config$baseline_methods <- paste(methods$label[methods$role == "baseline"], collapse = ",")
  config$ocs_lambda_group <- ocs_lambda_group
  config$ocs_lambda_mating <- ocs_lambda_mating

  list(
    summary = ng_bind_rows_fill(summaries),
    selections = ng_bind_rows_fill(selections),
    scores = scores,
    config = config
  )
}

ng_run_head_to_head_benchmark <- function(scenarios = c("compact_selfing", "cassava_diploid", "potato_tetraploid_stress"),
                                          parent_sizes = c(12L, 20L),
                                          reps = 1L,
                                          n_crosses = 4L,
                                          realized_progeny = 8L,
                                          seed = 1L,
                                          methods = NULL,
                                          baseline_methods = NULL,
                                          n_founders = NULL,
                                          n_chr = NULL,
                                          seg_sites = NULL,
                                          snp_per_chr = NULL,
                                          qtl_per_chr = NULL,
                                          genome_length_m = NULL,
                                          prediction_noise = 0.20,
                                          alphasimr_threads = 1L,
                                          ocs_lambda_group = 0.05,
                                          ocs_lambda_mating = 0,
                                          comparison_metrics = NULL,
                                          output_dir = NULL,
                                          prefix = "head_to_head_benchmark") {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    ng_stop("AlphaSimR is required for head-to-head benchmarking")
  }
  methods <- ng_head_to_head_normalize_methods(methods)
  if (is.null(baseline_methods)) baseline_methods <- methods$label[methods$role == "baseline"]
  baseline_methods <- unique(trimws(as.character(baseline_methods)))
  baseline_methods <- baseline_methods[nzchar(baseline_methods) & !is.na(baseline_methods)]
  missing_baseline <- setdiff(baseline_methods, methods$label)
  if (length(missing_baseline)) ng_stop("baseline method(s) not present in method registry: ", paste(missing_baseline, collapse = ", "))
  if (!length(baseline_methods)) ng_stop("head-to-head benchmark needs at least one baseline method")

  crop_scenarios <- ng_crop_genome_select(scenarios)
  parent_sizes <- unique(suppressWarnings(as.integer(parent_sizes)))
  parent_sizes <- parent_sizes[is.finite(parent_sizes)]
  if (!length(parent_sizes)) ng_stop("parent_sizes must contain at least one size")
  if (any(parent_sizes < 4L)) ng_stop("parent_sizes must all be at least 4")
  reps <- ng_multitrait_crop_validation_int(reps, 1L, "reps")
  n_crosses <- ng_multitrait_crop_validation_int(n_crosses, 1L, "n_crosses")
  realized_progeny <- ng_multitrait_crop_validation_int(realized_progeny, 1L, "realized_progeny")
  seed <- ng_multitrait_crop_validation_int(seed, 1L, "seed")
  alphasimr_threads <- ng_multitrait_crop_validation_int(alphasimr_threads, 1L, "alphasimr_threads")
  prediction_noise <- ng_multitrait_crop_validation_num(prediction_noise, 0.20, "prediction_noise", min_value = 0)
  ocs_lambda_group <- ng_multitrait_crop_validation_num(ocs_lambda_group, 0.05, "ocs_lambda_group", min_value = 0)
  ocs_lambda_mating <- ng_multitrait_crop_validation_num(ocs_lambda_mating, 0, "ocs_lambda_mating", min_value = 0)

  cases <- list()
  k <- 0L
  for (s in seq_len(nrow(crop_scenarios))) {
    scenario_name <- as.character(crop_scenarios$scenario[[s]])
    for (n_parents in parent_sizes) {
      for (rep in seq_len(reps)) {
        run_seed <- ng_multitrait_crop_validation_seed(seed, scenario_name, n_parents, rep)
        k <- k + 1L
        cases[[k]] <- ng_head_to_head_run_case(
          scenario_name = scenario_name,
          n_parents = n_parents,
          rep = rep,
          run_seed = run_seed,
          methods = methods,
          n_crosses = n_crosses,
          realized_progeny = realized_progeny,
          n_founders = n_founders,
          n_chr = n_chr,
          seg_sites = seg_sites,
          snp_per_chr = snp_per_chr,
          qtl_per_chr = qtl_per_chr,
          genome_length_m = genome_length_m,
          prediction_noise = prediction_noise,
          alphasimr_threads = alphasimr_threads,
          ocs_lambda_group = ocs_lambda_group,
          ocs_lambda_mating = ocs_lambda_mating
        )
      }
    }
  }

  summary <- ng_bind_rows_fill(lapply(cases, `[[`, "summary"))
  selections <- ng_bind_rows_fill(lapply(cases, `[[`, "selections"))
  scores <- ng_bind_rows_fill(lapply(cases, `[[`, "scores"))
  config <- ng_bind_rows_fill(lapply(cases, `[[`, "config"))
  rownames(summary) <- NULL
  rownames(selections) <- NULL
  rownames(scores) <- NULL
  rownames(config) <- NULL

  comparisons <- ng_head_to_head_compare_summary(
    summary = summary,
    baseline_methods = baseline_methods,
    metrics = comparison_metrics
  )
  winner_summary <- ng_multitrait_crop_validation_winner_summary(summary)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
    write.csv(selections, file.path(output_dir, paste0(prefix, "_selections.csv")), row.names = FALSE)
    write.csv(scores, file.path(output_dir, paste0(prefix, "_scores.csv")), row.names = FALSE)
    write.csv(comparisons, file.path(output_dir, paste0(prefix, "_comparisons.csv")), row.names = FALSE)
    write.csv(winner_summary, file.path(output_dir, paste0(prefix, "_winner_summary.csv")), row.names = FALSE)
    write.csv(config, file.path(output_dir, paste0(prefix, "_config.csv")), row.names = FALSE)
    write.csv(methods, file.path(output_dir, paste0(prefix, "_method_registry.csv")), row.names = FALSE)
  }

  list(
    summary = summary,
    selections = selections,
    scores = scores,
    comparisons = comparisons,
    winner_summary = winner_summary,
    config = config,
    method_registry = methods,
    scenarios = crop_scenarios,
    parent_sizes = parent_sizes,
    reps = reps,
    n_crosses = n_crosses,
    realized_progeny = realized_progeny,
    seed = seed,
    baseline_methods = baseline_methods
  )
}
