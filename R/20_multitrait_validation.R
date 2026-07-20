ng_multitrait_validation_default_methods <- function() {
  c("auto", "weighted", "economic_index", "desired_gain", "threshold")
}

ng_multitrait_validation_traits <- function() {
  ng_multitrait_spec(
    trait = c("yield", "disease", "quality"),
    column = c("pred_yield", "pred_disease", "pred_quality"),
    direction = c("maximize", "minimize", "maximize"),
    weight = c(1, 4, 2),
    economic_weight = c(1, 4, 2),
    desired_change = c(4, 25, 2),
    max_value = c(NA, 40, NA),
    threshold_weight = c(1, 4, 1)
  )
}

ng_multitrait_validation_realized_cols <- function(traits) {
  traits <- ng_multitrait_spec(traits)
  stats::setNames(paste0("realized_", traits$trait), traits$trait)
}

ng_multitrait_validation_scenario <- function(n_parents = 20L,
                                              seed = 1L,
                                              prediction_noise = 0.25) {
  n_parents <- suppressWarnings(as.integer(n_parents[[1]]))
  if (!is.finite(n_parents) || n_parents < 4L) ng_stop("n_parents must be at least 4")
  prediction_noise <- suppressWarnings(as.numeric(prediction_noise[[1]]))
  if (!is.finite(prediction_noise) || prediction_noise < 0) prediction_noise <- 0.25
  set.seed(as.integer(seed[[1]]))

  parent_id <- sprintf("P%02d", seq_len(n_parents))
  latent <- sort(stats::rnorm(n_parents))
  parent_yield <- 100 + 8 * latent + stats::rnorm(n_parents, sd = 2.0)
  parent_disease <- 45 + 13 * latent + stats::rnorm(n_parents, sd = 4.0)
  parent_quality <- 35 - 2.5 * latent + stats::rnorm(n_parents, sd = 1.5)
  names(latent) <- names(parent_yield) <- names(parent_disease) <- names(parent_quality) <- parent_id

  pairs <- ng_make_pairs(parent_id, include_self = FALSE)
  p1 <- pairs$parent1
  p2 <- pairs$parent2
  latent_gap <- abs(latent[p1] - latent[p2])
  realized_yield <- 0.5 * (parent_yield[p1] + parent_yield[p2]) +
    1.5 * latent_gap + stats::rnorm(nrow(pairs), sd = 2.0)
  realized_disease <- 0.5 * (parent_disease[p1] + parent_disease[p2]) +
    0.8 * latent_gap + stats::rnorm(nrow(pairs), sd = 3.0)
  realized_quality <- 0.5 * (parent_quality[p1] + parent_quality[p2]) +
    0.4 * latent_gap + stats::rnorm(nrow(pairs), sd = 1.2)

  scores <- pairs
  scores$pred_yield <- realized_yield + stats::rnorm(nrow(pairs), sd = 2.5 + 4 * prediction_noise)
  scores$pred_disease <- realized_disease + stats::rnorm(nrow(pairs), sd = 3.5 + 5 * prediction_noise)
  scores$pred_quality <- realized_quality + stats::rnorm(nrow(pairs), sd = 1.5 + 2 * prediction_noise)
  scores$realized_yield <- realized_yield
  scores$realized_disease <- realized_disease
  scores$realized_quality <- realized_quality
  scores$pair_kinship <- 0.02 + 0.18 * exp(-latent_gap)

  traits <- ng_multitrait_validation_traits()
  list(
    scores = scores,
    traits = traits,
    realized_cols = ng_multitrait_validation_realized_cols(traits),
    parent_values = data.frame(
      parent = parent_id,
      latent = latent,
      yield = parent_yield,
      disease = parent_disease,
      quality = parent_quality,
      stringsAsFactors = FALSE
    )
  )
}

ng_multitrait_validation_realized_index <- function(scores,
                                                    traits,
                                                    realized_cols = NULL,
                                                    reference_scores = NULL) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (is.null(reference_scores)) reference_scores <- scores
  reference_scores <- as.data.frame(reference_scores, stringsAsFactors = FALSE)
  traits <- ng_multitrait_spec(traits)
  if (is.null(realized_cols)) realized_cols <- ng_multitrait_validation_realized_cols(traits)
  realized_cols <- realized_cols[traits$trait]
  missing <- setdiff(realized_cols, names(scores))
  if (length(missing)) ng_stop("scores missing realized columns: ", paste(missing, collapse = ", "))
  missing_ref <- setdiff(realized_cols, names(reference_scores))
  if (length(missing_ref)) ng_stop("reference_scores missing realized columns: ", paste(missing_ref, collapse = ", "))
  weights <- ng_multitrait_economic_weights(traits)
  z <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  colnames(z) <- traits$trait
  for (i in seq_len(nrow(traits))) {
    x <- suppressWarnings(as.numeric(scores[[realized_cols[[i]]]]))
    ref <- suppressWarnings(as.numeric(reference_scores[[realized_cols[[i]]]]))
    oriented <- if (identical(traits$direction[[i]], "maximize")) x else -x
    ref_oriented <- if (identical(traits$direction[[i]], "maximize")) ref else -ref
    center <- mean(ref_oriented, na.rm = TRUE)
    scale <- stats::sd(ref_oriented, na.rm = TRUE)
    if (!is.finite(center)) center <- 0
    if (!is.finite(scale) || scale <= 0) scale <- 1
    z[, i] <- (oriented - center) / scale
    z[!is.finite(z[, i]), i] <- 0
  }
  as.numeric(z %*% weights)
}

ng_multitrait_validation_coefficients <- function(plan_summary, traits) {
  traits <- ng_multitrait_spec(traits)
  empty <- stats::setNames(rep(NA_real_, nrow(traits)), paste0("coef_", traits$trait))
  coef <- NULL
  if (!is.null(plan_summary$multitrait_economic_index_coefficients)) {
    coef <- plan_summary$multitrait_economic_index_coefficients
  } else if (!is.null(plan_summary$multitrait_desired_gain_coefficients)) {
    coef <- plan_summary$multitrait_desired_gain_coefficients
  }
  if (is.null(coef)) return(empty)
  coef <- coef[traits$trait]
  out <- as.numeric(coef)
  names(out) <- paste0("coef_", traits$trait)
  out
}

ng_multitrait_validation_evaluate_plan <- function(plan,
                                                   method,
                                                   traits,
                                                   realized_cols = NULL,
                                                   reference_scores = NULL) {
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  traits <- ng_multitrait_spec(traits)
  if (is.null(realized_cols)) realized_cols <- ng_multitrait_validation_realized_cols(traits)
  realized_cols <- realized_cols[traits$trait]
  plan$realized_index <- ng_multitrait_validation_realized_index(
    plan, traits, realized_cols, reference_scores = reference_scores
  )
  parents <- sort(unique(c(plan$parent1, plan$parent2)))
  counts <- ng_parent_counts(plan, parents)
  plan_summary <- attr(plan, "summary")
  if (is.null(plan_summary)) plan_summary <- list()
  coef <- ng_multitrait_validation_coefficients(plan_summary, traits)

  out <- data.frame(
    method = method,
    method_family = ng_multitrait_method_family(method),
    selected_crosses = nrow(plan),
    mean_predicted_score = mean(plan$multi_trait_score, na.rm = TRUE),
    mean_realized_index = mean(plan$realized_index, na.rm = TRUE),
    mean_realized_yield = mean(plan[[realized_cols[["yield"]]]], na.rm = TRUE),
    mean_realized_disease = mean(plan[[realized_cols[["disease"]]]], na.rm = TRUE),
    mean_realized_quality = mean(plan[[realized_cols[["quality"]]]], na.rm = TRUE),
    unique_parents = sum(counts > 0),
    max_parent_use = max(counts),
    stringsAsFactors = FALSE
  )
  if (!is.null(plan_summary$mean_pair_kinship)) {
    out$mean_pair_kinship <- plan_summary$mean_pair_kinship
  } else if ("pair_kinship" %in% names(plan)) {
    out$mean_pair_kinship <- mean(plan$pair_kinship, na.rm = TRUE)
  }
  if (!is.null(plan_summary$group_coancestry)) out$group_coancestry <- plan_summary$group_coancestry
  if (!is.null(plan_summary$lambda_group)) out$lambda_group <- plan_summary$lambda_group
  if (!is.null(plan_summary$lambda_mating)) out$lambda_mating <- plan_summary$lambda_mating
  if (!is.null(plan_summary$lambda_parent_use)) out$lambda_parent_use <- plan_summary$lambda_parent_use
  if (!is.null(plan_summary$objective)) out$ocs_objective <- plan_summary$objective
  cbind(out, as.data.frame(as.list(coef), check.names = FALSE))
}

ng_run_multitrait_validation <- function(scores = NULL,
                                         traits = NULL,
                                         realized_cols = NULL,
                                         n_parents = 20L,
                                         n_crosses = 6L,
                                         seed = 1L,
                                         methods = ng_multitrait_validation_default_methods(),
                                         allocator = c("topn", "ocs"),
                                         parent_kinship = NULL,
                                         ocs_lambda_group = 0.05,
                                         ocs_lambda_mating = 0,
                                         output_dir = NULL,
                                         prefix = "multitrait_validation") {
  allocator <- match.arg(allocator)
  if (is.null(scores) || is.null(traits)) {
    scenario <- ng_multitrait_validation_scenario(n_parents = n_parents, seed = seed)
    if (is.null(scores)) scores <- scenario$scores
    if (is.null(traits)) traits <- scenario$traits
    if (is.null(realized_cols)) realized_cols <- scenario$realized_cols
  }
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  traits <- ng_multitrait_spec(traits)
  if (is.null(realized_cols)) realized_cols <- ng_multitrait_validation_realized_cols(traits)
  methods <- unique(trimws(tolower(as.character(methods))))
  methods <- methods[nzchar(methods) & !is.na(methods)]
  if (!length(methods)) ng_stop("methods must contain at least one method")
  n_crosses <- suppressWarnings(as.integer(n_crosses[[1]]))
  if (!is.finite(n_crosses) || n_crosses < 1L) ng_stop("n_crosses must be positive")

  parent_ids <- sort(unique(c(scores$parent1, scores$parent2)))
  if (is.null(parent_kinship)) {
    parent_kinship <- diag(length(parent_ids))
    rownames(parent_kinship) <- colnames(parent_kinship) <- parent_ids
  } else {
    parent_kinship <- as.matrix(parent_kinship)
    storage.mode(parent_kinship) <- "double"
    if (is.null(rownames(parent_kinship)) || is.null(colnames(parent_kinship))) {
      ng_stop("parent_kinship must have row and column names")
    }
    missing_k <- setdiff(parent_ids, intersect(rownames(parent_kinship), colnames(parent_kinship)))
    if (length(missing_k)) ng_stop("parent_kinship missing parents: ", paste(missing_k, collapse = ", "))
    parent_kinship <- parent_kinship[parent_ids, parent_ids, drop = FALSE]
  }

  selected <- vector("list", length(methods))
  summaries <- vector("list", length(methods))
  names(selected) <- names(summaries) <- methods
  for (method in methods) {
    plan <- if (identical(allocator, "topn")) {
      ng_multitrait_select_topn(scores, traits, n_crosses = n_crosses, method = method)
    } else {
      ng_optimize_multitrait_mating_plan(
        scores = scores,
        traits = traits,
        n_crosses = n_crosses,
        parent_kinship = parent_kinship,
        multitrait_method = method,
        optimizer_method = "greedy_local",
        max_crosses_per_parent = max(2L, ceiling(n_crosses / 2)),
        lambda_group = ocs_lambda_group,
        lambda_mating = ocs_lambda_mating
      )
    }
    plan <- as.data.frame(plan, stringsAsFactors = FALSE)
    plan$realized_index <- ng_multitrait_validation_realized_index(
      plan, traits, realized_cols, reference_scores = scores
    )
    plan$method <- method
    plan$selection_rank <- seq_len(nrow(plan))
    selected[[method]] <- plan
    summaries[[method]] <- ng_multitrait_validation_evaluate_plan(
      plan, method, traits, realized_cols, reference_scores = scores
    )
  }

  summary <- do.call(rbind, summaries)
  rownames(summary) <- NULL
  selections <- do.call(rbind, selected)
  rownames(selections) <- NULL

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
    write.csv(selections, file.path(output_dir, paste0(prefix, "_selections.csv")), row.names = FALSE)
    write.csv(scores, file.path(output_dir, paste0(prefix, "_scores.csv")), row.names = FALSE)
  }

  list(
    summary = summary,
    selections = selections,
    scores = scores,
    traits = traits,
    realized_cols = realized_cols,
    allocator = allocator
  )
}

ng_multitrait_validation_winner_summary <- function(summary,
                                                    metrics = NULL) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  required <- c("method", "n_parents", "rep")
  missing <- setdiff(required, names(summary))
  if (length(missing)) ng_stop("summary missing columns: ", paste(missing, collapse = ", "))

  default_metrics <- c(
    "mean_realized_index",
    "mean_realized_yield",
    "mean_realized_disease",
    "mean_realized_quality",
    "unique_parents",
    "max_parent_use"
  )
  directions <- c(
    mean_realized_index = "maximize",
    mean_realized_yield = "maximize",
    mean_realized_disease = "minimize",
    mean_realized_quality = "maximize",
    unique_parents = "maximize",
    max_parent_use = "minimize"
  )
  if (is.null(metrics)) metrics <- default_metrics
  metrics <- unique(trimws(as.character(metrics)))
  metrics <- metrics[nzchar(metrics) & !is.na(metrics)]
  missing_metrics <- setdiff(metrics, names(summary))
  if (length(missing_metrics)) ng_stop("summary missing metric columns: ", paste(missing_metrics, collapse = ", "))
  missing_direction <- setdiff(metrics, names(directions))
  if (length(missing_direction)) ng_stop("unknown metric direction: ", paste(missing_direction, collapse = ", "))

  n_parents <- sort(unique(suppressWarnings(as.integer(summary$n_parents))))
  n_parents <- n_parents[is.finite(n_parents)]
  out <- list()
  k <- 0L
  for (n in n_parents) {
    by_parent <- summary[summary$n_parents == n, , drop = FALSE]
    methods <- sort(unique(as.character(by_parent$method)))
    for (metric in metrics) {
      method_values <- stats::setNames(rep(NA_real_, length(methods)), methods)
      method_reps <- stats::setNames(rep(0L, length(methods)), methods)
      for (method in methods) {
        rows <- by_parent[by_parent$method == method, , drop = FALSE]
        values <- suppressWarnings(as.numeric(rows[[metric]]))
        method_values[[method]] <- mean(values, na.rm = TRUE)
        method_reps[[method]] <- length(unique(rows$rep[is.finite(suppressWarnings(as.integer(rows$rep)))]))
      }
      finite <- is.finite(method_values)
      if (!any(finite)) next
      best <- if (identical(directions[[metric]], "minimize")) {
        min(method_values[finite])
      } else {
        max(method_values[finite])
      }
      best_methods <- names(method_values)[finite & method_values == best]
      best_method <- sort(best_methods)[[1]]
      k <- k + 1L
      out[[k]] <- data.frame(
        n_parents = n,
        metric = metric,
        direction = directions[[metric]],
        method = best_method,
        tied_methods = paste(sort(best_methods), collapse = ","),
        tied_method_count = length(best_methods),
        value = unname(method_values[[best_method]]),
        reps = unname(method_reps[[best_method]]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) {
    return(data.frame(
      n_parents = integer(),
      metric = character(),
      direction = character(),
      method = character(),
      tied_methods = character(),
      tied_method_count = integer(),
      value = numeric(),
      reps = integer(),
      stringsAsFactors = FALSE
    ))
  }
  combined <- do.call(rbind, out)
  rownames(combined) <- NULL
  combined
}

ng_run_multitrait_validation_grid <- function(parent_sizes = c(20L, 40L, 60L, 80L),
                                              reps = 3L,
                                              n_crosses = 6L,
                                              seed = 1L,
                                              methods = ng_multitrait_validation_default_methods(),
                                              allocator = c("topn", "ocs"),
                                              ocs_lambda_group = 0.05,
                                              ocs_lambda_mating = 0,
                                              output_dir = NULL,
                                              prefix = "multitrait_grid") {
  allocator <- match.arg(allocator)
  parent_sizes <- unique(suppressWarnings(as.integer(parent_sizes)))
  parent_sizes <- parent_sizes[is.finite(parent_sizes)]
  if (!length(parent_sizes)) ng_stop("parent_sizes must contain at least one size")
  if (any(parent_sizes < 4L)) ng_stop("parent_sizes must all be at least 4")
  reps <- suppressWarnings(as.integer(reps[[1]]))
  if (!is.finite(reps) || reps < 1L) ng_stop("reps must be positive")
  seed <- suppressWarnings(as.integer(seed[[1]]))
  if (!is.finite(seed)) seed <- 1L

  summaries <- list()
  selections <- list()
  scores_out <- list()
  configs <- list()
  k <- 0L
  for (n_parents in parent_sizes) {
    for (rep in seq_len(reps)) {
      run_seed <- seed + (as.integer(n_parents) * 1000L) + rep
      run <- ng_run_multitrait_validation(
        n_parents = n_parents,
        n_crosses = n_crosses,
        seed = run_seed,
        methods = methods,
        allocator = allocator,
        ocs_lambda_group = ocs_lambda_group,
        ocs_lambda_mating = ocs_lambda_mating
      )
      k <- k + 1L
      run$summary$n_parents <- n_parents
      run$summary$rep <- rep
      run$summary$seed <- run_seed
      run$summary$allocator <- allocator
      run$selections$n_parents <- n_parents
      run$selections$rep <- rep
      run$selections$seed <- run_seed
      run$selections$allocator <- allocator
      summaries[[k]] <- run$summary
      selections[[k]] <- run$selections
      score_rows <- run$scores
      score_rows$n_parents <- n_parents
      score_rows$rep <- rep
      score_rows$seed <- run_seed
      score_rows$allocator <- allocator
      scores_out[[k]] <- score_rows
      configs[[k]] <- data.frame(
        n_parents = n_parents,
        rep = rep,
        seed = run_seed,
        allocator = allocator,
        methods = paste(methods, collapse = ","),
        n_crosses = n_crosses,
        ocs_lambda_group = ocs_lambda_group,
        ocs_lambda_mating = ocs_lambda_mating,
        stringsAsFactors = FALSE
      )
    }
  }

  summary <- do.call(rbind, summaries)
  selections <- do.call(rbind, selections)
  scores <- do.call(rbind, scores_out)
  config <- do.call(rbind, configs)
  rownames(summary) <- NULL
  rownames(selections) <- NULL
  rownames(scores) <- NULL
  rownames(config) <- NULL
  winner_summary <- ng_multitrait_validation_winner_summary(summary)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
    write.csv(selections, file.path(output_dir, paste0(prefix, "_selections.csv")), row.names = FALSE)
    write.csv(scores, file.path(output_dir, paste0(prefix, "_scores.csv")), row.names = FALSE)
    write.csv(winner_summary, file.path(output_dir, paste0(prefix, "_winner_summary.csv")), row.names = FALSE)
    write.csv(config, file.path(output_dir, paste0(prefix, "_config.csv")), row.names = FALSE)
  }

  list(
    summary = summary,
    selections = selections,
    scores = scores,
    winner_summary = winner_summary,
    config = config,
    parent_sizes = parent_sizes,
    reps = reps,
    n_crosses = n_crosses,
    seed = seed,
    allocator = allocator
  )
}
