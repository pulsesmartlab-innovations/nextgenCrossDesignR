ng_first_numeric_col <- function(data, candidates, object_name = "data") {
  hit <- candidates[candidates %in% names(data)]
  for (col in hit) {
    x <- suppressWarnings(as.numeric(data[[col]]))
    if (any(is.finite(x))) return(col)
  }
  NULL
}

ng_criterion_vector <- function(criterion,
                                ids = NULL,
                                id_col = NULL,
                                criterion_col = NULL,
                                object_name = "criterion") {
  if (is.null(criterion)) {
    if (is.null(ids)) ng_stop(object_name, " is required")
    out <- stats::setNames(rep(0, length(ids)), ids)
    return(out)
  }
  if (is.data.frame(criterion)) {
    if (is.null(id_col)) {
      id_col <- ng_first_available_col(criterion, c("Genotype", "Id", "ID", "id", "parent"), character())
    }
    if (is.null(criterion_col)) {
      criterion_col <- ng_first_numeric_col(criterion, c("Criterion", "criterion", "Y", "score", "gebv", "value"))
    }
    if (is.null(id_col) || !(id_col %in% names(criterion))) {
      ng_stop(object_name, " data frame must contain an ID column")
    }
    if (is.null(criterion_col) || !(criterion_col %in% names(criterion))) {
      ng_stop(object_name, " data frame must contain a numeric criterion column")
    }
    out <- suppressWarnings(as.numeric(criterion[[criterion_col]]))
    names(out) <- as.character(criterion[[id_col]])
  } else {
    out <- suppressWarnings(as.numeric(criterion))
    if (is.null(names(out))) {
      if (is.null(ids) || length(out) != length(ids)) {
        ng_stop(object_name, " vector must be named or have the same length as ids")
      }
      names(out) <- ids
    }
  }
  out <- out[!is.na(names(out)) & nzchar(names(out))]
  if (any(duplicated(names(out)))) {
    keep <- !duplicated(names(out))
    out <- out[keep]
  }
  if (!is.null(ids)) {
    miss <- setdiff(ids, names(out))
    if (length(miss)) ng_stop(object_name, " is missing ", length(miss), " id(s)")
    out <- out[ids]
  }
  out
}

ng_popvar_style_scores <- function(scores,
                                   selection_prop = 0.10,
                                   mean_col = NULL,
                                   var_col = NULL,
                                   out_prefix = "popvar_style") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  required <- c("parent1", "parent2")
  missing <- setdiff(required, names(scores))
  if (length(missing)) ng_stop("scores missing columns: ", paste(missing, collapse = ", "))
  if (is.null(mean_col)) {
    mean_col <- ng_first_numeric_col(
      scores,
      c("cross_mean_blend", "cross_mean_gebv", "cross_mean", "multi_trait_score", "trait_index"),
      "scores"
    )
  }
  if (is.null(var_col)) {
    var_col <- ng_first_numeric_col(
      scores,
      c("pmv_cal", "vpm_cal",
        "pmv", "vpm"),
      "scores"
    )
  }
  if (is.null(mean_col)) ng_stop("No usable mean column found for PopVar-style scoring")
  mu <- suppressWarnings(as.numeric(scores[[mean_col]]))
  var_g <- if (is.null(var_col)) rep(0, nrow(scores)) else suppressWarnings(as.numeric(scores[[var_col]]))
  var_g[!is.finite(var_g)] <- NA_real_
  var_g <- pmax(var_g, 0, na.rm = FALSE)
  intensity <- ng_selection_intensity(selection_prop)
  uc <- mu + intensity * sqrt(pmax(var_g, 0))

  scores[[paste0(out_prefix, "_mu")]] <- mu
  scores[[paste0(out_prefix, "_varG")]] <- var_g
  scores[[paste0(out_prefix, "_uc")]] <- uc
  scores[[paste0(out_prefix, "_status")]] <- ifelse(is.finite(uc), "native_proxy", "unavailable")
  attr(scores, paste0(out_prefix, "_metadata")) <- list(
    implementation = "native_proxy",
    mean_col = mean_col,
    var_col = if (is.null(var_col)) NA_character_ else var_col,
    selection_prop = selection_prop
  )
  scores
}

ng_apply_popvar_native_proxy <- function(scores,
                                         tail_p = 0.10,
                                         mean_col = NULL,
                                         var_col = NULL,
                                         status = "native_proxy") {
  prox <- ng_popvar_style_scores(
    scores = scores,
    selection_prop = tail_p,
    mean_col = mean_col,
    var_col = var_col,
    out_prefix = ".ng_popvar_proxy"
  )
  scores$popvar_mu <- prox$.ng_popvar_proxy_mu
  scores$popvar_varG <- prox$.ng_popvar_proxy_varG
  scores$popvar_musp_low <- scores$popvar_mu
  scores$popvar_musp_high <- scores$popvar_uc <- prox$.ng_popvar_proxy_uc
  scores$popvar_status <- ifelse(is.finite(scores$popvar_uc), status, "unavailable")
  attr(scores, "popvar_proxy_metadata") <- attr(prox, ".ng_popvar_proxy_metadata")
  scores
}

ng_apply_simplemating_native_proxy <- function(scores,
                                               prop_sel = 0.10,
                                               mean_col = NULL,
                                               var_col = NULL,
                                               status = "native_proxy") {
  prox <- ng_popvar_style_scores(
    scores = scores,
    selection_prop = prop_sel,
    mean_col = mean_col,
    var_col = var_col,
    out_prefix = ".ng_simple_proxy"
  )
  scores$simple_mpv <- prox$.ng_simple_proxy_mu
  scores$simple_usefa_mean <- prox$.ng_simple_proxy_mu
  scores$simple_usefa_var <- prox$.ng_simple_proxy_varG
  scores$simple_usefa_sd <- sqrt(pmax(scores$simple_usefa_var, 0))
  scores$simple_usefa <- prox$.ng_simple_proxy_uc
  scores$simple_status <- ifelse(is.finite(scores$simple_usefa), status, "unavailable")
  attr(scores, "simplemating_proxy_metadata") <- attr(prox, ".ng_simple_proxy_metadata")
  scores
}

ng_relationship_components <- function(parent_kinship, threshold = 0.5) {
  K <- as.matrix(parent_kinship)
  storage.mode(K) <- "double"
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    ng_stop("parent_kinship must have row and column names")
  }
  ids <- rownames(K)
  K <- K[ids, ids, drop = FALSE]
  linked <- K >= threshold
  linked[!is.finite(K)] <- FALSE
  diag(linked) <- TRUE
  seen <- stats::setNames(rep(FALSE, length(ids)), ids)
  comps <- list()
  while (any(!seen)) {
    start <- ids[which(!seen)[[1L]]]
    queue <- start
    seen[start] <- TRUE
    members <- character()
    while (length(queue)) {
      current <- queue[[1L]]
      queue <- queue[-1L]
      members <- c(members, current)
      next_ids <- ids[linked[current, ids] & !seen[ids]]
      if (length(next_ids)) {
        seen[next_ids] <- TRUE
        queue <- c(queue, next_ids)
      }
    }
    comps[[length(comps) + 1L]] <- members
  }
  comps
}

ng_simplemating_relate_thinning <- function(parent_kinship,
                                            criterion,
                                            threshold = 0.5,
                                            max_per_cluster = 2L,
                                            id_col = NULL,
                                            criterion_col = NULL) {
  if (!is.finite(threshold)) ng_stop("threshold must be finite")
  max_per_cluster <- as.integer(max_per_cluster)
  if (!is.finite(max_per_cluster) || max_per_cluster < 1L) {
    ng_stop("max_per_cluster must be a positive integer")
  }
  K <- as.matrix(parent_kinship)
  ids <- rownames(K)
  if (is.null(ids)) ng_stop("parent_kinship must have row names")
  crit <- ng_criterion_vector(
    criterion,
    ids = ids,
    id_col = id_col,
    criterion_col = criterion_col,
    object_name = "criterion"
  )
  comps <- ng_relationship_components(K, threshold = threshold)
  keep <- unlist(lapply(comps, function(members) {
    values <- crit[members]
    members[order(values, members, decreasing = c(TRUE, FALSE), method = "radix")][seq_len(min(max_per_cluster, length(members)))]
  }), use.names = FALSE)
  keep <- unique(keep)
  attr(keep, "summary") <- data.frame(
    n_input = length(ids),
    n_kept = length(keep),
    threshold = threshold,
    max_per_cluster = max_per_cluster,
    n_clusters = length(comps),
    stringsAsFactors = FALSE
  )
  keep
}

ng_unordered_pair_key <- function(parent1, parent2) {
  p1 <- as.character(parent1)
  p2 <- as.character(parent2)
  lo <- ifelse(p1 <= p2, p1, p2)
  hi <- ifelse(p1 <= p2, p2, p1)
  paste(lo, hi, sep = "\r")
}

ng_simplemating_past_thinning <- function(candidate_pairs,
                                          past_crosses,
                                          parent1_col = "parent1",
                                          parent2_col = "parent2",
                                          past_parent1_col = NULL,
                                          past_parent2_col = NULL,
                                          reciprocal = TRUE) {
  pairs <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE)
  past <- as.data.frame(past_crosses, stringsAsFactors = FALSE)
  if (!(parent1_col %in% names(pairs)) && "Parent1" %in% names(pairs)) parent1_col <- "Parent1"
  if (!(parent2_col %in% names(pairs)) && "Parent2" %in% names(pairs)) parent2_col <- "Parent2"
  missing <- setdiff(c(parent1_col, parent2_col), names(pairs))
  if (length(missing)) ng_stop("candidate_pairs missing columns: ", paste(missing, collapse = ", "))
  if (is.null(past_parent1_col)) past_parent1_col <- if ("parent1" %in% names(past)) "parent1" else names(past)[[1L]]
  if (is.null(past_parent2_col)) past_parent2_col <- if ("parent2" %in% names(past)) "parent2" else names(past)[[2L]]
  missing_past <- setdiff(c(past_parent1_col, past_parent2_col), names(past))
  if (length(missing_past)) ng_stop("past_crosses missing columns: ", paste(missing_past, collapse = ", "))

  pair_key <- if (isTRUE(reciprocal)) {
    ng_unordered_pair_key(pairs[[parent1_col]], pairs[[parent2_col]])
  } else {
    ng_pair_key(pairs[[parent1_col]], pairs[[parent2_col]])
  }
  past_key <- if (isTRUE(reciprocal)) {
    ng_unordered_pair_key(past[[past_parent1_col]], past[[past_parent2_col]])
  } else {
    ng_pair_key(past[[past_parent1_col]], past[[past_parent2_col]])
  }
  keep <- !(pair_key %in% past_key)
  out <- pairs[keep, , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "n_removed") <- sum(!keep)
  out
}

ng_simplemating_build_crosses <- function(moms,
                                          dads,
                                          keep = NULL,
                                          criterion = NULL,
                                          parent_kinship = NULL,
                                          include_self = FALSE,
                                          max_pair_kinship = Inf,
                                          id_col = NULL,
                                          criterion_col = NULL) {
  moms <- unique(as.character(moms))
  dads <- unique(as.character(dads))
  if (!length(moms) || !length(dads)) ng_stop("moms and dads must be non-empty")
  if (!is.null(keep)) {
    keep <- unique(as.character(keep))
    moms <- intersect(moms, keep)
    dads <- intersect(dads, keep)
  }
  if (!length(moms) || !length(dads)) ng_stop("No moms/dads remain after keep filtering")
  ids <- sort(unique(c(moms, dads)))
  crit <- ng_criterion_vector(
    criterion,
    ids = ids,
    id_col = id_col,
    criterion_col = criterion_col,
    object_name = "criterion"
  )
  grid <- expand.grid(parent1 = moms, parent2 = dads, stringsAsFactors = FALSE)
  if (!isTRUE(include_self)) {
    grid <- grid[grid$parent1 != grid$parent2, , drop = FALSE]
  }
  if (!nrow(grid)) ng_stop("No feasible crosses after self-cross filtering")
  key <- ng_unordered_pair_key(grid$parent1, grid$parent2)
  grid <- grid[!duplicated(key), , drop = FALSE]
  rownames(grid) <- NULL
  grid$criterion_mean <- 0.5 * (crit[grid$parent1] + crit[grid$parent2])
  grid$Y <- grid$criterion_mean
  if (!is.null(parent_kinship)) {
    K <- as.matrix(parent_kinship)
    if (is.null(rownames(K)) || is.null(colnames(K))) ng_stop("parent_kinship must have row and column names")
    missing <- setdiff(unique(c(grid$parent1, grid$parent2)), rownames(K))
    if (length(missing)) ng_stop("parent_kinship missing parent IDs: ", paste(missing, collapse = ", "))
    grid$pair_relationship <- as.numeric(K[cbind(grid$parent1, grid$parent2)])
    grid$pair_kinship <- grid$pair_relationship / 2
  } else {
    grid$pair_relationship <- 0
    grid$pair_kinship <- 0
  }
  # SimpleMating-compatible K remains on relationship scale.
  grid$K <- grid$pair_relationship
  if (is.finite(max_pair_kinship)) {
    grid <- grid[is.finite(grid$pair_kinship) & grid$pair_kinship <= max_pair_kinship, , drop = FALSE]
    rownames(grid) <- NULL
  }
  grid
}

ng_simplemating_style_select <- function(scores,
                                                  score_col = NULL,
                                                  n_crosses,
                                                  parent_kinship = NULL,
                                                  max_crosses_per_parent = 4L,
                                                  min_unique_parents = NULL,
                                                  culling_pairwise_k = NULL,
                                                  method = "auto",
                                                  local_iter = 2000L) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (is.null(score_col)) {
    score_col <- ng_first_numeric_col(scores, c("Y", "criterion_mean", "simple_usefa", "usefulness_pmv_gebv", "usefulness_vpm_gebv", "cross_mean"), "scores")
  }
  if (is.null(score_col) || !(score_col %in% names(scores))) ng_stop("No usable score_col found")
  if (!("pair_kinship" %in% names(scores))) {
    if ("K" %in% names(scores)) {
      scores$pair_relationship <- suppressWarnings(as.numeric(scores$K))
      scores$pair_kinship <- scores$pair_relationship / 2
    } else {
      scores$pair_relationship <- 0
      scores$pair_kinship <- 0
    }
  }
  # SimpleMating's K threshold is relationship-scale; the native optimizer's
  # max_pair_kinship is coancestry-scale.
  max_pair_kinship <- if (is.null(culling_pairwise_k)) Inf else as.numeric(culling_pairwise_k) / 2
  plan <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = n_crosses,
    gain_col = score_col,
    parent_kinship = parent_kinship,
    max_crosses_per_parent = max_crosses_per_parent,
    min_unique_parents = min_unique_parents,
    max_pair_kinship = max_pair_kinship,
    lambda_group = 0,
    lambda_parent_use = 0,
    method = method,
    local_iter = local_iter
  )
  s <- attr(plan, "summary")
  s$simplemating_style <- "native_proxy"
  s$simplemating_score_col <- score_col
  s$simplemating_culling_pairwise_k <- if (is.null(culling_pairwise_k)) Inf else as.numeric(culling_pairwise_k)
  s$style_proxy <- TRUE
  s$style_proxy_note <- "Lambda-penalized OCS approximation of SimpleMating selectCrosses; not the exact SimpleMating algorithm."
  attr(plan, "summary") <- s
  attr(plan, "style_proxy") <- TRUE
  plan
}

# Back-compat alias: emits a deprecation warning so downstream code can update.
ng_simplemating_select_crosses_native <- function(...) {
  .Deprecated("ng_simplemating_style_select")
  ng_simplemating_style_select(...)
}

ng_alphamate_style_lambda_grid <- function(lambda_group = NULL) {
  if (!is.null(lambda_group)) return(as.numeric(lambda_group))
  unique(c(0, 10 ^ seq(-4, 4, length.out = 17)))
}

ng_alphamate_style_frontier <- function(scores,
                                        criterion_col,
                                        n_crosses,
                                        parent_kinship,
                                        max_contributions = NULL,
                                        lambda_grid = NULL,
                                        method = "auto",
                                        local_iter = 2000L) {
  if (is.null(lambda_grid)) lambda_grid <- ng_alphamate_style_lambda_grid()
  lambda_grid <- unique(as.numeric(lambda_grid))
  lambda_grid <- lambda_grid[is.finite(lambda_grid) & lambda_grid >= 0]
  if (!length(lambda_grid)) ng_stop("lambda_grid must contain at least one non-negative finite value")

  plans <- list()
  rows <- list()
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[[i]]
    plan <- tryCatch(
      ng_optimize_mating_plan(
        scores = scores,
        n_crosses = n_crosses,
        gain_col = criterion_col,
        parent_kinship = parent_kinship,
        max_crosses_per_parent = if (is.null(max_contributions)) n_crosses else max_contributions,
        lambda_group = lambda,
        lambda_parent_use = 0,
        method = method,
        local_iter = local_iter
      ),
      error = function(e) e
    )
    if (inherits(plan, "error")) next
    s <- attr(plan, "summary")
    key <- paste(sort(ng_unordered_pair_key(plan$parent1, plan$parent2)), collapse = "|")
    rows[[length(rows) + 1L]] <- data.frame(
      plan_index = length(plans) + 1L,
      lambda_group = lambda,
      plan_key = key,
      mean_gain = s$mean_gain,
      group_relationship = s$group_relationship,
      group_coancestry = s$group_coancestry,
      mean_pair_kinship = s$mean_pair_kinship,
      unique_parents = s$unique_parents,
      max_parent_use = s$max_parent_use,
      stringsAsFactors = FALSE
    )
    plans[[length(plans) + 1L]] <- plan
  }
  if (!length(plans)) ng_stop("No feasible AlphaMate-style native plans across lambda grid")
  frontier <- do.call(rbind, rows)
  keep <- !duplicated(frontier$plan_key)
  frontier <- frontier[keep, , drop = FALSE]
  plans <- plans[frontier$plan_index]
  frontier$plan_index <- seq_len(nrow(frontier))
  list(frontier = frontier, plans = plans)
}

# Map a diversity-emphasis value to the row index of a gain-vs-coancestry frontier
# data.frame. `emphasis` runs 0 = maximum genetic gain .. 100 = minimum parental
# coancestry (maximum diversity). The frontier must carry numeric `mean_gain` and
# `group_coancestry` columns, as produced by both ng_pareto_mate_allocation() (native)
# and ng_alphamate_style_frontier() (proxy). mode = "max_gain"/"min_coancestry" pick the
# extreme points directly; mode = "target" scalarizes normalized gain vs diversity at the
# requested emphasis. Shared by ng_select_by_strategy() and the AlphaMate proxy so the two
# agree exactly on how an emphasis maps to a frontier point.
ng_choose_frontier_point <- function(frontier,
                                     emphasis = 45,
                                     mode = c("target", "max_gain", "min_coancestry")) {
  mode <- match.arg(mode)
  frontier <- as.data.frame(frontier, stringsAsFactors = FALSE)
  if (!nrow(frontier)) ng_stop("frontier is empty")
  gain <- suppressWarnings(as.numeric(frontier$mean_gain))
  coanc <- suppressWarnings(as.numeric(frontier$group_coancestry))
  if (identical(mode, "max_gain")) {
    return(which.max(gain))
  }
  if (identical(mode, "min_coancestry")) {
    return(order(coanc, -gain, na.last = TRUE)[[1L]])
  }
  target <- suppressWarnings(as.numeric(emphasis))
  if (!is.finite(target)) target <- 45
  target <- max(0, min(100, target)) / 100
  gain_norm <- ng_scale01(gain, bigger_is_better = TRUE)
  diversity_norm <- ng_scale01(coanc, bigger_is_better = FALSE)
  score <- (1 - target) * gain_norm + target * diversity_norm
  order(score, gain, -coanc, decreasing = TRUE, na.last = TRUE)[[1L]]
}

ng_alphamate_style_choose_frontier_plan <- function(frontier, mode, target_degree) {
  if (!nrow(frontier)) ng_stop("AlphaMate-style frontier is empty")
  emphasis_mode <- switch(mode,
    ModeMaxCriterion = "max_gain",
    ModeMinCoancestry = "min_coancestry",
    "target")
  ng_choose_frontier_point(frontier, emphasis = target_degree, mode = emphasis_mode)
}

ng_alphamate_style_select <- function(scores,
                                      criterion_col = "cross_mean",
                                      n_crosses,
                                      parent_kinship,
                                      mode = c("ModeOptTarget1", "ModeMaxCriterion", "ModeMinCoancestry"),
                                      target_degree = 45,
                                      max_contributions = NULL,
                                      lambda_group = NULL,
                                      lambda_grid = NULL,
                                      method = "auto",
                                      local_iter = 2000L) {
  mode <- match.arg(mode)
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!(criterion_col %in% names(scores))) ng_stop("scores missing ", criterion_col)
  if (!is.null(lambda_group)) lambda_grid <- ng_alphamate_style_lambda_grid(lambda_group)
  frontier <- ng_alphamate_style_frontier(
    scores = scores,
    criterion_col = criterion_col,
    n_crosses = n_crosses,
    parent_kinship = parent_kinship,
    max_contributions = max_contributions,
    lambda_grid = lambda_grid,
    method = method,
    local_iter = local_iter
  )
  selected_idx <- ng_alphamate_style_choose_frontier_plan(
    frontier = frontier$frontier,
    mode = mode,
    target_degree = target_degree
  )
  plan <- frontier$plans[[selected_idx]]
  s <- attr(plan, "summary")
  s$alphamate_style <- "native_proxy"
  s$alphamate_mode <- mode
  s$alphamate_target_degree <- as.numeric(target_degree)
  s$alphamate_criterion_col <- criterion_col
  s$alphamate_selection_rule <- if (is.null(lambda_group)) "target_frontier" else "fixed_lambda"
  s$alphamate_frontier_n <- nrow(frontier$frontier)
  s$alphamate_selected_lambda_group <- frontier$frontier$lambda_group[[selected_idx]]
  # Approximate "achieved degree" from the selected plan's gain position on the
  # frontier: degree = 100 * (1 - gain_norm), where gain_norm in [0, 1] is the
  # normalized mean_gain across the frontier. Because the proxy scalarizes
  # gain vs coancestry rather than solving the constraint exactly, this is a
  # diagnostic, not a guarantee.
  gain_vec <- suppressWarnings(as.numeric(frontier$frontier$mean_gain))
  gain_rng <- range(gain_vec, na.rm = TRUE, finite = TRUE)
  achieved_degree <- if (all(is.finite(gain_rng)) && diff(gain_rng) > 0) {
    100 * (1 - (gain_vec[[selected_idx]] - gain_rng[[1L]]) / diff(gain_rng))
  } else {
    NA_real_
  }
  s$alphamate_achieved_degree <- achieved_degree
  s$alphamate_degree_gap <- if (is.finite(achieved_degree)) {
    as.numeric(target_degree) - achieved_degree
  } else {
    NA_real_
  }
  s$alphamate_frontier <- frontier$frontier
  s$style_proxy <- TRUE
  s$style_proxy_note <- "Scalarized lambda-grid approximation of the AlphaMate gain-diversity frontier; the achieved coancestry degree is NOT guaranteed to equal target_degree."
  attr(plan, "summary") <- s
  attr(plan, "style_proxy") <- TRUE
  plan$alphamate_mode <- mode
  plan$alphamate_target_degree <- as.numeric(target_degree)
  plan
}
