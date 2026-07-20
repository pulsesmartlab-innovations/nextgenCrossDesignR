ng_optimize_mating_plan <- function(scores,
                                    n_crosses,
                                    gain_col = "usefulness_pmv_gebv",
                                    parent_kinship = NULL,
                                    max_crosses_per_parent = 4,
                                    min_crosses_per_parent = 0,
                                    min_unique_parents = NULL,
                                    max_pair_kinship = Inf,
                                    committed_crosses = NULL,
                                    parent_group = NULL,
                                    group_permission = NULL,
                                    group_quota = NULL,
                                    cost_col = NULL,
                                    budget = Inf,
                                    lambda_cost = 0,
                                    logistic_col = NULL,
                                    lambda_logistic = 0,
                                    lambda_group = 0,
                                    lambda_mating = 0,
                                    lambda_progeny_inbreeding = 0,
                                    lambda_parent_use = 0,
                                    lambda_parent_use_mode = c("absolute", "adaptive"),
                                    method = c("auto", "greedy_local", "repair_local", "mip_linear", "mip_contribution", "evolution"),
                                    strategy = NULL,
                                    diversity_emphasis = NULL,
                                    target_coancestry = NULL,
                                    local_iter = 2000,
                                    ocs_iter = 5L,
                                    mip_time_limit = 10,
                                    mip_max_binary_vars = 2e5,
                                    evol_solutions = 100L,
                                    evol_iterations = 200L,
                                    evol_stop = 40L,
                                    evol_seed = NULL) {
  method <- match.arg(method)
  lambda_parent_use_mode <- match.arg(lambda_parent_use_mode)
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  required <- c("parent1", "parent2", gain_col)
  miss <- setdiff(required, names(scores))
  if (length(miss)) ng_stop("scores missing columns: ", paste(miss, collapse = ", "))
  if (!("pair_kinship" %in% names(scores))) scores$pair_kinship <- 0
  scores <- scores[is.finite(scores[[gain_col]]) & scores$pair_kinship <= max_pair_kinship, , drop = FALSE]
  # Mating-group legality: drop candidate crosses whose parents' groups are not an
  # allowed group x group pairing, so every downstream allocator only sees legal matings.
  scores <- ng_filter_group_permission(scores, parent_group, group_permission)
  if (nrow(scores) < n_crosses) ng_stop("Not enough feasible candidate crosses")
  parents <- sort(unique(c(scores$parent1, scores$parent2)))
  if (is.null(parent_kinship)) {
    parent_kinship <- diag(length(parents))
    rownames(parent_kinship) <- colnames(parent_kinship) <- parents
  } else {
    parent_kinship <- parent_kinship[parents, parents, drop = FALSE]
  }
  # Gain-diversity balancing dispatch: when the caller asks for a strategy label or a
  # diversity_emphasis instead of a raw lambda_group, hand off to the frontier-navigation
  # layer, which sweeps the native frontier and picks/refines the lambda_group that lands
  # at the requested gain-vs-coancestry balance. The inner sweep calls back into this
  # function with lambda_group set and strategy = NULL, so there is no unbounded recursion.
  # target_coancestry is an ABSOLUTE constraint (constrained OCS) and takes precedence over
  # the relative emphasis dial; both pick a point on the same frontier (no new penalty).
  frontier_forward <- list(
    max_crosses_per_parent = max_crosses_per_parent,
    min_crosses_per_parent = min_crosses_per_parent,
    min_unique_parents = min_unique_parents, max_pair_kinship = max_pair_kinship,
    committed_crosses = committed_crosses, parent_group = parent_group,
    group_permission = group_permission, group_quota = group_quota,
    cost_col = cost_col, budget = budget, lambda_cost = lambda_cost,
    logistic_col = logistic_col, lambda_logistic = lambda_logistic,
    lambda_mating = lambda_mating, lambda_progeny_inbreeding = lambda_progeny_inbreeding,
    lambda_parent_use = lambda_parent_use,
    lambda_parent_use_mode = lambda_parent_use_mode, local_iter = local_iter)
  if (!is.null(target_coancestry)) {
    return(do.call(ng_select_by_target_coancestry, c(list(
      scores = scores, n_crosses = n_crosses, parent_kinship = parent_kinship, gain_col = gain_col,
      target_coancestry = target_coancestry,
      method = if (identical(method, "auto")) "greedy_local" else method), frontier_forward)))
  }
  if (!is.null(strategy) || !is.null(diversity_emphasis)) {
    return(do.call(ng_select_by_strategy, c(list(
      scores = scores, n_crosses = n_crosses, parent_kinship = parent_kinship, gain_col = gain_col,
      strategy = strategy, diversity_emphasis = diversity_emphasis,
      method = if (identical(method, "auto")) "greedy_local" else method), frontier_forward)))
  }
  # .linear_gain folds the per-cross penalties into one value every allocator maximizes;
  # the group-coancestry (lambda_group) and parent-use penalties act on the whole plan.
  #
  # IMPORTANT — the two per-cross relatedness penalties act on the SAME axis (parent-pair
  # relatedness = immediate progeny inbreeding) and ADD; they are not independent:
  #   * lambda_mating * pair_kinship            : two-sided linear on the VanRaden-G
  #       relationship; penalizes related pairs AND rewards complementary (negative-G) pairs.
  #   * lambda_progeny_inbreeding * expected_progeny_inbreeding, where
  #       expected_progeny_inbreeding = max(0, pair_kinship / 2) : one-sided inbreeding
  #       penalty (an inbreeding coefficient is >= 0 by definition), neutral on
  #       complementary pairs.
  # For related pairs (pair_kinship > 0) the two are perfectly collinear (differ only by
  # the 1/2 scale), so setting BOTH double-emphasizes relatedness. Prefer ONE: lambda_mating
  # if you also want to actively favor complementary matings, lambda_progeny_inbreeding for
  # a pure inbreeding-avoidance penalty in interpretable coancestry units. lambda_group is
  # the separate, population-level (plan-wide) diversity/ΔF control.
  epi <- if ("expected_progeny_inbreeding" %in% names(scores)) {
    as.numeric(scores$expected_progeny_inbreeding)
  } else {
    pmax(0, as.numeric(scores$pair_kinship) / 2)
  }
  epi[!is.finite(epi)] <- 0
  if (!is.finite(lambda_progeny_inbreeding) || lambda_progeny_inbreeding < 0) lambda_progeny_inbreeding <- 0
  if (lambda_mating > 0 && lambda_progeny_inbreeding > 0 &&
      !isTRUE(getOption("ngcd.warned_relatedness_overlap"))) {
    warning("lambda_mating and lambda_progeny_inbreeding both penalize parent-pair ",
            "relatedness (immediate progeny inbreeding) and their effects add; prefer ",
            "setting one. lambda_group controls population-level coancestry separately.",
            call. = FALSE)
    options(ngcd.warned_relatedness_overlap = TRUE)
  }
  scores$.linear_gain <- scores[[gain_col]] - lambda_mating * scores$pair_kinship -
    lambda_progeny_inbreeding * epi
  # Optional per-cross logistic / operational factors folded into the same per-cross
  # value: a soft cost penalty (lambda_cost on cost_col) and a soft logistic/geographic
  # penalty (lambda_logistic on logistic_col, e.g. distance or reproductive difficulty).
  # cost_col also drives the hard `budget` cap applied to the finished plan below.
  if (!is.null(cost_col) && cost_col %in% names(scores) &&
      is.finite(lambda_cost) && lambda_cost != 0) {
    cst <- as.numeric(scores[[cost_col]]); cst[!is.finite(cst)] <- 0
    scores$.linear_gain <- scores$.linear_gain - lambda_cost * cst
  }
  if (!is.null(logistic_col) && logistic_col %in% names(scores) &&
      is.finite(lambda_logistic) && lambda_logistic != 0) {
    lg <- as.numeric(scores[[logistic_col]]); lg[!is.finite(lg)] <- 0
    scores$.linear_gain <- scores$.linear_gain - lambda_logistic * lg
  }
  score_scale <- ng_gain_scale(scores$.linear_gain, n_crosses)
  lambda_parent_use_input <- lambda_parent_use
  if (lambda_parent_use_mode == "adaptive") {
    lambda_parent_use <- lambda_parent_use_input * n_crosses * score_scale
  }
  if (!is.finite(lambda_parent_use) || lambda_parent_use < 0) lambda_parent_use <- 0
  # A finite min_unique_parents is a HARD constraint that the plain linear MIP (ng_mip_linear)
  # cannot express. Route such problems to the exact contribution MIP, which honors it (and with
  # lambda = 0 reduces to pure gain maximization subject to the min-unique / max-use constraints).
  min_unique_is_set <- !is.null(min_unique_parents) && (
    (is.character(min_unique_parents) && tolower(trimws(min_unique_parents[[1L]])) == "auto") ||
    (is.numeric(min_unique_parents) && is.finite(min_unique_parents[[1L]]) && min_unique_parents[[1L]] >= 1))
  if (method == "auto") {
    method <- if ((lambda_group > 0 || lambda_parent_use > 0 || min_unique_is_set) &&
                  requireNamespace("lpSolve", quietly = TRUE)) {
      "mip_contribution"
    } else if (lambda_group == 0 && requireNamespace("lpSolve", quietly = TRUE)) {
      "mip_linear"
    } else {
      "greedy_local"
    }
  }
  # MIP paths are size/time-guarded and NEVER hang the caller: on an oversized
  # problem, an infeasible solve, or a time budget overrun they return NULL and we
  # fall back to greedy_local (recorded in the summary as mip_fallback) instead of
  # erroring. lpSolve's branch-and-bound can blow up on large or degenerate/tied
  # score landscapes, so the size guard (deterministic) is the primary protection.
  mip_fallback <- FALSE
  mip_fallback_reason <- ""
  plan <- NULL
  if (method == "mip_contribution") {
    plan <- ng_mip_contribution(
      scores = scores, n_crosses = n_crosses, parents = parents, parent_kinship = parent_kinship,
      max_crosses_per_parent = max_crosses_per_parent, min_unique_parents = min_unique_parents,
      lambda_parent_use = lambda_parent_use, lambda_group = lambda_group, ocs_iter = ocs_iter,
      mip_time_limit = mip_time_limit, mip_max_binary_vars = mip_max_binary_vars
    )
    if (is.null(plan)) {
      mip_fallback <- TRUE
      mip_fallback_reason <- "mip_contribution unavailable/oversized/timed-out; used greedy_local"
      method <- "greedy_local"
    }
  } else if (method == "mip_linear") {
    plan <- if (min_unique_is_set) {
      # ng_mip_linear cannot express min_unique_parents; use the exact contribution MIP.
      ng_mip_contribution(
        scores = scores, n_crosses = n_crosses, parents = parents, parent_kinship = parent_kinship,
        max_crosses_per_parent = max_crosses_per_parent, min_unique_parents = min_unique_parents,
        lambda_parent_use = lambda_parent_use, lambda_group = lambda_group, ocs_iter = ocs_iter,
        mip_time_limit = mip_time_limit, mip_max_binary_vars = mip_max_binary_vars)
    } else {
      ng_mip_linear(scores, n_crosses, parents, max_crosses_per_parent,
                    mip_max_binary_vars = mip_max_binary_vars)
    }
    if (is.null(plan)) {
      mip_fallback <- TRUE
      mip_fallback_reason <- "mip_linear unavailable/oversized; used greedy_local"
      method <- "greedy_local"
    }
  }
  if (is.null(plan)) {
    if (mip_fallback) warning(mip_fallback_reason, call. = FALSE)
    if (method == "repair_local") {
      plan <- ng_repair_local(scores, n_crosses, parents, parent_kinship, max_crosses_per_parent,
                              min_unique_parents, lambda_group, local_iter,
                              lambda_parent_use = lambda_parent_use)
    } else if (method == "evolution") {
      plan <- ng_evolutionary_mate_allocation(
        scores, n_crosses, parents, parent_kinship, max_crosses_per_parent,
        min_unique_parents = min_unique_parents, lambda_group = lambda_group,
        lambda_parent_use = lambda_parent_use,
        local_iter = min(as.integer(local_iter), 400L),
        evol_solutions = evol_solutions, evol_iterations = evol_iterations,
        evol_stop = evol_stop, seed = evol_seed)
    } else {
      plan <- ng_greedy_local(scores, n_crosses, parents, parent_kinship, max_crosses_per_parent,
                              min_unique_parents, lambda_group, local_iter,
                              lambda_parent_use = lambda_parent_use)
    }
  }
  # Breeder mating constraints, applied repair-stage to the built plan (integer cross
  # indices) so they work with every allocator including MIP. Committed matings are
  # locked and protected from the subsequent min-use / quota drops.
  committed_idx <- ng_resolve_committed_crosses(scores, committed_crosses)
  if (length(committed_idx)) {
    plan <- ng_apply_committed_crosses(scores, plan, committed_idx, parents,
                                       max_crosses_per_parent, n_crosses)
  }
  min_use <- suppressWarnings(as.numeric(min_crosses_per_parent))
  if (is.finite(min_use) && min_use > 1) {
    plan <- ng_enforce_min_use_if_used(scores, plan, parents, max_crosses_per_parent,
                                       min_use, n_crosses, protected = committed_idx)
  }
  if (!is.null(group_quota) && length(group_quota)) {
    plan <- ng_enforce_group_quota(scores, plan, parents, parent_group, group_quota,
                                   max_crosses_per_parent, n_crosses, protected = committed_idx)
  }
  if (!is.null(cost_col) && cost_col %in% names(scores) && is.finite(budget)) {
    plan <- ng_enforce_budget(scores, plan, cost_col, budget, parents,
                              max_crosses_per_parent, n_crosses, protected = committed_idx)
  }
  constraints_reduced_plan <- length(plan) < n_crosses
  if (constraints_reduced_plan) {
    warning(sprintf(
      "Breeder constraints (min-use / group quota / budget) reduced the plan to %d of %d crosses; relax the constraints or add candidates for a full-size plan.",
      length(plan), n_crosses), call. = FALSE)
  }
  out <- ng_plan_summary(plan, scores, gain_col, parent_kinship, lambda_group, lambda_mating,
                         lambda_parent_use = lambda_parent_use,
                         lambda_parent_use_input = lambda_parent_use_input,
                         lambda_parent_use_mode = lambda_parent_use_mode,
                         lambda_progeny_inbreeding = lambda_progeny_inbreeding,
                         score_scale = score_scale)
  s <- attr(out, "summary")
  s$mip_fallback <- mip_fallback
  s$mip_fallback_reason <- mip_fallback_reason
  s$n_committed <- length(committed_idx)
  s$min_use_if_used <- if (is.finite(min_use) && min_use > 1) as.integer(min_use) else NA_integer_
  s$constraints_reduced_plan <- constraints_reduced_plan
  # Flag the (additive, same-axis) overlap so callers/reports can see it was in effect.
  s$relatedness_penalty_overlap <- lambda_mating > 0 && lambda_progeny_inbreeding > 0
  if (!is.null(cost_col) && cost_col %in% names(scores)) {
    plan_cost <- sum(as.numeric(scores[[cost_col]][plan]), na.rm = TRUE)
    s$total_cost <- plan_cost
    s$budget <- budget
    s$over_budget <- is.finite(budget) && plan_cost > budget + 1e-9
  }
  attr(out, "summary") <- s
  out
}

ng_gain_scale <- function(x, n_crosses, pool_multiplier = 10L, min_pool = 100L) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(1)
  pool_n <- min(length(x), max(as.integer(min_pool), as.integer(pool_multiplier * n_crosses)))
  pool <- sort(x, decreasing = TRUE)[seq_len(pool_n)]
  scale <- stats::IQR(pool, na.rm = TRUE) / 1.349
  if (!is.finite(scale) || scale <= 0) scale <- stats::mad(pool, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(pool, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- max(abs(stats::median(pool, na.rm = TRUE)), 1)
  scale
}

ng_rank_normalize <- function(x, bigger_is_better = TRUE) {
  x <- as.numeric(x)
  out <- rep(0, length(x))
  ok <- is.finite(x)
  if (sum(ok) <= 1L) return(out)
  value <- if (isTRUE(bigger_is_better)) x[ok] else -x[ok]
  r <- rank(value, ties.method = "average")
  p <- (r - 0.5) / length(r)
  out[ok] <- stats::qnorm(p)
  finite <- is.finite(out)
  if (any(finite)) out[!finite] <- min(out[finite], na.rm = TRUE) - 1
  out <- ng_standardize(out)
  out[!is.finite(out)] <- 0
  out
}

ng_first_available_col <- function(scores, preferred, fallback) {
  candidates <- unique(c(preferred, fallback))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% names(scores)]
  if (!length(hit)) NULL else hit[[1L]]
}

ng_default_min_unique_parents <- function(scores, n_crosses) {
  parents <- unique(c(as.character(scores$parent1), as.character(scores$parent2)))
  max_unique <- min(length(parents), 2L * as.integer(n_crosses))
  if (!is.finite(max_unique) || max_unique <= 0L) return(NULL)
  target <- max(as.integer(n_crosses) + 2L, ceiling(0.60 * max_unique))
  as.integer(min(max_unique, target))
}

ng_add_balanced_usefulness_score <- function(scores,
                                             gain_col = "usefulness_vpm_gebv",
                                             diversity_col = "parent_distance",
                                             pair_kinship_col = "pair_kinship",
                                             gain_weight = 1.0,
                                             diversity_weight = 0.30,
                                             pair_kinship_weight = 0.10,
                                             strict_gain_col = TRUE,
                                             out_col = ".balanced_usefulness_gain") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  requested_gain_col <- gain_col
  if (isTRUE(strict_gain_col) && length(gain_col) == 1L && !is.na(gain_col) && nzchar(gain_col)) {
    if (!(gain_col %in% names(scores))) ng_stop("scores missing balanced-usefulness gain column: ", gain_col)
  } else {
    gain_col <- ng_first_available_col(
      scores,
      gain_col,
      c("usefulness_vpm_gebv", "etk_vpm_gebv_cal", "usefulness_vpm",
        "etk_vpm_cal", "popvar_uc", "simple_usefa")
    )
  }
  if (is.null(gain_col)) ng_stop("No usable balanced-usefulness gain column found")
  diversity_col <- ng_first_available_col(
    scores,
    diversity_col,
    c("parent_distance_cal", "parent_distance", "vpm_cal", "vpm")
  )
  if (is.null(diversity_col)) diversity_weight <- 0
  if (!(pair_kinship_col %in% names(scores))) pair_kinship_weight <- 0

  gain_z <- ng_rank_normalize(scores[[gain_col]], bigger_is_better = TRUE)
  diversity_z <- if (diversity_weight > 0) {
    ng_rank_normalize(scores[[diversity_col]], bigger_is_better = TRUE)
  } else {
    rep(0, nrow(scores))
  }
  kinship_z <- if (pair_kinship_weight > 0) {
    ng_rank_normalize(scores[[pair_kinship_col]], bigger_is_better = TRUE)
  } else {
    rep(0, nrow(scores))
  }

  scores[[out_col]] <- gain_weight * gain_z +
    diversity_weight * diversity_z -
    pair_kinship_weight * kinship_z
  bad_gain <- !is.finite(as.numeric(scores[[gain_col]]))
  scores[[out_col]][bad_gain] <- -Inf
  attr(scores, "balanced_usefulness") <- list(
    requested_gain_col = requested_gain_col,
    gain_col = gain_col,
    diversity_col = diversity_col,
    pair_kinship_col = if (pair_kinship_weight > 0) pair_kinship_col else NA_character_,
    gain_weight = gain_weight,
    diversity_weight = diversity_weight,
    pair_kinship_weight = pair_kinship_weight,
    out_col = out_col
  )
  scores
}

ng_optimize_balanced_usefulness <- function(scores,
                                            n_crosses,
                                            gain_col = "usefulness_vpm_gebv",
                                            diversity_col = "parent_distance",
                                            parent_kinship = NULL,
                                            max_crosses_per_parent = 10L,
                                            min_unique_parents = NULL,
                                            max_pair_kinship = Inf,
                                            lambda_group = 1.0,
                                            lambda_mating = 0.0,
                                            lambda_parent_use = 2.0,
                                            lambda_parent_use_mode = c("absolute", "adaptive"),
                                            gain_weight = 1.0,
                                            diversity_weight = 0.30,
                                            pair_kinship_weight = 0.10,
                                            strict_gain_col = TRUE,
                                            method = "mip_contribution",
                                            local_iter = 2000,
                                            ocs_iter = 5L,
                                            out_col = ".balanced_usefulness_gain") {
  lambda_parent_use_mode <- match.arg(lambda_parent_use_mode)
  scores <- ng_add_balanced_usefulness_score(
    scores = scores,
    gain_col = gain_col,
    diversity_col = diversity_col,
    gain_weight = gain_weight,
    diversity_weight = diversity_weight,
    pair_kinship_weight = pair_kinship_weight,
    strict_gain_col = strict_gain_col,
    out_col = out_col
  )
  meta <- attr(scores, "balanced_usefulness")
  auto_min_unique <- is.character(min_unique_parents) &&
    length(min_unique_parents) == 1L &&
    tolower(trimws(min_unique_parents)) == "auto"
  if (auto_min_unique) {
    min_unique_parents <- ng_default_min_unique_parents(scores, n_crosses)
  } else if (is.null(min_unique_parents) || !is.finite(min_unique_parents)) {
    min_unique_parents <- NULL
  } else {
    min_unique_parents <- as.integer(min_unique_parents)
  }
  max_unique <- min(length(unique(c(scores$parent1, scores$parent2))), 2L * as.integer(n_crosses))
  if (!is.null(min_unique_parents)) {
    min_unique_parents <- max(1L, min(as.integer(max_unique), as.integer(min_unique_parents)))
  }

  min_unique_options <- list(min_unique_parents)
  if (!is.null(min_unique_parents) && min_unique_parents > n_crosses + 1L) {
    min_unique_options[[length(min_unique_options) + 1L]] <-
      max(as.integer(n_crosses) + 1L, floor(0.85 * min_unique_parents))
  }
  min_unique_options <- c(min_unique_options, list(NULL))

  last_error <- NULL
  relaxation_log <- character(0)
  requested_repr <- if (is.null(min_unique_parents)) "NULL" else as.character(min_unique_parents)
  for (option_idx in seq_along(min_unique_options)) {
    min_unique <- min_unique_options[[option_idx]]
    plan <- tryCatch(
      ng_optimize_mating_plan(
        scores = scores,
        n_crosses = n_crosses,
        gain_col = out_col,
        parent_kinship = parent_kinship,
        max_crosses_per_parent = max_crosses_per_parent,
        min_unique_parents = min_unique,
        max_pair_kinship = max_pair_kinship,
        lambda_group = lambda_group,
        lambda_mating = lambda_mating,
        lambda_parent_use = lambda_parent_use,
        lambda_parent_use_mode = lambda_parent_use_mode,
        method = method,
        local_iter = local_iter,
        ocs_iter = ocs_iter
      ),
      error = function(e) e
    )
    if (!inherits(plan, "error")) {
      if (option_idx > 1L) {
        used_repr <- if (is.null(min_unique)) "NULL (constraint dropped)" else as.character(min_unique)
        warning(sprintf(
          "ng_optimize_balanced_usefulness: min_unique_parents relaxed from %s to %s after %d failed attempt(s) (last error: %s). Reported as balanced_min_unique_used in the plan summary.",
          requested_repr, used_repr, option_idx - 1L,
          if (length(relaxation_log)) relaxation_log[[length(relaxation_log)]] else "n/a"
        ), call. = FALSE)
      }
      summary <- attr(plan, "summary")
      summary$balanced_gain_col <- meta$gain_col
      summary$balanced_diversity_col <- meta$diversity_col
      summary$balanced_gain_weight <- meta$gain_weight
      summary$balanced_diversity_weight <- meta$diversity_weight
      summary$balanced_pair_kinship_weight <- meta$pair_kinship_weight
      summary$balanced_min_unique_requested <- if (is.null(min_unique_parents)) NA_integer_ else min_unique_parents
      summary$balanced_min_unique_used <- if (is.null(min_unique)) NA_integer_ else min_unique
      summary$balanced_min_unique_relaxation_attempts <- option_idx - 1L
      attr(plan, "summary") <- summary
      return(plan)
    }
    relaxation_log <- c(relaxation_log, conditionMessage(plan))
    last_error <- plan
  }
  msg <- if (inherits(last_error, "error")) conditionMessage(last_error) else "unknown optimizer error"
  ng_stop("Balanced usefulness optimizer failed: ", msg)
}

ng_mip_linear <- function(scores, n_crosses, parents, max_crosses_per_parent,
                          mip_max_binary_vars = 2e5) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) return(NULL)
  n <- nrow(scores)
  # Size guard: n binary variables (one per candidate cross). Oversized MIPs route
  # to greedy instead of risking an lpSolve branch-and-bound blow-up.
  if (is.finite(mip_max_binary_vars) && n > mip_max_binary_vars) return(NULL)
  parent_mat <- matrix(0, nrow = length(parents), ncol = n)
  rownames(parent_mat) <- parents
  for (k in seq_len(n)) {
    parent_mat[scores$parent1[k], k] <- parent_mat[scores$parent1[k], k] + 1
    parent_mat[scores$parent2[k], k] <- parent_mat[scores$parent2[k], k] + 1
  }
  const_mat <- rbind(rep(1, n), parent_mat)
  const_dir <- c("=", rep("<=", length(parents)))
  const_rhs <- c(n_crosses, rep(max_crosses_per_parent, length(parents)))
  sol <- lpSolve::lp(
    direction = "max",
    objective.in = scores$.linear_gain,
    const.mat = const_mat,
    const.dir = const_dir,
    const.rhs = const_rhs,
    all.bin = TRUE
  )
  if (sol$status != 0) return(NULL)
  which(sol$solution > 0.5)
}

ng_mip_contribution <- function(scores,
                                n_crosses,
                                parents,
                                parent_kinship,
                                max_crosses_per_parent = n_crosses,
                                min_unique_parents = NULL,
                                lambda_parent_use = 0,
                                lambda_group = 0,
                                ocs_iter = 5L,
                                mip_time_limit = 10,
                                mip_max_binary_vars = 2e5) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) return(NULL)
  n <- nrow(scores)
  p <- length(parents)
  if (!is.finite(max_crosses_per_parent)) max_crosses_per_parent <- n_crosses
  cap <- max(1L, min(as.integer(max_crosses_per_parent), as.integer(n_crosses)))
  # Size guard: n candidate-cross binaries + p*cap increment binaries. Oversized
  # MIPs route to greedy (return NULL) rather than risk an lpSolve blow-up.
  if (is.finite(mip_max_binary_vars) && (n + p * cap) > mip_max_binary_vars) return(NULL)
  mip_start_time <- proc.time()[["elapsed"]]
  parent_mat <- matrix(0, nrow = p, ncol = n, dimnames = list(parents, NULL))
  for (k in seq_len(n)) {
    parent_mat[scores$parent1[k], k] <- parent_mat[scores$parent1[k], k] + 1
    parent_mat[scores$parent2[k], k] <- parent_mat[scores$parent2[k], k] + 1
  }

  inc_parent <- rep(seq_len(p), each = cap)
  inc_level <- rep(seq_len(cap), times = p)
  n_inc <- length(inc_parent)
  inc_mat <- matrix(0, nrow = p, ncol = n_inc)
  inc_mat[cbind(inc_parent, seq_len(n_inc))] <- 1
  total_parent_slots <- 2 * n_crosses
  increment_penalty <- ((2 * inc_level - 1) / (total_parent_slots ^ 2)) * lambda_parent_use

  build_constraints <- function() {
    edge_zero <- matrix(0, nrow = p, ncol = n)
    inc_zero <- matrix(0, nrow = p, ncol = n_inc)
    const_mat <- rbind(
      c(rep(1, n), rep(0, n_inc)),
      cbind(parent_mat, -inc_mat)
    )
    const_dir <- c("=", rep("=", p))
    const_rhs <- c(n_crosses, rep(0, p))

    if (!is.null(min_unique_parents)) {
      first_increment <- integer(n_inc)
      first_increment[inc_level == 1L] <- 1L
      const_mat <- rbind(const_mat, c(rep(0, n), first_increment))
      const_dir <- c(const_dir, ">=")
      const_rhs <- c(const_rhs, min_unique_parents)
    }

    if (cap > 1L) {
      monotone <- matrix(0, nrow = p * (cap - 1L), ncol = n + n_inc)
      row <- 0L
      for (parent_idx in seq_len(p)) {
        parent_cols <- which(inc_parent == parent_idx)
        for (level in 2:cap) {
          row <- row + 1L
          monotone[row, n + parent_cols[level]] <- 1
          monotone[row, n + parent_cols[level - 1L]] <- -1
        }
      }
      const_mat <- rbind(const_mat, monotone)
      const_dir <- c(const_dir, rep("<=", nrow(monotone)))
      const_rhs <- c(const_rhs, rep(0, nrow(monotone)))
    }
    list(mat = const_mat, dir = const_dir, rhs = const_rhs)
  }

  constraints <- build_constraints()
  solve_once <- function(parent_weight) {
    edge_penalty <- parent_weight[scores$parent1] + parent_weight[scores$parent2]
    objective <- c(scores$.linear_gain - edge_penalty, -increment_penalty)
    sol <- lpSolve::lp(
      direction = "max",
      objective.in = objective,
      const.mat = constraints$mat,
      const.dir = constraints$dir,
      const.rhs = constraints$rhs,
      all.bin = TRUE
    )
    if (sol$status != 0) return(NULL)
    which(sol$solution[seq_len(n)] > 0.5)
  }

  parent_weight <- setNames(numeric(p), parents)
  best <- NULL
  best_obj <- -Inf
  last_key <- NULL
  for (iter in seq_len(max(1L, as.integer(ocs_iter)))) {
    # Time budget across the dual-scaling iterations: return best-so-far if the
    # elapsed wall-clock exceeds mip_time_limit (checked between lpSolve calls).
    if (iter > 1L && is.finite(mip_time_limit) &&
        (proc.time()[["elapsed"]] - mip_start_time) > mip_time_limit) {
      return(best)
    }
    selected <- solve_once(parent_weight)
    if (is.null(selected) || length(selected) != n_crosses) return(best)
    obj <- ng_plan_objective_contribution(
      scores = scores,
      selected = selected,
      parent_kinship = parent_kinship,
      lambda_group = lambda_group,
      lambda_parent_use = lambda_parent_use
    )
    if (obj > best_obj + 1e-10) {
      best <- selected
      best_obj <- obj
    }
    key <- paste(sort(selected), collapse = ",")
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    contribution <- counts / sum(counts)
    kc <- as.numeric(parent_kinship[parents, parents, drop = FALSE] %*% contribution)
    parent_weight <- setNames(lambda_group * 2 * kc / sum(counts), parents)
    if (identical(key, last_key) || lambda_group <= 0) break
    last_key <- key
  }
  best
}

# Actively enforce a minimum number of distinct parents in the plan (greedy/repair
# only WARNED before, while the MIP path enforced it). Bounded swap loop: repeatedly
# swap in the highest-gain unselected cross that introduces the most currently-unused
# parents, dropping the lowest-gain selected cross that keeps the add capacity-feasible.
# Stops (best-effort) when the target is unreachable for this candidate pool.
ng_enforce_min_unique_parents <- function(scores, selected, parents,
                                          max_crosses_per_parent, min_unique_parents, n_crosses) {
  if (is.null(min_unique_parents) || !is.finite(min_unique_parents) || min_unique_parents <= 0) return(selected)
  target <- min(as.integer(min_unique_parents), min(2L * n_crosses, length(parents)))
  bump <- function(cnt, pp, delta) { for (q in pp) cnt[q] <- cnt[q] + delta; cnt }
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  cur_unique <- sum(counts > 0)
  guard <- 0L
  while (cur_unique < target && guard < 6L * n_crosses) {
    guard <- guard + 1L
    unused <- names(counts)[counts == 0]
    if (!length(unused)) break
    not_sel <- setdiff(seq_len(nrow(scores)), selected)
    intro <- vapply(not_sel, function(i) sum(c(scores$parent1[i], scores$parent2[i]) %in% unused), integer(1))
    cand <- not_sel[intro > 0]
    if (!length(cand)) break
    cand <- cand[order(-intro[match(cand, not_sel)], -scores$.linear_gain[cand])]  # most unused, then gain
    drop_order <- selected[order(scores$.linear_gain[selected])]                   # drop lowest gain first
    chosen <- NULL
    for (add_idx in cand) {
      add_p <- c(scores$parent1[add_idx], scores$parent2[add_idx])
      for (drop_idx in drop_order) {
        c2 <- bump(counts, c(scores$parent1[drop_idx], scores$parent2[drop_idx]), -1L)
        feasible <- if (add_p[1] == add_p[2]) c2[add_p[1]] + 2L <= max_crosses_per_parent
                    else c2[add_p[1]] + 1L <= max_crosses_per_parent && c2[add_p[2]] + 1L <= max_crosses_per_parent
        if (!isTRUE(feasible)) next
        c3 <- bump(c2, add_p, 1L)
        new_unique <- sum(c3 > 0)
        if (new_unique > cur_unique) {  # only accept strictly-improving swaps -> monotone -> terminates
          chosen <- list(add = add_idx, drop = drop_idx, counts = c3, unique = new_unique)
          break
        }
      }
      if (!is.null(chosen)) break
    }
    if (is.null(chosen)) break
    selected <- c(setdiff(selected, chosen$drop), chosen$add)
    counts <- chosen$counts
    cur_unique <- chosen$unique
  }
  selected
}

ng_greedy_local <- function(scores,
                            n_crosses,
                            parents,
                            parent_kinship,
                            max_crosses_per_parent,
                            min_unique_parents,
                            lambda_group,
                            local_iter = 2000,
                            lambda_parent_use = 0) {
  ord <- order(scores$.linear_gain, decreasing = TRUE)
  selected <- integer(0)
  counts <- setNames(integer(length(parents)), parents)
  for (idx in ord) {
    p <- c(scores$parent1[idx], scores$parent2[idx])
    if (any(counts[p] >= max_crosses_per_parent)) next
    selected <- c(selected, idx)
    counts[p] <- counts[p] + 1L
    if (length(selected) == n_crosses) break
  }
  if (length(selected) < n_crosses) ng_stop("Greedy allocator could not build a feasible plan")
  if (local_iter > 0) {
    selected <- ng_local_swap(scores, selected, parents, parent_kinship, max_crosses_per_parent, lambda_group, local_iter, lambda_parent_use = lambda_parent_use)
  }
  selected <- ng_enforce_min_unique_parents(scores, selected, parents, max_crosses_per_parent, min_unique_parents, n_crosses)
  if (!is.null(min_unique_parents)) {
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    if (sum(counts > 0) < min(min_unique_parents, min(2L * n_crosses, length(parents)))) {
      warning("min_unique_parents is infeasible for this candidate pool; returning best-effort plan", call. = FALSE)
    }
  }
  selected
}

ng_repair_local <- function(scores,
                            n_crosses,
                            parents,
                            parent_kinship,
                            max_crosses_per_parent,
                            min_unique_parents,
                            lambda_group,
                            local_iter = 2000,
                            lambda_parent_use = 0) {
  ord <- order(scores$.linear_gain, decreasing = TRUE)
  selected <- head(ord, n_crosses)
  selected <- ng_repair_capacity(scores, selected, parents, max_crosses_per_parent, n_crosses)
  if (local_iter > 0) {
    selected <- ng_local_swap(scores, selected, parents, parent_kinship, max_crosses_per_parent, lambda_group, local_iter, lambda_parent_use = lambda_parent_use)
    selected <- ng_upgrade_repair(scores, selected, parents, parent_kinship, max_crosses_per_parent, lambda_group, local_iter, lambda_parent_use = lambda_parent_use)
  }
  selected <- ng_enforce_min_unique_parents(scores, selected, parents, max_crosses_per_parent, min_unique_parents, n_crosses)
  if (!is.null(min_unique_parents)) {
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    if (sum(counts > 0) < min(min_unique_parents, min(2L * n_crosses, length(parents)))) {
      warning("min_unique_parents is infeasible for this candidate pool; returning best-effort plan", call. = FALSE)
    }
  }
  selected
}

ng_repair_capacity <- function(scores, selected, parents, max_crosses_per_parent, n_crosses) {
  selected <- unique(selected)
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  while (length(selected) && any(counts > max_crosses_per_parent)) {
    over <- names(counts)[counts > max_crosses_per_parent]
    touches_over <- vapply(selected, function(idx) {
      any(c(scores$parent1[idx], scores$parent2[idx]) %in% over)
    }, logical(1))
    candidates <- selected[touches_over]
    if (!length(candidates)) break
    relief <- vapply(candidates, function(idx) {
      p <- c(scores$parent1[idx], scores$parent2[idx])
      sum(pmax(counts[p] - max_crosses_per_parent, 0))
    }, numeric(1))
    drop_score <- scores$.linear_gain[candidates] / pmax(relief, 1)
    drop_idx <- candidates[which.min(drop_score)]
    selected <- setdiff(selected, drop_idx)
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  }
  ng_fill_best(scores, selected, parents, max_crosses_per_parent, n_crosses)
}

ng_fill_best <- function(scores, selected, parents, max_crosses_per_parent, n_crosses,
                         blocked_parents = integer(0), blocked_pair_keys = character(0),
                         pair_keys = NULL) {
  selected <- unique(selected)
  if (length(selected) >= n_crosses) return(selected[seq_len(n_crosses)])
  # Integer-indexed fill (no per-element character/name hash lookups): map parents to
  # integer indices once and keep a plain integer counts vector. This is on the repair
  # hot path (every greedy/repair/evolution repair) so it materially speeds all paths.
  selected_flag <- logical(nrow(scores)); selected_flag[selected] <- TRUE
  p1i <- match(as.character(scores$parent1), parents)
  p2i <- match(as.character(scores$parent2), parents)
  counts <- tabulate(c(p1i[selected], p2i[selected]), nbins = length(parents))
  ord <- order(scores$.linear_gain, decreasing = TRUE)
  ord <- ord[!selected_flag[ord]]
  need <- n_crosses - length(selected)
  has_block_parent <- length(blocked_parents) > 0L
  has_block_pair <- length(blocked_pair_keys) > 0L && !is.null(pair_keys)
  for (idx in ord) {
    if (need <= 0L) break
    a <- p1i[idx]; b <- p2i[idx]
    if (counts[a] >= max_crosses_per_parent || counts[b] >= max_crosses_per_parent) next
    if (has_block_parent && (a %in% blocked_parents || b %in% blocked_parents)) next
    if (has_block_pair && (pair_keys[[idx]] %in% blocked_pair_keys)) next
    selected <- c(selected, idx)
    counts[a] <- counts[a] + 1L; counts[b] <- counts[b] + 1L
    need <- need - 1L
  }
  if (length(selected) < n_crosses) {
    # When blocking is active (min-use bans a parent, or a group quota blocks a full
    # pair), a full-size plan may be genuinely unreachable. We honour the constraint and
    # return the shorter feasible plan rather than silently re-adding a blocked cross;
    # ng_optimize_mating_plan warns the caller that constraints reduced the plan. Without
    # blocking, an under-filled plan is a real infeasibility and stays an error.
    if (!has_block_parent && !has_block_pair) {
      ng_stop("Allocator could not repair/fill a feasible plan")
    }
  }
  selected
}

ng_local_swap <- function(scores, selected, parents, parent_kinship, max_crosses_per_parent, lambda_group, local_iter, use_cpp = TRUE, lambda_parent_use = 0) {
  # The hot path: thousands of ng_plan_objective() calls per ng_optimize_mating_plan
  # invocation. The C++ kernel fuses the swap-search loop with the K-quadratic
  # form, eliminating the O(p^2) R-level matrix multiply overhead per candidate.
  # Empirically ~30-100x speedup vs the R reference on n=80 fixtures (10-20s ->
  # 100-300ms).
  # The C++ kernel now models the FULL objective (gain - lambda_group * coancestry -
  # lambda_parent_use * sum(contribution^2)), so it is used for both penalties; the
  # pure-R reference below is retained only for when the compiled kernel is
  # unavailable (NGCD_SKIP_CPP / no Rcpp).
  if (isTRUE(use_cpp) &&
      exists("ng_local_swap_cpp", mode = "function", inherits = TRUE)) {
    # Map character parent IDs to 0-based row indices of parent_kinship so the C++
    # kernel can index the kinship matrix directly without character lookups.
    parent_idx <- setNames(seq_along(parents) - 1L, parents)
    p1_zero <- as.integer(parent_idx[scores$parent1])
    p2_zero <- as.integer(parent_idx[scores$parent2])
    if (anyNA(p1_zero) || anyNA(p2_zero)) {
      ng_stop("ng_local_swap: some parent IDs in scores are not in parents[]")
    }
    K_mat <- as.matrix(parent_kinship[parents, parents, drop = FALSE])
    storage.mode(K_mat) <- "double"
    out <- ng_local_swap_cpp(
      linear_gain = as.numeric(scores$.linear_gain),
      pair_p1_zero = p1_zero,
      pair_p2_zero = p2_zero,
      parent_kinship = K_mat,
      selected_zero = as.integer(selected - 1L),
      max_per_parent = as.integer(max_crosses_per_parent),
      lambda_group = as.numeric(lambda_group),
      lambda_parent_use = as.numeric(lambda_parent_use),
      local_iter = as.integer(local_iter),
      off_pool_size = 1000L
    )
    return(as.integer(out) + 1L)
  }
  selected_flag <- rep(FALSE, nrow(scores))
  selected_flag[selected] <- TRUE
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  best_obj <- ng_plan_objective(scores, selected, parent_kinship, lambda_group, lambda_parent_use)
  for (iter in seq_len(local_iter)) {
    improved <- FALSE
    off <- which(!selected_flag)
    on <- which(selected_flag)
    off <- off[order(scores$.linear_gain[off], decreasing = TRUE)]
    on <- on[order(scores$.linear_gain[on], decreasing = FALSE)]
    for (drop_idx in on) {
      drop_p <- c(scores$parent1[drop_idx], scores$parent2[drop_idx])
      counts2 <- counts
      counts2[drop_p] <- counts2[drop_p] - 1L
      for (add_idx in head(off, min(1000L, length(off)))) {
        add_p <- c(scores$parent1[add_idx], scores$parent2[add_idx])
        if (any(counts2[add_p] >= max_crosses_per_parent)) next
        candidate <- c(setdiff(selected, drop_idx), add_idx)
        obj <- ng_plan_objective(scores, candidate, parent_kinship, lambda_group, lambda_parent_use)
        if (obj > best_obj + 1e-10) {
          selected_flag[drop_idx] <- FALSE
          selected_flag[add_idx] <- TRUE
          selected <- candidate
          counts <- counts2
          counts[add_p] <- counts[add_p] + 1L
          best_obj <- obj
          improved <- TRUE
          break
        }
      }
      if (improved) break
    }
    if (!improved) break
  }
  selected
}

ng_upgrade_repair <- function(scores,
                              selected,
                              parents,
                              parent_kinship,
                              max_crosses_per_parent,
                              lambda_group,
                              local_iter = 2000,
                              add_pool = 1000L,
                              lambda_parent_use = 0) {
  selected <- unique(selected)
  selected_flag <- rep(FALSE, nrow(scores))
  selected_flag[selected] <- TRUE
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  best_obj <- ng_plan_objective(scores, selected, parent_kinship, lambda_group, lambda_parent_use)
  n_crosses <- length(selected)
  for (iter in seq_len(local_iter)) {
    improved <- FALSE
    off <- which(!selected_flag)
    off <- head(off[order(scores$.linear_gain[off], decreasing = TRUE)], add_pool)
    for (add_idx in off) {
      add_p <- c(scores$parent1[add_idx], scores$parent2[add_idx])
      counts_after <- counts
      counts_after[add_p] <- counts_after[add_p] + 1L
      over <- names(counts_after)[counts_after > max_crosses_per_parent]
      over <- intersect(over, add_p)
      if (!length(over)) next
      drop_set <- ng_min_loss_drop_set(scores, selected, over)
      if (!length(drop_set)) next
      candidate <- c(setdiff(selected, drop_set), add_idx)
      candidate <- ng_repair_capacity(scores, candidate, parents, max_crosses_per_parent, n_crosses)
      obj <- ng_plan_objective(scores, candidate, parent_kinship, lambda_group, lambda_parent_use)
      if (obj > best_obj + 1e-10) {
        selected <- candidate
        selected_flag <- rep(FALSE, nrow(scores))
        selected_flag[selected] <- TRUE
        counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
        best_obj <- obj
        improved <- TRUE
        break
      }
    }
    if (!improved) break
  }
  selected
}

ng_min_loss_drop_set <- function(scores, selected, over_parents) {
  over_parents <- unique(over_parents)
  touches <- lapply(over_parents, function(p) {
    selected[scores$parent1[selected] == p | scores$parent2[selected] == p]
  })
  if (any(lengths(touches) == 0)) return(integer(0))
  if (length(touches) == 1L) {
    cand <- touches[[1]]
    return(cand[which.min(scores$.linear_gain[cand])])
  }
  best <- integer(0)
  best_loss <- Inf
  for (a in touches[[1]]) {
    for (b in touches[[2]]) {
      drop <- unique(c(a, b))
      loss <- sum(scores$.linear_gain[drop])
      if (loss < best_loss) {
        best_loss <- loss
        best <- drop
      }
    }
  }
  best
}

ng_parent_counts <- function(plan, parents) {
  # Vectorized: match parent IDs to the parents[] index and tabulate in C, instead of
  # an R per-row loop. This is on the optimizer hot path (called by every repair,
  # enforcement, and fitness evaluation) so the speedup benefits all allocators.
  idx <- match(c(as.character(plan$parent1), as.character(plan$parent2)), parents)
  setNames(as.integer(tabulate(idx, nbins = length(parents))), parents)
}

# Group-coancestry term c'Kc of the parent contribution vector c (sums to 1). NOTE
# on units: parent_kinship is typically a VanRaden genomic relationship matrix
# (G ~ numerator relationship A ~ 2 x kinship coefficient), so this returns c'Gc ~
# 2 x Meuwissen group coancestry (i.e. a mean group RELATIONSHIP, not a coancestry
# coefficient in [0,1]). The optimizer ranking is unaffected (the 2x is absorbed into
# lambda_group), but do not interpret the reported value as a coancestry/inbreeding
# rate or compare it to a Delta-F target without dividing by 2.
ng_group_coancestry <- function(counts, parent_kinship) {
  total <- sum(counts)
  if (total <= 0) return(NA_real_)
  cvec <- counts[rownames(parent_kinship)] / total
  as.numeric(crossprod(cvec, parent_kinship %*% cvec))
}

ng_plan_objective <- function(scores, selected, parent_kinship, lambda_group, lambda_parent_use = 0) {
  base <- sum(scores$.linear_gain[selected])
  if (lambda_group <= 0 && lambda_parent_use <= 0) return(base)
  parents <- rownames(parent_kinship)
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  penalty <- 0
  if (lambda_group > 0) penalty <- penalty + lambda_group * ng_group_coancestry(counts, parent_kinship)
  if (lambda_parent_use > 0) {
    contribution <- counts / sum(counts)
    penalty <- penalty + lambda_parent_use * sum(contribution * contribution)
  }
  base - penalty
}

ng_plan_objective_contribution <- function(scores,
                                           selected,
                                           parent_kinship,
                                           lambda_group = 0,
                                           lambda_parent_use = 0) {
  base <- sum(scores$.linear_gain[selected])
  parents <- rownames(parent_kinship)
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  contribution <- counts / sum(counts)
  parent_use_sq <- sum(contribution * contribution)
  group <- ng_group_coancestry(counts, parent_kinship)
  base - lambda_parent_use * parent_use_sq - lambda_group * group
}

ng_plan_summary <- function(selected, scores, gain_col, parent_kinship, lambda_group, lambda_mating,
                            lambda_parent_use = 0,
                            lambda_parent_use_input = lambda_parent_use,
                            lambda_parent_use_mode = "absolute",
                            lambda_progeny_inbreeding = 0,
                            score_scale = NA_real_) {
  plan <- scores[selected, , drop = FALSE]
  parents <- rownames(parent_kinship)
  counts <- ng_parent_counts(plan, parents)
  contribution <- counts / sum(counts)
  progeny_f <- if ("expected_progeny_inbreeding" %in% names(plan)) {
    as.numeric(plan$expected_progeny_inbreeding)
  } else {
    as.numeric(plan$pair_kinship) / 2
  }
  attr(plan, "summary") <- list(
    n_crosses = nrow(plan),
    gain_col = gain_col,
    total_gain = sum(plan[[gain_col]], na.rm = TRUE),
    mean_gain = mean(plan[[gain_col]], na.rm = TRUE),
    mean_pair_kinship = mean(plan$pair_kinship, na.rm = TRUE),
    mean_progeny_inbreeding = mean(progeny_f, na.rm = TRUE),
    max_progeny_inbreeding = suppressWarnings(max(progeny_f, na.rm = TRUE)),
    lambda_progeny_inbreeding = lambda_progeny_inbreeding,
    group_coancestry = ng_group_coancestry(counts, parent_kinship),
    parent_use_sq = sum(contribution * contribution),
    unique_parents = sum(counts > 0),
    max_parent_use = max(counts),
    lambda_group = lambda_group,
    lambda_mating = lambda_mating,
    lambda_parent_use = lambda_parent_use,
    lambda_parent_use_input = lambda_parent_use_input,
    lambda_parent_use_mode = lambda_parent_use_mode,
    score_scale = score_scale,
    parent_counts = counts[counts > 0]
  )
  plan
}

ng_allocate_family_sizes <- function(plan,
                                     total_progeny,
                                     min_progeny = 1L,
                                     max_progeny = Inf,
                                     value_col = NULL,
                                     power = 1,
                                     method = c("score_weighted", "marginal_topk"),
                                     mean_col = "cross_mean",
                                     var_col = NULL,
                                     selected_top_n = NULL,
                                     selected_top_prop = 0.10) {
  method <- match.arg(method)
  if (method == "marginal_topk") {
    return(ng_allocate_family_sizes_marginal(
      plan = plan,
      total_progeny = total_progeny,
      min_progeny = min_progeny,
      max_progeny = max_progeny,
      mean_col = mean_col,
      var_col = var_col,
      selected_top_n = selected_top_n,
      selected_top_prop = selected_top_prop
    ))
  }
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  n <- nrow(plan)
  if (!n) return(plan)
  total_progeny <- as.integer(total_progeny)
  min_progeny <- as.integer(min_progeny)
  if (!is.finite(total_progeny) || total_progeny < n) ng_stop("total_progeny must be at least the number of crosses")
  if (!is.finite(min_progeny) || min_progeny < 1L) min_progeny <- 1L
  if (!is.finite(max_progeny)) max_progeny <- total_progeny
  max_progeny <- max(min_progeny, as.integer(max_progeny))
  if (total_progeny < n * min_progeny) ng_stop("total_progeny is too small for min_progeny")
  if (total_progeny > n * max_progeny) ng_stop("total_progeny is too large for max_progeny")

  size <- rep(min_progeny, n)
  extra <- total_progeny - sum(size)
  if (extra > 0) {
    value <- if (!is.null(value_col) && value_col %in% names(plan)) plan[[value_col]] else seq_len(n)
    value <- as.numeric(value)
    finite_value <- is.finite(value)
    if (!any(finite_value)) {
      value <- rep(1, n)
    } else {
      value[!finite_value] <- min(value[finite_value], na.rm = TRUE)
    }
    priority <- value - min(value, na.rm = TRUE)
    if (!any(priority > 0, na.rm = TRUE)) priority <- rep(1, n)
    priority <- pmax(priority, 0) ^ max(0, power)
    if (!any(priority > 0)) priority <- rep(1, n)
    room <- max_progeny - size
    while (extra > 0 && any(room > 0)) {
      p <- priority
      p[room <= 0] <- 0
      if (!any(p > 0)) p <- as.numeric(room > 0)
      quota <- extra * p / sum(p)
      add <- pmin(room, floor(quota))
      if (sum(add) == 0L) {
        frac <- quota - floor(quota)
        eligible <- which(room > 0)
        pick <- eligible[which.max(frac[eligible])]
        add[pick] <- 1L
      }
      size <- size + add
      room <- max_progeny - size
      extra <- total_progeny - sum(size)
    }
  }
  plan$n_progeny <- as.integer(size)
  attr(plan, "family_size_summary") <- list(
    total_progeny = sum(size),
    min_progeny = min(size),
    max_progeny = max(size),
    mean_progeny = mean(size),
    value_col = value_col,
    power = power
  )
  plan
}

ng_allocate_family_sizes_marginal <- function(plan,
                                              total_progeny,
                                              min_progeny = 1L,
                                              max_progeny = Inf,
                                              mean_col = "cross_mean",
                                              var_col,
                                              selected_top_n = NULL,
                                              selected_top_prop = 0.10) {
  plan <- as.data.frame(plan, stringsAsFactors = FALSE)
  n <- nrow(plan)
  if (!n) return(plan)
  if (missing(var_col) || is.null(var_col) || !(var_col %in% names(plan))) {
    ng_stop("marginal family allocation requires a valid variance column")
  }
  if (!(mean_col %in% names(plan))) ng_stop("plan missing mean_col: ", mean_col)
  total_progeny <- as.integer(total_progeny)
  min_progeny <- as.integer(min_progeny)
  if (!is.finite(total_progeny) || total_progeny < n) ng_stop("total_progeny must be at least the number of crosses")
  if (!is.finite(min_progeny) || min_progeny < 1L) min_progeny <- 1L
  if (!is.finite(max_progeny)) max_progeny <- total_progeny
  max_progeny <- max(min_progeny, as.integer(max_progeny))
  if (total_progeny < n * min_progeny) ng_stop("total_progeny is too small for min_progeny")
  if (total_progeny > n * max_progeny) ng_stop("total_progeny is too large for max_progeny")

  mu <- as.numeric(plan[[mean_col]])
  v <- pmax(as.numeric(plan[[var_col]]), 1e-8)
  ok_mu <- is.finite(mu)
  if (!any(ok_mu)) mu <- rep(0, n) else mu[!ok_mu] <- min(mu[ok_mu], na.rm = TRUE)
  v[!is.finite(v) | v <= 0] <- 1e-8
  sd <- sqrt(v)
  size <- rep(min_progeny, n)
  if (is.null(selected_top_n) || !is.finite(selected_top_n)) {
    selected_top_n <- max(1L, ceiling(total_progeny * selected_top_prop))
  }
  selected_top_n <- max(1L, min(total_progeny - 1L, as.integer(selected_top_n)))

  extra <- total_progeny - sum(size)
  last_threshold <- NA_real_
  if (extra > 0) {
    room <- max_progeny - size
    tie_break <- rank(-mu, ties.method = "first") * .Machine$double.eps
    for (step in seq_len(extra)) {
      threshold <- ng_family_selection_threshold(mu, sd, size, selected_top_n)
      last_threshold <- threshold
      marginal <- ng_expected_excess_above_threshold(mu, sd, threshold)
      marginal[room <= 0] <- -Inf
      if (!any(is.finite(marginal))) break
      pick <- which.max(marginal - tie_break)
      size[pick] <- size[pick] + 1L
      room[pick] <- room[pick] - 1L
    }
  }
  if (!is.finite(last_threshold)) {
    last_threshold <- ng_family_selection_threshold(mu, sd, size, selected_top_n)
  }
  plan$n_progeny <- as.integer(size)
  attr(plan, "family_size_summary") <- list(
    total_progeny = sum(size),
    min_progeny = min(size),
    max_progeny = max(size),
    mean_progeny = mean(size),
    method = "marginal_topk",
    mean_col = mean_col,
    var_col = var_col,
    selected_top_n = selected_top_n,
    threshold = last_threshold
  )
  plan
}

ng_family_selection_threshold <- function(mu, sd, n_progeny, selected_top_n) {
  mu <- as.numeric(mu)
  sd <- pmax(as.numeric(sd), 1e-6)
  n_progeny <- pmax(as.numeric(n_progeny), 0)
  target <- max(1e-6, min(sum(n_progeny) - 1e-6, as.numeric(selected_top_n)))
  f <- function(t) {
    sum(n_progeny * stats::pnorm((t - mu) / sd, lower.tail = FALSE)) - target
  }
  lo <- min(mu - 10 * sd, na.rm = TRUE)
  hi <- max(mu + 10 * sd, na.rm = TRUE)
  for (i in seq_len(20L)) {
    if (f(lo) >= 0) break
    lo <- lo - 2 * max(sd, na.rm = TRUE)
  }
  for (i in seq_len(20L)) {
    if (f(hi) <= 0) break
    hi <- hi + 2 * max(sd, na.rm = TRUE)
  }
  if (f(lo) < 0) return(lo)
  if (f(hi) > 0) return(hi)
  stats::uniroot(f, lower = lo, upper = hi, tol = 1e-6)$root
}

ng_expected_excess_above_threshold <- function(mu, sd, threshold) {
  sd <- pmax(as.numeric(sd), 1e-6)
  z <- (threshold - mu) / sd
  pmax(0, sd * stats::dnorm(z) + (mu - threshold) * stats::pnorm(z, lower.tail = FALSE))
}

ng_pareto_mate_allocation <- function(scores,
                                      n_crosses,
                                      gain_col = "usefulness_pmv_gebv",
                                      parent_kinship,
                                      lambdas = 10 ^ seq(-2, 3, length.out = 10),
                                      method = "greedy_local",
                                      ...) {
  plans <- vector("list", length(lambdas))
  for (i in seq_along(lambdas)) {
    plans[[i]] <- ng_optimize_mating_plan(
      scores = scores,
      n_crosses = n_crosses,
      gain_col = gain_col,
      parent_kinship = parent_kinship,
      lambda_group = lambdas[i],
      method = method,
      ...
    )
  }
  names(plans) <- paste0("lambda_", format(lambdas, scientific = TRUE))
  frontier <- do.call(rbind, lapply(seq_along(plans), function(i) {
    s <- attr(plans[[i]], "summary")
    data.frame(
      lambda_group = lambdas[i],
      total_gain = s$total_gain,
      mean_gain = s$mean_gain,
      group_coancestry = s$group_coancestry,
      mean_pair_kinship = s$mean_pair_kinship,
      unique_parents = s$unique_parents,
      max_parent_use = s$max_parent_use
    )
  }))
  list(frontier = frontier, plans = plans)
}
