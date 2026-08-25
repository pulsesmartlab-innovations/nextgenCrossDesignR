# Gain-diversity balancing: pick an operating point on the mate-allocation frontier.
#
# Instead of asking the breeder for a raw coancestry penalty weight (lambda_group), this
# layer exposes a single intuitive control -- a strategy label ("high_gain" | "balanced"
# | "diversity") or a numeric `diversity_emphasis` on a 0 (all genetic gain) .. 100 (all
# diversity) scale -- and translates it into the lambda_group that lands the plan at (or
# near) the requested point on the native gain-vs-coancestry frontier. The single-dial
# navigation of the gain/diversity trade-off follows the tactical-mate-selection ideas in
# Kinghorn (2011), Genet. Sel. Evol. 43:4; the implementation here is our own.
#
# We reuse the native frontier sweep (ng_pareto_mate_allocation) and the shared
# emphasis->frontier-point mapper (ng_choose_frontier_point, R/27), then bisection-refine
# lambda_group between the two bracketing grid points so the achieved emphasis lands
# within `tol` of the target when the grid is too coarse. This is an ergonomics layer
# over the same objective, not a new allocator; claims stay evidence-scoped.

# Strategy label -> diversity emphasis on the 0 (all gain) .. 100 (all diversity) scale.
ng_strategy_emphasis <- function(strategy = c("balanced", "high_gain", "diversity")) {
  strategy <- match.arg(strategy)
  switch(strategy,
    high_gain = 15,
    balanced = 45,
    diversity = 75)
}

# Achieved diversity emphasis (0..100) implied by a plan's realized mean gain, given the
# gain envelope [min_gain, max_gain] observed across the frontier sweep. Emphasis grows as
# realized gain falls from the max-gain end toward the min-gain (max-diversity) end.
ng_frontier_achieved_emphasis <- function(mean_gain, gain_range) {
  gain_range <- suppressWarnings(as.numeric(gain_range))
  if (length(gain_range) < 2L || !all(is.finite(gain_range)) || diff(range(gain_range)) <= 0) {
    return(NA_real_)
  }
  lo <- min(gain_range)
  span <- diff(range(gain_range))
  100 * (1 - (as.numeric(mean_gain) - lo) / span)
}

# Find the lambda_group bracket [lo, hi] on the (geometric) sweep grid whose achieved
# emphases straddle the target, so bisection has something to refine within.
ng_strategy_lambda_bracket <- function(frontier, lambdas, emphasis, gain_range) {
  eff <- vapply(seq_len(nrow(frontier)),
                function(i) ng_frontier_achieved_emphasis(frontier$mean_gain[[i]], gain_range),
                numeric(1))
  ord <- order(lambdas)
  lam <- lambdas[ord]
  eff <- eff[ord]
  ok <- is.finite(eff)
  lam <- lam[ok]
  eff <- eff[ok]
  if (!length(lam)) return(list(lo = NA_real_, hi = NA_real_))
  below <- which(eff <= emphasis)
  above <- which(eff >= emphasis)
  lo <- if (length(below)) lam[max(below)] else min(lam)
  hi <- if (length(above)) lam[min(above)] else max(lam)
  if (hi < lo) {
    tmp <- lo; lo <- hi; hi <- tmp
  }
  list(lo = lo, hi = hi)
}

ng_select_by_strategy <- function(scores,
                                  n_crosses,
                                  parent_kinship,
                                  gain_col = "usefulness_pmv_gebv",
                                  strategy = NULL,
                                  diversity_emphasis = NULL,
                                  lambdas = 10 ^ seq(-2, 3, length.out = 10),
                                  refine_iter = 12L,
                                  tol = 2,
                                  method = "greedy_local",
                                  ...) {
  if (is.null(diversity_emphasis)) {
    diversity_emphasis <- ng_strategy_emphasis(if (is.null(strategy)) "balanced" else strategy)
  }
  diversity_emphasis <- suppressWarnings(as.numeric(diversity_emphasis))
  if (!is.finite(diversity_emphasis)) diversity_emphasis <- 45
  diversity_emphasis <- max(0, min(100, diversity_emphasis))

  sweep <- ng_pareto_mate_allocation(
    scores = scores,
    n_crosses = n_crosses,
    gain_col = gain_col,
    parent_kinship = parent_kinship,
    lambdas = lambdas,
    method = method,
    ...
  )
  frontier <- sweep$frontier
  gain_range <- range(suppressWarnings(as.numeric(frontier$mean_gain)),
                      na.rm = TRUE, finite = TRUE)

  idx <- ng_choose_frontier_point(frontier, emphasis = diversity_emphasis, mode = "target")
  best_plan <- sweep$plans[[idx]]
  best_lambda <- frontier$lambda_group[[idx]]
  best_emphasis <- ng_frontier_achieved_emphasis(frontier$mean_gain[[idx]], gain_range)
  refine_evals <- 0L

  # Bisection refinement: higher lambda_group -> more diversity emphasis -> higher
  # achieved emphasis, so achieved_emphasis(lambda) is monotone non-decreasing and we can
  # bisect the bracketing grid interval on the geometric lambda scale.
  bracket <- ng_strategy_lambda_bracket(frontier, lambdas, diversity_emphasis, gain_range)
  lo <- bracket$lo
  hi <- bracket$hi
  if (is.finite(lo) && is.finite(hi) && hi > lo && is.finite(best_emphasis)) {
    for (k in seq_len(as.integer(refine_iter))) {
      if (abs(best_emphasis - diversity_emphasis) <= tol) break
      mid <- sqrt(lo * hi)
      cand <- ng_optimize_mating_plan(
        scores = scores, n_crosses = n_crosses, gain_col = gain_col,
        parent_kinship = parent_kinship, lambda_group = mid, method = method, ...)
      refine_evals <- refine_evals + 1L
      cand_gain <- attr(cand, "summary")$mean_gain
      cand_emphasis <- ng_frontier_achieved_emphasis(cand_gain, gain_range)
      if (!is.finite(cand_emphasis)) break
      if (abs(cand_emphasis - diversity_emphasis) < abs(best_emphasis - diversity_emphasis)) {
        best_plan <- cand
        best_lambda <- mid
        best_emphasis <- cand_emphasis
      }
      if (cand_emphasis < diversity_emphasis) lo <- mid else hi <- mid
    }
  }

  s <- attr(best_plan, "summary")
  s$strategy <- if (is.null(strategy)) NA_character_ else strategy
  s$diversity_emphasis <- diversity_emphasis
  s$achieved_emphasis <- best_emphasis
  s$emphasis_gap <- diversity_emphasis - best_emphasis
  s$selected_lambda_group <- best_lambda
  s$strategy_refine_evals <- refine_evals
  s$frontier <- frontier
  attr(best_plan, "summary") <- s
  best_plan
}

# Constrained-OCS operating point: the classic Meuwissen formulation -- maximize genetic
# gain subject to a cap on parental (group) coancestry (equivalently a target rate of
# inbreeding). This is a SECOND way to specify the SAME frontier point that
# ng_select_by_strategy navigates by emphasis; it is NOT a new penalty axis. Instead of a
# relative emphasis, the breeder gives the absolute constraint their program actually
# reports -- a target group coancestry -- and we auto-solve the lambda_group that maximizes
# gain while staying at/under it.
#
# `target_coancestry` is on the SAME scale that ng_plan_summary reports
# `group_coancestry`: c'Gc/2 for a VanRaden genomic-relationship matrix G.
# `group_relationship = c'Gc` is reported separately.
#
# group_coancestry decreases monotonically as lambda_group increases, so the max-gain plan
# meeting the cap is the SMALLEST lambda whose coancestry <= target. We locate it on the
# native frontier sweep, then bisection-refine between the bracketing grid lambdas, always
# keeping the plan feasible.
ng_select_by_target_coancestry <- function(scores,
                                           n_crosses,
                                           parent_kinship,
                                           gain_col = "usefulness_pmv_gebv",
                                           target_coancestry,
                                           lambdas = 10 ^ seq(-2, 3, length.out = 10),
                                           refine_iter = 12L,
                                           method = "greedy_local",
                                           ...) {
  target_coancestry <- suppressWarnings(as.numeric(target_coancestry)[[1L]])
  if (!is.finite(target_coancestry)) ng_stop("target_coancestry must be a finite number")

  sweep <- ng_pareto_mate_allocation(
    scores = scores, n_crosses = n_crosses, gain_col = gain_col,
    parent_kinship = parent_kinship, lambdas = lambdas, method = method, ...)
  fr <- sweep$frontier
  ord <- order(fr$lambda_group)
  lam <- fr$lambda_group[ord]
  coanc <- suppressWarnings(as.numeric(fr$group_coancestry))[ord]
  plans <- sweep$plans[ord]
  eps <- 1e-12

  feasible <- is.finite(coanc) & (coanc <= target_coancestry + eps)
  if (!any(feasible)) {
    # Even the strongest diversity penalty on the sweep cannot reach the target: return
    # the minimum-coancestry plan (closest feasible) and warn.
    j <- which.min(coanc)
    plan <- plans[[j]]
    best_lambda <- lam[[j]]
    status <- "infeasible"
  } else {
    j <- min(which(feasible))            # smallest lambda meeting the cap == most gain
    plan <- plans[[j]]
    best_lambda <- lam[[j]]
    status <- if (j == 1L) "met_slack" else "met"
    if (identical(status, "met")) {
      # Refine between the last infeasible grid lambda (lam[j-1], coanc > target) and this
      # feasible one (lam[j]); push lambda down for more gain while staying feasible.
      lo <- lam[[j - 1L]]
      hi <- lam[[j]]
      for (k in seq_len(as.integer(refine_iter))) {
        mid <- sqrt(lo * hi)
        cand <- ng_optimize_mating_plan(
          scores = scores, n_crosses = n_crosses, gain_col = gain_col,
          parent_kinship = parent_kinship, lambda_group = mid, method = method, ...)
        cc <- attr(cand, "summary")$group_coancestry
        if (is.finite(cc) && cc <= target_coancestry + eps) {
          plan <- cand
          best_lambda <- mid
          hi <- mid
        } else {
          lo <- mid
        }
      }
    }
  }

  s <- attr(plan, "summary")
  s$target_coancestry <- target_coancestry
  s$achieved_coancestry <- s$group_coancestry
  s$coancestry_gap <- target_coancestry - s$group_coancestry
  s$target_coancestry_status <- status
  s$selected_lambda_group <- best_lambda
  s$frontier <- fr
  attr(plan, "summary") <- s
  if (identical(status, "infeasible")) {
    warning(sprintf(
      "target_coancestry %.4g is below the minimum achievable on the frontier (%.4g); returned the minimum-coancestry plan.",
      target_coancestry, s$group_coancestry), call. = FALSE)
  }
  plan
}
