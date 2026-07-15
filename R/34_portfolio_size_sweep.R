# v0.3.2 — Portfolio-size sweep + K-selection criteria
#
# Wraps ng_optimize_mating_plan() in a loop over K (the number of crosses)
# and returns the marginal-utility curve plus an estimated K under any of
# several criteria. Helps breeders pick K honestly instead of guessing.
#
# All real allocation logic lives in ng_optimize_mating_plan(); this file
# only iterates and post-processes. Four criteria are offered:
#   - "elbow_relative" (default): smallest K where the marginal gain
#       drops below `relative_threshold` * (first marginal in the sweep).
#       Simple and breeder-friendly. Scale-stable.
#   - "elbow_kneedle": Satopaa et al. (2011) normalized-curvature
#       algorithm. Robust on smooth concave curves; returns NA when the
#       curve is too noisy, too linear, or has the maximum at a boundary.
#   - "ne_target": smallest K with Ne_estimate >= ne_min (default 30).
#       Ne_estimate = 1 / (2 * group_coancestry), Falconer-Mackay
#       approximation. Sustainability-aware: keeps the program above an
#       effective-population-size floor.
#   - "coancestry_budget": largest K with group_coancestry <= coancestry_max
#       (default 0.05). Hard ceiling on relatedness pressure.
#
# attr(., "elbow_K") is the recommended K under whatever criterion was
# used; attr(., "criterion") records which criterion produced it.

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("K", "total_gain"))
}

ng_optimize_mating_plan_curve <- function(scores,
                                          K_range = NULL,
                                          ...,
                                          criterion = c("elbow_relative",
                                                         "elbow_kneedle",
                                                         "ne_target",
                                                         "coancestry_budget"),
                                          relative_threshold = 0.05,
                                          ne_min = 30,
                                          coancestry_max = 0.05) {
  criterion <- match.arg(criterion)
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (is.null(K_range)) {
    K_max <- min(20L, nrow(scores))
    if (K_max < 3L) ng_stop("scores has fewer than 3 candidate crosses; cannot sweep")
    K_range <- seq.int(3L, K_max)
  }
  K_range <- as.integer(K_range)
  if (any(K_range < 1L)) ng_stop("K_range entries must be >= 1")
  if (anyDuplicated(K_range)) ng_stop("K_range must not contain duplicates")
  K_range <- sort(K_range)

  rows <- vector("list", length(K_range))
  for (i in seq_along(K_range)) {
    K <- K_range[[i]]
    plan <- ng_optimize_mating_plan(scores = scores, n_crosses = K, ...)
    plan_summary <- attr(plan, "summary")
    rows[[i]] <- data.frame(
      K = K,
      total_gain        = plan_summary$total_gain,
      mean_gain         = plan_summary$mean_gain,
      group_coancestry  = plan_summary$group_coancestry,
      unique_parents    = plan_summary$unique_parents,
      stringsAsFactors  = FALSE
    )
  }
  curve <- do.call(rbind, rows)
  rownames(curve) <- NULL

  curve$marginal_gain <- c(NA_real_, diff(curve$total_gain))
  base_marg <- curve$marginal_gain[2L]
  curve$relative_marginal <- if (is.finite(base_marg) && base_marg > 0) {
    curve$marginal_gain / base_marg
  } else {
    rep(NA_real_, nrow(curve))
  }
  # Ne_estimate: Falconer-Mackay approximation. Inf when coancestry is 0.
  curve$Ne_estimate <- ifelse(
    is.finite(curve$group_coancestry) & curve$group_coancestry > 0,
    1 / (2 * curve$group_coancestry),
    Inf
  )

  elbow_K <- switch(criterion,
    elbow_relative    = .ng_elbow_relative_marginal(curve, relative_threshold),
    elbow_kneedle     = .ng_elbow_kneedle(curve),
    ne_target         = .ng_pick_ne_target(curve, ne_min),
    coancestry_budget = .ng_pick_coancestry_budget(curve, coancestry_max)
  )
  attr(curve, "elbow_K")  <- elbow_K
  attr(curve, "criterion") <- criterion
  if (criterion == "elbow_relative")
    attr(curve, "relative_threshold") <- relative_threshold
  if (criterion == "ne_target")
    attr(curve, "ne_min") <- ne_min
  if (criterion == "coancestry_budget")
    attr(curve, "coancestry_max") <- coancestry_max
  curve
}

.ng_elbow_relative_marginal <- function(curve, threshold) {
  rel <- curve$relative_marginal
  candidates <- which(is.finite(rel) & rel < threshold)
  if (!length(candidates)) return(NA_integer_)
  as.integer(curve$K[candidates[[1L]]])
}

# Kneedle elbow detection (Satopaa et al. 2011). For a concave-down
# diminishing-returns curve, normalize x and y to [0, 1], compute the
# difference d = y_norm - x_norm, and return the K where d is maximized.
# Returns NA_integer_ when:
#   - fewer than 4 points (cannot define a knee reliably)
#   - the maximum is at the curve boundary (no interior knee)
#   - the curve is nearly linear (max(d) below a small tolerance)
.ng_elbow_kneedle <- function(curve) {
  if (nrow(curve) < 4L) return(NA_integer_)
  x <- as.numeric(curve$K)
  y <- as.numeric(curve$total_gain)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NA_integer_)
  x <- x[ok]; y <- y[ok]; K <- curve$K[ok]
  x_norm <- (x - min(x)) / max(diff(range(x)), .Machine$double.eps)
  y_norm <- (y - min(y)) / max(diff(range(y)), .Machine$double.eps)
  d <- y_norm - x_norm
  if (max(d, na.rm = TRUE) < 1e-3) return(NA_integer_)
  idx <- which.max(d)
  if (idx == 1L || idx == length(d)) return(NA_integer_)
  as.integer(K[[idx]])
}

# Smallest K (in the sweep) at which Ne_estimate >= ne_min.
# Returns NA_integer_ if no K reaches the floor.
.ng_pick_ne_target <- function(curve, ne_min) {
  ok <- (is.finite(curve$Ne_estimate) | is.infinite(curve$Ne_estimate)) &
        curve$Ne_estimate >= ne_min
  candidates <- which(ok)
  if (!length(candidates)) return(NA_integer_)
  as.integer(curve$K[candidates[[1L]]])
}

# Largest K (in the sweep) with group_coancestry <= coancestry_max.
.ng_pick_coancestry_budget <- function(curve, coancestry_max) {
  ok <- is.finite(curve$group_coancestry) &
        curve$group_coancestry <= coancestry_max
  candidates <- which(ok)
  if (!length(candidates)) return(NA_integer_)
  as.integer(curve$K[candidates[[length(candidates)]]])
}

# Render the diminishing-returns curve with the detected K marked.
# Pure ggplot2 helper; lazy-loaded so the package doesn't hard-depend on
# ggplot2 at install time.
ng_plot_diminishing_returns <- function(curve, show_elbow = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    ng_stop("ggplot2 is required for ng_plot_diminishing_returns; ",
            "install.packages(\"ggplot2\") or build the plot yourself ",
            "from the curve data frame.")
  }
  elbow_K <- attr(curve, "elbow_K")
  criterion <- if (is.null(attr(curve, "criterion"))) "elbow_relative" else attr(curve, "criterion")
  p <- ggplot2::ggplot(curve, ggplot2::aes(x = K, y = total_gain)) +
    ggplot2::geom_line(color = "#2b7d6b", linewidth = 0.6) +
    ggplot2::geom_point(color = "#2b7d6b", size = 2.2) +
    ggplot2::labs(
      title = "Diminishing returns on portfolio size",
      subtitle = sprintf("Total selection gain vs number of crosses; criterion = %s",
                         criterion),
      x = "K (number of crosses)",
      y = "Total selection gain"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  if (isTRUE(show_elbow) && length(elbow_K) == 1L && is.finite(elbow_K)) {
    p <- p +
      ggplot2::geom_vline(xintercept = elbow_K, color = "#b85c38",
                           linetype = "dashed", linewidth = 0.5) +
      ggplot2::annotate("text", x = elbow_K, y = max(curve$total_gain, na.rm = TRUE),
                        label = sprintf("K=%d", elbow_K),
                        vjust = -0.4, hjust = -0.1,
                        color = "#b85c38", size = 3.5)
  }
  p
}
