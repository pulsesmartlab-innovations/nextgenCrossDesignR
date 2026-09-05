ng_priority_validate_breaks <- function(breaks) {
  breaks <- suppressWarnings(as.numeric(breaks))
  if (length(breaks) != 4L || any(!is.finite(breaks))) {
    ng_stop("breaks must contain four cumulative proportions, e.g. c(0.10, 0.35, 0.70, 1)")
  }
  if (any(breaks <= 0) || any(breaks > 1) || is.unsorted(breaks, strictly = TRUE)) {
    ng_stop("breaks must be strictly increasing values in (0, 1]")
  }
  if (abs(utils::tail(breaks, 1L) - 1) > 1e-8) ng_stop("last priority break must be 1")
  breaks
}

ng_priority_component <- function(x, bigger_is_better = TRUE) {
  if (length(x) == 0L) return(numeric())
  ng_rank_normalize(as.numeric(x), bigger_is_better = bigger_is_better)
}

ng_rank_cross_priority <- function(crosses,
                                   score_col = "multi_trait_score",
                                   pair_kinship_col = "pair_kinship",
                                   threshold_violation_col = "multi_trait_threshold_violation",
                                   check_violation_col = "check_violation",
                                   score_weight = 1,
                                   kinship_weight = 0.15,
                                   threshold_weight = 1,
                                   check_weight = 0,
                                   breaks = c(0.10, 0.35, 0.70, 1.00),
                                   labels = c("highly_priority", "priority", "medium_priority", "low_priority"),
                                   sort = TRUE) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  if (!nrow(crosses)) return(crosses)
  if (!(score_col %in% names(crosses))) ng_stop("crosses missing score_col: ", score_col)
  breaks <- ng_priority_validate_breaks(breaks)
  labels <- as.character(labels)
  if (length(labels) != 4L || any(!nzchar(labels) | is.na(labels))) {
    ng_stop("labels must contain four non-empty priority labels")
  }
  score_weight <- suppressWarnings(as.numeric(score_weight[[1L]]))
  kinship_weight <- suppressWarnings(as.numeric(kinship_weight[[1L]]))
  threshold_weight <- suppressWarnings(as.numeric(threshold_weight[[1L]]))
  check_weight <- suppressWarnings(as.numeric(check_weight[[1L]]))
  if (!is.finite(score_weight) || score_weight <= 0) ng_stop("score_weight must be positive")
  if (!is.finite(kinship_weight) || kinship_weight < 0) ng_stop("kinship_weight must be >= 0")
  if (!is.finite(threshold_weight) || threshold_weight < 0) ng_stop("threshold_weight must be >= 0")
  if (!is.finite(check_weight) || check_weight < 0) ng_stop("check_weight must be >= 0")

  score <- suppressWarnings(as.numeric(crosses[[score_col]]))
  if (!any(is.finite(score))) ng_stop("score_col contains no finite values")
  score_component <- ng_priority_component(score, bigger_is_better = TRUE)
  kinship_component <- rep(0, nrow(crosses))
  threshold_component <- rep(0, nrow(crosses))
  check_component <- rep(0, nrow(crosses))

  has_kinship <- nzchar(pair_kinship_col) && pair_kinship_col %in% names(crosses)
  if (isTRUE(has_kinship)) {
    kinship_component <- ng_priority_component(crosses[[pair_kinship_col]], bigger_is_better = FALSE)
  }
  has_threshold <- nzchar(threshold_violation_col) && threshold_violation_col %in% names(crosses)
  if (isTRUE(has_threshold)) {
    threshold_component <- ng_priority_component(crosses[[threshold_violation_col]], bigger_is_better = FALSE)
  }
  # The check is a REFERENCE, not a selection criterion (package owner's rule): it is one
  # optional, weighted component of the same shape as kinship/threshold, defaulting to weight 0
  # so an existing run's tiers never move unless the breeder opts in. check_violation itself is
  # the dimensionless, NA-aware wrong-side COUNT ng_attach_check_reference() computes (R/51);
  # fewer violations is better, exactly like threshold_violation_col. It never vetoes a cross --
  # it only shifts the blended index, same as threshold_component.
  #
  # Gated on check_weight > 0, NOT merely on column presence: check_violation exists on every
  # check-enabled run regardless of check_weight (D-series feature), but this component (and
  # priority_rule's naming of it) must report EXACTLY as if the column were absent whenever the
  # caller has not opted in -- otherwise a check-enabled run at the default check_weight = 0 would
  # report a different priority_check_component / priority_rule than a run with no checks at all,
  # even though priority_index/rank/tier are unchanged. That distinction (a check run vs a
  # no-check run differing ONLY in check_violation and the columns ng_attach_check_reference()
  # itself adds) is exactly what check_reference_invariant.R holds the whole feature to.
  has_check <- isTRUE(check_weight > 0) &&
    nzchar(check_violation_col) && check_violation_col %in% names(crosses)
  if (isTRUE(has_check)) {
    check_component <- ng_priority_component(crosses[[check_violation_col]], bigger_is_better = FALSE)
  }

  index <- score_weight * score_component +
    kinship_weight * kinship_component +
    threshold_weight * threshold_component +
    check_weight * check_component
  index[!is.finite(index)] <- min(index[is.finite(index)], na.rm = TRUE) - 1
  tie_kinship <- if (isTRUE(has_kinship)) suppressWarnings(as.numeric(crosses[[pair_kinship_col]])) else rep(0, nrow(crosses))
  tie_threshold <- if (isTRUE(has_threshold)) suppressWarnings(as.numeric(crosses[[threshold_violation_col]])) else rep(0, nrow(crosses))
  ord <- order(-index, -score, tie_threshold, tie_kinship, na.last = TRUE)
  rank <- integer(nrow(crosses))
  rank[ord] <- seq_len(nrow(crosses))

  cutpoints <- pmax(1L, ceiling(breaks * nrow(crosses)))
  cutpoints[[4L]] <- nrow(crosses)
  tier <- character(nrow(crosses))
  ranked_position <- rank
  tier[ranked_position <= cutpoints[[1L]]] <- labels[[1L]]
  tier[ranked_position > cutpoints[[1L]] & ranked_position <= cutpoints[[2L]]] <- labels[[2L]]
  tier[ranked_position > cutpoints[[2L]] & ranked_position <= cutpoints[[3L]]] <- labels[[3L]]
  tier[ranked_position > cutpoints[[3L]]] <- labels[[4L]]

  out <- crosses
  out$priority_index <- as.numeric(index)
  out$priority_rank <- as.integer(rank)
  out$priority_percentile <- if (nrow(out) == 1L) 1 else 1 - ((rank - 1) / (nrow(out) - 1))
  out$priority_tier <- factor(tier, levels = labels, ordered = TRUE)
  out$priority_score_component <- as.numeric(score_component)
  out$priority_kinship_component <- as.numeric(kinship_component)
  out$priority_threshold_component <- as.numeric(threshold_component)
  out$priority_check_component <- as.numeric(check_component)
  out$priority_rule <- sprintf(
    "score_col=%s; score_weight=%.4g; kinship_col=%s; kinship_weight=%.4g; threshold_col=%s; threshold_weight=%.4g; check_col=%s; check_weight=%.4g; breaks=%s",
    score_col,
    score_weight,
    if (isTRUE(has_kinship)) pair_kinship_col else "none",
    kinship_weight,
    if (isTRUE(has_threshold)) threshold_violation_col else "none",
    threshold_weight,
    if (isTRUE(has_check)) check_violation_col else "none",
    check_weight,
    paste(format(breaks, trim = TRUE), collapse = "/")
  )
  if (isTRUE(sort)) out <- out[order(out$priority_rank), , drop = FALSE]
  rownames(out) <- NULL
  out
}

ng_cross_priority_summary <- function(crosses,
                                      tier_col = "priority_tier",
                                      score_col = "multi_trait_score",
                                      kinship_col = "pair_kinship") {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  if (!(tier_col %in% names(crosses))) ng_stop("crosses missing tier_col: ", tier_col)
  tiers <- unique(as.character(crosses[[tier_col]]))
  rows <- lapply(tiers, function(tier) {
    idx <- as.character(crosses[[tier_col]]) == tier
    score <- if (score_col %in% names(crosses)) suppressWarnings(as.numeric(crosses[[score_col]][idx])) else NA_real_
    kinship <- if (kinship_col %in% names(crosses)) suppressWarnings(as.numeric(crosses[[kinship_col]][idx])) else NA_real_
    data.frame(
      priority_tier = tier,
      crosses = sum(idx),
      mean_score = mean(score, na.rm = TRUE),
      min_score = min(score, na.rm = TRUE),
      max_score = max(score, na.rm = TRUE),
      mean_pair_kinship = mean(kinship, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$priority_tier <- factor(out$priority_tier, levels = levels(crosses[[tier_col]]), ordered = TRUE)
  out <- out[order(out$priority_tier), , drop = FALSE]
  rownames(out) <- NULL
  out
}
