ng_expected_topk_intensity <- function(n_progeny, top_k = NULL, top_prop = 0.10) {
  n <- as.integer(n_progeny)
  if (!is.finite(n) || n < 1) ng_stop("n_progeny must be positive")
  if (is.null(top_k)) top_k <- max(1L, ceiling(n * top_prop))
  k <- max(1L, min(n, as.integer(top_k)))
  ranks <- (n - k + 1L):n
  mean(stats::qnorm((ranks - 0.375) / (n + 0.25)))
}

ng_fit_variance_calibrator <- function(x,
                                       y,
                                       min_n = 20L,
                                       min_var = 1e-8,
                                       nonnegative_slope = TRUE,
                                       weighting = c("inv_pred_sq", "ols")) {
  # Realized within-family variance from finite progeny has chi-square noise:
  # Var(s^2_realized) is approximately proportional to (true_var)^2. OLS gives
  # high-variance observations excessive influence. Default to weighted least
  # squares with weights w_i = 1 / max(pred_i, eps)^2, which is the GLS weight
  # under the chi-square variance assumption (and unbiased for the slope).
  weighting <- match.arg(weighting)
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y) & x >= 0 & y >= 0
  x <- x[ok]
  y <- y[ok]
  if (length(x) < min_n || stats::sd(x) <= 0 || stats::sd(y) <= 0) {
    return(list(
      method = "identity",
      weighting = weighting,
      n = length(x),
      intercept = 0,
      slope = 1,
      cor = if (length(x) >= 3) suppressWarnings(stats::cor(x, y)) else NA_real_,
      rmse = if (length(x)) sqrt(mean((y - x)^2)) else NA_real_
    ))
  }

  weights <- if (identical(weighting, "inv_pred_sq")) {
    1 / pmax(x, max(min_var, stats::quantile(x[x > 0], 0.05, na.rm = TRUE)))^2
  } else {
    rep(1, length(x))
  }
  fit <- stats::lm(y ~ x, weights = weights)
  intercept <- unname(stats::coef(fit)[1])
  slope <- unname(stats::coef(fit)[2])
  if (!is.finite(intercept)) intercept <- 0
  if (!is.finite(slope)) slope <- 1
  method_used <- "wls"
  if (isTRUE(nonnegative_slope) && slope < 0) {
    # Negative slope often signals miscalibrated predictor scale rather than a
    # true negative relationship. Fall back to a non-negative-slope fit by
    # constraining slope >= 0 via a single-parameter regression through the
    # origin (slope = sum(w * x * y) / sum(w * x^2)) and re-using the
    # weighted intercept estimate.
    nnls_slope <- sum(weights * x * y) / sum(weights * x * x)
    if (!is.finite(nnls_slope) || nnls_slope < 0) nnls_slope <- 0
    intercept <- 0
    slope <- nnls_slope
    method_used <- "wls_nnls"
  }
  pred <- pmax(min_var, intercept + slope * x)
  list(
    method = method_used,
    weighting = weighting,
    n = length(x),
    intercept = intercept,
    slope = slope,
    cor = suppressWarnings(stats::cor(x, y)),
    rmse = sqrt(mean((y - pred)^2))
  )
}

ng_predict_variance_calibrator <- function(calibrator, x, min_var = 1e-8) {
  x <- as.numeric(x)
  pmax(min_var, calibrator$intercept + calibrator$slope * x)
}

ng_fit_family_variance_calibrators <- function(history,
                                               pred_map = c(
                                                 var_simple = "pred_var_simple",
                                                 dh_recomb_var = "pred_dh_recomb_var",
                                                 dh_pmv_var = "pred_dh_pmv_var"
                                               ),
                                               realized_col = "realized_var",
                                               min_n = 20L) {
  if (is.null(history) || !nrow(history) || !(realized_col %in% names(history))) {
    out <- lapply(names(pred_map), function(nm) {
      list(method = "identity", n = 0L, intercept = 0, slope = 1, cor = NA_real_, rmse = NA_real_)
    })
    names(out) <- names(pred_map)
    return(out)
  }
  out <- lapply(names(pred_map), function(score_col) {
    hist_col <- pred_map[[score_col]]
    if (!(hist_col %in% names(history))) {
      return(list(method = "identity", n = 0L, intercept = 0, slope = 1, cor = NA_real_, rmse = NA_real_))
    }
    ng_fit_variance_calibrator(history[[hist_col]], history[[realized_col]], min_n = min_n)
  })
  names(out) <- names(pred_map)
  out
}

ng_apply_family_variance_calibrators <- function(scores,
                                                 calibrators,
                                                 n_progeny,
                                                 top_k = NULL,
                                                 top_prop = 0.10,
                                                 mean_col = "cross_mean") {
  if (!(mean_col %in% names(scores))) ng_stop("scores missing mean_col: ", mean_col)
  intensity <- ng_expected_topk_intensity(n_progeny, top_k = top_k, top_prop = top_prop)
  for (score_col in names(calibrators)) {
    if (!(score_col %in% names(scores))) next
    cal_col <- paste0(score_col, "_cal")
    etk_col <- paste0("etk_", score_col, "_cal")
    scores[[cal_col]] <- ng_predict_variance_calibrator(calibrators[[score_col]], scores[[score_col]])
    scores[[etk_col]] <- scores[[mean_col]] + intensity * sqrt(pmax(scores[[cal_col]], 0))
    if ("cross_mean_gebv" %in% names(scores)) {
      etk_gebv_col <- paste0("etk_", score_col, "_gebv_cal")
      scores[[etk_gebv_col]] <- scores$cross_mean_gebv + intensity * sqrt(pmax(scores[[cal_col]], 0))
    }
    if ("cross_mean_adjusted_pheno" %in% names(scores)) {
      etk_adj_col <- paste0("etk_", score_col, "_adj_cal")
      scores[[etk_adj_col]] <- scores$cross_mean_adjusted_pheno + intensity * sqrt(pmax(scores[[cal_col]], 0))
    }
    if ("cross_mean_blend" %in% names(scores)) {
      etk_blend_col <- paste0("etk_", score_col, "_blend_cal")
      scores[[etk_blend_col]] <- scores$cross_mean_blend + intensity * sqrt(pmax(scores[[cal_col]], 0))
    }
  }
  attr(scores, "variance_calibrators") <- calibrators
  attr(scores, "topk_intensity") <- intensity
  scores
}

ng_calibrator_summary <- function(calibrators) {
  out <- lapply(names(calibrators), function(nm) {
    c <- calibrators[[nm]]
    data.frame(
      variance_metric = nm,
      calibration_method = c$method,
      calibration_n = c$n,
      calibration_intercept = c$intercept,
      calibration_slope = c$slope,
      calibration_cor = c$cor,
      calibration_rmse = c$rmse,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

ng_rank_normal_score <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) < 2L) {
    out[ok] <- 0
    return(out)
  }
  r <- rank(x[ok], ties.method = "average", na.last = "keep")
  p <- (r - 0.5) / length(r)
  p <- pmin(1 - 1e-6, pmax(1e-6, p))
  z <- stats::qnorm(p)
  z <- ng_standardize(z)
  out[ok] <- z
  out
}

ng_default_portfolio_weights <- function(scores, score_cols) {
  rel <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))

  mean_w <- c(
    gebv = max(0.05, rel),
    blend = 0.60,
    adj = max(0.05, 1 - rel)
  )
  mean_w <- mean_w / sum(mean_w)
  # Two correct-units variance choices: VPM (recomb) and PMV (pmv).
  var_w <- c(
    recomb = 0.45 + 0.10 * rel,
    pmv    = 0.55 - 0.10 * rel
  )
  var_w <- var_w / sum(var_w)

  out <- numeric(length(score_cols))
  names(out) <- score_cols
  for (col in score_cols) {
    m <- if (grepl("_gebv_", col, fixed = TRUE)) {
      "gebv"
    } else if (grepl("_adj_", col, fixed = TRUE)) {
      "adj"
    } else {
      "blend"
    }
    v <- if (grepl("dh_pmv_var", col, fixed = TRUE)) "pmv" else "recomb"
    out[col] <- mean_w[[m]] * var_w[[v]]
  }
  out / sum(out)
}

ng_history_portfolio_weights <- function(history,
                                         score_cols,
                                         target_col = "realized_top10",
                                         min_n = 20L) {
  if (is.null(history) || !nrow(history) || !(target_col %in% names(history))) {
    return(rep(NA_real_, length(score_cols)))
  }
  weights <- numeric(length(score_cols))
  names(weights) <- score_cols
  for (col in score_cols) {
    hist_col <- paste0("pred_", col)
    if (!(hist_col %in% names(history))) {
      weights[col] <- NA_real_
      next
    }
    x <- as.numeric(history[[hist_col]])
    y <- as.numeric(history[[target_col]])
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < min_n || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) {
      weights[col] <- NA_real_
      next
    }
    rho <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    weights[col] <- if (is.finite(rho)) pmax(0, rho)^2 * sqrt(sum(ok) / (sum(ok) + 25)) else NA_real_
  }
  weights
}

ng_add_local_portfolio_scores <- function(scores,
                                          history = NULL,
                                          target_col = "realized_top10",
                                          min_history_n = 20L,
                                          score_cols = c(
                                            "etk_dh_recomb_var_gebv_cal",
                                            "etk_dh_recomb_var_adj_cal",
                                            "etk_dh_recomb_var_blend_cal",
                                            "etk_dh_pmv_var_gebv_cal",
                                            "etk_dh_pmv_var_adj_cal",
                                            "etk_dh_pmv_var_blend_cal"
                                          ),
                                          history_weight = NULL) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  score_cols <- intersect(score_cols, names(scores))
  if (!length(score_cols)) {
    scores$ng_portfolio_score <- NA_real_
    return(scores)
  }

  default_w <- ng_default_portfolio_weights(scores, score_cols)
  hist_w <- ng_history_portfolio_weights(history, score_cols, target_col = target_col, min_n = min_history_n)
  if (any(is.finite(hist_w) & hist_w > 0)) {
    hist_w[!is.finite(hist_w) | hist_w < 0] <- 0
    hist_w <- hist_w / sum(hist_w)
    hist_n <- if (!is.null(history) && target_col %in% names(history)) {
      sum(is.finite(as.numeric(history[[target_col]])))
    } else {
      0L
    }
    if (is.null(history_weight)) {
      history_weight <- min(0.80, hist_n / (hist_n + 60))
    }
    history_weight <- max(0, min(0.80, history_weight))
    weights <- (1 - history_weight) * default_w[names(hist_w)] + history_weight * hist_w
  } else {
    weights <- default_w
    history_weight <- 0
  }
  weights <- weights / sum(weights)

  z <- lapply(names(weights), function(col) ng_rank_normal_score(scores[[col]]))
  z <- do.call(cbind, z)
  consensus <- as.numeric(z %*% weights)
  champion <- apply(z, 1, max, na.rm = TRUE)
  champion[!is.finite(champion)] <- consensus[!is.finite(champion)]
  raw <- 0.70 * consensus + 0.30 * champion
  raw <- ng_standardize(raw)
  fallback <- rowMeans(z, na.rm = TRUE)
  raw[!is.finite(raw)] <- fallback[!is.finite(raw)]
  var_terms <- matrix(0, nrow = nrow(scores), ncol = length(weights))
  colnames(var_terms) <- names(weights)
  for (col in names(weights)) {
    var_col <- if (grepl("dh_pmv_var", col, fixed = TRUE)) {
      "dh_pmv_var_cal"
    } else {
      "dh_recomb_var_cal"
    }
    if (var_col %in% names(scores)) {
      var_terms[, col] <- pmax(as.numeric(scores[[var_col]]), 0)
    }
  }
  scores$ng_portfolio_score <- raw
  scores$ng_portfolio_var <- as.numeric(var_terms %*% weights)
  scores$ng_portfolio_reliability <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  scores$ng_portfolio_history_weight <- history_weight
  attr(scores, "local_portfolio_weights") <- weights
  scores
}

ng_adaptive_stack_parent_use <- function(n_parents,
                                         n_crosses,
                                         reliability,
                                         min_value = 0.75,
                                         max_value = 2.25) {
  n_parents <- as.numeric(n_parents)
  n_crosses <- as.numeric(n_crosses)
  reliability <- as.numeric(reliability)
  if (!is.finite(n_parents) || n_parents <= 0) n_parents <- 2 * n_crosses
  if (!is.finite(n_crosses) || n_crosses <= 0) n_crosses <- 10
  if (!is.finite(reliability)) reliability <- 0.35
  reliability <- max(0, min(1, reliability))
  parent_pressure <- max(0, min(1, n_crosses / n_parents))
  value <- 0.75 + 2.0 * parent_pressure + 0.5 * (1 - reliability)
  max(min_value, min(max_value, value))
}

ng_adaptive_stack_fallback_weight <- function(reliability,
                                              history_n = 0L,
                                              min_history_n = 20L,
                                              low_reliability = 0.55,
                                              very_low_reliability = 0.30,
                                              max_weight = 0.35) {
  reliability <- as.numeric(reliability)
  if (!is.finite(reliability)) reliability <- 0.35
  reliability <- max(0, min(1, reliability))
  history_n <- as.numeric(history_n)
  if (!is.finite(history_n) || history_n < 0) history_n <- 0
  min_history_n <- max(1, as.numeric(min_history_n))
  low_reliability <- max(very_low_reliability + 1e-6, as.numeric(low_reliability))
  max_weight <- max(0, min(0.75, as.numeric(max_weight)))

  reliability_gap <- max(0, (low_reliability - reliability) /
                           (low_reliability - very_low_reliability))
  history_gap <- max(0, 1 - history_n / min_history_n)
  value <- max_weight * reliability_gap * (0.70 + 0.30 * history_gap)
  max(0, min(max_weight, value))
}

ng_adaptive_stack_prior_weights <- function(scores,
                                            score_cols,
                                            n_parents = NULL,
                                            n_crosses = NULL) {
  rel <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))
  if (is.null(n_parents)) n_parents <- length(unique(c(scores$parent1, scores$parent2)))
  if (is.null(n_crosses)) n_crosses <- max(1L, ceiling(0.10 * n_parents))
  pressure <- max(0, min(1, as.numeric(n_crosses) / max(1, as.numeric(n_parents))))

  out <- numeric(length(score_cols))
  names(out) <- score_cols
  fallback_boost <- ng_adaptive_stack_fallback_weight(
    reliability = rel,
    history_n = 0L,
    min_history_n = 20L,
    max_weight = 0.35
  )
  for (col in score_cols) {
    out[col] <- if (col %in% c("uc_recomb_gebv", "etk_dh_recomb_var_gebv_cal")) {
      0.55 + 0.25 * rel
    } else if (grepl("recomb_var_blend", col, fixed = TRUE) ||
               grepl("uc_recomb_blend", col, fixed = TRUE)) {
      0.22 + 0.18 * (1 - rel)
    } else if (grepl("dh_pmv_var", col, fixed = TRUE) ||
               col %in% c("uc_dh_gebv", "uc_dh_blend")) {
      0.20 + 0.15 * rel
    } else if (grepl("var_simple", col, fixed = TRUE)) {
      0.05 + 0.20 * (1 - rel) + 0.15 * pressure + fallback_boost
    } else if (grepl("portfolio", col, fixed = TRUE)) {
      0.12
    } else {
      0.08
    }
  }
  out[!is.finite(out) | out < 0] <- 0
  if (!any(out > 0)) out[] <- 1
  out / sum(out)
}

ng_variance_col_for_score <- function(score_col) {
  if (grepl("dh_pmv_var", score_col, fixed = TRUE) ||
      score_col %in% c("uc_dh_gebv", "uc_dh_blend", "uc_dh_adj", "uc_dh")) {
    "dh_pmv_var_cal"
  } else if (grepl("var_simple", score_col, fixed = TRUE)) {
    "var_simple_cal"
  } else if (grepl("portfolio", score_col, fixed = TRUE)) {
    "ng_portfolio_var"
  } else {
    "dh_recomb_var_cal"
  }
}

ng_add_adaptive_stack_scores <- function(scores,
                                         history = NULL,
                                         target_col = "realized_top10",
                                         min_history_n = 20L,
                                         score_cols = c(
                                           "uc_recomb_gebv",
                                           "etk_dh_recomb_var_gebv_cal",
                                           "etk_dh_recomb_var_blend_cal",
                                           "etk_dh_pmv_var_gebv_cal",
                                           "etk_dh_pmv_var_blend_cal",
                                           "etk_var_simple_cal"
                                         ),
                                         n_parents = NULL,
                                         n_crosses = NULL,
                                         history_weight = NULL,
                                         champion_weight = NULL,
                                         fallback_col = "etk_var_simple_cal",
                                         fallback_weight = NULL,
                                         fallback_max_weight = 0.35) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  score_cols <- intersect(unique(score_cols), names(scores))
  if (!length(score_cols)) {
    scores$ng_adaptive_score <- NA_real_
    scores$ng_adaptive_var <- NA_real_
    return(scores)
  }
  if (is.null(n_parents)) n_parents <- length(unique(c(scores$parent1, scores$parent2)))
  if (is.null(n_crosses)) n_crosses <- max(1L, ceiling(0.10 * n_parents))
  rel <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))

  prior_w <- ng_adaptive_stack_prior_weights(
    scores = scores,
    score_cols = score_cols,
    n_parents = n_parents,
    n_crosses = n_crosses
  )
  hist_w <- ng_history_portfolio_weights(
    history = history,
    score_cols = score_cols,
    target_col = target_col,
    min_n = min_history_n
  )
  hist_n <- if (!is.null(history) && target_col %in% names(history)) {
    sum(is.finite(as.numeric(history[[target_col]])))
  } else {
    0L
  }
  if (any(is.finite(hist_w) & hist_w > 0)) {
    hist_w[!is.finite(hist_w) | hist_w < 0] <- 0
    hist_w <- hist_w / sum(hist_w)
    if (is.null(history_weight) || !is.finite(history_weight)) {
      history_weight <- min(0.65, hist_n / (hist_n + 80))
    }
    history_weight <- max(0, min(0.65, history_weight))
    weights <- (1 - history_weight) * prior_w[names(hist_w)] + history_weight * hist_w
  } else {
    weights <- prior_w
    history_weight <- 0
  }
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) weights[] <- 1
  weights <- weights / sum(weights)

  z <- lapply(names(weights), function(col) ng_rank_normal_score(scores[[col]]))
  z <- do.call(cbind, z)
  colnames(z) <- names(weights)
  consensus <- as.numeric(z %*% weights)
  champion <- apply(z, 1, max, na.rm = TRUE)
  champion[!is.finite(champion)] <- consensus[!is.finite(champion)]
  if (is.null(champion_weight) || !is.finite(champion_weight)) {
    champion_weight <- 0.10 + 0.10 * rel
  }
  champion_weight <- max(0, min(0.30, champion_weight))
  raw <- (1 - champion_weight) * consensus + champion_weight * champion
  raw <- ng_standardize(raw)
  fallback <- rowMeans(z, na.rm = TRUE)
  raw[!is.finite(raw)] <- fallback[!is.finite(raw)]

  guarded <- raw
  fallback_col <- as.character(fallback_col[[1]])
  if (nzchar(fallback_col) && fallback_col %in% names(scores)) {
    if (is.null(fallback_weight) || !is.finite(fallback_weight)) {
      fallback_weight <- ng_adaptive_stack_fallback_weight(
        reliability = rel,
        history_n = hist_n,
        min_history_n = min_history_n,
        max_weight = fallback_max_weight
      )
    }
    fallback_weight <- max(0, min(fallback_max_weight, as.numeric(fallback_weight)))
    fallback_z <- ng_rank_normal_score(scores[[fallback_col]])
    ok <- is.finite(fallback_z)
    guarded[ok] <- (1 - fallback_weight) * raw[ok] + fallback_weight * fallback_z[ok]
    guarded <- ng_standardize(guarded)
  } else {
    fallback_weight <- 0
  }

  var_terms <- matrix(0, nrow = nrow(scores), ncol = length(weights))
  colnames(var_terms) <- names(weights)
  for (col in names(weights)) {
    var_col <- ng_variance_col_for_score(col)
    if (var_col %in% names(scores)) {
      var_terms[, col] <- pmax(as.numeric(scores[[var_col]]), 0)
    }
  }
  scores$ng_adaptive_score_raw <- raw
  scores$ng_adaptive_score <- guarded
  scores$ng_adaptive_var <- as.numeric(var_terms %*% weights)
  scores$ng_adaptive_reliability <- rel
  scores$ng_adaptive_history_weight <- history_weight
  scores$ng_adaptive_fallback_weight <- fallback_weight
  scores$ng_adaptive_fallback_col <- fallback_col
  scores$ng_adaptive_parent_use_input <- ng_adaptive_stack_parent_use(
    n_parents = n_parents,
    n_crosses = n_crosses,
    reliability = rel
  )
  attr(scores, "adaptive_stack") <- list(
    weights = weights,
    score_cols = names(weights),
    reliability = rel,
    history_weight = history_weight,
    history_n = hist_n,
    champion_weight = champion_weight,
    fallback_col = fallback_col,
    fallback_weight = fallback_weight,
    n_parents = n_parents,
    n_crosses = n_crosses,
    parent_use_input = scores$ng_adaptive_parent_use_input[[1]]
  )
  scores
}

ng_meta_portfolio_prior_weights <- function(scores,
                                            score_cols,
                                            n_parents = NULL,
                                            n_crosses = NULL) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  rel <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))
  if (is.null(n_parents)) n_parents <- length(unique(c(scores$parent1, scores$parent2)))
  if (is.null(n_crosses)) n_crosses <- max(1L, ceiling(0.10 * n_parents))
  n_parents <- max(2, as.numeric(n_parents))
  n_crosses <- max(1, as.numeric(n_crosses))
  pressure <- max(0, min(1, n_crosses / n_parents))
  middle_parent <- exp(-((n_parents - 50) / 22)^2)
  large_parent <- 1 / (1 + exp(-(n_parents - 65) / 8))

  out <- numeric(length(score_cols))
  names(out) <- score_cols
  for (col in score_cols) {
    out[col] <- if (identical(col, "ng_adaptive_score")) {
      0.28 + 0.12 * pressure + 0.10 * (1 - rel)
    } else if (col %in% c("uc_recomb_gebv", "etk_dh_recomb_var_gebv_cal")) {
      0.22 + 0.30 * rel + 0.12 * large_parent
    } else if (grepl("dh_pmv_var", col, fixed = TRUE) ||
               col %in% c("uc_dh_gebv", "uc_dh_blend")) {
      0.18 + 0.18 * rel + 0.20 * pressure
    } else if (col %in% c("popvar_uc", "simple_usefa")) {
      0.18 + 0.12 * rel + 0.18 * middle_parent
    } else if (grepl("var_simple", col, fixed = TRUE)) {
      0.10 + 0.35 * (1 - rel)
    } else {
      0.08
    }
  }
  out[!is.finite(out) | out < 0] <- 0
  if (!any(out > 0)) out[] <- 1
  out / sum(out)
}

ng_meta_weighted_rank_consensus <- function(z, weights) {
  z <- as.matrix(z)
  weights <- as.numeric(weights)
  names(weights) <- colnames(z)
  ok <- is.finite(z)
  zw <- z
  zw[!ok] <- 0
  denom <- as.numeric(ok %*% weights)
  out <- as.numeric(zw %*% weights)
  good <- denom > 0
  out[good] <- out[good] / denom[good]
  out[!good] <- NA_real_
  out
}

ng_meta_score_history_columns <- function(score_cols) {
  setNames(paste0("pred_", score_cols), score_cols)
}

ng_meta_method_patterns <- function(score_cols) {
  out <- lapply(score_cols, function(col) {
    if (identical(col, "ng_adaptive_score")) {
      # Adaptive score is already a composite of base metrics. Let cross-level
      # calibration weight it, but do not let method history self-reinforce it.
      NA_character_
    } else if (col %in% c("uc_recomb_gebv", "etk_dh_recomb_var_gebv_cal")) {
      "^ng_recomb_gebv_ocs"
    } else if (grepl("dh_pmv_var", col, fixed = TRUE) ||
               col %in% c("uc_dh_gebv", "uc_dh_blend")) {
      "^ng_pmv_blend_balanced_ocs|^ng_pmv_gebv_ocs"
    } else if (grepl("var_simple", col, fixed = TRUE)) {
      "^var_simple_ocs|^var_simple_select|^var_simple_allocator|^var_simple_mip|^var_simple_repair"
    } else if (identical(col, "popvar_uc")) {
      "^popvar_uc_ocs|^popvar_uc_select"
    } else if (identical(col, "simple_usefa")) {
      "^simple_usefa_ocs|^simple_usefa_select"
    } else {
      paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", col))
    }
  })
  names(out) <- score_cols
  out
}

ng_meta_method_history_weights <- function(history,
                                           score_cols,
                                           target_col = "realized_top10",
                                           min_n = 10L,
                                           recent_cycles = 3L,
                                           temperature = 1.25) {
  out <- rep(NA_real_, length(score_cols))
  names(out) <- score_cols
  attr(out, "n") <- 0L
  if (is.null(history) || !nrow(history) ||
      !(target_col %in% names(history)) ||
      !("method" %in% names(history))) {
    return(out)
  }
  h <- as.data.frame(history, stringsAsFactors = FALSE)
  if ("cycle" %in% names(h)) {
    cyc <- suppressWarnings(as.integer(h$cycle))
    if (any(is.finite(cyc))) {
      keep_cycles <- sort(unique(cyc[is.finite(cyc)]), decreasing = TRUE)
      keep_cycles <- keep_cycles[seq_len(min(length(keep_cycles), max(1L, as.integer(recent_cycles))))]
      h <- h[cyc %in% keep_cycles, , drop = FALSE]
    }
  }
  patterns <- ng_meta_method_patterns(score_cols)
  perf <- rep(NA_real_, length(score_cols))
  names(perf) <- score_cols
  y <- as.numeric(h[[target_col]])
  method_chr <- as.character(h$method)
  usable_patterns <- patterns[vapply(patterns, function(pattern) {
    is.character(pattern) && nzchar(pattern) && !is.na(pattern)
  }, logical(1))]
  hit_map <- if (length(usable_patterns)) {
    lapply(usable_patterns, function(pattern) {
      grepl(pattern, method_chr, perl = TRUE)
    })
  } else {
    list()
  }
  any_hit <- if (length(hit_map)) Reduce(`|`, hit_map) else rep(FALSE, nrow(h))
  n_by_col <- rep(0L, length(score_cols))
  names(n_by_col) <- score_cols
  for (col in names(hit_map)) {
    n_by_col[col] <- sum(hit_map[[col]] & is.finite(y))
  }
  attr(out, "n") <- sum(any_hit & is.finite(y))
  if (!length(hit_map) || !any(n_by_col >= min_n)) return(out)

  cycle_id <- if ("cycle" %in% names(h)) {
    cyc <- suppressWarnings(as.integer(h$cycle))
    ifelse(is.finite(cyc), cyc, 1L)
  } else {
    rep(1L, nrow(h))
  }
  cycle_levels <- sort(unique(cycle_id[is.finite(cycle_id)]))
  if (!length(cycle_levels)) cycle_levels <- 1L
  cycle_scores <- matrix(NA_real_, nrow = length(cycle_levels), ncol = length(score_cols))
  colnames(cycle_scores) <- score_cols
  rownames(cycle_scores) <- as.character(cycle_levels)
  for (i in seq_along(cycle_levels)) {
    in_cycle <- cycle_id == cycle_levels[[i]]
    for (col in names(hit_map)) {
      if (n_by_col[[col]] < min_n) next
      ok <- in_cycle & hit_map[[col]] & is.finite(y)
      if (any(ok)) cycle_scores[i, col] <- mean(y[ok], na.rm = TRUE)
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
  for (col in score_cols) {
    ok <- is.finite(cycle_scores[, col])
    if (any(ok)) {
      perf[col] <- stats::weighted.mean(cycle_scores[ok, col], recency_w[ok], na.rm = TRUE)
    }
  }
  ok <- is.finite(perf)
  if (sum(ok) < 2L || stats::sd(perf[ok]) <= 0) return(out)
  z <- ng_standardize(perf[ok])
  w <- exp(max(0.05, as.numeric(temperature)) * z)
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) return(out)
  out[ok] <- w / sum(w)
  out
}

ng_add_meta_portfolio_scores <- function(scores,
                                         history = NULL,
                                         method_history = history,
                                         target_col = "realized_top10",
                                         min_history_n = 20L,
                                         method_min_history_n = 10L,
                                         method_recent_cycles = 3L,
                                         score_cols = c(
                                           "ng_adaptive_score",
                                           "uc_recomb_gebv",
                                           "etk_dh_pmv_var_blend_cal",
                                           "etk_var_simple_cal",
                                           "popvar_uc",
                                           "simple_usefa"
                                         ),
                                         n_parents = NULL,
                                         n_crosses = NULL,
                                         history_weight = NULL,
                                         champion_weight = NULL) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  score_cols <- intersect(unique(score_cols), names(scores))
  score_cols <- score_cols[vapply(score_cols, function(col) {
    any(is.finite(as.numeric(scores[[col]])))
  }, logical(1))]
  if (!length(score_cols)) {
    scores$ng_meta_score <- NA_real_
    scores$ng_meta_var <- NA_real_
    return(scores)
  }
  if (is.null(n_parents)) n_parents <- length(unique(c(scores$parent1, scores$parent2)))
  if (is.null(n_crosses)) n_crosses <- max(1L, ceiling(0.10 * n_parents))
  rel <- suppressWarnings(stats::median(as.numeric(scores$effect_reliability), na.rm = TRUE))
  if (!is.finite(rel)) rel <- 0.35
  rel <- max(0, min(1, rel))

  prior_w <- ng_meta_portfolio_prior_weights(scores, score_cols, n_parents, n_crosses)
  hist_w <- ng_history_portfolio_weights(
    history = history,
    score_cols = score_cols,
    target_col = target_col,
    min_n = min_history_n
  )
  method_w <- ng_meta_method_history_weights(
    history = method_history,
    score_cols = score_cols,
    target_col = target_col,
    min_n = min(max(1L, as.integer(method_min_history_n)), max(3L, as.integer(n_crosses))),
    recent_cycles = method_recent_cycles
  )
  hist_n <- if (!is.null(history) && target_col %in% names(history)) {
    sum(is.finite(as.numeric(history[[target_col]])))
  } else {
    0L
  }
  method_hist_n <- attr(method_w, "n")
  if (!is.finite(method_hist_n)) method_hist_n <- 0L
  if (any(is.finite(hist_w) & hist_w > 0)) {
    hist_w[!is.finite(hist_w) | hist_w < 0] <- 0
    hist_w <- hist_w / sum(hist_w)
    if (is.null(history_weight) || !is.finite(history_weight)) {
      history_weight <- min(0.40, hist_n / (hist_n + 120))
    }
    history_weight <- max(0, min(0.40, history_weight))
    weights <- (1 - history_weight) * prior_w[names(hist_w)] + history_weight * hist_w
  } else {
    weights <- prior_w
    history_weight <- 0
  }
  method_history_weight <- 0
  if (any(is.finite(method_w) & method_w > 0)) {
    method_w[!is.finite(method_w) | method_w < 0] <- 0
    method_w <- method_w / sum(method_w)
    method_history_weight <- min(0.65, method_hist_n / (method_hist_n + 80))
    method_history_weight <- max(0, min(0.65, method_history_weight))
    weights <- (1 - method_history_weight) * weights[names(method_w)] +
      method_history_weight * method_w
  }
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) weights[] <- 1
  weights <- weights / sum(weights)
  leader_weights <- if (any(is.finite(method_w) & method_w > 0)) method_w else weights
  leader_weights[!is.finite(leader_weights) | leader_weights < 0] <- 0
  if (!any(leader_weights > 0)) leader_weights <- weights
  leader_col <- names(which.max(leader_weights))
  leader_weight <- if (length(leader_col)) leader_weights[[leader_col]] else NA_real_

  z <- lapply(names(weights), function(col) ng_rank_normal_score(scores[[col]]))
  z <- do.call(cbind, z)
  colnames(z) <- names(weights)
  consensus <- ng_meta_weighted_rank_consensus(z, weights)
  champion <- apply(z, 1, max, na.rm = TRUE)
  champion[!is.finite(champion)] <- consensus[!is.finite(champion)]
  if (is.null(champion_weight) || !is.finite(champion_weight)) {
    winner <- if (length(weights)) max(weights) else 0
    champion_weight <- 0.10 + 0.20 * max(0, winner - (1 / length(weights)))
  }
  champion_weight <- max(0, min(0.30, champion_weight))
  raw <- (1 - champion_weight) * consensus + champion_weight * champion
  raw <- ng_standardize(raw)
  fallback <- rowMeans(z, na.rm = TRUE)
  raw[!is.finite(raw)] <- fallback[!is.finite(raw)]

  var_terms <- matrix(0, nrow = nrow(scores), ncol = length(weights))
  colnames(var_terms) <- names(weights)
  for (col in names(weights)) {
    var_col <- if (identical(col, "popvar_uc")) {
      "popvar_varG"
    } else if (identical(col, "simple_usefa")) {
      "simple_usefa_var"
    } else if (identical(col, "ng_adaptive_score")) {
      "ng_adaptive_var"
    } else {
      ng_variance_col_for_score(col)
    }
    if (var_col %in% names(scores)) {
      var_terms[, col] <- pmax(as.numeric(scores[[var_col]]), 0)
    }
  }
  scores$ng_meta_score <- raw
  scores$ng_meta_var <- as.numeric(var_terms %*% weights)
  scores$ng_meta_reliability <- rel
  scores$ng_meta_history_weight <- history_weight
  scores$ng_meta_method_history_weight <- method_history_weight
  scores$ng_meta_leader_col <- leader_col
  scores$ng_meta_leader_weight <- leader_weight
  scores$ng_meta_leader_score <- if (!is.null(leader_col) && leader_col %in% names(scores)) {
    as.numeric(scores[[leader_col]])
  } else {
    scores$ng_meta_score
  }
  leader_var_col <- if (identical(leader_col, "popvar_uc")) {
    "popvar_varG"
  } else if (identical(leader_col, "simple_usefa")) {
    "simple_usefa_var"
  } else if (identical(leader_col, "ng_adaptive_score")) {
    "ng_adaptive_var"
  } else {
    ng_variance_col_for_score(leader_col)
  }
  scores$ng_meta_leader_var <- if (!is.null(leader_var_col) && leader_var_col %in% names(scores)) {
    pmax(as.numeric(scores[[leader_var_col]]), 0)
  } else {
    scores$ng_meta_var
  }
  scores$ng_meta_champion_weight <- champion_weight
  scores$ng_meta_parent_use_input <- ng_adaptive_stack_parent_use(
    n_parents = n_parents,
    n_crosses = n_crosses,
    reliability = rel
  )
  attr(scores, "meta_portfolio") <- list(
    weights = weights,
    score_cols = names(weights),
    reliability = rel,
    history_weight = history_weight,
    history_n = hist_n,
    method_history_weight = method_history_weight,
    method_history_n = method_hist_n,
    method_weights = method_w,
    leader_col = leader_col,
    leader_weight = leader_weight,
    leader_var_col = leader_var_col,
    champion_weight = champion_weight,
    n_parents = n_parents,
    n_crosses = n_crosses,
    parent_use_input = scores$ng_meta_parent_use_input[[1]]
  )
  scores
}
