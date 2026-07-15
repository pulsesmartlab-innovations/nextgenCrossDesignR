ng_union_screen <- function(scores,
                            keep = 500,
                            metric_cols = c("uc_dh_gebv", "uc_dh", "mpv", "dh_pmv_var", "var_simple"),
                            min_keep_per_metric = NULL) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  metric_cols <- intersect(metric_cols, names(scores))
  if (!length(metric_cols)) ng_stop("No metric_cols found in scores")
  if (is.null(min_keep_per_metric)) min_keep_per_metric <- max(1, ceiling(keep / length(metric_cols)))
  key <- ng_pair_key(scores$parent1, scores$parent2)
  chosen <- character(0)
  for (col in metric_cols) {
    ord <- order(scores[[col]], decreasing = TRUE, na.last = NA)
    chosen <- union(chosen, key[head(ord, min_keep_per_metric)])
  }
  if (length(chosen) < keep) {
    ord <- order(scores[[metric_cols[1]]], decreasing = TRUE, na.last = NA)
    chosen <- union(chosen, key[ord])
  }
  scores[key %in% head(chosen, keep), , drop = FALSE]
}

# Fisher-z 95% CI for a Pearson correlation given n observations. Returns
# NA for n < 4 or for r with |r| >= 1 (degenerate cases).
ng_pearson_fisher_ci <- function(r, n, level = 0.95) {
  if (!is.finite(r) || !is.finite(n) || n < 4L || abs(r) >= 1 - 1e-12) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- atanh(r)
  se <- 1 / sqrt(n - 3L)
  zq <- stats::qnorm(0.5 + level / 2)
  c(lower = tanh(z - zq * se), upper = tanh(z + zq * se))
}

# Heteroscedasticity-consistent (HC3) standard error for the slope of an
# OLS simple regression y ~ x. Robust to non-constant residual variance.
ng_hc3_slope_se <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  X <- cbind(1, x)
  fit <- tryCatch(stats::lm.fit(X, y), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  resid <- fit$residuals
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NA_real_)
  h <- pmin(rowSums((X %*% XtX_inv) * X), 1 - 1e-8)
  w <- resid / (1 - h)
  meat <- crossprod(X * w)
  V <- XtX_inv %*% meat %*% XtX_inv
  sqrt(V[2L, 2L])
}

ng_validate_metric_calibration <- function(scores,
                                           realized,
                                           pred_cols = c("var_simple", "dh_recomb_var", "dh_pmv_var"),
                                           realized_var_col = "realized_var",
                                           realized_mean_col = "realized_mean",
                                           ci_level = 0.95) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  realized <- as.data.frame(realized, stringsAsFactors = FALSE)
  scores$.key <- ng_pair_key(scores$parent1, scores$parent2)
  realized$.key <- ng_pair_key(realized$parent1, realized$parent2)
  dat <- merge(scores, realized, by = ".key", suffixes = c("_pred", "_real"))
  pred_cols <- intersect(pred_cols, names(dat))
  if (!length(pred_cols)) ng_stop("No pred_cols found after merge")
  if (!(realized_var_col %in% names(dat))) ng_stop("realized variance column not found: ", realized_var_col)
  out <- lapply(pred_cols, function(col) {
    x <- dat[[col]]
    y <- dat[[realized_var_col]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) {
      return(data.frame(metric = col, n = sum(ok),
                        cor = NA_real_, cor_lower = NA_real_, cor_upper = NA_real_,
                        spearman = NA_real_,
                        slope = NA_real_, slope_se_hc3 = NA_real_,
                        intercept = NA_real_, rmse = NA_real_,
                        stringsAsFactors = FALSE))
    }
    fit <- stats::lm(y[ok] ~ x[ok])
    pearson_r <- suppressWarnings(stats::cor(x[ok], y[ok]))
    spearman_r <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    ci <- ng_pearson_fisher_ci(pearson_r, sum(ok), level = ci_level)
    slope_se <- ng_hc3_slope_se(x[ok], y[ok])
    pred <- stats::fitted(fit)
    data.frame(
      metric = col,
      n = sum(ok),
      cor = pearson_r,
      cor_lower = unname(ci["lower"]),
      cor_upper = unname(ci["upper"]),
      spearman = spearman_r,
      slope = unname(stats::coef(fit)[2]),
      slope_se_hc3 = slope_se,
      intercept = unname(stats::coef(fit)[1]),
      rmse = sqrt(mean((y[ok] - pred)^2)),
      stringsAsFactors = FALSE
    )
  })
  calib <- do.call(rbind, out)
  attr(calib, "ci_level") <- ci_level
  if (realized_mean_col %in% names(dat) && "cross_mean" %in% names(dat)) {
    ok <- is.finite(dat$cross_mean) & is.finite(dat[[realized_mean_col]])
    if (sum(ok) >= 3) {
      mean_r <- suppressWarnings(stats::cor(dat$cross_mean[ok], dat[[realized_mean_col]][ok]))
      mean_ci <- ng_pearson_fisher_ci(mean_r, sum(ok), level = ci_level)
      attr(calib, "mean_correlation") <- mean_r
      attr(calib, "mean_correlation_lower") <- unname(mean_ci["lower"])
      attr(calib, "mean_correlation_upper") <- unname(mean_ci["upper"])
    } else {
      attr(calib, "mean_correlation") <- NA_real_
    }
  }
  calib
}

ng_compare_selected_families <- function(selected_by_method,
                                         realized,
                                         family_value_cols = c("realized_mean", "realized_top10", "realized_max")) {
  realized <- as.data.frame(realized, stringsAsFactors = FALSE)
  realized$.key <- ng_pair_key(realized$parent1, realized$parent2)
  out <- lapply(names(selected_by_method), function(method) {
    sel <- selected_by_method[[method]]
    sel <- as.data.frame(sel, stringsAsFactors = FALSE)
    sel$.key <- ng_pair_key(sel$parent1, sel$parent2)
    dat <- merge(sel, realized, by = ".key", suffixes = c("_sel", "_real"))
    vals <- lapply(intersect(family_value_cols, names(dat)), function(col) {
      data.frame(metric = col, mean = mean(dat[[col]], na.rm = TRUE),
                 max = max(dat[[col]], na.rm = TRUE), stringsAsFactors = FALSE)
    })
    tab <- if (length(vals)) do.call(rbind, vals) else data.frame()
    if (!nrow(tab)) {
      tab <- data.frame(metric = NA_character_, mean = NA_real_, max = NA_real_)
    }
    tab$method <- method
    tab$n_selected <- nrow(dat)
    tab
  })
  out <- do.call(rbind, out)
  out[, c("method", "n_selected", "metric", "mean", "max")]
}
