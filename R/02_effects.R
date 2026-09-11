ng_fit_ridge_effects <- function(geno,
                                 y,
                                 ids = rownames(geno),
                                 lambda = NULL,
                                 lambda_grid = 10 ^ seq(-2, 5, length.out = 20),
                                 h2_prior = NULL,
                                 kfold = 5,
                                 seed = 1L,
                                 return_beta_cov_full = FALSE) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  y <- ng_match_vector(y, ids, "y")
  ok <- is.finite(y)
  if (sum(ok) < 2L) ng_stop("ridge fitting requires at least two finite phenotypes")
  X_raw <- geno[ok, , drop = FALSE]
  y <- y[ok]
  marker_mean <- colMeans(X_raw, na.rm = TRUE)
  marker_mean[!is.finite(marker_mean)] <- 0
  X <- sweep(X_raw, 2, marker_mean, "-")
  X[!is.finite(X)] <- 0
  yc <- y - mean(y)

  lambda_grid <- suppressWarnings(as.numeric(lambda_grid))
  lambda_grid <- unique(lambda_grid[is.finite(lambda_grid) & lambda_grid > 0])
  if (!length(lambda_grid)) ng_stop("lambda_grid must contain at least one finite positive value")
  if (is.null(lambda)) {
    if (!is.null(h2_prior)) {
      h2_prior <- max(min(as.numeric(h2_prior), 0.99), 0.01)
      lambda0 <- ncol(X) * (1 - h2_prior) / h2_prior
      lambda_grid <- unique(sort(c(lambda0, lambda0 * 10 ^ seq(-2, 2, length.out = 17))))
    }
    lambda <- ng_choose_ridge_lambda(X_raw, y, lambda_grid, kfold = kfold, seed = seed)
  }
  lambda <- suppressWarnings(as.numeric(lambda))
  if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    ng_stop("lambda must be one finite positive value")
  }

  A <- tcrossprod(X) + diag(lambda, nrow(X))
  alpha <- as.numeric(solve(A, yc))
  beta <- as.numeric(crossprod(X, alpha))
  names(beta) <- colnames(geno)
  fitted <- mean(y) + as.numeric(X %*% beta)
  resid <- y - fitted
  # Posterior marker-effect variance for ridge in dual form via Woodbury:
  #   (X'X + lambda I)^{-1} = (1/lambda) (I - X' A^{-1} X),  A = XX' + lambda I.
  # So diag((X'X + lambda I)^{-1})_k = (1/lambda) (1 - x_k' A^{-1} x_k).
  AinvX <- solve(A, X)
  diag_xtAinvx <- pmax(colSums(X * AinvX), 0)
  df_eff <- sum(diag_xtAinvx)
  # Effective residual df for the centered ridge smoother plus one exact,
  # unpenalized intercept: n - 1 - 2 tr(H) + tr(H'H).
  H <- tcrossprod(X, AinvX)
  df_residual <- length(y) - 1 - 2 * df_eff + sum(H * H)
  sigma_e2 <- sum(resid^2) / max(1, df_residual)
  beta_var <- sigma_e2 * pmax(1 - diag_xtAinvx, 0) / lambda
  names(beta_var) <- names(beta)
  # Full posterior covariance of marker effects, opt-in because the m x m dense
  # matrix is O(m^2) memory (~300 MB at m = 6000). Same dual Woodbury identity
  # used for the diagonal: (X'X + lambda I)^{-1} = (1/lambda) (I_m - X' A^{-1} X),
  # so Sigma_beta = (sigma_e2 / lambda) (I_m - X' AinvX). Symmetrize defensively
  # against floating-point asymmetry and clip the diagonal at 0 so it is never
  # less than zero on a finite-precision machine.
  beta_cov_full <- NULL
  if (isTRUE(return_beta_cov_full)) {
    m <- ncol(X)
    if (m > 6000L) {
      warning(sprintf(
        "return_beta_cov_full = TRUE on %d markers; dense Sigma_beta is ~%.2f GB. Consider LD pruning before requesting the full posterior covariance.",
        m, (m * m * 8) / (1024 ^ 3)
      ), call. = FALSE)
    }
    XtAinvX <- crossprod(X, AinvX)
    Sigma <- (-XtAinvX)
    diag(Sigma) <- diag(Sigma) + 1
    Sigma <- 0.5 * (Sigma + t(Sigma))
    Sigma <- Sigma * (sigma_e2 / lambda)
    diag(Sigma) <- pmax(diag(Sigma), 0)
    rownames(Sigma) <- names(beta)
    colnames(Sigma) <- names(beta)
    beta_cov_full <- Sigma
  }
  sst <- sum((y - mean(y))^2)
  in_sample_predictive_r2 <- if (is.finite(sst) && sst > 0) {
    1 - sum((y - fitted)^2) / sst
  } else {
    NA_real_
  }
  in_sample_predictive_correlation <- suppressWarnings(
    stats::cor(y, fitted, use = "complete.obs")
  )
  cv <- ng_ridge_cv_predict(X_raw, y, lambda = lambda, kfold = kfold, seed = seed)
  cv_predictive_r2 <- if (!is.null(cv)) {
    keep_cv <- is.finite(y) & is.finite(cv)
    sst_cv <- sum((y[keep_cv] - mean(y[keep_cv]))^2)
    if (sum(keep_cv) >= 2L && is.finite(sst_cv) && sst_cv > 0) {
      1 - sum((y[keep_cv] - cv[keep_cv])^2) / sst_cv
    } else {
      NA_real_
    }
  } else NA_real_
  cv_predictive_correlation <- if (!is.null(cv)) suppressWarnings(
    stats::cor(y, cv, use = "complete.obs")
  ) else NA_real_

  list(
    method = "ridge_dual",
    beta = beta,
    beta_var = beta_var,
    beta_cov_full = beta_cov_full,
    intercept = mean(y) - sum(marker_mean * beta),
    marker_mean = marker_mean,
    lambda = lambda,
    sigma_e2 = sigma_e2,
    fitted_ids = ids[ok],
    fitted = setNames(fitted, ids[ok]),
    # Legacy names are deliberately missing: phenotype predictive diagnostics
    # are not accuracy^2 against true breeding value and not PEV reliability.
    reliability = NA_real_,
    in_sample_reliability = NA_real_,
    reliability_is_calibrated = FALSE,
    cv_predictive_r2 = cv_predictive_r2,
    cv_predictive_correlation = cv_predictive_correlation,
    in_sample_predictive_r2 = in_sample_predictive_r2,
    in_sample_predictive_correlation = in_sample_predictive_correlation,
    residual_df = df_residual
  )
}

ng_choose_ridge_lambda <- function(X, y, lambda_grid, kfold = 5, seed = 1L) {
  n <- nrow(X)
  if (n < 10 || kfold < 2) return(stats::median(lambda_grid))
  kfold <- min(as.integer(kfold), n)
  set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))
  pred <- matrix(NA_real_, nrow = n, ncol = length(lambda_grid))
  for (fold in seq_len(kfold)) {
    tr <- folds != fold
    te <- !tr
    marker_mean <- colMeans(X[tr, , drop = FALSE], na.rm = TRUE)
    marker_mean[!is.finite(marker_mean)] <- 0
    Xtr <- sweep(X[tr, , drop = FALSE], 2L, marker_mean, "-")
    Xte <- sweep(X[te, , drop = FALSE], 2L, marker_mean, "-")
    Xtr[!is.finite(Xtr)] <- 0
    Xte[!is.finite(Xte)] <- 0
    y_mean <- mean(y[tr])
    ytr <- y[tr] - y_mean
    eig <- eigen(tcrossprod(Xtr), symmetric = TRUE)
    vals <- pmax(eig$values, 0)
    Uy <- as.numeric(crossprod(eig$vectors, ytr))
    KteU <- tcrossprod(Xte, Xtr) %*% eig$vectors
    for (g in seq_along(lambda_grid)) {
      lam <- lambda_grid[g]
      pred[te, g] <- y_mean + as.numeric(KteU %*% (Uy / (vals + lam)))
    }
  }
  rmse <- sqrt(colMeans((y - pred)^2, na.rm = TRUE))
  lambda_grid[which.min(rmse)]
}

ng_ridge_cv_predict <- function(X, y, lambda, kfold = 5, seed = 1L) {
  n <- nrow(X)
  if (n < 10 || kfold < 2) return(NULL)
  kfold <- min(kfold, n)
  set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))
  pred <- rep(NA_real_, n)
  for (fold in seq_len(kfold)) {
    tr <- folds != fold
    te <- !tr
    marker_mean <- colMeans(X[tr, , drop = FALSE], na.rm = TRUE)
    marker_mean[!is.finite(marker_mean)] <- 0
    Xtr <- sweep(X[tr, , drop = FALSE], 2L, marker_mean, "-")
    Xte <- sweep(X[te, , drop = FALSE], 2L, marker_mean, "-")
    Xtr[!is.finite(Xtr)] <- 0
    Xte[!is.finite(Xte)] <- 0
    y_mean <- mean(y[tr])
    ytr <- y[tr] - y_mean
    A <- tcrossprod(Xtr) + diag(lambda, sum(tr))
    alpha <- as.numeric(solve(A, ytr))
    beta <- as.numeric(crossprod(Xtr, alpha))
    pred[te] <- y_mean + as.numeric(Xte %*% beta)
  }
  pred
}

ng_predict_gebv <- function(geno, effects) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  beta <- effects$beta[colnames(geno)]
  beta[!is.finite(beta)] <- 0
  as.numeric(effects$intercept + geno %*% beta)
}

ng_choose_mean_source <- function(geno,
                                  effects,
                                  adjusted_pheno = NULL,
                                  blue = NULL,
                                  blup = NULL,
                                  ids = rownames(geno),
                                  min_reliability = 0.35,
                                  min_cv_predictive_r2 = 0.35) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  gebv <- setNames(ng_predict_gebv(geno[ids, , drop = FALSE], effects), ids)
  rel <- suppressWarnings(as.numeric(effects$reliability))
  if (length(rel) != 1L || !is.finite(rel)) rel <- NA_real_
  rel_flag <- effects$reliability_is_calibrated
  # A numeric field named `reliability` is not self-authenticating. Require the
  # producer to state explicitly that it is calibrated to prediction-error
  # variance / true breeding-value accuracy before it can gate mean selection.
  rel_calibrated <- isTRUE(rel_flag) && is.finite(rel)
  predictive_r2 <- suppressWarnings(as.numeric(effects$cv_predictive_r2))
  if (length(predictive_r2) != 1L || !is.finite(predictive_r2)) predictive_r2 <- NA_real_
  if (rel_calibrated && is.finite(rel) && rel >= min_reliability) {
    return(list(source = "GEBV", value = gebv, reliability = rel,
                reliability_is_calibrated = TRUE,
                cv_predictive_r2 = predictive_r2,
                mean_source_criterion = "calibrated_reliability"))
  }
  # No calibrated reliability. `ng_fit_ridge_effects` deliberately reports none,
  # because a phenotype cross-validation statistic is not accuracy^2 against true
  # breeding value -- so gating ONLY on reliability made this branch unreachable
  # and silently forced the phenotype mid-parent on every run, however well the
  # markers predicted. Fall back to the statistic that IS available and IS
  # honestly named: out-of-sample predictive R^2 against the phenotype.
  #
  # The two thresholds are NOT interchangeable. cv_predictive_r2 is bounded above
  # by heritability, so for the same genomic model it sits BELOW a true
  # reliability; sharing the 0.35 default therefore demands more of the markers
  # here, not less. Lower it deliberately (and record why) rather than assuming
  # the numbers mean the same thing.
  if (is.finite(predictive_r2) && is.finite(min_cv_predictive_r2) &&
      predictive_r2 >= min_cv_predictive_r2) {
    return(list(source = "GEBV", value = gebv,
                reliability = if (rel_calibrated) rel else NA_real_,
                reliability_is_calibrated = rel_calibrated,
                cv_predictive_r2 = predictive_r2,
                mean_source_criterion = "cv_predictive_r2"))
  }
  candidates <- list(BLUP = blup, BLUE = blue, adjusted_pheno = adjusted_pheno)
  for (nm in names(candidates)) {
    if (!is.null(candidates[[nm]])) {
      v <- ng_match_vector(candidates[[nm]], ids, nm)
      return(list(source = nm, value = setNames(v, ids),
                  reliability = if (rel_calibrated) rel else NA_real_,
                  reliability_is_calibrated = rel_calibrated,
                  cv_predictive_r2 = predictive_r2,
                  mean_source_criterion = "below_threshold"))
    }
  }
  list(source = if (rel_calibrated) "GEBV_low_reliability" else "GEBV_uncalibrated",
       value = gebv, reliability = if (rel_calibrated) rel else NA_real_,
       reliability_is_calibrated = rel_calibrated,
       cv_predictive_r2 = predictive_r2,
       mean_source_criterion = "no_phenotype_fallback")
}
