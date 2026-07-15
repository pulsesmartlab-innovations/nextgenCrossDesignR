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
  X <- geno[ok, , drop = FALSE]
  y <- y[ok]
  marker_mean <- colMeans(X, na.rm = TRUE)
  X <- sweep(X, 2, marker_mean, "-")
  X[!is.finite(X)] <- 0
  yc <- y - mean(y)

  if (is.null(lambda)) {
    if (!is.null(h2_prior)) {
      h2_prior <- max(min(as.numeric(h2_prior), 0.99), 0.01)
      lambda0 <- ncol(X) * (1 - h2_prior) / h2_prior
      lambda_grid <- unique(sort(c(lambda0, lambda0 * 10 ^ seq(-2, 2, length.out = 17))))
    }
    lambda <- ng_choose_ridge_lambda(X, yc, lambda_grid, kfold = kfold, seed = seed)
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
  sigma_e2 <- sum(resid^2) / max(1, length(y) - df_eff)
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
  in_sample_reliability <- suppressWarnings(stats::cor(y, fitted, use = "complete.obs")^2)
  if (!is.finite(in_sample_reliability)) in_sample_reliability <- 0
  cv <- ng_ridge_cv_predict(X, yc, lambda = lambda, kfold = kfold, seed = seed)
  reliability <- if (!is.null(cv)) {
    rel <- suppressWarnings(stats::cor(y, mean(y) + cv, use = "complete.obs")^2)
    if (is.finite(rel)) rel else 0
  } else {
    in_sample_reliability
  }

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
    reliability = reliability,
    in_sample_reliability = in_sample_reliability
  )
}

ng_choose_ridge_lambda <- function(X, y, lambda_grid, kfold = 5, seed = 1L) {
  n <- nrow(X)
  if (n < 10 || kfold < 2) return(stats::median(lambda_grid))
  set.seed(seed)
  folds <- sample(rep(seq_len(kfold), length.out = n))
  pred <- matrix(NA_real_, nrow = n, ncol = length(lambda_grid))
  for (fold in seq_len(kfold)) {
    tr <- folds != fold
    te <- !tr
    Xtr <- X[tr, , drop = FALSE]
    Xte <- X[te, , drop = FALSE]
    eig <- eigen(tcrossprod(Xtr), symmetric = TRUE)
    vals <- pmax(eig$values, 0)
    Uy <- as.numeric(crossprod(eig$vectors, y[tr]))
    KteU <- tcrossprod(Xte, Xtr) %*% eig$vectors
    for (g in seq_along(lambda_grid)) {
      lam <- lambda_grid[g]
      pred[te, g] <- as.numeric(KteU %*% (Uy / (vals + lam)))
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
    A <- tcrossprod(X[tr, , drop = FALSE]) + diag(lambda, sum(tr))
    alpha <- as.numeric(solve(A, y[tr]))
    beta <- as.numeric(crossprod(X[tr, , drop = FALSE], alpha))
    pred[te] <- as.numeric(X[te, , drop = FALSE] %*% beta)
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
                                  min_reliability = 0.35) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  gebv <- setNames(ng_predict_gebv(geno[ids, , drop = FALSE], effects), ids)
  rel <- effects$reliability
  if (is.finite(rel) && rel >= min_reliability) {
    return(list(source = "GEBV", value = gebv, reliability = rel))
  }
  candidates <- list(BLUP = blup, BLUE = blue, adjusted_pheno = adjusted_pheno)
  for (nm in names(candidates)) {
    if (!is.null(candidates[[nm]])) {
      v <- ng_match_vector(candidates[[nm]], ids, nm)
      return(list(source = nm, value = setNames(v, ids), reliability = rel))
    }
  }
  list(source = "GEBV_low_reliability", value = gebv, reliability = rel)
}
