# Built-in estimator for the additive genetic covariance matrix G and the
# phenotypic covariance matrix P across traits. Both are needed by
# `ng_add_multitrait_score(method = "economic_index" | "desired_gain")` to drive
# Pesek-Baker and Smith-Hazel selection-index solves. The estimators are
# zero-external-dependency by default (two-stage ridge); they opportunistically
# use `sommer::mmer` when available for a fully multivariate REML fit.
#
# Statistical justification (two_stage_ridge):
#   For trait t, fit ridge marker effects beta_hat_t on the same n x m geno used
#   downstream. With centered genotypes X_c (mean 2*p subtracted) and
#   K = X_c X_c' / denom where denom = sum_k 2 p_k (1 - p_k) (VanRaden), the
#   additive-genetic covariance of breeding values g = X_c beta is
#       Var(g) = X_c X_c' * Var(beta) = K * denom * Var(beta_k)  (iid markers).
#
#   We construct G_hat from two parts that are individually robust:
#
#   (a) Diagonals via GBLUP lambda-inversion (method-of-moments). Under the
#       GBLUP model the optimal ridge penalty is
#           lambda_t = sigma_e^2_t / (sigma_g^2_t / denom),
#       so the genetic variance is identifiable from the chosen lambda and the
#       fitted residual variance as
#           sigma_g^2_t = sigma_e^2_t * denom / lambda_t.
#       This is a direct REML-style estimate that does not suffer from the
#       ridge-attenuation bias of inner-product plug-ins (the beta_hat scale
#       drops out entirely).
#
#   (b) Off-diagonals via the sample correlation of the trait marker-effect
#       vectors. Under independent-ridge fits with the SAME design X,
#           Cor(beta_hat_t, beta_hat_s) ~ Cor(beta_true_t, beta_true_s)
#       up to identical multiplicative shrinkage per trait, which cancels in
#       the Pearson correlation. The estimator is biased toward 0 by the
#       trait-specific noise in beta_hat, but it preserves the sign and the
#       relative ordering of the genetic correlations, which is what
#       downstream Smith-Hazel / Pesek-Baker solves use.
#
#   G_hat = D %*% R_beta %*% D with D = sqrt(diag(sigma_g^2)). The result is
#   on the same scale as a GBLUP variance component fit against the VanRaden
#   K, so it is directly comparable to ng_parent_kinship() downstream.
#
#   The G is then projected to the nearest positive-semidefinite matrix using
#   Matrix::nearPD() if available, or a documented eigen-clip fallback.
#
# PSD projection:
#   `Matrix::nearPD()` if installed; otherwise a simple eigen-clip that floors
#   negative eigenvalues at `eps_eig * mean(positive eigenvalues)`. nearPD is
#   preferred because it gives the closest-in-Frobenius PSD matrix; the
#   eigen-clip is a documented fallback so the function works without Matrix.

# --- helpers ----------------------------------------------------------------

ng_genetic_cov_eigen_clip <- function(M, eps_eig = 1e-6) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  vals <- ev$values
  pos <- vals[vals > 0]
  floor_val <- if (length(pos)) eps_eig * mean(pos) else eps_eig
  vals[vals < floor_val] <- floor_val
  ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
}

ng_genetic_cov_project_psd <- function(M, eps_eig = 1e-6) {
  M <- (M + t(M)) / 2
  if (requireNamespace("Matrix", quietly = TRUE)) {
    pd <- tryCatch(Matrix::nearPD(M, base.matrix = TRUE, keepDiag = FALSE,
                                  ensureSymmetry = TRUE),
                   error = function(e) NULL)
    if (!is.null(pd)) return(as.matrix(pd$mat))
  }
  ng_genetic_cov_eigen_clip(M, eps_eig = eps_eig)
}

ng_genetic_cov_to_correlation <- function(M) {
  d <- sqrt(pmax(diag(M), 0))
  d[!is.finite(d) | d <= 0] <- 1
  R <- M / outer(d, d)
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

# --- two-stage ridge engine -------------------------------------------------

ng_genetic_cov_two_stage_ridge <- function(geno, Y, ridge_lambda = NULL,
                                            kfold = 5L, seed = 1L,
                                            return_diagnostics = FALSE) {
  m <- ncol(geno)
  t <- ncol(Y)
  trait_names <- colnames(Y)

  p <- colMeans(geno, na.rm = TRUE) / 2
  denom <- sum(2 * p * (1 - p), na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) denom <- m

  Beta <- matrix(0, nrow = m, ncol = t)
  colnames(Beta) <- trait_names
  rownames(Beta) <- colnames(geno)
  reliabilities <- rep(NA_real_, t)
  lambdas <- rep(NA_real_, t)
  sigma_e2s <- rep(NA_real_, t)
  sigma_g2s <- rep(NA_real_, t)
  names(reliabilities) <- names(lambdas) <- names(sigma_e2s) <- names(sigma_g2s) <- trait_names

  for (j in seq_len(t)) {
    y_j <- Y[, j]
    fit_j <- ng_fit_ridge_effects(
      geno = geno,
      y = setNames(y_j, rownames(geno)),
      ids = rownames(geno),
      lambda = ridge_lambda,
      kfold = kfold,
      seed = seed + j
    )
    Beta[, j] <- fit_j$beta
    reliabilities[[j]] <- fit_j$reliability
    lambdas[[j]] <- fit_j$lambda
    sigma_e2s[[j]] <- fit_j$sigma_e2
    # GBLUP lambda-inversion: sigma_g^2 = sigma_e^2 * denom / lambda. Floor at
    # a tiny positive value so the projection step still produces a valid PSD G.
    sigma_g2_j <- fit_j$sigma_e2 * denom / max(fit_j$lambda, 1e-12)
    if (!is.finite(sigma_g2_j) || sigma_g2_j <= 0) sigma_g2_j <- 1e-8
    sigma_g2s[[j]] <- sigma_g2_j
  }

  # Off-diagonals: sample Pearson correlation of the m-vectors of beta_hats.
  # Pearson correlation is invariant to per-trait multiplicative shrinkage, so
  # this is robust to the ridge attenuation that wrecks naive inner-product
  # plug-ins under small n.
  if (m > 1L) {
    R_beta <- suppressWarnings(stats::cor(Beta))
    R_beta[!is.finite(R_beta)] <- 0
    diag(R_beta) <- 1
  } else {
    R_beta <- diag(t)
  }

  D <- sqrt(sigma_g2s)
  G_raw <- diag(D, t) %*% R_beta %*% diag(D, t)
  G_raw <- (G_raw + t(G_raw)) / 2
  dimnames(G_raw) <- list(trait_names, trait_names)

  G_hat <- ng_genetic_cov_project_psd(G_raw)
  dimnames(G_hat) <- list(trait_names, trait_names)

  out <- list(G_hat = G_hat,
              reliabilities = reliabilities,
              lambdas = lambdas,
              sigma_e2 = sigma_e2s,
              sigma_g2 = sigma_g2s,
              R_beta = R_beta,
              denom = denom,
              G_raw = G_raw,
              beta_hat = Beta)
  if (!return_diagnostics) {
    out$beta_hat <- NULL
    out$R_beta <- NULL
    out$G_raw <- NULL
  }
  out
}

# --- sommer REML engine -----------------------------------------------------

ng_genetic_cov_sommer_remml <- function(geno, Y) {
  if (!requireNamespace("sommer", quietly = TRUE)) {
    return(list(error = "sommer not installed"))
  }
  K <- ng_parent_kinship(geno)
  ids <- rownames(geno)
  trait_names <- colnames(Y)
  df <- as.data.frame(Y)
  df$id <- factor(ids, levels = ids)
  rownames(K) <- colnames(K) <- ids
  rhs <- paste(trait_names, collapse = ", ")
  form <- stats::as.formula(paste0("cbind(", rhs, ") ~ 1"))
  fit <- tryCatch({
    sommer::mmer(
      fixed = form,
      random = stats::as.formula("~ sommer::vsr(id, Gu = K, Gtc = sommer::unsm(length(trait_names)))"),
      rcov = stats::as.formula("~ sommer::vsr(units, Gtc = sommer::unsm(length(trait_names)))"),
      data = df,
      verbose = FALSE,
      tolParInv = 1e-3,
      date.warning = FALSE
    )
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(error = conditionMessage(fit)))
  }
  # sommer stores variance components in fit$sigma; the first element corresponds
  # to the random GRM term.
  G_hat <- tryCatch({
    sigma <- fit$sigma
    g_term <- sigma[[1L]]
    g_term <- as.matrix(g_term)
    g_term
  }, error = function(e) NULL)
  if (is.null(G_hat) || !all(dim(G_hat) == c(length(trait_names), length(trait_names)))) {
    return(list(error = "could not extract G from sommer fit"))
  }
  dimnames(G_hat) <- list(trait_names, trait_names)
  G_hat <- (G_hat + t(G_hat)) / 2
  list(G_hat = ng_genetic_cov_project_psd(G_hat))
}

# --- public estimator -------------------------------------------------------

# Estimate additive genetic covariance G across traits from a training set.
#
# `method = "auto"` chooses `sommer_remml` when the `sommer` package is
# installed, and otherwise falls back to `two_stage_ridge`. If the sommer fit
# errors (singular, non-convergence, etc.) we emit a `warning()` and fall back
# to `two_stage_ridge` so callers always receive a valid G_hat.
ng_estimate_genetic_covariance <- function(geno,
                                           Y,
                                           method = c("auto", "two_stage_ridge", "sommer_remml"),
                                           ridge_lambda = NULL,
                                           kfold = 5L,
                                           seed = 1L,
                                           return_diagnostics = FALSE) {
  method <- match.arg(method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  if (is.null(rownames(geno))) ng_stop("geno must have row names")
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  if (is.null(rownames(Y))) {
    if (nrow(Y) != nrow(geno)) ng_stop("Y must have row names or have nrow(Y) == nrow(geno)")
    rownames(Y) <- rownames(geno)
  }
  ids <- intersect(rownames(geno), rownames(Y))
  if (!length(ids)) ng_stop("no overlapping ids between geno and Y")
  geno <- ng_check_same_ids(geno, ids, "geno")
  Y <- Y[ids, , drop = FALSE]
  ok <- stats::complete.cases(Y)
  if (sum(ok) < 5L) ng_stop("need at least 5 complete-case rows to estimate G")
  geno <- geno[ok, , drop = FALSE]
  Y <- Y[ok, , drop = FALSE]
  storage.mode(Y) <- "double"

  resolved_method <- method
  if (method == "auto") {
    resolved_method <- if (requireNamespace("sommer", quietly = TRUE)) "sommer_remml" else "two_stage_ridge"
  }

  diag_list <- list()
  if (resolved_method == "sommer_remml") {
    sommer_out <- ng_genetic_cov_sommer_remml(geno, Y)
    if (!is.null(sommer_out$error)) {
      warning("sommer REML failed (", sommer_out$error, "); falling back to two_stage_ridge")
      resolved_method <- "two_stage_ridge"
    } else {
      G_hat <- sommer_out$G_hat
      diag_list$engine <- "sommer_remml"
    }
  }
  if (resolved_method == "two_stage_ridge") {
    ts <- ng_genetic_cov_two_stage_ridge(geno, Y, ridge_lambda = ridge_lambda,
                                          kfold = kfold, seed = seed,
                                          return_diagnostics = return_diagnostics)
    G_hat <- ts$G_hat
    diag_list$engine <- "two_stage_ridge"
    diag_list$reliabilities <- ts$reliabilities
    diag_list$lambdas <- ts$lambdas
    diag_list$sigma_e2 <- ts$sigma_e2
    diag_list$sigma_g2 <- ts$sigma_g2
    diag_list$denom <- ts$denom
    if (return_diagnostics) {
      diag_list$beta_hat <- ts$beta_hat
      diag_list$R_beta <- ts$R_beta
      diag_list$G_raw <- ts$G_raw
    }
  }

  dimnames(G_hat) <- list(colnames(Y), colnames(Y))
  attr(G_hat, "genetic_correlation") <- ng_genetic_cov_to_correlation(G_hat)
  attr(G_hat, "n_used") <- nrow(Y)
  attr(G_hat, "method") <- diag_list$engine
  attr(G_hat, "requested_method") <- method
  if (return_diagnostics) attr(G_hat, "diagnostics") <- diag_list
  G_hat
}

# --- phenotypic covariance --------------------------------------------------

# Ledoit-Wolf (2004) analytical optimal shrinkage of the sample covariance
# matrix toward a diagonal target. The closed-form intensity minimizes the
# expected Frobenius distance to the population covariance. Implemented from
# the original formulas (Ledoit & Wolf 2004, J. Multivariate Anal., eq. (14))
# so we do not need `corpcor`.
ng_genetic_cov_ledoit_wolf <- function(Y) {
  Y <- as.matrix(Y)
  Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  n <- nrow(Y)
  p <- ncol(Y)
  if (n < 2L) return(list(S_star = diag(p), intensity = 1, S = diag(p), T = diag(p)))
  mu <- colMeans(Y)
  Yc <- sweep(Y, 2L, mu, "-")
  # Sample covariance (1/n, not 1/(n-1), to match Ledoit-Wolf's derivation).
  S <- crossprod(Yc) / n
  # Diagonal target T (preserves the sample variances; assumes zero off-diag).
  T_target <- diag(diag(S), p)
  # pi-hat: sum_{i,j} Var(S_{ij}).
  Yc2 <- Yc * Yc
  pi_mat <- crossprod(Yc2) / n - S * S
  pi_hat <- sum(pi_mat)
  # rho-hat: contribution from the diagonal entries of T (off-diagonals of T
  # are zero, so only diagonal terms contribute).
  rho_hat <- sum(diag(pi_mat))
  # gamma-hat: ||S - T||_F^2.
  diff_ST <- S - T_target
  gamma_hat <- sum(diff_ST * diff_ST)
  # Optimal intensity in (0, 1).
  if (!is.finite(gamma_hat) || gamma_hat <= 0) {
    intensity <- 1
  } else {
    kappa <- (pi_hat - rho_hat) / gamma_hat
    intensity <- max(0, min(1, kappa / n))
  }
  S_star <- intensity * T_target + (1 - intensity) * S
  list(S_star = S_star, intensity = intensity, S = S, T = T_target)
}

ng_estimate_phenotypic_covariance <- function(Y, shrinkage = c("none", "auto")) {
  shrinkage <- match.arg(shrinkage)
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  trait_names <- colnames(Y)
  if (shrinkage == "none") {
    P <- stats::cov(Y, use = "pairwise.complete.obs")
    P[!is.finite(P)] <- 0
    P <- (P + t(P)) / 2
    P <- ng_genetic_cov_project_psd(P)
    dimnames(P) <- list(trait_names, trait_names)
    attr(P, "shrinkage") <- "none"
    attr(P, "intensity") <- 0
    return(P)
  }
  # Ledoit-Wolf auto shrinkage toward a diagonal target.
  lw <- ng_genetic_cov_ledoit_wolf(Y)
  P <- lw$S_star
  P <- (P + t(P)) / 2
  P <- ng_genetic_cov_project_psd(P)
  dimnames(P) <- list(trait_names, trait_names)
  attr(P, "shrinkage") <- "ledoit_wolf"
  attr(P, "intensity") <- lw$intensity
  attr(P, "sample_covariance") <- lw$S
  P
}

# --- posterior G across traits ---------------------------------------------

# Posterior over the additive genetic covariance matrix G across traits.
#
# Two modes:
#
#   method = "beta_posterior" (default; fast, coherent with D1).
#     For each trait, fit per-trait ridge once and draw S samples of beta_t
#     from the BCM closed-form posterior conditional on (sigma_e2_t, lambda_t).
#     For draw s, construct G_s = D R_beta^s D where D = sqrt(sigma_g2) is
#     fixed across draws (lambda-inversion is hyperparameter-conditional) and
#     R_beta^s = cor(beta_1^s, ..., beta_t^s) recomputed per draw. This is
#     the natural Bayesian companion to ng_posterior_cross_predict() since it
#     consumes the same BCM samples; the diagonals are fixed and only the
#     correlation matrix carries posterior uncertainty.
#
#   method = "parametric_bootstrap" (slower; covers hyperparameter uncertainty).
#     Resample Gaussian residuals around each per-trait ridge fit, regenerate
#     y_t^b = fitted_t + epsilon_t^b, and rerun ng_genetic_cov_two_stage_ridge
#     end-to-end (CV-tuning lambda included). This integrates over the joint
#     uncertainty in (beta, sigma_e2, lambda) at the cost of refitting B times.
#
# Returns a t x t x n_draws array with attribute `G_mean` (posterior mean
# matrix), `method`, `n_draws`. Off-diagonal credible intervals can be read
# off via apply(G_draws, c(1, 2), quantile, probs = c(.025, .975)).
ng_posterior_genetic_covariance <- function(geno,
                                            Y,
                                            n_draws = 100L,
                                            method = c("beta_posterior", "parametric_bootstrap"),
                                            ridge_lambda = NULL,
                                            kfold = 5L,
                                            seed = 1L) {
  method <- match.arg(method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  if (nrow(Y) != nrow(geno)) ng_stop("Y and geno must have the same number of rows")
  ids <- rownames(geno)
  trait_names <- colnames(Y)
  n_traits <- ncol(Y)
  m <- ncol(geno)
  p <- colMeans(geno, na.rm = TRUE) / 2
  denom <- sum(2 * p * (1 - p), na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) denom <- m

  G_draws <- array(0, dim = c(n_traits, n_traits, n_draws),
                   dimnames = list(trait_names, trait_names, NULL))

  if (identical(method, "beta_posterior")) {
    sigma_g2 <- numeric(n_traits); names(sigma_g2) <- trait_names
    beta_draw_list <- vector("list", n_traits)
    for (j in seq_len(n_traits)) {
      y_j <- as.numeric(Y[, j])
      names(y_j) <- ids
      post_j <- ng_fit_ridge_effects_posterior(
        geno = geno, y = y_j, ids = ids,
        lambda = ridge_lambda, kfold = kfold,
        n_draws = n_draws, method = "closed_form",
        seed = seed + j
      )
      sigma_g2[[j]] <- post_j$fit$sigma_e2 * denom / max(post_j$fit$lambda, 1e-12)
      if (!is.finite(sigma_g2[[j]]) || sigma_g2[[j]] <= 0) sigma_g2[[j]] <- 1e-8
      beta_draw_list[[j]] <- post_j$beta_draws
    }
    D <- sqrt(sigma_g2)
    for (b in seq_len(n_draws)) {
      B_b <- do.call(cbind, lapply(beta_draw_list, function(M) M[, b]))
      if (m > 1L) {
        R_beta_b <- suppressWarnings(stats::cor(B_b))
        R_beta_b[!is.finite(R_beta_b)] <- 0
        diag(R_beta_b) <- 1
      } else {
        R_beta_b <- diag(n_traits)
      }
      G_b <- diag(D, n_traits) %*% R_beta_b %*% diag(D, n_traits)
      G_b <- (G_b + t(G_b)) / 2
      G_b <- ng_genetic_cov_project_psd(G_b)
      dimnames(G_b) <- list(trait_names, trait_names)
      G_draws[, , b] <- G_b
    }
  } else {
    # parametric_bootstrap: resample residuals around each per-trait ridge fit.
    fits <- vector("list", n_traits)
    fitted_list <- vector("list", n_traits)
    ok_list <- vector("list", n_traits)
    for (j in seq_len(n_traits)) {
      y_j <- as.numeric(Y[, j])
      names(y_j) <- ids
      fit_j <- ng_fit_ridge_effects(
        geno = geno, y = y_j, ids = ids,
        lambda = ridge_lambda, kfold = kfold, seed = seed + j
      )
      fits[[j]] <- fit_j
      ok_j <- is.finite(y_j)
      ok_list[[j]] <- ok_j
      fitted_list[[j]] <- as.numeric(ng_predict_gebv(geno[ok_j, , drop = FALSE], fit_j))
    }
    set.seed(seed)
    for (b in seq_len(n_draws)) {
      Y_b <- Y
      for (j in seq_len(n_traits)) {
        eps_j <- stats::rnorm(sum(ok_list[[j]]), sd = sqrt(fits[[j]]$sigma_e2))
        Y_b[ok_list[[j]], j] <- fitted_list[[j]] + eps_j
      }
      out_b <- ng_genetic_cov_two_stage_ridge(
        geno = geno, Y = Y_b,
        ridge_lambda = ridge_lambda, kfold = kfold, seed = seed + b
      )
      G_draws[, , b] <- out_b$G_hat
    }
  }

  G_mean <- apply(G_draws, c(1L, 2L), mean)
  dimnames(G_mean) <- list(trait_names, trait_names)
  attr(G_draws, "G_mean") <- G_mean
  attr(G_draws, "method") <- method
  attr(G_draws, "n_draws") <- n_draws
  attr(G_draws, "n_used") <- nrow(Y)
  G_draws
}
