# Posterior-aware cross prediction (v0.2.0)
#
# Most genomic-prediction tools for cross design (PopVar, SimpleMating,
# AlphaMate, genomicMateSelectR) treat marker-effect estimates as fixed and
# report point estimates of cross mean, variance, and usefulness. With the
# corrected ridge posterior variance (T1.3 in v0.1.0), `nextgenCrossDesign`
# can ship a posterior-aware layer that propagates marker-effect uncertainty
# through to per-cross credible intervals, the probability of producing a
# progeny that clears a hard threshold tau, and posterior rank stability for
# OCS. This file implements:
#
#   - ng_sample_ridge_posterior_bcm()   closed-form posterior beta draws via
#                                       Bhattacharya-Chakraborty-Mallick (2016).
#   - ng_sample_ridge_posterior_mcmc()  Gibbs sampler over (beta, sigma_e2,
#                                       sigma_beta2) for full hyperparameter
#                                       posterior coverage. Slower; opt-in.
#   - ng_fit_ridge_effects_posterior()  user-facing wrapper that fits ridge
#                                       and returns S draws of beta plus
#                                       (sigma_e2, lambda) draws.
#   - ng_p_superior_progeny()           closed-form P(max of k progeny >= tau)
#                                       under the Gaussian family assumption.
#                                       Standalone helper; D4.
#   - ng_posterior_cross_predict()      per-cross usefulness credible
#                                       interval, P(superior progeny) CI,
#                                       and posterior_topn_prob for OCS.
#   - ng_optimize_robust_mating_plan()  posterior-aware OCS that maximizes a
#                                       robustness quantile of posterior
#                                       usefulness rather than the point
#                                       estimate. Reuses ng_optimize_mating_plan.

# ---- Posterior beta draws --------------------------------------------------

# Closed-form sampler from N(beta_hat, sigma_e2 * (X'X + lambda I)^{-1}) via
# Bhattacharya-Chakraborty-Mallick (2016). Conditional on (sigma_e2, lambda)
# being known. The dual factorization A = XX' + lambda I is computed once,
# then each draw costs O(n^2) back-substitution + O(np) for the X'w multiply.
#
# Algorithm (BCM, with prior beta ~ N(0, sigma_e2 / lambda * I_m)):
#   1. theta ~ N(0, (sigma_e2 / lambda) I_m)
#   2. eta   ~ N(0, sigma_e2 I_n)
#   3. nu    = X theta + eta
#   4. w     = (XX' + lambda I)^{-1} (yc - nu)
#   5. beta_draw = theta + X' w
#
# Centering: this returns draws of (beta - 0) given yc = y - mean(y); the
# overall intercept shift is handled by the caller, identical to the way
# ng_fit_ridge_effects() does.
ng_sample_ridge_posterior_bcm <- function(X, yc, sigma_e2, lambda, n_draws,
                                          seed = 1L, use_cpp = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  if (!is.finite(sigma_e2) || sigma_e2 <= 0) ng_stop("sigma_e2 must be > 0")
  if (!is.finite(lambda) || lambda <= 0) ng_stop("lambda must be > 0 for the BCM sampler")
  A <- tcrossprod(X) + diag(lambda, n)
  # Cholesky factorize once and reuse for all draws.
  A_chol <- tryCatch(chol(A), error = function(e) NULL)
  # Posterior reporting must not perturb the caller's RNG stream: doing so can
  # indirectly change a later unseeded stochastic allocation even though none
  # of these posterior columns belongs to the default allocation objective.
  ng_with_rng_seed(seed, {
    if (isTRUE(use_cpp) && !is.null(A_chol) &&
        exists("ng_bcm_posterior_sampler_cpp", mode = "function", inherits = TRUE)) {
      # C++ kernel: per-draw rnorm via R::rnorm (no callback overhead) +
      # hand-rolled forward/back triangular solve on the upper Cholesky factor +
      # cache-friendly X' w multiplication. Reproduces stats::rnorm under the
      # same seed for numerical equivalence with the R reference.
      ng_bcm_posterior_sampler_cpp(
        X_centered = X, yc = as.numeric(yc),
        A_chol_upper = A_chol,
        sigma_e2 = sigma_e2, lambda = lambda,
        n_draws = as.integer(n_draws), seed = as.integer(seed)
      )
    } else {
      solve_A <- if (is.null(A_chol)) {
        function(b) as.numeric(solve(A, b))
      } else {
        function(b) as.numeric(backsolve(A_chol, backsolve(A_chol, b, transpose = TRUE)))
      }
      beta_draws <- matrix(0, nrow = p, ncol = n_draws)
      sd_prior <- sqrt(sigma_e2 / lambda)
      sd_noise <- sqrt(sigma_e2)
      for (s in seq_len(n_draws)) {
        theta <- stats::rnorm(p, sd = sd_prior)
        eta   <- stats::rnorm(n, sd = sd_noise)
        nu    <- as.numeric(X %*% theta + eta)
        w     <- solve_A(yc - nu)
        beta_draws[, s] <- theta + as.numeric(crossprod(X, w))
      }
      rownames(beta_draws) <- colnames(X)
      beta_draws
    }
  })
}

# Gibbs sampler over (beta, sigma_e2, sigma_beta2). Vague inverse-gamma
# hyperpriors on the variance components: sigma_e2 ~ IG(a_e, b_e),
# sigma_beta2 ~ IG(a_b, b_b). lambda is treated as a derived quantity
# lambda = sigma_e2 / sigma_beta2. Each iteration:
#   1. beta | y, sigma_e2, lambda via BCM (uses current lambda)
#   2. sigma_e2 | y, beta ~ IG(a_e + n/2, b_e + 0.5 ||y - Xbeta||^2)
#   3. sigma_beta2 | beta ~ IG(a_b + m/2, b_b + 0.5 ||beta||^2)
#   4. lambda = sigma_e2 / sigma_beta2
#
# Standard reference: Park & Casella (2008) for Bayesian Lasso; the ridge
# variant is straightforward and is what BGLR / sommer ultimately implement
# under the hood. We use vague hyperpriors so the posterior is data-driven.
ng_sample_ridge_posterior_mcmc <- function(X, yc, n_draws,
                                           burnin = 500L,
                                           thin = 1L,
                                           a_e = 1e-3, b_e = 1e-3,
                                           a_b = 1e-3, b_b = 1e-3,
                                           lambda_init = NULL,
                                           sigma_e2_init = NULL,
                                           seed = 1L) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  yc <- as.numeric(yc)
  ng_with_rng_seed(seed, {
  total_iter <- as.integer(burnin) + as.integer(n_draws) * as.integer(max(1L, thin))
  if (is.null(sigma_e2_init) || !is.finite(sigma_e2_init) || sigma_e2_init <= 0) {
    sigma_e2_init <- max(stats::var(yc, na.rm = TRUE), 1e-4)
  }
  if (is.null(lambda_init) || !is.finite(lambda_init) || lambda_init <= 0) {
    lambda_init <- p
  }
  sigma_e2 <- sigma_e2_init
  sigma_beta2 <- sigma_e2_init / lambda_init
  beta_draws <- matrix(0, nrow = p, ncol = n_draws)
  sigma_e2_draws <- numeric(n_draws)
  sigma_beta2_draws <- numeric(n_draws)
  store_idx <- 0L
  XtX_diag_cache <- NULL  # not needed; each iter rebuilds A
  for (it in seq_len(total_iter)) {
    lambda <- sigma_e2 / sigma_beta2
    # Step 1: beta draw via BCM (one draw at the current hyperparameter state).
    A <- tcrossprod(X) + diag(lambda, n)
    A_chol <- tryCatch(chol(A), error = function(e) NULL)
    solve_A <- if (is.null(A_chol)) {
      function(b) as.numeric(solve(A, b))
    } else {
      function(b) as.numeric(backsolve(A_chol, backsolve(A_chol, b, transpose = TRUE)))
    }
    sd_prior <- sqrt(sigma_e2 / lambda)
    sd_noise <- sqrt(sigma_e2)
    theta <- stats::rnorm(p, sd = sd_prior)
    eta   <- stats::rnorm(n, sd = sd_noise)
    nu    <- as.numeric(X %*% theta + eta)
    w     <- solve_A(yc - nu)
    beta <- theta + as.numeric(crossprod(X, w))

    # Step 2: sigma_e2 from inverse-gamma full conditional.
    resid <- as.numeric(yc - X %*% beta)
    rss <- sum(resid * resid)
    shape_e <- a_e + n / 2
    rate_e  <- b_e + 0.5 * rss
    sigma_e2 <- 1 / stats::rgamma(1L, shape = shape_e, rate = rate_e)
    if (!is.finite(sigma_e2) || sigma_e2 <= 0) sigma_e2 <- max(rss / max(1, n - 1), 1e-6)

    # Step 3: sigma_beta2 from inverse-gamma full conditional.
    bss <- sum(beta * beta)
    shape_b <- a_b + p / 2
    rate_b  <- b_b + 0.5 * bss
    sigma_beta2 <- 1 / stats::rgamma(1L, shape = shape_b, rate = rate_b)
    if (!is.finite(sigma_beta2) || sigma_beta2 <= 0) sigma_beta2 <- max(bss / max(1, p - 1), 1e-8)

    if (it > burnin && ((it - burnin) %% max(1L, thin) == 0L)) {
      store_idx <- store_idx + 1L
      if (store_idx > n_draws) break
      beta_draws[, store_idx] <- beta
      sigma_e2_draws[store_idx] <- sigma_e2
      sigma_beta2_draws[store_idx] <- sigma_beta2
    }
  }
  rownames(beta_draws) <- colnames(X)
  list(
    beta_draws = beta_draws,
    sigma_e2_draws = sigma_e2_draws,
    sigma_beta2_draws = sigma_beta2_draws,
    lambda_draws = sigma_e2_draws / sigma_beta2_draws,
    iter_total = total_iter,
    burnin = burnin,
    thin = thin
  )
  })
}

# User-facing wrapper. method = "closed_form" uses BCM at the cross-validated
# lambda (empirical Bayes; same defensibility as rrBLUP under EB). method =
# "mcmc" runs the Gibbs sampler with hyperparameter posterior coverage.
ng_fit_ridge_effects_posterior <- function(geno, y,
                                           ids = rownames(geno),
                                           lambda = NULL,
                                           h2_prior = NULL,
                                           kfold = 5L,
                                           n_draws = 500L,
                                           method = c("closed_form", "mcmc"),
                                           mcmc_burnin = 500L,
                                           mcmc_thin = 1L,
                                           seed = 1L) {
  method <- match.arg(method)
  positive_integer <- function(x, name, minimum = 1L) {
    z <- suppressWarnings(as.numeric(x))
    if (length(z) != 1L || !is.finite(z) || z < minimum ||
        abs(z - round(z)) > 1e-8) {
      ng_stop(name, " must be one integer >= ", minimum)
    }
    as.integer(round(z))
  }
  n_draws <- positive_integer(n_draws, "n_draws")
  mcmc_burnin <- positive_integer(mcmc_burnin, "mcmc_burnin", minimum = 0L)
  mcmc_thin <- positive_integer(mcmc_thin, "mcmc_thin")
  fit <- ng_fit_ridge_effects(geno, y, ids = ids, lambda = lambda,
                              h2_prior = h2_prior, kfold = kfold, seed = seed)
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  y_vec <- ng_match_vector(y, ids, "y")
  ok <- is.finite(y_vec)
  X <- geno[ok, , drop = FALSE]
  marker_mean <- colMeans(X, na.rm = TRUE)
  marker_mean[!is.finite(marker_mean)] <- 0
  X <- sweep(X, 2L, marker_mean, "-")
  X[!is.finite(X)] <- 0
  yc <- y_vec[ok] - mean(y_vec[ok])

  if (identical(method, "closed_form")) {
    beta_draws <- ng_sample_ridge_posterior_bcm(
      X = X, yc = yc, sigma_e2 = fit$sigma_e2, lambda = fit$lambda,
      n_draws = n_draws, seed = seed
    )
    sigma_e2_draws <- rep(fit$sigma_e2, n_draws)
    lambda_draws   <- rep(fit$lambda, n_draws)
    sigma_beta2_draws <- sigma_e2_draws / lambda_draws
  } else {
    mc <- ng_sample_ridge_posterior_mcmc(
      X = X, yc = yc, n_draws = n_draws,
      burnin = mcmc_burnin, thin = mcmc_thin,
      lambda_init = fit$lambda, sigma_e2_init = fit$sigma_e2, seed = seed
    )
    beta_draws        <- mc$beta_draws
    sigma_e2_draws    <- mc$sigma_e2_draws
    sigma_beta2_draws <- mc$sigma_beta2_draws
    lambda_draws      <- mc$lambda_draws
  }
  list(
    fit = fit,
    beta_draws = beta_draws,
    sigma_e2_draws = sigma_e2_draws,
    sigma_beta2_draws = sigma_beta2_draws,
    lambda_draws = lambda_draws,
    marker_mean = marker_mean,
    method = method,
    n_draws = n_draws,
    seed = seed
  )
}

# ---- P(superior progeny >= tau) helpers (D4) -------------------------------

# Under a Gaussian within-family model X_i ~ N(mu, sigma^2) (i = 1..k),
# the probability that the family produces at least one line >= tau is
#   1 - Phi((tau - mu) / sigma) ^ k
# This is closed form and stable provided sigma > 0. For sigma = 0 the
# answer is the indicator mu >= tau.
ng_p_superior_progeny <- function(mu, sigma, tau, k_progeny) {
  mu <- as.numeric(mu); sigma <- as.numeric(sigma)
  tau <- as.numeric(tau); k_progeny <- as.numeric(k_progeny)
  out_len <- max(length(mu), length(sigma), length(tau), length(k_progeny))
  input_lengths <- c(length(mu), length(sigma), length(tau), length(k_progeny))
  if (out_len < 1L || any(!(input_lengths %in% c(1L, out_len)))) {
    ng_stop("mu, sigma, tau, and k_progeny must have length 1 or a common output length")
  }
  mu <- rep(mu, length.out = out_len)
  sigma <- rep(sigma, length.out = out_len)
  tau <- rep(tau, length.out = out_len)
  k_progeny <- rep(k_progeny, length.out = out_len)
  if (any(!is.finite(mu))) ng_stop("mu must contain finite family means")
  if (any(!is.finite(sigma)) || any(sigma < 0)) {
    ng_stop("sigma must contain finite non-negative family SDs")
  }
  if (any(is.na(tau))) ng_stop("tau must not be missing")
  if (any(!is.finite(k_progeny)) || any(k_progeny < 1) ||
      any(abs(k_progeny - round(k_progeny)) > 1e-8)) {
    ng_stop("k_progeny must contain positive integers")
  }
  k_progeny <- round(k_progeny)
  z <- (tau - mu) / pmax(sigma, .Machine$double.eps)
  # log-space to avoid (~1)^large precision loss when mu >> tau.
  log_pnorm_below <- stats::pnorm(z, lower.tail = TRUE, log.p = TRUE)
  one_minus_p <- exp(k_progeny * log_pnorm_below)
  out <- 1 - one_minus_p
  out[sigma <= 0] <- as.numeric(mu[sigma <= 0] >= tau[sigma <= 0])
  pmax(pmin(out, 1), 0)
}

# Add P(superior progeny) columns to an existing score table without changing
# anything else. Useful when the user wants the threshold-clearing probability
# but is not ready for the full posterior pipeline.
ng_add_p_superior_progeny <- function(scores,
                                      tau_superior,
                                      k_progeny = 100L,
                                      mean_col = "cross_mean_blend",
                                      var_col = "pmv",
                                      out_col = "p_superior_progeny") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!(mean_col %in% names(scores))) ng_stop("scores missing mean_col: ", mean_col)
  if (!(var_col %in% names(scores))) ng_stop("scores missing var_col: ", var_col)
  mu <- as.numeric(scores[[mean_col]])
  v <- suppressWarnings(as.numeric(scores[[var_col]]))
  if (any(!is.finite(mu))) ng_stop("mean_col must contain finite family means")
  if (any(!is.finite(v)) || any(v < 0)) {
    ng_stop("var_col must contain finite non-negative variances")
  }
  sigma <- sqrt(v)
  scores[[out_col]] <- ng_p_superior_progeny(mu, sigma, tau_superior, k_progeny)
  attr(scores, "p_superior_progeny") <- list(
    tau = tau_superior, k_progeny = k_progeny,
    mean_col = mean_col, var_col = var_col, out_col = out_col
  )
  scores
}

# ---- Posterior cross prediction --------------------------------------------

# Given posterior draws of beta (and matching sigma_e2/lambda draws), compute
# per-cross posterior summaries by running the cross-scoring kernel once per
# draw. Returns the point-estimate score table from ng_score_crosses() plus
# the following extra columns (all guarded with NA when not computable):
#
#   - pmv_post_mean, _lower, _upper   posterior CI on PMV
#   - usefulness_pmv_gebv_post_mean, _lower, _upper   posterior CI on usefulness
#   - p_superior_progeny_post_mean, _lower, _upper
#                                            posterior CI on P(max progeny >= tau)
#   - posterior_topn_prob_<N>                fraction of draws in which the
#                                            cross is in the top-N by usefulness
#
# Uses the existing ng_score_crosses() pipeline per draw, swapping in the
# posterior beta. The C++ kernel makes this practical for Haldane DH at
# realistic marker counts (~minutes for 80 parents, 5K markers, 500 draws).
# For Kosambi / RIL the dense path is much slower; the function emits a
# warning and recommends reducing n_draws.
ng_posterior_cross_predict <- function(geno,
                                       posterior_effects,
                                       marker_map = NULL,
                                       ids = rownames(geno),
                                       pairs = NULL,
                                       adjusted_pheno = NULL,
                                       blue = NULL,
                                       blup = NULL,
                                       include_self = FALSE,
                                       target = c("DH", "RIL"),
                                       parent_type = c("inbred", "dh", "ril"),
                                       selection_prop = 0.10,
                                       min_effect_reliability = 0.35,
                                       recomb_model = c("haldane", "kosambi"),
                                       window_cm = Inf,
                                       use_cpp = TRUE,
                                       assume_inbred = NULL,
                                       ci_level = 0.95,
                                       gain_col = "usefulness_pmv_gebv",
                                       var_col = "pmv",
                                       value_fun = NULL,
                                       tau_superior = NULL,
                                       k_progeny = 100L,
                                       top_n_targets = c(10L, 20L, 50L)) {
  target <- match.arg(target)
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  recomb_model <- match.arg(recomb_model)
  if (!is.list(posterior_effects) || is.null(posterior_effects$beta_draws)) {
    ng_stop("posterior_effects must be the list returned by ng_fit_ridge_effects_posterior()")
  }
  if (!identical(recomb_model, "haldane") || !identical(target, "DH")) {
    # The banded kernel can take this off the O(M^2) hot path when the user
    # supplies a finite window_cm small enough relative to mean chromosome
    # length (see ng_banded_kernel_preferred()). The warning is only
    # informative when we know we will still build the full m x m R per draw —
    # i.e., window_cm = Inf or the heuristic still picks the dense path.
    map_for_heuristic <- ng_prepare_marker_map(
      marker_map, marker_ids = colnames(geno), model = recomb_model
    )
    will_use_dense <- !is.finite(window_cm) ||
      !ng_banded_kernel_preferred(map_for_heuristic, window_cm)
    if (will_use_dense) {
      warning(sprintf(
        "Posterior cross prediction uses the dense O(M^2) recombination path under recomb_model = '%s', target = '%s' (window_cm = %s). Per-draw cost is ~O(pairs * M^2); consider reducing n_draws, LD-pruning markers, or passing a finite window_cm small relative to chromosome length.",
        recomb_model, target,
        if (is.finite(window_cm)) format(window_cm) else "Inf"
      ), call. = FALSE)
    }
  }
  beta_draws <- posterior_effects$beta_draws
  sigma_e2_draws <- posterior_effects$sigma_e2_draws
  lambda_draws <- posterior_effects$lambda_draws
  ci_level <- suppressWarnings(as.numeric(ci_level))
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    ng_stop("ci_level must be one finite probability in (0, 1)")
  }
  S <- ncol(beta_draws)
  if (S < 1L) ng_stop("posterior_effects$beta_draws has no draws")
  if (length(sigma_e2_draws) != S || length(lambda_draws) != S ||
      any(!is.finite(sigma_e2_draws)) || any(sigma_e2_draws <= 0) ||
      any(!is.finite(lambda_draws)) || any(lambda_draws <= 0)) {
    ng_stop("posterior sigma_e2_draws and lambda_draws must be positive, finite, and match beta_draws")
  }
  fit <- posterior_effects$fit

  # Point-estimate base table — gives us the cross identifiers plus the
  # mean-source / effect_reliability columns, all of which depend on the fit
  # not on the draw (they are computed from y, not beta).
  base <- ng_score_crosses(
    geno = geno, effects = fit, marker_map = marker_map, ids = ids,
    pairs = pairs, adjusted_pheno = adjusted_pheno, blue = blue, blup = blup,
    include_self = include_self, target = target,
    selection_prop = selection_prop,
    min_effect_reliability = min_effect_reliability,
    recomb_model = recomb_model, window_cm = window_cm, use_cpp = use_cpp,
    parent_type = parent_type
  )
  n_pairs <- nrow(base)
  if (!(gain_col %in% names(base))) ng_stop("base scores missing gain_col: ", gain_col)
  if (!(var_col %in% names(base))) ng_stop("base scores missing var_col: ", var_col)

  pmv_mat <- matrix(NA_real_, nrow = n_pairs, ncol = S)
  uc_mat  <- matrix(NA_real_, nrow = n_pairs, ncol = S)
  mu_mat <- if (!is.null(tau_superior)) matrix(NA_real_, nrow = n_pairs, ncol = S) else NULL
  if (!is.null(tau_superior) && !(var_col %in% c("vpm", "pmv"))) {
    ng_stop("tau_superior requires var_col = 'vpm' or 'pmv', a within-family genetic variance")
  }
  i_intensity <- ng_selection_intensity(selection_prop)
  # cross_mean_gebv changes with beta_draws because GEBV depends on beta.
  pair_p1 <- match(base$parent1, ids)
  pair_p2 <- match(base$parent2, ids)
  geno_mat <- ng_as_numeric_matrix(geno, "geno")
  intercept_y <- if (is.null(fit$intercept)) 0 else fit$intercept
  marker_mean_y <- if (is.null(fit$marker_mean)) colMeans(geno_mat) else fit$marker_mean

  for (s in seq_len(S)) {
    effects_s <- list(
      beta = beta_draws[, s],
      beta_var = setNames(rep(0, nrow(beta_draws)), rownames(beta_draws)),
      sigma_e2 = sigma_e2_draws[s],
      lambda = lambda_draws[s],
      intercept = mean(stats::na.omit(intercept_y)) -
                  sum(marker_mean_y * beta_draws[, s] - marker_mean_y * fit$beta),
      marker_mean = marker_mean_y,
      reliability = fit$reliability,
      cv_predictive_r2 = fit$cv_predictive_r2,
      reliability_is_calibrated = FALSE
    )
    # Centering shift: effects_s$intercept above keeps fitted ~ original.
    scored <- ng_score_crosses(
      geno = geno, effects = effects_s, marker_map = marker_map, ids = ids,
      pairs = base[, c("parent1", "parent2")],
      adjusted_pheno = adjusted_pheno, blue = blue, blup = blup,
      target = target, selection_prop = selection_prop,
      min_effect_reliability = min_effect_reliability,
      recomb_model = recomb_model, window_cm = window_cm, use_cpp = use_cpp,
      parent_type = parent_type
    )
    pmv_mat[, s] <- scored[[var_col]]
    if (!is.null(mu_mat)) mu_mat[, s] <- scored$cross_mean_gebv
    # `value_fun` lets a caller summarize the metric it actually RANKS on rather than the
    # hardcoded usefulness column. Without it, a posterior interval computed on usefulness would
    # be presented as the uncertainty of a `mean` / `var_complex` run -- the wrong quantity.
    uc_mat[, s]  <- if (is.null(value_fun)) scored[[gain_col]] else as.numeric(value_fun(scored))
  }

  alpha <- (1 - ci_level) / 2
  q_lower <- alpha
  q_upper <- 1 - alpha
  row_quantile <- function(M, q) {
    apply(M, 1L, function(row) {
      finite_row <- row[is.finite(row)]
      if (!length(finite_row)) return(NA_real_)
      unname(stats::quantile(finite_row, probs = q, names = FALSE))
    })
  }
  base$pmv_post_mean  <- rowMeans(pmv_mat, na.rm = TRUE)
  base$pmv_post_lower <- row_quantile(pmv_mat, q_lower)
  base$pmv_post_upper <- row_quantile(pmv_mat, q_upper)
  base[[paste0(gain_col, "_post_mean")]]  <- rowMeans(uc_mat, na.rm = TRUE)
  base[[paste0(gain_col, "_post_lower")]] <- row_quantile(uc_mat, q_lower)
  base[[paste0(gain_col, "_post_upper")]] <- row_quantile(uc_mat, q_upper)
  # Absolute posterior SD of the ranked value, on its native scale. This is the merit-DECOUPLED
  # spread the risk layer consumes; never a CV (usefulness can be ~0 or negative, so a ratio is
  # sign-ill-defined and would re-conflate merit with uncertainty).
  base$ranked_value_post_sd <- apply(uc_mat, 1L, function(r) {
    fr <- r[is.finite(r)]; if (length(fr) < 2L) NA_real_ else stats::sd(fr)
  })
  base$ranked_value_post_mean <- rowMeans(uc_mat, na.rm = TRUE)

  if (!is.null(tau_superior)) {
    tau <- as.numeric(tau_superior)
    k <- as.numeric(k_progeny)
    # Use the actual per-draw genetic family mean. Reconstructing mu from the
    # ranked value is invalid whenever gain_col/value_fun is mean-only,
    # threshold-penalized, multi-trait, or otherwise not exactly mu + i*sigma.
    pmv_pos <- pmax(pmv_mat, 0)
    sigma_mat <- sqrt(pmv_pos)
    p_mat <- matrix(NA_real_, nrow = n_pairs, ncol = S)
    for (s in seq_len(S)) {
      p_mat[, s] <- ng_p_superior_progeny(mu_mat[, s], sigma_mat[, s], tau, k)
    }
    base$p_superior_progeny_post_mean  <- rowMeans(p_mat, na.rm = TRUE)
    base$p_superior_progeny_post_lower <- row_quantile(p_mat, q_lower)
    base$p_superior_progeny_post_upper <- row_quantile(p_mat, q_upper)
    attr(base, "p_superior_progeny") <- list(
      tau = tau, k_progeny = k,
      mean_basis = "posterior_draw_cross_mean_gebv",
      variance_basis = var_col
    )
  }

  # Posterior top-N stability: for each N in top_n_targets, count fraction of
  # draws in which the cross is in the top-N by usefulness.
  for (N in as.integer(top_n_targets)) {
    if (!is.finite(N) || N < 1L || N >= n_pairs) next
    in_topn <- apply(uc_mat, 2L, function(col) {
      th <- sort(col, decreasing = TRUE, na.last = NA)[N]
      as.integer(col >= th & is.finite(col))
    })
    base[[paste0("posterior_topn_prob_", N)]] <- rowMeans(in_topn, na.rm = TRUE)
  }

  attr(base, "posterior") <- list(
    n_draws = S, method = posterior_effects$method,
    gain_col = gain_col, var_col = var_col,
    ranked_value = if (is.null(value_fun)) gain_col else "value_fun",
    ci_level = ci_level, selection_intensity = i_intensity,
    top_n_targets = as.integer(top_n_targets)
  )
  base
}

# ---- Posterior-aware OCS (D1c) --------------------------------------------

# Optimize the mating plan against a robustness quantile of the posterior
# usefulness instead of the point estimate. The NULL default uses the exact
# empirical lower credible bound already cached by ng_posterior_cross_predict().
# A different quantile requires a matching CI level in that prediction call;
# reconstructing it under a normal approximation is available only by explicit
# opt-in. Setting `objective = "posterior_topn_prob"` instead maximizes the sum
# of marginal per-cross top-N inclusion probabilities: the expected overlap
# with the posterior top-N set, not a joint probability for the whole plan.
ng_optimize_robust_mating_plan <- function(posterior_scores,
                                           n_crosses,
                                           parent_kinship = NULL,
                                           gain_col = "usefulness_pmv_gebv",
                                           robustness_quantile = NULL,
                                           allow_normal_approximation = FALSE,
                                           objective = c("posterior_quantile", "posterior_topn_prob"),
                                           top_n_target = NULL,
                                           max_crosses_per_parent = 4,
                                           min_crosses_per_parent = 0,
                                           min_unique_parents = NULL,
                                           max_pair_kinship = Inf,
                                           lambda_group = 0,
                                           lambda_mating = 0,
                                           lambda_parent_use = 0,
                                           lambda_parent_use_mode = c("absolute", "adaptive"),
                                           method = c("auto", "greedy_local", "repair_local", "mip_linear", "mip_contribution"),
                                           local_iter = 2000,
                                           ocs_iter = 5L) {
  objective <- match.arg(objective)
  if (!is.logical(allow_normal_approximation) ||
      length(allow_normal_approximation) != 1L || is.na(allow_normal_approximation)) {
    ng_stop("allow_normal_approximation must be TRUE or FALSE")
  }
  posterior_scores <- as.data.frame(posterior_scores, stringsAsFactors = FALSE)
  posterior_meta <- attr(posterior_scores, "posterior", exact = TRUE)
  posterior_ci <- suppressWarnings(as.numeric(posterior_meta$ci_level))
  if (is.null(posterior_meta) || length(posterior_ci) != 1L || !is.finite(posterior_ci) ||
      posterior_ci <= 0 || posterior_ci >= 1) {
    ng_stop("posterior_scores must retain posterior metadata from ng_posterior_cross_predict()")
  }
  posterior_meta$ci_level <- posterior_ci
  robust_col <- ".robust_gain"
  quantile_approximation <- FALSE
  if (identical(objective, "posterior_quantile")) {
    if (is.null(robustness_quantile)) {
      # Exact-by-default: use the lower tail already computed from the actual
      # posterior draws at ng_posterior_cross_predict()'s CI level.
      robustness_quantile <- (1 - posterior_meta$ci_level) / 2
    }
    robustness_quantile <- suppressWarnings(as.numeric(robustness_quantile))
    if (length(robustness_quantile) != 1L || !is.finite(robustness_quantile) ||
        robustness_quantile <= 0 || robustness_quantile >= 1) {
      ng_stop("robustness_quantile must be one finite probability in (0, 1)")
    }
    lower_col <- paste0(gain_col, "_post_lower")
    upper_col <- paste0(gain_col, "_post_upper")
    tail_prob <- (1 - posterior_meta$ci_level) / 2
    if (abs(robustness_quantile - tail_prob) <= 1e-12 &&
        lower_col %in% names(posterior_scores)) {
      # Reuse the cached lower CI when the user requested the same quantile.
      posterior_scores[[robust_col]] <- posterior_scores[[lower_col]]
    } else if (abs(robustness_quantile - (1 - tail_prob)) <= 1e-12 &&
               upper_col %in% names(posterior_scores)) {
      # The corresponding upper empirical quantile is cached as well.
      posterior_scores[[robust_col]] <- posterior_scores[[upper_col]]
    } else {
      # User asked for a different quantile; we need raw draws to recompute.
      # If not available on the table, fit a symmetric normal scale to the
      # cached central interval ONLY after explicit opt-in.
      mean_col <- paste0(gain_col, "_post_mean")
      if (all(c(mean_col, lower_col, upper_col) %in% names(posterior_scores))) {
        if (!isTRUE(allow_normal_approximation)) {
          matching_ci <- abs(1 - 2 * robustness_quantile)
          ng_stop(
            "Requested robustness_quantile is not an empirical tail quantile cached in ",
            "posterior_scores. Re-run ng_posterior_cross_predict() with ci_level = ",
            if (matching_ci > 0) format(matching_ci, digits = 8) else
              "a value whose empirical tail is the desired probability",
            " so that quantile is computed from draws, or explicitly set ",
            "allow_normal_approximation = TRUE."
          )
        }
        warning(
          "Robustness quantile is being reconstructed by a normal approximation; ",
          "posterior merit can be skewed. The plan summary records this approximation.",
          call. = FALSE
        )
        quantile_approximation <- TRUE
        z_target <- stats::qnorm(robustness_quantile)
        z_upper  <- stats::qnorm(1 - tail_prob)
        spread   <- (posterior_scores[[upper_col]] - posterior_scores[[lower_col]]) /
          (2 * z_upper)
        posterior_scores[[robust_col]] <- posterior_scores[[mean_col]] + z_target * spread
      } else {
        ng_stop("posterior_quantile objective requires *_post_* columns from ng_posterior_cross_predict()")
      }
    }
  } else {
    if (is.null(top_n_target)) ng_stop("posterior_topn_prob objective requires top_n_target")
    top_n_num <- suppressWarnings(as.numeric(top_n_target))
    if (length(top_n_num) != 1L || !is.finite(top_n_num) || top_n_num < 1 ||
        abs(top_n_num - round(top_n_num)) > 1e-8) {
      ng_stop("top_n_target must be one positive integer")
    }
    top_n_target <- as.integer(round(top_n_num))
    target_col <- paste0("posterior_topn_prob_", top_n_target)
    if (!(target_col %in% names(posterior_scores))) {
      ng_stop("posterior_scores missing column: ", target_col,
              ". Re-run ng_posterior_cross_predict() with top_n_targets including ",
              top_n_target)
    }
    topn_prob <- suppressWarnings(as.numeric(posterior_scores[[target_col]]))
    if (any(!is.finite(topn_prob)) || any(topn_prob < 0 | topn_prob > 1)) {
      ng_stop(target_col, " must contain finite probabilities in [0, 1]")
    }
    posterior_scores[[robust_col]] <- topn_prob
  }
  plan <- ng_optimize_mating_plan(
    scores = posterior_scores, n_crosses = n_crosses,
    gain_col = robust_col, parent_kinship = parent_kinship,
    max_crosses_per_parent = max_crosses_per_parent,
    min_crosses_per_parent = min_crosses_per_parent,
    min_unique_parents = min_unique_parents,
    max_pair_kinship = max_pair_kinship,
    lambda_group = lambda_group,
    lambda_mating = lambda_mating,
    lambda_parent_use = lambda_parent_use,
    lambda_parent_use_mode = lambda_parent_use_mode,
    method = method, local_iter = local_iter, ocs_iter = ocs_iter
  )
  s <- attr(plan, "summary")
  s$robust_objective <- objective
  s$robustness_quantile <- if (identical(objective, "posterior_quantile")) robustness_quantile else NA_real_
  s$robustness_quantile_is_normal_approximation <- quantile_approximation
  s$robust_top_n_target <- if (identical(objective, "posterior_topn_prob")) as.integer(top_n_target) else NA_integer_
  s$robust_gain_col <- gain_col
  attr(plan, "summary") <- s
  plan
}
