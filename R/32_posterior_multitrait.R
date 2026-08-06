# End-to-end multi-trait posterior cross prediction (v0.3.0 preview)
#
# Joins three v0.2.0 components in a single call:
#   - D1: per-trait posterior over marker effects (ng_fit_ridge_effects_posterior)
#   - D3: posterior over the genetic covariance G (ng_posterior_genetic_covariance)
#   - v0.1.0 multi-trait index solver (ng_add_multitrait_score, with G fed in)
#
# Returns a per-cross posterior on the multi-trait selection index along with
# rank-stability probabilities. This is the layer that lets a breeder ask:
# "for a multi-trait objective (yield up, disease down, lodging down), what's
# the credible interval on this candidate cross's index score, and how often
# does it land in the top-N across the joint posterior?"
#
# Computational design:
#   For each posterior draw s, we DON'T re-run the full ng_score_crosses
#   kernel per trait — that would be O(pairs * m^2) per draw for Kosambi/RIL.
#   Instead, we compute per-parent GEBVs once per draw per trait using the
#   posterior beta draws, then derive mid-parent values per cross in O(n_pairs).
#   The multi-trait index solve uses the per-draw G and is closed-form.
#   Total cost: O(n * m * t * n_draws) for GEBV updates + O(n_pairs * t * n_draws)
#   for mid-parent + O(t^3 * n_draws) for index solves. Tractable.
#
#   The default `value_mode = "mean"` uses cross_mean_gebv per trait — this is
#   the standard input to Smith-Hazel / Pesek-Baker indices, ignores within-
#   family variance, and is fast (~seconds for n=80, m=5000, t=3, S=200).
#
#   `value_mode = "usefulness"` opts in to a fuller PMV-aware calculation that
#   uses `cross_mean + i * sigma_t` per trait per draw, where sigma_t is the
#   per-trait posterior PMV under the BCM beta draws. This path requires
#   per-draw cross-scoring and is much slower; document the cost and warn.

ng_posterior_multitrait_cross_predict <- function(geno,
                                                  Y,
                                                  traits,
                                                  marker_map = NULL,
                                                  ids = rownames(geno),
                                                  pairs = NULL,
                                                  include_self = FALSE,
                                                  n_draws = 200L,
                                                  posterior_method = c("closed_form", "mcmc"),
                                                  genetic_covariance_method = c("beta_posterior", "parametric_bootstrap"),
                                                  phenotypic_covariance = NULL,
                                                  ridge_lambda = NULL,
                                                  kfold = 5L,
                                                  index_method = c("economic_index", "desired_gain", "auto", "weighted", "threshold"),
                                                  value_mode = c("mean", "usefulness"),
                                                  selection_prop = 0.10,
                                                  target = c("DH", "RIL"),
                                                  recomb_model = c("haldane", "kosambi"),
                                                  window_cm = Inf,
                                                  use_cpp = TRUE,
                                                  parent_type = c("inbred", "dh", "ril"),
                                                  assume_inbred = NULL,
                                                  tau_lower_vec = NULL,
                                                  tau_upper_vec = NULL,
                                                  threshold_k_progeny = 100L,
                                                  threshold_G_hat = NULL,
                                                  ci_level = 0.95,
                                                  top_n_targets = c(10L, 20L),
                                                  seed = 1L) {
  posterior_method <- match.arg(posterior_method)
  genetic_covariance_method <- match.arg(genetic_covariance_method)
  index_method <- match.arg(index_method)
  value_mode <- match.arg(value_mode)
  target <- match.arg(target)
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  recomb_model <- match.arg(recomb_model)

  # ---- Validate multivariate threshold arguments ---------------------------
  do_threshold <- !is.null(tau_lower_vec) || !is.null(tau_upper_vec)
  if (do_threshold) {
    if (is.null(tau_lower_vec) || is.null(tau_upper_vec)) {
      ng_stop("both tau_lower_vec and tau_upper_vec must be supplied together; ",
              "use -Inf / +Inf for unconstrained sides")
    }
    tau_lower_vec <- as.numeric(tau_lower_vec)
    tau_upper_vec <- as.numeric(tau_upper_vec)
    threshold_k_progeny <- as.integer(threshold_k_progeny)
    if (threshold_k_progeny < 1L)
      ng_stop("threshold_k_progeny must be a positive integer")
    if (!is.null(threshold_G_hat))
      threshold_G_hat <- as.matrix(threshold_G_hat)
    # Length validation deferred until n_traits is known (below).
  }

  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  if (nrow(Y) != nrow(geno)) ng_stop("Y must have one row per individual in geno")
  traits <- ng_multitrait_spec(traits)
  trait_cols <- traits$column
  trait_names <- traits$trait
  if (length(missing_traits <- setdiff(trait_cols, colnames(Y)))) {
    ng_stop("Y missing trait columns: ", paste(missing_traits, collapse = ", "))
  }
  Y <- Y[, trait_cols, drop = FALSE]
  colnames(Y) <- trait_cols
  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = include_self)
  pairs$parent1 <- as.character(pairs$parent1)
  pairs$parent2 <- as.character(pairs$parent2)
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  n_pairs <- nrow(pairs)
  n_traits <- length(trait_names)

  # ---- Validate tau vec lengths now that n_traits is known -----------------
  if (do_threshold) {
    if (length(tau_lower_vec) != n_traits)
      ng_stop("tau_lower_vec must have one entry per trait (length ", n_traits, ")")
    if (length(tau_upper_vec) != n_traits)
      ng_stop("tau_upper_vec must have one entry per trait (length ", n_traits, ")")
    if (!is.null(threshold_G_hat)) {
      if (!all(dim(threshold_G_hat) == n_traits))
        ng_stop("threshold_G_hat must be a ", n_traits, " x ", n_traits, " matrix")
    }
  }

  # ---- Fit per-trait posterior over beta -----------------------------------
  posteriors <- vector("list", n_traits)
  fits <- vector("list", n_traits)
  for (j in seq_len(n_traits)) {
    y_j <- as.numeric(Y[, j])
    names(y_j) <- ids
    post_j <- ng_fit_ridge_effects_posterior(
      geno = geno, y = y_j, ids = ids,
      lambda = ridge_lambda, kfold = kfold,
      n_draws = n_draws, method = posterior_method,
      seed = seed + j
    )
    posteriors[[j]] <- post_j
    fits[[j]] <- post_j$fit
  }
  names(posteriors) <- trait_names
  names(fits) <- trait_names

  # ---- Posterior G (S draws matched to the beta draws when possible) -------
  G_draws <- ng_posterior_genetic_covariance(
    geno = geno, Y = Y, n_draws = n_draws,
    method = genetic_covariance_method,
    ridge_lambda = ridge_lambda, kfold = kfold,
    seed = seed
  )
  # When genetic_covariance_method = "beta_posterior", G_draws shares the BCM
  # seed structure but uses seed + j internally per trait. Slices are still
  # valid independent draws even if not 1:1 matched to the trait posteriors.

  # ---- Phenotypic covariance (Smith-Hazel companion) -----------------------
  if (is.null(phenotypic_covariance)) {
    P_hat <- ng_estimate_phenotypic_covariance(Y, shrinkage = "auto")
  } else {
    P_hat <- as.matrix(phenotypic_covariance)
    dimnames(P_hat) <- list(trait_names, trait_names)
  }

  # ---- Pre-compute per-parent intercept and marker_mean -------------------
  intercepts <- vapply(fits, function(f) f$intercept, numeric(1L))
  marker_means <- vapply(fits, function(f) sum(f$marker_mean * f$beta), numeric(1L))

  # ---- Optional: per-draw within-family PMV (usefulness mode) -------------
  use_pmv <- identical(value_mode, "usefulness")
  if (use_pmv) {
    warning(
      "value_mode = 'usefulness' requires per-trait, per-draw cross-scoring (O(n_traits * n_draws * n_pairs * m)). For 80 parents x 5K markers x 3 traits x 200 draws expect ~5-15 minutes on one core; consider reducing n_draws or LD-pruning markers.",
      call. = FALSE
    )
    intensity <- ng_selection_intensity(selection_prop)
  }
  if (do_threshold && isFALSE(use_pmv)) {
    warning(
      "tau_lower_vec / tau_upper_vec require per-trait PMV computation; ",
      "forcing PMV scoring this run. Expect runtime closer to value_mode = ",
      "'usefulness' than 'mean'.",
      call. = FALSE
    )
  }

  # ---- Per-draw cross-level multi-trait index ------------------------------
  # Cache the rank-normalization and trait centering across draws to keep
  # the inner loop cheap. We rebuild the per-cross trait matrix per draw and
  # call ng_add_multitrait_score with G_s.
  index_mat <- matrix(NA_real_, nrow = n_pairs, ncol = n_draws)
  # Allocate threshold probability accumulator when tau args are supplied.
  if (do_threshold) {
    p_mt_mat <- matrix(NA_real_, nrow = n_pairs, ncol = n_draws)
  }
  ci_alpha <- (1 - ci_level) / 2

  for (s in seq_len(n_draws)) {
    # Per-trait posterior cross prediction.
    trait_cross_value <- matrix(NA_real_, nrow = n_pairs, ncol = n_traits)
    colnames(trait_cross_value) <- trait_cols
    # Per-trait per-pair PMV storage (needed by threshold path even in mean mode).
    if (do_threshold) {
      trait_cross_var <- matrix(NA_real_, nrow = n_pairs, ncol = n_traits)
    }
    for (j in seq_len(n_traits)) {
      beta_js <- posteriors[[j]]$beta_draws[, s]
      # GEBV per parent under draw s: intercept_j shifted to match the centered
      # marker-mean used during fitting (see ng_score_crosses commentary).
      gebv_js <- intercepts[[j]] + sum(fits[[j]]$marker_mean * fits[[j]]$beta) -
        sum(fits[[j]]$marker_mean * beta_js) + as.numeric(geno %*% beta_js)
      mp_js <- 0.5 * (gebv_js[p1] + gebv_js[p2])
      if (use_pmv || do_threshold) {
        # Per-pair PMV under draw s: VPM via the existing recursion.
        scored_pair <- ng_dh_recomb_variance_pairs(
          geno = geno, beta = beta_js, beta_var = rep(0, ncol(geno)),
          marker_map = ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model),
          ids = ids, pairs = pairs, window_cm = window_cm,
          use_cpp = use_cpp, recomb_model = recomb_model, target = target
        )
        pmv_js <- pmax(scored_pair$pmv, 0)
        if (do_threshold) trait_cross_var[, j] <- pmv_js
        if (use_pmv) {
          sigma_js <- sqrt(pmv_js)
          trait_cross_value[, j] <- mp_js + intensity * sigma_js
        } else {
          trait_cross_value[, j] <- mp_js
        }
      } else {
        trait_cross_value[, j] <- mp_js
      }
    }
    # Per-draw multi-trait index.
    scores_s <- data.frame(
      parent1 = pairs$parent1, parent2 = pairs$parent2,
      stringsAsFactors = FALSE
    )
    for (j in seq_len(n_traits)) {
      scores_s[[trait_cols[[j]]]] <- trait_cross_value[, j]
    }
    G_s <- G_draws[, , s]
    scored_s <- ng_add_multitrait_score(
      scores = scores_s, traits = traits, method = index_method,
      phenotypic_covariance = P_hat, genetic_covariance = G_s,
      threshold_penalty_autoscale = TRUE
    )
    index_mat[, s] <- as.numeric(scored_s$multi_trait_score)

    # Per-draw multivariate threshold probability (one value per cross).
    if (do_threshold) {
      for (i in seq_len(n_pairs)) {
        mu_is <- trait_cross_value[i, ]
        var_is <- trait_cross_var[i, ]
        Sigma_i <- ng_build_cross_trait_covariance(
          per_trait_var = var_is, G_hat = threshold_G_hat
        )
        p_mt_mat[i, s] <- ng_p_superior_progeny_multitrait(
          mu = mu_is, Sigma_c = Sigma_i,
          tau_lower = tau_lower_vec, tau_upper = tau_upper_vec,
          k_progeny = threshold_k_progeny
        )
      }
    }
  }

  # ---- Aggregate across draws ---------------------------------------------
  q_lower <- ci_alpha
  q_upper <- 1 - ci_alpha
  row_quantile <- function(M, q) {
    apply(M, 1L, function(row) {
      finite_row <- row[is.finite(row)]
      if (!length(finite_row)) return(NA_real_)
      unname(stats::quantile(finite_row, probs = q, names = FALSE))
    })
  }

  out <- pairs
  out$multi_trait_score_post_mean  <- rowMeans(index_mat, na.rm = TRUE)
  out$multi_trait_score_post_lower <- row_quantile(index_mat, q_lower)
  out$multi_trait_score_post_upper <- row_quantile(index_mat, q_upper)

  for (N in as.integer(top_n_targets)) {
    if (!is.finite(N) || N < 1L || N >= n_pairs) next
    in_topn <- apply(index_mat, 2L, function(col) {
      th <- sort(col, decreasing = TRUE, na.last = NA)[N]
      as.integer(col >= th & is.finite(col))
    })
    out[[paste0("multitrait_posterior_topn_prob_", N)]] <- rowMeans(in_topn, na.rm = TRUE)
  }

  # ---- Multivariate threshold posterior summary ----------------------------
  if (do_threshold) {
    out$p_superior_progeny_mt_post_mean  <- rowMeans(p_mt_mat, na.rm = TRUE)
    out$p_superior_progeny_mt_post_lower <- row_quantile(p_mt_mat, q_lower)
    out$p_superior_progeny_mt_post_upper <- row_quantile(p_mt_mat, q_upper)
  }

  attr(out, "posterior_multitrait") <- list(
    n_draws = n_draws,
    posterior_method = posterior_method,
    genetic_covariance_method = genetic_covariance_method,
    index_method = index_method,
    value_mode = value_mode,
    ci_level = ci_level,
    top_n_targets = as.integer(top_n_targets),
    G_mean = attr(G_draws, "G_mean"),
    P_hat = P_hat,
    tau_lower_vec = tau_lower_vec,
    tau_upper_vec = tau_upper_vec,
    threshold_k_progeny = if (do_threshold) threshold_k_progeny else NULL,
    threshold_G_hat = threshold_G_hat
  )
  out
}
