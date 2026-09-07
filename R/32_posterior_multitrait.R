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
#
# I2 fix: when tau_lower_vec/tau_upper_vec are supplied (do_threshold = TRUE), the per-draw,
# per-pair loop below calls ng_p_superior_progeny_multitrait() (R/33_threshold_probability_
# multitrait.R), which calls mvtnorm::pmvnorm(). At dimension >= 3 traits pmvnorm() switches
# from a deterministic closed form to the randomised GenzBretz lattice rule, which both DRAWS
# FROM and ADVANCES the ambient .Random.seed -- once per (draw, pair) call, so up to
# n_draws * n_pairs times per invocation. Left unguarded, this function would (a) perturb
# whatever the caller's ambient RNG stream happens to be at >= 3 traits, and (b) not be
# reproducible run to run on its own account, since its result would depend on that ambient
# state too. Guarded below by wrapping the ENTIRE per-draw loop -- not a per-call or per-draw
# reset, which would restart pmvnorm's lattice-shift sequence from the same point on every
# call/draw instead of letting it evolve continuously -- in a single ng_with_rng_seed(mt_pmvnorm_seed,
# ...) (R/00_utils.R) call that spans every pmvnorm() invocation this function makes. Nothing
# else in the loop consumes the ambient RNG (the posterior beta draws are already deterministic,
# precomputed columns of ng_fit_ridge_effects_posterior()'s beta_draws matrix, indexed -- not
# resampled -- per s), so this also makes p_superior_progeny_mt_post_* reproducible run to run at
# >= 3 traits, which it was not before. This DOES change the posterior threshold probabilities
# numerically at >= 3 traits vs. the unfixed code (the MC lattice noise it now integrates from is
# different) -- see the fix report for the measured size of that shift. At <= 2 traits pmvnorm is
# deterministic, so the seed is inert and output is bit-identical before/after.
NG_POSTERIOR_MT_THRESHOLD_PMVNORM_SEED <- 20260906L

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
                                                  genetic_covariance = NULL,
                                                  genetic_covariance_draws = NULL,
                                                  ridge_lambda = NULL,
                                                  kfold = 5L,
                                                  index_method = c("auto", "economic_index", "desired_gain", "weighted", "threshold"),
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
                                                  robustness_quantile = NULL,
                                                  top_n_targets = c(10L, 20L),
                                                  seed = 1L) {
  posterior_method <- match.arg(posterior_method)
  genetic_covariance_method <- match.arg(genetic_covariance_method)
  index_method <- match.arg(index_method)
  value_mode <- match.arg(value_mode)
  target <- match.arg(target)
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  recomb_model <- match.arg(recomb_model)

  # `ci_level` and `robustness_quantile` are SEPARATE controls, exactly as in
  # ng_posterior_cross_predict() (R/30): ci_level owns the reported credible interval
  # (_post_lower / _post_upper), robustness_quantile owns the tail an allocation is optimised
  # against. Both q and 1 - q are cached from the SAME draws, so a robust allocation on the
  # index is exact at that quantile with no extra sampling and no normal approximation, and the
  # reported interval is untouched.
  ci_level <- suppressWarnings(as.numeric(ci_level))
  if (length(ci_level) != 1L || !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1) {
    ng_stop("ci_level must be one finite probability in (0, 1)")
  }
  robust_probs <- numeric(0)
  if (!is.null(robustness_quantile)) {
    robust_probs <- suppressWarnings(as.numeric(robustness_quantile))
    if (!length(robust_probs) || any(!is.finite(robust_probs)) ||
        any(robust_probs <= 0) || any(robust_probs >= 1)) {
      ng_stop("robustness_quantile must be finite probabilities in (0, 1)")
    }
  }

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
  # Canonicalize marker order BEFORE fitting/sampling. This is required not
  # only by the linkage recursion later: posterior draws consume random
  # normals in column order, so leaving an arbitrary input order would give
  # different finite-draw summaries for the same named marker data. Sorting
  # genotype and map together makes the public result invariant to file-column
  # order under a fixed seed.
  if (!is.null(marker_map)) {
    prepared_map <- ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model)
    canonical <- ng_sort_by_map(
      geno, stats::setNames(rep(0, ncol(geno)), colnames(geno)),
      rep(0, ncol(geno)), marker_map = prepared_map
    )
    geno <- canonical$geno
    marker_map <- canonical$marker_map
  }
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
    # SEED (0.28.0): identity-derived, not position-derived. `seed + j` keyed each trait's
    # posterior stream to its ROW POSITION in the caller's traits table, so permuting that
    # table -- i.e. reordering the breeder's direction file -- silently moved every trait onto
    # a different draw stream and changed the index posterior it feeds. The streams must stay
    # DISTINCT per trait (shared innovations would manufacture cross-trait correlation in
    # exactly the index this function builds), so the key is the trait's name, not its row.
    post_j <- ng_fit_ridge_effects_posterior(
      geno = geno, y = y_j, ids = ids,
      lambda = ridge_lambda, kfold = kfold,
      n_draws = n_draws, method = posterior_method,
      seed = ng_trait_rng_seed(seed, trait_names[[j]])
    )
    posteriors[[j]] <- post_j
    fits[[j]] <- post_j$fit
  }
  names(posteriors) <- trait_names
  names(fits) <- trait_names

  # A formal economic/desired-gain index needs a defensible G. Independently
  # refitting univariate ridge models and correlating their marker effects is
  # not a multivariate variance-component posterior, so this routine no longer
  # manufactures G draws internally. Supply a REML/Bayesian G or matched draws.
  has_desired <- any(is.finite(suppressWarnings(as.numeric(traits$desired_change))) &
                       suppressWarnings(as.numeric(traits$desired_change)) > 0)
  has_economic <- any(is.finite(suppressWarnings(as.numeric(traits$economic_weight))) &
                        suppressWarnings(as.numeric(traits$economic_weight)) > 0)
  has_weight <- any(is.finite(suppressWarnings(as.numeric(traits$weight))) &
                      suppressWarnings(as.numeric(traits$weight)) > 0)
  resolved_index_method <- index_method
  if (identical(index_method, "auto")) {
    resolved_index_method <- if (has_desired) "desired_gain" else if (has_economic) {
      "economic_index"
    } else if (has_weight) "weighted" else "auto"
  }
  formal_index <- resolved_index_method %in% c("economic_index", "desired_gain")
  if (formal_index && is.null(genetic_covariance) && is.null(genetic_covariance_draws)) {
    ng_stop("posterior economic/desired-gain prediction requires genetic_covariance or ",
            "matched genetic_covariance_draws from a multivariate model")
  }
  if (!is.null(genetic_covariance_draws)) {
    genetic_covariance_draws <- as.array(genetic_covariance_draws)
    if (length(dim(genetic_covariance_draws)) != 3L ||
        !all(dim(genetic_covariance_draws)[1:2] == n_traits) ||
        dim(genetic_covariance_draws)[3] != n_draws) {
      ng_stop("genetic_covariance_draws must be n_traits x n_traits x n_draws")
    }
  }

  # ---- Phenotypic covariance (Smith-Hazel companion) -----------------------
  if (is.null(phenotypic_covariance)) {
    P_hat <- ng_estimate_phenotypic_covariance(Y, shrinkage = "auto")
  } else {
    # 0.27.0: this used to STAMP dimnames positionally, which silently relabelled a
    # user matrix whose rows/columns were in a different order (and defeated the
    # by-name alignment downstream). Align by name instead; an unlabelled matrix is
    # still read positionally in traits$trait order. See ng_multitrait_align_cov().
    P_hat <- ng_multitrait_align_cov(phenotypic_covariance, trait_names,
                                     "phenotypic_covariance")
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

  # I2 fix: see the header comment above ng_posterior_multitrait_cross_predict() for why this
  # scopes the WHOLE loop (every pmvnorm() call this function makes) under one fixed seed rather
  # than resetting per draw or per pair. mt_pmvnorm_seed is NULL (a documented no-op for
  # ng_with_rng_seed()) unless do_threshold is TRUE, since no pmvnorm() call happens otherwise.
  mt_pmvnorm_seed <- if (do_threshold) NG_POSTERIOR_MT_THRESHOLD_PMVNORM_SEED else NULL
  ng_with_rng_seed(mt_pmvnorm_seed, {
  for (s in seq_len(n_draws)) {
    # Per-trait posterior cross prediction.
    trait_cross_value <- matrix(NA_real_, nrow = n_pairs, ncol = n_traits)
    trait_cross_mean <- matrix(NA_real_, nrow = n_pairs, ncol = n_traits)
    colnames(trait_cross_value) <- trait_cols
    colnames(trait_cross_mean) <- trait_cols
    beta_s <- matrix(NA_real_, nrow = ncol(geno), ncol = n_traits,
                     dimnames = list(colnames(geno), trait_names))
    # Per-trait per-pair PMV storage (needed by threshold path even in mean mode).
    if (do_threshold) {
      trait_cross_var <- matrix(NA_real_, nrow = n_pairs, ncol = n_traits)
    }
    for (j in seq_len(n_traits)) {
      beta_js <- posteriors[[j]]$beta_draws[, s]
      beta_s[, j] <- beta_js
      # GEBV per parent under draw s: intercept_j shifted to match the centered
      # marker-mean used during fitting (see ng_score_crosses commentary).
      gebv_js <- intercepts[[j]] + sum(fits[[j]]$marker_mean * fits[[j]]$beta) -
        sum(fits[[j]]$marker_mean * beta_js) + as.numeric(geno %*% beta_js)
      mp_js <- 0.5 * (gebv_js[p1] + gebv_js[p2])
      trait_cross_mean[, j] <- mp_js
      if (use_pmv || do_threshold) {
        # Per-pair PMV under draw s: VPM via the existing recursion.
        # Prepare AND sort together. ng_prepare_marker_map() aligns map rows to
        # the caller's marker columns but deliberately preserves that order;
        # the chromosome recursion requires chromosome/position order.
        map_js <- ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model)
        sorted_js <- ng_sort_by_map(
          geno, beta_js, rep(0, ncol(geno)), marker_map = map_js
        )
        scored_pair <- ng_dh_recomb_variance_pairs(
          geno = sorted_js$geno, beta = sorted_js$effects,
          beta_var = sorted_js$beta_var, marker_map = sorted_js$marker_map,
          ids = ids, pairs = pairs, window_cm = window_cm,
          use_cpp = use_cpp, recomb_model = recomb_model, target = target
        )
        pmv_js <- pmax(scored_pair$pmv, 0)
        if (do_threshold) trait_cross_var[, j] <- pmv_js
        if (use_pmv) {
          sigma_js <- sqrt(pmv_js)
          # DIRECTION SIGN (0.26.0 fix). The usefulness criterion is
          # mean + sign * i * sigma with sign = +1 for a maximize trait and -1 for a minimize
          # trait, and the sign applies ONLY to the i*sigma term -- never to the mean, which is
          # already in the trait's own units and keeps its own orientation. Without the sign a
          # minimize trait (disease, lodging) was scored at mean + i*sigma, i.e. the UNFAVOURABLE
          # tail of the family, so within-family variance was charged as a liability instead of
          # credited as an opportunity, and then value_z applied -1 on top of that. This matches
          # ng_run_cp_trait_value() (R/39_cross_prediction_runner.R), which is the wired
          # single-trait convention; the two must not disagree.
          sign_j <- if (identical(traits$direction[[j]], "maximize")) 1 else -1
          trait_cross_value[, j] <- mp_js + sign_j * intensity * sigma_js
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
    G_s <- if (!is.null(genetic_covariance_draws)) {
      genetic_covariance_draws[, , s]
    } else {
      genetic_covariance
    }
    scored_s <- ng_add_multitrait_score(
      scores = scores_s, traits = traits, method = resolved_index_method,
      phenotypic_covariance = P_hat, genetic_covariance = G_s,
      threshold_penalty_autoscale = TRUE
    )
    index_mat[, s] <- as.numeric(scored_s$multi_trait_score)

    # Per-draw multivariate threshold probability (one value per cross).
    if (do_threshold) {
      exact_sigma <- NULL
      if (is.null(threshold_G_hat)) {
        wf <- ng_cross_trait_within_family_cov(
          geno = geno, betas = beta_s, marker_map = marker_map,
          ids = ids, pairs = pairs, target = target,
          recomb_model = recomb_model, window_cm = window_cm
        )
        exact_sigma <- ng_multitrait_exact_sigma_list(
          scores = pairs, trait_order = trait_names, cross_trait_cov = wf
        )
      }
      for (i in seq_len(n_pairs)) {
        mu_is <- trait_cross_mean[i, ]
        var_is <- trait_cross_var[i, ]
        Sigma_i <- if (!is.null(exact_sigma)) exact_sigma[[i]] else
          ng_build_cross_trait_covariance(per_trait_var = var_is, G_hat = threshold_G_hat)
        p_mt_mat[i, s] <- ng_p_superior_progeny_multitrait(
          mu = mu_is, Sigma_c = Sigma_i,
          tau_lower = tau_lower_vec, tau_upper = tau_upper_vec,
          k_progeny = threshold_k_progeny
        )
      }
    }
  }
  })

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
  # PER-DRAW RE-STANDARDISATION CAVEAT -- read before displaying this interval to a breeder.
  # ng_add_multitrait_score() is called INSIDE the draw loop above (that is what makes the
  # index a coherent joint draw rather than a recombination of per-trait marginals), and it
  # re-derives its own rank-normalisation / IQR centring from the rows it is handed. The index
  # is therefore re-standardised within EVERY draw, so any component of posterior uncertainty
  # that shifts or rescales the whole candidate pool together is removed before the quantile is
  # taken. Consequences, all deliberate and out of scope to "fix" here (removing the
  # re-standardisation is a design change to the index itself, not a bug fix):
  #   * multi_trait_score_post_lower/_upper are an interval on a PER-DRAW RELATIVE index, not
  #     on any fixed-scale quantity. They are narrower than a genuine index credible interval.
  #   * They are not on the same scale as, and must not be differenced against, the
  #     point-estimate multi_trait_score a runner call reports.
  #   * For RANK-stability (multitrait_posterior_topn_prob_*) and for ordering crosses by a
  #     conservative tail this is harmless, and arguably what is wanted.
  # Recorded in the metadata below as `index_rescaling` / `index_rescaling_note` so a frontend
  # can badge it rather than presenting it as a fixed-scale credible interval.
  out$multi_trait_score_post_mean  <- rowMeans(index_mat, na.rm = TRUE)
  out$multi_trait_score_post_lower <- row_quantile(index_mat, q_lower)
  out$multi_trait_score_post_upper <- row_quantile(index_mat, q_upper)
  # Absolute posterior SD of the index, mirroring ranked_value_post_sd in R/30. This is the
  # merit-DECOUPLED spread the cross-priority risk layer consumes; never a CV (the index can be
  # ~0 or negative, so a ratio is sign-ill-defined and would re-conflate merit with uncertainty).
  out$multi_trait_score_post_sd <- apply(index_mat, 1L, function(r) {
    fr <- r[is.finite(r)]; if (length(fr) < 2L) NA_real_ else stats::sd(fr)
  })

  # Extra empirical quantile(s) of the index, taken from the SAME index_mat draws that produced
  # the credible interval above -- the treatment 0.25.0 gave ng_posterior_cross_predict(). Both
  # q and 1 - q are cached so ng_optimize_robust_mating_plan() is exact in either orientation
  # without borrowing ci_level and without a normal approximation.
  quantile_probs <- if (length(robust_probs)) sort(unique(c(robust_probs, 1 - robust_probs))) else numeric(0)
  posterior_quantiles <- NULL
  if (length(quantile_probs)) {
    quantile_cols <- vapply(quantile_probs,
                            function(p) ng_posterior_quantile_col("multi_trait_score", p),
                            character(1L))
    for (jj in seq_along(quantile_probs)) {
      out[[quantile_cols[[jj]]]] <- row_quantile(index_mat, quantile_probs[[jj]])
    }
    posterior_quantiles <- data.frame(prob = quantile_probs, column = quantile_cols,
                                      stringsAsFactors = FALSE)
  }

  # INDEX ORIENTATION -- hard-coded, not an argument, and deliberately so.
  # Unlike the per-trait ranked value in R/30 (whose orientation depends on the trait's
  # breeding direction and on trait_value_metric), multi_trait_score is direction-normalised
  # to HIGHER = BETTER for every index method, and that is a stated design invariant of
  # ng_add_multitrait_score() (see the "Design invariant (deliberate; do not normalize this
  # away)" comment in R/19_multi_trait_selection.R). Both combination paths apply the trait
  # direction UPSTREAM of the combination: the rank path negates a minimize trait inside
  # ng_rank_normalize(bigger_is_better = FALSE) and combines with weights >= 0 summing to 1;
  # the solved path builds value_z = sign * (x - center) / scale and combines with b = P^-1 G a
  # / b = G^-1 d whose a, d >= 0 are enforced upstream. A `direction` argument here could only
  # ever be "maximize", so plumbing one would be a meaningless knob that a caller could set
  # wrong. The top-N sort below and the "posterior" metadata both fix it at "maximize".
  index_direction <- "maximize"
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
    genetic_covariance_method = if (!is.null(genetic_covariance_draws)) "user_matched_draws" else
      if (!is.null(genetic_covariance)) "user_fixed" else NA_character_,
    index_method = resolved_index_method,
    index_method_requested = index_method,
    value_mode = value_mode,
    ci_level = ci_level,
    direction = index_direction,
    robustness_quantile = if (length(robust_probs)) robust_probs else NA_real_,
    posterior_quantiles = posterior_quantiles,
    # Surfaced so the frontend can badge the interval instead of presenting it as a
    # fixed-scale credible interval; see the block comment above where the columns are built.
    index_rescaling = "per_draw_restandardized",
    index_rescaling_note = paste0(
      "multi_trait_score is re-standardized within every posterior draw (the index is ",
      "recomputed inside the draw loop by ng_add_multitrait_score(), which re-derives its ",
      "rank-normalization / IQR centring from the rows it is handed). ",
      "multi_trait_score_post_lower/_upper are therefore an interval on a PER-DRAW RELATIVE ",
      "index: pool-wide shifts and rescalings are removed before the quantile is taken, the ",
      "interval is narrower than a genuine index credible interval, and it is not on the same ",
      "scale as the point-estimate multi_trait_score. Rank-stability ",
      "(multitrait_posterior_topn_prob_*) and conservative-tail ORDERING are unaffected."),
    top_n_targets = as.integer(top_n_targets),
    G_mean = if (!is.null(genetic_covariance_draws)) apply(genetic_covariance_draws, c(1, 2), mean) else genetic_covariance,
    P_hat = P_hat,
    tau_lower_vec = tau_lower_vec,
    tau_upper_vec = tau_upper_vec,
    threshold_k_progeny = if (do_threshold) threshold_k_progeny else NULL,
    threshold_G_hat = threshold_G_hat
  )
  # Second, NARROWER attribute in the shape ng_optimize_robust_mating_plan() (R/30) reads, so
  # a robust allocation can be run on the INDEX the multi-trait plan actually ranks on instead
  # of on one trait's per-trait posterior. gain_col is "multi_trait_score"; direction is fixed
  # at "maximize" (see the index-orientation comment above); posterior_quantiles points at the
  # exact empirical quantile columns cached from these same draws. Additive -- the richer
  # "posterior_multitrait" attribute above is unchanged and remains the authoritative record.
  attr(out, "posterior") <- list(
    n_draws = n_draws,
    method = posterior_method,
    gain_col = "multi_trait_score",
    var_col = NA_character_,
    ranked_value = "multi_trait_score",
    ci_level = ci_level,
    top_n_targets = as.integer(top_n_targets),
    direction = index_direction,
    robustness_quantile = if (length(robust_probs)) robust_probs else NA_real_,
    posterior_quantiles = posterior_quantiles,
    index_rescaling = "per_draw_restandardized"
  )
  out
}
