ng_score_crosses <- function(geno,
                             effects,
                             marker_map = NULL,
                             ids = rownames(geno),
                             pairs = NULL,
                             adjusted_pheno = NULL,
                             blue = NULL,
                             blup = NULL,
                             include_self = FALSE,
                             target = c("DH", "RIL"),
                             selection_prop = 0.10,
                             min_effect_reliability = 0.35,
                             recomb_model = c("haldane", "kosambi"),
                             window_cm = Inf,
                             use_cpp = TRUE,
                             assume_inbred = TRUE,
                             inbred_tolerance = 0.05,
                             inbred_marker_fraction = 0.02,
                             posterior_cov_full = NULL,
                             parent_kinship = NULL,
                             grm_method = c("vanraden", "yang")) {
  target <- match.arg(target)
  recomb_model <- match.arg(recomb_model)
  grm_method <- match.arg(grm_method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  # The DH/RIL recombination-variance formula assumes both parents are fully
  # inbred (marker dosage in {0, 2}). When parents carry residual heterozygosity
  # the formula Cov(y_k, y_l) = d_k d_l (1 - 2 r_kl) no longer holds because the
  # F1 may be homozygous at some heterozygous loci. Refuse to silently produce
  # wrong PMV: detect non-inbred dosages and either error out or downgrade to
  # the relationship-distance baseline.
  inbred_audit <- ng_audit_inbred_dosage(
    geno, tolerance = inbred_tolerance,
    fraction_tolerance = inbred_marker_fraction
  )
  if (isTRUE(assume_inbred) && length(inbred_audit$violators)) {
    msg <- sprintf(
      "%d / %d parents exceed the residual-heterozygosity tolerance (max het-marker frac = %.3f, examples: %s). The DH/RIL recombination-variance kernel assumes inbred parents; set assume_inbred = FALSE and supply phased haplotypes to ng_exact_gms_additive_var() for outbred parents.",
      length(inbred_audit$violators), nrow(geno), inbred_audit$max_fraction,
      paste(head(inbred_audit$violators, 4L), collapse = ", ")
    )
    ng_stop(msg)
  } else if (!isTRUE(assume_inbred) && length(inbred_audit$violators)) {
    warning(sprintf(
      "assume_inbred = FALSE: %d parents are not fully inbred. vpm / pmv below treat parents as inbreds and will be biased; use ng_exact_gms_additive_var() with phased haplotypes for the exact outbred path.",
      length(inbred_audit$violators)
    ), call. = FALSE)
  }
  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = include_self)
  pairs$parent1 <- as.character(pairs$parent1)
  pairs$parent2 <- as.character(pairs$parent2)

  mean_source <- ng_choose_mean_source(
    geno = geno,
    effects = effects,
    adjusted_pheno = adjusted_pheno,
    blue = blue,
    blup = blup,
    ids = ids,
    min_reliability = min_effect_reliability
  )
  # Reuse a caller-supplied kinship matrix when given (avoids recomputing the O(n^2*m)
  # VanRaden G on identical genotypes across a per-trait scoring loop); else compute it.
  K <- if (is.null(parent_kinship)) {
    ng_parent_kinship(geno, method = grm_method)
  } else {
    pk <- as.matrix(parent_kinship)
    if (!is.null(rownames(pk)) && all(ids %in% rownames(pk))) pk[ids, ids, drop = FALSE] else pk
  }
  rel_var <- ng_pair_relationship_variance(pairs, K)

  beta <- effects$beta[colnames(geno)]
  beta[!is.finite(beta)] <- 0
  beta_var <- effects$beta_var[colnames(geno)]
  beta_var[!is.finite(beta_var)] <- 0
  marker_map <- ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model)
  sorted <- ng_sort_by_map(geno, beta, beta_var, marker_map)

  # Optional opt-in: full off-diagonal posterior PMV. Routes to a dense O(m^2)
  # kernel that takes the m x m Sigma_beta supplied by the caller (typically
  # from ng_fit_ridge_effects(return_beta_cov_full = TRUE)) and returns an extra
  # pmv_full_posterior column alongside the legacy diagonal one.
  if (!is.null(posterior_cov_full)) {
    if (!is.matrix(posterior_cov_full) ||
        nrow(posterior_cov_full) != ncol(posterior_cov_full)) {
      ng_stop("posterior_cov_full must be a square m x m matrix")
    }
    if (is.null(rownames(posterior_cov_full)) || is.null(colnames(posterior_cov_full))) {
      if (nrow(posterior_cov_full) != ncol(sorted$geno)) {
        ng_stop("posterior_cov_full has no dimnames and does not match colnames(geno) in length")
      }
      rownames(posterior_cov_full) <- colnames(geno)
      colnames(posterior_cov_full) <- colnames(geno)
    }
    missing_markers <- setdiff(sorted$marker_map$marker, rownames(posterior_cov_full))
    if (length(missing_markers)) {
      ng_stop(sprintf(
        "posterior_cov_full is missing rows/cols for %d marker(s): %s",
        length(missing_markers),
        paste(head(missing_markers, 4L), collapse = ", ")
      ))
    }
    dh <- ng_dh_recomb_variance_pairs_full_posterior(
      geno = sorted$geno,
      beta = sorted$effects,
      beta_cov_full = posterior_cov_full,
      marker_map = sorted$marker_map,
      ids = ids,
      pairs = pairs,
      target = target,
      recomb_model = recomb_model,
      window_cm = window_cm
    )
  } else {
    dh <- ng_dh_recomb_variance_pairs(
      geno = sorted$geno,
      beta = sorted$effects,
      beta_var = sorted$beta_var,
      marker_map = sorted$marker_map,
      ids = ids,
      pairs = pairs,
      window_cm = window_cm,
      use_cpp = use_cpp,
      target = target,
      recomb_model = recomb_model
    )
  }

  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  parent_mean <- 0.5 * (mean_source$value[p1] + mean_source$value[p2])
  adjusted_mean <- rep(NA_real_, nrow(pairs))
  if (!is.null(adjusted_pheno)) {
    adj <- ng_match_vector(adjusted_pheno, ids, "adjusted_pheno")
    adjusted_mean <- 0.5 * (adj[p1] + adj[p2])
  }
  gebv <- setNames(ng_predict_gebv(geno, effects), ids)
  mpv_gebv <- 0.5 * (gebv[p1] + gebv[p2])
  i <- ng_selection_intensity(selection_prop)
  dh_var <- pmax(dh$vpm, 0)
  dh_pmv <- pmax(dh$pmv, 0)
  effect_rel <- max(0, min(1, mean_source$reliability))
  # Blend mean sources with a smooth weight in the marker-effect reliability.
  # This mixes BLUEs/adjusted phenotypes (no marker-effect uncertainty) with
  # mid-parental GEBVs (genomic prediction), all in genetic-value units.
  blended_mean <- effect_rel * mpv_gebv + (1 - effect_rel) * adjusted_mean
  blended_mean[!is.finite(blended_mean)] <- parent_mean[!is.finite(blended_mean)]

  out <- pairs
  out$mean_source <- mean_source$source
  out$effect_reliability <- effect_rel
  out$progeny_target <- target
  out$mid_parent_value <- as.numeric(mpv_gebv)
  out$cross_mean <- as.numeric(parent_mean)
  out$cross_mean_gebv <- as.numeric(mpv_gebv)
  out$cross_mean_adj <- as.numeric(adjusted_mean)
  out$cross_mean_blend <- as.numeric(blended_mean)
  # Relationship-distance metric (NOT a variance). Reported separately so the
  # OCS optimizer can use it as a diversity penalty; never combined with PMV.
  out$parent_distance <- rel_var$parent_distance
  out$pair_kinship <- rel_var$pair_kinship
  # Expected inbreeding of the immediate progeny of each cross = coancestry of the two
  # parents (Meuwissen OCS: F_progeny = f(p1, p2)). pair_kinship is read off a VanRaden
  # genomic relationship matrix G, and G_XY ~ numerator relationship ~ 2 x coancestry
  # for non-inbred parents, so F_progeny ~ pair_kinship / 2. This is the pedigree/
  # relationship-scale progeny inbreeding managed in tactical mate selection; it is distinct
  # from the eventual homozygosity of a finished DH/RIL line (which tends to 1
  # regardless of the parents). Reported so the optimizer and reports can treat progeny
  # inbreeding as a first-class quantity alongside parental (group) coancestry.
  out$expected_progeny_inbreeding <- pmax(0, as.numeric(rel_var$pair_kinship) / 2)
  # Within-family genetic-value variances under the requested progeny target
  # (DH or RIL), in (genetic value)^2 units. pmv adds the diagonal
  # posterior marker-effect uncertainty to vpm (VPM). When the
  # caller supplies posterior_cov_full, pmv_full_posterior also
  # carries the full off-diagonal correction d'(R o Sigma_beta) d
  # (genomicMateSelectR formulation); NA otherwise.
  out$vpm <- dh_var
  out$pmv <- dh_pmv
  if (!is.null(dh$pmv_full_posterior)) {
    out$pmv_full_posterior <- pmax(dh$pmv_full_posterior, 0)
  } else {
    out$pmv_full_posterior <- NA_real_
  }
  # Usefulness criteria are mu + i * sigma in genetic-value units. Two
  # variance choices (VPM = vpm; PMV = pmv) crossed with four
  # mean sources. Pre-v0.1.0 columns usefulness_pmv_scaled, uc_gated, uc_hybrid (and
  # their _gebv/_adj/_blend variants) and dh_recomb_rel_var, dh_pmv_scaled_var,
  # gated_var, hybrid_var were removed because they mixed genetic-value
  # variance units with the parent_distance relatedness distance.
  out$usefulness_vpm       <- out$cross_mean                   + i * sqrt(out$vpm)
  out$usefulness_pmv           <- out$cross_mean                   + i * sqrt(out$pmv)
  out$usefulness_vpm_gebv  <- out$cross_mean_gebv              + i * sqrt(out$vpm)
  out$usefulness_pmv_gebv      <- out$cross_mean_gebv              + i * sqrt(out$pmv)
  out$usefulness_vpm_adj   <- out$cross_mean_adj    + i * sqrt(out$vpm)
  out$usefulness_pmv_adj       <- out$cross_mean_adj    + i * sqrt(out$pmv)
  out$usefulness_vpm_blend <- out$cross_mean_blend             + i * sqrt(out$vpm)
  out$usefulness_pmv_blend     <- out$cross_mean_blend             + i * sqrt(out$pmv)
  out$rank_score <- if (effect_rel >= min_effect_reliability) out$usefulness_pmv_gebv else out$usefulness_pmv_blend
  attr(out, "progeny_target") <- target
  attr(out, "recomb_model") <- recomb_model
  out
}

ng_pair_relationship_variance <- function(pairs, K) {
  p1 <- match(pairs$parent1, rownames(K))
  p2 <- match(pairs$parent2, rownames(K))
  v <- 0.5 * (diag(K)[p1] + diag(K)[p2]) - K[cbind(p1, p2)]
  data.frame(
    parent_distance = pmax(as.numeric(v), 0),
    pair_kinship = as.numeric(K[cbind(p1, p2)])
  )
}

# Distribution of expected progeny inbreeding across a set of crosses (a
# progeny-inbreeding histogram). `x` may be a scored candidate table or a
# selected mating plan (uses expected_progeny_inbreeding, else falls back to
# pair_kinship / 2), or a plain numeric vector of progeny-inbreeding values. Returns a
# per-bin data.frame (bin bounds, count, proportion) with mean/max/n as attributes;
# ng_write_progeny_inbreeding_histogram_json() serializes it for the frontend.
ng_progeny_inbreeding_histogram <- function(x, breaks = 20L, value_range = NULL) {
  f <- if (is.data.frame(x)) {
    if ("expected_progeny_inbreeding" %in% names(x)) {
      as.numeric(x$expected_progeny_inbreeding)
    } else if ("pair_kinship" %in% names(x)) {
      as.numeric(x$pair_kinship) / 2
    } else {
      ng_stop("x must have an expected_progeny_inbreeding or pair_kinship column")
    }
  } else {
    as.numeric(x)
  }
  f <- f[is.finite(f)]
  if (!length(f)) ng_stop("no finite progeny-inbreeding values to bin")
  rng <- if (is.null(value_range)) range(f) else suppressWarnings(as.numeric(value_range))
  if (length(rng) < 2L || !all(is.finite(rng)) || diff(range(rng)) <= 0) {
    rng <- f[[1L]] + c(-0.5, 0.5)
  }
  rng <- range(rng)
  brk <- seq(rng[[1L]], rng[[2L]], length.out = as.integer(breaks) + 1L)
  clamped <- pmin(pmax(f, rng[[1L]]), rng[[2L]])
  h <- graphics::hist(clamped, breaks = brk, plot = FALSE)
  nb <- length(h$breaks)
  out <- data.frame(
    bin_lower = h$breaks[-nb],
    bin_upper = h$breaks[-1L],
    bin_mid = h$mids,
    count = as.integer(h$counts),
    proportion = h$counts / sum(h$counts)
  )
  attr(out, "mean_progeny_inbreeding") <- mean(f)
  attr(out, "max_progeny_inbreeding") <- max(f)
  attr(out, "n") <- length(f)
  out
}

ng_dh_recomb_variance_pairs <- function(geno,
                                        beta,
                                        beta_var,
                                        marker_map,
                                        ids = rownames(geno),
                                        pairs,
                                        window_cm = Inf,
                                        use_cpp = TRUE,
                                        recomb_model = NULL,
                                        target = c("DH", "RIL")) {
  target <- match.arg(target)
  if (is.null(recomb_model)) {
    recomb_model <- if ("recomb_model" %in% names(marker_map)) {
      unique(as.character(marker_map$recomb_model))[[1L]]
    } else {
      "haldane"
    }
  }
  recomb_model <- match.arg(recomb_model, c("haldane", "kosambi"))
  # The closed-form chromosome recursion (R and C++ kernels) relies on the
  # multiplicative Markov property of the Haldane DH decay exp(-2 d / 100):
  # decay(a + b) = decay(a) * decay(b). This property fails for:
  #   - Kosambi DH:  1 - tanh(2 d / 100)
  #   - Haldane RIL: exp(-2 d / 100) / (2 - exp(-2 d / 100))
  # Both go through the dense O(M^2) path.
  #
  # The closed-form recursion also cannot apply a finite `window_cm` cutoff: it telescopes the
  # whole chromosome, so a distance window would break the multiplicative property. The recursion
  # is therefore valid ONLY for the full-chromosome case (window_cm = Inf); a finite window must
  # go through the dense/banded path, which zeroes off-window R entries and honors the cutoff.
  use_dense <- !identical(recomb_model, "haldane") || !identical(target, "DH") ||
    is.finite(window_cm)
  if (use_dense) {
    # Optimisation: when the user passes a finite window_cm small relative to
    # the average chromosome length, switch from the O(m^2) dense kernel to the
    # O(m * k_window) banded kernel. The dense kernel still builds the full
    # m x m R matrix even when window_cm zeros most of it; the banded kernel
    # never materialises off-window entries.
    if (ng_banded_kernel_preferred(marker_map, window_cm)) {
      return(ng_dh_recomb_variance_pairs_banded(
        geno = geno, beta = beta, beta_var = beta_var,
        marker_map = marker_map, ids = ids, pairs = pairs,
        recomb_model = recomb_model, target = target, window_cm = window_cm
      ))
    }
    return(ng_dh_recomb_variance_pairs_dense(
      geno = geno, beta = beta, beta_var = beta_var,
      marker_map = marker_map, ids = ids, pairs = pairs,
      recomb_model = recomb_model, target = target, window_cm = window_cm
    ))
  }
  if (isTRUE(use_cpp) && exists("ng_dh_recomb_pairs_cpp", mode = "function", inherits = TRUE)) {
    return(ng_dh_recomb_pairs_cpp(
      geno = geno,
      beta = as.numeric(beta),
      beta_var = as.numeric(beta_var),
      chr = as.integer(marker_map$chr_index),
      pos_cm = as.numeric(marker_map$pos_cm),
      ids = as.character(ids),
      pair_parent1 = as.character(pairs$parent1),
      pair_parent2 = as.character(pairs$parent2),
      window_cm = as.numeric(window_cm)
    ))
  }
  ng_dh_recomb_variance_pairs_r(geno, beta, beta_var, marker_map, ids, pairs, window_cm)
}

ng_dh_recomb_variance_pairs_full_posterior <- function(geno,
                                                       beta,
                                                       beta_cov_full,
                                                       marker_map,
                                                       ids,
                                                       pairs,
                                                       target = c("DH", "RIL"),
                                                       recomb_model = c("haldane", "kosambi"),
                                                       window_cm = Inf,
                                                       use_cpp = TRUE) {
  # Full off-diagonal posterior PMV (genomicMateSelectR formulation):
  #   PMV = a' R a + d' (R o Sigma_beta) d
  # where a_k = d_k * beta_hat_k, d_k = 0.5 * (geno[p1, k] - geno[p2, k]) and
  # R[k, l] = 1 - 2 R_kl is the progeny-target recombination decay matrix from
  # ng_recomb_decay_matrix(). Falls back to the diagonal-only correction when
  # beta_cov_full is NULL or all-zero. The 1 - 2 R term in R[k, l] doubles as
  # the kernel for both VPM (a'Ra) and the off-diagonal PMV correction
  # (d'(R o Sigma) d).
  recomb_model <- match.arg(recomb_model)
  target <- match.arg(target)
  m <- length(beta)
  marker_ids <- marker_map$marker
  if (is.null(rownames(beta_cov_full)) || is.null(colnames(beta_cov_full))) {
    if (nrow(beta_cov_full) != m || ncol(beta_cov_full) != m) {
      ng_stop("beta_cov_full must be m x m matching marker_map$marker when dimnames are missing")
    }
    Sigma <- beta_cov_full
  } else {
    Sigma <- beta_cov_full[marker_ids, marker_ids, drop = FALSE]
  }
  if (m > 12000L) {
    warning(sprintf(
      "Full-posterior recombination-variance path on %d markers (~%.1f GB); consider LD pruning.",
      m, (m * m * 8) / (1024 ^ 3)
    ), call. = FALSE)
  }
  R <- ng_recomb_decay_matrix(marker_map, model = recomb_model, target = target)
  if (is.finite(window_cm) && window_cm > 0) {
    dist <- outer(marker_map$pos_cm, marker_map$pos_cm, function(a, b) abs(a - b))
    R[dist > window_cm] <- 0
  }
  # Diagonal-only correction (legacy column) uses diag(Sigma_beta) only.
  beta_var_diag <- as.numeric(pmax(diag(Sigma), 0))
  # Precompute the Hadamard product R o Sigma_beta once so the per-pair PMV
  # extra term reduces to a single d' M d quadratic form.
  R_had_Sigma <- R * Sigma
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  beta <- as.numeric(beta)
  # Per-pair quadratic forms a'Ra and d'(R o Sigma) d are O(m^2). For m above
  # ~1000 the R loop becomes the dominant per-grid cost; route to the C++
  # kernel (column-major, cache-friendly outer-over-k inner-over-i) when
  # available. Numerical equivalence to the R reference is at machine
  # precision (~1e-12); verified by tests/posterior_pmv_full.R.
  if (isTRUE(use_cpp) && exists("ng_dh_recomb_pairs_full_posterior_cpp", mode = "function", inherits = TRUE)) {
    return(ng_dh_recomb_pairs_full_posterior_cpp(
      geno = geno, beta = beta, beta_var_diag = beta_var_diag,
      R_decay = R, R_had_Sigma = R_had_Sigma,
      pair_p1_zero = as.integer(p1 - 1L),
      pair_p2_zero = as.integer(p2 - 1L)
    ))
  }
  out_v <- numeric(nrow(pairs))
  out_pmv_diag <- numeric(nrow(pairs))
  out_pmv_full <- numeric(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    d <- 0.5 * (geno[p1[r], ] - geno[p2[r], ])
    d[!is.finite(d)] <- 0
    a <- d * beta
    a[!is.finite(a)] <- 0
    Ra <- as.numeric(R %*% a)
    v <- sum(a * Ra)
    pmv_diag_extra <- sum(d * d * beta_var_diag)
    pmv_off_extra <- as.numeric(crossprod(d, R_had_Sigma %*% d))
    out_v[r]        <- max(0, v)
    out_pmv_diag[r] <- max(0, v + pmv_diag_extra)
    out_pmv_full[r] <- max(0, v + pmv_off_extra)
  }
  data.frame(
    vpm = out_v,
    pmv = out_pmv_diag,
    pmv_full_posterior = out_pmv_full
  )
}

# Heuristic: when the user-supplied window_cm is finite and small relative to
# the average chromosome length, the banded kernel materialises only the in-
# window upper-triangle entries instead of the full m x m R matrix and runs in
# O(m * k_window) per pair vs. O(m^2) per pair for the dense kernel.
#
# Threshold: window_cm / mean_chr_len_cm < 0.5 means the window covers less
# than half of an average chromosome, so the band is non-trivially narrower
# than the dense matrix. The 0.5 ratio is empirical: under it the banded path
# starts dominating once m runs into the thousands; above it the dense kernel
# wins because of vectorised matrix multiply throughput. Override via the
# NGCD_BANDED_RATIO environment variable for benchmarking.
ng_banded_kernel_preferred <- function(marker_map, window_cm,
                                       ratio_threshold = NULL) {
  if (!is.finite(window_cm) || window_cm <= 0) return(FALSE)
  if (is.null(ratio_threshold)) {
    env <- suppressWarnings(as.numeric(Sys.getenv("NGCD_BANDED_RATIO", unset = "")))
    ratio_threshold <- if (is.finite(env) && env > 0) env else 0.5
  }
  pos <- as.numeric(marker_map$pos_cm)
  chr <- as.character(marker_map$chr)
  if (!length(pos)) return(FALSE)
  chr_lengths <- vapply(split(pos, chr), function(p) {
    if (length(p) < 2L) 0 else max(p) - min(p)
  }, numeric(1L))
  chr_lengths <- chr_lengths[chr_lengths > 0]
  if (!length(chr_lengths)) return(FALSE)
  mean_chr_len <- mean(chr_lengths)
  if (!is.finite(mean_chr_len) || mean_chr_len <= 0) return(FALSE)
  (window_cm / mean_chr_len) < ratio_threshold
}

ng_dh_recomb_variance_pairs_banded <- function(geno,
                                               beta,
                                               beta_var,
                                               marker_map,
                                               ids,
                                               pairs,
                                               target = c("DH", "RIL"),
                                               recomb_model = c("haldane", "kosambi"),
                                               window_cm = Inf,
                                               max_window_warn_marker_count = 12000L,
                                               use_cpp = TRUE) {
  # Banded sparse-COO consumer of the Kosambi / RIL decay kernel. Equivalent to
  # the dense path: VPM = a' R a, PMV = VPM + diag-only marker-effect variance
  # term sum(d^2 * beta_var). The (i, j, x) triplets cover the upper triangle
  # plus the diagonal (x = 1); the off-diagonal contribution is scaled by 2 to
  # account for symmetry, matching a' R a exactly when R is symmetric.
  recomb_model <- match.arg(recomb_model)
  target <- match.arg(target)
  m <- length(beta)
  if (m > max_window_warn_marker_count && !is.finite(window_cm)) {
    warning(sprintf(
      "Banded recombination-variance path called with window_cm = Inf on %d markers (~%.1f GB if it falls back to dense); supply a finite window_cm to keep the kernel sparse.",
      m, (m * m * 8) / (1024 ^ 3)
    ), call. = FALSE)
  }
  band <- ng_recomb_decay_banded(marker_map = marker_map, model = recomb_model,
                                 target = target, window_cm = window_cm)
  i_idx <- band$i
  j_idx <- band$j
  x_val <- band$x
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  beta <- as.numeric(beta)
  beta_var <- as.numeric(beta_var)
  # Route the per-pair COO sum through the C++ kernel when available; it
  # walks the triplets in C with no R-loop overhead and applies the
  # symmetry doubling inline. Numerically equivalent to the R reference to
  # ~1e-12.
  if (isTRUE(use_cpp) && exists("ng_dh_recomb_pairs_banded_cpp", mode = "function", inherits = TRUE)) {
    return(ng_dh_recomb_pairs_banded_cpp(
      geno = geno, beta = beta, beta_var = beta_var,
      band_i_zero = as.integer(i_idx - 1L),
      band_j_zero = as.integer(j_idx - 1L),
      band_x = as.numeric(x_val),
      pair_p1_zero = as.integer(p1 - 1L),
      pair_p2_zero = as.integer(p2 - 1L)
    ))
  }
  is_diag <- (i_idx == j_idx)
  # Scale factor: diagonal entries contribute once (a_i^2 * x), off-diagonals
  # contribute twice because the banded representation stores only i < j.
  weight <- ifelse(is_diag, 1, 2)
  out_v <- numeric(nrow(pairs))
  out_pmv <- numeric(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    d <- 0.5 * (geno[p1[r], ] - geno[p2[r], ])
    d[!is.finite(d)] <- 0
    a <- d * beta
    a[!is.finite(a)] <- 0
    v <- sum(a[i_idx] * a[j_idx] * x_val * weight)
    pmv_extra <- sum(d * d * pmax(beta_var, 0), na.rm = TRUE)
    out_v[r]   <- max(0, v)
    out_pmv[r] <- max(0, v + pmv_extra)
  }
  data.frame(vpm = out_v, pmv = out_pmv)
}

ng_dh_recomb_variance_pairs_dense <- function(geno,
                                              beta,
                                              beta_var,
                                              marker_map,
                                              ids,
                                              pairs,
                                              recomb_model = c("haldane", "kosambi"),
                                              target = c("DH", "RIL"),
                                              window_cm = Inf) {
  recomb_model <- match.arg(recomb_model)
  target <- match.arg(target)
  m <- length(beta)
  if (m > 12000L) {
    warning(sprintf(
      "Dense recombination-variance path on %d markers (~%.1f GB); consider LD pruning or chromosome-wise scoring.",
      m, (m * m * 8) / (1024 ^ 3)
    ), call. = FALSE)
  }
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  beta <- as.numeric(beta)
  beta_var <- as.numeric(beta_var)
  # Batch all pairs into BLAS GEMMs instead of one R-level matvec per pair, and exploit
  # the block-diagonal structure of the recombination-decay matrix: cross-chromosome
  # correlation is exactly 0, so v_r = a_r' R a_r = sum over chromosome blocks of
  # a_{r,c}' R_c a_{r,c}. This is EXACT (identical to the full dense a'Ra) while
  # avoiding the full m x m matrix and cutting GEMM flops ~n_chromosomes-fold.
  # a_r = d_r * beta -> column r of A; v = sum_c colSums(A_c * (R_c %*% A_c)).
  chr <- marker_map$chr_index; if (is.null(chr)) chr <- marker_map$chr
  D <- 0.5 * (geno[p1, , drop = FALSE] - geno[p2, , drop = FALSE])  # n_pairs x m
  D[!is.finite(D)] <- 0
  A <- t(D) * beta                                                  # m x n_pairs (beta recycled down cols)
  A[!is.finite(A)] <- 0
  v <- numeric(nrow(pairs))
  for (cc in unique(chr)) {
    idx <- which(chr == cc)
    Rc <- ng_recomb_decay_matrix(marker_map[idx, , drop = FALSE], model = recomb_model, target = target)
    if (is.finite(window_cm) && window_cm > 0) {
      dc <- outer(marker_map$pos_cm[idx], marker_map$pos_cm[idx], function(a, b) abs(a - b))
      Rc[dc > window_cm] <- 0
    }
    Ac <- A[idx, , drop = FALSE]
    v <- v + colSums(Ac * (Rc %*% Ac))
  }
  bv <- pmax(beta_var, 0); bv[!is.finite(bv)] <- 0
  pmv_extra <- as.numeric((D * D) %*% bv)
  data.frame(vpm = pmax(0, v), pmv = pmax(0, v + pmv_extra))
}

ng_dh_recomb_variance_pairs_r <- function(geno, beta, beta_var, marker_map, ids, pairs, window_cm = Inf) {
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  chr <- marker_map$chr_index
  pos <- marker_map$pos_cm
  out_v <- numeric(nrow(pairs))
  out_pmv <- numeric(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    d <- 0.5 * (geno[p1[r], ] - geno[p2[r], ])
    a <- d * beta
    a2 <- d * d * (beta * beta + beta_var)
    v <- 0
    pmv <- 0
    carry <- 0
    prev_chr <- NA_integer_
    prev_pos <- NA_real_
    for (k in seq_along(a)) {
      if (!identical(chr[k], prev_chr)) {
        carry <- 0
      } else {
        carry <- carry * ng_haldane_decay(pos[k] - prev_pos)
      }
      if (is.finite(a[k]) && is.finite(a2[k])) {
        v <- v + a[k] * a[k] + 2 * a[k] * carry
        pmv <- pmv + a2[k] + 2 * a[k] * carry
        carry <- carry + a[k]
      }
      prev_chr <- chr[k]
      prev_pos <- pos[k]
    }
    out_v[r] <- max(0, v)
    out_pmv[r] <- max(0, pmv)
  }
  data.frame(vpm = out_v, pmv = out_pmv)
}

ng_exact_gms_additive_var <- function(parent1,
                                      parent2,
                                      haplo_mat,
                                      beta,
                                      recomb_decay_mat,
                                      beta_cov = NULL) {
  sire_ld <- ng_gametic_ld_parent(parent1, haplo_mat, recomb_decay_mat)
  dam_ld <- ng_gametic_ld_parent(parent2, haplo_mat, recomb_decay_mat)
  cross_ld <- sire_ld + dam_ld
  beta <- beta[colnames(cross_ld)]
  vpm <- as.numeric(crossprod(beta, cross_ld %*% beta))
  if (!is.null(beta_cov)) {
    beta_cov <- beta_cov[colnames(cross_ld), colnames(cross_ld), drop = FALSE]
    pmv <- vpm + sum(diag(cross_ld %*% beta_cov))
  } else {
    pmv <- vpm
  }
  c(VPM = vpm, PMV = pmv)
}

ng_gametic_ld_parent <- function(parent, haplo_mat, recomb_decay_mat) {
  rows <- paste0(parent, c("_HapA", "_HapB"))
  X <- haplo_mat[rows, colnames(recomb_decay_mat), drop = FALSE]
  p <- colMeans(X)
  recomb_decay_mat * ((0.5 * crossprod(X)) - tcrossprod(p))
}
