#' Resolve the parent line-type policy governing the heterozygosity audit
#'
#' `parent_type` is the canonical control:
#'   * "inbred" / "dh" -- fully fixed lines. A DH line is homozygous by
#'     construction and finished inbreds are effectively so; a heterozygous
#'     locus beyond the QC tolerance is a genotyping/data error and is a
#'     BLOCKER (do not proceed).
#'   * "ril" -- recombinant inbred lines may retain residual heterozygosity
#'     after finite selfing; scoring then requires complete phased parental
#'     haplotypes so the package does not silently use a biased inbred kernel.
#' The legacy boolean `assume_inbred` is deprecated: if supplied it overrides
#' `parent_type` (TRUE -> "inbred", FALSE -> "ril") with a one-time warning.
#' Returns a single canonical parent_type string.
ng_reconcile_parent_type <- function(parent_type = c("inbred", "dh", "ril"),
                                     assume_inbred = NULL) {
  parent_type <- match.arg(parent_type)
  if (!is.null(assume_inbred)) {
    mapped <- if (isTRUE(assume_inbred)) "inbred" else "ril"
    warning(sprintf(
      "`assume_inbred` is deprecated; use parent_type = 'inbred' / 'dh' / 'ril'. Mapping assume_inbred = %s to parent_type = '%s'.",
      as.character(assume_inbred), mapped), call. = FALSE)
    return(mapped)
  }
  parent_type
}

# TRUE when the parent line type is a fully fixed line (DH or finished inbred),
# for which residual heterozygosity is a data error rather than biology.
ng_parent_type_blocks_het <- function(parent_type) parent_type %in% c("inbred", "dh")

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
                             parent_type = c("inbred", "dh", "ril"),
                             selection_prop = 0.10,
                             min_effect_reliability = 0.35,
                             recomb_model = c("haldane", "kosambi"),
                             window_cm = Inf,
                             use_cpp = TRUE,
                             assume_inbred = NULL,
                             ploidy = 2L,
                             inbred_tolerance = 0.05,
                             inbred_marker_fraction = 0.02,
                             dh_marker_fraction = 0.005,
                             phased_haplotypes = NULL,
                             posterior_cov_full = NULL,
                             parent_kinship = NULL,
                             grm_method = c("vanraden", "yang")) {
  target <- match.arg(target)
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  block_het <- ng_parent_type_blocks_het(parent_type)
  recomb_model <- match.arg(recomb_model)
  grm_method <- match.arg(grm_method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  if (is.null(effects$beta) || is.null(names(effects$beta)) || anyDuplicated(names(effects$beta))) {
    ng_stop("effects$beta must be a uniquely named marker-effect vector")
  }
  missing_beta <- setdiff(colnames(geno), names(effects$beta))
  if (length(missing_beta)) {
    ng_stop("effects$beta is missing markers used by geno: ",
            paste(utils::head(missing_beta, 6L), collapse = ", "))
  }
  beta_input <- suppressWarnings(as.numeric(effects$beta[colnames(geno)]))
  if (any(!is.finite(beta_input))) ng_stop("effects$beta must be finite for every genotype marker")
  names(beta_input) <- colnames(geno)
  intercept_input <- suppressWarnings(as.numeric(effects$intercept))
  if (length(intercept_input) != 1L || !is.finite(intercept_input)) {
    ng_stop("effects$intercept must be one finite number")
  }
  if (is.null(effects$beta_var)) {
    beta_var_input <- stats::setNames(rep(0, ncol(geno)), colnames(geno))
  } else {
    if (is.null(names(effects$beta_var)) || anyDuplicated(names(effects$beta_var))) {
      ng_stop("effects$beta_var must be NULL or a uniquely named marker-variance vector")
    }
    missing_beta_var <- setdiff(colnames(geno), names(effects$beta_var))
    if (length(missing_beta_var)) {
      ng_stop("effects$beta_var is missing markers used by geno: ",
              paste(utils::head(missing_beta_var, 6L), collapse = ", "))
    }
    beta_var_input <- suppressWarnings(as.numeric(effects$beta_var[colnames(geno)]))
    if (any(!is.finite(beta_var_input)) || any(beta_var_input < 0)) {
      ng_stop("effects$beta_var must contain finite non-negative marginal variances")
    }
    names(beta_var_input) <- colnames(geno)
  }
  # The DH/RIL recombination-variance formula assumes both parents are fully
  # inbred (marker dosage in {0, 2}). When parents carry residual heterozygosity
  # the formula Cov(y_k, y_l) = d_k d_l (1 - 2 r_kl) no longer holds because the
  # F1 may be homozygous at some heterozygous loci. Refuse to silently produce
  # wrong PMV: detect non-inbred dosages and require the information needed by
  # the exact phased-parent correction.
  # A doubled haploid is 100% homozygous by construction, so DH material gets a
  # STRICT het-marker fraction floor (`dh_marker_fraction`, default 0.5%): true
  # het is a data error (wrong ploidy, contamination, a RIL mislabelled DH), and
  # the small floor only absorbs routine genotyping noise (~0.1-0.5% per-call
  # het error on real SNP/GBS panels) so clean DH data is not false-blocked
  # while contamination (>> the floor) still blocks. Finished inbred lines
  # ('inbred') keep the looser `inbred_marker_fraction` (default 2%); RILs are
  # not audited for a blocker. The per-dosage `inbred_tolerance` (numeric
  # rounding, e.g. 1.998 -> 2) applies in all cases. ploidy is threaded so the
  # audit does not misread a homozygous polyploid dosage as heterozygous.
  het_marker_fraction <- if (identical(parent_type, "dh")) dh_marker_fraction else inbred_marker_fraction
  inbred_audit <- ng_audit_inbred_dosage(
    geno, ploidy = ploidy, tolerance = inbred_tolerance,
    fraction_tolerance = het_marker_fraction
  )
  if (block_het && length(inbred_audit$violators)) {
    msg <- sprintf(
      "%d / %d parents carry heterozygous loci beyond tolerance (max het-marker frac = %.3f, examples: %s), but parent_type = '%s' declares fully fixed lines. A doubled-haploid (DH) line is homozygous by construction, so heterozygous loci in DH / fixed material indicate a genotyping or data error -- BLOCKED, do not proceed. If these are RILs (which legitimately retain residual heterozygosity at a few loci after finite selfing), set parent_type = 'ril' to proceed. With phased_haplotypes supplied, het-parent crosses then get the exact residual-het variance (ng_gms_additive_var_general); without phased haplotypes the a'Ra kernel treats parents as inbred and is biased low at the het loci.",
      length(inbred_audit$violators), nrow(geno), inbred_audit$max_fraction,
      paste(head(inbred_audit$violators, 4L), collapse = ", "), parent_type
    )
    ng_stop(msg)
  }
  has_residual_het <- !block_het && any(inbred_audit$fraction > 0)
  if (has_residual_het && is.null(phased_haplotypes)) {
    ng_stop("parent_type = 'ril' includes residual-heterozygous parents, but phased_haplotypes ",
            "was not supplied. The inbred a'Ra kernel is biased low for these parents; ",
            "supply complete phased haplotypes or use fully inbred parents.")
  }
  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = include_self)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("pairs must contain parent1 and parent2")
  }
  pairs$parent1 <- as.character(pairs$parent1)
  pairs$parent2 <- as.character(pairs$parent2)
  if (anyNA(pairs$parent1) || anyNA(pairs$parent2) ||
      any(!nzchar(trimws(pairs$parent1))) || any(!nzchar(trimws(pairs$parent2)))) {
    ng_stop("pairs must contain non-missing, non-blank parent IDs")
  }
  unknown_parents <- setdiff(unique(c(pairs$parent1, pairs$parent2)), ids)
  if (length(unknown_parents)) {
    ng_stop("pairs reference parents absent from ids: ", paste(unknown_parents, collapse = ", "))
  }
  if (anyDuplicated(ng_group_pair_key(pairs$parent1, pairs$parent2))) {
    ng_stop("pairs contains duplicate unordered parent pairs")
  }

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
    if (!is.numeric(pk) || nrow(pk) != ncol(pk) || is.null(rownames(pk)) ||
        is.null(colnames(pk)) || anyDuplicated(rownames(pk)) || anyDuplicated(colnames(pk)) ||
        !setequal(rownames(pk), colnames(pk))) {
      ng_stop("parent_kinship must be a uniquely named square additive-relationship matrix")
    }
    missing_k <- setdiff(ids, intersect(rownames(pk), colnames(pk)))
    if (length(missing_k)) {
      ng_stop("parent_kinship is missing scored parents: ", paste(missing_k, collapse = ", "))
    }
    pk <- pk[ids, ids, drop = FALSE]
    if (any(!is.finite(pk)) || max(abs(pk - t(pk))) > 1e-8) {
      ng_stop("parent_kinship must contain a finite symmetric additive-relationship matrix")
    }
    pk <- (pk + t(pk)) / 2
    ev <- eigen(pk, symmetric = TRUE, only.values = TRUE)$values
    tol <- 1e-8 * max(1, max(abs(ev)))
    if (min(ev) < -tol) ng_stop("parent_kinship must be positive semidefinite")
    pk
  }
  rel_var <- ng_pair_relationship_variance(pairs, K)

  beta <- beta_input
  beta_var <- beta_var_input
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
  dh_pmv_full <- if (!is.null(dh$pmv_full_posterior)) pmax(dh$pmv_full_posterior, 0) else NULL
  # Exact residual-het-parent correction: when the parents are declared RILs and
  # phased haplotypes are supplied, replace the (biased-low) inbred a'Ra vpm/pmv
  # of het-parent crosses with the exact phased-haplotype variance
  # (ng_gms_additive_var_general). Inbred-parent crosses and the no-phase / DH /
  # inbred paths are untouched.
  if (identical(parent_type, "ril") && !is.null(phased_haplotypes)) {
    corr <- ng_apply_het_parent_correction(
      dh_var, dh_pmv, dh_pmv_full, phased_haplotypes, sorted, pairs,
      target = target, recomb_model = recomb_model,
      beta_var = sorted$beta_var, beta_cov_full = posterior_cov_full)
    dh_var <- corr$vpm; dh_pmv <- corr$pmv; dh_pmv_full <- corr$pmv_full
  }
  effect_rel <- suppressWarnings(as.numeric(mean_source$reliability[[1L]]))
  if (!is.finite(effect_rel)) effect_rel <- NA_real_
  # Compatibility name only: this is now the midpoint of the explicitly chosen
  # mean source, not an uncalibrated phenotype/GEBV interpolation.
  blended_mean <- parent_mean

  out <- pairs
  out$mean_source <- mean_source$source
  out$effect_reliability <- effect_rel
  out$effect_reliability_is_calibrated <- isTRUE(mean_source$reliability_is_calibrated)
  out$effect_cv_predictive_r2 <- suppressWarnings(as.numeric((mean_source$cv_predictive_r2 %||% NA_real_)[[1L]]))
  out$progeny_target <- target
  out$mid_parent_value <- as.numeric(mpv_gebv)
  out$cross_mean <- as.numeric(parent_mean)
  out$cross_mean_gebv <- as.numeric(mpv_gebv)
  out$cross_mean_adj <- as.numeric(adjusted_mean)
  out$cross_mean_blend <- as.numeric(blended_mean)
  # Relationship-distance metric (NOT a variance). Reported separately so the
  # OCS optimizer can use it as a diversity penalty; never combined with PMV.
  out$parent_distance <- rel_var$parent_distance
  out$pair_relationship <- rel_var$pair_relationship
  out$pair_kinship <- rel_var$pair_kinship
  # Expected inbreeding of the immediate progeny of each cross = coancestry of the two
  # parents (Meuwissen OCS: F_progeny = f(p1, p2)). pair_relationship is read from
  # the VanRaden genomic relationship matrix; pair_kinship is relationship / 2. This is the pedigree/
  # relationship-scale progeny inbreeding managed in tactical mate selection; it is distinct
  # from the eventual homozygosity of a finished DH/RIL line (which tends to 1
  # regardless of the parents). Reported so the optimizer and reports can treat progeny
  # inbreeding as a first-class quantity alongside parental (group) coancestry.
  out$expected_progeny_inbreeding <- pmax(0, as.numeric(rel_var$pair_kinship))
  # Within-family genetic-value variances under the requested progeny target
  # (DH or RIL), in (genetic value)^2 units. pmv adds the diagonal
  # posterior marker-effect uncertainty to vpm (VPM). When the
  # caller supplies posterior_cov_full, pmv_full_posterior also
  # carries the full off-diagonal correction d'(R o Sigma_beta) d
  # (genomicMateSelectR formulation); NA otherwise.
  out$vpm <- dh_var
  out$pmv <- dh_pmv
  if (!is.null(dh_pmv_full)) {
    out$pmv_full_posterior <- dh_pmv_full
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
  out$rank_score <- out$usefulness_pmv_blend
  attr(out, "progeny_target") <- target
  attr(out, "recomb_model") <- recomb_model
  out
}

ng_pair_relationship_variance <- function(pairs, K) {
  p1 <- match(pairs$parent1, rownames(K))
  p2 <- match(pairs$parent2, rownames(K))
  v <- 0.5 * (diag(K)[p1] + diag(K)[p2]) - K[cbind(p1, p2)]
  pair_relationship <- as.numeric(K[cbind(p1, p2)])
  data.frame(
    parent_distance = pmax(as.numeric(v), 0),
    pair_relationship = pair_relationship,
    pair_kinship = pair_relationship / 2
  )
}

# Distribution of expected progeny inbreeding across a set of crosses (a
# progeny-inbreeding histogram). `x` may be a scored candidate table or a
# selected mating plan (uses expected_progeny_inbreeding, else falls back to
# pair_kinship), or a plain numeric vector of progeny-inbreeding values. Returns a
# per-bin data.frame (bin bounds, count, proportion) with mean/max/n as attributes;
# ng_write_progeny_inbreeding_histogram_json() serializes it for the frontend.
ng_progeny_inbreeding_histogram <- function(x, breaks = 20L, value_range = NULL) {
  f <- if (is.data.frame(x)) {
    if ("expected_progeny_inbreeding" %in% names(x)) {
      as.numeric(x$expected_progeny_inbreeding)
    } else if ("pair_kinship" %in% names(x)) {
      as.numeric(x$pair_kinship)
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
  # This is an internal, performance-oriented kernel, but malformed dimensions
  # must still fail in R rather than reaching native code. In particular, an
  # absent chr_index previously became integer(0) at the .Call boundary and the
  # C++ loop indexed beyond it. Public scoring prepares/sorts the map before this
  # call; direct validation/research callers must do the same.
  geno <- ng_as_numeric_matrix(geno, "geno")
  marker_map <- as.data.frame(marker_map, stringsAsFactors = FALSE)
  m <- ncol(geno)
  if (length(beta) != m || length(beta_var) != m) {
    ng_stop("beta and beta_var must each have one value per genotype marker")
  }
  if (any(!is.finite(as.numeric(beta))) ||
      any(!is.finite(as.numeric(beta_var))) || any(as.numeric(beta_var) < 0)) {
    ng_stop("beta must be finite and beta_var must be finite and non-negative")
  }
  required_map <- c("marker", "chr", "chr_index", "pos_cm")
  missing_map <- setdiff(required_map, names(marker_map))
  if (length(missing_map)) {
    ng_stop("marker_map must be prepared with ng_prepare_marker_map and contain: ",
            paste(required_map, collapse = ", "), "; missing: ",
            paste(missing_map, collapse = ", "))
  }
  if (nrow(marker_map) != m ||
      length(marker_map$chr_index) != m || length(marker_map$pos_cm) != m) {
    ng_stop("marker_map must have exactly one prepared row per genotype marker")
  }
  if (is.null(colnames(geno)) ||
      !identical(as.character(marker_map$marker), as.character(colnames(geno)))) {
    ng_stop("prepared marker_map rows must match genotype markers in the same order")
  }
  if (anyNA(marker_map$chr_index) || any(!is.finite(as.numeric(marker_map$pos_cm)))) {
    ng_stop("prepared marker_map chr_index and pos_cm must be complete and finite")
  }
  map_order <- order(marker_map$chr_index, marker_map$pos_cm, marker_map$marker)
  if (!identical(as.integer(map_order), seq_len(m))) {
    ng_stop("marker_map and genotype columns must be sorted by chromosome and genetic position")
  }
  ids <- as.character(ids)
  if (length(ids) != nrow(geno) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    ng_stop("ids must be unique, non-missing, and have one value per genotype row")
  }
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("pairs must contain parent1 and parent2")
  }
  pairs$parent1 <- as.character(pairs$parent1)
  pairs$parent2 <- as.character(pairs$parent2)
  if (anyNA(pairs$parent1) || anyNA(pairs$parent2) ||
      length(setdiff(unique(c(pairs$parent1, pairs$parent2)), ids))) {
    ng_stop("pairs must reference non-missing parent IDs present in ids")
  }
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

#' Exact within-cross additive variance for arbitrary phased parents
#'
#' Generalizes the inbred-parent kernel `a'Ra` to parents that carry residual
#' heterozygosity (real RILs), from the phased parental haplotypes. See
#' docs/design/residual-het-parent-variance.md for the derivation. Writing
#' `d1 = HapA - HapB`, `d2 = HapA - HapB` for the two parents (nonzero only at
#' heterozygous loci) and `R` the target recombination kernel:
#'   DH  : Var = a'Ra + 1/2 b'[ (1-r) o R o (d1 d1' + d2 d2') ] b, (1-r) = (1+R)/2
#'   RIL(inf): Var = a'R*a + 1/2 b'[ R* o (d1 d1' + d2 d2') ] b,  R* = (1-2r)/(1+2r)
#' with `a_k = 1/2 (x1_k - x2_k) b_k`. For inbred parents (d1 = d2 = 0) this
#' reduces exactly to `a'Ra`. PMV uses the same full kernel K (Var = b'Kb):
#' `PMV = b'Kb + trace(K Sigma_b)`. `recomb_decay_mat` must be the
#' target-appropriate kernel (DH: 1-2r; RIL: (1-2r)/(1+2r)).
ng_gms_additive_var_general <- function(parent1, parent2, haplo_mat, beta,
                                        recomb_decay_mat, beta_cov = NULL,
                                        target = c("DH", "RIL")) {
  target <- match.arg(target)
  markers <- colnames(recomb_decay_mat)
  get_pair <- function(p) {
    rows <- paste0(p, c("_HapA", "_HapB"))
    miss <- setdiff(rows, rownames(haplo_mat))
    if (length(miss)) ng_stop("haplo_mat is missing haplotypes: ", paste(miss, collapse = ", "))
    matrix(haplo_mat[rows, markers, drop = FALSE], nrow = 2L,
           dimnames = list(c("HapA", "HapB"), markers))
  }
  H1 <- get_pair(parent1); H2 <- get_pair(parent2)
  x1 <- H1["HapA", ] + H1["HapB", ]; x2 <- H2["HapA", ] + H2["HapB", ]
  d1 <- H1["HapA", ] - H1["HapB", ]; d2 <- H2["HapA", ] - H2["HapB", ]
  beta <- beta[markers]; beta[!is.finite(beta)] <- 0
  R <- recomb_decay_mat
  # Term-2 co-inheritance x within-parent kernel: DH (1-r)(1-2r) = (1+R)R/2; RIL R* itself.
  K2 <- if (identical(target, "DH")) 0.5 * (1 + R) * R else R
  # Full progeny genotypic covariance K such that Var(G) = beta' K beta.
  Kfull <- R * (0.25 * tcrossprod(x1 - x2)) + 0.5 * K2 * (tcrossprod(d1) + tcrossprod(d2))
  Kfull <- 0.5 * (Kfull + t(Kfull))
  vpm <- as.numeric(crossprod(beta, Kfull %*% beta))
  pmv <- vpm
  if (!is.null(beta_cov)) {
    bc <- beta_cov[markers, markers, drop = FALSE]
    pmv <- vpm + sum(Kfull * bc)
  }
  c(VPM = max(0, vpm), PMV = max(0, pmv))
}

# Overwrite vpm/pmv for crosses whose parents carry residual heterozygosity with
# the exact phased-haplotype formula. Only het-parent crosses are touched;
# inbred-parent crosses keep the a'Ra path byte-identical. Uses a DENSE
# recombination kernel (the het correction is exact only densely; the a'Ra
# window_cm truncation does not apply to it). Phase input is validated strictly
# against the dosage matrix; incomplete or inconsistent phase is an error.
ng_apply_het_parent_correction <- function(vpm, pmv, pmv_full,
                                           phased_haplotypes, sorted, pairs,
                                           target, recomb_model,
                                           beta_var, beta_cov_full = NULL) {
  if (!is.array(phased_haplotypes) || length(dim(phased_haplotypes)) != 3L ||
      dim(phased_haplotypes)[2L] != 2L || is.null(dimnames(phased_haplotypes)[[1L]]) ||
      is.null(dimnames(phased_haplotypes)[[3L]])) {
    ng_stop("phased_haplotypes must be a named 3D array [sample, 2 homologs, marker]")
  }
  markers <- sorted$marker_map$marker
  ph <- phased_haplotypes
  ids_ph <- dimnames(ph)[[1L]]
  marker_ph <- dimnames(ph)[[3L]]
  if (anyDuplicated(ids_ph) || anyDuplicated(marker_ph)) {
    ng_stop("phased_haplotypes sample and marker names must be unique")
  }
  missing_markers <- setdiff(markers, marker_ph)
  if (length(missing_markers)) {
    ng_stop("phased_haplotypes is missing scored markers: ",
            paste(utils::head(missing_markers, 6L), collapse = ", "))
  }
  required_parents <- unique(c(as.character(pairs$parent1), as.character(pairs$parent2)))
  missing_parents <- setdiff(required_parents, ids_ph)
  if (length(missing_parents)) {
    ng_stop("phased_haplotypes is missing scored parents: ",
            paste(utils::head(missing_parents, 6L), collapse = ", "))
  }
  ph <- ph[required_parents, , markers, drop = FALSE]
  if (any(!is.finite(ph)) || any(abs(ph - round(ph)) > 1e-8) || any(!(round(ph) %in% c(0, 1)))) {
    ng_stop("phased_haplotypes must contain finite allele states 0/1")
  }
  dosage_from_phase <- ph[, 1L, , drop = FALSE] + ph[, 2L, , drop = FALSE]
  dim(dosage_from_phase) <- c(length(required_parents), length(markers))
  dimnames(dosage_from_phase) <- list(required_parents, markers)
  observed_dosage <- sorted$geno[required_parents, markers, drop = FALSE]
  if (any(!is.finite(observed_dosage)) || any(abs(dosage_from_phase - observed_dosage) > 1e-8)) {
    ng_stop("phased_haplotypes allele sums must exactly match the scored genotype dosages")
  }
  ids_ph <- dimnames(ph)[[1L]]
  het_of <- function(id) id %in% ids_ph && any(ph[id, 1L, ] != ph[id, 2L, ])
  need <- required_parents
  hap_rows <- matrix(0, nrow = 2L * length(need), ncol = length(markers),
                     dimnames = list(as.vector(rbind(paste0(need, "_HapA"),
                                                     paste0(need, "_HapB"))), markers))
  for (k in seq_along(need)) {
    hap_rows[2L * k - 1L, ] <- ph[need[k], 1L, ]
    hap_rows[2L * k,      ] <- ph[need[k], 2L, ]
  }
  R <- ng_recomb_decay_matrix(sorted$marker_map, model = recomb_model, target = target)
  beta <- sorted$effects[markers]
  bc_diag <- diag(as.numeric(beta_var[markers]), length(markers)); dimnames(bc_diag) <- list(markers, markers)
  bc_full <- if (!is.null(beta_cov_full)) beta_cov_full[markers, markers, drop = FALSE] else NULL
  for (i in seq_len(nrow(pairs))) {
    p1 <- as.character(pairs$parent1[i]); p2 <- as.character(pairs$parent2[i])
    if ((het_of(p1) || het_of(p2)) && p1 %in% ids_ph && p2 %in% ids_ph) {
      r_d <- ng_gms_additive_var_general(p1, p2, hap_rows, beta, R, beta_cov = bc_diag, target = target)
      vpm[i] <- r_d[["VPM"]]; pmv[i] <- r_d[["PMV"]]
      if (!is.null(bc_full) && !is.null(pmv_full)) {
        pmv_full[i] <- ng_gms_additive_var_general(p1, p2, hap_rows, beta, R, beta_cov = bc_full, target = target)[["PMV"]]
      }
    }
  }
  list(vpm = vpm, pmv = pmv, pmv_full = pmv_full)
}
