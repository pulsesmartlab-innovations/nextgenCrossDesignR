# Risk (confidence) + portfolio (level x genetic-SD) decision layer for cross priority.
# Single-trait, posterior-OFF core. See docs/design/2026-07-25-cross-priority-*-design.md.

# Per-cross prediction-error variance of the mid-parent GEBV:
#   Var(1/2 (g_p1 + g_p2)) = 1/4 (xc_p1 + xc_p2)' Sigma_beta (xc_p1 + xc_p2)
# on centered genotypes xc = geno - marker_mean. Returns NA vector if beta_cov is NULL.
ng_midparent_pev <- function(geno, pairs, beta_cov, marker_mean = NULL) {
  geno <- as.matrix(geno)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (is.null(beta_cov)) return(rep(NA_real_, nrow(pairs)))
  if (is.null(marker_mean)) marker_mean <- colMeans(geno, na.rm = TRUE)
  Xc <- sweep(geno, 2L, marker_mean, "-")
  i1 <- match(as.character(pairs$parent1), rownames(geno))
  i2 <- match(as.character(pairs$parent2), rownames(geno))
  if (anyNA(i1) || anyNA(i2))
    ng_stop("ng_midparent_pev: pair parents not found in geno rownames")
  s <- Xc[i1, , drop = FALSE] + Xc[i2, , drop = FALSE]
  unname(0.25 * rowSums((s %*% beta_cov) * s))
}

# Resolve the posterior-OFF confidence from the mid-parent PEV.
#   cross_confidence in [0,1], higher = more trustworthy (1 - normalized spread)
#   risk_bin low/med/high = within-plan tertiles of the spread (low spread = low risk)
#   confidence_method: *_partial when the merit's sqrt(X) term is effect-based (mean PEV
#     does not cover the variance-estimation error, spec 5.4); plain otherwise.
ng_cross_confidence <- function(pev, effect_based_x = TRUE,
                                method_prefix = "midparent_pev", spread = NULL) {
  # `spread` supplies the uncertainty scale DIRECTLY (e.g. a posterior SD of the ranked value)
  # instead of deriving it as sqrt(pev). Everything downstream -- normalization, tertiles -- is
  # identical; only the source of the scale differs.
  if (!is.null(spread)) pev <- pmax(as.numeric(spread), 0)^2
  n <- length(pev)
  na_bin <- factor(rep(NA_character_, n), levels = c("low", "med", "high"), ordered = TRUE)
  if (is.null(pev) || !any(is.finite(pev)) ||
      length(unique(pev[is.finite(pev)])) < 2L) {
    return(list(cross_confidence = rep(NA_real_, n), relative_precision = rep(NA_real_, n),
                risk_bin = na_bin, precision_bin = na_bin,
                confidence_method = "relative_precision_unavailable",
                is_calibrated = FALSE))
  }
  spread <- sqrt(pmax(as.numeric(pev), 0))
  fin <- is.finite(spread)
  rng <- range(spread[fin])
  # Degenerate spread contains no relative-precision information.
  if (diff(rng) <= 0) {
    return(list(cross_confidence = rep(NA_real_, n), relative_precision = rep(NA_real_, n),
                risk_bin = na_bin, precision_bin = na_bin,
                confidence_method = "relative_precision_unavailable",
                is_calibrated = FALSE))
  }
  norm <- rep(NA_real_, n)
  norm[fin] <- (spread[fin] - rng[1L]) / (rng[2L] - rng[1L])
  conf <- 1 - norm
  qs <- stats::quantile(spread[fin], c(1/3, 2/3), names = FALSE)
  brk <- unique(c(-Inf, qs, Inf))
  bin <- if (length(brk) == 4L) {
    cut(spread, breaks = brk, labels = c("low", "med", "high"))
  } else {                                   # tied tertiles: rank-thirds fallback
    r <- rank(spread, ties.method = "first", na.last = "keep")
    cut(r, breaks = 3L, labels = c("low", "med", "high"))
  }
  method <- paste0("relative_", if (isTRUE(effect_based_x)) paste0(method_prefix, "_partial") else method_prefix)
  bin <- factor(bin, levels = c("low", "med", "high"), ordered = TRUE)
  list(cross_confidence = conf,
       relative_precision = conf,
       risk_bin = bin,
       precision_bin = bin,
       confidence_method = method,
       is_calibrated = FALSE)
}

# Classify crosses into a 2x2 quadrant (level x upside) based on median cuts.
#   factor with levels = c("breakthrough","workhorse","long_shot","deprioritize")
#   (hi level, hi upside), (hi level, lo upside), (lo level, hi upside), (lo level, lo upside)
ng_cross_portfolio_profile <- function(level, upside,
    labels = c("breakthrough", "workhorse", "long_shot", "deprioritize")) {
  level <- as.numeric(level); upside <- as.numeric(upside)
  hi_l <- level  >= stats::median(level,  na.rm = TRUE)
  hi_u <- upside >= stats::median(upside, na.rm = TRUE)
  out <- ifelse(hi_l & hi_u, labels[1L],
         ifelse(hi_l & !hi_u, labels[2L],
         ifelse(!hi_l & hi_u, labels[3L], labels[4L])))
  out[is.na(level) | is.na(upside)] <- NA
  factor(out, levels = labels)
}

# Attach the portfolio + risk columns to a scored/selected cross table.
#   level = the additive-genetic mid-parent mean column (cross_mean_gebv basis)
#   vpm   = the pure within-family genetic variance column (a'Ra); upside = sqrt(vpm)
ng_annotate_cross_priority <- function(crosses, level, vpm, pev = NULL,
                                       effect_based_x = TRUE,
                                       post_sd = NULL, prob_top_tier = NULL) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  crosses$cross_level  <- as.numeric(level)
  crosses$cross_upside <- sqrt(pmax(as.numeric(vpm), 0))
  # Posterior ON: confidence comes from the posterior SD of the RANKED value, which already
  # covers both the mean and the variance term of the merit -- so it carries no "_partial"
  # caveat, unlike the mid-parent PEV fallback (which covers only the mean).
  use_post <- !is.null(post_sd) && any(is.finite(post_sd))
  cf <- if (use_post) {
    ng_cross_confidence(NULL, effect_based_x = FALSE,
                        method_prefix = "posterior_ci", spread = post_sd)
  } else {
    ng_cross_confidence(pev, effect_based_x = effect_based_x)
  }
  crosses$cross_confidence  <- cf$cross_confidence
  crosses$relative_precision <- cf$relative_precision
  crosses$risk_bin          <- cf$risk_bin
  crosses$precision_bin     <- cf$precision_bin
  crosses$confidence_method <- cf$confidence_method
  crosses$cross_confidence_is_calibrated <- cf$is_calibrated
  crosses$portfolio_profile <- ng_cross_portfolio_profile(crosses$cross_level,
                                                          crosses$cross_upside)
  crosses$portfolio_basis <- "single_trait"
  # prob_top_tier is a JOINT merit x uncertainty quantity ("is this cross genuinely top-N?"),
  # deliberately reported as its own continuous column and NEVER binned into risk_bin or
  # relabelled as confidence -- the two answer different questions.
  crosses$prob_top_tier <- if (is.null(prob_top_tier)) NA_real_ else as.numeric(prob_top_tier)
  crosses
}

# How literally the level/upside decomposition should be read, given the index method.
#   "single_trait"           - level is the trait's own mid-parent GEBV; exact.
#   "linear_index"           - economic_index / desired_gain: coefficients were SOLVED as a
#                              linear index, so level = b'm and upside = sqrt(b'Sb) are the
#                              index's own quantities, mapped back from value_z to raw units.
#   "linearized_rank_index"  - auto / weighted / threshold: multi_trait_score combines
#                              RANK-normalized traits, so no linear b exists in genetic units.
#                              The trait weights are reinterpreted as standardized-unit index
#                              coefficients to define the axes. The quadrant is internally
#                              consistent but is NOT a decomposition of multi_trait_score --
#                              surface it as indicative and never as the ranked merit.
#                              (Deliberate override of the 2026-07-25 design's "omit for
#                              non-linear methods" rule; see that doc's revision note.)
ng_multitrait_portfolio_basis <- function(method) {
  fam <- ng_multitrait_method_family(method)
  if (fam %in% c("economic_index", "desired_gain")) "linear_index" else "linearized_rank_index"
}

# ---- Multi-trait (selection-index) level / upside / risk ----------------------------
#
# For a multi-trait run the portfolio axes live on the SELECTION INDEX, not on any single
# trait. Both axes are derived from ONE raw-unit weight vector w so the level x upside
# quadrant stays coherent (mixing a standardized level with a raw-unit upside would make
# the median split meaningless):
#
#   cross_level  = w' m_i           m_i = per-trait mid-parent GEBV of cross i (raw units)
#   cross_upside = sqrt(w' S_i w)   S_i = within-family cross-trait covariance of cross i
#
# The index is solved on the oriented-standardized scale value_z_k = s_k (x_k - c_k)/scale_k
# (s_k = +1 maximize, -1 minimize; see ng_add_multitrait_score), so the coefficients must be
# pulled back to raw trait units: w_k = coef_k * s_k / scale_k. cross_level is then the index
# of mid-parent GEBVs up to an additive constant -- irrelevant to a median-split quadrant.
#
# NOTE (deliberate): the `weighted` family combines RANK-normalized traits in multi_trait_score,
# which has no linear raw-unit representation. cross_level uses the linear standardized form for
# every method so that level and upside remain consistent. cross_level is therefore a genetic
# quantity on the index, not a monotone transform of multi_trait_score.

# Resolve the raw-unit index weights and the index level from per-trait mid-parent GEBVs.
#   mean_gebv   - n x T numeric matrix (raw trait units), columns in trait order
#   directions  - length-T "maximize"/"minimize" (anything ng_multitrait_direction accepts)
#   coefficients- length-T index coefficients on the oriented-standardized scale
#   center/scale- REQUIRED whenever more than one table is annotated in the same run. The
#     per-trait scale is an IQR of the rows it is handed, so letting each table derive its own
#     would give the selected plan and the candidate pool DIFFERENT index weights -- two
#     incomparable sets of axes for the same run. Resolve the basis once on the candidate pool
#     (ng_multitrait_index_reference) and pass it to every table annotated afterwards.
ng_multitrait_index_basis <- function(mean_gebv, directions, coefficients,
                                      center = NULL, scale = NULL) {
  M <- as.matrix(mean_gebv)
  storage.mode(M) <- "double"
  t_n <- ncol(M)
  if (t_n < 1L) ng_stop("ng_multitrait_index_basis: need at least one trait column")
  directions <- ng_multitrait_direction(directions)
  if (length(directions) != t_n)
    ng_stop("ng_multitrait_index_basis: directions must have one entry per trait column")
  coef <- suppressWarnings(as.numeric(coefficients))
  if (!is.null(names(coefficients)) && !is.null(colnames(M)) &&
      all(colnames(M) %in% names(coefficients))) {
    coef <- suppressWarnings(as.numeric(coefficients[colnames(M)]))
  }
  if (length(coef) != t_n)
    ng_stop("ng_multitrait_index_basis: coefficients must have one entry per trait column")
  coef[!is.finite(coef)] <- 0
  if (!any(coef != 0)) ng_stop("ng_multitrait_index_basis: all index coefficients are zero")

  sgn <- ifelse(directions == "maximize", 1, -1)
  derive <- is.null(center) || is.null(scale)
  if (!derive && (length(center) != t_n || length(scale) != t_n))
    ng_stop("ng_multitrait_index_basis: center/scale must have one entry per trait column")
  if (derive) { scale <- rep(1, t_n); center <- rep(0, t_n) }
  Z <- matrix(0, nrow(M), t_n)
  for (k in seq_len(t_n)) {
    x <- M[, k]
    oriented <- sgn[[k]] * x
    if (derive) {
      scale[[k]] <- ng_multitrait_value_scale(x)
      fin <- oriented[is.finite(oriented)]
      center[[k]] <- if (length(fin)) stats::median(fin, na.rm = TRUE) else 0
    }
    Z[, k] <- (oriented - center[[k]]) / scale[[k]]
  }
  w <- coef * sgn / scale                       # raw-unit weights: level = w'm + const
  names(w) <- colnames(M)
  # A cross with a missing mid-parent GEBV for any component trait has no index level.
  level <- as.numeric(Z %*% coef)
  level[!stats::complete.cases(Z)] <- NA_real_
  list(level = level, w = w, coefficients = coef, sign = sgn,
       center = center, scale = scale)
}

# Index upside sqrt(w' S w) from the EXACT within-family cross-trait covariance columns
# (wf_var_<t> / wf_cov_<t>_<s>, produced by ng_cross_trait_within_family_cov and carried on
# the cross table).
#
# When `vpm` is supplied, S is rescaled so its diagonal reproduces the AUTHORITATIVE per-trait
# variance the run actually reported (S' = D^1/2 C D^1/2, C the correlation implied by a'Ra).
# This keeps two invariants that raw a'Ra alone would break:
#   (a) T == 1 collapses to exactly sqrt(vpm)  -> continuity with the single-trait axis;
#   (b) het-parent (parent_type = "ril" + phased haplotypes) corrected variances are honored,
#       since ng_cross_trait_within_family_cov applies the inbred kernel to the diagonal.
# Assemble the rescaled index weights W and the S %*% W product from the within-family
# covariance columns, per cross. The index variance is then rowSums(W * SW) = w'Sw. Both the
# upside and the per-trait variance attribution read this, so they are guaranteed to be
# describing the same S -- an attribution that did not sum back to the reported upside would be
# worse than no attribution at all.
ng_mt_index_sw <- function(crosses, trait_order, w, vpm = NULL) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  n <- nrow(crosses); t_n <- length(trait_order)
  var_col <- paste0("wf_var_", trait_order)
  missing_var <- var_col[!(var_col %in% names(crosses))]
  if (length(missing_var))
    ng_stop("ng_mt_index_sw: missing within-family variance column(s): ",
            paste(missing_var, collapse = ", "))
  V <- matrix(vapply(var_col, function(cc) suppressWarnings(as.numeric(crosses[[cc]])),
                     numeric(n)), nrow = n, ncol = t_n)

  # Row-wise weight rescaling is algebraically identical to rescaling S itself, and avoids
  # materializing one T x T matrix per cross.
  W <- matrix(rep(as.numeric(w), each = n), nrow = n, ncol = t_n)
  if (!is.null(vpm)) {
    P <- matrix(suppressWarnings(as.numeric(as.matrix(vpm))), nrow = n, ncol = t_n)
    ratio <- ifelse(is.finite(V) & V > 0 & is.finite(P) & P >= 0, sqrt(P / V), 1)
    W <- W * ratio
  }
  SW <- W * V                                   # diagonal contribution
  if (t_n > 1L) {
    for (a in seq_len(t_n - 1L)) for (b in seq(a + 1L, t_n)) {
      nm1 <- paste0("wf_cov_", trait_order[a], "_", trait_order[b])
      nm2 <- paste0("wf_cov_", trait_order[b], "_", trait_order[a])
      cc <- if (nm1 %in% names(crosses)) nm1 else if (nm2 %in% names(crosses)) nm2 else
        ng_stop("ng_mt_index_sw: missing covariance column for ",
                trait_order[a], " / ", trait_order[b])
      cv <- suppressWarnings(as.numeric(crosses[[cc]]))
      SW[, a] <- SW[, a] + W[, b] * cv
      SW[, b] <- SW[, b] + W[, a] * cv
    }
  }
  list(W = W, SW = SW, total = rowSums(W * SW))
}

ng_multitrait_index_upside <- function(crosses, trait_order, w, vpm = NULL) {
  if (!nrow(as.data.frame(crosses))) return(numeric(0))
  sqrt(pmax(ng_mt_index_sw(crosses, trait_order, w, vpm)$total, 0))
}

# Per-trait share of the index VARIANCE: w_k (Sw)_k / w'Sw. Shares sum to 1 by construction.
#
# A share can be NEGATIVE, and that is the informative case, not an error: a trait that is
# antagonistic to the rest of the index within the family REMOVES spread from the index (its
# covariance terms are negative), so it shows as a negative contribution. Read as "protein
# narrows the index spread by 15%", which is exactly the yield/protein intuition.
ng_multitrait_variance_shares <- function(crosses, trait_order, w, vpm = NULL) {
  s <- ng_mt_index_sw(crosses, trait_order, w, vpm)
  contrib <- s$W * s$SW
  denom <- ifelse(is.finite(s$total) & s$total != 0, s$total, NA_real_)
  out <- contrib / denom
  colnames(out) <- trait_order
  out
}

# Per-trait share of the index PREDICTION-ERROR variance: w_k^2 PEV_k / sum_j w_j^2 PEV_j.
# This is the exact decomposition of the block-diagonal index PEV, so shares sum to 1 and are
# always in [0, 1] (unlike the variance shares above, PEV has no cross-trait terms to go
# negative).
ng_multitrait_pev_shares <- function(pev, w, trait_order = NULL) {
  P <- as.matrix(pev)
  storage.mode(P) <- "double"
  if (ncol(P) != length(w))
    ng_stop("ng_multitrait_pev_shares: pev must have one column per trait weight")
  contrib <- sweep(P, 2L, as.numeric(w)^2, "*")
  tot <- rowSums(contrib)
  out <- contrib / ifelse(is.finite(tot) & tot > 0, tot, NA_real_)
  if (!is.null(trait_order)) colnames(out) <- trait_order
  out
}

# Which component trait is driving this cross's index uncertainty, and by how much.
#
# DELIBERATELY an attribution, NOT an adjustment to cross_confidence. sum_k w_k^2 PEV_k is
# already the correct propagation of marker-effect estimation error into the index: a lightly
# weighted trait contributes little because it genuinely moves the index little, and a lightly
# weighted trait with a huge PEV still dominates the sum on its own. Flooring confidence by the
# worst component would double-count uncertainty the weighting has already handled, and would
# make risk_bin stop meaning "how well is this cross's INDEX estimated".
#
# The related breeder question -- "is this cross unacceptable on a minor trait?" -- is a
# THRESHOLD question, and the package answers it separately and better via trait_checks (R/44)
# and the multi-trait min_value/max_value thresholds. What was actually missing is the ability
# to answer "high risk -- because of WHICH trait", which is what this provides.
ng_multitrait_risk_driver <- function(shares, trait_order) {
  S <- as.matrix(shares)
  n <- nrow(S)
  trait <- rep(NA_character_, n); share <- rep(NA_real_, n)
  ok <- which(rowSums(is.finite(S)) > 0L)
  if (length(ok)) {
    idx <- max.col(replace(S[ok, , drop = FALSE], !is.finite(S[ok, , drop = FALSE]), -Inf),
                   ties.method = "first")
    trait[ok] <- trait_order[idx]
    share[ok] <- S[cbind(ok, idx)]
  }
  list(trait = trait, share = share)
}

# Index mid-parent PEV from the per-trait mid-parent PEVs.
#
# APPROXIMATION (document it wherever the number is shown): each trait is fitted by an
# INDEPENDENT univariate ridge, so no cross-trait covariance of the marker-effect estimation
# errors exists to draw on. This treats Cov(beta_hat_t, beta_hat_s) = 0 for t != s, i.e.
#   PEV(index) = sum_k w_k^2 PEV_k
# which is exact under independent fits and biased LOW when the trait fits share information
# (correlated residuals on the same X). It is a confidence RANKING input, not a calibrated
# interval, and ng_cross_confidence only ever uses its within-plan ordering.
ng_multitrait_index_pev <- function(pev, w) {
  P <- as.matrix(pev)
  storage.mode(P) <- "double"
  if (ncol(P) != length(w))
    ng_stop("ng_multitrait_index_pev: pev must have one column per trait weight")
  as.numeric(rowSums(sweep(P, 2L, as.numeric(w)^2, "*")))
}

# Attach the portfolio + risk columns for a MULTI-TRAIT run. Mirrors
# ng_annotate_cross_priority column-for-column so the frontend panel is identical; the
# difference is entirely in how level/upside/pev are resolved (index, not single trait).
# Resolve the index basis ONCE, on the reference table (the full candidate pool), so every
# table annotated in the run shares one set of axes. Returns the object to hand to
# ng_annotate_cross_priority_multitrait(basis = ).
ng_multitrait_index_reference <- function(crosses, trait_order, mean_gebv_cols,
                                          directions, coefficients) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  n <- nrow(crosses)
  if (!all(mean_gebv_cols %in% names(crosses)))
    ng_stop("ng_multitrait_index_reference: missing mid-parent GEBV column(s): ",
            paste(setdiff(mean_gebv_cols, names(crosses)), collapse = ", "))
  M <- matrix(vapply(mean_gebv_cols, function(cc)
    suppressWarnings(as.numeric(crosses[[cc]])), numeric(n)), nrow = n,
    dimnames = list(NULL, trait_order))
  ng_multitrait_index_basis(M, directions, coefficients)
}

ng_annotate_cross_priority_multitrait <- function(crosses, trait_order, mean_gebv_cols,
                                                  vpm_cols, directions, coefficients,
                                                  pev_cols = NULL, effect_based_x = TRUE,
                                                  index_method = NA_character_,
                                                  basis = NULL) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  n <- nrow(crosses)
  need <- c(mean_gebv_cols, vpm_cols)
  if (!all(need %in% names(crosses)))
    ng_stop("ng_annotate_cross_priority_multitrait: missing per-trait column(s): ",
            paste(setdiff(need, names(crosses)), collapse = ", "))
  grab <- function(cols) matrix(vapply(cols, function(cc)
    suppressWarnings(as.numeric(crosses[[cc]])), numeric(n)), nrow = n,
    dimnames = list(NULL, trait_order))

  basis <- ng_multitrait_index_basis(grab(mean_gebv_cols), directions, coefficients,
                                     center = basis$center, scale = basis$scale)
  crosses$cross_level  <- basis$level
  crosses$cross_upside <- ng_multitrait_index_upside(crosses, trait_order, basis$w,
                                                     vpm = grab(vpm_cols))

  pev <- rep(NA_real_, n)
  risk_driver <- list(trait = rep(NA_character_, n), share = rep(NA_real_, n))
  if (!is.null(pev_cols) && all(pev_cols %in% names(crosses))) {
    pev_mat <- grab(pev_cols)
    pev <- ng_multitrait_index_pev(pev_mat, basis$w)
    # Attribution, not adjustment -- see ng_multitrait_risk_driver.
    risk_driver <- ng_multitrait_risk_driver(
      ng_multitrait_pev_shares(pev_mat, basis$w, trait_order), trait_order)
  }
  cf <- ng_cross_confidence(pev, effect_based_x = effect_based_x,
                            method_prefix = "midparent_pev_index")
  crosses$cross_confidence  <- cf$cross_confidence
  crosses$relative_precision <- cf$relative_precision
  crosses$risk_bin          <- cf$risk_bin
  crosses$precision_bin     <- cf$precision_bin
  crosses$confidence_method <- cf$confidence_method
  crosses$cross_confidence_is_calibrated <- cf$is_calibrated
  # "This cross is high risk -- because of protein (62% of the index PEV)."
  crosses$risk_driver_trait <- risk_driver$trait
  crosses$risk_driver_share <- risk_driver$share
  crosses$portfolio_profile <- ng_cross_portfolio_profile(crosses$cross_level,
                                                          crosses$cross_upside)
  # Provenance for the frontend: which weights produced these axes, in what units, and how
  # literally the decomposition may be read (see ng_multitrait_portfolio_basis).
  crosses$portfolio_basis <- if (is.na(index_method)) NA_character_ else
    ng_multitrait_portfolio_basis(index_method)
  crosses$cross_level_rule <- sprintf(
    "index_method=%s; portfolio_basis=%s; basis=mid_parent_gebv; weights=%s; units=index",
    if (is.na(index_method)) "unknown" else index_method,
    if (is.na(index_method)) "unknown" else ng_multitrait_portfolio_basis(index_method),
    paste(sprintf("%s=%.4g", trait_order, basis$w), collapse = "/")
  )
  attr(crosses, "index_basis") <- basis
  crosses
}

# Tier x profile cross-tab summary: aggregate counts and mean metrics.
#   Input: crosses with priority_tier and portfolio_profile columns.
#   Output: data.frame with one row per (tier, profile) combination found.
ng_cross_portfolio_summary <- function(crosses,
                                       tier_col = "priority_tier",
                                       profile_col = "portfolio_profile") {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  if (!all(c(tier_col, profile_col) %in% names(crosses)))
    ng_stop("ng_cross_portfolio_summary: missing tier/profile column")
  empty_out <- function() {
    data.frame(priority_tier = character(0), portfolio_profile = character(0),
               n = integer(0), mean_level = numeric(0), mean_upside = numeric(0),
               mean_confidence = numeric(0), stringsAsFactors = FALSE)
  }
  if (nrow(crosses) == 0L) return(empty_out())
  # NA tier/profile is a legitimate (unclassified) group, not a row to drop.
  tier    <- as.character(crosses[[tier_col]])
  profile <- as.character(crosses[[profile_col]])
  tier[is.na(tier)]       <- "(unclassified)"
  profile[is.na(profile)] <- "(unclassified)"
  key <- interaction(tier, profile, drop = TRUE, sep = "\r")
  parts <- split(seq_len(nrow(crosses)), key)
  if (length(parts) == 0L) return(empty_out())
  safe_mean <- function(x) {
    m <- mean(x, na.rm = TRUE)
    if (!is.finite(m)) NA_real_ else m
  }
  rows <- lapply(parts, function(idx) {
    kv <- strsplit(as.character(key[idx[1L]]), "\r", fixed = TRUE)[[1L]]
    data.frame(priority_tier = kv[1L], portfolio_profile = kv[2L], n = length(idx),
               mean_level = safe_mean(crosses$cross_level[idx]),
               mean_upside = safe_mean(crosses$cross_upside[idx]),
               mean_confidence = safe_mean(crosses$cross_confidence[idx]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
