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
ng_cross_confidence <- function(pev, effect_based_x = TRUE) {
  n <- length(pev)
  na_bin <- factor(rep(NA_character_, n), levels = c("low", "med", "high"), ordered = TRUE)
  if (is.null(pev) || !any(is.finite(pev)) ||
      length(unique(pev[is.finite(pev)])) < 2L) {
    return(list(cross_confidence = rep(NA_real_, n), risk_bin = na_bin,
                confidence_method = "reliability"))
  }
  spread <- sqrt(pmax(as.numeric(pev), 0))
  fin <- is.finite(spread)
  rng <- range(spread[fin])
  # Degenerate spread (all clamped to same value): fallback to reliability
  if (diff(rng) <= 0) {
    return(list(cross_confidence = rep(NA_real_, n), risk_bin = na_bin,
                confidence_method = "reliability"))
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
  method <- if (isTRUE(effect_based_x)) "midparent_pev_partial" else "midparent_pev"
  list(cross_confidence = conf,
       risk_bin = factor(bin, levels = c("low", "med", "high"), ordered = TRUE),
       confidence_method = method)
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
                                       effect_based_x = TRUE) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  crosses$cross_level  <- as.numeric(level)
  crosses$cross_upside <- sqrt(pmax(as.numeric(vpm), 0))
  cf <- ng_cross_confidence(pev, effect_based_x = effect_based_x)
  crosses$cross_confidence  <- cf$cross_confidence
  crosses$risk_bin          <- cf$risk_bin
  crosses$confidence_method <- cf$confidence_method
  crosses$portfolio_profile <- ng_cross_portfolio_profile(crosses$cross_level,
                                                          crosses$cross_upside)
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
