# v0.3.1 — Multivariate threshold probability Pr(G_c in A)
#
# Multi-trait extension of ng_p_superior_progeny(). A breeder can ask: what
# is the probability that at least one of k progeny falls inside a
# rectangular target region A = prod_t [tau_t^lo, tau_t^hi] across all
# traits at once? Closed form under the multivariate-normal progeny
# assumption (a decision approximation, not a biological claim; OK for
# polygenic scores per Allier et al. 2019).
#
# Formula:
#   p_one = Pr(one progeny in A) = pmvnorm(lower, upper, mu, Sigma_c)
#   p_max = 1 - (1 - p_one)^k          (order-statistic for max of k)
# Log-space via log1p(-p_one) to match ng_p_superior_progeny()'s
# underflow guard at R/30_posterior_prediction.R:249-261.

.ng_pmvnorm <- function(lower, upper, mu, Sigma) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    ng_stop("mvtnorm is required for multivariate threshold probabilities; ",
            "install.packages(\"mvtnorm\") or use the single-trait path.")
  }
  as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper,
                              mean = mu, sigma = Sigma))
}

ng_p_superior_progeny_multitrait <- function(mu, Sigma_c, tau_lower, tau_upper,
                                              k_progeny) {
  mu <- as.numeric(mu)
  tau_lower <- as.numeric(tau_lower)
  tau_upper <- as.numeric(tau_upper)
  k_progeny <- as.numeric(k_progeny)
  t <- length(mu)
  if (t == 0L) ng_stop("mu must have length >= 1")
  if (length(tau_lower) != t || length(tau_upper) != t) {
    ng_stop("tau_lower and tau_upper must have length(mu) = ", t)
  }
  if (!is.matrix(Sigma_c) || nrow(Sigma_c) != t || ncol(Sigma_c) != t) {
    ng_stop("Sigma_c must be a t x t matrix matching length(mu) = ", t)
  }
  if (!is.finite(k_progeny) || k_progeny < 1) {
    ng_stop("k_progeny must be >= 1")
  }
  if (any(tau_upper <= tau_lower)) return(0)
  diag_var <- diag(Sigma_c)
  if (any(diag_var <= 0)) {
    zero_traits <- which(diag_var <= 0)
    ok <- all(mu[zero_traits] >= tau_lower[zero_traits] &
              mu[zero_traits] <= tau_upper[zero_traits])
    if (!ok) return(0)
    keep <- setdiff(seq_len(t), zero_traits)
    # all traits are deterministic (Sigma diagonal = 0) and all fall inside
    # their intervals (verified above) -> probability is 1.
    if (!length(keep)) return(1)
    return(ng_p_superior_progeny_multitrait(
      mu = mu[keep],
      Sigma_c = Sigma_c[keep, keep, drop = FALSE],
      tau_lower = tau_lower[keep],
      tau_upper = tau_upper[keep],
      k_progeny = k_progeny
    ))
  }
  p_one <- .ng_pmvnorm(tau_lower, tau_upper, mu, Sigma_c)
  p_one <- max(0, min(1, p_one))
  log_complement <- log1p(-p_one)
  one_minus_p_max <- exp(k_progeny * log_complement)
  out <- 1 - one_minus_p_max
  max(0, min(1, out))
}

# Build a cross-level trait covariance Sigma_c from per-trait within-family
# variances and an optional t x t genetic correlation matrix G_hat.
#   - per_trait_var: length-t vector, typically the per-trait pmv
#     for the cross.
#   - G_hat: optional t x t matrix. May be a covariance or correlation
#     matrix; we coerce to correlation form via cov2cor() so the user can
#     pass the raw output of ng_estimate_genetic_covariance().
# Returns Sigma_c = D R D where D = diag(sqrt(per_trait_var)) and R is the
# correlation form of G_hat, projected to nearest PSD by eigen-clip.
ng_build_cross_trait_covariance <- function(per_trait_var, G_hat = NULL) {
  v <- as.numeric(per_trait_var)
  t <- length(v)
  if (t == 0L) ng_stop("per_trait_var must have length >= 1")
  if (any(!is.finite(v))) ng_stop("per_trait_var contains non-finite entries")
  if (any(v < 0)) ng_stop("per_trait_var contains negative values; variances must be >= 0")
  if (is.null(G_hat)) return(diag(pmax(v, 0), nrow = t))
  G_hat <- as.matrix(G_hat)
  if (nrow(G_hat) != t || ncol(G_hat) != t) {
    ng_stop("G_hat must be ", t, " x ", t, " to match per_trait_var")
  }
  G_diag <- diag(G_hat)
  if (any(G_diag <= 0)) ng_stop("G_hat must have positive diagonal")
  R <- stats::cov2cor(G_hat)
  # Postcondition: a valid covariance matrix must produce correlations in
  # [-1, 1]. cov2cor() doesn't enforce this when G_hat itself is invalid
  # (off-diagonals exceeding the geometric mean of diagonals). Reject
  # explicitly so we don't silently corrupt the Sigma_c diagonal during
  # the D R D construction.
  off_diag <- R[lower.tri(R)]
  if (any(abs(off_diag) > 1 + 1e-8)) {
    ng_stop("G_hat is not a valid covariance matrix: cov2cor produced ",
            "off-diagonal correlations with |R[i,j]| > 1 (max = ",
            sprintf("%.4f", max(abs(off_diag))), "). Ensure |G[i,j]| <= ",
            "sqrt(G[i,i] * G[j,j]) for all i != j.")
  }
  D <- sqrt(pmax(v, 0))
  Sigma <- outer(D, D) * R
  Sigma <- (Sigma + t(Sigma)) / 2
  eig <- eigen(Sigma, symmetric = TRUE)
  eig$values <- pmax(eig$values, 1e-10)
  eig$vectors %*% diag(eig$values, nrow = t) %*% t(eig$vectors)
}

# EXACT within-family cross-trait covariance for each cross: the recombination-aware two-trait
# generalization of the single-trait DH/RIL variance a'Ra,
#   Cov(trait_t, trait_s | cross i x j) = a_t' R a_s,  a_{t,k} = 0.5*(x_ik - x_jk) * beta_{t,k},
# with R the progeny-target recombination-decay matrix (same kernel ng_score_crosses uses for the
# per-trait variance). The diagonal (t == s) reproduces vpm exactly. Use this instead of
# the POPULATION genetic correlation (ng_build_cross_trait_covariance with G_hat) when the true
# within-family cross-trait covariance is wanted for the multi-trait threshold probability.
#
# `betas` is an m x T matrix of per-trait marker effects (rows named by marker to align with
# colnames(geno); columns named by trait). Returns a data.frame: parent1, parent2, and one column
# per trait pair -- wf_var_<trait> (diagonal) and wf_cov_<t>_<s> (t before s) -- plus a
# `trait_names` attribute.
ng_cross_trait_within_family_cov <- function(geno, betas, marker_map,
                                             ids = rownames(geno), pairs = NULL,
                                             target = c("DH", "RIL"),
                                             recomb_model = c("haldane", "kosambi"),
                                             window_cm = Inf) {
  target <- match.arg(target); recomb_model <- match.arg(recomb_model)
  geno <- as.matrix(geno); storage.mode(geno) <- "double"
  if (is.null(rownames(geno))) rownames(geno) <- ids
  betas <- as.matrix(betas)
  if (is.null(colnames(betas))) colnames(betas) <- paste0("trait", seq_len(ncol(betas)))
  if (!is.null(rownames(betas))) {
    common <- intersect(colnames(geno), rownames(betas))
    if (length(common) < 2L) ng_stop("betas and geno share fewer than 2 markers by name")
    geno <- geno[, common, drop = FALSE]; betas <- betas[common, , drop = FALSE]
  } else if (nrow(betas) != ncol(geno)) {
    ng_stop("betas must have one row per genotype marker (or be named by marker)")
  }
  mm <- ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model)
  ord <- order(mm$chr_index, mm$pos_cm, mm$marker)
  mm <- mm[ord, , drop = FALSE]; geno <- geno[, ord, drop = FALSE]; betas <- betas[ord, , drop = FALSE]

  if (is.null(pairs)) pairs <- ng_make_pairs(rownames(geno))
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  p1 <- match(as.character(pairs$parent1), rownames(geno))
  p2 <- match(as.character(pairs$parent2), rownames(geno))
  if (anyNA(p1) || anyNA(p2)) ng_stop("pairs reference ids not present in geno")

  D <- 0.5 * (geno[p1, , drop = FALSE] - geno[p2, , drop = FALSE])   # n_pairs x m contrasts
  D[!is.finite(D)] <- 0
  Tn <- ncol(betas); trait_names <- colnames(betas); npairs <- nrow(pairs)
  A <- lapply(seq_len(Tn), function(j) { Aj <- t(D) * betas[, j]; Aj[!is.finite(Aj)] <- 0; Aj })  # m x n_pairs

  colpair <- which(upper.tri(matrix(0, Tn, Tn), diag = TRUE), arr.ind = TRUE)  # (ti, si), ti <= si
  cov_ts <- matrix(0, npairs, nrow(colpair))
  chr <- mm$chr_index
  for (cc in unique(chr)) {                                   # cross-chromosome R is exactly 0
    idx <- which(chr == cc)
    Rc <- ng_recomb_decay_matrix(mm[idx, , drop = FALSE], model = recomb_model, target = target)
    if (is.finite(window_cm) && window_cm > 0) {
      dcm <- outer(mm$pos_cm[idx], mm$pos_cm[idx], function(a, b) abs(a - b))
      Rc[dcm > window_cm] <- 0
    }
    RcA <- lapply(A, function(Aj) Rc %*% Aj[idx, , drop = FALSE])
    for (cp in seq_len(nrow(colpair))) {
      ti <- colpair[cp, 1L]; si <- colpair[cp, 2L]
      cov_ts[, cp] <- cov_ts[, cp] + colSums(A[[ti]][idx, , drop = FALSE] * RcA[[si]])
    }
  }

  out <- data.frame(parent1 = as.character(pairs$parent1),
                    parent2 = as.character(pairs$parent2), stringsAsFactors = FALSE)
  for (cp in seq_len(nrow(colpair))) {
    ti <- colpair[cp, 1L]; si <- colpair[cp, 2L]
    nm <- if (ti == si) paste0("wf_var_", trait_names[ti]) else
      paste0("wf_cov_", trait_names[ti], "_", trait_names[si])
    out[[nm]] <- cov_ts[, cp]
  }
  attr(out, "trait_names") <- trait_names
  out
}

# Assemble a per-cross list of T x T within-family covariance matrices from the EXACT cross-trait
# covariance table (ng_cross_trait_within_family_cov), matched to `scores` rows by unordered parent
# pair, in the trait order `trait_order`. Each matrix is symmetrized and PSD-projected.
ng_multitrait_exact_sigma_list <- function(scores, trait_order, cross_trait_cov) {
  ctc <- as.data.frame(cross_trait_cov, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(scores)) ||
      !all(c("parent1", "parent2") %in% names(ctc))) {
    ng_stop("cross_trait_cov and scores both need parent1/parent2 to match crosses")
  }
  pkey <- function(df) paste(pmin(as.character(df$parent1), as.character(df$parent2)),
                             pmax(as.character(df$parent1), as.character(df$parent2)), sep = "||")
  idx <- match(pkey(scores), pkey(ctc))
  if (anyNA(idx)) ng_stop("cross_trait_cov is missing some crosses present in scores")
  t_n <- length(trait_order)
  cov_at <- function(row, a, b) {
    if (a == b) return(as.numeric(row[[paste0("wf_var_", trait_order[a])]]))
    nm1 <- paste0("wf_cov_", trait_order[a], "_", trait_order[b])
    nm2 <- paste0("wf_cov_", trait_order[b], "_", trait_order[a])
    if (nm1 %in% names(row)) as.numeric(row[[nm1]]) else if (nm2 %in% names(row)) as.numeric(row[[nm2]]) else
      ng_stop("cross_trait_cov missing covariance column for traits ",
              trait_order[a], " / ", trait_order[b])
  }
  lapply(idx, function(r) {
    row <- ctc[r, , drop = FALSE]
    S <- matrix(0, t_n, t_n)
    for (a in seq_len(t_n)) for (b in seq_len(t_n)) S[a, b] <- cov_at(row, a, b)
    S <- (S + t(S)) / 2
    eig <- eigen(S, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    eig$vectors %*% diag(eig$values, nrow = t_n) %*% t(eig$vectors)
  })
}

# Add a P(superior progeny in A) column to an existing multi-trait score
# table. Per row, builds Sigma_c from per-trait variances (plus optional
# G_hat for off-diagonals) and calls ng_p_superior_progeny_multitrait().
#
# trait_specs is a data.frame with columns:
#   trait     - trait label
#   mean_col  - score-table column with per-trait predicted mean
#   var_col   - score-table column with per-trait variance (pmv-like)
#
# Order of rows in trait_specs determines the order of tau_lower / tau_upper
# and (if supplied) of G_hat.
#
# cross_trait_cov (optional): the EXACT within-family cross-trait covariance table from
# ng_cross_trait_within_family_cov(). When supplied, Sigma_c per cross is taken directly from it
# (recombination-aware a_t' R a_s), overriding both var_col and the G_hat population-correlation
# proxy -- the statistically correct within-family covariance for the progeny distribution.
ng_add_p_superior_progeny_multitrait <- function(scores, trait_specs,
                                                  tau_lower, tau_upper,
                                                  k_progeny = 100L,
                                                  G_hat = NULL,
                                                  cross_trait_cov = NULL,
                                                  out_col = "p_superior_progeny_mt") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  trait_specs <- as.data.frame(trait_specs, stringsAsFactors = FALSE)
  required <- c("trait", "mean_col", "var_col")
  missing <- setdiff(required, names(trait_specs))
  if (length(missing)) ng_stop("trait_specs missing columns: ",
                                paste(missing, collapse = ", "))
  t_n <- nrow(trait_specs)
  if (t_n < 1L) ng_stop("trait_specs must have at least one row")
  if (length(tau_lower) != t_n || length(tau_upper) != t_n) {
    ng_stop("tau_lower and tau_upper must have one entry per trait_specs row")
  }
  needed_cols <- c(trait_specs$mean_col, trait_specs$var_col)
  missing_cols <- setdiff(needed_cols, names(scores))
  if (length(missing_cols)) ng_stop("scores missing columns referenced by trait_specs: ",
                                     paste(missing_cols, collapse = ", "))
  mean_mat <- as.matrix(scores[, trait_specs$mean_col, drop = FALSE])
  var_mat  <- as.matrix(scores[, trait_specs$var_col, drop = FALSE])
  storage.mode(mean_mat) <- "double"
  storage.mode(var_mat)  <- "double"
  exact_sigma <- if (!is.null(cross_trait_cov)) {
    ng_multitrait_exact_sigma_list(scores, trait_specs$trait, cross_trait_cov)
  } else NULL
  out <- numeric(nrow(scores))
  for (i in seq_len(nrow(scores))) {
    Sigma_i <- if (!is.null(exact_sigma)) exact_sigma[[i]] else
      ng_build_cross_trait_covariance(per_trait_var = var_mat[i, ], G_hat = G_hat)
    out[[i]] <- ng_p_superior_progeny_multitrait(
      mu = mean_mat[i, ], Sigma_c = Sigma_i,
      tau_lower = tau_lower, tau_upper = tau_upper, k_progeny = k_progeny
    )
  }
  scores[[out_col]] <- out
  attr(scores, "p_superior_progeny_mt") <- list(
    traits = trait_specs$trait, mean_cols = trait_specs$mean_col,
    var_cols = trait_specs$var_col, tau_lower = tau_lower,
    tau_upper = tau_upper, k_progeny = k_progeny,
    G_hat_supplied = !is.null(G_hat),
    exact_cross_trait_cov = !is.null(cross_trait_cov), out_col = out_col
  )
  scores
}
