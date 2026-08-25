# ---------------------------------------------------------------------------
# UCPC selected-fraction kernel (Allier et al. 2019 baseline), Phase 1.
#
# Allier et al. 2019 (Theor Appl Genet) define the "usefulness criterion
# parental contribution" (UCPC): given a truncation-selected fraction of a
# biparental progeny population, what fraction of the selected genome comes
# from each parent (c_p1, c_p2), and what expected heterozygosity does the
# selected fraction retain (He^(i))? This file implements that baseline for
# a single cross (ng_ucpc_cross()) and for a crossing plan
# (ng_ucpc_selected_diversity()), for BOTH DH and RIL targets.
#
# Allier et al. (2019) explicitly present DH-k and RIL-k generations. This
# implementation currently supports DH from the F1 and infinite-selfing RIL;
# the derivation is target-agnostic in structure -- at each locus a fully
# homozygous line carries P1's or P2's allele with indicator g_j in
# {-1, +1}, and for a biparental population with no prior selection
# E[g_j] = 0, Var(g_j) = 1 for BOTH DH and RIL. Only the between-locus
# covariance Cov(g_j, g_l) = R_jl changes between targets: DH -> 1 - 2r;
# RIL (infinite selfing) -> (1 - 2r) / (1 + 2r) (RILs accumulate more
# recombination than a single DH meiosis, so linked loci covary less).
# The truncation-selection response, He^(i), and post-selection
# contributions therefore hold for both targets by using the target's R.
#
# Scope: target = "RIL" here means the standard infinite-selfing (fully
# homozygous) line, matching the package's default ril_mode = "infinite".
# Finite-selfing RIL residual heterozygosity (ril_mode = "finite") breaks
# the g_j in {-1, +1} assumption underlying this kernel and is out of scope
# for this phase; callers requesting target = "RIL" get the infinite-
# selfing kernel unconditionally (there is no ril_mode argument here yet).
#
# Reliability-weighted / Sigma_beta extensions (Phase 2): implemented below
# as the explicitly experimental reliability = "experimental_diagonal" branch of ng_ucpc_cross(), backed by
# ng_ucpc_locus_reliability(). This is an experimental shrinkage of the MERIT
# (uc) only. It does not propagate effect uncertainty through the truncation
# response or selected-fraction diversity calculation, so it must not be read
# as a calibrated reliability model. Only the heuristic per-marker factor
# r_j = beta_hat_j^2 / (beta_hat_j^2 +
# max(beta_var_j, 0)) discounts each locus's contribution to the merit sum.
# At beta_var = 0 (or absent), r_j = 1 for every locus and RA-UCPC reduces
# to Phase 1's UCPC exactly.
#
# The full off-diagonal Sigma_beta form and OCS-optimizer integration
# (Phase 3) are deliberately NOT implemented in this file.
# ---------------------------------------------------------------------------

# Internal experimental kernel; deliberately not exported until its selected-
# fraction approximation and reliability extension are externally calibrated.
ng_ucpc_locus_response <- function(d, beta, R) {
  d <- as.numeric(d)
  beta <- as.numeric(beta)
  R <- as.matrix(R)
  if (!length(d) || length(beta) != length(d) ||
      nrow(R) != length(d) || ncol(R) != length(d)) {
    ng_stop("d, beta, and the square response kernel R must have matching dimensions")
  }
  if (any(!is.finite(d)) || any(!is.finite(beta)) || any(!is.finite(R))) {
    ng_stop("d, beta, and R must contain only finite values")
  }
  if (max(abs(R - t(R))) > 1e-8) ng_stop("R must be symmetric")
  a <- d * beta
  as.numeric(R %*% a)
}

# Internal experimental diagnostic.
ng_ucpc_locus_reliability <- function(beta, beta_var) {
  nm <- names(beta)
  beta <- as.numeric(beta)
  if (any(!is.finite(beta))) ng_stop("beta must contain only finite values")
  if (is.null(beta_var)) {
    beta_var <- rep(0, length(beta))
  } else {
    beta_var <- as.numeric(beta_var)
  }
  if (length(beta_var) != length(beta)) {
    ng_stop("beta_var must be the same length as beta")
  }
  if (any(!is.finite(beta_var)) || any(beta_var < 0)) {
    ng_stop("beta_var must contain finite non-negative marginal variances")
  }
  beta2 <- beta^2
  denom <- beta2 + beta_var
  r <- ifelse(denom > 0, beta2 / denom, 1)
  r <- pmin(pmax(r, 0), 1)
  names(r) <- nm
  r
}

# Internal experimental UCPC scorer (see scope statement above).
ng_ucpc_cross <- function(geno, effects, marker_map, ids, pair,
                          prop = 0.1, h = 1,
                          model = c("haldane", "kosambi"),
                          target = c("DH", "RIL"),
                          reliability = c("none", "experimental_diagonal")) {
  model <- match.arg(model)
  target <- match.arg(target)
  reliability <- as.character(reliability[[1L]])
  if (identical(reliability, "diagonal")) {
    warning("reliability = 'diagonal' is deprecated; using experimental_diagonal",
            call. = FALSE)
    reliability <- "experimental_diagonal"
  }
  if (!(reliability %in% c("none", "experimental_diagonal"))) {
    ng_stop("reliability must be 'none' or 'experimental_diagonal'")
  }
  h <- suppressWarnings(as.numeric(h))
  if (length(h) != 1L || !is.finite(h) || h < 0 || h > 1) {
    ng_stop("h must be one finite selection accuracy in [0, 1]")
  }
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  pair <- as.character(pair)
  if (length(pair) != 2L) ng_stop("pair must be length-2 (parent1, parent2)")
  geno_ids <- ng_check_same_ids(geno, ids, "geno")
  row1 <- match(pair[[1L]], ids)
  row2 <- match(pair[[2L]], ids)
  if (is.na(row1) || is.na(row2)) ng_stop("pair ids must both be present in ids")

  marker_map <- ng_prepare_marker_map(marker_map, colnames(geno_ids), model = model)
  x1 <- geno_ids[row1, ]
  x2 <- geno_ids[row2, ]
  if (any(!is.finite(x1)) || any(!is.finite(x2))) {
    ng_stop("UCPC parent genotypes must be complete; impute or remove missing markers before scoring")
  }
  if (any(!(x1 %in% c(0, 2))) || any(!(x2 %in% c(0, 2)))) {
    ng_stop("UCPC currently requires fully inbred diploid parents coded exactly 0/2; ",
            "heterozygous or dosage-uncertain parents violate its origin-indicator derivation")
  }
  d <- 0.5 * (x1 - x2)
  if (is.null(effects$beta) || is.null(names(effects$beta)) || anyDuplicated(names(effects$beta))) {
    ng_stop("effects$beta must be a uniquely named marker-effect vector")
  }
  missing_beta <- setdiff(colnames(geno_ids), names(effects$beta))
  if (length(missing_beta)) {
    ng_stop("effects$beta is missing markers used by geno: ", paste(head(missing_beta, 5L), collapse = ", "))
  }
  beta <- suppressWarnings(as.numeric(effects$beta[colnames(geno_ids)]))
  names(beta) <- colnames(geno_ids)
  if (any(!is.finite(beta))) ng_stop("effects$beta must be finite for every genotype marker")
  intercept <- if (is.null(effects$intercept)) 0 else suppressWarnings(as.numeric(effects$intercept))
  if (length(intercept) != 1L || !is.finite(intercept)) {
    ng_stop("effects$intercept must be one finite number")
  }

  R <- ng_recomb_decay_matrix(marker_map, model = model, target = target)
  s <- ng_ucpc_locus_response(d, beta, R)
  a <- d * beta
  vpm <- sum(a * s)
  vpm_tol <- 1e-10 * max(1, sum(abs(a * s)))
  if (!is.finite(vpm) || vpm < -vpm_tol) {
    ng_stop("UCPC response kernel produced a negative variance; check marker map and inputs")
  }
  vpm <- max(vpm, 0)
  sigma <- sqrt(vpm)
  i <- ng_selection_intensity(prop)

  if (sigma > 0) {
    p_sel_raw <- 0.5 + (i * h * s) / (2 * sigma)
  } else {
    p_sel_raw <- rep(0.5, length(s))
  }
  p_sel <- pmin(pmax(p_sel_raw, 0), 1)
  names(p_sel) <- colnames(geno_ids)
  clipped_fraction <- mean(abs(p_sel - p_sel_raw) > 0)
  if (clipped_fraction > 0) {
    warning(sprintf(
      "UCPC normal-response approximation exceeded [0,1] at %.1f%% of markers; probabilities were clipped and the clipped fraction is reported",
      100 * clipped_fraction
    ), call. = FALSE)
  }

  c_p1 <- mean(p_sel)
  c_p2 <- 1 - c_p1
  # He^(i): expected heterozygosity of the selected fraction, counted only
  # over loci that segregate WITHIN this cross (d_j != 0); loci fixed within
  # the cross contribute 0, not the population-level 0.5.
  seg <- d != 0
  he_locus <- seg * 2 * p_sel * (1 - p_sel)
  he_sum <- sum(he_locus)
  he_i <- mean(he_locus)

  effects_mu <- effects
  effects_mu$intercept <- intercept
  gebv_pair <- ng_predict_gebv(geno_ids[c(row1, row2), , drop = FALSE], effects_mu)
  mu <- 0.5 * sum(gebv_pair)

  # The MERIT uc is the only quantity reliability affects. UCPC (Phase 1):
  # uc = mu + i*h*sigma (equivalently mu + i*h*sum(a*s)/sigma, since
  # sum(a*s) == vpm == sigma^2). RA-UCPC replaces the sum with a per-marker
  # reliability-weighted sum(r*a*s); r = 1 everywhere reduces exactly to
  # UCPC (the beta_var = 0 / absent equivalence).
  if (reliability == "experimental_diagonal") {
    beta_var <- effects$beta_var
    if (is.null(beta_var)) {
      beta_var <- rep(0, length(beta))
      names(beta_var) <- names(beta)
    } else {
      if (is.null(names(beta_var)) || anyDuplicated(names(beta_var))) {
        ng_stop("experimental reliability requires effects$beta_var to be uniquely named by marker")
      }
      missing_var <- setdiff(colnames(geno_ids), names(beta_var))
      if (length(missing_var)) ng_stop("effects$beta_var is missing one or more genotype markers")
      beta_var <- beta_var[colnames(geno_ids)]
    }
    r <- ng_ucpc_locus_reliability(beta, beta_var)
    uc <- if (sigma > 0) mu + i * h * sum(r * a * s) / sigma else mu
  } else {
    r <- rep(1, length(beta))
    names(r) <- names(beta)
    uc <- mu + i * h * sigma
  }

  list(
    c_p1 = c_p1,
    c_p2 = c_p2,
    p_sel = p_sel,
    p_sel_unclipped = p_sel_raw,
    p_sel_clipped_fraction = clipped_fraction,
    he_i = he_i,
    he_sum = he_sum,
    uc = uc,
    vpm = vpm,
    mu = mu,
    sigma = sigma,
    target = target,
    reliability = reliability,
    r = r
  )
}

# Internal plan-level selected-fraction summary.
ng_ucpc_selected_diversity <- function(crosses, parent_kinship, weight = NULL) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE)
  n <- nrow(crosses)
  if (n == 0L) ng_stop("crosses must have at least one row")
  required <- c("parent1", "parent2", "c_p1", "c_p2", "he_i")
  missing_cols <- setdiff(required, names(crosses))
  if (length(missing_cols)) ng_stop("crosses missing columns: ", paste(missing_cols, collapse = ", "))
  if (is.null(weight)) {
    weight <- if ("weight" %in% names(crosses)) as.numeric(crosses$weight) else rep(1 / n, n)
  }
  weight <- as.numeric(weight)
  if (length(weight) != n) ng_stop("weight must have length nrow(crosses)")
  if (any(!is.finite(weight)) || any(weight < 0)) {
    ng_stop("weight must contain finite non-negative crossing proportions")
  }
  w_total <- sum(weight)
  if (!is.finite(w_total) || w_total <= 0) ng_stop("weight must sum to a positive finite value")
  w <- weight / w_total

  parent1 <- as.character(crosses$parent1)
  parent2 <- as.character(crosses$parent2)
  if (anyNA(parent1) || anyNA(parent2) || any(!nzchar(parent1)) || any(!nzchar(parent2))) {
    ng_stop("cross parent IDs must be non-missing and non-blank")
  }
  cp1 <- suppressWarnings(as.numeric(crosses$c_p1))
  cp2 <- suppressWarnings(as.numeric(crosses$c_p2))
  he <- suppressWarnings(as.numeric(crosses$he_i))
  if (any(!is.finite(cp1)) || any(!is.finite(cp2)) ||
      any(cp1 < 0 | cp1 > 1) || any(cp2 < 0 | cp2 > 1) ||
      any(abs(cp1 + cp2 - 1) > 1e-8)) {
    ng_stop("c_p1 and c_p2 must be finite probabilities in [0,1] that sum to one per cross")
  }
  if (any(!is.finite(he)) || any(he < 0 | he > 0.5 + 1e-8)) {
    ng_stop("he_i must be finite and in [0, 0.5]")
  }
  parents <- unique(c(parent1, parent2))
  c_vec <- setNames(numeric(length(parents)), parents)
  for (r in seq_len(n)) {
    c_vec[parent1[r]] <- c_vec[parent1[r]] + w[r] * cp1[r]
    c_vec[parent2[r]] <- c_vec[parent2[r]] + w[r] * cp2[r]
  }
  c_total <- sum(c_vec)
  if (is.finite(c_total) && c_total > 0) c_vec <- c_vec / c_total

  K <- as.matrix(parent_kinship)
  if (!is.numeric(K) || nrow(K) != ncol(K) || is.null(rownames(K)) || is.null(colnames(K)) ||
      anyDuplicated(rownames(K)) || anyDuplicated(colnames(K)) ||
      !setequal(rownames(K), colnames(K))) {
    ng_stop("parent_kinship must be a uniquely named square additive-relationship matrix")
  }
  missing_parents <- setdiff(names(c_vec), intersect(rownames(K), colnames(K)))
  if (length(missing_parents)) ng_stop("parent_kinship is missing one or more cross parents")
  K <- K[names(c_vec), names(c_vec), drop = FALSE]
  if (any(!is.finite(K)) || max(abs(K - t(K))) > 1e-8) {
    ng_stop("parent_kinship must contain a finite symmetric additive-relationship matrix")
  }
  K <- (K + t(K)) / 2
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  tol <- 1e-8 * max(1, max(abs(ev)))
  if (min(ev) < -tol) ng_stop("parent_kinship must be positive semidefinite")
  group_relationship <- as.numeric(crossprod(c_vec, K %*% c_vec))
  group_coancestry <- group_relationship / 2
  He_selected <- sum(w * he)

  list(c_vec = c_vec,
       group_relationship = group_relationship,
       group_coancestry = group_coancestry,
       # Legacy UCPC name retained; this value is on additive-relationship,
       # not kinship/coancestry, scale.
       D_selected = group_relationship,
       He_selected = He_selected)
}
