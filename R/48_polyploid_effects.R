# Additive + DOMINANCE marker-effect estimation and value prediction for polyploid dosage data.
#
# Outbred / clonally-propagated crops (cassava, potato, ...) select on TOTAL GENOTYPIC value
# (breeding value + dominance deviation), because dominance is transmitted intact when a clone is
# propagated; and mate performance depends on dominance (mid-parent heterosis). The additive-only
# ridge does not account for this. Here we fit both components jointly on designs constructed to be
# ORTHOGONAL under Hardy-Weinberg, so the estimates partition cleanly:
#   additive design   W = M - ploidy * p                        (average-effect basis)
#   dominance design  D = (H - E[H]) - b (M - ploidy p),        (dominance-deviation basis)
#                     H = d(ploidy - d),  b = (ploidy - 1)(1 - 2p)
# breeding value  BV = intercept + W beta_add ;  genotypic value  GV = BV + D beta_dom.
#
# The b term is NOT cosmetic. A plain mean-centered H is strongly correlated with W:
#   Cov(M, H) = ploidy (ploidy - 1) p q (1 - 2p)
# which vanishes ONLY at p = 0.5. At a low-MAF marker the two columns run to r ~ 0.96 (measured,
# ploidy 4, p = 0.1), so ridge cannot separate them and the additive/dominance split for those
# markers is decided by the penalty rather than by the data. Regressing H on W removes exactly that
# dependence.
#
# b is the OBSERVED per-marker regression coefficient Cov(W_k, H_k)/Var(W_k) in the training
# genotypes, not the Hardy-Weinberg value (ploidy - 1)(1 - 2p). This is the "statistical"
# parameterization (Alvarez-Castro & Carlborg 2007; Vitezica et al. 2013): breeding populations are
# selected, related and frequently out of HWE, and using the observed coefficient makes the two
# designs EXACTLY orthogonal in the data actually being fitted rather than only in expectation.
# Under HWE it converges to (ploidy - 1)(1 - 2p), and at ploidy 2 that limit reproduces the standard
# orthogonal dominance coding (-2p^2, 2pq, -2q^2) of Vitezica et al. (2013, Genetics 195:1223).
#
# The change is an exact reparameterization of the SAME model space --
# a W + d (H - E[H]) == (a + d b) W + d D -- so an unpenalized fit is invariant and only the
# ridge-implied A/D split moves. It is the split that was previously arbitrary, which is the point.

# Fit additive (and optionally digenic dominance) marker effects from dosage + a phenotype.
ng_polyploid_fit_effects <- function(dosage,
                                     y,
                                     ploidy = 2L,
                                     model = c("additive_dominance", "additive"),
                                     min_maf = 0,
                                     impute_missing = TRUE,
                                     seed = 1L,
                                     allow_experimental_dominance = FALSE) {
  model <- match.arg(model)
  if (identical(model, "additive_dominance") && !isTRUE(allow_experimental_dominance)) {
    ng_stop("polyploid additive+dominance fitting is experimental because it currently uses one ",
            "ridge penalty for both variance components; set allow_experimental_dominance = TRUE ",
            "only for research diagnostics")
  }
  if (identical(model, "additive_dominance")) {
    warning("experimental polyploid dominance fit: additive and dominance components share one ridge penalty",
            call. = FALSE)
  }
  prep <- ng_polyploid_prep_dosage(dosage, ploidy, min_maf = min_maf, impute_missing = impute_missing)
  M <- prep$M; ids <- prep$ids; p <- prep$p; ploidy <- prep$ploidy
  markers <- colnames(M)
  yv <- suppressWarnings(as.numeric(y))
  names(yv) <- if (!is.null(names(y))) names(y) else ids
  yv <- yv[ids]
  if (all(!is.finite(yv))) ng_stop("phenotype y has no finite values for the dosage samples")

  W <- sweep(M, 2L, ploidy * p, "-")                 # additive (freq-centered) design
  hbar <- NULL; b_orth <- NULL; m_add <- ncol(W)
  if (identical(model, "additive_dominance")) {
    H <- M * (ploidy - M)
    hbar <- vapply(seq_len(ncol(H)), function(j) {
      observed <- !prep$missing[, j]
      mean(H[observed, j])
    }, numeric(1L))
    # Observed per-marker regression of H on W (statistical parameterization; see header).
    # W is mean-centered by construction (p = colMeans(M)/ploidy), so b is a ratio of sums.
    # Monomorphic / invariant markers regress to 0, leaving D as plain mean-centered H.
    Hc <- sweep(H, 2L, hbar, "-")
    Wc <- sweep(W, 2L, colMeans(W), "-")
    ss <- colSums(Wc * Wc)
    b_orth <- ifelse(ss > 0, colSums(Wc * Hc) / ifelse(ss > 0, ss, 1), 0)
    b_orth[!is.finite(b_orth)] <- 0
    D <- Hc - sweep(W, 2L, b_orth, "*")              # orthogonal dominance design
    # Mean-imputed additive dosage is W=0, but its nonlinear H value is not the
    # conditional expectation of the centered dominance covariate. A missing
    # call carries no marker-level dominance information and must contribute 0.
    if (any(prep$missing)) D[prep$missing] <- 0
    X <- cbind(W, D)
  } else {
    X <- W
  }
  colnames(X) <- paste0("f", seq_len(ncol(X))); rownames(X) <- ids
  fit <- ng_fit_ridge_effects(X, yv, seed = as.integer(seed))
  beta <- as.numeric(fit$beta)
  structure(list(
    beta_add = beta[seq_len(m_add)],
    beta_dom = if (identical(model, "additive_dominance")) beta[(m_add + 1L):(2L * m_add)] else NULL,
    intercept = fit$intercept, allele_freq = p, hbar = hbar, b_orth = b_orth, markers = markers,
    ploidy = ploidy, model = model,
    cv_predictive_r2 = fit$cv_predictive_r2,
    reliability = NA_real_, reliability_is_calibrated = FALSE
  ), class = "ng_polyploid_effects")
}

# Predict breeding value (additive) or genotypic value (additive + dominance) for dosage samples.
# For CLONAL selection use type = "genotypic"; for expected progeny/breeding gain use "breeding".
ng_polyploid_predict_value <- function(fit, dosage, type = c("genotypic", "breeding")) {
  if (!inherits(fit, "ng_polyploid_effects")) ng_stop("fit must come from ng_polyploid_fit_effects")
  type <- match.arg(type)
  M <- ng_polyploid_as_dosage_matrix(dosage, ploidy = fit$ploidy, name = "dosage")
  miss <- setdiff(fit$markers, colnames(M))
  if (length(miss)) ng_stop("dosage is missing ", length(miss), " markers the model was fit on")
  M <- M[, fit$markers, drop = FALSE]
  W <- sweep(M, 2L, fit$ploidy * fit$allele_freq, "-")
  bv <- as.numeric(fit$intercept + W %*% fit$beta_add)
  names(bv) <- rownames(M)
  if (identical(type, "breeding") || is.null(fit$beta_dom)) return(bv)
  H <- M * (fit$ploidy - M)
  # Same orthogonal dominance basis the effects were fitted on (see header).
  D <- sweep(H, 2L, fit$hbar, "-") - sweep(W, 2L, fit$b_orth, "*")
  gv <- as.numeric(bv + D %*% fit$beta_dom); names(gv) <- rownames(M)
  gv
}
