# Additive + DOMINANCE marker-effect estimation and value prediction for polyploid dosage data.
#
# Outbred / clonally-propagated crops (cassava, potato, ...) select on TOTAL GENOTYPIC value
# (breeding value + dominance deviation), because dominance is transmitted intact when a clone is
# propagated; and mate performance depends on dominance (mid-parent heterosis). The additive-only
# ridge does not account for this. Here we fit both components jointly on frequency-centered
# designs that are orthogonal under random mating, so the estimates partition cleanly:
#   additive design   W = M - ploidy * p            (breeding-value / average-effect basis)
#   dominance design  D = H - E[H], H = d*(ploidy-d) (digenic dominance-deviation basis)
# breeding value  BV = intercept + W beta_add ;  genotypic value  GV = BV + D beta_dom.

# Fit additive (and optionally digenic dominance) marker effects from dosage + a phenotype.
ng_polyploid_fit_effects <- function(dosage,
                                     y,
                                     ploidy = 2L,
                                     model = c("additive_dominance", "additive"),
                                     min_maf = 0,
                                     impute_missing = TRUE,
                                     seed = 1L) {
  model <- match.arg(model)
  prep <- ng_polyploid_prep_dosage(dosage, ploidy, min_maf = min_maf, impute_missing = impute_missing)
  M <- prep$M; ids <- prep$ids; p <- prep$p; ploidy <- prep$ploidy
  markers <- colnames(M)
  yv <- suppressWarnings(as.numeric(y))
  names(yv) <- if (!is.null(names(y))) names(y) else ids
  yv <- yv[ids]
  if (all(!is.finite(yv))) ng_stop("phenotype y has no finite values for the dosage samples")

  W <- sweep(M, 2L, ploidy * p, "-")                 # additive (freq-centered) design
  hbar <- NULL; m_add <- ncol(W)
  if (identical(model, "additive_dominance")) {
    H <- M * (ploidy - M); hbar <- colMeans(H)
    D <- sweep(H, 2L, hbar, "-")                     # dominance (mean-centered) design
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
    intercept = fit$intercept, allele_freq = p, hbar = hbar, markers = markers,
    ploidy = ploidy, model = model, reliability = fit$reliability
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
  D <- sweep(H, 2L, fit$hbar, "-")
  gv <- as.numeric(bv + D %*% fit$beta_dom); names(gv) <- rownames(M)
  gv
}
