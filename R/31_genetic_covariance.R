# Estimators for the additive genetic covariance matrix G and phenotypic
# covariance matrix P across traits. Both are needed by formal Smith-Hazel and
# Pesek-Baker index calculations. `method = "auto"` requires the joint
# multivariate REML implementation in sommer. The separate-trait ridge method
# below is retained only as an explicitly requested diagnostic heuristic.
#
# Statistical justification (two_stage_ridge):
#   For trait t, fit ridge marker effects beta_hat_t on the same n x m geno used
#   downstream. With centered genotypes X_c (mean 2*p subtracted) and
#   K = X_c X_c' / denom where denom = sum_k 2 p_k (1 - p_k) (VanRaden), the
#   additive-genetic covariance of breeding values g = X_c beta is
#       Var(g) = X_c X_c' * Var(beta) = K * denom * Var(beta_k)  (iid markers).
#
#   The diagnostic heuristic constructs G_hat from two plug-in parts:
#
#   (a) Diagonals via GBLUP lambda-inversion (method-of-moments). Under the
#       GBLUP model the optimal ridge penalty is
#           lambda_t = sigma_e^2_t / (sigma_g^2_t / denom),
#       so the genetic variance is identifiable from the chosen lambda and the
#       fitted residual variance as
#           sigma_g^2_t = sigma_e^2_t * denom / lambda_t.
#       This identity is model-dependent and lambda is selected by prediction,
#       not estimated jointly as a variance-component ratio. It is therefore
#       not a REML estimate.
#
#   (b) Off-diagonals via the sample correlation of the trait marker-effect
#       vectors. Correlated residuals and trait-specific shrinkage can both
#       contaminate this correlation, so it must not be treated as a formal
#       genetic correlation.
#
#   G_hat = D %*% R_beta %*% D with D = sqrt(diag(sigma_g^2)). The result is
#   on the same scale as a GBLUP variance component fit against the VanRaden
#   K only under the working ridge assumptions above.
#
#   The G is then projected to the nearest positive-semidefinite matrix using
#   Matrix::nearPD() if available, or a documented eigen-clip fallback.
#
# PSD projection:
#   `Matrix::nearPD()` if installed; otherwise a simple eigen-clip that floors
#   negative eigenvalues at `eps_eig * mean(positive eigenvalues)`. nearPD is
#   preferred because it gives the closest-in-Frobenius PSD matrix; the
#   eigen-clip is a documented fallback so the function works without Matrix.

# --- helpers ----------------------------------------------------------------

ng_genetic_cov_eigen_clip <- function(M, eps_eig = 1e-6) {
  M <- (M + t(M)) / 2
  ev <- eigen(M, symmetric = TRUE)
  vals <- ev$values
  pos <- vals[vals > 0]
  floor_val <- if (length(pos)) eps_eig * mean(pos) else eps_eig
  vals[vals < floor_val] <- floor_val
  ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
}

ng_genetic_cov_project_psd <- function(M, eps_eig = 1e-6) {
  M <- (M + t(M)) / 2
  if (requireNamespace("Matrix", quietly = TRUE)) {
    pd <- tryCatch(Matrix::nearPD(M, base.matrix = TRUE, keepDiag = FALSE,
                                  ensureSymmetry = TRUE),
                   error = function(e) NULL)
    if (!is.null(pd)) return(as.matrix(pd$mat))
  }
  ng_genetic_cov_eigen_clip(M, eps_eig = eps_eig)
}

ng_genetic_cov_to_correlation <- function(M) {
  d <- sqrt(pmax(diag(M), 0))
  d[!is.finite(d) | d <= 0] <- 1
  R <- M / outer(d, d)
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  R
}

# --- two-stage ridge engine -------------------------------------------------

ng_genetic_cov_two_stage_ridge <- function(geno, Y, ridge_lambda = NULL,
                                            kfold = 5L, seed = 1L,
                                            return_diagnostics = FALSE) {
  m <- ncol(geno)
  t <- ncol(Y)
  trait_names <- colnames(Y)

  p <- colMeans(geno, na.rm = TRUE) / 2
  denom <- sum(2 * p * (1 - p), na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) denom <- m

  Beta <- matrix(0, nrow = m, ncol = t)
  colnames(Beta) <- trait_names
  rownames(Beta) <- colnames(geno)
  reliabilities <- rep(NA_real_, t)
  cv_predictive_r2 <- rep(NA_real_, t)
  lambdas <- rep(NA_real_, t)
  sigma_e2s <- rep(NA_real_, t)
  sigma_g2s <- rep(NA_real_, t)
  names(reliabilities) <- names(cv_predictive_r2) <- names(lambdas) <-
    names(sigma_e2s) <- names(sigma_g2s) <- trait_names

  for (j in seq_len(t)) {
    y_j <- Y[, j]
    # SEED (0.29.0): the base seed itself, NOT `seed + j`.
    #
    # This seed reaches ng_choose_ridge_lambda() -> set.seed(seed); sample(...),
    # i.e. it picks the k-fold CV partition behind ridge-lambda selection. Keying
    # it to `j` -- the trait's COLUMN POSITION in Y -- made every element of G_hat
    # depend on the order the caller happened to bind the phenotype columns in: a
    # different split can select a different lambda, which changes beta_j, which
    # changes both sigma_g2_j (via the lambda inversion below) and the Pearson
    # correlation of the beta vectors that supplies every off-diagonal.
    #
    # The fix matches ng_cp__stage_predict() (R/39, 0.28.0): all traits share ONE
    # fold partition. Folds are a nuisance parameter of lambda selection, not a
    # source of innovation -- given lambda, beta_j is a deterministic function of
    # y_j alone -- so a shared split cannot couple the traits, and it makes the
    # per-trait cv_predictive_r2 comparisons paired.
    fit_j <- ng_fit_ridge_effects(
      geno = geno,
      y = setNames(y_j, rownames(geno)),
      ids = rownames(geno),
      lambda = ridge_lambda,
      kfold = kfold,
      seed = seed
    )
    Beta[, j] <- fit_j$beta
    cv_predictive_r2[[j]] <- fit_j$cv_predictive_r2
    lambdas[[j]] <- fit_j$lambda
    sigma_e2s[[j]] <- fit_j$sigma_e2
    # GBLUP lambda-inversion: sigma_g^2 = sigma_e^2 * denom / lambda. Floor at
    # a tiny positive value so the projection step still produces a valid PSD G.
    sigma_g2_j <- fit_j$sigma_e2 * denom / max(fit_j$lambda, 1e-12)
    if (!is.finite(sigma_g2_j) || sigma_g2_j <= 0) sigma_g2_j <- 1e-8
    sigma_g2s[[j]] <- sigma_g2_j
  }

  # Off-diagonals: sample Pearson correlation of the m-vectors of beta_hats.
  # Pearson correlation is invariant to per-trait multiplicative shrinkage, so
  # this is robust to the ridge attenuation that wrecks naive inner-product
  # plug-ins under small n.
  if (m > 1L) {
    R_beta <- suppressWarnings(stats::cor(Beta))
    R_beta[!is.finite(R_beta)] <- 0
    diag(R_beta) <- 1
  } else {
    R_beta <- diag(t)
  }

  D <- sqrt(sigma_g2s)
  G_raw <- diag(D, t) %*% R_beta %*% diag(D, t)
  G_raw <- (G_raw + t(G_raw)) / 2
  dimnames(G_raw) <- list(trait_names, trait_names)

  G_hat <- ng_genetic_cov_project_psd(G_raw)
  dimnames(G_hat) <- list(trait_names, trait_names)

  out <- list(G_hat = G_hat,
              reliabilities = reliabilities,
              cv_predictive_r2 = cv_predictive_r2,
              lambdas = lambdas,
              sigma_e2 = sigma_e2s,
              sigma_g2 = sigma_g2s,
              R_beta = R_beta,
              denom = denom,
              G_raw = G_raw,
              beta_hat = Beta)
  if (!return_diagnostics) {
    out$beta_hat <- NULL
    out$R_beta <- NULL
    out$G_raw <- NULL
  }
  out
}

# --- sommer REML engine -----------------------------------------------------

ng_genetic_cov_sommer_remml <- function(geno, Y) {
  if (!requireNamespace("sommer", quietly = TRUE)) {
    return(list(error = "sommer not installed"))
  }
  K <- ng_parent_kinship(geno)
  ids <- rownames(geno)
  trait_names <- colnames(Y)
  df <- as.data.frame(Y)
  df$id <- factor(ids, levels = ids)
  rownames(K) <- colnames(K) <- ids
  rhs <- paste(trait_names, collapse = ", ")
  form <- stats::as.formula(paste0("cbind(", rhs, ") ~ 1"))
  fit <- tryCatch({
    sommer::mmer(
      fixed = form,
      random = stats::as.formula("~ sommer::vsr(id, Gu = K, Gtc = sommer::unsm(length(trait_names)))"),
      rcov = stats::as.formula("~ sommer::vsr(units, Gtc = sommer::unsm(length(trait_names)))"),
      data = df,
      verbose = FALSE,
      tolParInv = 1e-3,
      date.warning = FALSE
    )
  }, error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(error = conditionMessage(fit)))
  }
  # sommer stores variance components in fit$sigma; the first element corresponds
  # to the random GRM term.
  G_hat <- tryCatch({
    sigma <- fit$sigma
    g_term <- sigma[[1L]]
    g_term <- as.matrix(g_term)
    g_term
  }, error = function(e) NULL)
  if (is.null(G_hat) || !all(dim(G_hat) == c(length(trait_names), length(trait_names)))) {
    return(list(error = "could not extract G from sommer fit"))
  }
  dimnames(G_hat) <- list(trait_names, trait_names)
  G_hat <- (G_hat + t(G_hat)) / 2
  # Residual (units) covariance: needed to imply per-trait h2 for the
  # already-shrunk-input guard in ng_estimate_genetic_covariance(). sommer puts
  # the rcov term last in fit$sigma. Missing/odd shapes are not fatal -- the
  # guard simply reports NA heritability for that engine.
  R_hat <- tryCatch({
    sigma <- fit$sigma
    r <- as.matrix(sigma[[length(sigma)]])
    if (all(dim(r) == c(length(trait_names), length(trait_names)))) {
      dimnames(r) <- list(trait_names, trait_names)
      (r + t(r)) / 2
    } else {
      NULL
    }
  }, error = function(e) NULL)
  list(G_hat = ng_genetic_cov_project_psd(G_hat), R_hat = R_hat)
}

# --- already-shrunk-input guard ---------------------------------------------
#
# Per-trait REML heritability of Y against the GRM, computed ONLY as an input
# sanity check (it is never used to build G-hat). This is the EMMA-style
# eigen-decomposition profile REML for
#     y = 1 mu + g + e,   Var(y) = sigma_g2 K + sigma_e2 I,
# profiled over delta = sigma_e2 / sigma_g2 on one eigen(K). It is deliberately
# engine-independent: the two_stage_ridge sigma_e2 is an in-sample ridge
# residual at a CV-selected lambda, which does NOT collapse to zero on shrunk
# input, so it cannot serve as the trigger.
#
# The guard is skipped (all NA, with a stated reason) for very small n, where
# the estimate is too erratic to act on, and for very large n, where the n x n
# eigen decomposition would stop being cheap.
ng_genetic_cov_grm_reml_h2 <- function(geno, Y, min_n = 100L, max_n = 2000L) {
  trait_names <- colnames(Y)
  n_traits <- ncol(Y)
  n <- nrow(Y)
  na_vec <- function() stats::setNames(rep(NA_real_, n_traits), trait_names)
  out <- list(h2 = na_vec(), sigma_g2 = na_vec(), sigma_e2 = na_vec(),
              method = NA_character_, note = NA_character_)
  if (n < min_n) {
    out$note <- paste0("skipped: n = ", n, " < ", min_n,
                       " complete-case rows (REML heritability too erratic to act on)")
    return(out)
  }
  if (n > max_n) {
    out$note <- paste0("skipped: n = ", n, " > ", max_n,
                       " (n x n eigen decomposition is no longer a cheap guard)")
    return(out)
  }
  K <- tryCatch(ng_parent_kinship(geno), error = function(e) NULL)
  ev <- if (is.null(K)) NULL else tryCatch(eigen(K, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) {
    out$note <- "skipped: the GRM eigen decomposition failed"
    return(out)
  }
  U <- ev$vectors
  d <- pmax(ev$values, 1e-10)
  k_bar <- mean(d)
  U1 <- crossprod(U, rep(1, n))
  for (j in seq_len(n_traits)) {
    Uy <- crossprod(U, as.numeric(Y[, j]))
    neg_reml <- function(log_delta) {
      delta <- exp(log_delta)
      w <- 1 / (d + delta)
      A <- sum(U1 * w * U1)
      mu <- sum(U1 * w * Uy) / A
      r <- Uy - mu * U1
      s2 <- sum(w * r * r) / (n - 1)
      if (!is.finite(s2) || s2 <= 0) return(Inf)
      0.5 * ((n - 1) * log(s2) + sum(log(d + delta)) + log(A))
    }
    opt <- tryCatch(stats::optimize(neg_reml, c(log(1e-9), log(1e9))), error = function(e) NULL)
    if (is.null(opt)) next
    delta <- exp(opt$minimum)
    w <- 1 / (d + delta)
    A <- sum(U1 * w * U1)
    mu <- sum(U1 * w * Uy) / A
    r <- Uy - mu * U1
    s2 <- sum(w * r * r) / (n - 1)
    out$sigma_g2[[j]] <- s2 * k_bar
    out$sigma_e2[[j]] <- s2 * delta
    out$h2[[j]] <- k_bar / (k_bar + delta)
  }
  out$method <- "grm_profile_reml"
  out
}

# Cheap, non-invasive check for Y that is not phenotypes/BLUEs. Both BLUPs and
# GEBVs have had the residual removed before they reach this function, so a
# variance-component fit sees (almost) no residual and the implied heritability
# collapses onto 1. That is the signature this looks for. It WARNS and never
# errors: a user may knowingly supply near-h2 = 1 data (e.g. a simulated
# noiseless trait) and still want the covariance structure.
ng_genetic_cov_shrunk_input_guard <- function(h2, residual_variance,
                                              observed_variance, trait_names,
                                              h2_tol = 0.99, resid_frac_tol = 0.01) {
  resid_frac <- residual_variance / observed_variance
  names(h2) <- names(resid_frac) <- trait_names
  flagged <- (is.finite(h2) & h2 >= h2_tol) |
    (is.finite(resid_frac) & resid_frac <= resid_frac_tol)
  flagged[is.na(flagged)] <- FALSE
  if (any(flagged)) {
    bad <- trait_names[flagged]
    warning(
      "ng_estimate_genetic_covariance(): Y appears to contain ALREADY-SHRUNK predictions ",
      "(BLUPs or GEBVs) rather than phenotypes or BLUEs for trait(s): ",
      paste(bad, collapse = ", "), ". ",
      sprintf("Implied heritability %s and residual/observed variance ratio %s. ",
              paste(sprintf("%s=%.4f", bad, h2[flagged]), collapse = ", "),
              paste(sprintf("%s=%.2e", bad, resid_frac[flagged]), collapse = ", ")),
      "Var(BLUP) = sigma2_g - PEV, so REML sees too little genetic variance and almost no ",
      "residual: sigma2_g is biased LOW and h2 is pushed toward 1. GEBVs are worse and ",
      "circular -- they are a linear function of the same markers the GRM is built from, so ",
      "G-hat becomes the covariance of PREDICTIONS, understating sigma2_g by roughly the ",
      "reliability. Remedies: deregress BLUPs before use (Garrick-Taylor-Dekkers 2009), or ",
      "correct genetic correlations by dividing by sqrt(r_t * r_s), or supply a validated G ",
      "directly. Not an error: G-hat is returned, and is biased low.",
      call. = FALSE)
  }
  list(implied_heritability = h2, residual_fraction = resid_frac,
       suspected = any(flagged), traits = trait_names[flagged])
}

# --- public estimator -------------------------------------------------------

# Estimate additive genetic covariance G across traits from a training set.
#
# `method = "auto"` requires the multivariate REML engine. The two-stage ridge
# estimator is retained only as an explicit diagnostic heuristic; correlated
# residuals can induce correlated univariate marker-effect estimates and hence
# masquerade as genetic covariance.
#
# =====================================================================
# WHAT `Y` MUST CONTAIN  (read this before supplying anything but raw
# phenotypes -- BLUEs, BLUPs and GEBVs are NOT interchangeable here)
# =====================================================================
#
# `Y` is an n x t matrix of PHENOTYPIC observations, one column per trait, with
# row names matching `geno`. Both engines below are variance-component fits:
# they need input that still carries residual variation, because that is exactly
# what separates sigma2_g from sigma2_e.
#
#   BLUEs  -- VALID. This is the standard two-stage genomic analysis: stage-one
#             BLUEs enter stage two as the response, and REML splits the genetic
#             variance from the BLUE estimation error. On unbalanced data the
#             BLUEs ideally enter WEIGHTED by their inverse squared standard
#             errors; this function fits them unweighted, so treat unbalanced
#             stage-one designs as an approximation.
#
#   BLUPs  -- INVALID AS SUPPLIED. A BLUP is already shrunk toward the mean:
#             Var(BLUP) = sigma2_g - PEV. REML therefore sees too little genetic
#             variance AND almost no residual, so sigma2_g is biased LOW while
#             the implied h2 is pushed toward 1. Genetic CORRELATIONS are also
#             distorted whenever reliabilities differ across traits, because
#             each column is shrunk by a different factor.
#             Remedy: DEREGRESS the BLUPs before use (Garrick, Taylor & Dekkers
#             2009, Genet. Sel. Evol. 41:55), or correct the estimated genetic
#             correlation between traits t and s by dividing by sqrt(r_t * r_s)
#             with r the respective reliabilities, or supply a validated G to
#             the index functions directly and skip this estimator.
#
#   GEBVs  -- WORST, AND CIRCULAR. GEBVs are a linear function of the very
#             marker matrix the GRM is built from, so they lie in its span. The
#             residual goes to ~0, h2 -> 1, and G-hat becomes the covariance of
#             the PREDICTIONS rather than of the true breeding values --
#             understating sigma2_g by roughly the reliability. Do not use.
#
# A cheap guard below WARNS (never errors) when the REML residual variance
# against the GRM is at or near zero, or the implied h2 is at or near 1, because
# that is the signature of already-shrunk input. It is an ENGINE-INDEPENDENT
# profile REML computed only for this check and never used to build G-hat; see
# ng_genetic_cov_grm_reml_h2() and ng_genetic_cov_shrunk_input_guard(). It is
# skipped, with a stated reason in `implied_heritability_note`, for n < 100
# (too erratic) and n > 2000 (no longer cheap).
#
# Provenance attributes on the returned G-hat (for a caller to report):
#   method, requested_method, formal_variance_component_estimate,
#   genetic_correlation, n_used, n_markers, traits, y_input_contract,
#   genetic_variance, residual_variance          (from the fitting engine),
#   implied_heritability, implied_heritability_method, implied_heritability_note,
#   reml_genetic_variance, reml_residual_variance (from the guard),
#   shrunk_input_suspected, shrunk_input_traits
#   (and `diagnostics` when return_diagnostics = TRUE).
ng_estimate_genetic_covariance <- function(geno,
                                           Y,
                                           method = c("auto", "two_stage_ridge", "sommer_remml"),
                                           ridge_lambda = NULL,
                                           kfold = 5L,
                                           seed = 1L,
                                           return_diagnostics = FALSE) {
  method <- match.arg(method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  if (is.null(rownames(geno))) ng_stop("geno must have row names")
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  if (is.null(rownames(Y))) {
    if (nrow(Y) != nrow(geno)) ng_stop("Y must have row names or have nrow(Y) == nrow(geno)")
    rownames(Y) <- rownames(geno)
  }
  ids <- intersect(rownames(geno), rownames(Y))
  if (!length(ids)) ng_stop("no overlapping ids between geno and Y")
  geno <- ng_check_same_ids(geno, ids, "geno")
  Y <- Y[ids, , drop = FALSE]
  ok <- stats::complete.cases(Y)
  if (sum(ok) < 5L) ng_stop("need at least 5 complete-case rows to estimate G")
  geno <- geno[ok, , drop = FALSE]
  Y <- Y[ok, , drop = FALSE]
  storage.mode(Y) <- "double"

  resolved_method <- method
  if (method == "auto") {
    if (!requireNamespace("sommer", quietly = TRUE)) {
      ng_stop("method = 'auto' requires sommer for multivariate REML; install sommer, ",
              "supply a validated G, or request two_stage_ridge explicitly as a heuristic")
    }
    resolved_method <- "sommer_remml"
  }

  diag_list <- list()
  n_traits <- ncol(Y)
  genetic_variance <- stats::setNames(rep(NA_real_, n_traits), colnames(Y))
  residual_variance <- stats::setNames(rep(NA_real_, n_traits), colnames(Y))
  observed_variance <- stats::setNames(
    apply(Y, 2L, function(col) stats::var(col, na.rm = TRUE)), colnames(Y))
  if (resolved_method == "sommer_remml") {
    sommer_out <- ng_genetic_cov_sommer_remml(geno, Y)
    if (!is.null(sommer_out$error)) {
      ng_stop("sommer multivariate REML failed: ", sommer_out$error,
              ". No heuristic fallback was used.")
    } else {
      G_hat <- sommer_out$G_hat
      diag_list$engine <- "sommer_remml"
      genetic_variance[] <- diag(as.matrix(G_hat))
      if (!is.null(sommer_out$R_hat)) residual_variance[] <- diag(sommer_out$R_hat)
      diag_list$residual_covariance <- sommer_out$R_hat
    }
  }
  if (resolved_method == "two_stage_ridge") {
    warning("two_stage_ridge is a heuristic diagnostic, not a multivariate variance-component estimator; do not use it as formal G when residual traits may be correlated",
            call. = FALSE)
    ts <- ng_genetic_cov_two_stage_ridge(geno, Y, ridge_lambda = ridge_lambda,
                                          kfold = kfold, seed = seed,
                                          return_diagnostics = return_diagnostics)
    G_hat <- ts$G_hat
    diag_list$engine <- "two_stage_ridge"
    diag_list$reliabilities <- ts$reliabilities
    diag_list$cv_predictive_r2 <- ts$cv_predictive_r2
    diag_list$lambdas <- ts$lambdas
    diag_list$sigma_e2 <- ts$sigma_e2
    diag_list$sigma_g2 <- ts$sigma_g2
    diag_list$denom <- ts$denom
    genetic_variance[] <- as.numeric(ts$sigma_g2)
    residual_variance[] <- as.numeric(ts$sigma_e2)
    if (return_diagnostics) {
      diag_list$beta_hat <- ts$beta_hat
      diag_list$R_beta <- ts$R_beta
      diag_list$G_raw <- ts$G_raw
    }
  }

  dimnames(G_hat) <- list(colnames(Y), colnames(Y))
  # Guard: warn (never error) when Y looks like already-shrunk predictions.
  reml_h2 <- ng_genetic_cov_grm_reml_h2(geno, Y)
  guard <- ng_genetic_cov_shrunk_input_guard(
    h2 = reml_h2$h2,
    residual_variance = reml_h2$sigma_e2,
    observed_variance = observed_variance,
    trait_names = colnames(Y))
  attr(G_hat, "genetic_correlation") <- ng_genetic_cov_to_correlation(G_hat)
  attr(G_hat, "n_used") <- nrow(Y)
  attr(G_hat, "method") <- diag_list$engine
  attr(G_hat, "requested_method") <- method
  attr(G_hat, "formal_variance_component_estimate") <- identical(diag_list$engine, "sommer_remml")
  # Provenance a caller needs to report what this G-hat actually is.
  attr(G_hat, "traits") <- colnames(Y)
  attr(G_hat, "n_markers") <- ncol(geno)
  attr(G_hat, "y_input_contract") <-
    "Y must be phenotypes or BLUEs; BLUPs must be deregressed first and GEBVs must not be used"
  attr(G_hat, "genetic_variance") <- genetic_variance
  attr(G_hat, "residual_variance") <- residual_variance
  # Guard quantities (independent of the engine that built G-hat).
  attr(G_hat, "implied_heritability") <- guard$implied_heritability
  attr(G_hat, "implied_heritability_method") <- reml_h2$method
  attr(G_hat, "implied_heritability_note") <- reml_h2$note
  attr(G_hat, "reml_genetic_variance") <- reml_h2$sigma_g2
  attr(G_hat, "reml_residual_variance") <- reml_h2$sigma_e2
  attr(G_hat, "shrunk_input_suspected") <- guard$suspected
  attr(G_hat, "shrunk_input_traits") <- guard$traits
  if (return_diagnostics) attr(G_hat, "diagnostics") <- diag_list
  G_hat
}

# --- phenotypic covariance --------------------------------------------------

# Ledoit-Wolf (2004) analytical optimal shrinkage of the sample covariance
# matrix toward a diagonal target. The closed-form intensity minimizes the
# expected Frobenius distance to the population covariance. Implemented from
# the original formulas (Ledoit & Wolf 2004, J. Multivariate Anal., eq. (14))
# so we do not need `corpcor`.
ng_genetic_cov_ledoit_wolf <- function(Y) {
  Y <- as.matrix(Y)
  Y <- Y[stats::complete.cases(Y), , drop = FALSE]
  n <- nrow(Y)
  p <- ncol(Y)
  if (n < 2L) return(list(S_star = diag(p), intensity = 1, S = diag(p), T = diag(p)))
  mu <- colMeans(Y)
  Yc <- sweep(Y, 2L, mu, "-")
  # Sample covariance (1/n, not 1/(n-1), to match Ledoit-Wolf's derivation).
  S <- crossprod(Yc) / n
  # Diagonal target T (preserves the sample variances; assumes zero off-diag).
  T_target <- diag(diag(S), p)
  # pi-hat: sum_{i,j} Var(S_{ij}).
  Yc2 <- Yc * Yc
  pi_mat <- crossprod(Yc2) / n - S * S
  pi_hat <- sum(pi_mat)
  # rho-hat: contribution from the diagonal entries of T (off-diagonals of T
  # are zero, so only diagonal terms contribute).
  rho_hat <- sum(diag(pi_mat))
  # gamma-hat: ||S - T||_F^2.
  diff_ST <- S - T_target
  gamma_hat <- sum(diff_ST * diff_ST)
  # Optimal intensity in (0, 1).
  if (!is.finite(gamma_hat) || gamma_hat <= 0) {
    intensity <- 1
  } else {
    kappa <- (pi_hat - rho_hat) / gamma_hat
    intensity <- max(0, min(1, kappa / n))
  }
  S_star <- intensity * T_target + (1 - intensity) * S
  list(S_star = S_star, intensity = intensity, S = S, T = T_target)
}

ng_estimate_phenotypic_covariance <- function(Y, shrinkage = c("none", "auto")) {
  shrinkage <- match.arg(shrinkage)
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  trait_names <- colnames(Y)
  if (shrinkage == "none") {
    P <- stats::cov(Y, use = "pairwise.complete.obs")
    P[!is.finite(P)] <- 0
    P <- (P + t(P)) / 2
    P <- ng_genetic_cov_project_psd(P)
    dimnames(P) <- list(trait_names, trait_names)
    attr(P, "shrinkage") <- "none"
    attr(P, "intensity") <- 0
    return(P)
  }
  # Ledoit-Wolf auto shrinkage toward a diagonal target.
  lw <- ng_genetic_cov_ledoit_wolf(Y)
  P <- lw$S_star
  P <- (P + t(P)) / 2
  P <- ng_genetic_cov_project_psd(P)
  dimnames(P) <- list(trait_names, trait_names)
  attr(P, "shrinkage") <- "ledoit_wolf"
  attr(P, "intensity") <- lw$intensity
  attr(P, "sample_covariance") <- lw$S
  P
}

# --- posterior G across traits ---------------------------------------------

# Posterior over the additive genetic covariance matrix G across traits.
#
# Two modes:
#
#   method = "beta_posterior" (fast diagnostic approximation).
#     For each trait, fit per-trait ridge once and draw S samples of beta_t
#     from the BCM closed-form posterior conditional on (sigma_e2_t, lambda_t).
#     For draw s, construct G_s = D R_beta^s D where D = sqrt(sigma_g2) is
#     fixed across draws (lambda-inversion is hyperparameter-conditional) and
#     R_beta^s = cor(beta_1^s, ..., beta_t^s) recomputed per draw. This is
#     This is not a joint multi-trait posterior: the diagonals are fixed and
#     only the plug-in marker-effect correlation carries uncertainty.
#
#   method = "parametric_bootstrap" (slower diagnostic approximation).
#     Resample Gaussian residuals around each per-trait ridge fit, regenerate
#     y_t^b = fitted_t + epsilon_t^b, and rerun ng_genetic_cov_two_stage_ridge
#     end-to-end (CV-tuning lambda included). It propagates uncertainty within
#     the same separate-trait heuristic; it does not create a joint multivariate
#     variance-component model.
#
# Returns a t x t x n_draws array with attribute `G_mean` (posterior mean
# matrix), `method`, `n_draws`. Off-diagonal credible intervals can be read
# off via apply(G_draws, c(1, 2), quantile, probs = c(.025, .975)).
ng_posterior_genetic_covariance <- function(geno,
                                            Y,
                                            n_draws = 100L,
                                            method = c("beta_posterior", "parametric_bootstrap"),
                                            ridge_lambda = NULL,
                                            kfold = 5L,
                                            seed = 1L,
                                            allow_heuristic = FALSE) {
  method <- match.arg(method)
  if (!isTRUE(allow_heuristic)) {
    ng_stop("ng_posterior_genetic_covariance currently uses univariate ridge heuristics, ",
            "not a joint multivariate posterior; set allow_heuristic = TRUE only for diagnostics")
  }
  warning("returning heuristic covariance draws from separate univariate ridge fits",
          call. = FALSE)
  geno <- ng_as_numeric_matrix(geno, "geno")
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) ng_stop("Y must have column names (one per trait)")
  if (nrow(Y) != nrow(geno)) ng_stop("Y and geno must have the same number of rows")
  ids <- rownames(geno)
  trait_names <- colnames(Y)
  n_traits <- ncol(Y)
  m <- ncol(geno)
  p <- colMeans(geno, na.rm = TRUE) / 2
  denom <- sum(2 * p * (1 - p), na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0) denom <- m

  G_draws <- array(0, dim = c(n_traits, n_traits, n_draws),
                   dimnames = list(trait_names, trait_names, NULL))

  if (identical(method, "beta_posterior")) {
    sigma_g2 <- numeric(n_traits); names(sigma_g2) <- trait_names
    beta_draw_list <- vector("list", n_traits)
    for (j in seq_len(n_traits)) {
      y_j <- as.numeric(Y[, j])
      names(y_j) <- ids
      # SEED (0.29.0). ng_fit_ridge_effects_posterior() uses ONE `seed` for two
      # different jobs: the k-fold CV partition that selects lambda, and the
      # posterior innovations. Those jobs need opposite treatments (see
      # ng_genetic_cov_two_stage_ridge() above), so they are separated here by
      # selecting lambda first and passing it in, which removes lambda selection
      # from the seed's remit entirely.
      #
      #   lambda -- chosen on the SHARED fold partition (the base `seed`), which
      #     is both the order-invariant choice and the one that preserves this
      #     function's documented invariant: because sigma_g2 = sigma_e2 * denom /
      #     lambda is hyperparameter-conditional and fixed across draws, the
      #     beta_posterior diagonals equal the point ng_estimate_genetic_covariance()
      #     diagonals EXACTLY. Letting the posterior pick lambda on a different
      #     partition from the point estimator broke that by 1.3 in testing.
      #   innovations -- DISTINCT per trait and keyed on the trait NAME. A shared
      #     stream would give every trait identical draws and manufacture
      #     cross-trait correlation in exactly the R_beta_b computed below.
      #     `seed + j` kept them distinct but keyed to the trait's COLUMN
      #     POSITION, so permuting Y moved every trait onto a different stream.
      lambda_j <- ridge_lambda
      if (is.null(lambda_j)) {
        lambda_j <- ng_fit_ridge_effects(
          geno = geno, y = y_j, ids = ids,
          lambda = NULL, kfold = kfold, seed = seed
        )$lambda
      }
      post_j <- ng_fit_ridge_effects_posterior(
        geno = geno, y = y_j, ids = ids,
        lambda = lambda_j, kfold = kfold,
        n_draws = n_draws, method = "closed_form",
        seed = ng_trait_rng_seed(seed, trait_names[[j]])
      )
      sigma_g2[[j]] <- post_j$fit$sigma_e2 * denom / max(post_j$fit$lambda, 1e-12)
      if (!is.finite(sigma_g2[[j]]) || sigma_g2[[j]] <= 0) sigma_g2[[j]] <- 1e-8
      beta_draw_list[[j]] <- post_j$beta_draws
    }
    D <- sqrt(sigma_g2)
    for (b in seq_len(n_draws)) {
      B_b <- do.call(cbind, lapply(beta_draw_list, function(M) M[, b]))
      if (m > 1L) {
        R_beta_b <- suppressWarnings(stats::cor(B_b))
        R_beta_b[!is.finite(R_beta_b)] <- 0
        diag(R_beta_b) <- 1
      } else {
        R_beta_b <- diag(n_traits)
      }
      G_b <- diag(D, n_traits) %*% R_beta_b %*% diag(D, n_traits)
      G_b <- (G_b + t(G_b)) / 2
      G_b <- ng_genetic_cov_project_psd(G_b)
      dimnames(G_b) <- list(trait_names, trait_names)
      G_draws[, , b] <- G_b
    }
  } else {
    # parametric_bootstrap: resample residuals around each per-trait ridge fit.
    fits <- vector("list", n_traits)
    fitted_list <- vector("list", n_traits)
    ok_list <- vector("list", n_traits)
    for (j in seq_len(n_traits)) {
      y_j <- as.numeric(Y[, j])
      names(y_j) <- ids
      # SEED (0.29.0): shared fold partition, as in ng_genetic_cov_two_stage_ridge()
      # above -- this is the lambda-CV seed, not an innovation stream.
      fit_j <- ng_fit_ridge_effects(
        geno = geno, y = y_j, ids = ids,
        lambda = ridge_lambda, kfold = kfold, seed = seed
      )
      fits[[j]] <- fit_j
      ok_j <- is.finite(y_j)
      ok_list[[j]] <- ok_j
      fitted_list[[j]] <- as.numeric(ng_predict_gebv(geno[ok_j, , drop = FALSE], fit_j))
    }
    # RESIDUAL DRAWS (0.29.0): one identity-derived stream per (trait, draw).
    #
    # The previous code did a single set.seed(seed) outside both loops and then
    # consumed one shared stream in COLUMN ORDER, so trait j's bootstrap
    # residuals were whatever the stream happened to hold after traits 1..j-1 had
    # taken theirs. Permuting the columns of Y therefore handed every trait a
    # different residual sample -- the same position dependence as the seeds
    # above, just expressed through stream position instead of a seed offset.
    #
    # Each (trait, draw) now gets its own stream keyed on the trait NAME. `salt`
    # multiplies the draw index by a large prime so that the (name hash + salt)
    # sums of two different (trait, draw) pairs cannot coincide except by an
    # astronomically improbable hash collision, which would otherwise give two
    # traits identical residuals in some draw and manufacture correlation.
    # ng_with_rng_seed() restores the caller's global RNG state afterwards.
    salt_stride <- 1000003
    for (b in seq_len(n_draws)) {
      Y_b <- Y
      for (j in seq_len(n_traits)) {
        eps_j <- ng_with_rng_seed(
          ng_trait_rng_seed(seed, trait_names[[j]], salt = b * salt_stride),
          stats::rnorm(sum(ok_list[[j]]), sd = sqrt(fits[[j]]$sigma_e2)))
        Y_b[ok_list[[j]], j] <- fitted_list[[j]] + eps_j
      }
      out_b <- ng_genetic_cov_two_stage_ridge(
        geno = geno, Y = Y_b,
        ridge_lambda = ridge_lambda, kfold = kfold, seed = seed + b
      )
      G_draws[, , b] <- out_b$G_hat
    }
  }

  G_mean <- apply(G_draws, c(1L, 2L), mean)
  dimnames(G_mean) <- list(trait_names, trait_names)
  attr(G_draws, "G_mean") <- G_mean
  attr(G_draws, "method") <- method
  attr(G_draws, "n_draws") <- n_draws
  attr(G_draws, "n_used") <- nrow(Y)
  G_draws
}
