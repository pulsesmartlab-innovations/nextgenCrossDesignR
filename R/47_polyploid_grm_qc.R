# Correct polyploid genomic relationship (GRM) and quality control for allele-dosage data.
#
# The earlier ng_poly4x_parent_relationship centers dosage on the ploidy MIDPOINT (ploidy/2),
# i.e. it implicitly assumes every allele frequency is 0.5 -- a diploid-oriented shortcut. A
# correct additive GRM centers on the ACTUAL allele frequency and scales by the genotypic
# variance under random mating, generalizing VanRaden (2008) to arbitrary ploidy (Endelman 2018;
# AGHmatrix / Kerr 2012 "VanRaden"): for dosage M in 0..ploidy with allele frequency p,
#   W = M - ploidy * p ,   G = W W' / ( ploidy * sum_k p_k (1 - p_k) ).
# This reduces exactly to VanRaden's diploid G when ploidy = 2, and gives diagonal ~ 1 + F.

# Per-marker allele frequency from dosage (mean dosage / ploidy), ignoring missing values.
ng_polyploid_allele_freq <- function(dosage, ploidy) {
  colMeans(dosage, na.rm = TRUE) / ploidy
}

# Shared GRM preamble: validate ploidy + dosage range, impute missing to expected dosage, drop
# monomorphic / below-MAF markers. Returns the cleaned dosage (M), kept allele frequencies (p),
# sample ids, and integer ploidy. Used by both the additive and dominance GRMs.
ng_polyploid_prep_dosage <- function(dosage, ploidy, min_maf = 0, impute_missing = TRUE) {
  ploidy <- suppressWarnings(as.numeric(ploidy)[[1]])
  if (!is.finite(ploidy) || ploidy < 2 || abs(ploidy - round(ploidy)) > 1e-8) {
    ng_stop("ploidy must be an integer >= 2")
  }
  ploidy <- as.integer(ploidy)
  M <- ng_as_numeric_matrix(dosage, "dosage")
  ids <- rownames(M)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) ng_stop("dosage must have sample (row) names")
  if (anyDuplicated(ids)) ng_stop("dosage sample (row) names must be unique")
  obs <- is.finite(M)
  bad <- obs & (M < 0 | M > ploidy | abs(M - round(M)) > 1e-8)
  if (any(bad)) {
    idx <- which(bad, arr.ind = TRUE)[1, ]
    ng_stop("dosage outside 0..", ploidy, " at row ", ids[[idx[[1]]]], ", column ", colnames(M)[[idx[[2]]]])
  }
  p <- ng_polyploid_allele_freq(M, ploidy)
  if (impute_missing && any(!obs)) {
    na_idx <- which(!obs, arr.ind = TRUE)
    M[na_idx] <- (ploidy * p)[na_idx[, 2]]
  } else if (any(!obs)) {
    ng_stop("dosage has missing values; set impute_missing = TRUE or clean first (ng_polyploid_qc)")
  }
  maf <- pmin(p, 1 - p)
  keep <- is.finite(maf) & maf >= min_maf & p > 0 & p < 1
  if (!any(keep)) ng_stop("no polymorphic markers pass min_maf = ", min_maf)
  list(M = M[, keep, drop = FALSE], p = p[keep], ids = ids, ploidy = ploidy)
}

# Correct polyploid additive GRM, generalizing two standard estimators to arbitrary ploidy.
# Both center on the ACTUAL allele frequency (W = M - ploidy * p):
#   method = "vanraden" (default): a single overall scaling, G = W W' / (ploidy * sum_k p_k(1-p_k))
#     -- VanRaden (2008), weights rare-allele markers more; reduces exactly to the diploid G.
#   method = "yang": per-marker standardization to unit variance (GCTA / Yang et al. 2010),
#     Z_ik = W_ik / sqrt(ploidy * p_k(1-p_k)), G = Z Z' / m -- weights every marker equally.
# Missing dosages are imputed to the marker's expected dosage (ploidy * p); monomorphic and
# (optionally) low-MAF markers are dropped (no relationship information; avoids divide-by-zero).
ng_polyploid_grm <- function(dosage,
                             ploidy = 2L,
                             method = c("vanraden", "yang"),
                             min_maf = 0,
                             impute_missing = TRUE,
                             return_freq = FALSE) {
  method <- match.arg(method)
  prep <- ng_polyploid_prep_dosage(dosage, ploidy, min_maf = min_maf, impute_missing = impute_missing)
  Mk <- prep$M; ids <- prep$ids; pk <- prep$p; ploidy <- prep$ploidy
  W <- sweep(Mk, 2L, ploidy * pk, "-")
  if (identical(method, "yang")) {
    Z <- sweep(W, 2L, sqrt(ploidy * pk * (1 - pk)), "/")
    G <- tcrossprod(Z) / ncol(Z)
  } else {
    denom <- ploidy * sum(pk * (1 - pk))
    if (!is.finite(denom) || denom <= 0) ng_stop("GRM scaling denominator is not positive")
    G <- tcrossprod(W) / denom
  }
  dimnames(G) <- list(ids, ids)
  attr(G, "ploidy") <- ploidy
  attr(G, "n_markers") <- ncol(Mk)
  attr(G, "method") <- paste0(method, "_polyploid")
  if (isTRUE(return_freq)) attr(G, "allele_freq") <- pk
  G
}

# Polyploid DOMINANCE (digenic) relationship matrix -- needed when both additive and dominance
# effects drive parent selection and crossing (e.g. cassava and other outbred clonal crops). The
# digenic dominance covariate at a marker is the number of heterozygous allele pairs, H = d*(ploidy
# - d) (for diploid this is the classic heterozygote indicator scaled), centered on its mean. As
# with the additive GRM: "vanraden" uses one overall scaling; "yang" standardizes each marker.
# The crossproduct is BLAS-backed; for very large marker sets the fused C++ kernel avoids
# materializing the design matrix (used when compiled and use_cpp = TRUE).
ng_polyploid_dominance_grm <- function(dosage,
                                       ploidy = 2L,
                                       method = c("vanraden", "yang"),
                                       min_maf = 0,
                                       impute_missing = TRUE,
                                       return_freq = FALSE,
                                       use_cpp = TRUE) {
  method <- match.arg(method)
  prep <- ng_polyploid_prep_dosage(dosage, ploidy, min_maf = min_maf, impute_missing = impute_missing)
  M <- prep$M; ids <- prep$ids; pk <- prep$p; ploidy <- prep$ploidy
  H <- M * (ploidy - M)                              # digenic heterozygosity covariate
  hbar <- colMeans(H)
  D <- sweep(H, 2L, hbar, "-")                       # centered dominance design
  if (identical(method, "yang")) {
    sdv <- sqrt(apply(H, 2L, stats::var)); sdv[!is.finite(sdv) | sdv <= 0] <- 1
    Z <- sweep(D, 2L, sdv, "/")
    G <- tcrossprod(Z) / ncol(Z)
  } else {
    denom <- sum(apply(H, 2L, stats::var))
    if (!is.finite(denom) || denom <= 0) ng_stop("dominance GRM scaling denominator is not positive")
    G <- tcrossprod(D) / denom
  }
  dimnames(G) <- list(ids, ids)
  attr(G, "ploidy") <- ploidy
  attr(G, "n_markers") <- ncol(M)
  attr(G, "method") <- paste0(method, "_dominance_polyploid")
  attr(G, "component") <- "dominance"
  if (isTRUE(return_freq)) attr(G, "allele_freq") <- pk
  G
}

# Ploidy-aware QC for a dosage matrix (samples x markers, dosage 0..ploidy). Flags out-of-range /
# non-integer dosages, per-marker and per-sample missingness, monomorphic and low-MAF markers, and
# duplicate sample/marker IDs, and returns a cleaned matrix ready for GRM / cross scoring. Unlike
# the diploid preflight this validates against 0..ploidy and computes MAF as mean-dosage/ploidy.
ng_polyploid_qc <- function(dosage,
                            ploidy = 2L,
                            max_missing_marker = 0.20,
                            max_missing_sample = 0.20,
                            min_maf = 0.0,
                            drop_monomorphic = TRUE) {
  ploidy <- suppressWarnings(as.numeric(ploidy)[[1]])
  if (!is.finite(ploidy) || ploidy < 2 || abs(ploidy - round(ploidy)) > 1e-8) {
    ng_stop("ploidy must be an integer >= 2")
  }
  ploidy <- as.integer(ploidy)
  M <- ng_as_numeric_matrix(dosage, "dosage")
  ids <- rownames(M); mk <- colnames(M)
  if (is.null(ids)) ids <- sprintf("sample%03d", seq_len(nrow(M)))
  if (is.null(mk)) mk <- sprintf("marker%03d", seq_len(ncol(M)))

  # out-of-range / non-integer (observed) dosages -> treated as missing after being flagged
  obs <- is.finite(M)
  out_of_range <- obs & (M < 0 | M > ploidy | abs(M - round(M)) > 1e-8)
  n_out_of_range <- sum(out_of_range)
  M[out_of_range] <- NA_real_
  miss <- !is.finite(M)

  marker_missing <- colMeans(miss)
  sample_missing <- rowMeans(miss)
  p <- ng_polyploid_allele_freq(M, ploidy)
  maf <- pmin(p, 1 - p)
  monomorphic <- !is.finite(maf) | p <= 0 | p >= 1

  drop_marker <- marker_missing > max_missing_marker |
    (is.finite(maf) & maf < min_maf) |
    (isTRUE(drop_monomorphic) & monomorphic)
  drop_sample <- sample_missing > max_missing_sample
  dup_samples <- unique(ids[duplicated(ids)])
  dup_markers <- unique(mk[duplicated(mk)])

  keep_m <- !drop_marker; keep_s <- !drop_sample
  clean <- M[keep_s, keep_m, drop = FALSE]

  marker_report <- data.frame(marker = mk, missing = marker_missing, allele_freq = p,
                              maf = maf, monomorphic = monomorphic, dropped = drop_marker,
                              stringsAsFactors = FALSE)
  sample_report <- data.frame(sample = ids, missing = sample_missing, dropped = drop_sample,
                              stringsAsFactors = FALSE)
  summary <- list(
    ploidy = ploidy, n_samples = nrow(M), n_markers = ncol(M),
    n_out_of_range = n_out_of_range, overall_missing = mean(miss),
    n_markers_dropped = sum(drop_marker), n_samples_dropped = sum(drop_sample),
    n_monomorphic = sum(monomorphic), n_duplicate_samples = length(dup_samples),
    n_duplicate_markers = length(dup_markers),
    n_samples_clean = nrow(clean), n_markers_clean = ncol(clean)
  )
  pass <- n_out_of_range == 0L && !length(dup_samples) && ncol(clean) > 0L && nrow(clean) > 1L
  list(pass = pass, summary = summary, marker_report = marker_report,
       sample_report = sample_report, duplicate_samples = dup_samples,
       duplicate_markers = dup_markers, clean = clean)
}
