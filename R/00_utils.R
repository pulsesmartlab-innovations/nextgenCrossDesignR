ng_stop <- function(..., call. = FALSE) stop(paste0(...), call. = call.)

ng_as_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(rownames(x))) ng_stop(name, " must have row names")
  if (is.null(colnames(x))) colnames(x) <- paste0("m", seq_len(ncol(x)))
  x
}

ng_check_same_ids <- function(x, ids, object_name = "matrix") {
  if (is.null(rownames(x))) ng_stop(object_name, " must have row names")
  miss <- setdiff(ids, rownames(x))
  if (length(miss)) ng_stop(object_name, " is missing ", length(miss), " requested ids")
  x[ids, , drop = FALSE]
}

ng_make_pairs <- function(ids, include_self = FALSE) {
  ids <- as.character(ids)
  out_n <- if (include_self) length(ids) * (length(ids) + 1) / 2 else length(ids) * (length(ids) - 1) / 2
  p1 <- character(out_n)
  p2 <- character(out_n)
  k <- 1L
  for (i in seq_along(ids)) {
    j0 <- if (include_self) i else i + 1L
    if (j0 <= length(ids)) {
      for (j in j0:length(ids)) {
        p1[k] <- ids[i]
        p2[k] <- ids[j]
        k <- k + 1L
      }
    }
  }
  data.frame(parent1 = p1, parent2 = p2, stringsAsFactors = FALSE)
}

ng_selection_intensity <- function(prop_selected, n_progeny = NULL) {
  p <- as.numeric(prop_selected)
  if (!is.finite(p) || p <= 0 || p >= 1) ng_stop("prop_selected must be in (0, 1)")
  if (is.null(n_progeny) || !is.finite(n_progeny[[1L]]) || n_progeny[[1L]] <= 0) {
    z <- stats::qnorm(1 - p)
    return(stats::dnorm(z) / p)
  }
  # Finite-population correction. Use the expected mean of the top k order
  # statistics of an N(0,1) sample of size n, via Blom (1958) plotting
  # positions: E[Z_{(j:n)}] ~= Phi^{-1}((j - 3/8) / (n + 1/4)). This is
  # deterministic, fast, and converges to phi(z)/p as n -> Inf.
  n <- max(1L, as.integer(round(n_progeny[[1L]])))
  k <- max(1L, as.integer(round(n * p)))
  k <- min(k, n)
  j <- seq.int(n - k + 1L, n)
  mean(stats::qnorm((j - 0.375) / (n + 0.25)))
}

ng_scale01 <- function(x, bigger_is_better = TRUE) {
  x <- as.numeric(x)
  r <- range(x[is.finite(x)], na.rm = TRUE)
  if (!all(is.finite(r)) || abs(diff(r)) < .Machine$double.eps) return(rep(0.5, length(x)))
  y <- (x - r[1]) / diff(r)
  if (!bigger_is_better) y <- 1 - y
  y
}

ng_standardize <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

ng_match_vector <- function(x, ids, name = deparse(substitute(x))) {
  if (is.null(names(x))) {
    if (length(x) != length(ids)) ng_stop(name, " must be named or have length ids")
    names(x) <- ids
  }
  miss <- setdiff(ids, names(x))
  if (length(miss)) ng_stop(name, " is missing ", length(miss), " ids")
  as.numeric(x[ids])
}

ng_parent_kinship <- function(geno, method = c("vanraden", "yang")) {
  method <- match.arg(method)
  geno <- ng_as_numeric_matrix(geno, "geno")
  p <- colMeans(geno, na.rm = TRUE) / 2
  X <- sweep(geno, 2, 2 * p, "-")
  X[!is.finite(X)] <- 0
  if (identical(method, "yang")) {
    # Yang / GCTA: standardize each marker to unit variance (equal weight), then average.
    s <- sqrt(2 * p * (1 - p)); s[!is.finite(s) | s <= 0] <- Inf   # monomorphic -> X/Inf = 0
    Z <- sweep(X, 2, s, "/"); Z[!is.finite(Z)] <- 0
    denom <- sum(is.finite(s)); if (denom <= 0) denom <- ncol(geno)
    K <- tcrossprod(Z) / denom
  } else {
    # VanRaden (default): single overall allele-frequency scaling.
    denom <- sum(2 * p * (1 - p), na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- ncol(geno)
    K <- tcrossprod(X) / denom
  }
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  attr(K, "method") <- method
  attr(K, "relationship_scale") <- "additive_relationship"
  attr(K, "coancestry_divisor") <- 2
  K
}

ng_pair_key <- function(parent1, parent2) paste(parent1, parent2, sep = "||")

# Evaluate `expr` under a fixed RNG seed, restoring the caller's global RNG state
# afterwards so seeded internal randomness (optimizers, benchmark fixtures) is
# reproducible without perturbing the surrounding stream. `expr` is a promise and is
# only forced AFTER the seed is set. A non-finite/NULL seed evaluates `expr` as-is.
ng_with_rng_seed <- function(seed, expr) {
  s <- suppressWarnings(as.numeric(seed[[1L]]))
  if (!length(s) || !is.finite(s)) return(force(expr))
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(as.integer(s))
  force(expr)
}

# Per-parent diagnostics for the inbred-parent assumption used by the
# DH/RIL recombination-variance kernel. Dosage is considered inbred at a
# marker when min(d, ploidy - d) <= tolerance (e.g. 0 or 2 +/- noise for
# diploids). Returns the row indices and IDs of parents whose fraction of
# heterozygous-looking markers exceeds `fraction_tolerance`, plus the worst
# observed fraction.
ng_audit_inbred_dosage <- function(geno, ploidy = 2,
                                   tolerance = 0.05,
                                   fraction_tolerance = 0.02) {
  if (!is.matrix(geno)) geno <- as.matrix(geno)
  d <- as.numeric(geno)
  dim(d) <- dim(geno)
  het_dist <- pmin(abs(d), abs(ploidy - d))
  het_dist[!is.finite(het_dist)] <- 0
  het_marker <- het_dist > tolerance
  frac <- rowMeans(het_marker, na.rm = TRUE)
  frac[!is.finite(frac)] <- 0
  bad <- which(frac > fraction_tolerance)
  ids <- if (!is.null(rownames(geno))) rownames(geno) else as.character(seq_len(nrow(geno)))
  list(
    violators = ids[bad],
    violators_index = bad,
    fraction = frac,
    max_fraction = if (length(frac)) max(frac, na.rm = TRUE) else 0,
    tolerance = tolerance,
    fraction_tolerance = fraction_tolerance
  )
}

# Optional external packages (SimpleMating, genomicMateSelectR) are used only for
# opt-in benchmark comparisons and are intentionally NOT declared dependencies:
# they are not on CRAN, and the package is fully functional without them (callers
# fall back to native methods). These helpers reach them through variable
# arguments so R CMD check does not treat them as undeclared imports, and so pak
# never tries to resolve them during dependency setup.
ng_has_optional_pkg <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

ng_optional_pkg_fun <- function(pkg, fun) {
  get(fun, envir = asNamespace(pkg), mode = "function")
}
