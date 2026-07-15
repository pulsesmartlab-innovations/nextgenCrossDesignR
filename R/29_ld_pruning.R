ng_ld_backend_available <- function() {
  exists("ng_ld_prune_graph_cpp", mode = "function", inherits = TRUE)
}

ng_ld_validate_inputs <- function(geno,
                                  window,
                                  r2_threshold,
                                  maf_threshold,
                                  ploidy,
                                  object_name = "geno") {
  geno <- ng_as_numeric_matrix(geno, object_name)
  if (any(duplicated(colnames(geno)))) {
    ng_stop(object_name, " marker names must be unique before LD pruning")
  }
  if (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 1) {
    ng_stop("window must be a single positive integer")
  }
  if (!is.numeric(r2_threshold) || length(r2_threshold) != 1L || !is.finite(r2_threshold) ||
      r2_threshold < 0 || r2_threshold > 1) {
    ng_stop("r2_threshold must be in [0, 1]")
  }
  if (!is.numeric(maf_threshold) || length(maf_threshold) != 1L || !is.finite(maf_threshold) ||
      maf_threshold < 0 || maf_threshold > 0.5) {
    ng_stop("maf_threshold must be in [0, 0.5]")
  }
  if (!is.numeric(ploidy) || length(ploidy) != 1L || !is.finite(ploidy) || ploidy <= 0) {
    ng_stop("ploidy must be a single positive number")
  }
  geno
}

ng_ld_marker_maf <- function(geno, ploidy = 2) {
  allele_frequency <- colMeans(geno, na.rm = TRUE) / ploidy
  maf <- pmin(allele_frequency, 1 - allele_frequency)
  maf[!is.finite(maf)] <- NA_real_
  maf
}

ng_ld_pair_r2 <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  n <- sum(keep)
  if (n < 2L) return(NA_real_)
  x <- x[keep]
  y <- y[keep]
  sx <- sum(x)
  sy <- sum(y)
  sxx <- sum(x * x)
  syy <- sum(y * y)
  sxy <- sum(x * y)
  denom <- (n * sxx - sx * sx) * (n * syy - sy * sy)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  r2 <- (n * sxy - sx * sy)^2 / denom
  if (is.finite(r2)) min(1, max(0, r2)) else NA_real_
}

ng_ld_prune_graph_r <- function(geno,
                                window = 100L,
                                r2_threshold = 0.9,
                                maf_threshold = 0.01,
                                ploidy = 2) {
  m <- ncol(geno)
  if (m < 1L) return(integer())
  maf <- ng_ld_marker_maf(geno, ploidy = ploidy)
  pass_maf <- is.finite(maf) & maf >= maf_threshold
  if (!any(pass_maf)) return(integer())
  if (m == 1L) return(if (pass_maf[[1L]]) 1L else integer())

  parent <- seq_len(m)
  rank <- integer(m)
  find_root <- function(x) {
    while (!identical(parent[[x]], x)) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  unite <- function(a, b) {
    ra <- find_root(a)
    rb <- find_root(b)
    if (identical(ra, rb)) return(invisible(FALSE))
    if (rank[[ra]] < rank[[rb]]) {
      tmp <- ra
      ra <- rb
      rb <- tmp
    }
    parent[[rb]] <<- ra
    if (identical(rank[[ra]], rank[[rb]])) rank[[ra]] <<- rank[[ra]] + 1L
    invisible(TRUE)
  }

  window <- as.integer(window)
  for (i in seq_len(m - 1L)) {
    if (!pass_maf[[i]]) next
    j_end <- min(m, i + window)
    if (j_end <= i) next
    for (j in seq.int(i + 1L, j_end)) {
      if (!pass_maf[[j]]) next
      r2 <- ng_ld_pair_r2(geno[, i], geno[, j])
      if (is.finite(r2) && r2 > r2_threshold) unite(i, j)
    }
  }

  roots <- vapply(seq_len(m), function(i) if (pass_maf[[i]]) find_root(i) else NA_integer_, integer(1L))
  best <- rep(NA_integer_, m)
  for (i in seq_len(m)) {
    if (!pass_maf[[i]]) next
    root <- roots[[i]]
    current <- best[[root]]
    if (is.na(current) || maf[[i]] > maf[[current]] ||
        (identical(maf[[i]], maf[[current]]) && i < current)) {
      best[[root]] <- i
    }
  }
  keep <- best[is.finite(best)]
  keep <- sort(as.integer(keep))
  keep
}

ng_ld_prune_indices <- function(geno,
                                window = 100L,
                                r2_threshold = 0.9,
                                maf_threshold = 0.01,
                                ploidy = 2,
                                backend = c("auto", "cpp", "r")) {
  backend <- match.arg(backend)
  use_cpp <- backend %in% c("auto", "cpp") && ng_ld_backend_available()
  force_cpp <- identical(backend, "cpp")
  if (isTRUE(use_cpp)) {
    out <- tryCatch(
      ng_ld_prune_graph_cpp(
        geno = geno,
        window = as.integer(window),
        r2_threshold = as.numeric(r2_threshold),
        maf_threshold = as.numeric(maf_threshold),
        ploidy = as.numeric(ploidy)
      ),
      error = function(e) {
        if (isTRUE(force_cpp)) stop(e)
        NULL
      }
    )
    if (!is.null(out)) {
      out <- as.integer(out)
      attr(out, "backend") <- "cpp"
      return(out)
    }
  } else if (isTRUE(force_cpp)) {
    ng_stop("C++ LD pruning backend is not available. Load the package with use_cpp = TRUE or install from source.")
  }
  out <- ng_ld_prune_graph_r(
    geno = geno,
    window = window,
    r2_threshold = r2_threshold,
    maf_threshold = maf_threshold,
    ploidy = ploidy
  )
  attr(out, "backend") <- "r"
  out
}

ng_ld_prune_markers <- function(geno,
                                window = 100L,
                                r2_threshold = 0.9,
                                maf_threshold = 0.01,
                                ploidy = 2,
                                backend = c("auto", "cpp", "r")) {
  backend <- match.arg(backend)
  geno <- ng_ld_validate_inputs(
    geno = geno,
    window = window,
    r2_threshold = r2_threshold,
    maf_threshold = maf_threshold,
    ploidy = ploidy
  )
  maf <- ng_ld_marker_maf(geno, ploidy = ploidy)
  pass_maf <- is.finite(maf) & maf >= maf_threshold
  keep_idx <- ng_ld_prune_indices(
    geno = geno,
    window = window,
    r2_threshold = r2_threshold,
    maf_threshold = maf_threshold,
    ploidy = ploidy,
    backend = backend
  )
  backend_used <- attr(keep_idx, "backend", exact = TRUE)
  if (is.null(backend_used)) backend_used <- "r"
  keep_markers <- colnames(geno)[keep_idx]
  report <- data.frame(
    backend = backend_used,
    markers_before = as.integer(ncol(geno)),
    markers_after = as.integer(length(keep_markers)),
    low_maf_removed = as.integer(sum(!pass_maf)),
    ld_redundant_removed = as.integer(sum(pass_maf) - length(keep_markers)),
    window = as.integer(window),
    r2_threshold = as.numeric(r2_threshold),
    maf_threshold = as.numeric(maf_threshold),
    ploidy = as.numeric(ploidy),
    stringsAsFactors = FALSE
  )
  list(
    keep_markers = keep_markers,
    drop_markers = setdiff(colnames(geno), keep_markers),
    marker_maf = stats::setNames(maf, colnames(geno)),
    report = report
  )
}

ng_ld_prune_geno <- function(geno,
                             window = 100L,
                             r2_threshold = 0.9,
                             maf_threshold = 0.01,
                             ploidy = 2,
                             backend = c("auto", "cpp", "r")) {
  backend <- match.arg(backend)
  geno <- ng_ld_validate_inputs(
    geno = geno,
    window = window,
    r2_threshold = r2_threshold,
    maf_threshold = maf_threshold,
    ploidy = ploidy
  )
  pruned <- ng_ld_prune_markers(
    geno = geno,
    window = window,
    r2_threshold = r2_threshold,
    maf_threshold = maf_threshold,
    ploidy = ploidy,
    backend = backend
  )
  out <- geno[, pruned$keep_markers, drop = FALSE]
  attr(out, "ld_pruning_report") <- pruned$report
  attr(out, "ld_pruning_drop_markers") <- pruned$drop_markers
  out
}
