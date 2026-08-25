ng_duplicate_recode_to_012 <- function(geno) {
  geno <- as.matrix(geno)
  storage.mode(geno) <- "double"
  vals <- unique(as.numeric(geno[is.finite(geno)]))
  # A -1 observation disambiguates centered {-1,0,1} coding from ordinary
  # binary 0/1 data. Inspect the full matrix so a late dosage-2 call is not
  # converted to the invalid value 3.
  if (length(vals) && any(vals == -1) && all(vals %in% c(-1, 0, 1))) geno <- geno + 1
  geno
}

ng_duplicate_filter_markers <- function(geno, maf_min = 0, max_missing_prop = 1) {
  geno <- as.matrix(geno)
  storage.mode(geno) <- "double"
  miss_prop <- colMeans(!is.finite(geno))
  p <- colMeans(geno, na.rm = TRUE) / 2
  p[is.nan(p)] <- NA_real_
  maf <- pmin(p, 1 - p)
  keep <- miss_prop <= max_missing_prop & is.finite(maf) & maf > 0 & maf >= maf_min
  list(
    geno = geno[, keep, drop = FALSE],
    marker_summary = data.frame(
      metric = c("markers_input", "markers_used", "markers_removed"),
      value = c(ncol(geno), sum(keep), sum(!keep)),
      stringsAsFactors = FALSE
    )
  )
}

ng_duplicate_ibs01 <- function(geno) {
  obs <- is.finite(geno)
  g0 <- geno == 0 & obs
  g1 <- geno == 1 & obs
  g2 <- geno == 2 & obs
  compared <- tcrossprod(obs)
  ibs2 <- tcrossprod(g0) + tcrossprod(g1) + tcrossprod(g2)
  ibs1 <- tcrossprod(g0, g1) + tcrossprod(g1, g0) + tcrossprod(g1, g2) + tcrossprod(g2, g1)
  sim <- (2 * ibs2 + ibs1) / (2 * pmax(compared, 1))
  sim[compared <= 0] <- NA_real_
  diag(sim) <- 1
  list(similarity = sim, compared_markers = compared)
}

ng_duplicate_pairs_from_similarity <- function(similarity,
                                               compared_markers,
                                               ids,
                                               threshold,
                                               min_compared_markers) {
  n <- nrow(similarity)
  if (n < 2L) {
    return(data.frame(
      parent1 = character(),
      parent2 = character(),
      similarity = numeric(),
      markers_compared = integer(),
      stringsAsFactors = FALSE
    ))
  }
  ut <- which(upper.tri(similarity), arr.ind = TRUE)
  sim <- as.numeric(similarity[ut])
  compared <- as.integer(compared_markers[ut])
  keep <- is.finite(sim) & sim >= threshold & compared >= min_compared_markers
  if (!any(keep)) {
    return(data.frame(
      parent1 = character(),
      parent2 = character(),
      similarity = numeric(),
      markers_compared = integer(),
      stringsAsFactors = FALSE
    ))
  }
  out <- data.frame(
    parent1 = ids[ut[keep, 1L]],
    parent2 = ids[ut[keep, 2L]],
    similarity = round(sim[keep], 6),
    markers_compared = compared[keep],
    stringsAsFactors = FALSE
  )
  out[order(-out$similarity, -out$markers_compared, out$parent1, out$parent2), , drop = FALSE]
}

ng_duplicate_clusters <- function(pairs) {
  if (!nrow(pairs)) return(list())
  ids <- unique(c(as.character(pairs$parent1), as.character(pairs$parent2)))
  parent <- stats::setNames(ids, ids)
  find <- function(x) {
    while (!identical(parent[[x]], x)) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  union <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (!identical(ra, rb)) parent[[rb]] <<- ra
  }
  for (i in seq_len(nrow(pairs))) union(as.character(pairs$parent1[[i]]), as.character(pairs$parent2[[i]]))
  root <- vapply(ids, find, character(1L))
  clusters <- split(ids, root)
  clusters <- lapply(clusters, sort)
  clusters[order(lengths(clusters), decreasing = TRUE)]
}

ng_duplicate_nearest_neighbor <- function(similarity, compared_markers, ids) {
  if (length(ids) < 2L) {
    return(data.frame(parent = ids, nearest_parent = NA_character_, similarity = NA_real_,
                      markers_compared = NA_integer_, stringsAsFactors = FALSE))
  }
  sim <- similarity
  diag(sim) <- -Inf
  idx <- max.col(sim, ties.method = "first")
  data.frame(
    parent = ids,
    nearest_parent = ids[idx],
    similarity = round(sim[cbind(seq_along(ids), idx)], 6),
    markers_compared = as.integer(compared_markers[cbind(seq_along(ids), idx)]),
    stringsAsFactors = FALSE
  )
}

ng_putative_duplicate_summary <- function(pairs, clusters, n_parents, n_markers_used) {
  dup_ids <- if (nrow(pairs)) unique(c(as.character(pairs$parent1), as.character(pairs$parent2))) else character()
  cluster_sizes <- if (length(clusters)) lengths(clusters) else integer()
  data.frame(
    metric = c(
      "parents_checked",
      "markers_used",
      "duplicate_pairs",
      "putative_duplicate_parents",
      "duplicate_clusters",
      "largest_cluster_size",
      "duplicate_excess_if_keep_one_per_cluster"
    ),
    value = c(
      n_parents,
      n_markers_used,
      nrow(pairs),
      length(dup_ids),
      length(clusters),
      if (length(cluster_sizes)) max(cluster_sizes) else 0L,
      if (length(cluster_sizes)) sum(pmax(cluster_sizes - 1L, 0L)) else 0L
    ),
    stringsAsFactors = FALSE
  )
}

ng_detect_putative_duplicates <- function(geno,
                                          ids = rownames(geno),
                                          maf_min = 0,
                                          max_missing_prop = 1,
                                          duplicate_threshold = 0.995,
                                          min_compared_markers = 100L,
                                          return_similarity = TRUE) {
  geno <- as.matrix(geno)
  if (is.null(ids)) ids <- rownames(geno)
  if (is.null(ids)) ids <- sprintf("id_%05d", seq_len(nrow(geno)))
  ids <- trimws(as.character(ids))
  if (length(ids) != nrow(geno)) ng_stop("ids length must equal nrow(geno)")
  if (any(!nzchar(ids) | is.na(ids))) ng_stop("ids must not contain missing or blank values")
  if (anyDuplicated(ids)) ng_stop("ids must be unique before putative duplicate detection")
  duplicate_threshold <- suppressWarnings(as.numeric(duplicate_threshold[[1L]]))
  if (!is.finite(duplicate_threshold) || duplicate_threshold <= 0 || duplicate_threshold > 1) {
    ng_stop("duplicate_threshold must be in (0, 1]")
  }
  maf_min <- suppressWarnings(as.numeric(maf_min[[1L]]))
  max_missing_prop <- suppressWarnings(as.numeric(max_missing_prop[[1L]]))
  min_compared_markers <- suppressWarnings(as.integer(min_compared_markers[[1L]]))
  if (!is.finite(maf_min) || maf_min < 0 || maf_min > 0.5) ng_stop("maf_min must be in [0, 0.5]")
  if (!is.finite(max_missing_prop) || max_missing_prop < 0 || max_missing_prop > 1) {
    ng_stop("max_missing_prop must be in [0, 1]")
  }
  if (!is.finite(min_compared_markers) || min_compared_markers < 1L) min_compared_markers <- 1L

  geno <- ng_duplicate_recode_to_012(geno)
  filtered <- ng_duplicate_filter_markers(geno, maf_min = maf_min, max_missing_prop = max_missing_prop)
  if (ncol(filtered$geno) < 1L) ng_stop("No polymorphic markers remain after duplicate-detection filtering")
  rownames(filtered$geno) <- ids
  ibs <- ng_duplicate_ibs01(filtered$geno)
  rownames(ibs$similarity) <- colnames(ibs$similarity) <- ids
  rownames(ibs$compared_markers) <- colnames(ibs$compared_markers) <- ids
  pairs <- ng_duplicate_pairs_from_similarity(
    similarity = ibs$similarity,
    compared_markers = ibs$compared_markers,
    ids = ids,
    threshold = duplicate_threshold,
    min_compared_markers = min_compared_markers
  )
  clusters <- ng_duplicate_clusters(pairs)
  out <- list(
    pairs = pairs,
    clusters = clusters,
    summary = ng_putative_duplicate_summary(pairs, clusters, nrow(geno), ncol(filtered$geno)),
    nearest_neighbor = ng_duplicate_nearest_neighbor(ibs$similarity, ibs$compared_markers, ids),
    marker_summary = filtered$marker_summary,
    params = list(
      method = "IBS01",
      duplicate_threshold = duplicate_threshold,
      maf_min = maf_min,
      max_missing_prop = max_missing_prop,
      min_compared_markers = min_compared_markers
    )
  )
  if (isTRUE(return_similarity)) {
    out$similarity <- ibs$similarity
    out$compared_markers <- ibs$compared_markers
  }
  class(out) <- c("ng_putative_duplicates", "list")
  out
}

ng_plot_open_device <- function(output_path, width, height, res = 150) {
  ext <- tolower(tools::file_ext(output_path))
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (identical(ext, "pdf")) {
    grDevices::pdf(output_path, width = width, height = height)
  } else {
    grDevices::png(output_path, width = width * res, height = height * res, res = res)
  }
  invisible(grDevices::dev.cur())
}

# Close a device opened by ng_plot_open_device, flushing the file to disk. Safe to
# call more than once (no-op if already closed) so it can be used BOTH explicitly
# -- before verifying the output path exists -- and from on.exit() as an error
# safety net. This matters because png()/pdf() may not materialize the output file
# until the device is closed (observed on macOS); normalizePath(mustWork = TRUE)
# run while the device is still open would then spuriously fail.
ng_plot_close_device <- function(dev_no) {
  if (!is.null(dev_no) && dev_no %in% grDevices::dev.list()) grDevices::dev.off(dev_no)
  invisible(NULL)
}

ng_plot_putative_duplicates <- function(x,
                                        output_path = NULL,
                                        top_k_groups = 10L,
                                        max_ids = 250L,
                                        title = NULL) {
  if (!inherits(x, "ng_putative_duplicates")) ng_stop("x must be returned by ng_detect_putative_duplicates()")
  if (is.null(output_path)) {
    output_path <- tempfile("ng_putative_duplicates_", fileext = ".png")
  }
  dev_no <- ng_plot_open_device(output_path, width = 9, height = 7, res = 160)
  on.exit(ng_plot_close_device(dev_no), add = TRUE)

  title <- if (is.null(title)) "Putative duplicate genotype similarity" else as.character(title[[1L]])
  old_par <- graphics::par(no.readonly = TRUE)
  # Restore par only if a device is still open; once we close our owned device the
  # restore would otherwise spawn a stray default device (Rplots.pdf).
  on.exit(if (length(grDevices::dev.list())) graphics::par(old_par), add = TRUE)
  if (!nrow(x$pairs) || is.null(x$similarity)) {
    graphics::plot.new()
    graphics::title(main = title)
    graphics::text(0.5, 0.58, "No putative duplicate pairs exceeded the threshold", cex = 1.05)
    graphics::text(
      0.5, 0.48,
      sprintf("threshold = %.4f; markers used = %s",
              x$params$duplicate_threshold,
              x$summary$value[x$summary$metric == "markers_used"]),
      cex = 0.9
    )
    ng_plot_close_device(dev_no)
    return(invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE)))
  }

  groups <- x$clusters
  groups <- groups[order(lengths(groups), decreasing = TRUE)]
  groups <- utils::head(groups, top_k_groups)
  ids <- unique(unlist(groups, use.names = FALSE))
  if (length(ids) > max_ids) {
    degree <- table(c(as.character(x$pairs$parent1), as.character(x$pairs$parent2)))
    ids <- names(sort(degree[ids], decreasing = TRUE))[seq_len(max_ids)]
  }
  cluster_id <- rep(NA_integer_, length(ids))
  names(cluster_id) <- ids
  for (i in seq_along(groups)) cluster_id[intersect(groups[[i]], ids)] <- i
  ids <- ids[order(cluster_id[ids], ids)]
  sim <- x$similarity[ids, ids, drop = FALSE]

  graphics::par(mar = c(7, 7, 4, 6))
  cols <- grDevices::colorRampPalette(c("#f7fbff", "#fdd66f", "#c53b2c"))(100)
  graphics::image(
    x = seq_len(ncol(sim)),
    y = seq_len(nrow(sim)),
    z = t(sim[nrow(sim):1L, , drop = FALSE]),
    col = cols,
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = title
  )
  graphics::axis(3, at = seq_len(ncol(sim)), labels = colnames(sim), las = 2, cex.axis = 0.55)
  graphics::axis(2, at = seq_len(nrow(sim)), labels = rev(rownames(sim)), las = 2, cex.axis = 0.55)
  graphics::box()
  graphics::mtext(sprintf("IBS similarity, threshold >= %.4f", x$params$duplicate_threshold),
                  side = 1, line = 5, cex = 0.85)
  ng_plot_close_device(dev_no)
  invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
}
