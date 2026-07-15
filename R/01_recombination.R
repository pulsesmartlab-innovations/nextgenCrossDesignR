ng_haldane_decay <- function(dist_cm) {
  exp(-2 * abs(dist_cm) / 100)
}

ng_kosambi_decay <- function(dist_cm) {
  1 - tanh(2 * abs(dist_cm) / 100)
}

# Single-meiosis recombination fraction r from the chosen mapping function.
# Preserves matrix dimensions if `dist_cm` is a matrix (used by
# ng_recomb_decay_matrix()).
ng_meiosis_r <- function(dist_cm, model = c("haldane", "kosambi")) {
  model <- match.arg(model)
  d <- abs(dist_cm)
  storage.mode(d) <- "double"
  if (model == "haldane") 0.5 * (1 - exp(-2 * d / 100)) else 0.5 * tanh(2 * d / 100)
}

# Equilibrium between-marker recombination fraction R in F2-derived RILs by
# selfing (Haldane & Waddington 1931): R = 2r / (1 + 2r). Independent of the
# mapping function chosen for r itself.
ng_ril_R_from_r <- function(r) {
  r <- pmax(pmin(as.numeric(r), 0.5), 0)
  (2 * r) / (1 + 2 * r)
}

# Covariance kernel (1 - 2 R) for the target generation. For DH from F1 of two
# inbred parents this is (1 - 2 r); for selfed RILs it is (1 - 2 r) / (1 + 2 r).
# Note: the RIL kernel is NOT multiplicative in distance under Haldane (it is
# multiplicative under Kosambi, reducing to exp(-4 d / 100)).
ng_progeny_decay <- function(dist_cm, model = c("haldane", "kosambi"),
                             target = c("DH", "RIL")) {
  model <- match.arg(model)
  target <- match.arg(target)
  r <- ng_meiosis_r(dist_cm, model = model)
  if (identical(target, "DH")) 1 - 2 * r else (1 - 2 * r) / (1 + 2 * r)
}

ng_prepare_marker_map <- function(marker_map, marker_ids, model = c("haldane", "kosambi")) {
  model <- match.arg(model)
  marker_ids <- as.character(marker_ids)
  if (is.null(marker_map)) {
    marker_map <- data.frame(
      marker = marker_ids,
      chr = 1L,
      pos_cm = seq_along(marker_ids) - 1,
      stringsAsFactors = FALSE
    )
  } else if (is.vector(marker_map) && !is.data.frame(marker_map)) {
    if (is.null(names(marker_map))) names(marker_map) <- marker_ids
    marker_map <- data.frame(
      marker = names(marker_map),
      chr = 1L,
      pos_cm = as.numeric(marker_map),
      stringsAsFactors = FALSE
    )
  } else {
    marker_map <- as.data.frame(marker_map, stringsAsFactors = FALSE)
  }
  if (!("marker" %in% names(marker_map))) {
    marker_map$marker <- if (!is.null(rownames(marker_map))) rownames(marker_map) else marker_ids
  }
  if (!("chr" %in% names(marker_map))) marker_map$chr <- 1L
  if (!("pos_cm" %in% names(marker_map))) {
    pos_col <- intersect(c("pos", "position", "cm", "cM", "genetic_pos"), names(marker_map))
    if (!length(pos_col)) ng_stop("marker_map must contain pos_cm or a recognized position column")
    marker_map$pos_cm <- marker_map[[pos_col[1]]]
  }
  marker_map$marker <- as.character(marker_map$marker)
  marker_map$chr <- as.character(marker_map$chr)
  marker_map$pos_cm <- as.numeric(marker_map$pos_cm)
  marker_map <- marker_map[match(marker_ids, marker_map$marker), , drop = FALSE]
  if (anyNA(marker_map$marker)) ng_stop("marker_map is missing markers used in genotype/effect matrix")
  marker_map$chr_index <- as.integer(factor(marker_map$chr, levels = unique(marker_map$chr)))
  marker_map$recomb_model <- model
  marker_map
}

ng_recomb_decay_matrix <- function(marker_map, model = c("haldane", "kosambi"),
                                   target = c("DH", "RIL")) {
  # m x m matrix R with R[k, l] = (1 - 2 R_kl) for the requested progeny target
  # and mapping function (Haldane or Kosambi). Off-chromosome entries are 0 and
  # the diagonal is 1 by convention. Used in two places:
  #   1. Kosambi DH and any RIL request fall back to the dense O(m^2) path in
  #      R/03_metrics.R (the closed-form chromosome recursion only holds for
  #      multiplicative Haldane DH decay).
  #   2. Full off-diagonal posterior PMV in
  #      ng_dh_recomb_variance_pairs_full_posterior() uses R element-wise
  #      multiplied with Sigma_beta to extend the diagonal-only marker-effect
  #      uncertainty correction to the genomicMateSelectR formulation
  #      PMV = a'Ra + d'(R o Sigma_beta) d.
  model <- match.arg(model)
  target <- match.arg(target)
  d <- outer(marker_map$pos_cm, marker_map$pos_cm, "-")
  same_chr <- outer(marker_map$chr, marker_map$chr, "==")
  decay <- ng_progeny_decay(d, model = model, target = target)
  decay[!same_chr] <- 0
  diag(decay) <- 1
  rownames(decay) <- marker_map$marker
  colnames(decay) <- marker_map$marker
  decay
}

ng_recomb_decay_banded <- function(marker_map, model = c("haldane", "kosambi"),
                                   target = c("DH", "RIL"),
                                   window_cm = Inf) {
  # Sparse banded (COO upper-triangle) representation of the m x m progeny-decay
  # matrix R, restricted to same-chromosome marker pairs with |pos_i - pos_j|
  # <= window_cm. Stores the diagonal (always 1) plus all in-window off-diagonal
  # entries with i < j. Off-diagonal symmetry is implicit: callers reconstruct
  # 2 * x for i != j to obtain a' R a.
  #
  # Returns a list with integer i, j, double x (length = m + N_band entries),
  # plus m, window_cm, model, target for diagnostics. Uses base R only — no
  # Matrix dependency — so the kernel works in barebones environments. Cost is
  # O(m * k_window) where k_window is the average number of markers within
  # window_cm on the same chromosome.
  model <- match.arg(model)
  target <- match.arg(target)
  m <- nrow(marker_map)
  chr_idx <- as.integer(marker_map$chr_index)
  if (anyNA(chr_idx)) {
    chr_idx <- as.integer(factor(marker_map$chr, levels = unique(marker_map$chr)))
  }
  pos <- as.numeric(marker_map$pos_cm)
  # Assumes markers have been sorted by ng_sort_by_map() so that all markers on
  # the same chromosome are contiguous and ordered by position. Walk each
  # chromosome with an expanding right pointer and emit (i, j, decay) triplets
  # for j > i with positions within the window.
  if (!is.finite(window_cm) || window_cm < 0) window_cm <- Inf
  # Diagonal entries (always exactly 1, regardless of model/target).
  diag_i <- seq_len(m)
  diag_j <- diag_i
  diag_x <- rep(1.0, m)
  # Pre-allocate triplet buffers in chunks to avoid quadratic append cost.
  chunk_size <- max(m, 1024L)
  off_i <- integer(0); off_j <- integer(0); off_x <- numeric(0)
  buf_i <- integer(chunk_size); buf_j <- integer(chunk_size); buf_x <- numeric(chunk_size)
  buf_n <- 0L
  flush_buf <- function() {
    if (buf_n == 0L) return(invisible())
    off_i <<- c(off_i, buf_i[seq_len(buf_n)])
    off_j <<- c(off_j, buf_j[seq_len(buf_n)])
    off_x <<- c(off_x, buf_x[seq_len(buf_n)])
    buf_n <<- 0L
  }
  push <- function(i_vec, j_vec, x_vec) {
    n <- length(i_vec)
    if (n == 0L) return(invisible())
    if (buf_n + n > length(buf_i)) {
      flush_buf()
      if (n > length(buf_i)) {
        buf_i <<- integer(max(n, chunk_size))
        buf_j <<- integer(length(buf_i))
        buf_x <<- numeric(length(buf_i))
      }
    }
    idx <- (buf_n + 1L):(buf_n + n)
    buf_i[idx] <<- i_vec
    buf_j[idx] <<- j_vec
    buf_x[idx] <<- x_vec
    buf_n <<- buf_n + n
  }
  if (m > 0L) {
    # Iterate over chromosomes via run-length encoding of the (sorted) chr_idx.
    chr_run <- rle(chr_idx)
    starts <- cumsum(c(1L, head(chr_run$lengths, -1L)))
    for (b in seq_along(chr_run$lengths)) {
      s <- starts[b]
      e <- s + chr_run$lengths[b] - 1L
      if (e <= s) next  # no off-diagonal entries on a single-marker chromosome
      r_ptr <- s
      for (i in s:(e - 1L)) {
        if (r_ptr < i + 1L) r_ptr <- i + 1L
        # Advance r_ptr while the next marker is still within window of i.
        while (r_ptr <= e && (pos[r_ptr] - pos[i]) <= window_cm) r_ptr <- r_ptr + 1L
        j_hi <- r_ptr - 1L
        if (j_hi >= i + 1L) {
          js <- (i + 1L):j_hi
          dvec <- pos[js] - pos[i]
          xvec <- ng_progeny_decay(dvec, model = model, target = target)
          push(rep.int(i, length(js)), js, xvec)
        }
      }
    }
    flush_buf()
  }
  list(
    i = c(diag_i, off_i),
    j = c(diag_j, off_j),
    x = c(diag_x, off_x),
    m = m,
    window_cm = window_cm,
    model = model,
    target = target,
    n_offdiag = length(off_i)
  )
}

ng_sort_by_map <- function(geno, effects, beta_var = NULL, marker_map) {
  ord <- order(marker_map$chr_index, marker_map$pos_cm, marker_map$marker)
  list(
    geno = geno[, ord, drop = FALSE],
    effects = effects[ord],
    beta_var = if (is.null(beta_var)) rep(0, length(effects))[ord] else beta_var[ord],
    marker_map = marker_map[ord, , drop = FALSE]
  )
}
