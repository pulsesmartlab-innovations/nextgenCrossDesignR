ng_poly4x_as_dosage_matrix <- function(geno, ploidy = 4L, name = deparse(substitute(geno))) {
  if (length(ploidy) != 1L) ng_stop("ploidy must be an integer >= 2")
  ploidy <- suppressWarnings(as.numeric(ploidy))
  if (!is.finite(ploidy) || ploidy < 2L || abs(ploidy - round(ploidy)) > 1e-8) {
    ng_stop("ploidy must be an integer >= 2")
  }
  ploidy <- as.integer(ploidy)
  geno <- ng_as_numeric_matrix(geno, name)
  if (anyDuplicated(rownames(geno))) ng_stop(name, " row names must be unique")
  bad <- !is.finite(geno) | geno < 0 | geno > ploidy | abs(geno - round(geno)) > 1e-8
  if (any(bad)) {
    idx <- which(bad, arr.ind = TRUE)[1, ]
    ng_stop(
      name, " has dosage outside 0..", ploidy,
      " at row ", rownames(geno)[idx[[1]]],
      ", column ", colnames(geno)[idx[[2]]]
    )
  }
  storage.mode(geno) <- "double"
  geno
}

ng_poly4x_additive_scaled <- function(geno, ploidy = 4L) {
  geno <- ng_poly4x_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  (geno - ploidy / 2) * (2 / ploidy)
}

ng_poly4x_digenic_scaled <- function(geno, ploidy = 4L) {
  geno <- ng_poly4x_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  geno * (ploidy - geno) * (2 / ploidy) ^ 2
}

ng_poly4x_parent_relationship <- function(geno, ploidy = 4L) {
  geno <- ng_poly4x_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  X <- ng_poly4x_additive_scaled(geno, ploidy = ploidy)
  denom <- sum(diag(stats::var(X)))
  if (!is.finite(denom) || denom <= 0) denom <- ncol(X)
  K <- tcrossprod(X) / denom
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  K
}

ng_poly4x_pair_coancestry <- function(parent_K, pairs) {
  parent_K <- as.matrix(parent_K)
  if (!is.numeric(parent_K)) ng_stop("parent_K must be numeric")
  if (nrow(parent_K) != ncol(parent_K)) ng_stop("parent_K must be a square matrix")
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("pairs must contain parent1 and parent2")
  }
  row_ids <- rownames(parent_K)
  col_ids <- colnames(parent_K)
  if (is.null(row_ids) || is.null(col_ids)) ng_stop("parent_K must have row and column dimnames")
  if (anyDuplicated(row_ids)) ng_stop("parent_K row IDs must be unique")
  if (anyDuplicated(col_ids)) ng_stop("parent_K column IDs must be unique")
  requested <- unique(c(as.character(pairs$parent1), as.character(pairs$parent2)))
  miss <- setdiff(requested, row_ids)
  if (length(miss)) ng_stop("parent_K missing parent IDs: ", paste(miss, collapse = ", "))
  miss <- setdiff(requested, col_ids)
  if (length(miss)) ng_stop("parent_K column IDs missing parent IDs: ", paste(miss, collapse = ", "))
  if (!setequal(row_ids, col_ids)) ng_stop("parent_K row and column IDs must describe the same parents")
  vapply(seq_len(nrow(pairs)), function(i) {
    parent_K[as.character(pairs$parent1[[i]]), as.character(pairs$parent2[[i]])]
  }, numeric(1))
}
