ng_polyploid_as_dosage_matrix <- function(geno, ploidy = 4L, name = deparse(substitute(geno))) {
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
  geno <- ng_polyploid_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  (geno - ploidy / 2) * (2 / ploidy)
}

ng_poly4x_digenic_scaled <- function(geno, ploidy = 4L) {
  geno <- ng_polyploid_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  geno * (ploidy - geno) * (2 / ploidy) ^ 2
}

# DEPRECATED -- do not use for coancestry, OCS, or any relatedness decision.
#
# This centers dosage on the ploidy MIDPOINT, which implicitly assumes every allele frequency is
# 0.5. It is therefore an allele-frequency-dominated IBS similarity, not an additive relationship
# matrix. Measured on 200 UNRELATED autotetraploids (independent draws from common frequencies,
# i.e. true relatedness ~ 0) it returns a mean off-diagonal of 1.47 (range 1.20..1.72) where the
# correct GRM returns -0.005 (range -0.15..0.15), and mean diagonal 2.47 versus 0.995. It is not
# merely rescaled: the RANK correlation of pair values against the correct matrix is only ~0.56,
# so any coancestry ORDERING taken from it is wrong, which is what an OCS penalty consumes.
#
# Use ng_polyploid_grm() (R/47; VanRaden/Yang generalized to ploidy) instead. Retained only so
# existing calls keep resolving, and to keep the regression test that pins the contrast between
# the two matrices.
ng_poly4x_parent_relationship <- function(geno, ploidy = 4L) {
  geno <- ng_polyploid_as_dosage_matrix(geno, ploidy = ploidy, name = "geno")
  X <- ng_poly4x_additive_scaled(geno, ploidy = ploidy)
  denom <- sum(diag(stats::var(X)))
  if (!is.finite(denom) || denom <= 0) denom <- ncol(X)
  K <- tcrossprod(X) / denom
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  K
}

ng_poly4x_pair_relationship <- function(parent_kinship, pairs) {
  parent_kinship <- as.matrix(parent_kinship)
  if (!is.numeric(parent_kinship)) ng_stop("parent_kinship must be numeric")
  if (nrow(parent_kinship) != ncol(parent_kinship)) ng_stop("parent_kinship must be a square matrix")
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("pairs must contain parent1 and parent2")
  }
  row_ids <- rownames(parent_kinship)
  col_ids <- colnames(parent_kinship)
  if (is.null(row_ids) || is.null(col_ids)) ng_stop("parent_kinship must have row and column dimnames")
  if (anyDuplicated(row_ids)) ng_stop("parent_kinship row IDs must be unique")
  if (anyDuplicated(col_ids)) ng_stop("parent_kinship column IDs must be unique")
  requested <- unique(c(as.character(pairs$parent1), as.character(pairs$parent2)))
  miss <- setdiff(requested, row_ids)
  if (length(miss)) ng_stop("parent_kinship missing parent IDs: ", paste(miss, collapse = ", "))
  miss <- setdiff(requested, col_ids)
  if (length(miss)) ng_stop("parent_kinship column IDs missing parent IDs: ", paste(miss, collapse = ", "))
  if (!setequal(row_ids, col_ids)) ng_stop("parent_kinship row and column IDs must describe the same parents")
  vapply(seq_len(nrow(pairs)), function(i) {
    parent_kinship[as.character(pairs$parent1[[i]]), as.character(pairs$parent2[[i]])]
  }, numeric(1))
}

ng_poly4x_pair_coancestry <- function(parent_kinship, pairs) {
  ng_poly4x_pair_relationship(parent_kinship, pairs) / 2
}
