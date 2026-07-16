# Small synthetic panel shared by the testthat suite.
#
# This suite is the fast gate that ships with the package and runs under
# R CMD check, so it exercises the exported API only (no ::: reach-ins) and
# stays deliberately small. The deep validation harness lives in the
# repository's tests/*.R scripts, which load the source tree via ng_load() and
# are not part of the installed package.

ng_test_panel <- function(n_parents = 20L, n_markers = 120L, n_chr = 3L, seed = 42L) {
  set.seed(seed)
  ids <- paste0("P", seq_len(n_parents))
  markers <- paste0("M", seq_len(n_markers))
  geno <- matrix(2L * rbinom(n_parents * n_markers, 1L, 0.45),
                 nrow = n_parents, ncol = n_markers,
                 dimnames = list(ids, markers))
  beta <- rnorm(n_markers, sd = 0.1)
  y <- as.numeric(geno %*% beta + rnorm(n_parents, sd = 1))
  names(y) <- ids
  per_chr <- ceiling(n_markers / n_chr)
  marker_map <- data.frame(
    marker = markers,
    chr = rep(seq_len(n_chr), each = per_chr)[seq_len(n_markers)],
    pos_cm = rep(seq(0, 100, length.out = per_chr), n_chr)[seq_len(n_markers)]
  )
  list(ids = ids, geno = geno, y = y, marker_map = marker_map)
}

# Scored crosses for the fixture panel, computed once per test file.
ng_test_scores <- function(p = ng_test_panel()) {
  fit <- ng_fit_ridge_effects(p$geno, p$y, ids = p$ids, kfold = 3L, seed = 1L)
  ng_score_crosses(p$geno, fit, marker_map = p$marker_map, ids = p$ids,
                   adjusted_pheno = p$y, selection_prop = 0.1, use_cpp = FALSE)
}
