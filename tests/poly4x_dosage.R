helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

dosage <- matrix(
  c(0, 1, 2, 3, 4,
    4, 3, 2, 1, 0),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("P1", "P2"), paste0("m", 1:5))
)

additive <- ng_poly4x_additive_scaled(dosage)
expected_additive <- matrix(
  c(-1.0, -0.5, 0.0, 0.5, 1.0,
     1.0,  0.5, 0.0, -0.5, -1.0),
  nrow = 2,
  byrow = TRUE,
  dimnames = dimnames(dosage)
)
stopifnot(isTRUE(all.equal(additive, expected_additive, tolerance = 1e-12)))

digenic <- ng_poly4x_digenic_scaled(dosage)
expected_digenic <- matrix(
  c(0.00, 0.75, 1.00, 0.75, 0.00,
    0.00, 0.75, 1.00, 0.75, 0.00),
  nrow = 2,
  byrow = TRUE,
  dimnames = dimnames(dosage)
)
stopifnot(isTRUE(all.equal(digenic, expected_digenic, tolerance = 1e-12)))

checked <- ng_poly4x_as_dosage_matrix(dosage)
stopifnot(identical(dimnames(checked), dimnames(dosage)))
stopifnot(storage.mode(checked) == "double")

bad <- dosage
bad[1, 1] <- 5
err <- tryCatch(ng_poly4x_as_dosage_matrix(bad), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("outside 0..4", conditionMessage(err), fixed = TRUE))

err <- tryCatch(ng_poly4x_as_dosage_matrix(dosage, ploidy = 4.9), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("ploidy must be an integer >= 2", conditionMessage(err), fixed = TRUE))

dup_dosage <- dosage
rownames(dup_dosage) <- c("P1", "P1")
err <- tryCatch(ng_poly4x_as_dosage_matrix(dup_dosage, name = "geno"), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("geno row names must be unique", conditionMessage(err), fixed = TRUE))

K <- ng_poly4x_parent_relationship(dosage)
expected_K <- matrix(
  c(0.5, -0.5,
    -0.5, 0.5),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(rownames(dosage), rownames(dosage))
)
stopifnot(is.matrix(K))
stopifnot(identical(rownames(K), rownames(dosage)))
stopifnot(identical(colnames(K), rownames(dosage)))
stopifnot(isTRUE(all.equal(K, expected_K, tolerance = 1e-12)))

pairs <- data.frame(parent1 = c("P1", "P1"), parent2 = c("P2", "P1"), stringsAsFactors = FALSE)
co <- ng_poly4x_pair_coancestry(K, pairs)
stopifnot(length(co) == 2L)
stopifnot(isTRUE(all.equal(co, K[cbind(pairs$parent1, pairs$parent2)], tolerance = 1e-12)))

K_dup <- K
rownames(K_dup) <- c("P1", "P1")
err <- tryCatch(ng_poly4x_pair_coancestry(K_dup, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K row IDs must be unique", conditionMessage(err), fixed = TRUE))

K_dup <- K
colnames(K_dup) <- c("P1", "P1")
err <- tryCatch(ng_poly4x_pair_coancestry(K_dup, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K column IDs must be unique", conditionMessage(err), fixed = TRUE))

K_no_dimnames <- unname(K)
err <- tryCatch(ng_poly4x_pair_coancestry(K_no_dimnames, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K must have row and column dimnames", conditionMessage(err), fixed = TRUE))

K_non_numeric <- matrix("x", nrow = 2, ncol = 2, dimnames = dimnames(K))
err <- tryCatch(ng_poly4x_pair_coancestry(K_non_numeric, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K must be numeric", conditionMessage(err), fixed = TRUE))

K_not_square <- K[, 1, drop = FALSE]
err <- tryCatch(ng_poly4x_pair_coancestry(K_not_square, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K must be a square matrix", conditionMessage(err), fixed = TRUE))

K_bad_cols <- K
colnames(K_bad_cols) <- c("P1", "PX")
err <- tryCatch(ng_poly4x_pair_coancestry(K_bad_cols, pairs), error = function(e) e)
stopifnot(inherits(err, "error"))
stopifnot(grepl("parent_K column IDs", conditionMessage(err), fixed = TRUE))

cat("poly4x dosage tests passed\n")
