ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- Task 1: check-genotype alignment ---------------------------------------
markers <- c("m1", "m2", "m3")
chk <- matrix(c(0, 2, 2,
                2, 0, 0), nrow = 2, byrow = TRUE,
              dimnames = list(c("CHK_A", "CHK_B"), c("m3", "m1", "m2")))

# columns are reordered to the parents' marker order, rownames preserved
al <- ng_align_check_geno(chk, markers)
stopifnot(identical(colnames(al), markers))
stopifnot(identical(rownames(al), c("CHK_A", "CHK_B")))
stopifnot(al["CHK_A", "m1"] == 2, al["CHK_A", "m3"] == 0)

# extra markers in the check file are dropped, not an error
chk_extra <- cbind(chk, m9 = c(1, 1))
al2 <- ng_align_check_geno(chk_extra, markers)
stopifnot(identical(colnames(al2), markers))

# a MISSING marker is a hard error naming the count
err <- tryCatch(ng_align_check_geno(chk[, c("m1", "m2"), drop = FALSE], markers),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("m3", err), grepl("1 of 3", err))

# unnamed columns are a hard error (cannot verify alignment)
bare <- matrix(0, nrow = 1, ncol = 3, dimnames = list("CHK_A", NULL))
err2 <- tryCatch(ng_align_check_geno(bare, markers), error = function(e) conditionMessage(e))
stopifnot(is.character(err2), grepl("column names", err2))

# dosages outside [0, ploidy] are a hard error (wrong coding / flipped file)
badcode <- matrix(c(0, 5, 1), nrow = 1, dimnames = list("CHK_A", markers))
err3 <- tryCatch(ng_align_check_geno(badcode, markers), error = function(e) conditionMessage(e))
stopifnot(is.character(err3), grepl("dosage", err3))

cat("task 1 ok\n")
