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

# --- Task 2: value resolution follows the run's mean source -----------------
eff <- list(beta = c(m1 = 1, m2 = 0.5, m3 = -2), intercept = 0)
al <- ng_align_check_geno(chk, markers)   # CHK_A = (m1=2, m2=2, m3=0); CHK_B = (0, 0, 2)

# GEBV source -> predicted from the check's own markers with the SAME effects
gv <- ng_check_reference_value("GEBV", al, eff)
stopifnot(is.numeric(gv), identical(names(gv), c("CHK_A", "CHK_B")))
stopifnot(abs(gv[["CHK_A"]] - (2 * 1 + 2 * 0.5 + 0 * -2)) < 1e-8)
stopifnot(abs(gv[["CHK_B"]] - (0 * 1 + 0 * 0.5 + 2 * -2)) < 1e-8)

# the low-reliability / uncalibrated GEBV variants are still GEBV
stopifnot(abs(ng_check_reference_value("GEBV_uncalibrated", al, eff)[["CHK_A"]] - gv[["CHK_A"]]) < 1e-8)
stopifnot(abs(ng_check_reference_value("GEBV_low_reliability", al, eff)[["CHK_A"]] - gv[["CHK_A"]]) < 1e-8)

# a phenotypic source reads the supplied record, NOT the markers
bv <- ng_check_reference_value("BLUE", al, eff,
                               check_records = list(BLUE = c(CHK_A = 11.5, CHK_B = 9.25)))
stopifnot(abs(bv[["CHK_A"]] - 11.5) < 1e-8, abs(bv[["CHK_B"]] - 9.25) < 1e-8)

# NOT-EVALUABLE: run source is BLUE, the check has no BLUE record -> NA, never a GEBV fallback
nv <- ng_check_reference_value("BLUE", al, eff, check_records = list(BLUE = c(CHK_A = 11.5)))
stopifnot(abs(nv[["CHK_A"]] - 11.5) < 1e-8, is.na(nv[["CHK_B"]]))

# no records at all on a phenotypic source -> all NA (still not a GEBV fallback)
allna <- ng_check_reference_value("adjusted_pheno", al, eff)
stopifnot(all(is.na(allna)))

cat("task 1 ok\n")
cat("task 2 ok\n")
