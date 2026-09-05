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

# --- Task 3: reference columns, no row loss ---------------------------------
scores <- data.frame(
  parent1     = c("P1", "P1", "P2"),
  parent2     = c("P2", "P3", "P3"),
  yield_mean  = c(10, 4, 7),      # check at 6 -> above, below, above
  yield_pmv_used = c(4, 4, 0),    # sd = 2, 2, 0
  matur_mean  = c(70, 80, 75),    # check at 75, DECREASE -> below is good
  matur_pmv_used = c(1, 1, 1),
  stringsAsFactors = FALSE)

spec <- ng_trait_check_spec(trait = c("yield", "matur"), check = c("CHK_A", "CHK_B"),
                            trait_direction = c(yield = "increase", matur = "decrease"))
cv <- list(yield = c(CHK_A = 6), matur = c(CHK_B = 75))

out <- ng_attach_check_reference(scores, spec, trait_values = NULL, check_values = cv,
                                 k_progeny = 50L)

# ROW COUNT IS UNCHANGED -- this is the whole point
stopifnot(nrow(out) == 3L)
stopifnot(identical(out$parent1, scores$parent1))

# check id and value carried through
stopifnot(all(out$yield_check_id == "CHK_A"), all(out$yield_check_value == 6))

# direction-aware margin: POSITIVE ALWAYS MEANS BETTER, both directions
stopifnot(abs(out$yield_vs_check - c(4, -2, 1)) < 1e-8)      # increase: mean - check
stopifnot(abs(out$matur_vs_check - c(5, -5, 0)) < 1e-8)      # decrease: check - mean

# ok flags; a tie is NOT a violation
stopifnot(identical(out$yield_check_ok, c(TRUE, FALSE, TRUE)))
stopifnot(identical(out$matur_check_ok, c(TRUE, FALSE, TRUE)))   # row 3 is an exact tie

# combined flag
stopifnot(identical(out$checks_all_ok, c(TRUE, FALSE, TRUE)))

# P(beat check), increase trait: 1 - pnorm((tau-mu)/sd)^k
expect1 <- 1 - stats::pnorm((6 - 10) / 2)^50
stopifnot(abs(out$yield_p_beat_check[[1L]] - expect1) < 1e-10)
# a BELOW-check cross can still have a high tail probability -- the reason this column exists
expect2 <- 1 - stats::pnorm((6 - 4) / 2)^50
stopifnot(abs(out$yield_p_beat_check[[2L]] - expect2) < 1e-10, out$yield_p_beat_check[[2L]] > 0.5)
# sd == 0 degenerates to the indicator
stopifnot(out$yield_p_beat_check[[3L]] == 1)

# P(beat check), DECREASE trait: mirrored, 1 - pnorm((mu-tau)/sd)^k
expect_d <- 1 - stats::pnorm((70 - 75) / 1)^50
stopifnot(abs(out$matur_p_beat_check[[1L]] - expect_d) < 1e-10)

# a NOT-EVALUABLE check value yields NA columns and never a FALSE flag
out_na <- ng_attach_check_reference(scores, spec, NULL,
                                    list(yield = c(CHK_A = NA_real_), matur = c(CHK_B = 75)),
                                    k_progeny = 50L)
stopifnot(all(is.na(out_na$yield_check_value)), all(is.na(out_na$yield_check_ok)))
stopifnot(identical(out_na$checks_all_ok, c(TRUE, FALSE, TRUE)))   # driven by matur alone

d <- attr(out, "check_reference_diagnostics")
stopifnot(is.list(d), d$n_wrong_side$yield == 1L, d$n_wrong_side$matur == 1L)
stopifnot(d$n_not_evaluable == 0L)

# A trait whose name is not a valid R name: the cross table sanitises it for column names
# (make.names + dots to underscores) while the spec keeps what the breeder typed. Columns are
# looked up and written by the KEY; diagnostics and reporting use the RAW name.
scores_sp <- data.frame(parent1 = "P1", parent2 = "P2",
                        `Days_to_flower_mean` = 70, `Days_to_flower_pmv_used` = 1,
                        check.names = FALSE, stringsAsFactors = FALSE)
spec_sp <- ng_trait_check_spec("Days to flower", "CHK_B",
                               trait_direction = c(`Days to flower` = "decrease"))
spec_sp$column_key <- "Days_to_flower"
out_sp <- ng_attach_check_reference(scores_sp, spec_sp, NULL,
                                    list(`Days to flower` = c(CHK_B = 75)), k_progeny = 50L)
stopifnot("Days_to_flower_check_value" %in% names(out_sp))   # written by key
stopifnot(out_sp$Days_to_flower_check_ok)                     # 70 < 75, decrease -> good
stopifnot(names(attr(out_sp, "check_reference_diagnostics")$n_wrong_side) == "Days to flower")

cat("task 3 ok\n")

# --- Task 4: tau bounds for the multi-trait joint probability ---------------
b <- ng_check_tau_bounds(spec, cv)
stopifnot(identical(names(b$tau_lower), c("yield", "matur")))
# increase trait -> a LOWER bound at the check; decrease trait -> an UPPER bound
stopifnot(b$tau_lower[["yield"]] == 6, is.infinite(b$tau_upper[["yield"]]), b$tau_upper[["yield"]] > 0)
stopifnot(b$tau_upper[["matur"]] == 75, is.infinite(b$tau_lower[["matur"]]), b$tau_lower[["matur"]] < 0)

# an unevaluable check widens to an unbounded side rather than excluding the trait
b_na <- ng_check_tau_bounds(spec, list(yield = c(CHK_A = NA_real_), matur = c(CHK_B = 75)))
stopifnot(is.infinite(b_na$tau_lower[["yield"]]), b_na$tau_lower[["yield"]] < 0)

# the joint column: P(a progeny beats EVERY check at once)
j <- ng_attach_joint_check_probability(scores, spec, cv, k_progeny = 50L)
stopifnot("p_beat_all_checks" %in% names(j), nrow(j) == 3L)
stopifnot(all(j$p_beat_all_checks >= 0 & j$p_beat_all_checks <= 1, na.rm = TRUE))
# beating both checks can never be more likely than beating either one alone
single <- ng_attach_check_reference(scores, spec, NULL, cv, k_progeny = 50L)
stopifnot(all(j$p_beat_all_checks <= single$yield_p_beat_check + 1e-8))
stopifnot(all(j$p_beat_all_checks <= single$matur_p_beat_check + 1e-8))

cat("task 4 ok\n")

# --- Task 5: the spec builder, moved and basis-free ------------------------
td <- c(yield = "increase", maturity = "decrease")
s <- ng_trait_check_spec(trait = c("yield", "maturity"), check = c("CkY", "CkM"),
                         trait_direction = td)
stopifnot(nrow(s) == 2L, !("basis" %in% names(s)))
stopifnot(s$reject_if[s$trait == "yield"] == "below")     # increase -> wrong side is below
stopifnot(s$reject_if[s$trait == "maturity"] == "above")  # decrease -> wrong side is above
s2 <- ng_trait_check_spec("protein", "CkP", direction = "above",
                          trait_direction = c(protein = "increase"))
stopifnot(s2$reject_if == "above")                        # explicit override wins
err_d <- tryCatch(ng_trait_check_spec("x", "C", direction = "sideways"),
                  error = function(e) conditionMessage(e))
stopifnot(is.character(err_d), grepl("direction", err_d))
err_dup <- tryCatch(ng_trait_check_spec(c("yield", "yield"), c("CkA", "CkB"),
                                        trait_direction = c(yield = "increase")),
                    error = function(e) conditionMessage(e))
stopifnot(is.character(err_dup), grepl("duplicate", err_dup))
# basis is GONE: passing it must be an unused-argument error, not silently ignored
err_b <- tryCatch(ng_trait_check_spec("yield", "CkY", basis = "phenotype",
                                      trait_direction = c(yield = "increase")),
                  error = function(e) conditionMessage(e))
stopifnot(is.character(err_b), grepl("unused argument", err_b))

cat("task 5 ok\n")

# --- Task 10: check_pheno -> check_records ----------------------------------
cp <- data.frame(NAME = c("CHK_A", "CHK_B"),
                 yield = c(11.0, 9.5), matur = c(74, 71),
                 stringsAsFactors = FALSE)
rec <- ng_check_records_from_pheno(cp, id_col = "NAME",
                                   trait_columns = c(yield = "yield", matur = "matur"),
                                   source = "adjusted_pheno")
stopifnot(is.list(rec), setequal(names(rec), c("yield", "matur")))
stopifnot(names(rec$yield) == "adjusted_pheno")
stopifnot(abs(rec$yield$adjusted_pheno[["CHK_A"]] - 11.0) < 1e-8)
stopifnot(abs(rec$matur$adjusted_pheno[["CHK_B"]] - 71) < 1e-8)

# the source key follows the RUN, not the input's name
rec_blue <- ng_check_records_from_pheno(cp, "NAME", c(yield = "yield"), source = "BLUE")
stopifnot(names(rec_blue$yield) == "BLUE")

# a trait column absent from check_pheno yields NA, never a silent drop
rec_miss <- ng_check_records_from_pheno(cp, "NAME", c(protein = "protein"),
                                        source = "adjusted_pheno")
stopifnot(all(is.na(rec_miss$protein$adjusted_pheno)))

# a GEBV source consults markers, so check_pheno is not used: converter returns NULL
stopifnot(is.null(ng_check_records_from_pheno(cp, "NAME", c(yield = "yield"), source = "GEBV")))

cat("task 10 ok\n")
