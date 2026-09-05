ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

crosses <- data.frame(
  parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
  yield_check_id = "CHK_A", yield_check_value = 6,
  yield_vs_check = c(4, -2), yield_check_ok = c(TRUE, FALSE),
  yield_p_beat_check = c(0.98, 0.71), checks_all_ok = c(TRUE, FALSE),
  stringsAsFactors = FALSE)

# n_candidates (66) is DELIBERATELY different from nrow(crosses) (2): n_wrong_side is counted
# once, over the full candidate pool, before allocation/selection narrows it down to whatever
# `crosses` table the workbook layer happens to pass in here (typically the much smaller
# selected plan). The denominator must come from diagnostics$n_candidates, NOT nrow(crosses) --
# if a future change reverted to nrow(crosses), this fixture's mismatch (2 vs 66) is what would
# catch it: the ratio would silently become "1 / 2" instead of the correct "1 / 66".
ref <- list(
  active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                      stringsAsFactors = FALSE),
  values = list(yield = c(CHK_A = 6)),
  source = c(yield = "GEBV"),
  # n_not_evaluable is PER TRAIT (a named list mirroring n_wrong_side), not a single scalar
  # aggregated across traits -- see R/51_check_reference.R::ng_attach_check_reference().
  diagnostics = list(n_wrong_side = list(yield = 1L), n_not_evaluable = list(yield = 0L),
                     n_candidates = 66L))

sh <- ng_cpw_checks_sheet(ref, crosses)
stopifnot(nrow(sh) == 1L)
stopifnot(sh$trait == "yield", sh$check_id == "CHK_A")
stopifnot(sh$direction == "want above")          # reject_if 'below' reads as 'want above'
stopifnot(sh$value == 6, sh$source == "GEBV")
# denominator is diagnostics$n_candidates (66), not nrow(crosses) (2, the selected-plan size)
stopifnot(sh$n_crosses_on_wrong_side == "1 / 66")

# a decrease trait reads the other way round
ref2 <- ref; ref2$active$reject_if <- "above"
stopifnot(ng_cpw_checks_sheet(ref2, crosses)$direction == "want below")

# --- C2 part 3: a check whose value never resolved renders "not evaluable", never "0 / N" ------
# "0 / 66" (0 crosses on the wrong side) is an affirmative claim -- it must be unreachable for a
# check that was never compared to anything. A never-evaluated check has n_wrong_side == 0 for
# the trivial reason that `ok` was NA for every cross (never FALSE), which is exactly what makes
# this case dangerous if not special-cased.
ref_na <- ref
ref_na$values$yield <- c(CHK_A = NA_real_)
ref_na$diagnostics$n_wrong_side$yield <- 0L
ref_na$diagnostics$n_not_evaluable$yield <- 66L
sh_na <- ng_cpw_checks_sheet(ref_na, crosses)
stopifnot(is.na(sh_na$value))
stopifnot(identical(sh_na$n_crosses_on_wrong_side, "not evaluable"))

# --- a check that DID resolve, but some individual crosses' own mean/variance did not, surfaces
# the partial gap alongside the wrong/total ratio rather than silently dropping it -------------
ref_partial <- ref
ref_partial$diagnostics$n_not_evaluable$yield <- 5L
sh_partial <- ng_cpw_checks_sheet(ref_partial, crosses)
stopifnot(identical(sh_partial$n_crosses_on_wrong_side, "1 / 66 (5 not evaluable)"))

cat("workbook checks sheet ok\n")
