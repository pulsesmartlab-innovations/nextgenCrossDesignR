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

ref <- list(
  active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                      stringsAsFactors = FALSE),
  values = list(yield = c(CHK_A = 6)),
  source = c(yield = "GEBV"),
  diagnostics = list(n_wrong_side = list(yield = 1L), n_not_evaluable = 0L, n_candidates = 2L))

sh <- ng_cpw_checks_sheet(ref, crosses)
stopifnot(nrow(sh) == 1L)
stopifnot(sh$trait == "yield", sh$check_id == "CHK_A")
stopifnot(sh$direction == "want above")          # reject_if 'below' reads as 'want above'
stopifnot(sh$value == 6, sh$source == "GEBV")
stopifnot(sh$n_crosses_on_wrong_side == "1 / 2")

# a decrease trait reads the other way round
ref2 <- ref; ref2$active$reject_if <- "above"
stopifnot(ng_cpw_checks_sheet(ref2, crosses)$direction == "want below")

cat("workbook checks sheet ok\n")
