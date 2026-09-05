ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# This test exercises the FULL runner (ng_run_cross_prediction -> ng_cp_pipeline stages ->
# ng_run_cp_output_files -> ng_write_cross_priority_workbook) with write_outputs = TRUE, twice:
# once with no checks at all (the most common real call pattern) and once with checks. Every
# other check-reference test either sets write_outputs = FALSE (so ng_run_cp_output_files never
# runs) or calls the workbook functions directly (bypassing the ctx/list2env staged-pipeline
# mechanism entirely) -- so none of them can catch a bug in how trait_check_reference is threaded
# from ctx into ng_run_cp_output_files.

set.seed(11)
n_p <- 12L; n_m <- 40L
ids <- paste0("P", seq_len(n_p))
markers <- paste0("m", seq_len(n_m))
gm <- matrix(2L * rbinom(n_p * n_m, 1, 0.4), nrow = n_p, dimnames = list(ids, markers))
genotype <- data.frame(id = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(marker = markers, chr = rep(1:4, length.out = n_m),
                         bp = rep(seq(0, 100, length.out = 10), 4)[seq_len(n_m)] * 1e6,
                         stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)
pheno <- data.frame(id = ids, yield = rnorm(n_p, 10, 2), stringsAsFactors = FALSE)

args <- list(genotype = genotype, phenotype = pheno, marker_map = marker_map,
             map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
             bp_per_cm = 1e6, id_col = "id",
             trait_direction = direction, n_crosses = 5L,
             write_outputs = TRUE, write_figures = FALSE, seed = 5L)

# (a) WITHOUT trait_checks -- must complete without error, and no Checks sheet.
out_dir_a <- tempfile("check_ref_outputs_a_")
res_a <- do.call(ng_run_cross_prediction, c(args, list(output_dir = out_dir_a)))
stopifnot(!is.null(res_a$output_files$workbook), file.exists(res_a$output_files$workbook))
sheets_a <- openxlsx::getSheetNames(res_a$output_files$workbook)
stopifnot(!("Checks" %in% sheets_a))

# (b) WITH trait_checks + check_geno + check_progeny_size -- must complete, and the workbook
# MUST contain a Checks sheet and the per-trait check columns on Selected_All.
chk <- matrix(2L * rbinom(2 * n_m, 1, 0.4), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), markers))
chk_records <- list(yield = list(adjusted_pheno = c(CHK_A = 11.5, CHK_B = 9.25)))
out_dir_b <- tempfile("check_ref_outputs_b_")
res_b <- do.call(ng_run_cross_prediction, c(args, list(
  output_dir = out_dir_b,
  check_geno = chk, check_progeny_size = 200L, check_records = chk_records,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(!is.null(res_b$output_files$workbook), file.exists(res_b$output_files$workbook))
sheets_b <- openxlsx::getSheetNames(res_b$output_files$workbook)
stopifnot("Checks" %in% sheets_b)

selected_b <- openxlsx::readWorkbook(res_b$output_files$workbook, sheet = "Selected_All", startRow = 3L)
stopifnot(all(c("yield_check_id", "yield_check_value", "yield_vs_check",
                "yield_check_ok", "yield_p_beat_check", "checks_all_ok") %in% names(selected_b)))

# M5: check_violation (also on Candidate_Crosses) and priority_check_component (Selected_All
# only -- priority ranking runs on the selected plan, not the full candidate table) must reach
# the workbook. Both exist on a check-enabled run regardless of priority_check_weight (default
# 0): check_violation is computed unconditionally by ng_attach_check_reference(), and
# ng_rank_cross_priority() always reports priority_check_component (0-weighted contribution)
# once ranking runs at all -- see tests/cross_priority.R.
stopifnot("check_violation" %in% names(selected_b))
stopifnot("priority_check_component" %in% names(selected_b))
candidate_b <- openxlsx::readWorkbook(res_b$output_files$workbook, sheet = "Candidate_Crosses", startRow = 3L)
stopifnot("check_violation" %in% names(candidate_b))

cat("check reference output-files wiring ok\n")

# (c) C2: check_geno supplied with NEITHER check_pheno NOR check_records (the vignette's
# original example, and the design's headline workflow) -- this fixture's ridge fit never
# calibrates reliability (R/02_effects.R), so mean_source resolves to the phenotypic
# "adjusted_pheno" for "yield", and with no phenotypic record supplied for the check, its value
# is NA for every cross. This must: (1) emit a diagnostic warning naming the trait and source,
# (2) NEVER report checks_all_ok = TRUE for a cross whose only check was never evaluated (it must
# be NA, not an affirmative pass), and (3) render "not evaluable" in the workbook, never "0 / N".
out_dir_c <- tempfile("check_ref_outputs_c_")
warned <- character(0)
res_c <- withCallingHandlers(
  do.call(ng_run_cross_prediction, c(args, list(
    output_dir = out_dir_c,
    check_geno = chk, check_progeny_size = 200L,
    trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(any(grepl("yield", warned) & grepl("not evaluable", warned)))
stopifnot(any(grepl("check_pheno", warned)))

ct_c <- res_c$candidate_crosses
stopifnot(all(is.na(ct_c$yield_check_value)))
stopifnot(all(is.na(ct_c$checks_all_ok)))   # never an affirmative TRUE for an unevaluated check

sheets_c <- openxlsx::getSheetNames(res_c$output_files$workbook)
stopifnot("Checks" %in% sheets_c)
checks_sheet_c <- openxlsx::readWorkbook(res_c$output_files$workbook, sheet = "Checks", startRow = 3L)
stopifnot(is.na(checks_sheet_c$value[[1L]]))
stopifnot(identical(checks_sheet_c$n_crosses_on_wrong_side[[1L]], "not evaluable"))

cat("check reference not-evaluable diagnostic + workbook wiring ok\n")

# --- M1: a default (single-check) run must still get a check visual -------------------------
# Previously the per-trait check panel was gated on nrow(trait_check_reference$active) > 1L
# (multi-trait only), and the main scatter's line is itself gated (D4, inside
# ng_check_line_value()) on trait_value_metric == "mean" -- so a run with exactly one active
# check and the package default trait_value_metric ("usefulness") got NEITHER visual. Reuse
# scenario (b)'s single yield check, this time with write_figures = TRUE.
out_dir_d <- tempfile("check_ref_outputs_d_")
args_d <- modifyList(args, list(write_figures = TRUE))
res_d <- do.call(ng_run_cross_prediction, c(args_d, list(
  output_dir = out_dir_d,
  check_geno = chk, check_progeny_size = 200L, check_records = chk_records,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(!is.null(res_d$output_files$check_panels_png))
stopifnot(file.exists(res_d$output_files$check_panels_png))

cat("M1 single-check panel ok\n")
