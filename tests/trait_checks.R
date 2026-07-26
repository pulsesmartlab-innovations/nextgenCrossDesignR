ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- spec: direction default from trait_direction, explicit override, basis, validation ---
td <- c(yield = "increase", maturity = "decrease")
s <- ng_trait_check_spec(trait = c("yield", "maturity"), check = c("CkY", "CkM"),
                         trait_direction = td)
stopifnot(nrow(s) == 2L)
stopifnot(s$reject_if[s$trait == "yield"] == "below")     # increase -> reject if below check
stopifnot(s$reject_if[s$trait == "maturity"] == "above")  # decrease -> reject if above check
stopifnot(all(s$basis == "gebv"))                          # default basis
# explicit override + phenotype basis
s2 <- ng_trait_check_spec("protein", "CkP", direction = "above", basis = "phenotype",
                          trait_direction = c(protein = "increase"))
stopifnot(s2$reject_if == "above", s2$basis == "phenotype")
# invalid direction / basis error
err <- tryCatch(ng_trait_check_spec("x", "C", direction = "sideways"),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("direction", err))
# duplicate-trait error (Fix 3): v1 allows only one check per trait
errdup <- tryCatch(ng_trait_check_spec(c("yield", "yield"), c("CkA", "CkB"),
                                       trait_direction = c(yield = "increase")),
                    error = function(e) conditionMessage(e))
stopifnot(is.character(errdup), grepl("duplicate", errdup), grepl("yield", errdup))
cat("ng_trait_check_spec test passed\n")

# --- apply: gebv-basis violation, not-evaluable (NA), flag + optional exclude, diagnostics ---
scores <- data.frame(parent1 = c("A","A","B"), parent2 = c("B","C","C"),
                     stringsAsFactors = FALSE)
# yield increase -> reject if mid-parent BELOW check CkY (gebv 5)
gy <- c(A = 8, B = 6, C = 2, CkY = 5)
sp <- ng_trait_check_spec("yield", "CkY", trait_direction = c(yield = "increase"))
tv <- list(yield = list(gebv = gy, phenotype = NULL))
out <- ng_apply_trait_checks(scores, sp, tv, exclude = FALSE)
# mid-parents: A×B=7 (ok), A×C=5 (==check, NOT a violation), B×C=4 (<5 -> violation)
stopifnot(isFALSE(out$yield_check_violation[1]), isFALSE(out$yield_check_violation[2]),
          isTRUE(out$yield_check_violation[3]))
stopifnot(out$threshold_ok[1], out$threshold_ok[2], isFALSE(out$threshold_ok[3]))
stopifnot(out$threshold_violation[3] == "yield", out$threshold_violation[1] == "")
d <- attr(out, "trait_check_diagnostics")
stopifnot(d$n_flagged == 1L)
# exclude drops the violator
ex <- ng_apply_trait_checks(scores, sp, tv, exclude = TRUE)
stopifnot(nrow(ex) == 2L, attr(ex, "trait_check_diagnostics")$n_excluded == 1L)
# not-evaluable: check missing on basis -> NA violation, not a fail, counted
gy2 <- c(A = 8, B = 6, C = 2)          # no CkY
out2 <- ng_apply_trait_checks(scores, sp, list(yield = list(gebv = gy2, phenotype = NULL)))
stopifnot(all(is.na(out2$yield_check_violation)), all(out2$threshold_ok))
stopifnot(attr(out2, "trait_check_diagnostics")$n_not_evaluable >= 1L)
cat("ng_apply_trait_checks test passed\n")

# --- multi-trait apply: cross violating one trait vs. cross violating both, comma-joined labels ---
scores_mt <- data.frame(parent1 = c("A", "A", "B"), parent2 = c("B", "C", "C"),
                        stringsAsFactors = FALSE)
# yield increase -> reject if mid-parent BELOW check (CkY gebv = 5)
gy_mt <- c(A = 8, B = 6, C = 2, CkY = 5)
# protein increase -> reject if mid-parent BELOW check (CkP gebv = 5)
gp_mt <- c(A = 8, B = 8, C = 1, CkP = 5)
sp_mt <- ng_trait_check_spec(c("yield", "protein"), c("CkY", "CkP"),
                             trait_direction = c(yield = "increase", protein = "increase"))
tv_mt <- list(yield = list(gebv = gy_mt, phenotype = NULL),
             protein = list(gebv = gp_mt, phenotype = NULL))
out_mt <- ng_apply_trait_checks(scores_mt, sp_mt, tv_mt, exclude = FALSE)
# mid-parents -- yield: A×B=7(pass) A×C=5(==thr,pass) B×C=4(violation)
#              -- protein: A×B=8(pass) A×C=4.5(violation) B×C=4.5(violation)
stopifnot(out_mt$threshold_ok[1], out_mt$threshold_violation[1] == "")              # A×B: passes both
stopifnot(isFALSE(out_mt$threshold_ok[2]), out_mt$threshold_violation[2] == "protein")            # A×C: protein only
stopifnot(isFALSE(out_mt$threshold_ok[3]), out_mt$threshold_violation[3] == "yield,protein")      # B×C: both
cat("ng_apply_trait_checks multi-trait test passed\n")

# --- e2e: runner wiring (Task 3) -- flag columns + exclude behavior through the full pipeline ---
set.seed(9)
n <- 16L; mk <- 60L; gid <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * mk, 1, 0.5), n, mk, dimnames = list(gid, sprintf("M%03d", seq_len(mk))))
y  <- as.numeric(gm %*% rnorm(mk, 0, 0.1)) + rnorm(n)
genotype  <- data.frame(NAME = gid, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = gid, yield = y, stringsAsFactors = FALSE)
runmm <- data.frame(SNP = colnames(gm), chr = rep(1:2, length.out = mk),
                    bp = rep(seq(0, 100, length.out = 30), 2)[seq_len(mk)] * 1e6)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase")
# Elite check: the single BEST phenotype in the population is unusable here -- it is
# mathematically guaranteed to exceed every pairwise mid-parent average (the average of any two
# values can never exceed the true max), so "reject if below the top individual" would exclude
# every candidate cross and there would be nothing left to allocate. Use a strong-but-not-extreme
# individual (5th highest phenotype) so exclusion trims the pool without emptying it.
check_id <- gid[order(y, decreasing = TRUE)[5]]
checks <- data.frame(trait = "yield", check = check_id, stringsAsFactors = FALSE)
run <- function(excl) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, n_crosses = 8L, max_crosses_per_parent = 3L, use_ocs = TRUE,
  write_outputs = FALSE, write_figures = FALSE, seed = 5L,
  trait_checks = checks, exclude_threshold_violators = excl)
r0 <- run(FALSE)
stopifnot(all(c("yield_check_violation","threshold_ok","threshold_violation") %in% names(r0$selected_crosses)))
stopifnot(!is.null(r0$trait_check_diagnostics))
# with an elite check, many candidate crosses fall below it -> exclusion shrinks (but doesn't
# empty) the candidate pool
rE <- run(TRUE)
stopifnot(nrow(rE$candidate_crosses) > 0L)
stopifnot(all(rE$candidate_crosses$threshold_ok))         # survivors all pass
stopifnot(rE$trait_check_diagnostics$n_excluded > 0L)     # the check did trim something
cat("trait checks e2e test passed\n")

# --- Fix 2: a check line that is not a genotyped candidate parent must error clearly ---
bad_check_id <- "NOT_A_PARENT"
bad_checks <- data.frame(trait = "yield", check = bad_check_id, stringsAsFactors = FALSE)
err_check <- tryCatch({
  ng_run_cross_prediction(
    phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
    id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
    map_pos_cm_divisor = 1e6, n_crosses = 8L, max_crosses_per_parent = 3L, use_ocs = TRUE,
    write_outputs = FALSE, write_figures = FALSE, seed = 5L,
    trait_checks = bad_checks, exclude_threshold_violators = FALSE)
  NULL
}, error = function(e) conditionMessage(e))
stopifnot(is.character(err_check), grepl(bad_check_id, err_check, fixed = TRUE),
          grepl("yield", err_check, fixed = TRUE))
cat("trait checks check-not-a-parent validation test passed\n")

# --- json export: trait_check_diagnostics round-trip ---
env <- list(selected_crosses = r0$selected_crosses,
            trait_check_diagnostics = r0$trait_check_diagnostics)
sanitize <- function(x) { if (is.data.frame(x)) { for (j in seq_along(x))
  if (is.numeric(x[[j]])) x[[j]][!is.finite(x[[j]])] <- NA; return(x) }
  if (is.list(x)) return(lapply(x, sanitize)); if (is.numeric(x)) x[!is.finite(x)] <- NA; x }
tf <- tempfile(fileext = ".json")
jsonlite::write_json(sanitize(env), tf, auto_unbox = TRUE, na = "null", null = "null",
                     dataframe = "rows", digits = 8)
back <- jsonlite::read_json(tf, simplifyVector = TRUE)
stopifnot("threshold_ok" %in% names(back$selected_crosses))
stopifnot(!is.null(back$trait_check_diagnostics))
cat("trait checks json test passed\n")

# --- Excel opt-in (Task 5): include_trait_gebv adds a per-trait mid-parent GEBV column block ---
out_dir_on <- tempfile("ngcd_wb_on_")
r_wb_on <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, n_crosses = 8L, max_crosses_per_parent = 3L, use_ocs = TRUE,
  write_outputs = TRUE, write_figures = FALSE, seed = 5L,
  output_dir = out_dir_on, include_trait_gebv = TRUE)
stopifnot(!is.null(r_wb_on$output_files$workbook), file.exists(r_wb_on$output_files$workbook))

out_dir_off <- tempfile("ngcd_wb_off_")
r_wb_off <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, n_crosses = 8L, max_crosses_per_parent = 3L, use_ocs = TRUE,
  write_outputs = TRUE, write_figures = FALSE, seed = 5L,
  output_dir = out_dir_off)                          # include_trait_gebv defaults to FALSE
stopifnot(!is.null(r_wb_off$output_files$workbook), file.exists(r_wb_off$output_files$workbook))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  # Sheets carry a title row (row 1) before the real header row (row 3, see ng_cpw_write_sheet).
  sheet_on <- openxlsx::read.xlsx(r_wb_on$output_files$workbook, sheet = "Selected_All", startRow = 3L)
  stopifnot(any(grepl("_mid_parent_gebv$", names(sheet_on))))
  cand_on <- openxlsx::read.xlsx(r_wb_on$output_files$workbook, sheet = "Candidate_Crosses", startRow = 3L)
  stopifnot(any(grepl("_mid_parent_gebv$", names(cand_on))))

  sheet_off <- openxlsx::read.xlsx(r_wb_off$output_files$workbook, sheet = "Selected_All", startRow = 3L)
  stopifnot(!any(grepl("_mid_parent_gebv$", names(sheet_off))))
  cand_off <- openxlsx::read.xlsx(r_wb_off$output_files$workbook, sheet = "Candidate_Crosses", startRow = 3L)
  stopifnot(!any(grepl("_mid_parent_gebv$", names(cand_off))))
}
cat("include_trait_gebv workbook test passed\n")
