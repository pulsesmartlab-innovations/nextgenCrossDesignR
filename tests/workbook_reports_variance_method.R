# The workbook must state which within-family variance produced its numbers.
#
# 0.31.0 added a "Cross-mean basis (this run)" row because a breeder could not tell
# whether a cross mean was genomic or phenotypic. The variance side of the same
# question was left unanswered, and it is the harder half:
#
#   * seven variance-ish columns ship on the cross table, four of them numerically
#     identical on a typical run (pmv, pmv_fast, pmv_used, var_complex), one all-NA
#     (pmv_full_posterior) and one ~30x smaller (vpm);
#   * uc_variance_source appeared nowhere in the workbook at all, so PMV-vs-VPM was
#     unreportable even though fast-vs-full PMV was reported;
#   * and the het-parent correction silently substitutes a DIFFERENT estimator for
#     a subset of rows, under a truncation regime (dense, window_cm inapplicable)
#     that differs from every other row in the same table.
#
# A spreadsheet outlives the session that produced it. If the variance method is
# not on the sheet, it is not reported -- and the reader is left comparing numbers
# that may not be comparable.
#
# The row is driven by data, not boilerplate: a VPM run must not describe itself as
# PMV, and the het-correction sentence must appear only when rows were actually
# corrected.
suppressPackageStartupMessages(library(openxlsx))
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

n <- 8L
crosses <- data.frame(
  parent1 = sprintf("P%02d", seq_len(n)),
  parent2 = sprintf("Q%02d", seq_len(n)),
  multi_trait_score = seq(0.9, 0.5, length.out = n),
  pair_kinship = seq(-0.2, 0.1, length.out = n),
  YIELD_value = seq(104, 97, length.out = n),
  YIELD_mean = seq(103.5, 96.9, length.out = n),
  YIELD_mean_gebv = seq(97.0, 97.3, length.out = n),
  stringsAsFactors = FALSE)
crosses <- ng_rank_cross_priority(crosses, breaks = c(0.25, 0.5, 0.75, 1),
                                  kinship_weight = 0.1)
td <- data.frame(trait = "YIELD", column = "YIELD_value", direction = "maximize",
                 weight = 1, stringsAsFactors = FALSE)

es <- function(var_col, het_n = 0L, requires = TRUE) {
  data.frame(trait = "YIELD", variance_column_used = var_col,
             trait_value_metric_resolved = "usefulness",
             uc_variance_source_resolved = var_col,
             n_crosses_het_corrected = het_n,
             beta_var_available = TRUE, pmv_is_degenerate = FALSE,
             stringsAsFactors = FALSE)
}

sheet_text <- function(effect_summary) {
  t <- ng_cross_priority_workbook_tables(
    crosses = crosses, scored = crosses, trait_directions = td,
    n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE,
    trait_value_metric = "usefulness",
    trait_mean_source = list(YIELD = "adjusted_pheno"),
    effect_summary = effect_summary)
  paste(unlist(t[["Scoring_Method"]]), collapse = " | ")
}

# ---- 1. the variance method is named, and the trait it applies to ----------
pmv <- sheet_text(es("pmv"))
stopifnot(grepl("variance", pmv, ignore.case = TRUE))
stopifnot(grepl("YIELD", pmv, fixed = TRUE))
stopifnot(grepl("pmv", pmv, fixed = TRUE))

# ---- 2. it is driven by the run, not boilerplate ---------------------------
# A VPM run must describe itself as VPM and must NOT claim PMV.
vpm <- sheet_text(es("vpm"))
stopifnot(grepl("vpm", vpm, fixed = TRUE))
stopifnot(!grepl("pmv", vpm, fixed = TRUE))

# ---- 3. marker dependence is stated ---------------------------------------
# The fact this whole release turns on: every within-family variance here is a
# quadratic form in beta-hat, so a reader must be told the variance inherits the
# marker model's quality even when the MEAN did not.
stopifnot(grepl("marker", pmv, ignore.case = TRUE))

# ---- 4. the het-correction sentence appears only when rows were corrected --
none <- sheet_text(es("pmv", het_n = 0L))
some <- sheet_text(es("pmv", het_n = 3L))
stopifnot(!grepl("phased", none, ignore.case = TRUE))
stopifnot(grepl("phased", some, ignore.case = TRUE))
stopifnot(grepl("3", some, fixed = TRUE))          # names the count
stopifnot(grepl("window_cm", some, fixed = TRUE))  # and that the window does not apply

# ---- 5. absent evidence produces no row, rather than an empty one ---------
# Same rule the mean-basis row follows: a fixed glossary would put claims on a
# sheet for a run that never made them.
bare <- ng_cross_priority_workbook_tables(
  crosses = crosses, scored = crosses, trait_directions = td,
  n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE)
bt <- paste(unlist(bare[["Scoring_Method"]]), collapse = " | ")
stopifnot(!grepl("Within-family variance method", bt, fixed = TRUE))

cat("workbook_reports_variance_method: PASS\n")
