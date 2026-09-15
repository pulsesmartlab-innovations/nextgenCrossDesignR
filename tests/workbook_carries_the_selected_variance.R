# The workbook must carry the variance the run actually selected.
#
# A breeder picks ONE within-family variance -- uc_variance_source = "pmv" or "vpm" --
# and that choice drives the merit the plan is ranked on. The delivered workbook did
# not carry its value. It carried:
#
#   <trait>_usefulness   the merit  (built on the SELECTED variance)
#   <trait>_mid_parent   the mean
#   cross_upside         sqrt(VPM)  -- the OTHER variance, always, regardless of choice
#
# So on a default run (PMV) the sheet showed a merit derived from PMV beside a spread
# derived from VPM, with neither raw value present to reconcile them. PMV and VPM are
# different estimands -- PMV = VPM + tr(R Sigma_beta) -- and their ratio varies per
# cross (1.3x to 3.2x on a 40-line panel), so they can order crosses differently. A
# reader could not tell which quantity the ranking rested on, or check it.
#
# 0.32.0 added a Scoring_Method row NAMING the estimator. This adds the value, so the
# name can actually be checked against a number on the same sheet.
suppressPackageStartupMessages(library(openxlsx))
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

n <- 8L
crosses <- data.frame(
  parent1 = sprintf("P%02d", seq_len(n)), parent2 = sprintf("Q%02d", seq_len(n)),
  multi_trait_score = seq(0.9, 0.5, length.out = n),
  pair_kinship = seq(-0.2, 0.1, length.out = n),
  YIELD_value = seq(104, 97, length.out = n),
  YIELD_mean  = seq(103.5, 96.9, length.out = n),
  YIELD_pmv   = seq(3.0, 1.2, length.out = n),
  YIELD_vpm   = seq(2.0, 0.4, length.out = n),      # deliberately DIFFERENT from pmv
  stringsAsFactors = FALSE)
crosses <- ng_rank_cross_priority(crosses, breaks = c(0.25, 0.5, 0.75, 1), kinship_weight = 0.1)
td <- data.frame(trait = "YIELD", column = "YIELD_value", direction = "maximize",
                 weight = 1, stringsAsFactors = FALSE)
es <- function(var_col) data.frame(
  trait = "YIELD", variance_column_used = var_col,
  trait_value_metric_resolved = "usefulness", uc_variance_source_resolved = var_col,
  n_crosses_het_corrected = 0L, beta_var_available = TRUE, pmv_is_degenerate = FALSE,
  stringsAsFactors = FALSE)

tabs <- function(effect_summary) ng_cross_priority_workbook_tables(
  crosses = crosses, scored = crosses, trait_directions = td,
  n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE,
  trait_value_metric = "usefulness", trait_mean_source = list(YIELD = "GEBV"),
  effect_summary = effect_summary)

# ---- 1. a PMV run puts the PMV numbers on the sheet ------------------------
tp <- tabs(es("pmv"))
sheet <- tp[["Selected_All"]]
vcol <- grep("variance", names(sheet), ignore.case = TRUE, value = TRUE)
if (!length(vcol)) {
  stop("the workbook carries no within-family variance column, so the value the ",
       "ranking rests on is absent from the delivered file", call. = FALSE)
}
stopifnot(length(vcol) == 1L)
got <- suppressWarnings(as.numeric(sheet[[vcol]]))
stopifnot(isTRUE(all.equal(sort(got), sort(crosses$YIELD_pmv))))

# ---- 2. a VPM run puts the VPM numbers there instead -----------------------
# Driven by the run, not boilerplate: the same fixture with a different selection
# must show different numbers.
tv <- tabs(es("vpm"))
gotv <- suppressWarnings(as.numeric(tv[["Selected_All"]][[vcol]]))
stopifnot(isTRUE(all.equal(sort(gotv), sort(crosses$YIELD_vpm))))
stopifnot(!isTRUE(all.equal(sort(got), sort(gotv))))   # the choice actually matters

# ---- 3. the column names the trait, so a multi-trait sheet stays readable --
stopifnot(grepl("YIELD", vcol, fixed = TRUE))

# ---- 4. and it agrees with the Scoring_Method row on the same sheet --------
# The row names the estimator; the column carries its value. If they can disagree,
# the reader is back to guessing.
sm <- paste(unlist(tp[["Scoring_Method"]]), collapse = " ")
stopifnot(grepl("pmv", sm, fixed = TRUE))
smv <- paste(unlist(tv[["Scoring_Method"]]), collapse = " ")
stopifnot(grepl("vpm", smv, fixed = TRUE), !grepl("pmv", smv, fixed = TRUE))

# ---- 5. absent evidence -> no column, rather than an empty one -------------
bare <- ng_cross_priority_workbook_tables(
  crosses = crosses, scored = crosses, trait_directions = td,
  n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE)
stopifnot(!length(grep("variance", names(bare[["Selected_All"]]), ignore.case = TRUE)))

cat("workbook_carries_the_selected_variance: PASS  (column:", vcol, ")\n")
