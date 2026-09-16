# The delivered sheet must say whether the safety was on.
#
# The workbench is about to expose `effect_gate` to breeders as a plain on/off control,
# deliberately, so a breeder can decide from the outset what they are willing to accept.
# That decision has to survive into the file someone else reads. Today it does not:
# effect_gate reaches result.json, but `grep effect_gate R/37_cross_priority_workbook.R`
# returns nothing, so a plan built on markers the engine would have REFUSED looks
# exactly like a plan built on markers it accepted.
#
# That is the unlabelled-number defect 0.31.0 and 0.32.0 spent two releases removing,
# in its most consequential form yet: not a mislabelled basis, but a silently disabled
# refusal.
#
# The row also carries the two things a breeder cannot reconstruct from the numbers:
#   * the threshold ON ITS OWN SCALE. cv_predictive_r2 is an out-of-fold R2 against
#     PHENOTYPE. Breeders reason in accuracy, and 0.35 R2 is r ~ 0.59 -- a demanding
#     bar, not a modest one. A bare "0.35" invites the opposite reading.
#   * which traits fell back, AND that their variance stayed marker-derived. The
#     fallback replaces the MEAN only; a mixed-basis plan is unreadable without that.
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
  stringsAsFactors = FALSE)
crosses <- ng_rank_cross_priority(crosses, breaks = c(0.25, 0.5, 0.75, 1), kinship_weight = 0.1)
td <- data.frame(trait = "YIELD", column = "YIELD_value", direction = "maximize",
                 weight = 1, stringsAsFactors = FALSE)

es <- function(gate = "on", thresh = 0.35, mean_source = "GEBV", traits = "YIELD",
               cv = 0.52) data.frame(
  trait = traits, variance_column_used = "pmv",
  trait_value_metric_resolved = "usefulness", uc_variance_source_resolved = "pmv",
  mean_source = mean_source, cv_predictive_r2 = cv,
  effect_gate = gate, min_cv_predictive_r2_applied = thresh,
  n_crosses_het_corrected = 0L, beta_var_available = TRUE, pmv_is_degenerate = FALSE,
  stringsAsFactors = FALSE)

sheet_text <- function(effect_summary) {
  t <- ng_cross_priority_workbook_tables(
    crosses = crosses, scored = crosses, trait_directions = td,
    n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE,
    trait_value_metric = "usefulness", trait_mean_source = list(YIELD = "GEBV"),
    effect_summary = effect_summary)
  paste(unlist(t[["Scoring_Method"]]), collapse = " | ")
}

# ---- 1. a gate-ON run says so, and states the threshold with its scale ----
on <- sheet_text(es(gate = "on", thresh = 0.35))
stopifnot(grepl("gate", on, ignore.case = TRUE))
stopifnot(grepl("0.35", on, fixed = TRUE))
# the correlation equivalent, so the bar is read on the scale breeders use
stopifnot(grepl("0.59", on, fixed = TRUE))
stopifnot(grepl("R2|R²", on))

# ---- 2. a gate-OFF run says the refusal was bypassed ----------------------
# This is the case the row exists for. It must be unmistakable, not a footnote.
off <- sheet_text(es(gate = "off"))
stopifnot(grepl("OFF", off))
stopifnot(grepl("refus", off, ignore.case = TRUE))
stopifnot(!identical(on, off))          # driven by the run, not boilerplate

# ---- 3. a fallback trait is named, with the asymmetry stated --------------
fb <- sheet_text(es(gate = "on", mean_source = "adjusted_pheno", cv = 0.11))
stopifnot(grepl("YIELD", fb, fixed = TRUE))
stopifnot(grepl("phenotyp", fb, ignore.case = TRUE))
# the variance did NOT fall back, and a mixed-basis plan is unreadable without this
stopifnot(grepl("variance", fb, ignore.case = TRUE))
stopifnot(grepl("marker", fb, ignore.case = TRUE))
# a run where nothing fell back must NOT claim a fallback
stopifnot(!grepl("fell back", on, ignore.case = TRUE))

# ---- 4. the threshold shown is the one that APPLIED ----------------------
# A breeder who lowered the bar must see the bar they set, not the package default.
lowered <- sheet_text(es(gate = "on", thresh = 0.10))
stopifnot(grepl("0.10|0.1\\b", lowered))
stopifnot(!grepl("0.35", lowered, fixed = TRUE))

# ---- 5. absent evidence -> no row, rather than an empty one --------------
bare <- ng_cross_priority_workbook_tables(
  crosses = crosses, scored = crosses, trait_directions = td,
  n_crosses_requested = n, block_size = 4L, include_trait_gebv = TRUE)
stopifnot(!grepl("gate", paste(unlist(bare[["Scoring_Method"]]), collapse = " "),
                 ignore.case = TRUE))

# ---- 6. an unreachable threshold does not print an impossible correlation --
# R2 <= 1 always, so sqrt() of a threshold above 1 yields "r = 1.05" -- a statistic
# that cannot exist. A threshold above 1 is a legitimate way to force phenotypic means
# (no model can clear it), so the row says that rather than inventing a number. Found
# by reading a real workbook, which is the only place it would have shown up.
imp <- sheet_text(es(gate = "on", thresh = 1.1, mean_source = "adjusted_pheno"))
stopifnot(!grepl("r = 1[.]", imp))
stopifnot(!grepl("r = [2-9]", imp))
stopifnot(grepl("cannot exceed 1", imp, fixed = TRUE))
stopifnot(grepl("by construction", imp, fixed = TRUE))
# and a reachable threshold still gets its correlation
stopifnot(grepl("r = 0.45", sheet_text(es(gate = "on", thresh = 0.20)), fixed = TRUE))

cat("workbook_reports_the_gate: PASS\n")
