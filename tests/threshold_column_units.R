helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# DEFECT 1 (see .superpowers/sdd/2026-09-03-check-reference-lines-backend/threshold-fix-report.md):
# a breeder's min_value/max_value are trait-unit thresholds, but the threshold comparison was
# tested against `traits$column` -- the RANKING metric, which under the default
# trait_value_metric = "usefulness" is `mean + i*sqrt(pmv)`, not the family mean. A cross whose
# predicted MEAN is below min_value can still clear the threshold because its inflated `_value`
# is above it.
#
# This case is built explicitly so it is impossible to satisfy by accident: row 1's mean (66) is
# below min_value (70), but its usefulness-inflated value (66 + 1.755*5 = 74.775) is comfortably
# above 70. Row 2 is a true negative control (mean 80, above min_value; value even higher) that
# must never violate, so a fix that over-corrects (e.g. always flags every row) is also caught.
scores <- data.frame(
  parent1 = c("A", "C"),
  parent2 = c("B", "D"),
  yield_value = c(66 + 1.755 * 5, 80 + 1.755 * 4),
  yield_mean = c(66, 80),
  stringsAsFactors = FALSE
)

traits <- data.frame(
  trait = "yield",
  column = "yield_value",
  threshold_column = "yield_mean",
  direction = "maximize",
  weight = 1,
  min_value = 70,
  max_value = NA_real_,
  stringsAsFactors = FALSE
)

scored <- ng_add_multitrait_score(scores, traits, method = "weighted")

violation_col <- "multi_trait_yield_violation"
stopifnot(violation_col %in% names(scored))

# Row 1: mean (66) is below min_value (70) -> MUST register a violation, even though the
# ranking value (74.775) clears the threshold. This is the exact scenario that silently passed
# before the fix (threshold compared against `column`, i.e. the usefulness value, instead of
# `threshold_column`, i.e. the mean).
stopifnot(scored[[violation_col]][[1L]] > 0)

# Row 2: mean (80) is above min_value (70) -> must NOT violate.
stopifnot(scored[[violation_col]][[2L]] == 0)

# Also exercise the wired-up caller path: ng_breeder_selection_objective +
# ng_score_breeder_objective must thread threshold_column through ng_multitrait_spec() end to
# end (this is the path R/39_cross_prediction_runner.R actually uses).
objective <- ng_breeder_selection_objective(
  trait = traits,
  method = "weighted",
  threshold_policy = "strict"
)
stopifnot(identical(objective$traits$threshold_column, "yield_mean"))
scored2 <- ng_score_breeder_objective(scores, objective)
stopifnot(scored2$multi_trait_score[[1L]] == -Inf)  # strict policy: violating row is excluded
stopifnot(is.finite(scored2$multi_trait_score[[2L]]))

# A trait spec that omits threshold_column entirely must keep behaving exactly as before: the
# comparison basis defaults to `column`. This is the caller-compatibility contract from the fix
# report -- direct callers who control both the spec and the scores table are unaffected.
traits_no_threshold_col <- traits
traits_no_threshold_col$threshold_column <- NULL
scored_default <- ng_add_multitrait_score(scores, traits_no_threshold_col, method = "weighted")
# With no threshold_column, the comparison falls back to `column` (yield_value), which never
# dips below 70 for either row here -- so neither row violates under this fallback spec.
stopifnot(scored_default[[violation_col]][[1L]] == 0)
stopifnot(scored_default[[violation_col]][[2L]] == 0)

# A threshold_column naming a column absent from `scores` is a programming error and must fail
# loudly, naming the trait and the missing column.
traits_bad <- traits
traits_bad$threshold_column <- "yield_mean_typo"
err <- tryCatch(
  ng_add_multitrait_score(scores, traits_bad, method = "weighted"),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err))
stopifnot(grepl("yield", err, fixed = TRUE))
stopifnot(grepl("yield_mean_typo", err, fixed = TRUE))

cat("threshold_column units tests passed\n")
