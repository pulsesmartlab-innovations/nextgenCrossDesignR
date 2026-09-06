helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# DEFECT 2 (see .superpowers/sdd/2026-09-03-check-reference-lines-backend/threshold-fix-report.md):
# ng_preflight_input_tables() (R/26_data_preflight.R) accepts min/minimum as aliases for
# min_value (and max/maximum for max_value) -- it even cross-validates min <= max and reports a
# positive "your thresholds are consistent" result. But ng_multitrait_spec() (R/19) recognised
# only the literal names "min_value"/"max_value": any other name was silently overwritten with
# NA. So a spec using min/max passed preflight and then had its thresholds silently discarded.
#
# This test proves a spec using min/max produces IDENTICAL violations to the same spec written
# with min_value/max_value.
scores <- data.frame(
  parent1 = c("A", "C"),
  parent2 = c("B", "D"),
  yield_value = c(66, 90),
  stringsAsFactors = FALSE
)

traits_canonical <- data.frame(
  trait = "yield",
  column = "yield_value",
  direction = "maximize",
  weight = 1,
  min_value = 70,
  max_value = NA_real_,
  stringsAsFactors = FALSE
)

traits_alias <- data.frame(
  trait = "yield",
  column = "yield_value",
  direction = "maximize",
  weight = 1,
  min = 70,
  max = NA_real_,
  stringsAsFactors = FALSE
)

scored_canonical <- ng_add_multitrait_score(scores, traits_canonical, method = "weighted")
scored_alias <- ng_add_multitrait_score(scores, traits_alias, method = "weighted")

violation_col <- "multi_trait_yield_violation"
stopifnot(violation_col %in% names(scored_canonical))
stopifnot(violation_col %in% names(scored_alias))

# Row 1 (66) is below min_value/min (70) -> both specs must register the same violation.
stopifnot(scored_canonical[[violation_col]][[1L]] > 0)
stopifnot(isTRUE(all.equal(scored_alias[[violation_col]], scored_canonical[[violation_col]])))

# Row 2 (90) is above the threshold in both specs -> no violation either way.
stopifnot(scored_canonical[[violation_col]][[2L]] == 0)
stopifnot(scored_alias[[violation_col]][[2L]] == 0)

# minimum/maximum spellings must also work.
traits_alias2 <- data.frame(
  trait = "yield",
  column = "yield_value",
  direction = "maximize",
  weight = 1,
  minimum = 70,
  maximum = NA_real_,
  stringsAsFactors = FALSE
)
scored_alias2 <- ng_add_multitrait_score(scores, traits_alias2, method = "weighted")
stopifnot(isTRUE(all.equal(scored_alias2[[violation_col]], scored_canonical[[violation_col]])))

# Collision: a spec somehow carrying both min_value and min must prefer the canonical min_value
# and warn, rather than silently picking one.
traits_collision <- data.frame(
  trait = "yield",
  column = "yield_value",
  direction = "maximize",
  weight = 1,
  min_value = 70,
  min = 10,        # would NOT trigger a violation for row 1 if wrongly preferred
  max_value = NA_real_,
  stringsAsFactors = FALSE
)
collision_warning <- NULL
scored_collision <- withCallingHandlers(
  ng_add_multitrait_score(scores, traits_collision, method = "weighted"),
  warning = function(w) {
    collision_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(!is.null(collision_warning))
stopifnot(grepl("min_value", collision_warning, fixed = TRUE))
stopifnot(grepl("min", collision_warning, fixed = TRUE))
# canonical min_value (70) wins, so row 1 still violates.
stopifnot(scored_collision[[violation_col]][[1L]] > 0)

# The alias must also survive ng_breeder_selection_objective(), the entry point R/39 uses.
objective_alias <- ng_breeder_selection_objective(trait = traits_alias, method = "weighted")
stopifnot(objective_alias$traits$min_value == 70)
stopifnot(is.na(objective_alias$traits$max_value))

cat("threshold alias columns tests passed\n")
