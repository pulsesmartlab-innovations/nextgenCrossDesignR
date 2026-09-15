# The frontend (nextgenCrossWorkbench) drives this package out-of-process and
# keys off these strings and argument names, so renaming any of them silently
# breaks a separately released consumer that this repository cannot see. These
# check the *installed* package, which is what the frontend actually loads.
#
# The wider "every documented parameter is a real argument" check needs
# docs/frontend/contracts/config_schema.json, which is not installed with the
# package; it lives in the repository harness (tests/contract_schema_drift.R).

test_that("the capability registry keeps its schema version", {
  reg <- ng_backend_capability_registry()
  expect_identical(reg$schema_version, "ng_backend_capabilities.v2")
})

test_that("the controls section enumerates dropdowns with valid defaults", {
  reg <- ng_backend_capability_registry()
  expect_true(is.list(reg$controls) && length(reg$controls) >= 15L)
  ids <- vapply(reg$controls, function(x) x$id, character(1))
  expect_true(all(c("trait_value_metric", "multi_trait_method", "optimizer",
                    "allocation_method", "progeny") %in% ids))
  # Type-aware. This demanded `choices` of EVERY control, which is exactly what pinned
  # the registry to enums: no numeric parameter could be declared, so the workbench had
  # to hardcode each one -- and ended up exposing the inert min_effect_reliability while
  # min_cv_predictive_r2, the threshold that governs mean selection, was undeclarable.
  for (ctl in reg$controls) {
    expect_true(all(c("id", "label", "group", "type", "default") %in% names(ctl)))
    if (identical(ctl$type, "enum")) {
      expect_true("choices" %in% names(ctl))
      vals <- vapply(ctl$choices, function(x) x$value, character(1))
      expect_true(length(vals) >= 2L && ctl$default %in% vals)
    } else if (identical(ctl$type, "number")) {
      expect_true(all(c("min", "max", "step") %in% names(ctl)))
      expect_true(is.numeric(ctl$default) && length(ctl$default) == 1L &&
                    is.finite(ctl$default))
      expect_true(is.finite(ctl$min) && is.finite(ctl$max) && ctl$min < ctl$max)
      expect_true(ctl$default >= ctl$min && ctl$default <= ctl$max)
    } else {
      fail(paste("unknown control type in the registry:", ctl$type))
    }
  }
})

test_that("the capability registry keeps the sections the frontend reads", {
  reg <- ng_backend_capability_registry()
  expect_true(all(c("data_qc", "breeding_systems", "method_families",
                    "workflows", "external_integrations", "navigation")
                  %in% names(reg)))
})

test_that("the run entrypoint keeps the arguments the frontend sends", {
  # Note the names: the contract uses map_file / direction_file, not
  # marker_map_file / trait_direction_file.
  fm <- names(formals(ng_run_cross_prediction))
  expect_true(all(c("genotype_file", "phenotype_file", "map_file",
                    "direction_file", "allocation_method", "optimizer",
                    "n_crosses", "max_crosses_per_parent", "prediction_mode") %in% fm))
})

test_that("the polyploid entrypoint the frontend routes to still exists", {
  # The workbench sends workflow = "polyploid_design" and calls this directly.
  expect_true(is.function(ng_polyploid_design_crosses))
})
