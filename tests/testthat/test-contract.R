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
  expect_identical(reg$schema_version, "ng_backend_capabilities.v1")
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
                    "n_crosses", "max_uses_per_parent", "prediction_mode") %in% fm))
})

test_that("the polyploid entrypoint the frontend routes to still exists", {
  # The workbench sends workflow = "polyploid_design" and calls this directly.
  expect_true(is.function(ng_design_crosses_poly))
})
