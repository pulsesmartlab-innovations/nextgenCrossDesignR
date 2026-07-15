root_candidates <- c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

registry <- ng_backend_capability_registry(
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)

stopifnot(identical(registry$schema_version, "ng_backend_capabilities.v1"))
stopifnot(identical(registry$generated_at, "2026-05-07T00:00:00Z"))
stopifnot(is.data.frame(registry$data_qc))
stopifnot(is.data.frame(registry$breeding_systems))
stopifnot(is.data.frame(registry$method_families))
stopifnot(is.data.frame(registry$workflows))
stopifnot(is.data.frame(registry$external_integrations))

qc_ids <- registry$data_qc$id
stopifnot(all(c(
  "duplicate_parent_ids",
  "duplicate_marker_ids",
  "duplicate_trait_names",
  "duplicate_candidate_crosses",
  "reciprocal_candidate_crosses",
  "genotype_phenotype_id_mismatch"
) %in% qc_ids))

systems <- registry$breeding_systems
stopifnot(all(c("diploid_dh", "diploid_ril", "autotetraploid_4x", "allopolyploid_subgenome") %in% systems$id))
stopifnot(identical(systems$default[systems$id == "diploid_dh"], TRUE))
stopifnot(identical(systems$backend_target[systems$id == "diploid_ril"], "RIL"))

methods <- registry$method_families
stopifnot(all(c(
  "multitrait_auto",
  "multitrait_weighted",
  "multitrait_economic_index",
  "multitrait_desired_gain",
  "multitrait_threshold",
  "dh_ril_pmv_scoring",
  "ocs_allocation",
  "external_baselines",
  "native_external_components",
  "poly4x_policy"
) %in% methods$id))
stopifnot(identical(methods$default[methods$id == "multitrait_auto"], TRUE))
stopifnot(grepl("covariance", methods$evidence[methods$id == "multitrait_economic_index"], ignore.case = TRUE))

workflow_ids <- registry$workflows$id
stopifnot(all(c(
  "export_backend_capabilities",
  "head_to_head",
  "multitrait_validation",
  "multitrait_validation_grid",
  "multitrait_crop_validation_grid",
  "crop_genome_scenarios",
  "family_calibration",
  "poly4x_benchmark",
  "poly4x_controlled_ocs",
  "render_visual_report",
  "crossing_plan_workbook",
  "export_dashboard_json"
) %in% workflow_ids))

external_ids <- registry$external_integrations$id
stopifnot(all(c("popvar", "simplemating", "alphamate") %in% external_ids))
stopifnot("native_fallback_function" %in% names(registry$external_integrations))
stopifnot(all(nzchar(registry$external_integrations$native_fallback_function)))
stopifnot(any(grepl("ng_popvar_style_scores", registry$external_integrations$native_fallback_function, fixed = TRUE)))
stopifnot(any(grepl("ng_alphamate_style_select", registry$external_integrations$native_fallback_function, fixed = TRUE)))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  message("jsonlite unavailable; skipping registry JSON export check")
  quit(status = 0)
}

tmp <- tempfile("ng_backend_capabilities_")
dir.create(tmp, recursive = TRUE)
out_path <- file.path(tmp, "backend_capabilities.json")
returned <- ng_write_backend_capability_registry_json(
  output_path = out_path,
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)
stopifnot(identical(normalizePath(returned, winslash = "/", mustWork = TRUE),
                    normalizePath(out_path, winslash = "/", mustWork = TRUE)))
payload <- jsonlite::fromJSON(out_path, simplifyVector = FALSE)
stopifnot(identical(payload$schema_version, "ng_backend_capabilities.v1"))
stopifnot(length(payload$data_qc) >= length(qc_ids))

cli_out <- file.path(tmp, "cli_backend_capabilities.json")
env <- c(NG_BACKEND_CAPABILITIES_OUT = cli_out)
old_env <- Sys.getenv(names(env), unset = NA_character_)
restore_env <- function() {
  for (name in names(old_env)) {
    if (is.na(old_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old_env[[name]], name)))
    }
  }
}
on.exit(restore_env(), add = TRUE)
do.call(Sys.setenv, as.list(env))
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
cli_output <- system2(
  rscript,
  normalizePath(file.path(root, "tools", "export_backend_capabilities_json.R"), winslash = "/", mustWork = FALSE),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cli_output, "status")
if (is.null(status)) status <- 0L
if (!identical(as.integer(status), 0L)) print(cli_output)
stopifnot(length(status) == 1L, is.finite(status), status == 0L)
stopifnot(file.exists(cli_out))

cat("backend capability registry tests passed\n")
