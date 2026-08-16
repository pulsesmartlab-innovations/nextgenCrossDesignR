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

stopifnot(identical(registry$schema_version, "ng_backend_capabilities.v2"))
stopifnot(identical(registry$generated_at, "2026-05-07T00:00:00Z"))

# controls section: enumerable UI controls the frontend renders dropdowns from
stopifnot(is.list(registry$controls), length(registry$controls) >= 15L)
for (ctl in registry$controls) {
  stopifnot(all(c("id", "label", "group", "type", "default", "choices") %in% names(ctl)))
  stopifnot(is.list(ctl$choices), length(ctl$choices) >= 2L)
  for (cho in ctl$choices) stopifnot(all(c("value", "label") %in% names(cho)))
}
ctl_ids <- vapply(registry$controls, function(x) x$id, character(1))
stopifnot(all(c("trait_value_metric", "multi_trait_method", "optimizer",
                "allocation_method", "progeny") %in% ctl_ids))
# the choice list must be exhaustive: multi_trait_method includes every backend method
mtm <- registry$controls[[which(ctl_ids == "multi_trait_method")]]
mtm_vals <- vapply(mtm$choices, function(x) x$value, character(1))
stopifnot(all(c("auto", "weighted", "economic_index", "desired_gain", "threshold") %in% mtm_vals))
# defaults must be valid choices
for (ctl in registry$controls) {
  vals <- vapply(ctl$choices, function(x) x$value, character(1))
  stopifnot(ctl$default %in% vals)
}
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
stopifnot(identical(payload$schema_version, "ng_backend_capabilities.v2"))
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

# --- cross-priority risk & portfolio family (0.19.0) ---------------------------------------
# The whole risk/portfolio layer was previously undiscoverable: the frontend reads capabilities
# from this registry and had no entry for it.
crp <- methods[methods$id == "cross_priority_risk_portfolio", , drop = FALSE]
stopifnot(nrow(crp) == 1L)
stopifnot(identical(crp$category, "scoring"), isTRUE(crp$default))
stopifnot(grepl("ng_annotate_cross_priority_multitrait", crp$backend_function, fixed = TRUE))
# The evidence has to carry the load-bearing frontend contract: the basis vocabulary, the badge
# obligation for rank indices, and the within-run-only caveat.
for (needle in c("portfolio_basis", "linearized_rank_index", "linear_index", "single_trait",
                 "badge", "within-run", "cross_upside", "risk_driver_trait")) {
  stopifnot(grepl(needle, crp$evidence, fixed = TRUE))
}
# Appending the row must not have shifted any other family's parallel fields.
stopifnot(identical(methods$label[methods$id == "multitrait_auto"], "Auto multi-trait default"))
stopifnot(identical(methods$category[methods$id == "ocs_allocation"], "allocation"))
stopifnot(!anyNA(methods$id), !anyNA(methods$label), !anyNA(methods$category),
          !anyNA(methods$default), !anyNA(methods$backend_function), !anyNA(methods$evidence))
cat("cross-priority risk/portfolio registry family test passed\n")
