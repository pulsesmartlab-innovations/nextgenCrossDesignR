# Headless JSON-in / JSON-out run wrapper for the frontend.
#
# A single entry point any frontend can spawn: it reads a JSON config whose keys mirror
# ng_run_cross_prediction() arguments (input file paths + parameters), runs the full pipeline,
# and writes a schema-versioned result.json (the frontend "result contract") plus the optional
# xlsx workbook / PNG figures. No R knowledge required on the frontend side.
#
# Usage:
#   NG_RUN_CONFIG=path/to/config.json NG_RUN_RESULT_OUT=path/to/result.json \
#     Rscript tools/run_cross_prediction_json.R
#   # or: Rscript tools/run_cross_prediction_json.R config.json result.json
#
# Config: a JSON object of ng_run_cross_prediction() arguments. Use *_file keys for inputs
# (phenotype_file/genotype_file/map_file/direction_file[/training_*_file]) so the frontend never
# ships in-memory matrices. Unknown keys are a hard error (never silently dropped). See
# docs/frontend/contracts/config_schema.json for the full parameter menu and
# docs/frontend/contracts/example_config.json for a runnable example.

# --- locate the package root and load it (pure-R; mirrors tools/export_*.R) -------------------
root_candidates <- unique(normalizePath(c(
  file.path(getwd(), "nextgen_cross_design"), getwd(),
  file.path("..", "nextgen_cross_design"), file.path("..")
), winslash = "/", mustWork = FALSE))
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
if (!length(root_hits)) stop("Could not locate nextgen_cross_design root", call. = FALSE)
root <- root_hits[[1L]]
source(file.path(root, "tools", "ng_project_libpath.R"))
ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = !identical(Sys.getenv("NGCD_SKIP_CPP"), "1"), verbose = FALSE)

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required", call. = FALSE)

# --- resolve config + output paths (env vars, then positional args) ---------------------------
args <- commandArgs(trailingOnly = TRUE)
config_path <- Sys.getenv("NG_RUN_CONFIG", unset = if (length(args) >= 1L) args[[1L]] else "")
result_path <- Sys.getenv("NG_RUN_RESULT_OUT",
                          unset = if (length(args) >= 2L) args[[2L]] else file.path(root, "results", "run_result.json"))
if (!nzchar(config_path)) stop("Set NG_RUN_CONFIG (or pass the config path as arg 1)", call. = FALSE)
if (!file.exists(config_path)) stop("Config file not found: ", config_path, call. = FALSE)

# --- parse + validate the config against the runner's formals ---------------------------------
cfg <- jsonlite::fromJSON(config_path, simplifyVector = TRUE, simplifyDataFrame = TRUE)
if (!is.list(cfg)) stop("Config must be a JSON object of ng_run_cross_prediction arguments", call. = FALSE)
valid_args <- names(formals(ng_run_cross_prediction))
unknown <- setdiff(names(cfg), valid_args)
if (length(unknown)) {
  stop("Unknown config key(s) (not ng_run_cross_prediction arguments): ",
       paste(unknown, collapse = ", "),
       ". See docs/frontend/contracts/config_schema.json.", call. = FALSE)
}
# A JSON object of scalars (e.g. {"yield":0.4,...}) parses to a named list; the runner wants a
# named numeric vector for trait_weights.
if (!is.null(cfg$trait_weights) && is.list(cfg$trait_weights)) {
  cfg$trait_weights <- unlist(cfg$trait_weights)
}

res <- do.call(ng_run_cross_prediction, cfg)

# --- assemble the schema-versioned result envelope (JSON-friendly subset) ----------------------
version <- tryCatch(read.dcf(file.path(root, "DESCRIPTION"), fields = "Version")[1L, 1L],
                    error = function(e) NA_character_)
generated_at <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")

audit <- res$input_match_audit
audit$marker_order <- NULL   # drop the (potentially huge) marker-name vector; marker_count kept

envelope <- list(
  schema          = "ng_run_result.v1",
  generated_at    = generated_at,
  package_version = version,
  prediction_mode = res$prediction_mode,
  settings          = res$settings,
  input_match_audit = audit,
  qc = list(status = res$qc$status, counts = res$qc$counts,
            tables = res$qc$tables, issues = res$qc$issues),
  effect_summary  = res$effect_summary,
  trait_direction = res$trait_direction,
  objective       = if (!is.null(res$objective)) res$objective$diagnostics else NULL,
  plan_summary    = res$plan_summary,
  constraint_diagnostics = res$constraint_diagnostics,
  priority_risk_diagnostics = res$priority_risk_diagnostics,
  trait_check_diagnostics = res$trait_check_diagnostics,
  candidate_crosses = res$candidate_crosses,
  selected_crosses  = res$selected_crosses,
  ld_pruning_report = res$ld_pruning_report,
  output_files      = res$output_files
)

# Non-finite (Inf/NaN) -> null so the JSON is valid and consumable everywhere.
sanitize <- function(x) {
  if (is.data.frame(x)) {
    for (j in seq_along(x)) if (is.numeric(x[[j]])) x[[j]][!is.finite(x[[j]])] <- NA
    return(x)
  }
  if (is.list(x)) return(lapply(x, sanitize))
  if (is.numeric(x)) x[!is.finite(x)] <- NA
  x
}
envelope <- sanitize(envelope)

dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(envelope, result_path, auto_unbox = TRUE, na = "null", null = "null",
                     dataframe = "rows", pretty = TRUE, digits = 8)
cat("Wrote run result: ", normalizePath(result_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (length(res$output_files)) {
  for (nm in names(res$output_files)) cat("  ", nm, ": ", res$output_files[[nm]], "\n", sep = "")
}
