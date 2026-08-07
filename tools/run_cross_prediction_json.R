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
# Staged pipeline: the frontend can drive ONE stage at a time by adding
# workflow="stage" + stage=<qc|predict|index|allocate|rank>. Split those two
# control keys off before validating the rest against ng_run_cross_prediction().
workflow <- if (!is.null(cfg$workflow)) as.character(cfg$workflow) else "full"
stage    <- if (!is.null(cfg$stage)) as.character(cfg$stage) else NULL
cfg$workflow <- NULL; cfg$stage <- NULL
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
if (identical(workflow, "stage")) {
  # ng_run_stage()/ng_cp__build_ctx() expect a COMPLETE config (every formal
  # resolved); the monolith path gets that for free from do.call. Fill any keys
  # the frontend did not send from ng_run_cross_prediction()'s own defaults.
  filler <- function() NULL
  formals(filler) <- formals(ng_run_cross_prediction)
  body(filler) <- quote(mget(names(formals())))
  cfg <- tryCatch(do.call(filler, cfg), error = function(e) cfg)
}

version <- tryCatch(read.dcf(file.path(root, "DESCRIPTION"), fields = "Version")[1L, 1L],
                    error = function(e) NA_character_)
generated_at <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")

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
write_result <- function(env) {
  env <- sanitize(env)
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(env, result_path, auto_unbox = TRUE, na = "null", null = "null",
                       dataframe = "rows", pretty = TRUE, digits = 8)
}

# Run the pipeline while (1) capturing every warning it emits (e.g. the
# residual-heterozygous-RIL advisory, the assume_inbred deprecation) so the
# frontend can surface them, and (2) turning a blocker -- a hard ng_stop such as
# a DH/inbred parent carrying heterozygosity, or any input error -- into a
# structured result.json (status = "error") instead of an uncaught crash, so the
# frontend can DISPLAY the message rather than only see a failed process.
run_warnings <- list()
staged <- identical(workflow, "stage")
res <- tryCatch(
  withCallingHandlers(
    if (staged) {
      if (is.null(stage)) ng_stop("workflow = 'stage' requires a 'stage' key")
      ng_run_stage(stage, dirname(result_path), cfg)   # one gated stage; persists to run_dir/artifacts
    } else {
      do.call(ng_run_cross_prediction, cfg)
    },
    warning = function(w) {
      run_warnings[[length(run_warnings) + 1L]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) structure(list(message = conditionMessage(e)), class = "ng_run_json_error")
)

if (inherits(res, "ng_run_json_error")) {
  write_result(list(
    schema          = "ng_run_result.v1",
    status          = "error",
    ok              = FALSE,               # frontend keys success on `ok`
    stage           = stage,
    error           = TRUE,
    error_message   = res$message,         # flat field the frontend renders
    generated_at    = generated_at,
    package_version = version,
    warnings        = run_warnings
  ))
  cat("Run blocked: ", res$message, "\n", sep = "", file = stderr())
  cat("Wrote error result: ", normalizePath(result_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  quit(status = 1L, save = "no")   # nonzero exit, but result.json carries the message
}

# --- assemble the schema-versioned result envelope (JSON-friendly subset) ----------------------
# The full monolith/rank result contract. `warnings` (residual-het RIL advisory,
# deprecations, ...) is carried for the frontend to display alongside the results.
full_envelope <- function(r) {
  audit <- r$input_match_audit
  audit$marker_order <- NULL   # drop the (potentially huge) marker-name vector; marker_count kept
  list(
    schema          = "ng_run_result.v1",
    status          = "ok",
    ok              = TRUE,               # frontend keys success on `ok`
    generated_at    = generated_at,
    package_version = version,
    prediction_mode = r$prediction_mode,
    settings          = r$settings,
    input_match_audit = audit,
    qc = list(status = r$qc$status, counts = r$qc$counts,
              tables = r$qc$tables, issues = r$qc$issues),
    effect_summary  = r$effect_summary,
    trait_direction = r$trait_direction,
    objective       = if (!is.null(r$objective)) r$objective$diagnostics else NULL,
    plan_summary    = r$plan_summary,
    constraint_diagnostics = r$constraint_diagnostics,
    priority_risk_diagnostics = r$priority_risk_diagnostics,
    trait_check_diagnostics = r$trait_check_diagnostics,
    candidate_crosses = r$candidate_crosses,
    selected_crosses  = r$selected_crosses,
    ld_pruning_report = r$ld_pruning_report,
    warnings          = run_warnings,
    output_files      = r$output_files
  )
}

if (staged) {
  # One gated stage. The rank stage carries the full assembled result (attr);
  # earlier stages return {status, stage, files}. Either way, surface `status`,
  # `stage`, and any `warnings` so the frontend's staged view can display them.
  full <- attr(res, "result")
  if (!is.null(full)) {
    envelope <- full_envelope(full)
    envelope$stage <- res$stage
    envelope$files <- res$files
  } else {
    envelope <- list(
      schema = "ng_run_result.v1", status = res$status, ok = TRUE, stage = res$stage,
      files = res$files, generated_at = generated_at, package_version = version,
      warnings = run_warnings)
  }
  write_result(envelope)
  cat("Wrote stage result (", res$stage, ", ", res$status, "): ",
      normalizePath(result_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  quit(status = 0L, save = "no")
}

envelope <- full_envelope(res)
write_result(envelope)
cat("Wrote run result: ", normalizePath(result_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (length(res$output_files)) {
  for (nm in names(res$output_files)) cat("  ", nm, ": ", res$output_files[[nm]], "\n", sep = "")
}
