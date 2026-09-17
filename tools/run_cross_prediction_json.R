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
# Resolution order matters: this searched getwd()/nextgen_cross_design BEFORE
# getwd(), so an untracked stale copy beside the real sources won and the headless
# frontend path ran an eleven-version-old package while every in-process test ran
# the current one. Shared with tests/helper_load.R so the two cannot diverge again.
local({
  here <- dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1L]),
    winslash = "/", mustWork = FALSE))
  finder <- file.path(here, "ng_find_package_root.R")
  source(if (file.exists(finder)) finder
         else file.path(getwd(), "tools", "ng_find_package_root.R"))
})
root <- ng_find_package_root(getwd())
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
# workflow="batch" runs many INDEPENDENT single-trait analyses at once, each writing its
# own result.json + workbook under `output_root`, and returns a manifest instead of a single
# result. Its control keys are split off here for the same reason workflow/stage are: the
# validation below rejects anything that is not a runner argument.
workflow <- if (!is.null(cfg$workflow)) as.character(cfg$workflow) else "full"
stage    <- if (!is.null(cfg$stage)) as.character(cfg$stage) else NULL
batch_jobs <- cfg$batch_jobs
batch_output_root <- if (!is.null(cfg$batch_output_root)) as.character(cfg$batch_output_root) else NULL
batch_workers <- if (!is.null(cfg$batch_workers)) as.integer(cfg$batch_workers) else NULL
batch_job_dir <- if (!is.null(cfg$job_dir)) as.character(cfg$job_dir) else NULL
batch_shared_dir <- if (!is.null(cfg$shared_dir)) as.character(cfg$shared_dir) else NULL
batch_resume <- isTRUE(cfg$resume)
cfg$workflow <- NULL; cfg$stage <- NULL
cfg$batch_jobs <- NULL; cfg$batch_output_root <- NULL; cfg$batch_workers <- NULL
cfg$job_dir <- NULL; cfg$shared_dir <- NULL; cfg$resume <- NULL

# Write the job record FIRST. The frontend spawns this process detached and then looks for
# job.json; if the backend cannot load or the config is bad, a job that never wrote a record
# is indistinguishable from one that was never launched, and "nothing happened" is the worst
# possible feedback for a run expected to take hours.
if (!is.null(batch_job_dir)) {
  ng_job_create(batch_job_dir, config = cfg, label = basename(batch_job_dir))
}

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

# The envelope, the non-finite sanitation and the writer options all live in the package
# now (R/54_batch_runner.R). A batch writes one result.json per job, and a second copy of
# these rules here would drift from the contract the frontend parses.
write_result <- function(env) ng_write_result_json(env, result_path)

# Run the pipeline while (1) capturing every warning it emits (e.g. the
# residual-heterozygous-RIL advisory, the assume_inbred deprecation) so the
# frontend can surface them, and (2) turning a blocker -- a hard ng_stop such as
# a DH/inbred parent carrying heterozygosity, or any input error -- into a
# structured result.json (status = "error") instead of an uncaught crash, so the
# frontend can DISPLAY the message rather than only see a failed process.
run_warnings <- list()
staged <- identical(workflow, "stage")
batched <- identical(workflow, "batch")
res <- tryCatch(
  withCallingHandlers(
    if (batched) {
      # Default the output root beside result.json, matching how run_dir is derived for the
      # staged workflow -- a caller that says nothing still gets its files somewhere sensible.
      root_dir <- if (!is.null(batch_output_root)) batch_output_root else
        file.path(dirname(result_path), "batch")
      ng_run_cross_prediction_batch(cfg, jobs = batch_jobs, output_root = root_dir,
                                    batch_workers = batch_workers,
                                    job_dir = batch_job_dir,
                                    shared_dir = batch_shared_dir,
                                    resume = batch_resume,
                                    generated_at = generated_at, package_version = version)
    } else if (staged) {
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
  if (!is.null(batch_job_dir) && file.exists(file.path(batch_job_dir, "job.json"))) {
    ng_job_mark(batch_job_dir, "failed", list(error_message = res$message))
  }
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

# --- assemble the schema-versioned result envelope --------------------------------------------
# ng_run_result_envelope() is the single definition, shared with the batch runner.
full_envelope <- function(r) {
  ng_run_result_envelope(r, warnings = run_warnings, generated_at = generated_at,
                         package_version = version)
}
if (batched) {
  # A batch has no single result -- each job wrote its own result.json where the caller can
  # read it. What goes here is the manifest: which traits ran, how they fared, and where
  # their files are. `ok` reflects the BATCH completing, not every job succeeding; a job that
  # failed is recorded with its message rather than taking the others down, so a reader must
  # look at the per-job statuses. That distinction is stated here because a caller keying
  # only on `ok` would otherwise believe seventeen plans exist when sixteen do.
  mf <- res$manifest
  write_result(list(
    schema = "ng_batch_result.v1",
    status = "ok", ok = TRUE,
    generated_at = generated_at, package_version = version,
    output_root = mf$output_root,
    n_jobs = mf$n_jobs,
    n_ok = sum(vapply(mf$jobs, function(j) identical(j$status, "ok"), logical(1))),
    workers = mf$workers, worker_basis = mf$worker_basis,
    manifest_path = res$manifest_path,
    jobs = mf$jobs,
    warnings = run_warnings))
  cat("Wrote batch manifest: ",
      normalizePath(result_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
  for (j in mf$jobs) cat("  ", j$id, ": ", j$status, "\n", sep = "")
  quit(status = 0L, save = "no")
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
