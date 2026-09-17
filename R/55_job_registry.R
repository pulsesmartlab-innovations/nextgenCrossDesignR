# A job is a batch that outlives the process that asked for it.
#
# 0.36.0 made a 17-trait batch computable. This file makes one DURABLE: submitted, left
# running, and readable afterwards by anything that can see the directory. Everything here
# is about a job as an object on disk -- its record, its per-trait progress, whether it is
# still alive -- and nothing here knows how a cross is scored.
#
# The whole design rests on one asymmetry: a worker can report what it did, but it cannot
# report that it died. So state is written optimistically by whoever is alive, and the
# absence of an update is what the reader interprets. See ng_job_status().

ng_job_trait_states <- c("running", "done", "error")

# Record one trait's state beside the results it writes.
#
# Called twice per trait: once on entry so a reader can distinguish a trait being worked on
# from one not yet started, and once on exit with the summary. The summary is SMALL on
# purpose -- an overview of seventeen traits must not mean parsing seventeen result.json
# files of ~125,000 rows each.
ng_job_trait_status_write <- function(job_dir, trait_id, state, fields = list()) {
  if (!state %in% ng_job_trait_states) {
    ng_stop("trait state must be one of ", paste(ng_job_trait_states, collapse = ", "),
            "; got '", state, "'. Job-level states (queued, finished, crashed) belong in ",
            "job.json, not beside a trait.")
  }
  dir <- file.path(job_dir, trait_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "status.json")
  now <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  # Preserve started_at across the second write: without it nothing can say how long a trait
  # took, and elapsed time is the first thing anyone asks of a job that is still running.
  prior <- if (file.exists(path)) {
    tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
  } else NULL
  rec <- c(list(
    schema = "ng_job_trait_status.v1",
    trait = trait_id,
    state = state,
    started_at = if (!is.null(prior$started_at)) prior$started_at else now,
    finished_at = if (identical(state, "running")) NULL else now
  ), fields)
  ng_write_json_atomic(rec, path)
}

# A job may be written as any of these. `crashed` is absent on purpose -- see ng_job_status().
ng_job_states <- c("queued", "running", "finished", "failed")

ng_job__now <- function() format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")

# Create the job record. Written BEFORE the work starts, so a job that dies during startup
# is still a job someone can find and read an error out of, rather than an empty directory.
ng_job_create <- function(job_dir, config = NULL, label = NULL, n_traits = NA_integer_) {
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(config)) ng_write_json_atomic(config, file.path(job_dir, "config.json"))
  rec <- list(
    schema = "ng_job.v1",
    id = basename(job_dir),
    label = if (is.null(label)) basename(job_dir) else label,
    created_at = ng_job__now(),
    state = "queued",
    n_traits = if (is.na(n_traits)) NA_integer_ else as.integer(n_traits),
    pid = NA_integer_
  )
  ng_write_json_atomic(rec, file.path(job_dir, "job.json"))
}

# Move the job to a new state, preserving what was already recorded.
ng_job_mark <- function(job_dir, state, fields = list()) {
  if (!state %in% ng_job_states) {
    ng_stop("job state must be one of ", paste(ng_job_states, collapse = ", "),
            "; got '", state, "'. 'crashed' is derived by ng_job_status() from a stale ",
            "heartbeat and is never written -- the process that would write it is the one ",
            "that died.")
  }
  path <- file.path(job_dir, "job.json")
  rec <- if (file.exists(path)) {
    tryCatch(as.list(jsonlite::fromJSON(path, simplifyVector = TRUE)),
             error = function(e) list())
  } else list()
  rec$schema <- "ng_job.v1"
  if (is.null(rec$id)) rec$id <- basename(job_dir)
  if (is.null(rec$created_at)) rec$created_at <- ng_job__now()
  rec$state <- state
  if (state %in% c("finished", "failed")) rec$finished_at <- ng_job__now()
  for (nm in names(fields)) rec[[nm]] <- fields[[nm]]
  ng_write_json_atomic(rec, path)
}

# Prove the process is still alive. Only the mtime matters; the file's content never does.
ng_job_heartbeat <- function(job_dir) {
  path <- file.path(job_dir, "heartbeat")
  if (!file.exists(path)) {
    dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    file.create(path)
  } else {
    Sys.setFileTime(path, Sys.time())
  }
  invisible(path)
}
