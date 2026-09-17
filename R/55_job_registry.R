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
