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

# Read a job's state, deriving what cannot be written.
#
# Two derivations live here and nowhere else:
#
#   crashed    -- the record says "running" but the heartbeat has gone quiet. Nothing writes
#                 this, because the writer would have to be the process that died.
#   incomplete -- a trait claiming "running" under a job that is not running. A worker that
#                 dies cannot write its own epitaph, and if the parent died too, nobody
#                 reconciles it.
#
# Keeping both in one function is the point. A frontend that recomputed them would eventually
# disagree with the backend about whether a job is alive, and the disagreement would surface
# as a breeder acting on a plan that was never finished.
ng_job_status <- function(job_dir, stale_after_sec = 120) {
  path <- file.path(job_dir, "job.json")
  if (!file.exists(path)) {
    ng_stop("no job record at ", job_dir,
            " -- a directory without job.json is not a job; if a process was launched for ",
            "it and never wrote one, it failed during startup and its log will say why.")
  }
  rec <- as.list(jsonlite::fromJSON(path, simplifyVector = TRUE))
  state <- rec$state %||% "queued"

  if (identical(state, "running")) {
    hb <- file.path(job_dir, "heartbeat")
    age <- if (file.exists(hb)) {
      as.numeric(difftime(Sys.time(), file.mtime(hb), units = "secs"))
    } else Inf
    if (is.finite(stale_after_sec) && age > stale_after_sec) state <- "crashed"
  }

  dirs <- list.dirs(job_dir, recursive = FALSE, full.names = TRUE)
  files <- file.path(dirs, "status.json")
  files <- files[file.exists(files)]
  live <- identical(state, "running")
  traits <- if (!length(files)) {
    data.frame(trait = character(0), state = character(0), cv_predictive_r2 = numeric(0),
               mean_source = character(0), effect_gate = character(0),
               n_selected = integer(0), error_message = character(0),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, lapply(files, function(f) {
      st <- tryCatch(as.list(jsonlite::fromJSON(f, simplifyVector = TRUE)),
                     error = function(e) list())
      ts <- st$state %||% NA_character_
      # The rule: a trait's state is only meaningful relative to its job's.
      if (identical(ts, "running") && !live) ts <- "incomplete"
      chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x)[[1L]]
      num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[[1L]]
      data.frame(trait = chr(st$trait), state = ts, cv_predictive_r2 = num(st$cv_predictive_r2),
                 mean_source = chr(st$mean_source), effect_gate = chr(st$effect_gate),
                 n_selected = as.integer(num(st$n_selected)),
                 error_message = chr(st$error_message), stringsAsFactors = FALSE)
    }))
  }
  traits <- traits[order(traits$trait), , drop = FALSE]
  rownames(traits) <- NULL

  count <- function(w) sum(traits$state == w, na.rm = TRUE)
  list(schema = "ng_job_status.v1", id = rec$id %||% basename(job_dir),
       label = rec$label %||% basename(job_dir), created_at = rec$created_at,
       state = state, pid = rec$pid,
       n_traits = if (is.null(rec$n_traits)) nrow(traits) else as.integer(rec$n_traits),
       n_done = count("done"), n_error = count("error"),
       n_running = count("running"), n_incomplete = count("incomplete"),
       traits = traits, path = job_dir)
}

# List the jobs in a directory, newest first.
#
# Tolerant by construction: the jobs directory is shared, and a stray folder, a dotfile or a
# half-made directory from an interrupted submit must not take the listing down. The listing
# is the only route a breeder has to their results, so it degrades by omitting what it cannot
# read rather than by failing.
ng_job_list <- function(jobs_dir, stale_after_sec = 120, limit = 100L) {
  empty <- data.frame(id = character(0), label = character(0), created_at = character(0),
                      state = character(0), n_traits = integer(0), n_done = integer(0),
                      n_error = integer(0), path = character(0), stringsAsFactors = FALSE)
  if (!length(jobs_dir) || is.na(jobs_dir) || !dir.exists(jobs_dir)) return(empty)
  dirs <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "job.json"))]
  if (!length(dirs)) return(empty)
  dirs <- dirs[order(file.mtime(dirs), decreasing = TRUE)]
  if (is.finite(limit) && length(dirs) > limit) dirs <- dirs[seq_len(limit)]
  rows <- lapply(dirs, function(d) {
    s <- tryCatch(ng_job_status(d, stale_after_sec = stale_after_sec), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    data.frame(id = s$id, label = s$label,
               created_at = if (is.null(s$created_at)) NA_character_ else s$created_at,
               state = s$state, n_traits = as.integer(s$n_traits),
               n_done = as.integer(s$n_done), n_error = as.integer(s$n_error),
               path = d, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out <- out[order(out$created_at, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
