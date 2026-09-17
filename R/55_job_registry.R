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
#
# Two-pass filtering for correctness and efficiency:
#   - Pass 1: read ONLY job.json (tiny) from each candidate, sort by created_at (not mtime,
#     which reflects the last state transition, not creation), and slice to limit. This avoids
#     expensive ng_job_status() calls on jobs that will be discarded.
#   - Pass 2: call ng_job_status() only on survivors to build the full row.
ng_job_list <- function(jobs_dir, stale_after_sec = 120, limit = 100L) {
  empty <- data.frame(id = character(0), label = character(0), created_at = character(0),
                      state = character(0), n_traits = integer(0), n_done = integer(0),
                      n_error = integer(0), path = character(0), stringsAsFactors = FALSE)
  if (!length(jobs_dir) || is.na(jobs_dir) || !dir.exists(jobs_dir)) return(empty)
  dirs <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "job.json"))]
  if (!length(dirs)) return(empty)

  # First pass: cheap. Read only job.json to get created_at.
  metadata <- lapply(dirs, function(d) {
    job_rec <- tryCatch(as.list(jsonlite::fromJSON(file.path(d, "job.json"),
                                                    simplifyVector = TRUE)),
                        error = function(e) NULL)
    if (is.null(job_rec)) return(NULL)
    data.frame(id = job_rec$id %||% basename(d),
               created_at = job_rec$created_at %||% NA_character_,
               path = d, stringsAsFactors = FALSE)
  })
  metadata <- do.call(rbind, metadata[!vapply(metadata, is.null, logical(1))])
  if (!nrow(metadata)) return(empty)

  # Sort by created_at (when the job was submitted), not mtime (when it last changed state).
  metadata <- metadata[order(metadata$created_at, decreasing = TRUE), , drop = FALSE]
  if (is.finite(limit) && nrow(metadata) > limit) {
    metadata <- metadata[seq_len(limit), , drop = FALSE]
  }

  # Second pass: expensive. Call ng_job_status() only on survivors.
  rows <- lapply(metadata$path, function(d) {
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

# The content key of the shared artefact for a configuration.
#
# Derived from ng_cp__batch_shared_keys -- the list 0.36.0 already refuses to let a job
# override, precisely because changing one of those settings would invalidate the shared
# work. That list therefore IS the definition of what the artefact depends on, and the
# existing test asserting it names only real runner arguments protects this cache too.
#
# File inputs are hashed by CONTENT, via tools::md5sum (base R, portable, no dependency).
# Keying on a path, or a path and a timestamp, would let edited data silently reuse an
# artefact built from the old data -- wrong numbers with nothing to signal them.
#
# In-memory inputs are hashed by SERIALISING them, never by summarising them. str() or
# dim() would be cheaper and would be a latent disaster: two different genotype matrices of
# the same shape would hash identically and silently reuse each other's artefact. A cache
# key must be a function of the content or it is not a cache key.
#
# Serialising a large matrix costs a temp file, paid once per batch against recomputing
# quality control. Callers submitting through the frontend materialise CSVs first anyway, so
# the path that matters takes the cheap md5sum-of-file branch.
ng_shared_artifact__digest <- function(x) {
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, compress = FALSE)   # uncompressed: deterministic bytes, and faster
  unname(tools::md5sum(f))
}

ng_shared_artifact_key <- function(config) {
  parts <- lapply(ng_cp__batch_shared_keys, function(k) {
    v <- config[[k]]
    if (is.null(v)) return(paste0(k, "=<null>"))
    if (is.character(v) && length(v) == 1L && !is.na(v) && file.exists(v)) {
      return(paste0(k, "=file:", unname(tools::md5sum(v))))
    }
    paste0(k, "=obj:", ng_shared_artifact__digest(v))
  })
  ng_shared_artifact__digest(paste(unlist(parts), collapse = "\n"))
}

ng_shared_artifact_dir <- function(shared_dir, key) {
  d <- file.path(shared_dir, key)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Record that a job depends on this artefact, so retention can refuse to delete it.
ng_shared_artifact_reference <- function(shared_dir, key, job_id) {
  d <- ng_shared_artifact_dir(shared_dir, key)
  path <- file.path(d, "meta.json")
  meta <- if (file.exists(path)) {
    tryCatch(as.list(jsonlite::fromJSON(path, simplifyVector = TRUE)), error = function(e) list())
  } else list()
  meta$schema <- "ng_shared_artifact.v1"
  meta$key <- key
  if (is.null(meta$created_at)) meta$created_at <- ng_job__now()
  meta$referenced_by <- unique(c(as.character(meta$referenced_by %||% character(0)),
                                 as.character(job_id)))
  ng_write_json_atomic(meta, path)
  invisible(path)
}

# Remove old jobs, and any shared artefact nothing references any more.
#
# This function deletes with unlink(recursive = TRUE, force = TRUE), so its refusals matter
# more than its removals:
#
#   * a RUNNING, QUEUED or CRASHED job is never removed, however old. Deleting a running one
#     would take hours of finished traits with it, and the process writing into it would carry
#     on writing into nothing.
#
#     "crashed" is protected too, and deliberately so. ng_job_status() derives "crashed" from
#     a stale heartbeat, but a stale heartbeat is exactly what a perfectly healthy SERIAL job
#     looks like for the entire duration of one trait: the serial dispatch path beats the
#     heart once per job launch and then blocks inside ng_cp__batch_run_one() for as long as
#     that trait takes -- hours, in this package -- because R is single-threaded and nothing
#     else can touch the heartbeat file while it blocks. A single-trait job therefore has no
#     "between traits" to beat in, goes stale after stale_after_sec, and would otherwise look
#     identical to a job whose process actually died. Raising stale_after_sec cannot fix this:
#     any fixed timeout is a guess that some real single-trait run will still exceed. So
#     "crashed" gets the same protection as "running", and this function accepts the two
#     states can be indistinguishable from disk alone. The asymmetry is intentional: treating
#     a genuinely dead job as alive costs disk space; treating a genuinely alive job as dead
#     destroys a breeder's work. Erring toward keeping a dead job is the only safe direction --
#     and a crashed job is also precisely the one holding partial results a breeder would want
#     to resume, so deleting it would be the worst available choice even when it really is dead.
#   * a shared artefact referenced by ANY surviving job is never removed. The artefact
#     outlives the batch that built it -- that is the point of sharing it -- so deletion is
#     driven by references, not by age.
#
# This is also why jobs_dir and shared_dir are siblings of the frontend's runs_dir rather
# than living inside it: ngcd_prune_runs() there enumerates list.dirs(runs_dir) and unlinks
# everything past keep_runs, and would happily take the whole job store with it.
ng_job_prune <- function(jobs_dir, shared_dir = NULL, keep = 20L, stale_after_sec = 120) {
  removed <- character(0)
  keep <- suppressWarnings(as.integer(keep))
  if (is.na(keep) || keep < 0L || !dir.exists(jobs_dir)) return(invisible(removed))

  lst <- ng_job_list(jobs_dir, stale_after_sec = stale_after_sec, limit = Inf)
  if (nrow(lst)) {
    protected <- lst$state %in% c("running", "queued", "crashed")
    candidates <- lst[!protected, , drop = FALSE]
    if (nrow(candidates) > keep) {
      doomed <- candidates$path[seq.int(keep + 1L, nrow(candidates))]
      unlink(doomed, recursive = TRUE, force = TRUE)
      removed <- c(removed, doomed)
    }
  }

  if (!is.null(shared_dir) && dir.exists(shared_dir)) {
    survivors <- ng_job_list(jobs_dir, stale_after_sec = stale_after_sec, limit = Inf)
    refs <- unique(unlist(lapply(survivors$path, function(d) {
      f <- file.path(d, "shared_ref")
      if (!file.exists(f)) return(character(0))
      # A shared_ref that cannot be read is treated as referencing nothing, not as referencing
      # everything. That is the cautious reading in ONE sense (it never keeps this function
      # from making progress on a batch it can otherwise clean up) but the risky one in
      # another: the alternative -- treating unreadable as "keep" -- would leak an artefact
      # forever the moment its reference file got corrupted, with no way to ever reclaim it.
      # Failing closed on the read, not on the retention decision, keeps that risk bounded.
      tryCatch(trimws(readLines(f, warn = FALSE)), error = function(e) character(0))
    })))
    for (d in list.dirs(shared_dir, recursive = FALSE, full.names = TRUE)) {
      if (!(basename(d) %in% refs)) {
        unlink(d, recursive = TRUE, force = TRUE)
        removed <- c(removed, d)
      }
    }
  }
  invisible(removed)
}

# ---------------------------------------------------------------------------
# What the shared artefact OWNS -- and, just as importantly, what it does not.
# ---------------------------------------------------------------------------
#
# The artefact used to be the whole `ctx`, saved with saveRDS(ctx) and restored with
# readRDS(). But ctx IS the config: ng_cp__build_ctx() starts `ctx <- config` and adds
# derived fields to it. The cache key deliberately covers only ng_cp__batch_shared_keys --
# the settings quality control and the predict prologue actually spend -- so restoring the
# whole ctx restored every NON-shared setting too. A second batch asking for nine crosses on
# the mid-parent mean silently got the first batch's three crosses on usefulness: the spec's
# own headline scenario producing wrong numbers with nothing to signal it, which is the same
# defect class as the 0.36.0 RNG-kind bug.
#
# So the artefact carries ONLY what it computed, and the caller's settings always come from
# the batch that is running now. The list below is exactly the set of fields the two
# artefact-producing functions write:
#
#   ng_cp__stage_qc()        -- its terminal ng_ctx_put() (R/39): phenotype_id_col_used,
#                               genotype_id_col_used, trait_spec, direction_columns,
#                               direction_canonical, marker_map_std, qc, geno, pheno, ids
#   ng_cp__predict_prologue() -- stored whole under ctx$predict_prologue
#
# Nothing else. If a future stage adds a derived field that belongs to the shared work, it
# must be added here deliberately; a field absent from this list is simply recomputed, which
# is the safe direction to fail.
ng_shared_artifact_schema <- "ng_shared_artifact_ctx.v1"

ng_shared_artifact_fields <- c(
  "qc", "geno", "pheno", "ids", "trait_spec", "marker_map_std",
  "direction_columns", "direction_canonical",
  "phenotype_id_col_used", "genotype_id_col_used",
  "predict_prologue")

ng_shared_artifact_write <- function(ctx, path) {
  fields <- ctx[intersect(ng_shared_artifact_fields, names(ctx))]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  # .part-then-rename for the same reason every JSON payload here uses it: a concurrent
  # batch may be reading this exact path while another writes it, and a truncated RDS is a
  # much worse failure than a recomputed one.
  tmp <- paste0(path, ".part")
  saveRDS(list(schema = ng_shared_artifact_schema, fields = fields), tmp)
  file.rename(tmp, path)
  path
}

# Returns the artefact-owned fields, or NULL when there is nothing trustworthy to return.
#
# An unrecognised schema is treated as a MISS rather than as an error or as data: a future
# release that changes what the artefact owns must not silently read an old one as though it
# still meant the same thing, and recomputing is always available as an answer.
ng_shared_artifact_read <- function(path) {
  if (!length(path) || !file.exists(path)) return(NULL)
  art <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.list(art) || !identical(art$schema, ng_shared_artifact_schema)) return(NULL)
  if (!is.list(art$fields) || !length(art$fields)) return(NULL)
  art$fields
}

# The one place a batch context is assembled, used by the parent and by every worker, so the
# two cannot drift: THIS batch's config decides everything, and the artefact supplies only
# the shared work laid over it.
ng_cp__batch_ctx <- function(config, fields) {
  ctx <- ng_cp__build_ctx(config)
  if (length(fields)) ctx <- do.call(ng_ctx_put, c(list(ctx), fields))
  ctx
}
