# Retention deletes with unlink(recursive = TRUE, force = TRUE). Everything it will not
# delete has to be asserted, because the cost of being wrong is a breeder's work.
#
# Two absolute constraints:
#   1. never delete a RUNNING job -- it would take hours of finished traits with it;
#   2. never delete a shared artefact any surviving job still references -- the artefact
#      outlives the batch that built it, which is the entire point of sharing it.
#
# The frontend's ngcd_prune_runs() is the cautionary example: it enumerates
# list.dirs(cfg$runs_dir) and unlinks everything past keep_runs, which is exactly why the
# jobs and shared stores are siblings of runs_dir rather than living inside it.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_pr_", as.integer(runif(1, 1, 1e6))))
jobs_dir <- file.path(root, "jobs"); shared_dir <- file.path(root, "shared")
dir.create(jobs_dir, recursive = TRUE); dir.create(shared_dir, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

mk <- function(id, state, key, write_ref = TRUE) {
  d <- file.path(jobs_dir, id)
  ng_job_create(d, config = NULL, label = id, n_traits = 1L)
  ng_job_mark(d, state, list(pid = Sys.getpid()))
  if (identical(state, "running")) ng_job_heartbeat(d)
  if (write_ref) writeLines(key, file.path(d, "shared_ref"))
  ng_shared_artifact_reference(shared_dir, key, id)
  saveRDS(list(1), file.path(shared_dir, key, "shared.rds"))
  Sys.sleep(1.1)   # distinct created_at, so ordering is by design rather than by tie-break
  d
}

# The serial dispatch path beats the heart once per job launch, then blocks inside
# ng_cp__batch_run_one() for the entire trait -- hours, in this package -- because R is
# single-threaded and nothing else can touch the heartbeat file while it blocks. A
# single-trait job therefore has no "between traits" to beat in: its heartbeat goes stale
# long before the job is actually dead, ng_job_status() reports it "crashed", and
# ng_job_prune() must protect "crashed" exactly as it protects "running" -- deleting a job
# because its heartbeat looks old is indistinguishable, from disk alone, from deleting one
# that is still computing. Made the OLDEST job here so it is also the first candidate a
# naive age-based prune would reach.
stale <- mk("stale", "running", "KEYCRASHED")
Sys.setFileTime(file.path(stale, "heartbeat"), Sys.time() - 200)   # > default stale_after_sec

old1 <- mk("old1", "finished", "KEYOLD")
old2 <- mk("old2", "finished", "KEYSHARED")
live <- mk("live", "running",  "KEYLIVE")     # oldest-but-one, and must survive anyway
new1 <- mk("new1", "finished", "KEYSHARED")   # shares KEYSHARED with old2

# A RUNNING job whose shared_ref file never made it to disk -- an interrupted submit, a
# corrupted file, the reader above deliberately failing closed. Its artefact must still
# survive, because meta.json$referenced_by names it and its directory is still there. Without
# that second source, unlink(recursive = TRUE) would take a live job's shared work.
noref <- mk("noref", "running", "KEYNOREF", write_ref = FALSE)

removed <- ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 1L)

# --- a crashed job -- heartbeat stale, record still "running" -- is never removed -------------
# It is the oldest job of all, and the one a naive prune would reach first.
stopifnot(dir.exists(stale))
stopifnot(!(stale %in% removed))
# --- a running job is never removed, however old ---------------------------------------------
stopifnot(dir.exists(live))
stopifnot(!(live %in% removed))
# --- the newest survives, as retention intends ------------------------------------------------
stopifnot(dir.exists(new1))
# --- the genuinely old, finished ones go -------------------------------------------------------
stopifnot(!dir.exists(old1))

# --- artefacts: referenced ones stay, orphans go ------------------------------------------------
# KEYSHARED is still referenced by new1, so it must survive old2's removal. Deleting it would
# make a surviving job unresumable, and its results unexplainable.
stopifnot(dir.exists(file.path(shared_dir, "KEYSHARED")))
stopifnot(dir.exists(file.path(shared_dir, "KEYLIVE")))
stopifnot(dir.exists(file.path(shared_dir, "KEYCRASHED")))
# KEYOLD's only referent is gone.
stopifnot(!dir.exists(file.path(shared_dir, "KEYOLD")))
# --- meta.json$referenced_by is READ, not merely written ----------------------------------
# noref has no shared_ref file at all, so the shared_ref pass sees nothing. The artefact
# survives only because meta.json names a job directory that still exists.
stopifnot(!file.exists(file.path(noref, "shared_ref")))
stopifnot(dir.exists(file.path(shared_dir, "KEYNOREF")))
# And the field is an ARRAY at length one, not a bare string, so one consumer shape serves
# both -- asserted against what landed on disk.
meta1 <- jsonlite::fromJSON(file.path(shared_dir, "KEYNOREF", "meta.json"),
                            simplifyVector = FALSE)
stopifnot(is.list(meta1$referenced_by), is.null(names(meta1$referenced_by)),
          length(meta1$referenced_by) == 1L)

# --- nothing to do is not an error ----------------------------------------------------------------
stopifnot(identical(length(ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 100L)), 0L))

# --- and an artefact whose referents are ALL gone is still collectable -----------------------
# meta.json only ever grows, so treating its list as proof of life would make every artefact
# immortal. Directory existence is the liveness signal: remove the job, prune again, and the
# artefact goes.
unlink(noref, recursive = TRUE, force = TRUE)
invisible(ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 0L))
stopifnot(!dir.exists(file.path(shared_dir, "KEYNOREF")))

cat("PASS: job_prune_protects_running_and_referenced\n")
