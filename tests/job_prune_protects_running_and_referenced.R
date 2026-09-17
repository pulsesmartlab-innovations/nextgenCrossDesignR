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

mk <- function(id, state, key) {
  d <- file.path(jobs_dir, id)
  ng_job_create(d, config = NULL, label = id, n_traits = 1L)
  ng_job_mark(d, state, list(pid = Sys.getpid()))
  if (identical(state, "running")) ng_job_heartbeat(d)
  writeLines(key, file.path(d, "shared_ref"))
  ng_shared_artifact_reference(shared_dir, key, id)
  saveRDS(list(1), file.path(shared_dir, key, "shared.rds"))
  Sys.sleep(0.3)
  d
}
old1 <- mk("old1", "finished", "KEYOLD")
old2 <- mk("old2", "finished", "KEYSHARED")
live <- mk("live", "running",  "KEYLIVE")     # oldest-but-one, and must survive anyway
new1 <- mk("new1", "finished", "KEYSHARED")   # shares KEYSHARED with old2

removed <- ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 1L)

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
# KEYOLD's only referent is gone.
stopifnot(!dir.exists(file.path(shared_dir, "KEYOLD")))

# --- nothing to do is not an error ----------------------------------------------------------------
stopifnot(identical(length(ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 100L)), 0L))

cat("PASS: job_prune_protects_running_and_referenced\n")
