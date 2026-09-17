# Job status is DERIVED, because the two states that matter most cannot be written.
#
# A worker writes "running" on entry. If it dies it cannot write its own epitaph. Normally
# the batch parent reconciles that trait to "error" -- but if the PARENT died too, nobody
# reconciles anything and the trait sits at "running" forever.
#
# So: a trait's state is only meaningful RELATIVE TO ITS JOB'S state. Under a crashed job,
# a trait still claiming "running" is reported as "incomplete". This rule lives here, in one
# tested place, precisely so a frontend never reimplements it and drifts.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

mk <- function(tag) {
  d <- file.path(tempdir(), paste0("ngcd_js_", tag, "_", as.integer(runif(1, 1, 1e6))))
  ng_job_create(d, config = NULL, label = tag, n_traits = 3L)
  d
}

# --- a healthy running job ----------------------------------------------------------------
d <- mk("live")
ng_job_mark(d, "running", list(pid = Sys.getpid()))
ng_job_heartbeat(d)
ng_job_trait_status_write(d, "A", "done", list(cv_predictive_r2 = 0.6, n_selected = 10L))
ng_job_trait_status_write(d, "B", "running")
s <- ng_job_status(d)
stopifnot(identical(s$state, "running"))
stopifnot(identical(as.integer(s$n_done), 1L))
stopifnot(identical(as.integer(s$n_running), 1L))
stopifnot(is.data.frame(s$traits), nrow(s$traits) == 2L)
stopifnot(all(c("trait", "state", "cv_predictive_r2") %in% names(s$traits)))
unlink(d, recursive = TRUE)

# --- a job whose process died --------------------------------------------------------------
d <- mk("dead")
ng_job_mark(d, "running", list(pid = 999999L))
ng_job_heartbeat(d)
Sys.setFileTime(file.path(d, "heartbeat"), Sys.time() - 3600)  # an hour of silence
ng_job_trait_status_write(d, "A", "done", list(cv_predictive_r2 = 0.6, n_selected = 10L))
ng_job_trait_status_write(d, "B", "running")                   # died mid-trait
s <- ng_job_status(d, stale_after_sec = 120)
stopifnot(identical(s$state, "crashed"))
# THE RULE: a trait cannot still be "running" under a job that is not.
b <- s$traits[s$traits$trait == "B", , drop = FALSE]
stopifnot(identical(as.character(b$state), "incomplete"))
# The finished trait keeps its result. A crash must not erase what completed.
a <- s$traits[s$traits$trait == "A", , drop = FALSE]
stopifnot(identical(as.character(a$state), "done"))
stopifnot(isTRUE(all.equal(as.numeric(a$cv_predictive_r2), 0.6)))
unlink(d, recursive = TRUE)

# --- a finished job is never reinterpreted ---------------------------------------------------
# Its heartbeat is stale by definition, because it stopped when the work did.
d <- mk("done")
ng_job_mark(d, "running", list(pid = Sys.getpid()))
ng_job_heartbeat(d)
Sys.setFileTime(file.path(d, "heartbeat"), Sys.time() - 86400)
ng_job_mark(d, "finished")
ng_job_trait_status_write(d, "A", "done", list(n_selected = 4L))
s <- ng_job_status(d)
stopifnot(identical(s$state, "finished"))
unlink(d, recursive = TRUE)

# --- a job that never started ----------------------------------------------------------------
d <- mk("queued")
s <- ng_job_status(d)
stopifnot(identical(s$state, "queued"))
stopifnot(identical(as.integer(s$n_done), 0L))
unlink(d, recursive = TRUE)

cat("PASS: job_status_derives_crashed_and_incomplete\n")
