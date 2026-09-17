# A job must say what it is before it starts, and prove it is alive while it runs.
#
# "Alive" cannot be asked of the process directly: checking a pid means `ps` on Unix and
# `tasklist` on Windows, and this package has a hard requirement to behave identically on
# both. So the running process leaves a heartbeat, and a reader interprets its staleness.
# The pid is recorded for a human debugging afterwards, never for the decision.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

job_dir <- file.path(tempdir(), paste0("ngcd_jr_", as.integer(runif(1, 1, 1e6))))
on.exit(unlink(job_dir, recursive = TRUE), add = TRUE)

# --- created before anything runs ---------------------------------------------------------
p <- ng_job_create(job_dir, config = list(seed = 7L), label = "17-trait run", n_traits = 17L)
stopifnot(file.exists(p))
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$schema, "ng_job.v1"))
stopifnot(identical(rec$state, "queued"))
stopifnot(identical(rec$label, "17-trait run"))
stopifnot(identical(as.integer(rec$n_traits), 17L))
stopifnot(is.character(rec$id), nzchar(rec$id))
stopifnot(is.character(rec$created_at), nzchar(rec$created_at))
created_before <- rec$created_at

# The config travels with the job, so a run is replayable from its own directory rather than
# from whatever the submitting session happened to still hold.
stopifnot(file.exists(file.path(job_dir, "config.json")))

# --- state transitions --------------------------------------------------------------------
# ng_job__now() stamps milliseconds, so a re-stamping mutant would almost always differ
# even without a pause. The sleep stays anyway: "almost always" is not an assertion, and a
# coarse filesystem clock or a virtualised timer can still collapse two calls made
# microseconds apart. Cheap insurance for a check whose whole value is being exact.
Sys.sleep(1.1)
ng_job_mark(job_dir, "running", list(pid = Sys.getpid()))
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$state, "running"))
stopifnot(identical(as.integer(rec$pid), Sys.getpid()))
# created_at must survive -- it is how a listing orders jobs. Equality against the value
# captured right after ng_job_create() is what makes this a real assertion: checking only
# is.character()/nzchar() would pass identically for an implementation that stamped a FRESH
# timestamp on every mark() call, since a fresh timestamp is also a non-empty string.
stopifnot(identical(rec$created_at, created_before))

ng_job_mark(job_dir, "finished")
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$state, "finished"))
stopifnot(is.character(rec$finished_at), nzchar(rec$finished_at))

# --- a bad state is refused ----------------------------------------------------------------
# "crashed" is DERIVED, never written: the process that would have written it is the one
# that died. Accepting it here would invite code to write a state it cannot honestly know.
bad <- tryCatch({ ng_job_mark(job_dir, "crashed"); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(bad), grepl("derived", bad, fixed = TRUE))

# --- heartbeat ------------------------------------------------------------------------------
hb <- file.path(job_dir, "heartbeat")
ng_job_heartbeat(job_dir)
stopifnot(file.exists(hb))
t1 <- file.mtime(hb)
Sys.sleep(1.1)
ng_job_heartbeat(job_dir)
stopifnot(file.mtime(hb) > t1)   # it must actually advance, not merely exist

cat("PASS: job_record_and_heartbeat\n")
