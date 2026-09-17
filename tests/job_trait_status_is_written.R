# A trait's state must be readable from disk WHILE the batch runs, and after a trait fails.
#
# The state cannot be inferred from result.json: that file appears only at the END, so a
# running trait would look identical to a queued one, and a FAILED trait writes no
# result.json at all -- so the one state a breeder most needs to see would be invisible.
#
# It is also deliberately SMALL. The overview table needs R2 and the gate verdict for every
# trait; reading those out of seventeen result.json files, each carrying ~125,000 candidate
# crosses, would make the list view as costly as the thing it exists to avoid opening.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

job_dir <- file.path(tempdir(), paste0("ngcd_ts_", as.integer(runif(1, 1, 1e6))))
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(job_dir, recursive = TRUE), add = TRUE)

# --- entry: the trait announces itself before doing any work -----------------------------
p <- ng_job_trait_status_write(job_dir, "YIELD", "running")
stopifnot(file.exists(p))
st <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(st$schema, "ng_job_trait_status.v1"))
stopifnot(identical(st$trait, "YIELD"))
stopifnot(identical(st$state, "running"))
stopifnot(is.character(st$started_at), nzchar(st$started_at))

# --- exit: the summary a reader needs, without opening result.json -----------------------
ng_job_trait_status_write(job_dir, "YIELD", "done", list(
  cv_predictive_r2 = 0.61, mean_source = "gebv", effect_gate = "on", n_selected = 20L))
st <- jsonlite::fromJSON(file.path(job_dir, "YIELD", "status.json"), simplifyVector = TRUE)
stopifnot(identical(st$state, "done"))
stopifnot(isTRUE(all.equal(st$cv_predictive_r2, 0.61)))
stopifnot(identical(st$mean_source, "gebv"), identical(st$effect_gate, "on"))
stopifnot(identical(as.integer(st$n_selected), 20L))
stopifnot(is.character(st$finished_at), nzchar(st$finished_at))
# started_at must SURVIVE the second write -- elapsed time is otherwise unknowable, and a
# reader cannot tell a trait that took ten seconds from one that took ten hours.
stopifnot(is.character(st$started_at), nzchar(st$started_at))

# --- failure carries its reason ----------------------------------------------------------
ng_job_trait_status_write(job_dir, "EMPTY", "running")
ng_job_trait_status_write(job_dir, "EMPTY", "error", list(error_message = "no finite values"))
st <- jsonlite::fromJSON(file.path(job_dir, "EMPTY", "status.json"), simplifyVector = TRUE)
stopifnot(identical(st$state, "error"))
stopifnot(grepl("no finite values", st$error_message, fixed = TRUE))

# --- a bad state is refused, not written --------------------------------------------------
# "queued" and "finished" belong to the JOB, not a trait. Accepting them here would let two
# vocabularies drift into one file and make the reader guess which it is looking at.
bad <- tryCatch({ ng_job_trait_status_write(job_dir, "X", "finished"); NULL },
                error = function(e) conditionMessage(e))
stopifnot(!is.null(bad), grepl("finished", bad, fixed = TRUE))

cat("PASS: job_trait_status_is_written\n")
