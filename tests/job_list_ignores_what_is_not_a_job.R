# A listing must survive the directory containing things that are not jobs.
#
# The jobs directory is a shared location on disk. Editors leave dotfiles, users drop
# folders in, an interrupted submit leaves a half-made directory with no job.json. A listing
# that errored on any of those would take the whole Jobs view down with it -- and the view
# is the only way a breeder reaches results at all.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_jl_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

for (nm in c("job_a", "job_b", "job_c")) {
  d <- file.path(root, nm)
  ng_job_create(d, config = NULL, label = paste("run", nm), n_traits = 2L)
  ng_job_mark(d, "finished")
  ng_job_trait_status_write(d, "A", "done", list(n_selected = 3L))
  Sys.sleep(1.1)   # distinct created_at, so ordering is testable rather than tied
}
# Things that are not jobs.
dir.create(file.path(root, "not_a_job"), showWarnings = FALSE)          # no job.json
dir.create(file.path(root, ".DS_Store_dir"), showWarnings = FALSE)
writeLines("x", file.path(root, "stray.txt"))                            # a file, not a dir

lst <- ng_job_list(root)
stopifnot(is.data.frame(lst))
stopifnot(identical(nrow(lst), 3L))                       # exactly the three real jobs
stopifnot(all(c("id", "label", "state", "n_done", "path") %in% names(lst)))
stopifnot(!("not_a_job" %in% lst$id), !("stray.txt" %in% lst$id))
# Newest first -- a breeder looking for what they just submitted should not have to scroll.
stopifnot(identical(lst$id, c("job_c", "job_b", "job_a")))
stopifnot(all(lst$state == "finished"))
stopifnot(all(as.integer(lst$n_done) == 1L))

# Truncation test: verify limit works on created_at, not mtime.
# Create a 4th job (newest), then exercise limit=2 twice.
d_d <- file.path(root, "job_d")
ng_job_create(d_d, config = NULL, label = "run job_d", n_traits = 2L)
ng_job_mark(d_d, "finished")
ng_job_trait_status_write(d_d, "A", "done", list(n_selected = 3L))

# With 4 jobs, limit=2 should return only the 2 newest (job_d, job_c).
lst_limited <- ng_job_list(root, limit = 2L)
stopifnot(identical(nrow(lst_limited), 2L))
stopifnot(identical(lst_limited$id, c("job_d", "job_c")))

# Now mark the OLDEST job (job_a) to bump its directory mtime.
# This tests that we sort by creation time, not mtime.
ng_job_mark(file.path(root, "job_a"), "finished")

# limit=2 should STILL return job_d and job_c, not the recently-touched job_a.
lst_after_mark <- ng_job_list(root, limit = 2L)
stopifnot(identical(nrow(lst_after_mark), 2L))
stopifnot(identical(lst_after_mark$id, c("job_d", "job_c")))
stopifnot(!("job_a" %in% lst_after_mark$id))

# An empty or absent directory is a normal state on a fresh install, not an error.
stopifnot(identical(nrow(ng_job_list(file.path(root, "nope"))), 0L))
empty <- file.path(root, "empty"); dir.create(empty)
stopifnot(identical(nrow(ng_job_list(empty)), 0L))

cat("PASS: job_list_ignores_what_is_not_a_job\n")
