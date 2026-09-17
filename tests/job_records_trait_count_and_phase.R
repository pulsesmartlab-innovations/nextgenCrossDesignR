# "k of n done" needs an n that does not move, and a phase that says what silence means.
#
# TWO defects, both on the REAL path rather than on a hand-made job record, which is why
# this test drives an actual batch instead of calling ng_job_mark() directly:
#
#   n_traits was never populated by ng_run_cross_prediction_batch(), so ng_job_status() fell
#   back to counting the status.json files that existed. The denominator then GREW as workers
#   started: a 17-trait job read "1 of 1 done", then "1 of 2", then "2 of 5". The parent is
#   the only process that knows length(jobs) before any worker has written anything.
#
#   phase was absent entirely. The parent beats the heartbeat once and then runs quality
#   control, the duplicate scan, LD pruning and the GRM before dispatching -- and R is
#   single-threaded, so nothing can touch the heartbeat during it. A reader with only an
#   mtime cannot tell that from a dead process. The phase says which it is without any timing
#   assumption at all.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_nt_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(5); n <- 16L; m <- 16L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
cfg <- list(
  genotype = data.frame(NAME = ids, g, check.names = FALSE), genotype_id_col = "NAME",
  phenotype = data.frame(NAME = ids, A = tr(1), B = tr(2), C = tr(3)),
  phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B", "C"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 45, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, effect_gate = "off", seed = 1L)

job_dir <- file.path(root, "jobs", "job_1")

# Observe the job mid-run, from inside the serial dispatch loop: after the FIRST trait has
# finished and before the other two have. This is the moment the overview actually renders,
# and the only moment at which a grown denominator is visible.
snaps <- list()
orig <- ng_cp__batch_run_one
assign("ng_cp__batch_run_one", function(...) {
  r <- orig(...)
  snaps[[length(snaps) + 1L]] <<- ng_job_status(job_dir)
  r
}, envir = .GlobalEnv)
on.exit(assign("ng_cp__batch_run_one", orig, envir = .GlobalEnv), add = TRUE)

b <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = job_dir, batch_workers = 1L, job_dir = job_dir,
  shared_dir = file.path(root, "shared")))
for (j in b$jobs) stopifnot(identical(j$status, "ok"))

stopifnot(length(snaps) == 3L)

# --- n_traits is the FULL count from the first observation onward ----------------------------
after_first <- snaps[[1L]]
stopifnot(identical(as.integer(after_first$n_traits), 3L))
stopifnot(as.integer(after_first$n_done) < 3L)      # genuinely mid-run
stopifnot(identical(as.integer(after_first$n_done), 1L))
# It does not move as workers start. A fallback that counted status files would have read
# 1, then 2, then 3 here.
stopifnot(identical(vapply(snaps, function(s) as.integer(s$n_traits), integer(1)),
                    c(3L, 3L, 3L)))
stopifnot(identical(vapply(snaps, function(s) as.integer(s$n_done), integer(1)),
                    c(1L, 2L, 3L)))

# --- phase distinguishes shared setup from scoring from finished ------------------------------
stopifnot(identical(as.character(after_first$phase), "scoring"))
final <- ng_job_status(job_dir)
stopifnot(identical(final$state, "finished"))
stopifnot(identical(as.integer(final$n_traits), 3L))
stopifnot(identical(as.character(final$phase), "complete"))

# The shared-setup phase is what a reader sees while the parent cannot beat its heartbeat.
# Asserted on a fresh job record rather than by racing a real one: ng_run_cross_prediction_batch()
# writes it as its first act, and the mid-run snapshots above prove it is later replaced.
d2 <- file.path(root, "jobs", "job_2")
invisible(ng_job_create(d2, config = NULL, label = "job_2"))
invisible(ng_job_mark(d2, "running", list(pid = Sys.getpid(), phase = "shared_setup")))
invisible(ng_job_heartbeat(d2))
s2 <- ng_job_status(d2)
stopifnot(identical(as.character(s2$phase), "shared_setup"))
# No traits yet, and no heartbeat it could possibly have beaten -- which is exactly why the
# default stale window is hours, not two minutes.
stopifnot(identical(s2$state, "running"))
stopifnot(ng_job_stale_after_default >= 3600)

# The listing carries it too -- that is the view a reader polls, and the one where a job
# still inside its shared setup would otherwise be misread as dead.
lst <- ng_job_list(dirname(job_dir))
stopifnot("phase" %in% names(lst))
stopifnot(identical(as.character(lst$phase[lst$id == "job_1"]), "complete"))
stopifnot(identical(as.character(lst$phase[lst$id == "job_2"]), "shared_setup"))

# A job record written before any phase was recorded still reads, as NA rather than an error:
# the registry must tolerate a directory written by an older release.
d3 <- file.path(root, "jobs", "job_3")
invisible(ng_job_create(d3, config = NULL, label = "job_3"))
stopifnot(is.na(ng_job_status(d3)$phase))

cat("PASS: job_records_trait_count_and_phase\n")
