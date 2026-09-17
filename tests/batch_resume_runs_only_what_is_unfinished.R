# Resuming must finish a job, not restart it.
#
# Failure isolation (0.36.0) means one bad trait does not cost the other sixteen. But
# "the manifest records that it failed" is only half an answer; the other half is being able
# to finish. A batch may have run for hours, so re-running the traits that already succeeded
# would be its own kind of loss.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_rs_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(77); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
# EMPTY has no usable data, so its job fails -- and only its job.
pheno <- data.frame(NAME = ids, A = tr(1), B = tr(2), EMPTY = rep(NA_real_, n))
cfg <- list(
  genotype = data.frame(NAME = ids, g, check.names = FALSE), genotype_id_col = "NAME",
  phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B", "EMPTY"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  # write_outputs = TRUE (not FALSE): FIX B needs an actual workbook path on disk for a
  # carried-forward entry to reproduce -- see the output_files assertions below.
  write_outputs = TRUE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L)

out <- file.path(root, "out")
b1 <- suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root = out, batch_workers = 1L))
st <- vapply(b1$jobs, function(j) j$status, character(1))
names(st) <- vapply(b1$jobs, function(j) j$id, character(1))
stopifnot(identical(unname(st[c("A", "B")]), c("ok", "ok")))
stopifnot(identical(unname(st["EMPTY"]), "error"))

# FIX B: a freshly-run trait's manifest entry carries output_files -- the workbook path a
# reader needs to find the deliverable. Record it so we can prove resume does not drop it.
a_entry_fresh <- Filter(function(j) identical(j$id, "A"), b1$jobs)[[1L]]
stopifnot(is.list(a_entry_fresh$output_files))
stopifnot(is.character(a_entry_fresh$output_files$workbook))
stopifnot(file.exists(a_entry_fresh$output_files$workbook))

# Record when A's result was written, so we can prove it is not rewritten.
a_mtime <- file.mtime(file.path(out, "A", "result.json"))
Sys.sleep(1.1)

# --- resume: only the unfinished trait runs -------------------------------------------------
ran <- character(0)
orig <- ng_cp__batch_run_one
assign("ng_cp__batch_run_one", function(job, ...) { ran <<- c(ran, job$id); orig(job, ...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__batch_run_one", orig, envir = .GlobalEnv), add = TRUE)

b2 <- suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root = out,
                                                     batch_workers = 1L, resume = TRUE))
stopifnot(identical(sort(ran), "EMPTY"))                 # A and B were not re-run
stopifnot(identical(file.mtime(file.path(out, "A", "result.json")), a_mtime))

# The manifest still describes ALL traits, not only the resumed one -- a resumed job is the
# same job, and a manifest listing one trait would misrepresent what the breeder has.
ids2 <- sort(vapply(b2$jobs, function(j) j$id, character(1)))
stopifnot(identical(ids2, c("A", "B", "EMPTY")))
carried <- Filter(function(j) identical(j$id, "A"), b2$jobs)[[1L]]
stopifnot(identical(carried$status, "ok"))

# FIX B: a carried-forward entry (built from status.json on resume, not from a fresh run)
# must still name the workbook. Before the fix, output_files was never recorded in
# status.json, so a resumed job silently dropped the deliverable path for every trait that
# did not re-run -- even though the file is right there on disk.
stopifnot(is.list(carried$output_files))
stopifnot(identical(carried$output_files$workbook, a_entry_fresh$output_files$workbook))

cat("PASS: batch_resume_runs_only_what_is_unfinished\n")
