# A cache hit must reuse the shared WORK, never the previous batch's SETTINGS.
#
# The shared artefact is keyed on ng_cp__batch_shared_keys -- deliberately only the settings
# quality control and the predict prologue actually spend. Everything else (n_crosses, the
# scoring metric, the acceptability bar, the seed, where outputs go) is per-batch and must
# come from the batch that is running NOW.
#
# The defect this guards: the artefact used to be the whole ctx, and ctx IS the config
# (ng_cp__build_ctx() starts with `ctx <- config`). saveRDS(ctx)/readRDS() therefore restored
# every non-shared setting from whichever batch built the artefact FIRST. A second batch
# asking for nine crosses scored on the mid-parent mean silently got the first batch's three
# crosses scored on usefulness. No error -- the spec's own headline scenario ("a second batch
# on the same data -- a different metric, a different acceptability bar") producing wrong
# numbers silently, which is the same defect class as the 0.36.0 RNG-kind bug.
#
# So this asserts BOTH halves at once: the artefact is reused (quality control runs once,
# counted the way tests/shared_artifact_is_computed_once_across_jobs.R counts it) AND the
# second batch's own settings are what produced its numbers.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_sl_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(77); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
gfile <- file.path(root, "geno.csv"); pfile <- file.path(root, "pheno.csv")
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2)), pfile, row.names = FALSE)

base_cfg <- list(
  genotype_file = gfile, genotype_id_col = "NAME",
  phenotype_file = pfile, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril",
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L, ld_pruning = TRUE, grm_method = "vanraden")

# Every difference between these two is a NON-shared setting: none of them appears in
# ng_cp__batch_shared_keys, so both batches must hit the same artefact.
cfg1 <- c(base_cfg, list(n_crosses = 3L, trait_value_metric = "usefulness",
                         uc_variance_source = "vpm", selection_prop = 0.10))
cfg2 <- c(base_cfg, list(n_crosses = 9L, trait_value_metric = "mean",
                         uc_variance_source = "vpm", selection_prop = 0.25))
stopifnot(identical(ng_shared_artifact_key(cfg1), ng_shared_artifact_key(cfg2)))

shared_dir <- file.path(root, "shared")
qc_calls <- 0L
orig <- ng_cp__stage_qc
assign("ng_cp__stage_qc", function(...) { qc_calls <<- qc_calls + 1L; orig(...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__stage_qc", orig, envir = .GlobalEnv), add = TRUE)

b1 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg1, output_root = file.path(root, "out1"), batch_workers = 1L, shared_dir = shared_dir))
stopifnot(identical(qc_calls, 1L))

b2 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg2, output_root = file.path(root, "out2"), batch_workers = 1L, shared_dir = shared_dir))
# THE REUSE HALF: the second batch ran no quality control at all.
stopifnot(identical(qc_calls, 1L))

read_job <- function(out, trait) {
  p <- file.path(root, out, trait, "result.json")
  stopifnot(file.exists(p))
  jsonlite::fromJSON(p, simplifyVector = TRUE)
}
r1 <- read_job("out1", "A")
r2 <- read_job("out2", "A")

# THE NO-LEAK HALF: batch 2's settings, not batch 1's.
#
# settings is the RESOLVED echo (R/39 ng_cp__assemble_result), so it reports what actually
# produced the numbers rather than what was typed.
stopifnot(identical(as.character(r1$settings$trait_value_metric), "usefulness"))
stopifnot(identical(as.character(r2$settings$trait_value_metric), "mean"))
stopifnot(identical(as.integer(r1$settings$n_crosses), 3L))
stopifnot(identical(as.integer(r2$settings$n_crosses), 9L))
stopifnot(identical(as.numeric(r2$settings$selection_prop), 0.25))

# The scored column, not just the echoed setting: variance_column_used names the column
# ng_run_cp_trait_value() actually READ. Under usefulness+vpm that is "vpm"; under the
# mid-parent mean there is no variance column at all. A leaked config would make batch 2
# report "vpm" here, because it would have been scored as a usefulness run.
vcol <- function(r) as.character(r$effect_summary$variance_column_used)[[1L]]
stopifnot(identical(vcol(r1), "vpm"))
stopifnot(is.na(vcol(r2)))

# And the deliverable itself: the plan batch 2 hands the breeder is nine crosses, not three.
stopifnot(identical(nrow(r1$selected_crosses), 3L))
stopifnot(identical(nrow(r2$selected_crosses), 9L))

# The per-trait status file -- what an overview reads -- agrees with the plan on disk.
st2 <- jsonlite::fromJSON(file.path(root, "out2", "A", "status.json"), simplifyVector = TRUE)
stopifnot(identical(as.integer(st2$n_selected), 9L))

# Still exactly one artefact: honouring the caller's settings must not have been achieved by
# quietly building a second copy of the shared work.
stopifnot(length(list.dirs(shared_dir, recursive = FALSE, full.names = TRUE)) == 1L)

cat("PASS: shared_artifact_does_not_leak_config\n")
