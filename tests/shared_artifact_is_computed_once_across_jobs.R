# The same data must be processed once, not once per job.
#
# Quality control, genotype cleaning, duplicate detection, LD PRUNING, the training-set
# alignment, the pair table and the GRM depend only on the genotypes, the map, and the
# settings that govern them -- never on which trait is being scored. 0.36.0 computes them
# once per batch; that is still once too often. A breeder who submits a second batch on the
# same data pays for all of it again, and LD pruning over a real marker panel is not cheap.
#
# The dangerous half of a cache is the KEY. A key that omitted an input a job depends on
# would silently reuse an artefact built from different data -- wrong numbers, no error.
# That is the same class of defect as the RNG-kind bug found in 0.36.0, which produced a
# different crossing plan from identical data and signalled nothing. So the safety
# assertions below matter more than the reuse assertion.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_sa_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(31); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
gfile <- file.path(root, "geno.csv"); pfile <- file.path(root, "pheno.csv")
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2)), pfile, row.names = FALSE)

cfg <- list(
  genotype_file = gfile, genotype_id_col = "NAME",
  phenotype_file = pfile, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L, ld_pruning = TRUE, grm_method = "vanraden")

k1 <- ng_shared_artifact_key(cfg)
stopifnot(is.character(k1), length(k1) == 1L, nzchar(k1))

# --- STABILITY: the same inputs give the same key -----------------------------------------
stopifnot(identical(ng_shared_artifact_key(cfg), k1))
# Settings that are NOT shared must not move it -- otherwise nothing would ever be reused,
# since per-trait settings are exactly what varies between jobs.
#
# This is a statement about the KEY ONLY. It says n_crosses does not invalidate the cache;
# it says nothing about the batch HONOURING n_crosses, and for a while the two came apart:
# the artefact was the whole ctx, so a cache hit restored the first batch's n_crosses along
# with its shared work. That half is asserted in tests/shared_artifact_does_not_leak_config.R,
# which runs a second batch with different non-shared settings and checks the plan it
# actually produced. Read the two together -- neither alone is the property that matters.
c2 <- cfg; c2$n_crosses <- 9L; c2$min_cv_predictive_r2 <- 0.5
stopifnot(identical(ng_shared_artifact_key(c2), k1))

# --- SAFETY: anything the artefact depends on must move it ---------------------------------
for (chg in list(
      function(x) { x$grm_method <- "yang"; x },
      function(x) { x$ld_pruning <- FALSE; x },
      function(x) { x$duplicate_threshold <- 0.97; x })) {
  stopifnot(!identical(ng_shared_artifact_key(chg(cfg)), k1))
}
# One byte of one input file must move it. Keying on a filename, or a name and a timestamp,
# would let edited data silently reuse an artefact built from the old data.
g2 <- g; g2[1, 1] <- if (g[1, 1] == 0L) 2L else 0L
write.csv(data.frame(NAME = ids, g2, check.names = FALSE), gfile, row.names = FALSE)
stopifnot(!identical(ng_shared_artifact_key(cfg), k1))
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
stopifnot(identical(ng_shared_artifact_key(cfg), k1))   # and restoring it restores the key

# --- SAFETY: an in-memory input must be hashed by CONTENT, not by shape --------------------
# Summarising with str() or dim() would be cheaper and would be a latent disaster: two
# different genotype matrices of the same shape would hash identically and silently reuse
# each other's artefact.
cm <- cfg
cm$marker_map$Position_cM[[1L]] <- cm$marker_map$Position_cM[[1L]] + 7   # same shape, new data
stopifnot(!identical(ng_shared_artifact_key(cm), k1))

# --- REUSE: a second batch on the same data does not redo the work -------------------------
shared_dir <- file.path(root, "shared")
qc_calls <- 0L
orig <- ng_cp__stage_qc
assign("ng_cp__stage_qc", function(...) { qc_calls <<- qc_calls + 1L; orig(...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__stage_qc", orig, envir = .GlobalEnv), add = TRUE)

b1 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = file.path(root, "out1"), batch_workers = 1L, shared_dir = shared_dir))
stopifnot(identical(qc_calls, 1L))

b2 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = file.path(root, "out2"), batch_workers = 1L, shared_dir = shared_dir))
stopifnot(identical(qc_calls, 1L))   # the second batch ran NO quality control

# The second batch still produced real results, not an empty shell.
for (j in b2$jobs) stopifnot(identical(j$status, "ok"))
stopifnot(file.exists(file.path(root, "out2", "A", "result.json")))

# --- the artefact records who is using it, so retention can be safe ------------------------
# Locate the artefact by listing shared_dir rather than assuming its key is k1: the runner
# keys on the RESOLVED config (after back-fill), while k1 above was computed directly on the
# raw cfg as a pure-function stability/safety check, so the two need not be the same string.
# Asserting there is exactly ONE artefact directory is itself a stronger check than looking
# up k1 would have been -- it proves the second batch did not create a second artefact.
artifact_dirs <- list.dirs(shared_dir, recursive = FALSE, full.names = TRUE)
stopifnot(length(artifact_dirs) == 1L)
meta <- jsonlite::fromJSON(file.path(artifact_dirs[[1L]], "meta.json"), simplifyVector = TRUE)
stopifnot(identical(meta$schema, "ng_shared_artifact.v1"))
stopifnot(length(meta$referenced_by) >= 1L)

cat("PASS: shared_artifact_is_computed_once_across_jobs\n")
