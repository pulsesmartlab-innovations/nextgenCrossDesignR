# traits_to_use must key the artefact even though a job may override it.
#
# ng_cp__batch_shared_keys and ng_cp__batch_artifact_keys (R/54) are deliberately different
# sets. traits_to_use is EXCLUDED from the first -- a per-trait job's whole purpose is to
# override it -- and INCLUDED in the second, because ng_cp__stage_qc() bakes it into
# trait_spec (ng_run_cp_trait_spec(), R/39), and trait_spec is a field the artefact persists
# (ng_shared_artifact_fields, R/55).
#
# The defect this prevents: two batches on the SAME genotypes differing only in
# traits_to_use would, before this fix, compute the SAME key -- ng_shared_artifact_key()
# iterated only ng_cp__batch_shared_keys. The second batch would then reuse the first
# batch's artefact and silently inherit the FIRST batch's trait_spec. Run traits A,B and
# then A,C, and the second batch would score B instead of C -- no error, no warning, a
# crossing plan for the wrong trait.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_satc_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(53); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
gfile <- file.path(root, "geno.csv"); pfile <- file.path(root, "pheno.csv")
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2), C = tr(3)), pfile, row.names = FALSE)

base_cfg <- list(
  genotype_file = gfile, genotype_id_col = "NAME",
  phenotype_file = pfile, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B", "C"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L, ld_pruning = TRUE, grm_method = "vanraden")

cfg1 <- c(base_cfg, list(traits_to_use = c("A", "B")))
cfg2 <- c(base_cfg, list(traits_to_use = c("A", "C")))
cfg3 <- cfg2   # identical to batch 2 -- the reuse case

# --- THE KEYS MUST DIFFER -------------------------------------------------------------
k1 <- ng_shared_artifact_key(cfg1)
k2 <- ng_shared_artifact_key(cfg2)
stopifnot(is.character(k1), nzchar(k1), is.character(k2), nzchar(k2))
stopifnot(!identical(k1, k2))

shared_dir <- file.path(root, "shared")
qc_calls <- 0L
orig <- ng_cp__stage_qc
assign("ng_cp__stage_qc", function(...) { qc_calls <<- qc_calls + 1L; orig(...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__stage_qc", orig, envir = .GlobalEnv), add = TRUE)

job_ids <- function(b) vapply(b$jobs, function(j) as.character(j$id), character(1))

b1 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg1, output_root = file.path(root, "out1"), batch_workers = 1L, shared_dir = shared_dir))
for (j in b1$jobs) if (!identical(j$status, "ok")) stop("batch1 ", j$id, ": ", j$error_message)
stopifnot(identical(qc_calls, 1L))
stopifnot(setequal(job_ids(b1), c("A", "B")))
stopifnot(file.exists(file.path(root, "out1", "A", "result.json")))
stopifnot(file.exists(file.path(root, "out1", "B", "result.json")))

# --- BATCH 2: same genotypes, different traits_to_use -----------------------------------
# A different key means quality control runs again -- this IS the fix. Before it, this
# second call would have hit batch 1's artefact and inherited batch 1's trait_spec (A, B).
b2 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg2, output_root = file.path(root, "out2"), batch_workers = 1L, shared_dir = shared_dir))
for (j in b2$jobs) if (!identical(j$status, "ok")) stop("batch2 ", j$id, ": ", j$error_message)
stopifnot(identical(qc_calls, 2L))   # QC recomputed -- a genuinely different artefact

# Batch 2 asked for A and C. Its own job manifest must name exactly those traits.
stopifnot(setequal(job_ids(b2), c("A", "C")))

# Batch 2 must have produced results for ITS OWN traits...
stopifnot(file.exists(file.path(root, "out2", "A", "result.json")))
stopifnot(file.exists(file.path(root, "out2", "C", "result.json")))
# ...and NOT for a trait only batch 1 requested. This is the exact failure mode: batch 2
# silently scoring trait B (batch 1's second trait) instead of its own trait C.
stopifnot(!dir.exists(file.path(root, "out2", "B")))

read_job <- function(out, trait) {
  p <- file.path(root, out, trait, "result.json")
  stopifnot(file.exists(p))
  jsonlite::fromJSON(p, simplifyVector = TRUE)
}
# Each job's result must be keyed to the trait its OWN directory name promises -- the
# trait_spec the run actually scored against (ng_cp__assemble_result() ->
# trait_direction <- trait_spec, R/39; trait_spec is the field the artefact key must
# cover, since it is what value_kind/direction/weight for a trait ultimately come from).
r2a <- read_job("out2", "A")
r2c <- read_job("out2", "C")
stopifnot(identical(as.character(r2a$trait_direction$trait), "A"))
stopifnot(identical(as.character(r2c$trait_direction$trait), "C"))
# And never the trait that only batch 1 asked for.
stopifnot(!("B" %in% as.character(r2a$trait_direction$trait)))
stopifnot(!("B" %in% as.character(r2c$trait_direction$trait)))

# --- BATCH 3: identical to batch 2 -- the reuse case must still work ---------------------
b3 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg3, output_root = file.path(root, "out3"), batch_workers = 1L, shared_dir = shared_dir))
for (j in b3$jobs) if (!identical(j$status, "ok")) stop("batch3 ", j$id, ": ", j$error_message)
stopifnot(identical(qc_calls, 2L))   # NO third quality-control run: the artefact was reused
stopifnot(setequal(job_ids(b3), c("A", "C")))
stopifnot(file.exists(file.path(root, "out3", "A", "result.json")))
stopifnot(file.exists(file.path(root, "out3", "C", "result.json")))

# Exactly two artefacts exist: one for {A,B}, one for {A,C} shared by batches 2 and 3.
artifact_dirs <- list.dirs(shared_dir, recursive = FALSE, full.names = TRUE)
stopifnot(length(artifact_dirs) == 2L)

cat("PASS: shared_artifact_key_covers_the_trait_set\n")
