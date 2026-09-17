# A field that is conceptually a LIST must serialise as a JSON array at every length.
#
# jsonlite's auto_unbox = TRUE renders a length-1 vector as a bare scalar and anything longer
# as an array. Every conceptually-plural field written by this package therefore changed shape
# with its own content: manifest.json's jobs[].traits was "YIELD" for a one-trait job and
# ["A","B"] for a two-trait one, and a result envelope's `warnings` was a string after one
# warning and an array after two. A consumer -- and the next round is a frontend reading
# exactly these files -- has to handle both shapes and breaks on whichever it did not test.
# Worse, the DEFAULT batch is one trait per job, so the scalar form is the common case and the
# array the surprise.
#
# This is the same defect the manifest's `jobs` field was fixed for in
# tests/batch_manifest_jobs_field_is_always_array.R, and it is asserted the same way: by
# parsing what LANDED ON DISK with simplifyVector = FALSE, where a JSON array is an unnamed R
# list and a JSON string is a length-1 character. An R-side assertion on the object before
# jsonlite touched it would miss the defect entirely, because the defect is jsonlite's.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_pj_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(101); n <- 16L; m <- 16L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
cfg <- list(
  genotype = data.frame(NAME = ids, g, check.names = FALSE), genotype_id_col = "NAME",
  phenotype = data.frame(NAME = ids, A = tr(1), B = tr(2)), phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 45, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, effect_gate = "off", seed = 1L)

# A JSON array parses to an UNNAMED R list under simplifyVector = FALSE; a JSON string does
# not parse to a list at all. That difference is the whole assertion.
is_json_array <- function(x) is.list(x) && is.null(names(x))

# --- jobs[].traits: one trait per job (the DEFAULT) and two traits in one job ----------------
out1 <- file.path(root, "one_each")
b1 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = out1, batch_workers = 1L, shared_dir = file.path(root, "shared")))
for (j in b1$jobs) stopifnot(identical(j$status, "ok"))
m1 <- jsonlite::fromJSON(file.path(out1, "manifest.json"), simplifyVector = FALSE)
stopifnot(length(m1$jobs) == 2L)
for (j in m1$jobs) {
  stopifnot(is_json_array(j$traits))     # ONE trait, and still an array
  stopifnot(length(j$traits) == 1L)
}

out2 <- file.path(root, "both_in_one")
b2 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, jobs = list(pair = list(id = "pair", traits = c("A", "B"))),
  output_root = out2, batch_workers = 1L, shared_dir = file.path(root, "shared")))
for (j in b2$jobs) stopifnot(identical(j$status, "ok"))
m2 <- jsonlite::fromJSON(file.path(out2, "manifest.json"), simplifyVector = FALSE)
stopifnot(length(m2$jobs) == 1L)
stopifnot(is_json_array(m2$jobs[[1L]]$traits), length(m2$jobs[[1L]]$traits) == 2L)

# --- a result envelope's `warnings`, at zero, one and two -------------------------------------
# Written through ng_write_result_json(), the same writer a job's result.json uses, so this
# asserts the contract the frontend parses rather than a private helper.
for (k in 0:2) {
  p <- file.path(root, paste0("env_", k, ".json"))
  env <- ng_run_result_envelope(list(qc = list(status = "pass")),
                                warnings = if (k == 0L) character(0) else paste0("w", seq_len(k)))
  ng_write_result_json(env, p)
  parsed <- jsonlite::fromJSON(p, simplifyVector = FALSE)
  stopifnot(is_json_array(parsed$warnings))
  stopifnot(length(parsed$warnings) == k)
}

# --- a job's own result.json, written by a real worker ------------------------------------------
rj <- jsonlite::fromJSON(file.path(out1, "A", "result.json"), simplifyVector = FALSE)
stopifnot(is_json_array(rj$warnings))

cat("PASS: batch_plural_json_fields_are_always_arrays\n")
