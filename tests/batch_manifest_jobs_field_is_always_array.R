# The manifest's `jobs` field must be a JSON ARRAY regardless of worker count.
#
# ng_cp__batch_dispatch() returns a NAMED list on its serial path (lapply() over a named
# `jobs` list) and an UNNAMED list on its mirai path (mirai_map()'s result is walked by
# seq_along(), which drops names). jsonlite then serialises manifest.json's `jobs` field as
# a JSON OBJECT keyed by trait at batch_workers = 1 and a JSON ARRAY at batch_workers = 2.
# Verified directly on 0.36.0. A consumer -- and the entire next round is a frontend reading
# this file -- would have to handle both shapes and would break on whichever it did not
# test. This asserts the shape AT THE SERIALISED JSON LEVEL, not merely on the R object's
# class: the defect is in what reaches disk, and an R-side-only assertion (e.g. is.null(
# names(x))) would miss jsonlite's own auto_unbox/dataframe rules re-adding structure on the
# way out.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(414); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
pheno <- data.frame(NAME = ids, A = tr(1), B = tr(2))
cfg <- list(
  genotype = data.frame(NAME = ids, g, check.names = FALSE), genotype_id_col = "NAME",
  phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 45, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L)

run_and_read_shape <- function(tag, workers) {
  root <- file.path(tempdir(), paste0("ngcd_mj_", tag))
  unlink(root, recursive = TRUE)
  b <- suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root = root,
                                                      batch_workers = workers))
  for (j in b$jobs) if (!identical(j$status, "ok")) stop(tag, "/", j$id, ": ", j$error_message)
  manifest_path <- file.path(root, "manifest.json")
  stopifnot(file.exists(manifest_path))
  # simplifyVector = FALSE preserves the JSON-level distinction: a JSON array parses to an
  # UNNAMED R list, a JSON object to a NAMED one. This is the check that reads what actually
  # landed on disk, not what the R object looked like before jsonlite touched it.
  parsed <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  is_array <- is.null(names(parsed$jobs))
  unlink(root, recursive = TRUE)
  is_array
}

array_at_1 <- run_and_read_shape("w1", 1L)
array_at_2 <- run_and_read_shape("w2", 2L)

stopifnot(isTRUE(array_at_1))
stopifnot(isTRUE(array_at_2))
stopifnot(identical(array_at_1, array_at_2))

cat("PASS: batch_manifest_jobs_field_is_always_array\n")
