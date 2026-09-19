# This test lives in tests/testthat/, not tests/, ON PURPOSE.
#
# Only tests/testthat/ runs against the INSTALLED package -- R CMD check installs the
# package first and then attaches it with library(), exactly as every production caller
# does. Everything under tests/*.R runs against the SOURCE TREE via ng_load(), which sources
# every internal function into globalenv() and so can never see a defect that only exists
# when the package is installed. ng_run_cross_prediction_batch(batch_workers >= 2) is exactly
# that defect: its mirai::everywhere() setup block called an internal function
# (ng_shared_artifact_read) that library() does not expose, so every worker daemon's setup
# silently failed and every job came back "no batch configuration reached the worker". The
# 183/183 deep harness was green throughout, because none of its eight batch tests could ever
# see this. This test is the one place that exercises the installed path with
# batch_workers >= 2 and would have caught the regression before release.
test_that("a batch with batch_workers = 2 runs against the installed package", {
  skip_if_not_installed("mirai")

  set.seed(1101)
  n <- 16L; m <- 16L
  geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
  ids <- paste0("P", seq_len(n))
  rownames(geno) <- ids
  colnames(geno) <- paste0("M", seq_len(m))
  trait_of <- function(s) {
    set.seed(s)
    b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
    as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.3)
  }
  traits <- c("YIELD", "PROTEIN")
  pheno <- data.frame(NAME = ids, YIELD = trait_of(11L), PROTEIN = trait_of(12L),
                      stringsAsFactors = FALSE)
  gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
  map <- data.frame(SNP_code = colnames(geno), Chromosome = rep(1:2, length.out = m),
                    Position_cM = rep(seq(0, 100, length.out = m / 2), times = 2),
                    stringsAsFactors = FALSE)
  dir <- data.frame(Trait = traits, Selection_direction = "increase",
                    stringsAsFactors = FALSE)

  config <- list(
    genotype = gdf, genotype_id_col = "NAME",
    phenotype = pheno, phenotype_id_col = "NAME",
    trait_direction = dir, direction_trait_col = "Trait",
    direction_column_col = "Trait", direction_direction_col = "Selection_direction",
    marker_map = map, map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_cm_col = "Position_cM", map_position_unit = "cM",
    progeny = "RIL", parent_type = "ril", n_crosses = 3L,
    write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
    min_cv_predictive_r2 = -1, seed = 20260918L)

  root <- file.path(tempdir(), paste0("ngcd_batch_installed_", as.integer(runif(1, 1, 1e6))))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  b <- suppressWarnings(ng_run_cross_prediction_batch(config, output_root = root,
                                                      batch_workers = 2L))

  expect_s3_class(b, "ng_cross_prediction_batch")
  expect_identical(length(b$jobs), 2L)
  statuses <- vapply(b$jobs, function(j) j$status, character(1))
  errs <- vapply(b$jobs, function(j) {
    if (is.null(j$error_message)) "" else j$error_message
  }, character(1))
  expect_true(all(statuses == "ok"),
              info = paste("job statuses:", paste(statuses, collapse = ", "),
                          "| errors:", paste(errs, collapse = " || ")))
})
