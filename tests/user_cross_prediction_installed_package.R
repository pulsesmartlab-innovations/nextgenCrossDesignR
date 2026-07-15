test_lib <- Sys.getenv("NG_TEST_INSTALLED_LIB", unset = "")
if (!nzchar(test_lib)) {
  cat("installed package user workflow test skipped; set NG_TEST_INSTALLED_LIB to run\n")
  quit(status = 0)
}

.libPaths(c(
  normalizePath(test_lib, winslash = "/", mustWork = TRUE),
  if (.Platform$OS.type == "windows") normalizePath(file.path(getwd(), "..", ".Rlib"), winslash = "/", mustWork = FALSE) else NULL,
  .libPaths()
))

library(nextgenCrossDesign)

tmp <- tempfile("ng_installed_user_workflow_")
dir.create(tmp, recursive = TRUE)

ids <- c("P01", "P02", "P02_copy", "P03", "P04", "P05", "P06")
geno <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 0, 2, 2, 0, 2),
  M02 = c(0, 2, 2, 0, 2, 0, 2),
  M03 = c(2, 0, 0, 2, 0, 2, 0),
  M04 = c(2, 2, 2, 0, 0, 2, 0),
  M05 = c(0, 0, 0, 2, 0, 2, 2),
  M06 = c(2, 0, 0, 0, 2, 2, 0),
  M07 = c(0, 2, 2, 2, 0, 0, 2),
  M08 = c(2, 2, 2, 0, 2, 0, 0),
  check.names = FALSE
)
phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 61, 61, 53, 67, 55, 64),
  disease = c(4.0, 3.1, 3.1, 5.2, 2.6, 4.6, 2.9),
  stringsAsFactors = FALSE
)
marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:8)),
  Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
  Position_BP = c(0, 1, 2, 3, 0, 1, 2, 3) * 1e6,
  stringsAsFactors = FALSE
)
direction <- data.frame(
  trait = c("yield", "disease"),
  column = c("yield", "disease"),
  direction = c("increase", "decrease"),
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(tmp, "pheno.csv")
genotype_file <- file.path(tmp, "geno.csv")
map_file <- file.path(tmp, "map.csv")
direction_file <- file.path(tmp, "direction.csv")
write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(geno, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(direction, direction_file, row.names = FALSE, quote = FALSE)

result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  id_col = "NAME",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_col = "Position_BP",
  map_pos_cm_divisor = 1e6,
  prediction_mode = "trait_by_trait",
  trait_value_metric = "uc",
  uc_variance_source = "pmv",
  multi_trait_method = "auto",
  progeny = "DH",
  duplicate_action = "remove",
  duplicate_threshold = 0.995,
  duplicate_maf_min = 0,
  duplicate_max_missing_prop = 1,
  duplicate_min_compared_markers = 8,
  n_crosses = 5,
  max_uses_per_parent = 3,
  optimizer = "auto",
  use_ocs = TRUE,
  output_dir = file.path(tmp, "outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,
  assume_inbred = TRUE,
  seed = 11
)

stopifnot(packageVersion("nextgenCrossDesign") >= "0.3.7")
stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(nrow(result$selected_crosses) == 5L)
stopifnot("putative_duplicate_genotypes_removed" %in% result$qc$issues$id)
stopifnot(file.exists(result$output_files$priority_score_vs_kinship_png))
stopifnot(file.info(result$output_files$priority_score_vs_kinship_png)$size > 1000)

cat("installed package user workflow test passed\n")
