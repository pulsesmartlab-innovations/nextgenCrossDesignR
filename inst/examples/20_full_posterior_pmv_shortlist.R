# Full-posterior PMV shortlist example.
#
# A practical workflow is to screen broadly with method_varPMV = "fast", then
# rerun a shortlist-sized plan with method_varPMV = "full_posterior" so the PMV
# score uses the full marker-effect covariance.

library(nextgenCrossDesign)

required_version <- "0.3.13"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop(
    "This example requires nextgenCrossDesign ", required_version,
    " or newer. Reinstall the current tarball and restart R.",
    call. = FALSE
  )
}

set.seed(20260623)

out_dir <- file.path(tempdir(), "nextgenCrossDesign_full_posterior_pmv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- paste0("P", sprintf("%02d", 1:8))

genotype <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 2, 2, 0, 2, 0, 2),
  M02 = c(0, 2, 0, 2, 0, 2, 2, 0),
  M03 = c(2, 0, 2, 0, 2, 0, 2, 0),
  M04 = c(2, 2, 0, 0, 2, 0, 0, 2),
  M05 = c(0, 0, 2, 0, 2, 2, 2, 0),
  M06 = c(2, 0, 0, 2, 2, 0, 0, 2),
  M07 = c(0, 2, 2, 2, 0, 0, 2, 0),
  M08 = c(2, 2, 2, 0, 2, 0, 0, 2),
  check.names = FALSE
)

phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 61, 53, 67, 55, 64, 60, 62),
  protein = c(11.0, 10.8, 12.1, 10.5, 11.7, 10.9, 11.2, 10.7),
  disease = c(4.0, 3.1, 5.2, 2.6, 4.6, 2.9, 3.4, 3.0),
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:8)),
  Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
  Position_BP = c(0, 1, 2, 4, 0, 1, 3, 6) * 1e6,
  stringsAsFactors = FALSE
)

trait_direction <- data.frame(
  Trait = c("yield", "protein", "disease"),
  PhenotypeColumn = c("yield", "protein", "disease"),
  Selection_direction = c("increase", "increase", "decrease"),
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(out_dir, "phenotype.csv")
genotype_file <- file.path(out_dir, "genotype.csv")
map_file <- file.path(out_dir, "map.csv")
direction_file <- file.path(out_dir, "trait_direction.csv")

write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(genotype, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(trait_direction, direction_file, row.names = FALSE, quote = FALSE)

shortlist_size <- 4

method_varPMV <- "fast"
fast_result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,

  phenotype_id_col = "NAME",
  genotype_id_col = "NAME",
  direction_trait_col = "Trait",
  direction_column_col = "PhenotypeColumn",
  direction_direction_col = "Selection_direction",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,

  prediction_mode = "trait_by_trait",
  multi_trait_method = "auto",
  trait_value_metric = "pmv",
  method_varPMV = method_varPMV,

  duplicate_action = "none",
  n_crosses = 6,
  max_crosses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",

  output_dir = file.path(out_dir, "fast_outputs"),
  write_outputs = FALSE,
  write_figures = FALSE,
  seed = 20260623
)

fast_shortlist <- head(
  fast_result$candidate_crosses[order(-fast_result$candidate_crosses$multi_trait_score), ],
  shortlist_size
)

method_varPMV <- "full_posterior"
full_result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,

  phenotype_id_col = "NAME",
  genotype_id_col = "NAME",
  direction_trait_col = "Trait",
  direction_column_col = "PhenotypeColumn",
  direction_direction_col = "Selection_direction",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,

  prediction_mode = "trait_by_trait",
  multi_trait_method = "auto",
  trait_value_metric = "pmv",
  method_varPMV = method_varPMV,

  duplicate_action = "none",
  n_crosses = shortlist_size,
  max_crosses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",

  output_dir = file.path(out_dir, "full_posterior_outputs"),
  write_outputs = FALSE,
  write_figures = FALSE,
  seed = 20260623
)

stopifnot(inherits(fast_result, "ng_cross_prediction_result"))
stopifnot(inherits(full_result, "ng_cross_prediction_result"))
stopifnot(identical(fast_result$settings$method_varPMV, "fast"))
stopifnot(identical(full_result$settings$method_varPMV, "full_posterior"))
stopifnot(all(is.finite(full_result$candidate_crosses$yield_pmv_full_posterior)))
stopifnot(max(abs(full_result$candidate_crosses$yield_pmv_used -
                    full_result$candidate_crosses$yield_pmv_full_posterior), na.rm = TRUE) < 1e-10)

comparison <- merge(
  fast_shortlist[, c("parent1", "parent2", "yield_pmv_fast", "yield_pmv_used")],
  full_result$candidate_crosses[, c("parent1", "parent2", "yield_pmv_full_posterior", "yield_pmv_used")],
  by = c("parent1", "parent2"),
  suffixes = c("_fast_run", "_full_run")
)

comparison
head(full_result$selected_crosses)

message("Full-posterior PMV shortlist example completed.")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
