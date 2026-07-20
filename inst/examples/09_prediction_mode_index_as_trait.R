# Prediction mode example 2: index_as_trait.
#
# Use this when the breeder already created a single selection-index column in
# the phenotype file and wants the package to predict crosses for that index.

library(nextgenCrossDesign)

required_version <- "0.3.12"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop("Install nextgenCrossDesign ", required_version, " or newer before running this example.", call. = FALSE)
}

set.seed(20260623)
out_dir <- file.path(tempdir(), "nextgenCrossDesign_prediction_index_as_trait")
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
phenotype$selection_index <- as.numeric(scale(phenotype$yield)) +
  0.5 * as.numeric(scale(phenotype$protein)) -
  0.8 * as.numeric(scale(phenotype$disease))
marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:8)),
  Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
  Position_BP = c(0, 1, 2, 4, 0, 1, 2, 4) * 1e6,
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(out_dir, "phenotype.csv")
genotype_file <- file.path(out_dir, "genotype.csv")
map_file <- file.path(out_dir, "map.csv")
write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(genotype, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)

# Parameters specific to prediction_mode = "index_as_trait".
prediction_mode <- "index_as_trait"
index_col <- "selection_index"
index_direction <- "increase"
direction_file <- NULL
multi_trait_method <- "auto"
trait_value_metric <- "var_complex"

parameter_notes <- data.frame(
  parameter = c("prediction_mode", "index_col", "index_direction", "direction_file", "multi_trait_method"),
  value = c(prediction_mode, index_col, index_direction, "NULL", multi_trait_method),
  when_to_change = c(
    "Use index_as_trait only when the phenotype file already contains one index column.",
    "Set to the phenotype column containing the user's precomputed index.",
    "Use increase when a larger index is better; use decrease when a smaller index is better.",
    "No direction file is required for index_as_trait.",
    "Usually auto because there is only one modeled response."
  ),
  stringsAsFactors = FALSE
)

result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = NULL,
  phenotype_id_col = "NAME",
  genotype_id_col = "NAME",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,
  prediction_mode = prediction_mode,
  index_col = index_col,
  index_direction = index_direction,
  multi_trait_method = multi_trait_method,
  trait_weights = NULL,
  trait_value_metric = trait_value_metric,
  uc_variance_source = "pmv",
  method_varPMV = "fast",
  ril_mode = "infinite",
  run_posterior_prediction = FALSE,
  posterior_method = "mcmc",
  n_iter = 5000,
  burn_in = 500,
  use_parallel = FALSE,
  progeny = "DH",
  recomb_model = "haldane",
  duplicate_action = "none",
  n_crosses = 4,
  max_crosses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",
  use_ocs = TRUE,
  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,
  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$prediction_mode, "index_as_trait"))
stopifnot(identical(result$trait_direction$column, index_col))
stopifnot("selection_index_value" %in% names(result$candidate_crosses))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))

parameter_notes
result$effect_summary
head(result$selected_crosses)

message("Prediction mode index_as_trait example completed.")
