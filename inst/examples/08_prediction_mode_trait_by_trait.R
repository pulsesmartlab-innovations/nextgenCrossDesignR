# Prediction mode example 1: trait_by_trait.
#
# Use this when the phenotype file has separate trait columns and the
# direction file tells the package which traits to increase or decrease.

library(nextgenCrossDesign)

required_version <- "0.3.12"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop("Install nextgenCrossDesign ", required_version, " or newer before running this example.", call. = FALSE)
}

set.seed(20260623)
out_dir <- file.path(tempdir(), "nextgenCrossDesign_prediction_trait_by_trait")
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
  Position_BP = c(0, 1, 2, 4, 0, 1, 2, 4) * 1e6,
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

# Parameters specific to prediction_mode = "trait_by_trait".
prediction_mode <- "trait_by_trait"
traits_to_use <- c("yield", "protein", "disease")
multi_trait_method <- "auto"
trait_weights <- NULL
trait_value_metric <- "var_complex"
threshold_policy <- "soft"

parameter_notes <- data.frame(
  parameter = c("prediction_mode", "traits_to_use", "direction_file", "multi_trait_method", "trait_value_metric"),
  value = c(prediction_mode, paste(traits_to_use, collapse = ", "), basename(direction_file),
            multi_trait_method, trait_value_metric),
  when_to_change = c(
    "Keep trait_by_trait when the phenotype file has separate traits.",
    "Set to NULL to use every trait in the direction file, or list selected traits.",
    "Required for trait_by_trait so the package knows trait columns and directions.",
    "Use auto unless the breeder provides weights, economic weights, or desired gains.",
    "Use var_complex as the practical default for multi-trait cross value."
  ),
  stringsAsFactors = FALSE
)

result <- ng_run_cross_prediction(
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
  prediction_mode = prediction_mode,
  traits_to_use = traits_to_use,
  multi_trait_method = multi_trait_method,
  trait_weights = trait_weights,
  threshold_policy = threshold_policy,
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
stopifnot(identical(result$prediction_mode, "trait_by_trait"))
stopifnot(identical(result$settings$multi_trait_method, "auto"))
stopifnot(all(paste0(traits_to_use, "_value") %in% names(result$candidate_crosses)))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))

parameter_notes
result$effect_summary
head(result$selected_crosses)

message("Prediction mode trait_by_trait example completed.")
