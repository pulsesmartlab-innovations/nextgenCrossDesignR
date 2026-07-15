# Multi-trait method example 4: desired_gain.
#
# Use desired_gain when the breeder has explicit response targets. Put the
# desired changes and economic weights in the direction file.

library(nextgenCrossDesign)

required_version <- "0.3.12"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop("Install nextgenCrossDesign ", required_version, " or newer before running this example.", call. = FALSE)
}

set.seed(20260623)
out_dir <- file.path(tempdir(), "nextgenCrossDesign_multitrait_desired_gain")
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
  desired_change = c(5.0, 0.4, 1.0),
  economic_weight = c(1.00, 0.45, 1.20),
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

# Parameters specific to multi_trait_method = "desired_gain".
multi_trait_method <- "desired_gain"
trait_weights <- NULL
threshold_policy <- "soft"

parameter_notes <- data.frame(
  parameter = c("multi_trait_method", "desired_change", "economic_weight", "trait_weights", "threshold_policy"),
  value = c(multi_trait_method,
            paste(trait_direction$Trait, trait_direction$desired_change, sep = "=", collapse = ", "),
            paste(trait_direction$Trait, trait_direction$economic_weight, sep = "=", collapse = ", "),
            "NULL", threshold_policy),
  when_to_change = c(
    "Use desired_gain when explicit response targets are available.",
    "Set in the direction file in the beneficial direction for each trait.",
    "Used with desired gains to choose index coefficients.",
    "Keep NULL when desired_change and economic_weight are supplied.",
    "Soft keeps near-miss crosses available but penalized."
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
  prediction_mode = "trait_by_trait",
  multi_trait_method = multi_trait_method,
  trait_weights = trait_weights,
  threshold_policy = threshold_policy,
  trait_value_metric = "var_complex",
  uc_variance_source = "pmv",
  method_varPMV = "fast",
  ril_mode = "infinite",
  run_posterior_prediction = FALSE,
  posterior_method = "mcmc",
  nIter = 5000,
  burnIn = 500,
  use_parallel = FALSE,
  progeny = "DH",
  recombination_model = "haldane",
  duplicate_action = "none",
  n_crosses = 4,
  max_uses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",
  use_ocs = TRUE,
  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,
  seed = 20260623
)

desired_gain_coefficients <- attr(result$candidate_crosses, "multi_trait")$desired_gain_coefficients
desired_gain_target <- attr(result$candidate_crosses, "multi_trait")$desired_gain_target

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$objective$method, "desired_gain"))
stopifnot(identical(result$objective$diagnostics$method_reason, "requested_desired_gain"))
stopifnot(!is.null(desired_gain_coefficients))
stopifnot(all(is.finite(desired_gain_coefficients)))
stopifnot(!is.null(desired_gain_target))

parameter_notes
desired_gain_coefficients
desired_gain_target
head(result$selected_crosses)

message("Multi-trait desired_gain example completed.")
