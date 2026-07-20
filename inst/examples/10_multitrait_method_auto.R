# Multi-trait method example 1: auto.
#
# Use auto when the user only supplies trait directions. The package gives
# each oriented trait equal normalized weight and records that choice in
# result$objective$diagnostics.
#
# WHICH multi_trait_method SHOULD I USE? (see docs/BACKEND_USER_GUIDE.md "Choosing A
# Multi-Trait Method"). The methods optimize DIFFERENT objectives, so the choice is driven by
# your inputs, not a performance race:
#   * only trait directions (+ rough thresholds)     -> "auto"  (RECOMMENDED default; this file)
#   * declared relative importance weights           -> "weighted"        (example 11)
#   * reliable economic weights + genetic cov matrix -> "economic_index"  (example 12)
#   * target genetic changes per trait               -> "desired_gain"    (example 13)
#   * hard/soft min/max constraints per trait        -> "threshold" (combinable with the above)
# economic_index is NOT the default: reliable economic weights are hard to get and users often
# substitute phenotypic for genetic correlations (violating Smith-Hazel) -- when in doubt, auto.

library(nextgenCrossDesign)

required_version <- "0.3.12"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop("Install nextgenCrossDesign ", required_version, " or newer before running this example.", call. = FALSE)
}

set.seed(20260623)
out_dir <- file.path(tempdir(), "nextgenCrossDesign_multitrait_auto")
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

# Parameters specific to multi_trait_method = "auto".
multi_trait_method <- "auto"
trait_weights <- NULL
threshold_policy <- "soft"
threshold_penalty_weight <- 1.0
threshold_penalty_autoscale <- TRUE

parameter_notes <- data.frame(
  parameter = c("multi_trait_method", "trait_weights", "threshold_policy", "threshold_penalty_weight"),
  value = c(multi_trait_method, "NULL", threshold_policy, threshold_penalty_weight),
  when_to_change = c(
    "Use auto when only trait directions are supplied.",
    "Keep NULL when the breeder has no defensible trait weights.",
    "Soft keeps near-miss crosses available but penalized.",
    "Increase when threshold misses should matter more in ranking."
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
  threshold_penalty_weight = threshold_penalty_weight,
  threshold_penalty_autoscale = threshold_penalty_autoscale,
  trait_value_metric = "var_complex",
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
stopifnot(identical(result$objective$method, "auto"))
stopifnot(identical(result$objective$diagnostics$method_reason, "equal_weight_rank_default"))
stopifnot(nrow(result$selected_crosses) == 4L)

parameter_notes
result$objective$diagnostics
head(result$selected_crosses)

message("Multi-trait auto example completed.")
