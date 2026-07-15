# Run nextgenCrossDesign with package OCS allocation.
#
# This script is self-contained. It creates small phenotype, genotype, map, and
# trait-direction files, then runs ng_run_cross_prediction() with
# allocation_method = "ocs". Use this when you want package-native constrained
# mate allocation rather than AlphaMate-style or external AlphaMate allocation.

library(nextgenCrossDesign)

required_args <- c(
  "allocation_method",
  "lambda_parent_use",
  "lambda_parent_use_mode",
  "min_unique_parents",
  "max_pair_kinship",
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "nIter",
  "burnIn",
  "use_parallel"
)
missing_args <- setdiff(required_args, names(formals(nextgenCrossDesign::ng_run_cross_prediction)))
if (length(missing_args)) {
  stop(
    "This example requires nextgenCrossDesign 0.3.11 or newer. Reinstall the current tarball, restart R, ",
    "and confirm packageVersion('nextgenCrossDesign') >= '0.3.11'. Missing arguments in your active package: ",
    paste(missing_args, collapse = ", "),
    call. = FALSE
  )
}

set.seed(20260623)

out_dir <- file.path(tempdir(), "nextgenCrossDesign_ocs_user_run")
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
  M09 = c(0, 2, 0, 0, 2, 2, 0, 2),
  M10 = c(2, 0, 2, 2, 0, 0, 2, 0),
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
  SNP_code = paste0("M", sprintf("%02d", 1:10)),
  Chromosome = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
  Position_BP = c(0, 1, 2, 3, 4, 0, 1, 2, 3, 4) * 1e6,
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

# OCS settings.
# For a real breeding run, values such as n_crosses <- 100 and
# min_unique_parents <- 50 are reasonable starting points when the candidate
# parent set is large enough. This small example uses smaller feasible values.
allocation_method <- "ocs"

n_crosses <- 5
max_uses_per_parent <- 3
min_unique_parents <- 6
max_pair_kinship <- Inf

optimizer <- "greedy_local"
# Other optimizer choices:
#   "auto"
#   "repair_local"
#   "mip_linear"
#   "mip_contribution"

use_ocs <- TRUE
lambda_group <- 0.05
lambda_mating <- 0.02
lambda_parent_use <- 0
lambda_parent_use_mode <- "absolute"
# Other choice:
#   "adaptive"

local_iter <- 2000
ocs_iter <- 5

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
  multi_trait_method = "auto",
  trait_weights = NULL,

  trait_value_metric = "var_complex",
  uc_variance_source = "pmv",
  selection_prop = 0.20,
  method_varPMV = "fast",
  ril_mode = "infinite",
  run_posterior_prediction = FALSE,
  posterior_method = "mcmc",
  nIter = 5000,
  burnIn = 500,
  use_parallel = FALSE,

  progeny = "DH",
  recombination_model = "haldane",
  assume_inbred = TRUE,

  duplicate_action = "none",
  n_crosses = n_crosses,
  max_uses_per_parent = max_uses_per_parent,
  min_unique_parents = min_unique_parents,
  max_pair_kinship = max_pair_kinship,
  optimizer = optimizer,
  allocation_method = allocation_method,
  use_ocs = use_ocs,
  lambda_group = lambda_group,
  lambda_mating = lambda_mating,
  lambda_parent_use = lambda_parent_use,
  lambda_parent_use_mode = lambda_parent_use_mode,
  local_iter = local_iter,
  ocs_iter = ocs_iter,

  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,

  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$settings$allocation_method, "ocs"))
stopifnot(nrow(result$selected_crosses) == n_crosses)
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))

result$input_match_audit
result$effect_summary
head(result$selected_crosses)
result$output_files

message("OCS nextgenCrossDesign user run completed.")
message("Output directory: ", normalizePath(file.path(out_dir, "outputs"), winslash = "/", mustWork = TRUE))
