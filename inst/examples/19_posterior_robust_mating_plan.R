# Posterior robust mating-plan example.
#
# This script turns on posterior cross prediction, then uses
# ng_optimize_robust_mating_plan() to select crosses using a pessimistic
# posterior usefulness quantile rather than a point estimate only.

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_posterior_robust_plan")
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

run_posterior_prediction <- TRUE
posterior_method <- "closed_form"
n_iter <- 12
burn_in <- 4

robustness_quantile <- 0.25
objective <- "posterior_quantile"
optimizer <- "greedy_local"

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
  trait_value_metric = "pmv",
  method_varPMV = "fast",
  run_posterior_prediction = run_posterior_prediction,
  posterior_method = posterior_method,
  n_iter = n_iter,
  burn_in = burn_in,

  duplicate_action = "none",
  n_crosses = 4,
  max_crosses_per_parent = 3,
  optimizer = optimizer,
  allocation_method = "ocs",
  use_ocs = TRUE,

  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = FALSE,
  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(length(result$posterior_predictions) > 0L)

posterior_scores <- result$posterior_predictions[[1L]]
stopifnot(all(c("usefulness_pmv_gebv_post_mean", "usefulness_pmv_gebv_post_lower", "usefulness_pmv_gebv_post_upper") %in%
                names(posterior_scores)))

parent_kinship <- ng_parent_kinship(result$cleaned_data$genotype)

robust_plan <- ng_optimize_robust_mating_plan(
  posterior_scores = posterior_scores,
  n_crosses = 4,
  parent_kinship = parent_kinship,
  gain_col = "usefulness_pmv_gebv",
  robustness_quantile = robustness_quantile,
  objective = objective,
  max_crosses_per_parent = 3,
  lambda_group = 0.05,
  lambda_mating = 0.02,
  method = optimizer
)

stopifnot(nrow(robust_plan) == 4L)
stopifnot(identical(attr(robust_plan, "summary")$robust_objective, objective))

attr(robust_plan, "summary")
robust_plan
head(posterior_scores)

message("Posterior robust mating-plan example completed.")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
