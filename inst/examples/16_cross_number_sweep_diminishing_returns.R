# Cross-number sweep and diminishing-returns example.
#
# Use this when a breeder does not want to guess the number of crosses. The
# script scores candidate crosses once, then asks the package to optimize
# portfolios for a range of K values and identify a practical K criterion.

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_cross_number_sweep")
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

max_crosses_per_parent <- 4
max_pair_kinship <- Inf
lambda_group <- 0.05
lambda_mating <- 0.02
lambda_parent_use <- 0
lambda_parent_use_mode <- "absolute"
optimizer <- "greedy_local"
local_iter <- 1000
ocs_iter <- 3

K_range <- 2:6
criterion <- "elbow_relative"
relative_threshold <- 0.05
ne_min <- 4
coancestry_max <- 0.25

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
  trait_value_metric = "var_complex",
  method_varPMV = "fast",

  duplicate_action = "none",
  n_crosses = max(K_range),
  max_crosses_per_parent = max_crosses_per_parent,
  optimizer = optimizer,
  allocation_method = "ocs",
  use_ocs = TRUE,
  lambda_group = lambda_group,
  lambda_mating = lambda_mating,
  lambda_parent_use = lambda_parent_use,
  lambda_parent_use_mode = lambda_parent_use_mode,

  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = FALSE,
  seed = 20260623
)

parent_kinship <- ng_parent_kinship(result$cleaned_data$genotype)

curve <- ng_optimize_mating_plan_curve(
  scores = result$candidate_crosses,
  K_range = K_range,
  gain_col = "multi_trait_score",
  parent_kinship = parent_kinship,
  max_crosses_per_parent = max_crosses_per_parent,
  max_pair_kinship = max_pair_kinship,
  lambda_group = lambda_group,
  lambda_mating = lambda_mating,
  lambda_parent_use = lambda_parent_use,
  lambda_parent_use_mode = lambda_parent_use_mode,
  method = optimizer,
  local_iter = local_iter,
  ocs_iter = ocs_iter,
  criterion = criterion,
  relative_threshold = relative_threshold,
  ne_min = ne_min,
  coancestry_max = coancestry_max
)

recommended_K <- attr(curve, "elbow_K")

plot_file <- NA_character_
if (requireNamespace("ggplot2", quietly = TRUE)) {
  curve_plot <- ng_plot_diminishing_returns(curve)
  plot_file <- file.path(out_dir, "diminishing_returns.png")
  ggplot2::ggsave(plot_file, curve_plot, width = 7, height = 4.5, dpi = 150)
  stopifnot(file.exists(plot_file))
}

stopifnot(nrow(curve) == length(K_range))
stopifnot(identical(attr(curve, "criterion"), criterion))
stopifnot(length(recommended_K) == 1L)

curve
recommended_K
plot_file

message("Cross-number sweep example completed.")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
