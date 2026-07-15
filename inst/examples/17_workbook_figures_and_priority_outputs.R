# Workbook, figures, and priority-tier outputs.
#
# This script demonstrates user-facing outputs from ng_run_cross_prediction():
# selected crosses, priority tiers, the priority-vs-kinship figure, and the
# workbook when openxlsx is installed.

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_workbook_figures")
output_dir <- file.path(out_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ids <- paste0("P", sprintf("%02d", 1:10))

genotype <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 2, 2, 0, 2, 0, 2, 0, 2),
  M02 = c(0, 2, 0, 2, 0, 2, 2, 0, 2, 0),
  M03 = c(2, 0, 2, 0, 2, 0, 2, 0, 2, 0),
  M04 = c(2, 2, 0, 0, 2, 0, 0, 2, 2, 0),
  M05 = c(0, 0, 2, 0, 2, 2, 2, 0, 0, 2),
  M06 = c(2, 0, 0, 2, 2, 0, 0, 2, 2, 0),
  M07 = c(0, 2, 2, 2, 0, 0, 2, 0, 0, 2),
  M08 = c(2, 2, 2, 0, 2, 0, 0, 2, 0, 2),
  M09 = c(0, 2, 0, 0, 2, 2, 0, 2, 2, 0),
  M10 = c(2, 0, 2, 2, 0, 0, 2, 0, 0, 2),
  check.names = FALSE
)

phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 61, 53, 67, 55, 64, 60, 62, 57, 66),
  protein = c(11.0, 10.8, 12.1, 10.5, 11.7, 10.9, 11.2, 10.7, 11.5, 10.6),
  disease = c(4.0, 3.1, 5.2, 2.6, 4.6, 2.9, 3.4, 3.0, 4.4, 2.7),
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

write_outputs <- requireNamespace("openxlsx", quietly = TRUE)
write_figures <- TRUE

priority_breaks <- c(0.25, 0.50, 0.75, 1.00)
priority_labels <- c(
  "highly_priority",
  "priority",
  "medium_priority",
  "low_priority"
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
  multi_trait_method = "auto",
  trait_value_metric = "var_complex",
  method_varPMV = "fast",

  duplicate_action = "none",
  n_crosses = 8,
  max_uses_per_parent = 4,
  optimizer = "greedy_local",
  allocation_method = "ocs",
  use_ocs = TRUE,

  priority_breaks = priority_breaks,
  priority_labels = priority_labels,

  output_dir = output_dir,
  output_file = "crossing_plan.xlsx",
  write_outputs = write_outputs,
  write_figures = write_figures,
  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(nrow(result$selected_crosses) == 8L)
stopifnot("priority_tier" %in% names(result$selected_crosses))
stopifnot(all(as.character(result$selected_crosses$priority_tier) %in% priority_labels))
stopifnot(file.exists(result$output_files$priority_score_vs_kinship_png))

if (isTRUE(write_outputs)) {
  stopifnot(file.exists(result$output_files$workbook))
}

result$output_files
table(result$selected_crosses$priority_tier)
head(result$selected_crosses)

message("Workbook, figures, and priority output example completed.")
message("Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
