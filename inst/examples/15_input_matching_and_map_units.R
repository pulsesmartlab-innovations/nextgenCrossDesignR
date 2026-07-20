# Input matching and map-unit example.
#
# This script uses nonstandard column names and deliberately shuffled
# phenotype and marker-map rows. The package matches parents by ID, aligns
# marker-map rows to genotype marker columns, and converts base-pair positions
# to the internal cM scale using bp_per_cm.

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_input_matching_example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- paste0("Line", sprintf("%02d", 1:8))
marker_ids <- paste0("SNP", sprintf("%02d", 1:10))

genotype <- data.frame(
  GenoID = ids,
  SNP01 = c(0, 0, 2, 2, 0, 2, 0, 2),
  SNP02 = c(0, 2, 0, 2, 0, 2, 2, 0),
  SNP03 = c(2, 0, 2, 0, 2, 0, 2, 0),
  SNP04 = c(2, 2, 0, 0, 2, 0, 0, 2),
  SNP05 = c(0, 0, 2, 0, 2, 2, 2, 0),
  SNP06 = c(2, 0, 0, 2, 2, 0, 0, 2),
  SNP07 = c(0, 2, 2, 2, 0, 0, 2, 0),
  SNP08 = c(2, 2, 2, 0, 2, 0, 0, 2),
  SNP09 = c(0, 2, 0, 0, 2, 2, 0, 2),
  SNP10 = c(2, 0, 2, 2, 0, 0, 2, 0),
  check.names = FALSE
)

phenotype_base <- data.frame(
  LineID = ids,
  Yield_kg = c(58, 61, 53, 67, 55, 64, 60, 62),
  Protein_pct = c(11.0, 10.8, 12.1, 10.5, 11.7, 10.9, 11.2, 10.7),
  Disease_score = c(4.0, 3.1, 5.2, 2.6, 4.6, 2.9, 3.4, 3.0),
  stringsAsFactors = FALSE
)
phenotype <- phenotype_base[c(4, 2, 8, 1, 6, 3, 5, 7), , drop = FALSE]

marker_map_base <- data.frame(
  MarkerName = marker_ids,
  ChrLabel = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
  PhysicalPositionBP = c(0, 1, 2, 3, 4, 0, 1, 2, 3, 4) * 1e6,
  stringsAsFactors = FALSE
)
marker_map <- marker_map_base[c(10, 1, 5, 8, 3, 6, 2, 9, 4, 7), , drop = FALSE]

trait_direction <- data.frame(
  TraitName = c("grain_yield", "protein", "disease"),
  PhenotypeColumn = c("Yield_kg", "Protein_pct", "Disease_score"),
  SelectionGoal = c("increase", "increase", "decrease"),
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(out_dir, "phenotype_shuffled.csv")
genotype_file <- file.path(out_dir, "genotype.csv")
map_file <- file.path(out_dir, "map_shuffled_bp.csv")
direction_file <- file.path(out_dir, "trait_direction_custom_names.csv")

write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(genotype, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(trait_direction, direction_file, row.names = FALSE, quote = FALSE)

phenotype_id_col <- "LineID"
genotype_id_col <- "GenoID"
direction_trait_col <- "TraitName"
direction_column_col <- "PhenotypeColumn"
direction_direction_col <- "SelectionGoal"
map_marker_col <- "MarkerName"
map_chr_col <- "ChrLabel"
map_pos_bp_col <- "PhysicalPositionBP"
map_position_unit <- "bp"
bp_per_cm <- 1e6

result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,

  phenotype_id_col = phenotype_id_col,
  genotype_id_col = genotype_id_col,
  direction_trait_col = direction_trait_col,
  direction_column_col = direction_column_col,
  direction_direction_col = direction_direction_col,
  map_marker_col = map_marker_col,
  map_chr_col = map_chr_col,
  map_pos_bp_col = map_pos_bp_col,
  map_position_unit = map_position_unit,
  bp_per_cm = bp_per_cm,

  prediction_mode = "trait_by_trait",
  traits_to_use = c("grain_yield", "protein", "disease"),
  multi_trait_method = "auto",
  trait_value_metric = "var_complex",
  method_varPMV = "fast",

  duplicate_action = "none",
  n_crosses = 5,
  max_crosses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",

  output_dir = file.path(out_dir, "outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,
  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$input_match_audit$phenotype_id_col, phenotype_id_col))
stopifnot(identical(result$input_match_audit$genotype_id_col, genotype_id_col))
stopifnot(identical(result$input_match_audit$map_marker_col, map_marker_col))
stopifnot(identical(result$input_match_audit$map_position_unit, map_position_unit))
stopifnot(identical(result$input_match_audit$bp_per_cm, bp_per_cm))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))
stopifnot(identical(result$cleaned_data$marker_map$pos_bp, marker_map_base$PhysicalPositionBP))

result$input_match_audit
head(result$selected_crosses)

message("Input matching and map-unit example completed.")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
