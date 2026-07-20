# QC duplicate removal and reporting.
#
# This script shows duplicate-genotype QC as part of the package workflow.
# It first reports and plots a deliberate duplicate, then runs the full
# cross-prediction workflow with duplicate_action = "remove".

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_qc_duplicate_example")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- c(paste0("P", sprintf("%02d", 1:8)), "P02_dup")

base_genotype <- data.frame(
  NAME = paste0("P", sprintf("%02d", 1:8)),
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
genotype <- rbind(
  base_genotype,
  transform(base_genotype[base_genotype$NAME == "P02", ], NAME = "P02_dup")
)
rownames(genotype) <- NULL

phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 61, 53, 67, 55, 64, 60, 62, 61),
  protein = c(11.0, 10.8, 12.1, 10.5, 11.7, 10.9, 11.2, 10.7, 10.8),
  disease = c(4.0, 3.1, 5.2, 2.6, 4.6, 2.9, 3.4, 3.0, 3.1),
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

phenotype_file <- file.path(out_dir, "phenotype_with_duplicate.csv")
genotype_file <- file.path(out_dir, "genotype_with_duplicate.csv")
map_file <- file.path(out_dir, "map.csv")
direction_file <- file.path(out_dir, "trait_direction.csv")

write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(genotype, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(trait_direction, direction_file, row.names = FALSE, quote = FALSE)

duplicate_action <- "remove"
duplicate_action_report <- "report"
duplicate_threshold <- 0.995
duplicate_maf_min <- 0
duplicate_max_missing_prop <- 0.40
duplicate_min_compared_markers <- 6

qc_report <- ng_preflight_input_tables(
  geno = genotype,
  phenotype = phenotype,
  trait_spec = trait_direction,
  marker_map = marker_map,
  ploidy = 2L,
  putative_duplicate_check = TRUE,
  putative_duplicate_action = duplicate_action_report,
  duplicate_threshold = duplicate_threshold,
  duplicate_maf_min = duplicate_maf_min,
  duplicate_max_missing_prop = duplicate_max_missing_prop,
  duplicate_min_compared_markers = duplicate_min_compared_markers,
  putative_duplicate_return_similarity = TRUE
)

stopifnot("putative_duplicate_genotypes" %in% qc_report$issues$id)
stopifnot(nrow(qc_report$putative_duplicates$pairs) >= 1L)

duplicate_plot <- file.path(out_dir, "putative_duplicate_qc.png")
ng_plot_putative_duplicates(
  qc_report$putative_duplicates,
  output_path = duplicate_plot,
  title = "Putative duplicate check before cross prediction"
)
stopifnot(file.exists(duplicate_plot))

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

  duplicate_action = duplicate_action,
  duplicate_threshold = duplicate_threshold,
  duplicate_maf_min = duplicate_maf_min,
  duplicate_max_missing_prop = duplicate_max_missing_prop,
  duplicate_min_compared_markers = duplicate_min_compared_markers,

  n_crosses = 5,
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
stopifnot("putative_duplicate_genotypes_removed" %in% result$qc$issues$id)
stopifnot(!("P02_dup" %in% rownames(result$cleaned_data$genotype)))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
stopifnot(file.exists(result$output_files$priority_score_vs_kinship_png))

result$qc$issues
result$qc$cleaning$putative_duplicates$removed_parents
head(result$selected_crosses)

message("QC duplicate removal and reporting example completed.")
message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = TRUE))
