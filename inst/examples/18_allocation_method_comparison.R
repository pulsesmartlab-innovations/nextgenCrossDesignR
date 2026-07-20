# Allocation-method comparison.
#
# This script runs the same cross-prediction input with package OCS,
# native AlphaMate-style allocation, and external AlphaMate when an executable
# is available through NG_ALPHAMATE_EXE.

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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_allocation_comparison")
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

alphamate_executable <- Sys.getenv("NG_ALPHAMATE_EXE", unset = "")
alphamate_runtime_path <- Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")

allocation_grid <- data.frame(
  run_name = c(
    "package_ocs",
    "native_alphamate_style",
    "external_alphamate"
  ),
  allocation_method = c(
    "ocs",
    "alphamate_style",
    "alphamate_executable"
  ),
  optimizer = c(
    "greedy_local",
    "greedy_local",
    "greedy_local"
  ),
  alphamate_mode = c(
    "ModeOptTarget1",
    "ModeMaxCriterion",
    "ModeOptTarget1"
  ),
  alphamate_target_degree = c(45, 45, 45),
  stringsAsFactors = FALSE
)

run_allocation <- function(row) {
  method_out <- file.path(out_dir, "outputs", row$run_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)

  if (identical(row$allocation_method, "alphamate_executable") &&
      !file.exists(alphamate_executable)) {
    return(data.frame(
      run_name = row$run_name,
      allocation_method = row$allocation_method,
      status = "skipped",
      reason = "Set NG_ALPHAMATE_EXE to run the external AlphaMate example.",
      selected_crosses = NA_integer_,
      top_parent1 = NA_character_,
      top_parent2 = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

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
    n_crosses = 5,
    max_crosses_per_parent = 3,
    optimizer = row$optimizer,
    allocation_method = row$allocation_method,
    use_ocs = TRUE,
    lambda_group = 0.05,
    lambda_mating = 0.02,

    alphamate_mode = row$alphamate_mode,
    alphamate_target_degree = row$alphamate_target_degree,
    alphamate_max_contributions = 3,
    alphamate_executable = if (nzchar(alphamate_executable)) alphamate_executable else NULL,
    alphamate_runtime_path = alphamate_runtime_path,

    output_dir = method_out,
    write_outputs = FALSE,
    write_figures = TRUE,
    seed = 20260623
  )

  stopifnot(inherits(result, "ng_cross_prediction_result"))
  stopifnot(identical(result$settings$allocation_method, row$allocation_method))
  stopifnot(nrow(result$selected_crosses) == 5L)

  top <- result$selected_crosses[1, , drop = FALSE]
  data.frame(
    run_name = row$run_name,
    allocation_method = row$allocation_method,
    status = "completed",
    reason = NA_character_,
    selected_crosses = nrow(result$selected_crosses),
    top_parent1 = top$parent1,
    top_parent2 = top$parent2,
    stringsAsFactors = FALSE
  )
}

allocation_results <- lapply(seq_len(nrow(allocation_grid)), function(i) {
  run_allocation(allocation_grid[i, , drop = FALSE])
})

allocation_summary <- do.call(rbind, allocation_results)
stopifnot(all(allocation_summary$status == "completed" | allocation_summary$status == "skipped"))
stopifnot(allocation_summary$status[allocation_summary$allocation_method == "ocs"] == "completed")
stopifnot(allocation_summary$status[allocation_summary$allocation_method == "alphamate_style"] == "completed")

summary_file <- file.path(out_dir, "allocation_method_summary.csv")
write.csv(allocation_summary, summary_file, row.names = FALSE)

allocation_summary

message("Allocation-method comparison example completed.")
message("Summary file: ", normalizePath(summary_file, winslash = "/", mustWork = TRUE))
