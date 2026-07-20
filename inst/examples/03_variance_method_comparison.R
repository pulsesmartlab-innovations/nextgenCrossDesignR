# Compare trait-value variance methods in one reproducible run.
#
# This example uses the same package workflow with several trait_value_metric
# choices. It is useful when a user wants to see how to switch between:
#   - var_complex
#   - uc with PMV
#   - uc with VPM
#   - pmv
#   - vpm
#   - le
#   - mean
#
# The script is self-contained. It creates small example input files and writes
# a summary table comparing the selected crosses from each method.

library(nextgenCrossDesign)

required_args <- c(
  "allocation_method",
  "alphamate_mode",
  "alphamate_target_degree",
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "n_iter",
  "burn_in",
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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_variance_method_comparison")
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

# # Balanced gain-diversity target mode
# alphamate_mode <- "ModeOptTarget1"
# 
# # Emphasize gain / criterion more strongly
# alphamate_mode <- "ModeMaxCriterion"
# 
# # Emphasize low coancestry / diversity more strongly
# alphamate_mode <- "ModeMinCoancestry"
# 
# So if the user wants to compare AlphaMate-style allocation behavior, a good example grid is:
#   
#   alphamate_modes <- c(
#     "ModeOptTarget1",
#     "ModeMaxCriterion",
#     "ModeMinCoancestry"
#   )
# 
# Important detail: alphamate_target_degree mainly matters for ModeOptTarget1. For ModeMaxCriterion and ModeMinCoancestry, the mode
# itself is already choosing the gain-heavy or diversity-heavy end of the frontier.

method_grid <- data.frame(
  run_name = c(
    "var_complex_alphamate_style",
    "uc_pmv_ocs",
    "uc_vpm_ocs",
    "pmv_ocs",
    "vpm_ocs",
    "var_simple_ocs",
    "mean_ocs"
  ),
  trait_value_metric = c(
    "var_complex",
    "usefulness",
    "usefulness",
    "pmv",
    "vpm",
    "le",
    "mean"
  ),
  uc_variance_source = c(
    "pmv",
    "pmv",
    "vpm",
    "pmv",
    "pmv",
    "pmv",
    "pmv"
  ),
  allocation_method = c(
    "alphamate_style",
    "ocs",
    "ocs",
    "ocs",
    "ocs",
    "ocs",
    "ocs"
  ),
  method_varPMV = c(
    "fast",
    "fast",
    "fast",
    "full_posterior",
    "fast",
    "fast",
    "fast"
  ),
  stringsAsFactors = FALSE
)

# To call the external AlphaMate executable, add a row with
# allocation_method = "alphamate_executable" and set alphamate_executable below.
alphamate_executable <- Sys.getenv("NG_ALPHAMATE_EXE", unset = "")
alphamate_runtime_path <- Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")

run_one_method <- function(row) {
  method_out <- file.path(out_dir, "outputs", row$run_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)

  ng_run_cross_prediction(
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

    trait_value_metric = row$trait_value_metric,
    uc_variance_source = row$uc_variance_source,
    selection_prop = 0.20,
    method_varPMV = row$method_varPMV,
    ril_mode = "infinite",
    run_posterior_prediction = FALSE,
    posterior_method = "mcmc",
    n_iter = 5000,
    burn_in = 500,
    use_parallel = FALSE,

    progeny = "DH",
    recomb_model = "haldane",
    assume_inbred = TRUE,

    duplicate_action = "none",
    n_crosses = 5,
    max_crosses_per_parent = 3,
    optimizer = "greedy_local",
    allocation_method = row$allocation_method,
    use_ocs = TRUE,
    lambda_group = 0.05,
    lambda_mating = 0.02,

    alphamate_mode = "ModeOptTarget1",
    alphamate_target_degree = 45,
    alphamate_max_contributions = 3,
    alphamate_executable = if (nzchar(alphamate_executable)) alphamate_executable else NULL,
    alphamate_runtime_path = alphamate_runtime_path,

    output_dir = method_out,
    write_outputs = FALSE,
    write_figures = TRUE,

    seed = 20260623
  )
}

results <- setNames(vector("list", nrow(method_grid)), method_grid$run_name)

for (i in seq_len(nrow(method_grid))) {
  row <- method_grid[i, , drop = FALSE]
  result <- run_one_method(row)

  stopifnot(inherits(result, "ng_cross_prediction_result"))
  stopifnot(identical(result$settings$trait_value_metric, row$trait_value_metric))
  stopifnot(identical(result$settings$allocation_method, row$allocation_method))
  stopifnot(nrow(result$selected_crosses) == 5L)
  stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
  stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))

  results[[row$run_name]] <- result
}

variance_summary <- do.call(rbind, lapply(names(results), function(run_name) {
  result <- results[[run_name]]
  top <- result$selected_crosses[1, , drop = FALSE]
  data.frame(
    run_name = run_name,
    trait_value_metric = result$settings$trait_value_metric,
    uc_variance_source = result$settings$uc_variance_source,
    method_varPMV = result$settings$method_varPMV,
    allocation_method = result$settings$allocation_method,
    selected_crosses = nrow(result$selected_crosses),
    top_parent1 = top$parent1,
    top_parent2 = top$parent2,
    top_multi_trait_score = top$multi_trait_score,
    top_pair_kinship = top$pair_kinship,
    output_dir = normalizePath(file.path(out_dir, "outputs", run_name), winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}))

summary_file <- file.path(out_dir, "variance_method_summary.csv")
write.csv(variance_summary, summary_file, row.names = FALSE)

variance_summary

message("Variance-method comparison completed.")
message("Summary file: ", normalizePath(summary_file, winslash = "/", mustWork = TRUE))
