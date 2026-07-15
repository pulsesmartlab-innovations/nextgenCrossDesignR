# Compare OCS optimizer choices and their parameters.
#
# This script is self-contained. It creates small phenotype, genotype, map, and
# trait-direction files, then runs ng_run_cross_prediction() with
# allocation_method = "ocs" while changing the optimizer row by row.
#
# WHICH optimizer SHOULD I USE? (see docs/BACKEND_USER_GUIDE.md "Choosing An Optimizer")
# The optimizer only SOLVES the allocation problem -- it does not change the metric or the
# objective, only how well/fast the mating plan is found. It does not affect prediction accuracy.
#   * optimizer = "evolution" -- RECOMMENDED default for real breeding programs. The native
#     memetic genetic algorithm: warm-started from greedy (never worse), competitive with the
#     exact MIP, and it scales to large candidate sets where MIP is size/time-guarded. Tune with
#     evol_solutions / evol_iterations / evol_stop. See example 26.
#   * "auto" -- routes to the exact MIP when the problem is small/tractable, else greedy_local.
#     Good when you want the exact single-decision objective and the candidate set is small.
#   * "mip_contribution" (alias "mip") -- exact-style OCS objective for ONE mating decision;
#     slowest, size/time-guarded. "mip_linear" is the no-coancestry-penalty variant (lambda_group=0).
#   * "greedy_local" -- fastest / deterministic diagnostic; behind on realized gain.
# Evidence (config-scoped, directional):
#   - Single-shot objective bake-off (tools/run_optimizer_benchmark.R): mip_contribution ranks
#     first on ACHIEVED OBJECTIVE, evolution a close second, both beating greedy/repair.
#   - Recurrent REALIZED gain (tools/run_method_optimizer_metric_study.R, 10 reps x 25 cycles,
#     metric held at uc/pmv): evolution ~= MIP (cycle-25 gain 7.20 +/- 0.30 vs 7.07 +/- 0.44,
#     paired diff not significant), both clearly ahead of greedy (6.78 +/- 0.30).
# Bottom line: use "evolution" for real/recurrent/large programs; "auto" or "mip" when you want
# the exact objective for a single small mating decision.

library(nextgenCrossDesign)

required_args <- c(
  "allocation_method",
  "lambda_parent_use",
  "lambda_parent_use_mode",
  "local_iter",
  "ocs_iter",
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

out_dir <- file.path(tempdir(), "nextgenCrossDesign_optimizer_parameter_guide")
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

# Shared OCS capacity settings.
allocation_method <- "ocs"
n_crosses <- 5
max_uses_per_parent <- 4
min_unique_parents <- 4
max_pair_kinship <- Inf
use_ocs <- TRUE

# Optimizer-specific guide.
# - auto chooses the best available optimizer from the settings and installed packages.
# - greedy_local and repair_local use local_iter for local search.
# - mip_linear uses lpSolve to solve a linear parent-capacity problem.
# - mip_contribution uses lpSolve plus lambda_parent_use, lambda_group, and ocs_iter.
# lpSolve is a package import, so these optimizers should be available after
# installing nextgenCrossDesign with its dependencies.
if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop(
    "lpSolve is required for mip_linear and mip_contribution. Reinstall ",
    "nextgenCrossDesign with dependencies before running this guide.",
    call. = FALSE
  )
}

optimizer_grid <- data.frame(
  run_name = c(
    "auto_default",
    "greedy_local_search",
    "repair_local_search",
    "mip_linear_exact",
    "mip_contribution_ocs"
  ),
  optimizer = c(
    "auto",
    "greedy_local",
    "repair_local",
    "mip_linear",
    "mip_contribution"
  ),
  lambda_group = c(
    0.05,
    0.05,
    0.05,
    0.00,
    0.05
  ),
  lambda_mating = c(
    0.02,
    0.02,
    0.02,
    0.02,
    0.02
  ),
  lambda_parent_use = c(
    0.25,
    0.00,
    0.00,
    0.00,
    0.25
  ),
  lambda_parent_use_mode = c(
    "adaptive",
    "absolute",
    "absolute",
    "absolute",
    "adaptive"
  ),
  local_iter = c(
    2000,
    2000,
    4000,
    2000,
    2000
  ),
  ocs_iter = c(
    5,
    5,
    5,
    1,
    5
  ),
  note = c(
    "auto dispatch; uses mip_contribution when MIP penalties are active",
    "fast local search; local_iter controls swap attempts",
    "local repair search; local_iter controls repair attempts",
    "exact linear MIP for score and parent caps; lambda_group and lambda_parent_use do not change the solve",
    "MIP contribution optimizer; uses parent-use penalty, group penalty, and ocs_iter"
  ),
  stringsAsFactors = FALSE
)

run_optimizer <- function(row) {
  method_out <- file.path(out_dir, "outputs", row$run_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)

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
    optimizer = row$optimizer,
    allocation_method = allocation_method,
    use_ocs = use_ocs,
    lambda_group = row$lambda_group,
    lambda_mating = row$lambda_mating,
    lambda_parent_use = row$lambda_parent_use,
    lambda_parent_use_mode = row$lambda_parent_use_mode,
    local_iter = row$local_iter,
    ocs_iter = row$ocs_iter,

    output_dir = method_out,
    write_outputs = FALSE,
    write_figures = TRUE,

    seed = 20260623
  )

  stopifnot(inherits(result, "ng_cross_prediction_result"))
  stopifnot(identical(result$settings$allocation_method, "ocs"))
  stopifnot(identical(result$settings$optimizer, row$optimizer))
  stopifnot(nrow(result$selected_crosses) == n_crosses)
  stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
  stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))

  top <- result$selected_crosses[1, , drop = FALSE]
  data.frame(
    run_name = row$run_name,
    optimizer = row$optimizer,
    status = "completed",
    optimizer_method = result$settings$optimizer_method,
    selected_crosses = nrow(result$selected_crosses),
    top_parent1 = top$parent1,
    top_parent2 = top$parent2,
    top_multi_trait_score = top$multi_trait_score,
    output_dir = normalizePath(method_out, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

optimizer_results <- lapply(seq_len(nrow(optimizer_grid)), function(i) {
  run_optimizer(optimizer_grid[i, , drop = FALSE])
})

optimizer_summary <- do.call(rbind, optimizer_results)
stopifnot(all(optimizer_summary$status == "completed"))
summary_file <- file.path(out_dir, "optimizer_parameter_summary.csv")
write.csv(optimizer_summary, summary_file, row.names = FALSE)

optimizer_grid
optimizer_summary

message("Optimizer parameter guide completed.")
message("Summary file: ", normalizePath(summary_file, winslash = "/", mustWork = TRUE))
