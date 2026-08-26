fail <- function(...) stop(sprintf(...), call. = FALSE)

description_path <- "DESCRIPTION"
vignette_path <- file.path("vignettes", "nextgenCrossDesign.Rmd")
example_path <- file.path("inst", "examples", "01_chronological_cross_prediction_workflow.R")
minimal_example_path <- file.path("inst", "examples", "02_minimal_reproducible_user_run.R")
variance_example_path <- file.path("inst", "examples", "03_variance_method_comparison.R")
external_alpha_example_path <- file.path("inst", "examples", "04_alphamate_executable_user_run.R")
ocs_example_path <- file.path("inst", "examples", "05_ocs_user_run.R")
optimizer_example_path <- file.path("inst", "examples", "06_optimizer_parameter_guide.R")
trait_value_parameter_example_path <- file.path("inst", "examples", "07_trait_value_metric_parameter_guide.R")
prediction_trait_example_path <- file.path("inst", "examples", "08_prediction_mode_trait_by_trait.R")
prediction_index_example_path <- file.path("inst", "examples", "09_prediction_mode_index_as_trait.R")
multitrait_auto_example_path <- file.path("inst", "examples", "10_multitrait_method_auto.R")
multitrait_weighted_example_path <- file.path("inst", "examples", "11_multitrait_method_weighted.R")
multitrait_economic_example_path <- file.path("inst", "examples", "12_multitrait_method_economic_index.R")
multitrait_desired_example_path <- file.path("inst", "examples", "13_multitrait_method_desired_gain.R")
qc_duplicate_example_path <- file.path("inst", "examples", "14_qc_duplicate_removal_and_reporting.R")
input_matching_example_path <- file.path("inst", "examples", "15_input_matching_and_map_units.R")
cross_number_example_path <- file.path("inst", "examples", "16_cross_number_sweep_diminishing_returns.R")
workbook_figures_example_path <- file.path("inst", "examples", "17_workbook_figures_and_priority_outputs.R")
allocation_comparison_example_path <- file.path("inst", "examples", "18_allocation_method_comparison.R")
posterior_robust_example_path <- file.path("inst", "examples", "19_posterior_robust_mating_plan.R")
full_posterior_pmv_example_path <- file.path("inst", "examples", "20_full_posterior_pmv_shortlist.R")
examples_readme_path <- file.path("inst", "examples", "README.md")

if (!file.exists(description_path)) fail("Missing DESCRIPTION")
if (!file.exists(vignette_path)) fail("Missing vignette: %s", vignette_path)
if (!file.exists(example_path)) fail("Missing chronological example: %s", example_path)
if (!file.exists(minimal_example_path)) fail("Missing minimal reproducible example: %s", minimal_example_path)
if (!file.exists(variance_example_path)) fail("Missing variance-method comparison example: %s", variance_example_path)
if (!file.exists(external_alpha_example_path)) fail("Missing AlphaMate executable example: %s", external_alpha_example_path)
if (!file.exists(ocs_example_path)) fail("Missing OCS example: %s", ocs_example_path)
if (!file.exists(optimizer_example_path)) fail("Missing optimizer parameter guide example: %s", optimizer_example_path)
if (!file.exists(trait_value_parameter_example_path)) fail("Missing trait-value parameter guide example: %s", trait_value_parameter_example_path)
if (!file.exists(prediction_trait_example_path)) fail("Missing trait-by-trait prediction-mode example: %s", prediction_trait_example_path)
if (!file.exists(prediction_index_example_path)) fail("Missing index-as-trait prediction-mode example: %s", prediction_index_example_path)
if (!file.exists(multitrait_auto_example_path)) fail("Missing auto multi-trait method example: %s", multitrait_auto_example_path)
if (!file.exists(multitrait_weighted_example_path)) fail("Missing weighted multi-trait method example: %s", multitrait_weighted_example_path)
if (!file.exists(multitrait_economic_example_path)) fail("Missing economic-index multi-trait method example: %s", multitrait_economic_example_path)
if (!file.exists(multitrait_desired_example_path)) fail("Missing desired-gain multi-trait method example: %s", multitrait_desired_example_path)
if (!file.exists(qc_duplicate_example_path)) fail("Missing QC duplicate example: %s", qc_duplicate_example_path)
if (!file.exists(input_matching_example_path)) fail("Missing input-matching example: %s", input_matching_example_path)
if (!file.exists(cross_number_example_path)) fail("Missing cross-number sweep example: %s", cross_number_example_path)
if (!file.exists(workbook_figures_example_path)) fail("Missing workbook/figures example: %s", workbook_figures_example_path)
if (!file.exists(allocation_comparison_example_path)) fail("Missing allocation comparison example: %s", allocation_comparison_example_path)
if (!file.exists(posterior_robust_example_path)) fail("Missing posterior robust mating example: %s", posterior_robust_example_path)
if (!file.exists(full_posterior_pmv_example_path)) fail("Missing full-posterior PMV shortlist example: %s", full_posterior_pmv_example_path)
# Examples for the mate-selection capability additions (gain-diversity balance,
# progeny-inbreeding management, breeder constraints, marker/lethal management, cost, and
# the evolution optimizer).
capability_example_paths <- file.path("inst", "examples", c(
  "21_gain_diversity_balance_and_target_coancestry.R",
  "22_progeny_inbreeding_management.R",
  "23_breeder_mating_constraints.R",
  "24_marker_steering_and_lethal_guarding.R",
  "25_cost_and_logistics_factors.R",
  "26_evolution_optimizer.R"))
for (capability_example in capability_example_paths) {
  if (!file.exists(capability_example)) fail("Missing capability example: %s", capability_example)
}
if (!file.exists(examples_readme_path)) fail("Missing examples README: %s", examples_readme_path)
if (file.exists(file.path("inst", "examples", ".Rhistory"))) fail("inst/examples contains .Rhistory")

vignette_text <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")
example_text <- paste(readLines(example_path, warn = FALSE), collapse = "\n")
minimal_example_text <- paste(readLines(minimal_example_path, warn = FALSE), collapse = "\n")
variance_example_text <- paste(readLines(variance_example_path, warn = FALSE), collapse = "\n")
external_alpha_example_text <- paste(readLines(external_alpha_example_path, warn = FALSE), collapse = "\n")
ocs_example_text <- paste(readLines(ocs_example_path, warn = FALSE), collapse = "\n")
optimizer_example_text <- paste(readLines(optimizer_example_path, warn = FALSE), collapse = "\n")
trait_value_parameter_example_text <- paste(readLines(trait_value_parameter_example_path, warn = FALSE), collapse = "\n")
prediction_trait_example_text <- paste(readLines(prediction_trait_example_path, warn = FALSE), collapse = "\n")
prediction_index_example_text <- paste(readLines(prediction_index_example_path, warn = FALSE), collapse = "\n")
multitrait_auto_example_text <- paste(readLines(multitrait_auto_example_path, warn = FALSE), collapse = "\n")
multitrait_weighted_example_text <- paste(readLines(multitrait_weighted_example_path, warn = FALSE), collapse = "\n")
multitrait_economic_example_text <- paste(readLines(multitrait_economic_example_path, warn = FALSE), collapse = "\n")
multitrait_desired_example_text <- paste(readLines(multitrait_desired_example_path, warn = FALSE), collapse = "\n")
qc_duplicate_example_text <- paste(readLines(qc_duplicate_example_path, warn = FALSE), collapse = "\n")
input_matching_example_text <- paste(readLines(input_matching_example_path, warn = FALSE), collapse = "\n")
cross_number_example_text <- paste(readLines(cross_number_example_path, warn = FALSE), collapse = "\n")
workbook_figures_example_text <- paste(readLines(workbook_figures_example_path, warn = FALSE), collapse = "\n")
allocation_comparison_example_text <- paste(readLines(allocation_comparison_example_path, warn = FALSE), collapse = "\n")
posterior_robust_example_text <- paste(readLines(posterior_robust_example_path, warn = FALSE), collapse = "\n")
full_posterior_pmv_example_text <- paste(readLines(full_posterior_pmv_example_path, warn = FALSE), collapse = "\n")
examples_readme_text <- paste(readLines(examples_readme_path, warn = FALSE), collapse = "\n")
runner_text <- paste(readLines(file.path("R", "39_cross_prediction_runner.R"), warn = FALSE), collapse = "\n")

description <- read.dcf(description_path)[1, ]
imports <- if ("Imports" %in% names(description)) {
  trimws(strsplit(description[["Imports"]], ",", fixed = TRUE)[[1]])
} else {
  character(0)
}
suggests <- if ("Suggests" %in% names(description)) {
  trimws(strsplit(description[["Suggests"]], ",", fixed = TRUE)[[1]])
} else {
  character(0)
}
if (!("lpSolve" %in% imports)) fail("lpSolve must be listed in DESCRIPTION Imports")
if ("lpSolve" %in% suggests) fail("lpSolve must not be left in DESCRIPTION Suggests")

required_runner_args <- c(
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "n_iter",
  "burn_in",
  "use_parallel"
)
has_runner_arg <- vapply(required_runner_args, function(arg) {
  grepl(paste0(arg, " ="), runner_text, fixed = TRUE)
}, logical(1L))
missing_runner_args <- required_runner_args[!has_runner_arg]
if (length(missing_runner_args)) {
  fail("ng_run_cross_prediction missing required exposed argument(s): %s",
       paste(missing_runner_args, collapse = ", "))
}

# The vignette must still walk a reader from installation to exported outputs in order.
# These are the workflow milestones of the topic-structured tutorial; the ordering check
# below is what actually guards the chronology, so this list only needs the load-bearing
# stops, not every heading (leaving room to add topic sections without touching the test).
required_sections <- c(
  "## Installation and loading",
  "## Your input data",
  "## Quick start: from data to a crossing plan",
  "## Reading the result",
  "## Data QC, duplicates, and LD pruning",
  "## The cross-scoring metric",
  "## The optimizer",
  "## The allocation method and core OCS constraints",
  "## Building a multi-trait objective",
  "## Priority ranking and tiers",
  "## Workbooks, figures, and exports",
  "## Reproducible full script",
  "## Public API map"
)

positions <- vapply(required_sections, function(section) {
  pos <- regexpr(section, vignette_text, fixed = TRUE)[[1]]
  if (pos < 0) fail("Vignette missing chronological section: %s", section)
  pos
}, integer(1))

if (is.unsorted(positions, strictly = TRUE)) {
  fail("Chronological workflow sections are not in the expected order")
}

required_example_tokens <- c(
  "variance_method <-",
  "variance_method <- \"var_complex\"",
  "allocation_method <-",
  "alphamate_mode <-",
  "alphamate_target_degree <-",
  "optimizer <-",
  "solver <-",
  "use_ocs <-",
  "max_crosses <-",
  "max_crosses_per_parent <-",
  "min_variance <-",
  "ng_preflight_input_tables(",
  "ng_score_crosses(",
  "ng_breeder_selection_objective(",
  "ng_optimize_breeder_selection_plan(",
  "ng_alphamate_style_select(",
  "ng_select_alphamate(",
  "ng_rank_cross_priority(",
  "ng_write_cross_priority_workbook("
)

for (token in required_example_tokens) {
  if (!grepl(token, example_text, fixed = TRUE)) {
    fail("Chronological example missing token: %s", token)
  }
}

required_minimal_tokens <- c(
  "trait_value_metric <- \"var_complex\"",
  "allocation_method <- \"alphamate_style\"",
  "direction_column_col = \"PhenotypeColumn\"",
  "map_position_unit = \"bp\"",
  "ng_run_cross_prediction(",
  "stopifnot(inherits(result, \"ng_cross_prediction_result\"))"
)

for (token in required_minimal_tokens) {
  if (!grepl(token, minimal_example_text, fixed = TRUE)) {
    fail("Minimal example missing token: %s", token)
  }
}

if (!grepl("02_minimal_reproducible_user_run.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the minimal reproducible example")
}

required_variance_tokens <- c(
  "method_grid <- data.frame(",
  "trait_value_metric = c(",
  "\"usefulness\"",
  "\"pmv\"",
  "\"vpm\"",
  "\"le\"",
  "\"mean\"",
  "\"var_complex\"",
  "allocation_method = c(",
  "\"ocs\"",
  "\"alphamate_style\"",
  "ng_run_cross_prediction(",
  "variance_summary <- do.call(rbind,"
)

for (token in required_variance_tokens) {
  if (!grepl(token, variance_example_text, fixed = TRUE)) {
    fail("Variance-method example missing token: %s", token)
  }
}

if (!grepl("03_variance_method_comparison.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the variance-method comparison example")
}

required_external_alpha_tokens <- c(
  "allocation_method <- \"alphamate_executable\"",
  "alphamate_executable <- Sys.getenv(\"NG_ALPHAMATE_EXE\"",
  "alphamate_runtime_path <- Sys.getenv(\"NG_ALPHAMATE_RUNTIME_PATH\"",
  "alphamate_workdir <-",
  "alphamate_keep_files <- TRUE",
  "alphamate_mode <- \"ModeOptTarget1\"",
  "alphamate_target_degree <- 45",
  "alphamate_max_contributions <-",
  "alphamate_number_of_parents <-",
  "alphamate_evol_solutions <-",
  "alphamate_evol_iterations <-",
  "alphamate_evol_stop <-",
  "alphamate_n_threads <-",
  "file.exists(alphamate_executable)",
  "ng_run_cross_prediction(",
  "allocation_method = allocation_method",
  "alphamate_executable = alphamate_executable",
  "stopifnot(identical(result$settings$allocation_method, \"alphamate_executable\"))"
)

for (token in required_external_alpha_tokens) {
  if (!grepl(token, external_alpha_example_text, fixed = TRUE)) {
    fail("AlphaMate executable example missing token: %s", token)
  }
}

if (!grepl("04_alphamate_executable_user_run.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the AlphaMate executable example")
}

required_ocs_tokens <- c(
  "allocation_method <- \"ocs\"",
  "n_crosses <-",
  "max_crosses_per_parent <-",
  "min_unique_parents <-",
  "max_pair_kinship <-",
  "optimizer <- \"greedy_local\"",
  "use_ocs <- TRUE",
  "lambda_group <- 0.05",
  "lambda_mating <- 0.02",
  "lambda_parent_use <- 0",
  "lambda_parent_use_mode <- \"absolute\"",
  "local_iter <-",
  "ocs_iter <-",
  "ng_run_cross_prediction(",
  "allocation_method = allocation_method",
  "use_ocs = use_ocs",
  "lambda_parent_use = lambda_parent_use",
  "lambda_parent_use_mode = lambda_parent_use_mode",
  "stopifnot(identical(result$settings$allocation_method, \"ocs\"))"
)

for (token in required_ocs_tokens) {
  if (!grepl(token, ocs_example_text, fixed = TRUE)) {
    fail("OCS example missing token: %s", token)
  }
}

if (!grepl("05_ocs_user_run.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the OCS example")
}

required_optimizer_tokens <- c(
  "optimizer_grid <- data.frame(",
  "\"auto\"",
  "\"greedy_local\"",
  "\"repair_local\"",
  "\"mip_linear\"",
  "\"mip_contribution\"",
  "local_iter = c(",
  "ocs_iter = c(",
  "lambda_parent_use = c(",
  "lambda_parent_use_mode = c(",
  "requireNamespace(\"lpSolve\", quietly = TRUE)",
  "run_optimizer <- function(row)",
  "ng_run_cross_prediction(",
  "optimizer = row$optimizer",
  "local_iter = row$local_iter",
  "ocs_iter = row$ocs_iter",
  "optimizer_summary <- do.call(rbind,",
  "stopifnot(all(optimizer_summary$status == \"completed\"))",
  "stopifnot(identical(result$settings$allocation_method, \"ocs\"))"
)

for (token in required_optimizer_tokens) {
  if (!grepl(token, optimizer_example_text, fixed = TRUE)) {
    fail("Optimizer parameter guide example missing token: %s", token)
  }
}

if (!grepl("06_optimizer_parameter_guide.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the optimizer parameter guide example")
}

required_trait_value_parameter_tokens <- c(
  "metric_grid <- data.frame(",
  "\"var_complex\"",
  "\"usefulness\"",
  "\"pmv\"",
  "\"vpm\"",
  "\"le\"",
  "\"mean\"",
  "method_varPMV <- \"fast\"",
  "method_varPMV_choices <- c(\"fast\", \"full_posterior\")",
  "ril_mode <- \"infinite\"",
  "run_posterior_prediction <- FALSE",
  "posterior_method <- \"mcmc\"",
  "posterior_method_choices <- c(\"closed_form\", \"mcmc\")",
  "n_iter <- 5000",
  "burn_in <- 500",
  "n_draws <- max(1L, n_iter - burn_in)",
  "use_parallel <- FALSE",
  "parameter_notes <- data.frame(",
  "ng_run_cross_prediction(",
  "method_varPMV = row$method_varPMV",
  "ril_mode = row$ril_mode",
  "run_posterior_prediction = row$run_posterior_prediction",
  "posterior_method = row$posterior_method",
  "n_iter = row$n_iter",
  "burn_in = row$burn_in",
  "use_parallel = row$use_parallel",
  "trait_value_summary <- do.call(rbind,",
  "stopifnot(all(trait_value_summary$status == \"completed\"))",
  "ng_fit_ridge_effects(",
  "return_beta_cov_full = TRUE",
  "posterior_cov_full = effects_full$beta_cov_full",
  "ng_fit_ridge_effects_posterior(",
  "ng_posterior_cross_predict("
)

for (token in required_trait_value_parameter_tokens) {
  if (!grepl(token, trait_value_parameter_example_text, fixed = TRUE)) {
    fail("Trait-value parameter guide example missing token: %s", token)
  }
}

if (!grepl("07_trait_value_metric_parameter_guide.R", examples_readme_text, fixed = TRUE)) {
  fail("Examples README does not mention the trait-value parameter guide example")
}

required_prediction_trait_tokens <- c(
  "prediction_mode <- \"trait_by_trait\"",
  "traits_to_use <- c(",
  "direction_file",
  "direction_column_col = \"PhenotypeColumn\"",
  "multi_trait_method <- \"auto\"",
  "trait_value_metric <- \"var_complex\"",
  "ng_run_cross_prediction(",
  "prediction_mode = prediction_mode",
  "traits_to_use = traits_to_use",
  "stopifnot(identical(result$prediction_mode, \"trait_by_trait\"))"
)
for (token in required_prediction_trait_tokens) {
  if (!grepl(token, prediction_trait_example_text, fixed = TRUE)) {
    fail("Trait-by-trait prediction-mode example missing token: %s", token)
  }
}

required_prediction_index_tokens <- c(
  "prediction_mode <- \"index_as_trait\"",
  "index_col <- \"selection_index\"",
  "index_direction <- \"increase\"",
  "direction_file = NULL",
  "multi_trait_method <- \"auto\"",
  "ng_run_cross_prediction(",
  "prediction_mode = prediction_mode",
  "index_col = index_col",
  "index_direction = index_direction",
  "stopifnot(identical(result$prediction_mode, \"index_as_trait\"))"
)
for (token in required_prediction_index_tokens) {
  if (!grepl(token, prediction_index_example_text, fixed = TRUE)) {
    fail("Index-as-trait prediction-mode example missing token: %s", token)
  }
}

required_multitrait_auto_tokens <- c(
  "multi_trait_method <- \"auto\"",
  "trait_weights <- NULL",
  "threshold_policy <- \"soft\"",
  "ng_run_cross_prediction(",
  "multi_trait_method = multi_trait_method",
  "trait_weights = trait_weights",
  "stopifnot(identical(result$objective$method, \"auto\"))",
  "stopifnot(identical(result$objective$diagnostics$method_reason, \"equal_weight_rank_default\"))"
)
for (token in required_multitrait_auto_tokens) {
  if (!grepl(token, multitrait_auto_example_text, fixed = TRUE)) {
    fail("Auto multi-trait method example missing token: %s", token)
  }
}

required_multitrait_weighted_tokens <- c(
  "multi_trait_method <- \"weighted\"",
  "trait_weights <- c(",
  "ng_run_cross_prediction(",
  "multi_trait_method = multi_trait_method",
  "trait_weights = trait_weights",
  "stopifnot(identical(result$objective$method, \"weighted\"))",
  "stopifnot(identical(result$objective$diagnostics$method_reason, \"requested_weighted\"))"
)
for (token in required_multitrait_weighted_tokens) {
  if (!grepl(token, multitrait_weighted_example_text, fixed = TRUE)) {
    fail("Weighted multi-trait method example missing token: %s", token)
  }
}

required_multitrait_economic_tokens <- c(
  "multi_trait_method <- \"economic_index\"",
  "economic_weight = c(",
  "ng_run_cross_prediction(",
  "multi_trait_method = multi_trait_method",
  "stopifnot(identical(result$objective$method, \"economic_index\"))",
  "economic_index_coefficients"
)
for (token in required_multitrait_economic_tokens) {
  if (!grepl(token, multitrait_economic_example_text, fixed = TRUE)) {
    fail("Economic-index multi-trait method example missing token: %s", token)
  }
}

required_multitrait_desired_tokens <- c(
  "multi_trait_method <- \"desired_gain\"",
  "desired_change = c(",
  "phenotypic_covariance <-",
  "genetic_covariance <-",
  "ng_run_cross_prediction(",
  "multi_trait_method = multi_trait_method",
  "stopifnot(identical(result$objective$method, \"desired_gain\"))",
  "desired_gain_coefficients"
)
for (token in required_multitrait_desired_tokens) {
  if (!grepl(token, multitrait_desired_example_text, fixed = TRUE)) {
    fail("Desired-gain multi-trait method example missing token: %s", token)
  }
}

for (script in c(
  "08_prediction_mode_trait_by_trait.R",
  "09_prediction_mode_index_as_trait.R",
  "10_multitrait_method_auto.R",
  "11_multitrait_method_weighted.R",
  "12_multitrait_method_economic_index.R",
  "13_multitrait_method_desired_gain.R",
  "14_qc_duplicate_removal_and_reporting.R",
  "15_input_matching_and_map_units.R",
  "16_cross_number_sweep_diminishing_returns.R",
  "17_workbook_figures_and_priority_outputs.R",
  "18_allocation_method_comparison.R",
  "19_posterior_robust_mating_plan.R",
  "20_full_posterior_pmv_shortlist.R"
)) {
  if (!grepl(script, examples_readme_text, fixed = TRUE)) {
    fail("Examples README does not mention %s", script)
  }
}

required_qc_duplicate_tokens <- c(
  "duplicate_action <- \"remove\"",
  "duplicate_action_report <- \"report\"",
  "duplicate_threshold <- 0.995",
  "ng_preflight_input_tables(",
  "putative_duplicate_action = duplicate_action_report",
  "ng_plot_putative_duplicates(",
  "ng_run_cross_prediction(",
  "duplicate_action = duplicate_action",
  "stopifnot(\"putative_duplicate_genotypes_removed\" %in% result$qc$issues$id)"
)
for (token in required_qc_duplicate_tokens) {
  if (!grepl(token, qc_duplicate_example_text, fixed = TRUE)) {
    fail("QC duplicate example missing token: %s", token)
  }
}

required_input_matching_tokens <- c(
  "phenotype_id_col <- \"LineID\"",
  "genotype_id_col <- \"GenoID\"",
  "direction_trait_col <- \"TraitName\"",
  "direction_column_col <- \"PhenotypeColumn\"",
  "direction_direction_col <- \"SelectionGoal\"",
  "map_marker_col <- \"MarkerName\"",
  "map_position_unit <- \"bp\"",
  "bp_per_cm <- 1e6",
  "ng_run_cross_prediction(",
  "stopifnot(identical(result$input_match_audit$phenotype_id_col, phenotype_id_col))",
  "stopifnot(identical(result$input_match_audit$genotype_id_col, genotype_id_col))",
  "stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))"
)
for (token in required_input_matching_tokens) {
  if (!grepl(token, input_matching_example_text, fixed = TRUE)) {
    fail("Input-matching example missing token: %s", token)
  }
}

required_cross_number_tokens <- c(
  "K_range <- 2:6",
  "criterion <- \"elbow_relative\"",
  "relative_threshold <- 0.05",
  "ne_min <- 4",
  "coancestry_max <-",
  "ng_optimize_mating_plan_curve(",
  "ng_plot_diminishing_returns(",
  "recommended_K <- attr(curve, \"elbow_K\")"
)
for (token in required_cross_number_tokens) {
  if (!grepl(token, cross_number_example_text, fixed = TRUE)) {
    fail("Cross-number sweep example missing token: %s", token)
  }
}

required_workbook_figures_tokens <- c(
  "write_outputs <- requireNamespace(\"openxlsx\", quietly = TRUE)",
  "write_figures <- TRUE",
  "priority_breaks <- c(0.25, 0.50, 0.75, 1.00)",
  "priority_labels <- c(",
  "ng_run_cross_prediction(",
  "output_dir = output_dir",
  "write_outputs = write_outputs",
  "write_figures = write_figures",
  "stopifnot(file.exists(result$output_files$priority_score_vs_kinship_png))"
)
for (token in required_workbook_figures_tokens) {
  if (!grepl(token, workbook_figures_example_text, fixed = TRUE)) {
    fail("Workbook/figures example missing token: %s", token)
  }
}

required_allocation_comparison_tokens <- c(
  "allocation_grid <- data.frame(",
  "\"ocs\"",
  "\"alphamate_style\"",
  "\"alphamate_executable\"",
  "NG_ALPHAMATE_EXE",
  "ng_run_cross_prediction(",
  "allocation_method = row$allocation_method",
  "allocation_summary <- do.call(rbind,",
  "stopifnot(all(allocation_summary$status == \"completed\" | allocation_summary$status == \"skipped\"))"
)
for (token in required_allocation_comparison_tokens) {
  if (!grepl(token, allocation_comparison_example_text, fixed = TRUE)) {
    fail("Allocation comparison example missing token: %s", token)
  }
}

required_posterior_robust_tokens <- c(
  "run_posterior_prediction <- TRUE",
  "posterior_method <- \"closed_form\"",
  "n_iter <- 12",
  "burn_in <- 4",
  "ng_run_cross_prediction(",
  "result$posterior_predictions",
  "ng_optimize_robust_mating_plan(",
  "robustness_quantile <- 0.25",
  "objective <- \"posterior_quantile\"",
  "stopifnot(inherits(result, \"ng_cross_prediction_result\"))"
)
for (token in required_posterior_robust_tokens) {
  if (!grepl(token, posterior_robust_example_text, fixed = TRUE)) {
    fail("Posterior robust mating example missing token: %s", token)
  }
}

required_full_posterior_pmv_tokens <- c(
  "method_varPMV <- \"fast\"",
  "shortlist_size <-",
  "method_varPMV <- \"full_posterior\"",
  "ng_run_cross_prediction(",
  "yield_pmv_fast",
  "yield_pmv_full_posterior",
  "yield_pmv_used",
  "stopifnot(all(is.finite(full_result$candidate_crosses$yield_pmv_full_posterior)))"
)
for (token in required_full_posterior_pmv_tokens) {
  if (!grepl(token, full_posterior_pmv_example_text, fixed = TRUE)) {
    fail("Full-posterior PMV shortlist example missing token: %s", token)
  }
}

forbidden_fake_api <- c("ng_predict_crosses(", "ng_rank_crosses(")
for (token in forbidden_fake_api) {
  if (grepl(token, vignette_text, fixed = TRUE) || grepl(token, example_text, fixed = TRUE)) {
    fail("Found non-exported/fake API name in user docs: %s", token)
  }
}

cat("chronological vignette example tests passed\n")
