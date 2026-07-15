script <- file.path("tools", "run_ril_breeding_program_benchmark.R")
if (!file.exists(script)) {
  stop("Missing RIL breeding-program benchmark script: ", script, call. = FALSE)
}

txt <- paste(readLines(script, warn = FALSE), collapse = "\n")

required_tokens <- c(
  "RIL",
  "run_ril_breeding_program_benchmark",
  "ril_benchmark_method_registry",
  "advance_ril_pipeline",
  "F2",
  "F3",
  "F4",
  "F5",
  "F6",
  "stage_heritabilities",
  "within_family_selection_prop",
  "between_family_selection_prop",
  "AlphaSimR",
  "PopVar",
  "SimpleMating",
  "AlphaMate",
  "install_required_benchmark_dependencies",
  "ensure_r_github_package",
  "NG_POPVAR_GITHUB_REPO",
  "NG_SIMPLEMATING_GITHUB_REPO",
  "UMN-BarleyOatSilphium/PopVar",
  "Resende-Lab/SimpleMating",
  "multi_trait_method",
  "trait_by_trait",
  "index_as_trait",
  "index_direction",
  "trait_direction",
  "Selection_direction",
  "economic_weight",
  "desired_change",
  "NG_RIL_TRAITS",
  "NG_RIL_TRAIT_DIRECTIONS",
  "NG_RIL_MULTI_TRAIT_METHOD",
  "NG_RIL_PREDICTION_MODE",
  "_trait_spec.csv",
  "SimpleMating single-index comparison",
  "trait_value_metric",
  "uc_variance_source",
  "method_varPMV",
  "allocation_method",
  "exact_external_status",
  "software_versions",
  "NG_RIL_",
  # Prediction-accuracy evaluation (predicted vs realized, both shared-set and
  # range-restricted): the package must be graded on accuracy, not just gain.
  "evaluate_shared_accuracy",
  "selected_cross_accuracy",
  "realize_cross_families",
  "prediction_accuracy",
  "acc_pearson_top10",
  "acc_var_pearson",
  # Per-cycle breeder metrics: genetic gain rising, diversity eroding.
  "genetic_gain_index",
  "expected_heterozygosity",
  "prop_polymorphic_markers",
  "mean_parent_relationship",
  "population_diversity",
  # Feed single-trait external tools the package's exact composite index.
  "package_composite_index",
  "external_training_index",
  # Complete PopVar/SimpleMating criterion coverage (no underselling the comparison).
  "popvar_mu_topn",
  "popvar_var_topn",
  "popvar_musp_topn",
  "simple_mpv_topn",
  "simple_mpv_select",
  # Replicates run in parallel; cycles stay sequential (build on each other).
  "run_single_rep",
  "parallel_lapply",
  "parallel_scope",
  "mclapply",
  # native evolutionary (memetic GA) allocator, head-to-head with OCS/MIP
  "var_complex_evolution",
  "evolution",
  "evol_solutions",
  "evol_iterations",
  # gain-diversity balancing dial exercised at three frontier points
  "strategy_high_gain",
  "strategy_balanced",
  "strategy_diversity",
  "strategy = sub",
  # the package driven through its top-level user API, like a user (and like SimpleMating)
  "package_user_api",
  "user_api",
  "ng_run_cross_prediction("
)

for (token in required_tokens) {
  if (!grepl(token, txt, fixed = TRUE)) {
    stop("RIL benchmark script is missing required token: ", token, call. = FALSE)
  }
}

current_example_methods <- c(
  "var_complex_ocs",
  "var_complex_alphamate_style",
  "uc_pmv_ocs",
  "uc_vpm_ocs",
  "pmv_ocs",
  "vpm_ocs",
  "var_simple_ocs",
  "mean_ocs"
)

for (method in current_example_methods) {
  if (!grepl(method, txt, fixed = TRUE)) {
    stop("RIL benchmark script is missing current package/example method: ", method, call. = FALSE)
  }
}

stale_default_methods <- c(
  "nextgen_ril_uc_ocs",
  "nextgen_ril_recomb_ocs",
  "nextgen_ril_var_topn"
)

for (method in stale_default_methods) {
  if (grepl(method, txt, fixed = TRUE)) {
    stop("RIL benchmark script still contains stale benchmark method label: ", method, call. = FALSE)
  }
}

cat("ril breeding-program benchmark contract passed\n")
