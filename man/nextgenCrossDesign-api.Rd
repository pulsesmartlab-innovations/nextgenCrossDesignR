\name{nextgenCrossDesign-api}
\alias{nextgenCrossDesign-api}
\alias{ng_evolutionary_mate_allocation}
\alias{ng_with_rng_seed}
\alias{ng_choose_frontier_point}
\alias{ng_select_by_strategy}
\alias{ng_select_by_target_coancestry}
\alias{ng_progeny_inbreeding_histogram}
\alias{ng_apply_marker_management}
\alias{ng_marker_target_spec}
\alias{ng_marker_target_scores}
\alias{ng_lethal_recessive_spec}
\alias{ng_lethal_recessive_cross_risk}
\alias{ng_trait_check_spec}
\alias{ng_apply_trait_checks}
\alias{ng_frontier_export_payload}
\alias{ng_write_frontier_json}
\alias{ng_write_progeny_inbreeding_histogram_json}
\alias{ng_run_comparison_payload}
\alias{ng_write_run_comparison_json}
\alias{ng_add_balanced_usefulness_score}
\alias{ng_add_external_baseline_scores}
\alias{ng_add_gms_vpm_scores}
\alias{ng_add_multitrait_score}
\alias{ng_add_popvar_scores}
\alias{ng_add_simplemating_scores}
\alias{ng_polyploid_score_crosses}
\alias{ng_polyploid_design_crosses}
\alias{ng_polyploid_grm}
\alias{ng_polyploid_dominance_grm}
\alias{ng_polyploid_qc}
\alias{ng_polyploid_fit_effects}
\alias{ng_polyploid_predict_value}
\alias{ng_polyploid_score_crosses_dominance}
\alias{ng_polyploid_progeny_moment_table}
\alias{ng_polyploid_allele_freq}
\alias{ng_alphamate_style_select}
\alias{ng_allocate_family_sizes}
\alias{ng_alphamate_default_executable}
\alias{ng_apply_family_variance_calibrators}
\alias{ng_backend_capability_registry}
\alias{ng_breeder_selection_objective}
\alias{ng_compare_selected_families}
\alias{ng_crop_aware_policy_select}
\alias{ng_crop_genome_env}
\alias{ng_crop_genome_scenarios}
\alias{ng_crop_genome_select}
\alias{ng_cross_priority_summary}
\alias{ng_cross_priority_workbook_tables}
\alias{ng_design_crosses}
\alias{ng_detect_putative_duplicates}
\alias{ng_build_cross_trait_covariance}
\alias{ng_cross_trait_within_family_cov}
\alias{ng_estimate_genetic_covariance}
\alias{ng_estimate_phenotypic_covariance}
\alias{ng_posterior_genetic_covariance}
\alias{ng_posterior_multitrait_cross_predict}
\alias{ng_expected_excess_above_threshold}
\alias{ng_external_methods_needed}
\alias{ng_family_metric_registry}
\alias{ng_family_metric_summary}
\alias{ng_family_selection_threshold}
\alias{ng_fit_family_variance_calibrators}
\alias{ng_fit_ridge_effects}
\alias{ng_fit_ridge_effects_posterior}
\alias{ng_frontier_default_policy}
\alias{ng_frontier_policy_from_spec}
\alias{ng_frontier_policy_select}
\alias{ng_head_to_head_compare_summary}
\alias{ng_head_to_head_dashboard_payload}
\alias{ng_head_to_head_default_methods}
\alias{ng_ld_prune_geno}
\alias{ng_ld_prune_markers}
\alias{ng_load}
\alias{ng_load_cpp}
\alias{ng_make_pairs}
\alias{ng_multitrait_crop_validation_scenario}
\alias{ng_multitrait_policy_select}
\alias{ng_multitrait_select_topn}
\alias{ng_multitrait_spec}
\alias{ng_multitrait_validation_default_methods}
\alias{ng_multitrait_validation_scenario}
\alias{ng_optimize_balanced_usefulness}
\alias{ng_optimize_mating_plan}
\alias{ng_optimize_mating_plan_curve}
\alias{ng_plot_diminishing_returns}
\alias{ng_plot_priority_score_vs_kinship}
\alias{ng_plot_putative_duplicates}
\alias{ng_optimize_breeder_selection_plan}
\alias{ng_optimize_multitrait_mating_plan}
\alias{ng_optimize_robust_mating_plan}
\alias{ng_pareto_mate_allocation}
\alias{ng_parent_kinship}
\alias{ng_p_superior_progeny}
\alias{ng_p_superior_progeny_multitrait}
\alias{ng_add_p_superior_progeny}
\alias{ng_add_p_superior_progeny_multitrait}
\alias{ng_posterior_cross_predict}
\alias{ng_sample_ridge_posterior_bcm}
\alias{ng_sample_ridge_posterior_mcmc}
\alias{ng_polyploid_as_dosage_matrix}
\alias{ng_polyploid_model_select}
\alias{ng_polyploid_policy}
\alias{ng_polyploid_policy_modes}
\alias{ng_polyploid_subgenome_as_dosage_list}
\alias{ng_polyploid_subgenome_grm}
\alias{ng_polyploid_subgenome_score_crosses}
\alias{ng_predict_gebv}
\alias{ng_preflight_input_tables}
\alias{ng_prepare_marker_map}
\alias{ng_run_alphamate}
\alias{ng_run_cross_prediction}
\alias{ng_run_head_to_head_benchmark}
\alias{ng_run_multitrait_crop_validation_grid}
\alias{ng_run_multitrait_validation}
\alias{ng_run_multitrait_validation_grid}
\alias{ng_rank_cross_priority}
\alias{ng_score_crosses}
\alias{ng_score_breeder_objective}
\alias{ng_select_alphamate}
\alias{ng_simplemating_build_crosses}
\alias{ng_simplemating_past_thinning}
\alias{ng_simplemating_relate_thinning}
\alias{ng_simplemating_select_crosses_native}
\alias{ng_simplemating_style_select}
\alias{ng_validate_metric_calibration}
\alias{ng_popvar_style_scores}
\alias{ng_write_backend_capability_registry_json}
\alias{ng_write_crossing_plan_workbook}
\alias{ng_write_cross_priority_workbook}
\alias{ng_write_data_preflight_json}
\alias{ng_write_head_to_head_dashboard_json}
\alias{ng_write_head_to_head_visual_report}
\title{Public NextGen Cross Design Backend API}
\description{
Public backend functions for cross scoring, mate allocation, multi-trait
objectives, validation, benchmark reporting, external baseline comparison,
data preflight, and scoped polyploid workflows.
}
\value{
Most functions return data frames, lists, or file paths containing cross scores,
mating plans, validation summaries, benchmark reports, or JSON-ready payloads.
}
\details{
This manual page documents the package's public API surface. Internal helpers
remain available within the namespace but are intentionally not exported.

Native external-tool-inspired components are integrated directly in the package
so the backend does not require PopVar, SimpleMating, or AlphaMate to operate.
Use \code{ng_run_cross_prediction(trait_value_metric = "var_complex")} for the
native complex usefulness metric, \code{ng_popvar_style_scores()} for
lower-level PopVar-style compatibility scoring,
\code{ng_simplemating_relate_thinning()}, \code{ng_simplemating_build_crosses()},
\code{ng_simplemating_past_thinning()}, and
\code{ng_simplemating_select_crosses_native()} for SimpleMating-style parent and
cross handling, and \code{ng_alphamate_style_select()} for AlphaMate-style
gain-diversity target allocation. Exact external calls remain optional through
the PopVar, SimpleMating, and AlphaMate integration wrappers, or through
\code{ng_run_cross_prediction(allocation_method = "alphamate_executable")} for
AlphaMate.

\code{ng_run_cross_prediction()} exposes the main trait-value and uncertainty
controls directly: \code{trait_value_metric}, \code{uc_variance_source},
\code{method_varPMV = "fast"} or \code{"full_posterior"}, \code{ril_mode =
"infinite"}, \code{run_posterior_prediction}, \code{posterior_method},
\code{n_iter}, \code{burn_in}, and \code{use_parallel}. When
\code{method_varPMV = "full_posterior"}, the wrapper fits and uses the full
marker-effect covariance in PMV scoring. When
\code{run_posterior_prediction = TRUE}, posterior cross-prediction tables are
returned in \code{result$posterior_predictions}.

Use \code{ng_write_crossing_plan_workbook()} to create breeder-facing Excel
workbooks with top parents, recommended crosses, full-sib screening, editable
criteria, and an embedded selection-index chart.

Use \code{ng_ld_prune_markers()} and \code{ng_ld_prune_geno()} to remove
low-MAF and high-LD redundant markers before formal model fitting or cross
design. The pruning workflow has a native C++ backend for large dosage matrices
and an R fallback.

Use \code{ng_preflight_input_tables()} before scoring. Set
\code{putative_duplicate_check = TRUE} to run the package IBS 0-1 screen for
near-identical genotype profiles. Set
\code{putative_duplicate_action = "remove"} to keep one parent per duplicate
cluster and return cleaned input tables in \code{qc$cleaned_tables}; the audit
metadata is returned in \code{qc$cleaning$putative_duplicates}. Use
\code{ng_write_data_preflight_json()} for compact UI or batch logs, and use
\code{ng_detect_putative_duplicates()} plus \code{ng_plot_putative_duplicates()}
when duplicate-pair tables or duplicate heatmaps are needed for QC review.

Use \code{ng_rank_cross_priority()} after selecting a practical crossing list
such as 100 crosses to add \code{priority_rank}, \code{priority_index}, and
breeder-facing \code{priority_tier} values. The default tiers are
\code{highly_priority}, \code{priority}, \code{medium_priority}, and
\code{low_priority}; use \code{ng_cross_priority_summary()} to summarize the
number of selected crosses in each tier.

Use \code{ng_write_cross_priority_workbook()} to create a breeder-facing Excel
workbook from a selected crossing plan. The workbook reports objective evidence
such as priority tier, score, kinship, threshold violation, parent-use QC, and
duplicate QC, while leaving \code{Breeder_Rationale},
\code{Breeder_Notes}, \code{Final_Decision}, and \code{Crossing_Status}
editable because final rationale differs by breeding program.

Use \code{ng_plot_priority_score_vs_kinship()} to show all candidate crosses in
score-versus-relatedness space and overlay the selected crossing plan colored by
priority tier.
}
\keyword{methods}
