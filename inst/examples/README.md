# nextgenCrossDesign Examples

These examples are intended for users who want to run the package directly from
R.

- `00_user_friendly_cross_prediction.R`: first script to adapt for real
  phenotype, genotype, marker-map, and trait-direction files.
- `01_chronological_cross_prediction_workflow.R`: same workflow in a more
  detailed order, exposing QC, method choice, optimizer, allocation method,
  AlphaMate settings, priority tiers, and cross-number sweep settings.
- `02_minimal_reproducible_user_run.R`: self-contained smoke test that creates
  small example files, runs `ng_run_cross_prediction()` with
  `trait_value_metric = "var_complex"` and
  `allocation_method = "alphamate_style"`, and checks the output.
- `03_variance_method_comparison.R`: self-contained method comparison that
  runs the same input files with `var_complex`, `uc`, `pmv`, `vpm`,
  `le`, and `mean` trait-value settings, then writes a concise
  summary table.
- `04_alphamate_executable_user_run.R`: self-contained external AlphaMate
  example. It shows the exact `allocation_method = "alphamate_executable"`
  call and all executable-specific settings, and runs when `NG_ALPHAMATE_EXE`
  points to an installed AlphaMate binary.
- `05_ocs_user_run.R`: self-contained package OCS example. It shows the exact
  `allocation_method = "ocs"` call and the OCS-specific capacity, kinship,
  optimizer, and penalty settings.
- `06_optimizer_parameter_guide.R`: self-contained optimizer guide for OCS
  allocation. It walks through `auto`, `greedy_local`, `repair_local`,
  `mip_linear`, and `mip_contribution`; `lpSolve` is installed as a package
  dependency for the MIP optimizers.
- `07_trait_value_metric_parameter_guide.R`: self-contained guide for
  `trait_value_metric`, `method_varPMV`, `ril_mode`,
  `run_posterior_prediction`, `posterior_method`, `n_iter`, `burn_in`, and
  `use_parallel`. It runs several metrics and includes a small posterior
  prediction row so the settings are exercised through the package wrapper.
- `08_prediction_mode_trait_by_trait.R`: first prediction-mode example. Use
  it when phenotype files contain separate trait columns and a direction file.
- `09_prediction_mode_index_as_trait.R`: second prediction-mode example. Use
  it when the phenotype file already has a single selection-index column.
- `10_multitrait_method_auto.R`: first multi-trait method example. Uses only
  trait directions and lets the package build equal normalized directional
  weights.
- `11_multitrait_method_weighted.R`: second multi-trait method example. Shows
  user-supplied relative trait weights through `trait_weights`.
- `12_multitrait_method_economic_index.R`: third multi-trait method example.
  Shows `multi_trait_method = "economic_index"` with `economic_weight` in the
  direction file.
- `13_multitrait_method_desired_gain.R`: fourth multi-trait method example.
  Shows `multi_trait_method = "desired_gain"` with `desired_change` and
  `economic_weight` in the direction file.
- `14_qc_duplicate_removal_and_reporting.R`: QC example for putative
  duplicate genotype reporting, duplicate heatmap output, and removal through
  `duplicate_action = "remove"` inside the package workflow.
- `15_input_matching_and_map_units.R`: input-matching example using custom
  phenotype ID, genotype ID, direction-file, and marker-map columns. It also
  shows `map_position_unit = "bp"` and `bp_per_cm`.
- `16_cross_number_sweep_diminishing_returns.R`: portfolio-size example using
  `ng_optimize_mating_plan_curve()` and `ng_plot_diminishing_returns()` to
  choose the number of crosses by diminishing returns or diversity criteria.
- `17_workbook_figures_and_priority_outputs.R`: output example showing
  priority tiers, `priority_score_vs_kinship.png`, and the crossing-plan
  workbook when `openxlsx` is installed.
- `18_allocation_method_comparison.R`: allocation comparison across
  `ocs`, `alphamate_style`, and optional `alphamate_executable`.
- `19_posterior_robust_mating_plan.R`: posterior-aware example using
  `run_posterior_prediction = TRUE` and
  `ng_optimize_robust_mating_plan()`.
- `20_full_posterior_pmv_shortlist.R`: shortlist workflow that screens with
  `method_varPMV = "fast"` and reruns with
  `method_varPMV = "full_posterior"`.
- `21_gain_diversity_balance_and_target_coancestry.R`: steer the gain-vs-diversity
  trade-off with `strategy` / `diversity_emphasis` (a relative dial) or
  `target_coancestry` (constrained OCS: max gain at/under an inbreeding target), and
  export the frontier as JSON for a decision UI.
- `22_progeny_inbreeding_management.R`: penalize progeny inbreeding with
  `lambda_progeny_inbreeding`, read mean/max progeny F from `plan_summary`, and export a
  progeny-inbreeding histogram (JSON).
- `23_breeder_mating_constraints.R`: `min_crosses_per_parent` (min-use-if-used),
  `committed_crosses` (locked matings), and `parent_group` + `group_permission` /
  `group_quota` (mating-group rules), applied in the repair stage.
- `24_marker_steering_and_lethal_guarding.R`: drop carrier x carrier matings with
  `lethal_spec`, and drive a favourable marker allele's frequency with
  `ng_marker_target_spec()` + `ng_apply_marker_management()`.
- `25_cost_and_logistics_factors.R`: fold per-cross cost (`cost_col` + `budget` /
  `lambda_cost`) and a logistic/distance penalty (`logistic_col` + `lambda_logistic`) into
  allocation through the breeder path `ng_run_cross_prediction()` by supplying `cross_cost`
  (a per-cross data frame joined onto the candidate crosses).
- `26_evolution_optimizer.R`: the native memetic genetic-algorithm allocator via
  `optimizer = "evolution"` (with `evol_solutions` / `evol_iterations` / `evol_stop`).
- `27_breeder_controls_combined.R`: gain metric + strategy dial (`strategy = "balanced"`) +
  cost/`budget` + management constraints (`min_crosses_per_parent`, `parent_group` /
  `group_permission`, `lambda_progeny_inbreeding`) all in one `ng_run_cross_prediction()` run.
- `28_polyploid_mate_design.R`: any-ploidy polyploid mate design in one call --
  `ng_polyploid_design_crosses()` takes an allele-dosage matrix (0..ploidy) + effects or phenotype,
  runs ploidy-aware QC (`ng_polyploid_qc`), scores (mid-parent GEBV + a correct allele-frequency
  polyploid GRM, `ng_polyploid_grm`), and runs the full native control suite (strategy dial,
  committed matings, ...).
- `29_polyploid_qc_and_grm.R`: ploidy-aware QC (`ng_polyploid_qc`) and the correct polyploid
  relationship matrices -- additive `ng_polyploid_grm(method = "vanraden" | "yang")` and digenic
  dominance `ng_polyploid_dominance_grm()`.
- `30_polyploid_additive_dominance_effects.R`: additive + dominance genomic prediction --
  `ng_polyploid_fit_effects(model = "additive_dominance")` and `ng_polyploid_predict_value(type =
  "genotypic" | "breeding")`. Select clones on genotypic value, parents on breeding value.
- `31_polyploid_dominance_crossing.R`: dominance-aware mate design for clonal/heterosis crops
  (cassava, sugarcane) -- `ng_polyploid_score_crosses_dominance()` (heterosis mean + within-family
  variance) and `ng_polyploid_design_crosses(dominance = TRUE, gain =, double_reduction =, grm_method =)`.
- `32_exact_cross_trait_covariance.R`: exact recombination-aware within-family cross-trait
  covariance `ng_cross_trait_within_family_cov()` (a_t' R a_s) feeding the multi-trait threshold
  probability via `ng_add_p_superior_progeny_multitrait(cross_trait_cov = ...)`, versus the
  population-correlation proxy.
- `33_marker_effect_training_set.R`: enlarge the marker-effect training set with extra
  genotyped+phenotyped individuals that are NOT candidate parents --
  `ng_run_cross_prediction(training_genotype =, training_phenotype =)`. The parents are the main
  tables; training-only individuals sharpen the ridge fit and never appear in the crossing plan
  (`input_match_audit$training_only_count` / `$effect_training_n`).

For real data, do not manually intersect phenotype/genotype IDs or marker names
outside the package. Supply the correct ID, direction-file, and marker-map
column names and let package QC enforce exact matching before analysis.
