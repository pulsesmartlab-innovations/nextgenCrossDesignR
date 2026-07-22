ng_backend_registry_timestamp <- function(generated_at = Sys.time()) {
  format(as.POSIXct(generated_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

# Enumerable user-facing controls: the authoritative choice lists + defaults a frontend
# should render its dropdowns from, so new/renamed methods appear without editing the UI.
# `group` maps to a workbench panel; `capability` links to method_families; `depends_on`
# gates visibility. Values mirror ng_run_cross_prediction()'s match.arg() choices.
ng_backend_controls <- function() {
  ch <- function(x) unname(Map(function(v, l) list(value = v, label = unname(l)), names(x), x))
  enum <- function(id, label, group, default, x, capability = NA_character_,
                   depends_on = NA_character_) {
    list(id = id, label = label, group = group, type = "enum", default = default,
         capability = capability, depends_on = depends_on, choices = ch(x))
  }
  list(
    enum("map_position_unit", "Map position unit", "data", "bp",
      c(bp = "Base pairs (bp)", cM = "centiMorgans (cM)")),
    enum("prediction_mode", "Prediction mode", "objective", "trait_by_trait",
      c(trait_by_trait = "Trait by trait", index_as_trait = "Index as a single trait")),
    enum("multi_trait_method", "Multi-trait method", "objective", "auto",
      c(auto = "Auto (equal weights from directions)", weighted = "Weighted index",
        economic_index = "Economic index", desired_gain = "Desired gain",
        threshold = "Threshold (soft by default)"), capability = "multitrait_auto"),
    enum("threshold_policy", "Threshold policy", "objective", "soft",
      c(soft = "Soft (penalty)", strict = "Strict (hard cutoff)"),
      depends_on = "multi_trait_method=threshold"),
    enum("trait_value_metric", "Cross-scoring metric", "scoring", "usefulness",
      c(usefulness = "Usefulness (mean + i * within-family SD)",
        pmv = "Usefulness via PMV variance", vpm = "Usefulness via VPM (recombination) variance",
        le = "Relationship-distance proxy (screening only)",
        var_complex = "var_complex (native PopVar-style usefulness)",
        mean = "Cross mean only"), capability = "dh_ril_pmv_scoring"),
    enum("uc_variance_source", "Usefulness variance source", "scoring", "pmv",
      c(pmv = "PMV", vpm = "VPM (recombination)", le = "Relationship distance"),
      depends_on = "trait_value_metric=usefulness"),
    enum("method_varPMV", "PMV method", "scoring", "fast",
      c(fast = "Fast (point estimate)", full_posterior = "Full posterior")),
    enum("progeny", "Progeny system", "scoring", "DH",
      c(DH = "DH (doubled haploid)", DHs = "DHs", RIL = "RIL", RILs = "RILs")),
    enum("ril_mode", "RIL selfing model", "scoring", "infinite",
      c(infinite = "Infinite (F-infinity)", finite = "Finite selfing"),
      depends_on = "progeny=RIL|RILs"),
    enum("recomb_model", "Recombination map function", "scoring", "haldane",
      c(haldane = "Haldane", kosambi = "Kosambi")),
    enum("grm_method", "GRM method", "scoring", "vanraden",
      c(vanraden = "VanRaden", yang = "Yang")),
    enum("duplicate_action", "Duplicate handling", "qc", "remove",
      c(remove = "Remove duplicates", report = "Report only", none = "Ignore")),
    enum("ld_backend", "LD-pruning backend", "qc", "auto",
      c(auto = "Auto", cpp = "C++", r = "R")),
    enum("optimizer", "Optimizer", "allocation", "auto",
      c(auto = "Auto", evolution = "Evolution (recommended)", greedy_local = "Greedy local",
        repair_local = "Repair local", mip_linear = "MIP (linear)",
        mip_contribution = "MIP (contribution)"), capability = "ocs_allocation"),
    enum("allocation_method", "Allocation method", "allocation", "ocs",
      c(ocs = "OCS (coancestry-constrained)", alphamate_style = "AlphaMate-style (native)",
        alphamate_executable = "AlphaMate executable (external)"), capability = "ocs_allocation"),
    enum("strategy", "Gain-diversity strategy", "allocation", "balanced",
      c(balanced = "Balanced", high_gain = "High gain", diversity = "High diversity"),
      capability = "gain_diversity_balance"),
    enum("lambda_parent_use_mode", "Parent-use penalty mode", "allocation", "absolute",
      c(absolute = "Absolute", adaptive = "Adaptive")),
    enum("cross_sweep_criterion", "Cross-number recommendation rule", "allocation", "elbow_relative",
      c(elbow_relative = "Diminishing returns (relative)", elbow_kneedle = "Diminishing returns (kneedle)",
        ne_target = "Effective population size target", coancestry_budget = "Coancestry budget")),
    enum("alphamate_mode", "AlphaMate mode", "allocation", "ModeOptTarget1",
      c(ModeOptTarget1 = "Optimize to a target degree", ModeMaxCriterion = "Maximize criterion",
        ModeMinCoancestry = "Minimize coancestry"), depends_on = "allocation_method=alphamate_style|alphamate_executable"),
    enum("posterior_method", "Posterior method", "advanced", "mcmc",
      c(mcmc = "MCMC", closed_form = "Closed form"), depends_on = "run_posterior_prediction=true")
  )
}

ng_backend_capability_registry <- function(generated_at = Sys.time()) {
  data_qc <- data.frame(
    id = c(
      "duplicate_parent_ids",
      "putative_duplicate_genotypes_removed",
      "duplicate_marker_ids",
      "duplicate_phenotype_rows",
      "duplicate_trait_names",
      "duplicate_candidate_crosses",
      "reciprocal_candidate_crosses",
      "genotype_phenotype_id_mismatch",
      "missing_marker_map_entries",
      "ld_redundant_markers",
      "invalid_polyploid_dosage",
      "parent_relationship_id_mismatch"
    ),
    label = c(
      "Duplicate parent IDs",
      "Putative duplicate genotype removal",
      "Duplicate marker IDs",
      "Duplicate phenotype rows",
      "Duplicate trait names",
      "Duplicate candidate crosses",
      "Reciprocal duplicate crosses",
      "Genotype/phenotype ID mismatch",
      "Missing marker-map entries",
      "LD-redundant markers",
      "Invalid polyploid dosage",
      "Parent relationship ID mismatch"
    ),
    severity = c("blocker", "warning", "blocker", "warning", "blocker", "blocker", "blocker", "blocker", "blocker", "warning", "blocker", "blocker"),
    backend_status = c(
      "implemented_core",
      "implemented_import_preflight",
      "implemented_map_polyploid",
      "implemented_import_preflight",
      "implemented_multitrait",
      "implemented_polyploid_pairs",
      "implemented_polyploid_pairs",
      "implemented_core",
      "implemented_marker_map",
      "implemented_native_cpp",
      "implemented_polyploid",
      "implemented_polyploid_relationship"
    ),
    backend_function = c(
      "ng_check_same_ids",
      "ng_preflight_input_tables; ng_detect_putative_duplicates",
      "ng_prepare_marker_map",
      "frontend_preflight_needed",
      "ng_multitrait_spec",
      "ng_poly_validate_candidate_pairs",
      "ng_poly_validate_candidate_pairs",
      "ng_check_same_ids",
      "ng_prepare_marker_map",
      "ng_ld_prune_markers; ng_ld_prune_geno; ng_run_cross_prediction(ld_pruning=); ng_design_crosses(ld_pruning=)",
      "ng_polyploid_as_dosage_matrix",
      "ng_polyploid_grm; ng_polyploid_subgenome_grm"
    ),
    ui_stage = rep("Data QC", 12L),
    description = c(
      "Parent/sample identifiers must be unique before scoring or allocation.",
      "Putative duplicate genotype profiles can be detected and redundant parents removed inside package QC before scoring.",
      "Marker identifiers must be unique and aligned with effects and maps.",
      "Phenotype imports should flag repeated observations for the same parent and trait.",
      "Multi-trait objectives require unique non-empty trait names.",
      "Candidate-pair tables must not repeat the same cross.",
      "Crosses A x B and B x A are the same mating plan entry for these workflows.",
      "Genotype, phenotype, and requested parent IDs must align before scoring.",
      "Every marker used in genotype/effect matrices must exist in the marker map.",
      "Optional native LD pruning removes low-MAF and high-LD redundant markers before formal analysis.",
      "Dosage matrices must stay within the declared ploidy range.",
      "Relationship matrix row and column IDs must be unique and aligned."
    ),
    stringsAsFactors = FALSE
  )

  breeding_systems <- data.frame(
    id = c(
      "diploid_dh",
      "diploid_ril",
      "diploidized_crop_stress",
      "autotetraploid_4x",
      "allopolyploid_subgenome",
      "complex_polyploid_guard"
    ),
    label = c(
      "Diploid DH",
      "Diploid RIL",
      "Diploidized crop stress",
      "Autotetraploid 4x",
      "Allopolyploid subgenome",
      "Complex polyploid guard"
    ),
    category = c("diploid", "diploid", "stress", "polyploid", "polyploid", "guarded"),
    backend_target = c("DH", "RIL", "DH", "4x_dosage", "disomic_subgenome", "unsupported_without_empirical_generator"),
    default = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    validation_scope = c(
      "validated_dh_ril",
      "validated_dh_ril",
      "diploidized_portability_screen",
      "true_autotetraploid_4x_experimental",
      "wheat_like_disomic_phase2a",
      "guarded_not_claimed"
    ),
    status = c("default", "supported", "stress_screen", "experimental", "experimental", "guarded"),
    backend_function = c(
      "ng_score_crosses(target = 'DH')",
      "ng_score_crosses(target = 'RIL')",
      "ng_crop_genome_select",
      "ng_polyploid_policy",
      "ng_poly_subgenome_policy",
      "ng_polyploid_model_select"
    ),
    description = c(
      "Primary validated inbred-line workflow for doubled-haploid family decisions.",
      "Primary diploid recombinant inbred line target using the same scoring API.",
      "Crop portability stress screen for wheat, potato, cassava, and sugarcane approximations.",
      "True AlphaSimR 4x dosage workflow for potato/cassava-like autotetraploid experiments.",
      "Named subgenome disomic dosage workflow for true allopolyploids; recombination-aware within-family variance (per-subgenome a'Ra) when a chromosome+cM map is supplied, else the linkage-equilibrium approximation.",
      "Complex or aneuploid polyploids are guarded until an empirical generator is supplied."
    ),
    stringsAsFactors = FALSE
  )

  method_families <- data.frame(
    id = c(
      "user_friendly_cross_prediction",
      "multitrait_auto",
      "multitrait_weighted",
      "multitrait_economic_index",
      "multitrait_desired_gain",
      "multitrait_threshold",
      "dh_ril_pmv_scoring",
      "frontier_policy",
      "ocs_allocation",
      "family_size_allocation",
      "external_baselines",
      "native_external_components",
      "crop_aware_policy",
      "poly4x_policy",
      "poly_subgenome_policy",
      "gain_diversity_balance",
      "progeny_inbreeding_management",
      "mating_constraints",
      "marker_steering",
      "lethal_allele_guarding",
      "cost_logistics",
      "polyploid_grm_qc",
      "polyploid_dominance_design",
      "training_set_augmentation"
    ),
    label = c(
      "User-friendly cross prediction workflow",
      "Auto multi-trait default",
      "Weighted index",
      "Economic index",
      "Desired-gain optimizer",
      "Threshold selection",
      "DH/RIL PMV and usefulness scoring",
      "DH/RIL frontier policy",
      "OCS constrained allocation",
      "Family-size allocation",
      "External baseline comparison",
      "Native external-tool components",
      "Crop-aware policy",
      "Autotetraploid 4x policy",
      "Allopolyploid subgenome policy",
      "Gain-diversity balance (strategy dial)",
      "Progeny-inbreeding management",
      "Breeder mating constraints",
      "Marker steering",
      "Lethal-allele guarding",
      "Cost / logistics factors",
      "Polyploid GRM (VanRaden/Yang) + ploidy-aware QC",
      "Polyploid dominance-aware mate design",
      "Marker-effect training-set augmentation"
    ),
    category = c(
      "workflow",
      "multi_trait",
      "multi_trait",
      "multi_trait",
      "multi_trait",
      "multi_trait",
      "scoring",
      "policy",
      "allocation",
      "allocation",
      "baseline",
      "baseline",
      "policy",
      "polyploid",
      "polyploid",
      "allocation",
      "scoring",
      "allocation",
      "scoring",
      "scoring",
      "allocation",
      "polyploid",
      "polyploid",
      "scoring"
    ),
    default = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    backend_function = c(
      "ng_run_cross_prediction",
      "ng_add_multitrait_score(method = 'auto')",
      "ng_add_multitrait_score(method = 'weighted')",
      "ng_add_multitrait_score(method = 'economic_index')",
      "ng_add_multitrait_score(method = 'desired_gain')",
      "ng_add_multitrait_score(method = 'threshold')",
      "ng_score_crosses(grm_method=); ng_parent_kinship(method=)",
      "ng_frontier_policy_ocs10_lps2",
      "ng_optimize_mating_plan",
      "ng_allocate_family_sizes",
      "ng_external_baseline_scores",
      "ng_run_cross_prediction(trait_value_metric = 'var_complex'); ng_simplemating_*; ng_alphamate_style_select",
      "ng_crop_genome_select",
      "ng_polyploid_policy",
      "ng_poly_subgenome_policy",
      "ng_optimize_mating_plan(strategy=, diversity_emphasis=, target_coancestry=); ng_select_by_strategy; ng_select_by_target_coancestry; ng_write_frontier_json",
      "ng_score_crosses (expected_progeny_inbreeding); ng_optimize_mating_plan(lambda_progeny_inbreeding=); ng_progeny_inbreeding_histogram",
      "ng_optimize_mating_plan(min_crosses_per_parent=, committed_crosses=, parent_group=, group_permission=, group_quota=)",
      "ng_marker_target_spec; ng_marker_target_scores; ng_apply_marker_management; ng_run_cross_prediction(marker_target_spec=, lambda_marker=); ng_design_crosses(marker_target_spec=, lambda_marker=)",
      "ng_lethal_recessive_spec; ng_lethal_recessive_cross_risk",
      "ng_optimize_mating_plan(cost_col=, budget=, lambda_cost=, logistic_col=, lambda_logistic=); ng_run_cross_prediction(cross_cost=, cost_col=, budget=, lambda_cost=, logistic_col=, lambda_logistic=)",
      "ng_polyploid_grm(method=); ng_polyploid_dominance_grm(method=); ng_polyploid_qc",
      "ng_polyploid_design_crosses(dominance=, gain=, double_reduction=, grm_method=); ng_polyploid_fit_effects; ng_polyploid_predict_value; ng_polyploid_score_crosses_dominance",
      "ng_run_cross_prediction(training_genotype=, training_phenotype=, training_genotype_file=, training_phenotype_file=)"
    ),
    evidence = c(
      "End-to-end breeder workflow for QC, duplicate removal, trait-by-trait marker effects, cross prediction, multi-trait scoring, OCS allocation, priority tiers, and output files.",
      "Rank-normalized default for unknown breeder weights with threshold penalties.",
      "Uses declared relative trait weights after orienting increase/decrease directions.",
      "Covariance-aware economic coefficients estimated from oriented/scaled trait covariance.",
      "Full economic desired-gain solver using desired changes, economic weights, stabilized covariance, and diagnostic responses.",
      "Soft thresholds by default with optional strict filtering.",
      "Diploid DH/RIL recombination-aware PMV, usefulness, kinship, and reliability diagnostics.",
      "Current validated DH/RIL parent-count decision rule from replicated evidence.",
      "Gain/diversity/parent-use constrained mating plans with greedy, repair, and MIP modes.",
      "Allocates progeny counts across selected families using score-weighted or marginal top-k logic.",
      "PopVar, SimpleMating, AlphaMate exact/proxy baseline surfaces with fallback status.",
      "Native var_complex usefulness, SimpleMating thinning/selection, and AlphaMate-style target-degree allocation, without requiring those tools.",
      "Routes crop/ploidy stress screens to the best validated available family.",
      "True 4x policy modes for gain, OCS, and diversity experiments.",
      "Wheat-like disomic subgenome gain, usefulness, diversity, and OCS policy.",
      "Single-dial gain-vs-diversity balance (0=gain..100=diversity) over the native frontier with bisection refinement and strategy presets.",
      "Expected progeny inbreeding as a first-class metric, optional objective emphasis distinct from parental coancestry, and a distribution histogram.",
      "Repair-stage min-use-if-used, committed/fixed matings, and mating-group permission matrix + quotas; returns a shorter plan with a warning when infeasible.",
      "Steer a nominated marker's allele frequency up/down or toward a target; ploidy-aware expected progeny frequency.",
      "Flag and optionally exclude carrier x carrier matings at deleterious recessive loci.",
      "Optional per-cross cost (soft penalty + hard budget) and logistic/geographic penalty folded into the allocation objective.",
      "Correct allele-frequency-based polyploid additive and dominance GRMs (VanRaden or Yang/GCTA, generalized to ploidy) plus ploidy-aware QC (0..ploidy range, missingness, MAF, monomorphic, duplicates); replaces the diploid-oriented kinship.",
      "Any-ploidy one-call mate design with optional additive+dominance modelling for clonal/heterosis crops (cassava, sugarcane, potato): genotypic-value parent selection, heterosis-inclusive cross mean + within-family additive+dominance variance, optional double reduction, C++-accelerated. Additive-only is the default.",
      "Enlarge the marker-effect training set with extra genotyped+phenotyped individuals that are NOT candidate parents (the parents are the main genotype/phenotype tables). Improves effect-based metrics (mean/usefulness/pmv/vpm/var_complex), not le. The result reports training_only_count, effect_training_n, and training_ids so the frontend can display 'trained on N, crossing K parents'."
    ),
    stringsAsFactors = FALSE
  )

  workflows <- data.frame(
    id = c(
      "user_friendly_cross_prediction",
      "export_backend_capabilities",
      "head_to_head",
      "export_dashboard_json",
      "render_visual_report",
      "crossing_plan_workbook",
      "multitrait_validation",
      "multitrait_validation_grid",
      "multitrait_crop_validation_grid",
      "crop_genome_scenarios",
      "family_calibration",
      "poly4x_benchmark",
      "poly4x_controlled_ocs",
      "alphasimr_benchmark"
    ),
    label = c(
      "User-friendly cross prediction",
      "Export backend capability registry",
      "Head-to-head benchmark",
      "Export dashboard JSON",
      "Render visual report",
      "Crossing-plan Excel workbook",
      "Multi-trait validation",
      "Multi-trait validation grid",
      "Multi-trait crop validation grid",
      "Crop-genome scenarios",
      "Family calibration benchmark",
      "Autotetraploid 4x benchmark",
      "Controlled 4x OCS comparison",
      "AlphaSimR DH/RIL benchmark"
    ),
    category = c("example", "metadata", "benchmark", "report", "report", "report", "validation", "validation", "validation", "crop", "calibration", "polyploid", "polyploid", "benchmark"),
    script = c(
      "inst/examples/00_user_friendly_cross_prediction.R",
      "tools/export_backend_capabilities_json.R",
      "tools/run_head_to_head_benchmark.R",
      "tools/export_head_to_head_dashboard_json.R",
      "tools/render_head_to_head_visual_report.R",
      "tools/export_crossing_plan_workbook.R",
      "tools/run_multitrait_validation.R",
      "tools/run_multitrait_validation_grid.R",
      "tools/run_multitrait_crop_validation_grid.R",
      "tools/run_crop_genome_scenarios.R",
      "tools/run_family_calibration_benchmark.R",
      "tools/run_poly4x_benchmark.R",
      "tools/run_poly4x_controlled_ocs.R",
      "tools/run_alphasimr_benchmark.R"
    ),
    job_type = c(
      "user_cross_prediction",
      "export_backend_capabilities",
      "head_to_head",
      "export_dashboard_json",
      "render_visual_report",
      "crossing_plan_workbook",
      "multitrait_validation",
      "multitrait_validation_grid",
      "multitrait_crop_validation_grid",
      "crop_genome_scenarios",
      "family_calibration",
      "poly4x_benchmark",
      "poly4x_controlled_ocs",
      "alphasimr_benchmark"
    ),
    default = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
    estimated_runtime = c("short", "short", "medium", "short", "short", "short", "short", "medium", "medium", "medium", "medium", "long", "long", "long"),
    description = c(
      "Runs the package-first user workflow from input files to selected priority-ranked crosses.",
      "Writes the frontend capability contract used by the workbench.",
      "Compares NextGen candidates against optional PopVar/SimpleMating/AlphaMate-style baselines.",
      "Converts head-to-head CSV artifacts into frontend dashboard JSON.",
      "Creates a standalone evidence-oriented HTML report.",
      "Creates a breeder-facing multi-sheet Excel crossing-plan workbook with family screening.",
      "Runs deterministic multi-trait method validation.",
      "Replicated multi-trait parent-size/method grid.",
      "Crop-aware multi-trait validation grid with DH progeny realization.",
      "Runs crop-genome portability stress screens.",
      "Calibrates predicted cross families against realized DH family outcomes.",
      "Runs true autotetraploid 4x benchmark scenarios.",
      "Compares controlled 4x OCS policy configurations.",
      "Runs the broader diploid DH/RIL AlphaSimR benchmark harness."
    ),
    stringsAsFactors = FALSE
  )

  external_integrations <- data.frame(
    id = c("popvar", "simplemating", "alphamate"),
    label = c("PopVar", "SimpleMating", "AlphaMate"),
    status = c("optional_package_or_proxy", "optional_package_or_proxy", "optional_external_executable"),
    backend_function = c("ng_add_popvar_scores", "ng_add_simplemating_scores", "ng_select_alphamate"),
    native_fallback_function = c(
      "ng_run_cross_prediction trait_value_metric='var_complex'; ng_popvar_style_scores",
      "ng_simplemating_* native compatibility helpers",
      "ng_alphamate_style_select"
    ),
    evidence = c(
      "Native var_complex is preferred for user workflows; exact PopVar-style outputs are optional when package support is available.",
      "Exact SimpleMating-style outputs when package support is available; style proxies are labeled.",
      "Official executable integration when configured, with target-degree AlphaMate methods; native alphamate_style is available without the executable."
    ),
    stringsAsFactors = FALSE
  )

  navigation <- data.frame(
    id = c("data", "traits", "scoring", "allocation", "validation", "benchmarks", "polyploid", "external-tools", "reports"),
    label = c("Data QC", "Traits", "Scoring", "Allocation", "Validation", "Benchmarks", "Polyploid", "External Tools", "Reports"),
    href = c("/data", "/traits", "/scoring", "/allocation", "/validation", "/benchmarks", "/polyploid", "/external-tools", "/runs"),
    category = c("readiness", "objective", "evidence", "decision", "validation", "benchmark", "model", "baseline", "report"),
    count_hint = c(nrow(data_qc), 5L, 6L, 4L, 5L, 4L, 3L, nrow(external_integrations), 3L),
    stringsAsFactors = FALSE
  )

  list(
    schema_version = "ng_backend_capabilities.v2",
    generated_at = ng_backend_registry_timestamp(generated_at),
    product_position = list(
      primary_scope = "validated diploid DH/RIL cross design with explicit experimental polyploid modules",
      default_breeding_system = "diploid_dh",
      default_multitrait_method = "auto",
      default_allocation = "ocs_allocation"
    ),
    data_qc = data_qc,
    breeding_systems = breeding_systems,
    method_families = method_families,
    workflows = workflows,
    external_integrations = external_integrations,
    navigation = navigation,
    controls = ng_backend_controls()
  )
}

ng_write_backend_capability_registry_json <- function(output_path = NULL,
                                                      generated_at = Sys.time()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to write backend capability registry JSON")
  }
  if (is.null(output_path) || !nzchar(output_path)) {
    output_path <- file.path("results", "backend_capabilities.json")
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  payload <- ng_backend_capability_registry(generated_at = generated_at)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
