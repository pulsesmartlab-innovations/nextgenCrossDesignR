# Point-and-Click Frontend Developer Guide

This guide is for the developer building a user-facing frontend around the
`nextgenCrossDesign` R package. The frontend should behave like a breeding
decision workbench: users upload data, choose methods and parameters, run the
backend, inspect QC/results, and export a crossing plan without writing R code.

The frontend must not reimplement quantitative genetics logic in JavaScript or
TypeScript. It should collect user choices, call the R backend, and display the
returned objects and files.

Current backend package:

```text
nextgenCrossDesign 0.3.13
dist/nextgenCrossDesign_0.3.13.tar.gz
```

## 1. Backend Installation Check

The frontend backend service should verify the package at startup.

```r
library(nextgenCrossDesign)

required_version <- "0.3.13"
if (utils::packageVersion("nextgenCrossDesign") < required_version) {
  stop(
    "Install nextgenCrossDesign ", required_version,
    " or newer before starting the frontend backend service.",
    call. = FALSE
  )
}

required_args <- c(
  "phenotype_file",
  "genotype_file",
  "map_file",
  "direction_file",
  "phenotype_id_col",
  "genotype_id_col",
  "direction_trait_col",
  "direction_column_col",
  "direction_direction_col",
  "map_marker_col",
  "map_chr_col",
  "map_pos_bp_col",
  "map_position_unit",
  "bp_per_cm",
  "prediction_mode",
  "traits_to_use",
  "index_col",
  "index_direction",
  "trait_value_metric",
  "uc_variance_source",
  "multi_trait_method",
  "trait_weights",
  "threshold_policy",
  "threshold_penalty_weight",
  "threshold_penalty_autoscale",
  "progeny",
  "recombination_model",
  "selection_prop",
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "nIter",
  "burnIn",
  "use_parallel",
  "duplicate_action",
  "duplicate_threshold",
  "duplicate_maf_min",
  "duplicate_max_missing_prop",
  "duplicate_min_compared_markers",
  "n_crosses",
  "max_uses_per_parent",
  "min_unique_parents",
  "max_pair_kinship",
  "optimizer",
  "allocation_method",
  "use_ocs",
  "lambda_group",
  "lambda_mating",
  "lambda_parent_use",
  "lambda_parent_use_mode",
  "local_iter",
  "ocs_iter",
  "alphamate_mode",
  "alphamate_target_degree",
  "alphamate_max_contributions",
  "alphamate_number_of_parents",
  "alphamate_lambda_group",
  "alphamate_lambda_grid",
  "alphamate_executable",
  "alphamate_runtime_path",
  "alphamate_workdir",
  "alphamate_keep_files",
  "alphamate_evol_solutions",
  "alphamate_evol_iterations",
  "alphamate_evol_stop",
  "alphamate_n_threads",
  "priority_breaks",
  "priority_labels",
  "priority_score_weight",
  "priority_kinship_weight",
  "priority_threshold_weight",
  "output_dir",
  "output_file",
  "write_outputs",
  "write_figures",
  "seed"
)

missing_args <- setdiff(
  required_args,
  names(formals(nextgenCrossDesign::ng_run_cross_prediction))
)
if (length(missing_args)) {
  stop(
    "Active package is missing required frontend arguments: ",
    paste(missing_args, collapse = ", "),
    call. = FALSE
  )
}
```

## 2. Backend Capability Registry

Use the package registry to populate frontend menus, capability cards, and
evidence labels.

```r
registry <- nextgenCrossDesign::ng_backend_capability_registry()
nextgenCrossDesign::ng_write_backend_capability_registry_json(
  "results/backend_capabilities.json"
)
```

The registry contains these UI groups:

- `data_qc`
- `breeding_systems`
- `method_families`
- `workflows`
- `external_integrations`
- `navigation`

The frontend should use this registry as the source of truth for labels,
module names, validation status, and exact/proxy/experimental/guarded badges.

## 3. Recommended Frontend Screens

Build the point-and-click workflow in this order.

### Screen 1. Project And File Setup

User provides:

```r
phenotype_file
genotype_file
map_file
direction_file
output_dir
output_file
```

Required UI behavior:

- Show file preview after upload.
- Detect columns from each file.
- Let user map required columns explicitly.
- Do not start analysis until required files and columns are mapped.

### Screen 2. Column Mapping And Input Matching

Expose:

```r
phenotype_id_col
genotype_id_col
direction_trait_col
direction_column_col
direction_direction_col
map_marker_col
map_chr_col
map_pos_bp_col
map_position_unit = "bp"
bp_per_cm = 1e6
```

Important frontend rule:

- `direction_trait_col` is the trait label shown to the user.
- `direction_column_col` is the phenotype column used for analysis.
- `direction_direction_col` is the increase/decrease direction.

The frontend should not manually reorder phenotype, genotype, or map files.
Pass the user-selected column names to the package and display:

```r
result$input_match_audit
```

### Screen 3. Data QC

Expose:

```r
duplicate_action = c("remove", "report", "none")
duplicate_threshold = 0.995
duplicate_maf_min = 0.01
duplicate_max_missing_prop = 0.40
duplicate_min_compared_markers = 100
```

Backend functions:

```r
qc <- ng_preflight_input_tables(
  geno = genotype_table,
  phenotype = phenotype_table,
  trait_spec = trait_direction_table,
  marker_map = marker_map_table,
  ploidy = 2,
  putative_duplicate_check = TRUE,
  putative_duplicate_action = duplicate_action,
  duplicate_threshold = duplicate_threshold,
  duplicate_maf_min = duplicate_maf_min,
  duplicate_max_missing_prop = duplicate_max_missing_prop,
  duplicate_min_compared_markers = duplicate_min_compared_markers
)
```

Display:

```r
qc$status
qc$issues
qc$cleaned_tables_summary
qc$putative_duplicates$pairs
qc$cleaning$putative_duplicates$kept_parents
qc$cleaning$putative_duplicates$removed_parents
```

Optional duplicate figure:

```r
ng_plot_putative_duplicates(
  qc$putative_duplicates,
  output_path = "putative_duplicate_qc.png"
)
```

Frontend blocking rules:

- Block analysis on `qc$status == "blocker"`.
- If duplicate genotypes are detected and `duplicate_action = "report"`, show
  the warning and require the user to decide whether to keep or remove them.
- If `duplicate_action = "remove"`, downstream analysis must use the
  backend-cleaned data returned by the package.

### Screen 4. Prediction Mode

Expose:

```r
prediction_mode = c("trait_by_trait", "index_as_trait")
traits_to_use
index_col
index_direction = c("increase", "decrease")
```

Use `trait_by_trait` when the phenotype file has separate trait columns and a
direction file.

Use `index_as_trait` only when the phenotype file already has a trusted
selection-index column.

### Screen 5. Multi-Trait Objective

Expose:

```r
multi_trait_method = c("auto", "weighted", "economic_index", "desired_gain")
trait_weights
threshold_policy = c("soft", "strict")
threshold_penalty_weight
threshold_penalty_autoscale
```

When to show each method:

- `auto`: default when user only provides trait directions.
- `weighted`: user has relative trait weights.
- `economic_index`: direction file includes economic weights.
- `desired_gain`: direction file includes desired changes and economic weights.

Do not hide trait direction. Disease, lodging, maturity, plant height risk, and
other risk traits must be allowed to decrease.

### Screen 6. Trait Value, Variance, And Usefulness

Expose:

```r
trait_value_metric = c(
  "var_complex",
  "uc",
  "pmv",
  "vpm",
  "var_simple",
  "mean"
)

uc_variance_source = c("pmv", "vpm", "var_simple")
selection_prop
method_varPMV = c("fast", "full_posterior")
```

Recommended default:

```r
trait_value_metric = "var_complex"
uc_variance_source = "pmv"
selection_prop = 0.10
method_varPMV = "fast"
```

User-facing notes:

- `var_complex`: practical default; native PopVar-inspired usefulness.
- `uc`: usefulness criterion using selected variance source.
- `pmv`: marker-effect uncertainty variance.
- `vpm`: recombination variance.
- `var_simple`: relationship-distance diversity proxy.
- `mean`: ranks cross mean only.
- `full_posterior`: best for shortlist reruns, not first broad screens.

### Screen 7. Breeding System And Recombination

Expose:

```r
progeny = c("DH", "RIL")
recombination_model = c("haldane", "kosambi")
ril_mode = "infinite"
assume_inbred = TRUE
```

Recommended default:

```r
progeny = "DH"
recombination_model = "haldane"
```

Use `RIL` only when the program is actually selecting recombinant inbred
families.

### Screen 8. Posterior Prediction

Expose:

```r
run_posterior_prediction = c(TRUE, FALSE)
posterior_method = c("closed_form", "mcmc")
nIter
burnIn
use_parallel
parallel_cores
```

Recommended defaults:

```r
run_posterior_prediction = FALSE
posterior_method = "closed_form"
nIter = 5000
burnIn = 500
use_parallel = FALSE
```

Display posterior outputs when enabled:

```r
result$posterior_predictions
```

For robust posterior allocation, use:

```r
robust_plan <- ng_optimize_robust_mating_plan(
  posterior_scores = result$posterior_predictions[[1L]],
  n_crosses = n_crosses,
  parent_K = ng_parent_kinship(result$cleaned_data$genotype),
  gain_col = "usefulness_pmv_gebv",
  robustness_quantile = 0.25,
  objective = "posterior_quantile",
  max_crosses_per_parent = max_uses_per_parent,
  method = optimizer
)
```

### Screen 9. Number Of Crosses And Portfolio Size

Expose:

```r
n_crosses
max_uses_per_parent
min_unique_parents
max_pair_kinship
```

Optional K-sweep controls:

```r
K_range
criterion = c(
  "elbow_relative",
  "elbow_kneedle",
  "ne_target",
  "coancestry_budget"
)
relative_threshold
ne_min
coancestry_max
```

Backend:

```r
curve <- ng_optimize_mating_plan_curve(
  scores = result$candidate_crosses,
  K_range = K_range,
  gain_col = "multi_trait_score",
  parent_K = ng_parent_kinship(result$cleaned_data$genotype),
  max_crosses_per_parent = max_uses_per_parent,
  method = optimizer,
  criterion = criterion,
  relative_threshold = relative_threshold,
  ne_min = ne_min,
  coancestry_max = coancestry_max
)

recommended_K <- attr(curve, "elbow_K")
```

Optional plot:

```r
plot <- ng_plot_diminishing_returns(curve)
```

### Screen 10. Allocation Method

Expose:

```r
allocation_method = c(
  "ocs",
  "alphamate_style",
  "alphamate_executable"
)

optimizer = c(
  "auto",
  "greedy_local",
  "repair_local",
  "mip_linear",
  "mip_contribution"
)
```

General OCS parameters:

```r
use_ocs
lambda_group
lambda_mating
lambda_parent_use
lambda_parent_use_mode = c("absolute", "adaptive")
local_iter
ocs_iter
```

Recommended default:

```r
allocation_method = "ocs"
optimizer = "auto"
use_ocs = TRUE
lambda_group = 0.05
lambda_mating = 0.02
lambda_parent_use = 0
lambda_parent_use_mode = "absolute"
```

### Screen 11. AlphaMate-Style Controls

Expose when `allocation_method = "alphamate_style"` or
`allocation_method = "alphamate_executable"`:

```r
alphamate_mode = c(
  "ModeOptTarget1",
  "ModeMaxCriterion",
  "ModeMinCoancestry"
)
alphamate_target_degree
alphamate_max_contributions
alphamate_number_of_parents
alphamate_lambda_group
alphamate_lambda_grid
```

Mode meaning:

- `ModeOptTarget1`: balanced gain-diversity frontier.
- `ModeMaxCriterion`: gain-heavy.
- `ModeMinCoancestry`: diversity-heavy.

Suggested target degree:

- `20`: gain-heavy.
- `45`: balanced default.
- `70`: diversity-heavy.

### Screen 12. External AlphaMate Executable

Expose only when `allocation_method = "alphamate_executable"`:

```r
alphamate_executable
alphamate_runtime_path
alphamate_workdir
alphamate_keep_files
alphamate_evol_solutions
alphamate_evol_iterations
alphamate_evol_stop
alphamate_n_threads
```

Frontend behavior:

- Check whether `alphamate_executable` exists.
- If not available, disable external AlphaMate run and suggest
  `allocation_method = "alphamate_style"`.
- Label external AlphaMate as exact executable mode.
- Label `alphamate_style` as native package-developed style mode.

### Screen 13. Priority Ranking

Expose:

```r
priority_breaks = c(0.10, 0.35, 0.70, 1.00)
priority_labels = c(
  "highly_priority",
  "priority",
  "medium_priority",
  "low_priority"
)
priority_score_weight
priority_kinship_weight
priority_threshold_weight
```

Default meaning for 100 selected crosses:

- 10 highly priority crosses.
- 25 priority crosses.
- 35 medium priority crosses.
- 30 low priority crosses.

Display:

```r
result$selected_crosses
table(result$selected_crosses$priority_tier)
```

### Screen 14. Outputs And Reports

Expose:

```r
write_outputs = c(TRUE, FALSE)
write_figures = c(TRUE, FALSE)
output_dir
output_file
```

Display:

```r
result$output_files
```

Expected files:

- `priority_score_vs_kinship.png`
- crossing-plan workbook when `write_outputs = TRUE` and `openxlsx` is
  installed

The frontend should show the figure and provide download links for workbook,
figures, logs, and run settings.

## 4. Main Backend Runner Template

> **Preferred: use the headless JSON runner.** Instead of hand-building the R call, spawn
> `tools/run_cross_prediction_json.R` with a JSON config (keys mirror `ng_run_cross_prediction`;
> unknown keys are a hard error) and read back an `ng_run_result.v1` JSON. See
> `docs/frontend/README.md` and `docs/frontend/contracts/` (`config_schema.json`, `example_config.json`,
> `example_result.json`). The template below shows the underlying call for reference; it also accepts
> the marker-effect **training set** (`training_genotype_file` / `training_phenotype_file` /
> `training_genotype_id_col` / `training_phenotype_id_col`) and `grm_method` = `vanraden` | `yang`.

The frontend should generate a run object and pass it into a backend R script
similar to this.

```r
library(nextgenCrossDesign)

run_cross_prediction_from_frontend <- function(params) {
  result <- ng_run_cross_prediction(
    phenotype_file = params$phenotype_file,
    genotype_file = params$genotype_file,
    map_file = params$map_file,
    direction_file = params$direction_file,

    phenotype_id_col = params$phenotype_id_col,
    genotype_id_col = params$genotype_id_col,
    direction_trait_col = params$direction_trait_col,
    direction_column_col = params$direction_column_col,
    direction_direction_col = params$direction_direction_col,
    map_marker_col = params$map_marker_col,
    map_chr_col = params$map_chr_col,
    map_pos_bp_col = params$map_pos_bp_col,
    map_position_unit = params$map_position_unit,
    bp_per_cm = params$bp_per_cm,

    prediction_mode = params$prediction_mode,
    traits_to_use = params$traits_to_use,
    index_col = params$index_col,
    index_direction = params$index_direction,

    trait_value_metric = params$trait_value_metric,
    uc_variance_source = params$uc_variance_source,
    multi_trait_method = params$multi_trait_method,
    trait_weights = params$trait_weights,
    threshold_policy = params$threshold_policy,
    threshold_penalty_weight = params$threshold_penalty_weight,
    threshold_penalty_autoscale = params$threshold_penalty_autoscale,

    progeny = params$progeny,
    recombination_model = params$recombination_model,
    selection_prop = params$selection_prop,
    method_varPMV = params$method_varPMV,
    ril_mode = params$ril_mode,
    run_posterior_prediction = params$run_posterior_prediction,
    posterior_method = params$posterior_method,
    nIter = params$nIter,
    burnIn = params$burnIn,
    use_parallel = params$use_parallel,
    parallel_cores = params$parallel_cores,

    duplicate_action = params$duplicate_action,
    duplicate_threshold = params$duplicate_threshold,
    duplicate_maf_min = params$duplicate_maf_min,
    duplicate_max_missing_prop = params$duplicate_max_missing_prop,
    duplicate_min_compared_markers = params$duplicate_min_compared_markers,

    n_crosses = params$n_crosses,
    max_uses_per_parent = params$max_uses_per_parent,
    min_unique_parents = params$min_unique_parents,
    max_pair_kinship = params$max_pair_kinship,
    optimizer = params$optimizer,
    allocation_method = params$allocation_method,
    use_ocs = params$use_ocs,
    lambda_group = params$lambda_group,
    lambda_mating = params$lambda_mating,
    lambda_parent_use = params$lambda_parent_use,
    lambda_parent_use_mode = params$lambda_parent_use_mode,
    local_iter = params$local_iter,
    ocs_iter = params$ocs_iter,

    alphamate_mode = params$alphamate_mode,
    alphamate_target_degree = params$alphamate_target_degree,
    alphamate_max_contributions = params$alphamate_max_contributions,
    alphamate_number_of_parents = params$alphamate_number_of_parents,
    alphamate_lambda_group = params$alphamate_lambda_group,
    alphamate_lambda_grid = params$alphamate_lambda_grid,
    alphamate_executable = params$alphamate_executable,
    alphamate_runtime_path = params$alphamate_runtime_path,
    alphamate_workdir = params$alphamate_workdir,
    alphamate_keep_files = params$alphamate_keep_files,
    alphamate_evol_solutions = params$alphamate_evol_solutions,
    alphamate_evol_iterations = params$alphamate_evol_iterations,
    alphamate_evol_stop = params$alphamate_evol_stop,
    alphamate_n_threads = params$alphamate_n_threads,

    priority_breaks = params$priority_breaks,
    priority_labels = params$priority_labels,
    priority_score_weight = params$priority_score_weight,
    priority_kinship_weight = params$priority_kinship_weight,
    priority_threshold_weight = params$priority_threshold_weight,

    output_dir = params$output_dir,
    output_file = params$output_file,
    write_outputs = params$write_outputs,
    write_figures = params$write_figures,
    assume_inbred = params$assume_inbred,
    use_cpp = TRUE,
    seed = params$seed
  )

  stopifnot(inherits(result, "ng_cross_prediction_result"))
  result
}
```

## 5. Frontend Result Contract

After a run completes, display these objects.

```r
result$prediction_mode
result$qc
result$cleaned_data
result$trait_direction
result$input_match_audit
result$effect_summary
result$marker_effects
result$trait_scores
result$posterior_effects
result$posterior_predictions
result$candidate_crosses
result$selected_crosses
result$objective
result$plan_summary
result$output_files
result$settings
```

Minimum user-facing result tabs:

- QC issues and duplicate-removal audit.
- Input matching audit.
- Marker-effect summary.
- Candidate cross scores.
- Selected priority-ranked crosses.
- Parent-use summary.
- Method/settings summary.
- Figures.
- Workbook/downloads.

**New result fields (this release).** `input_match_audit$training_only_count`,
`$effect_training_n`, and `$training_ids` distinguish the candidate parents (`matched_parent_count`,
which are crossed) from the effect-only training individuals (never crossed) — show "trained on N,
crossing K parents". `effect_summary$marker_effect_training_n` reports how many individuals fit each
trait's effects. The headless runner (`tools/run_cross_prediction_json.R`) serializes exactly this
result as `ng_run_result.v1` — code your views against `docs/frontend/contracts/example_result.json`.

## 6. Required Example Coverage

The frontend developer should use the installed examples as acceptance tests.
Every capability below should have a corresponding point-and-click path.

| Example | Frontend capability covered |
| --- | --- |
| `00_user_friendly_cross_prediction.R` | Main user workflow |
| `01_chronological_cross_prediction_workflow.R` | Full workflow order and lower-level functions |
| `02_minimal_reproducible_user_run.R` | Minimal smoke test |
| `03_variance_method_comparison.R` | Trait value metrics |
| `04_alphamate_executable_user_run.R` | External AlphaMate executable |
| `05_ocs_user_run.R` | OCS allocation |
| `06_optimizer_parameter_guide.R` | Optimizer choices |
| `07_trait_value_metric_parameter_guide.R` | PMV, posterior, RIL, parallel settings |
| `08_prediction_mode_trait_by_trait.R` | Trait-by-trait prediction |
| `09_prediction_mode_index_as_trait.R` | Index-as-trait prediction |
| `10_multitrait_method_auto.R` | Auto multi-trait objective |
| `11_multitrait_method_weighted.R` | User-weighted multi-trait objective |
| `12_multitrait_method_economic_index.R` | Economic-index objective |
| `13_multitrait_method_desired_gain.R` | Desired-gain objective |
| `14_qc_duplicate_removal_and_reporting.R` | Duplicate QC report/remove and duplicate figure |
| `15_input_matching_and_map_units.R` | ID matching, marker matching, BP map units |
| `16_cross_number_sweep_diminishing_returns.R` | Cross-number sweep |
| `17_workbook_figures_and_priority_outputs.R` | Workbook and figures |
| `18_allocation_method_comparison.R` | OCS, AlphaMate-style, external AlphaMate comparison |
| `19_posterior_robust_mating_plan.R` | Posterior robust mating plan |
| `20_full_posterior_pmv_shortlist.R` | Full-posterior PMV shortlist |
| `21_gain_diversity_balance_and_target_coancestry.R` | Strategy dial / diversity_emphasis / target_coancestry + frontier export |
| `22_progeny_inbreeding_management.R` | Progeny-inbreeding penalty + histogram export |
| `23_breeder_mating_constraints.R` | min_crosses_per_parent, committed_crosses, parent_group/permission/quota |
| `24_marker_steering_and_lethal_guarding.R` | marker_target_spec / lambda_marker, lethal_spec |
| `25_cost_and_logistics_factors.R` | cross_cost / cost_col / budget / logistic_col |
| `26_evolution_optimizer.R` | Native evolution optimizer |
| `27_breeder_controls_combined.R` | Several breeder controls in one run |
| `28_polyploid_mate_design.R` | Any-ploidy one-call mate design |
| `29_polyploid_qc_and_grm.R` | Polyploid QC + VanRaden/Yang GRM |
| `30_polyploid_additive_dominance_effects.R` | Additive+dominance genomic prediction |
| `31_polyploid_dominance_crossing.R` | Dominance-aware crossing (heterosis) |
| `32_exact_cross_trait_covariance.R` | Exact within-family cross-trait covariance → multi-trait threshold probability |
| `33_marker_effect_training_set.R` | Marker-effect training-set augmentation (extra effect-only individuals) |

## 7. Frontend Defaults

Use these defaults for the first production UI.

```r
prediction_mode <- "trait_by_trait"
multi_trait_method <- "auto"
trait_weights <- NULL
threshold_policy <- "soft"

trait_value_metric <- "var_complex"
uc_variance_source <- "pmv"
selection_prop <- 0.10
method_varPMV <- "fast"

progeny <- "DH"
recombination_model <- "haldane"
ril_mode <- "infinite"
assume_inbred <- TRUE

duplicate_action <- "remove"
duplicate_threshold <- 0.995
duplicate_maf_min <- 0.01
duplicate_max_missing_prop <- 0.40
duplicate_min_compared_markers <- 100

n_crosses <- 100
max_uses_per_parent <- 6
min_unique_parents <- NULL
max_pair_kinship <- Inf

allocation_method <- "ocs"
optimizer <- "auto"
use_ocs <- TRUE
lambda_group <- 0.05
lambda_mating <- 0.02
lambda_parent_use <- 0
lambda_parent_use_mode <- "absolute"
local_iter <- 2000
ocs_iter <- 5

alphamate_mode <- "ModeOptTarget1"
alphamate_target_degree <- 45
alphamate_max_contributions <- max_uses_per_parent

priority_breaks <- c(0.10, 0.35, 0.70, 1.00)
priority_labels <- c(
  "highly_priority",
  "priority",
  "medium_priority",
  "low_priority"
)
priority_score_weight <- 1.0
priority_kinship_weight <- 0.15
priority_threshold_weight <- 1.0

run_posterior_prediction <- FALSE
posterior_method <- "closed_form"
nIter <- 5000
burnIn <- 500
use_parallel <- FALSE

write_outputs <- TRUE
write_figures <- TRUE
```

## 8. Validation And Guardrails

The frontend is acceptable when all these conditions are true:

- Data QC runs before analysis.
- Blocker QC issues prevent run launch.
- Duplicate genotype report/remove is available inside the workflow.
- User can map phenotype ID, genotype ID, direction columns, and marker-map
  columns.
- Map positions support base-pair input with `map_position_unit = "bp"`.
- Trait directions are visible and required.
- `trait_by_trait` and `index_as_trait` are both available.
- `auto`, `weighted`, `economic_index`, and `desired_gain` are all available.
- `var_complex`, `uc`, `pmv`, `vpm`, `var_simple`, and `mean` are all
  available.
- `ocs`, `alphamate_style`, and guarded `alphamate_executable` are available.
- OCS parameters and AlphaMate-style parameters are separated in the UI.
- Cross-number sweep is available.
- Posterior prediction and robust posterior mating plan are available as
  advanced controls.
- Full-posterior PMV is available as a shortlist rerun.
- Priority tiers are displayed and editable.
- Workbook and figures are generated and downloadable.
- Experimental polyploid and exact external-tool features are labeled by
  evidence status from the backend registry.

## 9. Do Not Reimplement These In The Frontend

Do not reimplement:

- duplicate genotype detection or removal;
- marker map alignment;
- phenotype/genotype ID alignment;
- marker-effect estimation;
- cross scoring;
- multi-trait scoring;
- PMV, VPM, usefulness, or posterior prediction;
- OCS optimization;
- AlphaMate-style allocation;
- external AlphaMate execution;
- priority ranking;
- workbook writing;
- figure writing.

The frontend should be an excellent user interface over the package, not a
second genetics engine.
