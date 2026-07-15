# nextgenCrossDesign Backend User Guide

This guide explains how to use the `nextgenCrossDesign` R backend, what each
major function family does, and which defaults are recommended in common
breeding-program scenarios.

The backend is designed as a quantitative genetics package first. The frontend
can call it later, but every core workflow here can be run directly from R or
from the command-line wrappers in `tools/`.

## Current Backend Scope

The package supports:

- data preflight, duplicate checking, and putative duplicate genotype removal;
- native LD pruning for redundant marker removal before formal analysis;
- diploid DH/RIL cross scoring;
- single-trait and multi-trait cross selection;
- positive and negative trait directions;
- unknown-weight, weighted, economic-index, threshold, and desired-gain
  multi-trait objectives;
- mate allocation with gain, diversity, parent-use, and relationship control;
- native `var_complex`, SimpleMating-style, and AlphaMate-style components
  without requiring those tools;
- optional exact PopVar, SimpleMating, AlphaMate, and genomicMateSelectR
  integrations when available;
- crop and parent-size validation runners;
- experimental autotetraploid 4x and allopolyploid subgenome workflows;
- JSON/HTML/Excel outputs for frontend and breeder-facing reports.

The strongest validated production path is still **diploid DH/RIL cross design**.
Polyploid modules are implemented and tested, but should be treated as
experimental until validated for the target crop, ploidy model, marker density,
and breeding scheme.

## Installation And Loading

From the repository root:

```r
source("R/load.R")
ng_load(use_cpp = TRUE)
```

For package build verification:

```powershell
Rscript tools\check_r_package.R
```

The package check is staged because this repository also contains frontend code.
Use `tools/check_r_package.R` as the backend package gate.


## Recommended Defaults By Scenario

| Scenario | Recommended default | Why |
| --- | --- | --- |
| First-time data import | `ng_preflight_input_tables()` | Catch duplicate IDs, duplicate crosses, putative duplicate genotype profiles, marker-map issues, dosage errors, and relationship-matrix mismatches before scoring. |
| Breeder has phenotype, genotype, map, and trait-direction files | `ng_run_cross_prediction(prediction_mode = "trait_by_trait")` | Runs the practical end-to-end workflow: QC, duplicate removal, per-trait marker effects, cross prediction, multi-trait objective, OCS allocation, priority tiers, and optional workbook/figures. |
| High-density marker set before modeling | `ng_ld_prune_geno()`, `ng_design_crosses(ld_pruning = TRUE)`, or `ng_run_cross_prediction(ld_pruning = TRUE)` | Removes low-MAF and high-LD redundant markers with a native C++ backend for large matrices and an R fallback. |
| Diploid DH selection, single trait | `ng_design_crosses()` or `ng_score_crosses()` plus `ng_optimize_mating_plan()` | Uses ridge marker effects, cross means, usefulness, PMV-style variance, kinship, and constrained allocation. Choose the map function with `recomb_model = "haldane"` (default) or `"kosambi"`. |
| Diploid RIL selection | `ng_design_crosses(target = "RIL")`, or `ng_score_crosses(target = "RIL")` / `ng_run_cross_prediction(progeny = "RIL")` | Same API family; the target must match the progeny type so the within-family recombination variance is correct. |
| Choose the relationship matrix (GRM) | `grm_method = "vanraden"` (default) or `"yang"` on `ng_design_crosses()`, `ng_run_cross_prediction()`, `ng_score_crosses()`, and `ng_design_crosses_poly()` | VanRaden uses one overall allele-frequency scaling; Yang/GCTA standardizes each marker to unit variance. |
| Steer allocation toward a target allele | `marker_target_spec = ng_marker_target_spec(...)` with `lambda_marker` on `ng_design_crosses()` or `ng_run_cross_prediction()` | Blends a marker-target score into the merit the allocator optimizes; pair with `lethal_spec` to also drop carrier x carrier crosses. |
| Breeder has several traits but no reliable weights | `ng_add_multitrait_score(method = "auto")` | Rank-normalized default with soft threshold penalties. Safer than forcing arbitrary weights. |
| Breeder has relative economic weights | `method = "economic_index"` | Uses oriented/scaled trait covariance to build a covariance-aware economic index. |
| Breeder has desired response targets | `method = "desired_gain"` | Uses desired changes, economic weights, and stabilized covariance to estimate index coefficients. |
| Some traits must increase and others decrease | Always define `direction` in `ng_multitrait_spec()` | Prevents disease, lodging, height, or risk traits from being optimized in the wrong direction. |
| Must meet hard market/disease cutoffs | Start with soft thresholds; use `strict_thresholds = TRUE` only when necessary | Soft thresholds avoid empty feasible sets; strict thresholds are best for true non-negotiable cutoffs. |
| Need diversity and parent-use control | `ng_optimize_mating_plan(method = "auto", lambda_group > 0, lambda_parent_use > 0)` | Balances gain, group coancestry, and parent contribution. |
| Need a practical crossing list | Request about 100 crosses, then use `ng_rank_cross_priority()` | Produces breeder-friendly tiers: `highly_priority`, `priority`, `medium_priority`, and `low_priority`. |
| Want AlphaMate-like allocation without AlphaMate | `ng_alphamate_style_select()` or `ng_run_cross_prediction(allocation_method = "alphamate_style")` | Native target-degree gain-diversity frontier selection. |
| Want SimpleMating-like parent/cross handling | `ng_simplemating_relate_thinning()`, `ng_simplemating_build_crosses()`, `ng_simplemating_past_thinning()`, `ng_simplemating_select_crosses_native()` | Native thinning, past-cross removal, and constrained selection. |
| Want PopVar-like usefulness without PopVar | `ng_run_cross_prediction(trait_value_metric = "var_complex")` | Native mean plus complex within-family variance/usefulness scoring; no PopVar dependency. |
| Need breeder-facing Excel output | `ng_write_crossing_plan_workbook()` | Produces top parents, recommended crosses, criteria, family screening, and a chart. |
| Need frontend capability metadata | `ng_backend_capability_registry()` or `tools/export_backend_capabilities_json.R` | Gives the UI the backend method/navigation contract. |
| True autotetraploid 4x experiment | `ng_poly4x_*` functions and runners | Use as experimental validation or research workflow, not a universal claim. |
| Wheat-like allopolyploid/subgenome experiment | `ng_poly_subgenome_*` functions | Use when subgenomes can be represented as disomic dosage matrices. |

## Data Inputs

### Genotype Matrix

Use parents/lines as rows and markers as columns. Row names must be parent IDs.

```r
geno <- as.matrix(read.csv("geno.csv", row.names = 1, check.names = FALSE))
```

For diploid inbred DH/RIL workflows, marker dosage is usually `0`, `1`, or `2`.
Some exact PopVar calls require inbred `0/2` genotypes; native
`nextgenCrossDesign` scoring can handle broader numeric dosage inputs.

### Phenotype Or Adjusted Values

For direct marker-effect fitting:

```r
y <- read.csv("phenotypes.csv")
y_vec <- setNames(y$yield_blup, y$parent_id)
```

For production use, prefer adjusted phenotypes, BLUEs, or BLUPs from your trial
analysis. Raw plot-level data should be adjusted before cross design unless the
analysis model is part of your workflow.

### Marker Map

Recommended columns for user files:

```text
marker, chr, pos_bp
```

`marker` must match genotype marker column names. User-facing workflows should
provide base-pair positions for consistency, for example `Position_BP`; pass
that column as `map_pos_bp_col` with `map_position_unit = "bp"`. The package
keeps `pos_bp` in the cleaned marker map and derives `pos_cm` internally using
`bp_per_cm` because recombination-aware variance functions use cM distances.

### Candidate Crosses

Candidate-pair tables should contain:

```text
parent1, parent2
```

Self-crosses, exact duplicate crosses, and reciprocal duplicates are flagged by
the preflight checker. Near-identical genotype profiles with different sample
IDs can also be screened and removed inside package QC before scoring.

### Trait Specification

Multi-trait specs can be built directly:

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "lodging"),
  column = c("pred_yield", "pred_disease", "pred_lodging"),
  direction = c("maximize", "minimize", "minimize"),
  weight = c(1, 1, 1)
)
```

Or supplied as a data frame with columns:

```text
trait, column, direction, weight, min_value, max_value,
threshold_weight, desired_change, economic_weight
```

For the file-based user workflow, the direction file may use breeder-friendly
column names. Pass those names explicitly:

```r
direction_trait_col <- "Trait"                  # name shown in reports
direction_column_col <- "PhenotypeColumn"       # actual phenotype column
direction_direction_col <- "Selection_direction" # increase/decrease
```

Use `direction_column_col = NULL` when the trait names in the direction file are
also the phenotype column names.

## Preflight Data QC

Run preflight before scoring or allocation:

```r
qc <- ng_preflight_input_tables(
  geno = geno,
  phenotype = phenotype,
  candidate_pairs = candidate_pairs,
  trait_spec = traits,
  marker_map = marker_map,
  parent_K = parent_K,
  ploidy = 2,
  putative_duplicate_check = TRUE,
  putative_duplicate_action = "remove",
  duplicate_threshold = 0.995
)

qc$status
qc$issues
qc$cleaned_tables_summary
clean_geno <- qc$cleaned_tables$geno
clean_pheno <- qc$cleaned_tables$phenotype
```

When `putative_duplicate_check = TRUE`, the package computes an IBS 0-1
similarity screen after marker filtering and reports near-identical genotype
profiles even when the sample IDs are different. With
`putative_duplicate_action = "remove"`, QC keeps one parent per duplicate
cluster, removes redundant duplicate parents from genotype and phenotype input
tables, filters candidate pairs and relationship matrices when supplied, and
returns the cleaned inputs in `qc$cleaned_tables`. The removal audit is stored
in `qc$cleaning$putative_duplicates`, including kept parents, removed parents,
and row counts removed from each table.

Use `ng_detect_putative_duplicates()` directly when you need the full
duplicate-pair table, nearest-neighbor table, or similarity matrix for
reporting. Use `ng_plot_putative_duplicates()` to save a duplicate heatmap for
QC review.

Write JSON for a UI or batch log:

```r
ng_write_data_preflight_json(
  output_path = "results/data_preflight.json",
  geno = geno,
  phenotype = phenotype,
  candidate_pairs = candidate_pairs,
  trait_spec = traits,
  marker_map = marker_map,
  putative_duplicate_check = TRUE,
  putative_duplicate_action = "remove"
)
```

The JSON report includes issue summaries, duplicate-pair metadata, cleaned-table
summaries, and duplicate-removal audit tables. It intentionally omits full
cleaned genotype/phenotype matrices to keep UI and batch logs compact. Use
`ng_preflight_input_tables()` directly when downstream R analysis needs the
cleaned tables.

Command-line wrapper:

```powershell
$env:NG_PREFLIGHT_GENO="geno.csv"
$env:NG_PREFLIGHT_PHENOTYPE="phenotype.csv"
$env:NG_PREFLIGHT_TRAIT_SPEC="traits.csv"
$env:NG_PREFLIGHT_MARKER_MAP="marker_map.csv"
$env:NG_PREFLIGHT_PUTATIVE_DUPLICATES="true"
$env:NG_PREFLIGHT_DUPLICATE_ACTION="remove"
$env:NG_PREFLIGHT_OUT="results/data_preflight.json"
Rscript tools\run_data_preflight.R
```

Treat `status = "blocker"` as a stop condition. Fix IDs, duplicates, dosage,
or marker-map problems before running selection. If putative duplicate genotype
profiles are expected, prefer `putative_duplicate_action = "remove"` and run
scoring on `qc$cleaned_tables`; if both records are intentionally kept, document
that decision in the QC report.

## LD Pruning Before Formal Analysis

Use LD pruning when marker density is high or when redundant markers are likely
to overweight local chromosome regions. The native backend is loaded by
`ng_load(use_cpp = TRUE)` and falls back to R when C++ is unavailable.

```r
geno_pruned <- ng_ld_prune_geno(
  geno,
  window = 100,
  r2_threshold = 0.9,
  maf_threshold = 0.01,
  ploidy = 2,
  backend = "auto"
)

attr(geno_pruned, "ld_pruning_report")
```

For the full one-call workflow:

```r
result <- ng_design_crosses(
  geno = geno,
  y = y_vec,
  marker_map = marker_map,
  ld_pruning = TRUE,
  ld_backend = "auto"
)

result$ld_pruning_report
```

The algorithm first drops markers below `maf_threshold`, then scans markers
within `window`, groups marker pairs with dosage-correlation `r2` above
`r2_threshold`, and keeps one representative marker per LD group. The default
representative is the marker with the highest MAF, with original marker order as
the tie-breaker.

## User-Friendly Diploid DH/RIL Workflow

The recommended first workflow for a breeder is `ng_run_cross_prediction()`.
Use `prediction_mode = "trait_by_trait"` when the user provides multiple
traits and a direction file. Marker effects are estimated separately for every
trait, and the multi-trait objective is built after cross prediction.

```r
result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  phenotype_id_col = "NAME",
  genotype_id_col = "NAME",
  direction_trait_col = "Trait",
  direction_column_col = NULL,
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
  selection_prop = 0.10,
  progeny = "RILs",
  duplicate_action = "remove",
  n_crosses = 100,
  max_uses_per_parent = 20,
  optimizer = "auto",
  allocation_method = "alphamate_style",
  use_ocs = TRUE,
  alphamate_mode = "ModeOptTarget1",
  alphamate_target_degree = 45,
  alphamate_max_contributions = 20
)

result$input_match_audit
```

`ng_run_cross_prediction()` also surfaces the same scoring and mate-design controls
as the lower-level API, so the end-to-end runner reaches them without dropping to
`ng_score_crosses()` / `ng_design_crosses()`:

- `grm_method = "vanraden"` (default) or `"yang"` selects the relationship matrix.
- `recombination_model = "haldane"` (default) or `"kosambi"` sets the map function.
- `marker_target_spec = ng_marker_target_spec(...)` with `lambda_marker` steers
  allocation toward a target allele (the reported `candidate_crosses` gains a
  `marker_target_score` column); combine with `lethal_spec` to drop carrier crosses.
- `ld_pruning = TRUE` (with `ld_r2_threshold`, `ld_maf_threshold`, ...) prunes
  redundant markers before scoring and returns `result$ld_pruning_report`.

`trait_value_metric = "var_complex"` is the native PopVar-inspired option for
multi-trait user runs. It does not call PopVar. It uses oriented cross mean
plus or minus selection intensity times within-family standard deviation,
preferring PMV variance when available and falling back to recombination
variance or `var_simple`.

### Choosing A Trait-Value Metric

The metrics differ in how they score a cross:

- `mean` -- mid-parent GEBV only (expected progeny mean, no within-family variance).
- `vpm` / `uc` (with `uc_variance_source = "vpm"`) -- usefulness = mean +/- i * SD using
  the recombination-aware within-family variance.
- `pmv` / `var_complex` (default) -- usefulness using PMV, which adds marker-effect
  estimation uncertainty to the recombination variance.
- `var_simple` -- a relationship-DISTANCE metric (how unrelated the two parents are). It is a
  diversity/QC proxy, NOT a cross-merit predictor.

Evidence-based guidance (from the package's own simulation studies,
`tools/run_metric_merit_study.R` single-generation and
`tools/run_ril_breeding_program_benchmark.R` recurrent; config-scoped, directional):

- Default: **`var_complex`** (PMV-based usefulness). It is never much worse than the best
  and is clearly best for oligogenic traits with adequate training and heritability. `vpm`
  and `pmv` rank crosses almost identically to it, so treat them as interchangeable choices.
- Use **`mean`** when the trait is highly polygenic OR your training set is small / low
  heritability -- there the within-family variance term adds little and can add noise, so the
  cheaper mean-only metric is competitive.
- Do NOT select on **`var_simple`** as a merit metric -- it ranks crosses poorly
  single-generation. (Interestingly it does well over MANY cycles of recurrent selection,
  because favouring unrelated parents preserves genetic variance -- but that is a diversity
  effect, not merit.)
- For LONG-TERM recurrent selection (recycling parents over many cycles), pure usefulness
  metrics exhaust genetic variance and gain plateaus. Manage diversity EXPLICITLY instead --
  set a coancestry penalty (`lambda_group`), or use the strategy dial
  (`strategy` / `diversity_emphasis`) or an inbreeding target (`target_coancestry`), on top
  of a good merit metric like `var_complex`. See
  `inst/examples/21_gain_diversity_balance_and_target_coancestry.R`.

Prediction accuracy is similar across all these metrics; the differences are about how they
model within-family variance and manage diversity, not about prediction quality.

`allocation_method` controls the mating-plan selection layer:

```text
allocation_method = "ocs"                  package OCS-style optimizer
allocation_method = "alphamate_style"      native package-developed AlphaMate-style allocation
allocation_method = "alphamate_executable" call the external AlphaMate executable
```

### Choosing An Optimizer

Within `allocation_method = "ocs"`, the `optimizer` argument picks the engine that solves the
mate-allocation problem (all optimize the same OCS objective):

| optimizer | What it is | Speed | Solution quality |
|---|---|---|---|
| `auto` (default) | dispatches to `mip_contribution` when a diversity penalty is active and `lpSolve` is available, else `mip_linear` (no penalty) or `greedy_local`; falls back to `greedy_local` on oversized/timed-out MIP problems (`summary$mip_fallback`). | adaptive | best available |
| `mip_contribution` | exact-style MIP OCS with the contribution (coancestry) penalty | slowest; size/time-guarded | best achieved objective |
| `mip_linear` | MIP without the coancestry penalty (use when `lambda_group = 0`) | slow | exact for the gain-only problem |
| `greedy_local` | fast greedy + local swap | fastest | slightly below MIP |
| `repair_local` | top-N + capacity repair + local swap | fast | near greedy |
| `evolution` | native memetic genetic algorithm (aliases `ga`/`de`/`memetic`); warm-started from greedy with elitism, so never worse than greedy | moderate (`evol_*` tunable) | competitive with MIP |

Evidence (from the package's own studies): on a fixed-input optimizer bake-off
(`tools/run_optimizer_benchmark.R`), `mip_contribution` ranks first on achieved objective and
`evolution` is a close second, both beating `greedy_local`/`repair_local`; in the recurrent RIL
study `evolution` is within noise of OCS on realized gain.

Recommendation:

- **Leave `optimizer = "auto"` (recommended).** It already routes to the exact MIP path when
  that is best and degrades gracefully to `greedy_local` when the MIP would be too large or
  slow, so you get the best available solution without hanging.
- For a **large candidate set** where MIP is slow (or you want a better solution than greedy
  without the MIP cost), use **`"evolution"`** (tune with `evol_solutions` / `evol_iterations` /
  `evol_stop`).
- Use **`"greedy_local"`** when you need the fastest possible run or a deterministic diagnostic.
- The optimizer choice affects the *allocation* only; it does not change prediction accuracy.

Use `prediction_mode = "index_as_trait"` only when the phenotype table already
contains a trusted selection-index column. In that case marker effects are
estimated once for that index.

### Polyploid (Any-Ploidy) Mate Design In One Call

The mate allocator (`ng_optimize_mating_plan`) and all its controls operate on a candidate-cross
table + a kinship matrix and are **ploidy-agnostic**. So a polyploid breeder gets the whole control
suite — the usefulness/mean criteria, the strategy dial, `target_coancestry`, committed matings,
group permission/quotas, cost/budget/logistics, and the evolution optimizer — from a single entry,
`ng_design_crosses_poly()` (see `inst/examples/28_polyploid_mate_design.R`):

```r
plan <- ng_design_crosses_poly(
  dosage,                      # parents x markers, allele dosage 0..ploidy
  n_crosses = 15L, ploidy = 4L,
  phenotype = pheno,           # or effects = <per-marker additive effects>
  run_qc = TRUE,               # ploidy-aware QC (0..ploidy, missingness, MAF, monomorphic)
  strategy = "balanced",       # any native control forwards: target_coancestry, committed_crosses,
  committed_crosses = fixed)   #   group_permission/quota, cost_col/budget, method = "evolution", ...
```

It runs ploidy-aware QC (`ng_polyploid_qc`), scores analytically (correct allele-frequency-based
polyploid GRM, `ng_polyploid_grm` — VanRaden **or** Yang via `grm_method =`, *not* a diploid-style
shortcut), then runs `ng_optimize_mating_plan()` with everything you pass through `...`.
`ng_poly_score_crosses()` is the lower-level scorer, and `ng_polyploid_grm()` / `ng_polyploid_qc()`
are usable standalone. Selectable metric: `gain = "mean"` (mid-parent breeding value) or
`gain = "usefulness"` (mean + i·within-family SD; the additive segregation variance comes from the
progeny-moment table). For clonal/heterosis crops set `dominance = TRUE` to add mid-parent heterosis
and dominance segregation variance (and `double_reduction` for autopolyploids); additive-only is the
default.

## Basic Lower-Level Diploid DH/RIL Workflow

The shortest full workflow is `ng_design_crosses()`:

```r
source("R/load.R")
ng_load(use_cpp = TRUE)

result <- ng_design_crosses(
  geno = geno,
  y = y_vec,
  marker_map = marker_map,
  n_crosses = 100,
  max_crosses_per_parent = 6,
  selection_prop = 0.10,
  lambda_group = 0.05,
  lambda_mating = 0
)

scores <- result$scores
plan <- result$plan
summary <- result$plan_summary
```

Use this when you have one main target trait and want a complete default path:
fit effects, score all crosses, calculate parent kinship, and allocate a mating
plan.

For more control, run the steps explicitly:

```r
effects <- ng_fit_ridge_effects(geno, y_vec, ids = rownames(geno), seed = 1)

scores <- ng_score_crosses(
  geno = geno,
  effects = effects,
  marker_map = marker_map,
  ids = rownames(geno),
  adjusted_pheno = y_vec,
  selection_prop = 0.10,
  use_cpp = TRUE
)

parent_K <- ng_parent_kinship(geno)

plan <- ng_optimize_mating_plan(
  scores = scores,
  n_crosses = 100,
  gain_col = "uc_gated",
  parent_K = parent_K,
  max_crosses_per_parent = 6,
  lambda_group = 0.05,
  lambda_parent_use = 1.0,
  lambda_parent_use_mode = "adaptive",
  method = "auto"
)
```

Recommended single-trait defaults:

- `gain_col = "uc_gated"` for a conservative usefulness default.
- `selection_prop = 0.10` unless the program routinely selects a different
  top-family fraction.
- `max_crosses_per_parent = 3` to `5` for normal breeding use.
- `lambda_group = 0.03` to `0.10` when relatedness control matters.
- `lambda_parent_use_mode = "adaptive"` when parent-use penalty should scale
  with score spread and number of crosses.

## Multi-Trait Selection

### Choosing A Multi-Trait Method

`multi_trait_method` (in `ng_run_cross_prediction`, or `method` in
`ng_add_multitrait_score` / `ng_multitrait_spec`) selects how per-trait predictions are
combined into one selection value. IMPORTANT: unlike the trait-value metric and the optimizer
(which the package benchmarks head-to-head, because each answers the *same* question), the
multi-trait methods optimize DIFFERENT objectives -- balanced progress vs weighted importance
vs economic merit vs target gains vs constraints. There is no single "best"; the choice is
dictated by **what information you have and your breeding goal**, not by a performance race.

Pick by your inputs:

| You have… | Use | Why |
|---|---|---|
| only trait directions (+ maybe rough min/max thresholds) | **`auto`** (recommended default) | rank-normalizes each oriented trait and uses equal normalized weights; robust when weights/units are unknown, hard to misuse |
| declared relative importance weights per trait | **`weighted`** | applies your weights to rank-normalized, direction-oriented traits |
| reliable economic values AND a trustworthy genetic variance–covariance matrix | **`economic_index`** (Smith–Hazel) | maximizes aggregate economic merit `b = P^-1 a` |
| target genetic changes per trait (a desired-gain vector) | **`desired_gain`** (Pesek–Baker) | solves for the index that delivers the requested change profile |
| hard/soft minimum or maximum constraints per trait | **`threshold`** | keeps crosses within bounds (soft penalty by default); can combine with the above |

Recommendation: **default to `auto`.** Economic-index is deliberately NOT the default because
reliable economic weights are hard to obtain, and users often substitute *phenotypic* for
*genetic* correlations, which violates the Smith–Hazel assumptions and can mislead -- when in
doubt, the rank-normalized `auto` (a rank-summation-style index) is the safe choice. Use
`economic_index` / `desired_gain` only when you genuinely have trustworthy economic weights or
desired-gain targets. See the per-method sections below and examples
`10_multitrait_method_auto.R` … `13_multitrait_method_desired_gain.R`.

### Practical Breeder Objective

Use `ng_breeder_selection_objective()` for user-facing workflows. It requires
every trait direction by default, chooses the strongest practical method
supported by the supplied inputs, and keeps threshold handling explicit. The
result still uses the same multi-trait scoring and OCS allocation engine.

```r
objective <- ng_breeder_selection_objective(
  trait = c("yield", "disease", "lodging"),
  column = c("pred_yield", "pred_disease", "pred_lodging"),
  direction = c("maximize", "minimize", "minimize"),
  economic_weight = c(1, 4, 2),
  max_value = c(NA, 20, 15),
  threshold_policy = "soft"
)

plan <- ng_optimize_breeder_selection_plan(
  scores = scores,
  objective = objective,
  n_crosses = 100,
  parent_K = parent_K,
  optimizer_method = "greedy_local",
  max_crosses_per_parent = 6,
  lambda_group = 0.05
)
```

Use this path when the package is being driven by a breeder, UI, or scripted
decision workflow. Use the lower-level `ng_multitrait_spec()` path when you
need direct control over method selection or are writing validation harnesses.

### Unknown Weights: Recommended Default

Use `method = "auto"` when breeders know the directions and rough thresholds
but not exact economic weights.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "quality"),
  column = c("pred_yield", "pred_disease", "pred_quality"),
  direction = c("maximize", "minimize", "maximize"),
  min_value = c(NA, NA, 65),
  max_value = c(NA, 20, NA)
)

scored <- ng_add_multitrait_score(scores, traits, method = "auto")
```

`auto` rank-normalizes all oriented traits, assigns equal normalized weights
when weights are missing, and applies soft threshold penalties.

Recommended when:

- a breeder wants balanced progress;
- exact trait weights are unknown;
- traits are in different units;
- thresholds are useful but should not discard every cross.

### Known Relative Weights

Use `method = "weighted"` when weights are meaningful and you do not want a
covariance-derived index.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "protein"),
  direction = c("maximize", "minimize", "maximize"),
  weight = c(3, 2, 1)
)

scored <- ng_add_multitrait_score(scores, traits, method = "weighted")
```

Recommended when weights come from breeder consensus, market classes, or a
simple operational priority.

### Economic Index

Use `method = "economic_index"` when economic weights are known but desired
gains are not.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "quality"),
  direction = c("maximize", "minimize", "maximize"),
  economic_weight = c(1, 4, 2)
)

scored <- ng_add_multitrait_score(scores, traits, method = "economic_index")
```

This estimates the covariance among oriented/scaled traits and solves
stabilized index coefficients from economic weights.

Recommended when:

- trait covariance matters;
- improving one trait can harm another;
- disease, quality, or lodging has clear economic importance.

### Desired-Gain Optimizer

Use `method = "desired_gain"` when the breeder can state desired changes.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "quality"),
  direction = c("maximize", "minimize", "maximize"),
  desired_change = c(5, 30, 3),
  economic_weight = c(1, 4, 2)
)

scored <- ng_add_multitrait_score(scores, traits, method = "desired_gain")
```

This uses desired changes, economic weights, and stabilized covariance to
estimate the index coefficients. Diagnostics are stored in:

```r
attr(scored, "multi_trait")
```

Recommended when:

- program goals are explicit;
- the breeder wants controlled response, not only ranking;
- traits are antagonistic or have different economic consequences.

### Threshold Selection

Thresholds are soft by default:

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease"),
  direction = c("maximize", "minimize"),
  min_value = c(70, NA),
  max_value = c(NA, 15),
  threshold_weight = c(1, 3)
)

scored <- ng_add_multitrait_score(scores, traits, method = "threshold")
```

Only use strict thresholds when failure is unacceptable:

```r
scored <- ng_add_multitrait_score(
  scores, traits,
  method = "threshold",
  strict_thresholds = TRUE
)
```

Recommended strict-threshold examples:

- disease score must be below a release threshold;
- quality must meet a market class;
- height/lodging must stay below mechanical-harvest limits.

### Multi-Trait Mating Plan

```r
plan <- ng_optimize_multitrait_mating_plan(
  scores = scores,
  traits = traits,
  n_crosses = 100,
  parent_K = parent_K,
  multitrait_method = "auto",
  max_crosses_per_parent = 6,
  lambda_group = 0.05
)
```

Recommended default:

- `multitrait_method = "auto"` if weights are unknown;
- `economic_index` when economic weights are known;
- `desired_gain` when explicit desired changes are known;
- keep `lambda_group > 0` if maintaining diversity matters.

### Practical Cross Priority Tiers

Most production programs need a larger crossing menu than 20 crosses. A
practical workflow is to request about 100 selected crosses, control parent use
and relatedness during allocation, then rank the selected list into priority
tiers. The priority tier is not a replacement for the mating optimizer; it is a
decision layer for deciding which selected crosses should be made first if
capacity, seed, greenhouse space, or nursery slots become limiting.

```r
plan <- ng_rank_cross_priority(
  plan,
  score_col = "multi_trait_score",
  pair_kinship_col = "pair_kinship",
  threshold_violation_col = "multi_trait_threshold_violation",
  score_weight = 1,
  kinship_weight = 0.15,
  threshold_weight = 1,
  breaks = c(0.10, 0.35, 0.70, 1.00)
)

priority_summary <- ng_cross_priority_summary(plan)
```

With 100 selected crosses, the default breaks create:

- 10 `highly_priority` crosses;
- 25 `priority` crosses;
- 35 `medium_priority` crosses;
- 30 `low_priority` crosses.

Increase or decrease `n_crosses` for the breeding program's crossing capacity.
Change `breaks` to alter tier sizes. Increase `kinship_weight` when avoiding
close-parent crosses matters more, and increase `threshold_weight` when crosses
that violate disease, quality, lodging, or market thresholds should move down
the list.

For breeder-facing delivery, write an Excel workbook from the selected plan:

```r
ng_write_cross_priority_workbook(
  output_path = "results/breeder_crossing_plan.xlsx",
  crosses = plan,
  scored = scored,
  trait_directions = traits,
  parent_use = selected_parent_use,
  duplicate_pairs = putative_duplicate_pairs,
  n_crosses_requested = 100
)
```

The workbook reports objective evidence such as priority tier, score,
kinship, threshold violation, parent-use flags, duplicate-QC flags, and the
strongest favorable/risk trait evidence. It intentionally leaves
`Breeder_Rationale`, `Breeder_Notes`, `Final_Decision`, and
`Crossing_Status` blank because final rationale differs by breeder, market,
nursery capacity, and cycle.

## Native var_complex, SimpleMating, And AlphaMate-Style Components

These are native additions to `nextgenCrossDesign`. They do not require the old
tools.

### Native var_complex Usefulness

Use the one-call workflow when users want PopVar-like usefulness without adding
PopVar as a package dependency:

```r
result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  trait_value_metric = "var_complex"
)
```

`var_complex` is native code in this package. It does not call PopVar and it is
not a copy-paste dependency. It uses the package's cross means and the best
available within-family variance column for the trait-value calculation.

### AlphaMate-Style Target-Degree Allocation

Use this when you want AlphaMate-like allocation without the AlphaMate
executable:

```r
result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  allocation_method = "alphamate_style",
  alphamate_mode = "ModeOptTarget1",
  alphamate_target_degree = 45,
  alphamate_max_contributions = 20
)
```

The lower-level equivalent is:

```r
plan <- ng_alphamate_style_select(
  scores = scores,
  criterion_col = "multi_trait_score",
  n_crosses = 100,
  parent_K = parent_K,
  mode = "ModeOptTarget1",
  target_degree = 45,
  max_contributions = 6
)
```

Modes:

- `ModeMaxCriterion`: prioritizes gain.
- `ModeMinCoancestry`: prioritizes coancestry reduction.
- `ModeOptTarget1`: searches a gain-diversity frontier and chooses according
  to `target_degree`.

Recommended target-degree starting points:

- `target_degree = 20`: gain-heavy.
- `target_degree = 45`: balanced default.
- `target_degree = 70`: diversity-heavy.

This is a package-developed gain-diversity frontier allocator. It is not the
AlphaMate executable and does not claim to solve AlphaMate's target-degree
constraint exactly.

## Optional Exact External Tool Calls

Exact external tools are optional and should be used for benchmarking,
comparison, or breeder trust.

```r
scores <- ng_add_external_baseline_scores(
  scores = scores,
  geno = geno,
  effects = effects,
  marker_map = marker_map,
  methods = c("popvar_uc_ocs10_lps1", "simple_usefa_ocs10_lps1"),
  popvar_engine = "auto",
  simplemating_engine = "auto"
)
```

Engine options:

- `"native"`: always use native `nextgenCrossDesign` behavior.
- `"external"`: require exact external package and return unavailable/failed
  status if not possible.
- `"auto"`: try exact external behavior when possible, otherwise use native
  proxy fallback.

AlphaMate exact execution:

```r
plan <- ng_select_alphamate(
  scores = scores,
  criterion_col = "multi_trait_score",
  n_crosses = 100,
  parent_K = parent_K,
  executable = "external/AlphaMate/binaries/AlphaMate.exe",
  target_degree = 45
)
```

The one-call workflow can also call the executable directly:

```r
result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  allocation_method = "alphamate_executable",
  alphamate_executable = "external/AlphaMate/binaries/AlphaMate.exe",
  alphamate_runtime_path = Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = ""),
  alphamate_mode = "ModeOptTarget1",
  alphamate_target_degree = 45
)
```

Use exact external calls when:

- reproducing published or legacy workflows;
- validating native approximations;
- comparing against PopVar, SimpleMating, or AlphaMate for reporting.

For default production use, prefer native `nextgenCrossDesign` components.

## Breeder-Facing Excel Crossing Plan

Use this when the output needs to be understandable without opening R.

Required candidate columns:

- `line` or `id`;
- `cross` or `family`;
- optional `seed_color`;
- numeric traits such as `yield`, `harvestability`, `tsw`.

Criteria example:

```r
criteria <- data.frame(
  trait = c("yield", "harvestability", "tsw"),
  label = c("Yield (bu/a)", "Harvestability", "TSW"),
  weight = c(0.5, 0.3, 0.2),
  direction = c("increase", "decrease", "increase"),
  meaning = c(
    "Higher yield is better",
    "Lower score is better",
    "Higher seed weight is preferred"
  )
)
```

Write workbook:

```r
ng_write_crossing_plan_workbook(
  output_path = "results/crossing_plan.xlsx",
  candidates = candidates,
  criteria = criteria,
  top_n = 5,
  cross_n = 10,
  title = "Top 5 Crossing Candidates with Crossing Plan - 2026"
)
```

The workbook includes:

- `Top Parents`: selected unique-family parents and breeding-use notes;
- `Recommended Crosses`: ranked crosses, tiers, rationale, and full-sib checks;
- `Selection Index - All Lines`: all candidates and scores;
- `Criteria`: weights, directions, and interpretation;
- `Family Screening`: within-family ranks and full-sib exclusions;
- embedded selection-index chart.

Command-line wrapper:

```powershell
$env:NG_CROSSING_PLAN_CANDIDATES="candidates.csv"
$env:NG_CROSSING_PLAN_CRITERIA="criteria.csv"
$env:NG_CROSSING_PLAN_OUT="results/crossing_plan.xlsx"
$env:NG_CROSSING_PLAN_TOP_N="5"
$env:NG_CROSSING_PLAN_CROSS_N="10"
Rscript tools\export_crossing_plan_workbook.R
```

## Reports And Frontend Payloads

The full framework-agnostic frontend integration kit lives in **`docs/frontend/`** (start at
`docs/frontend/README.md`). The single headless entry point is the JSON-in / JSON-out run wrapper:

```bash
NG_RUN_CONFIG=my_config.json NG_RUN_RESULT_OUT=results/run_result.json \
  Rscript tools/run_cross_prediction_json.R
```

It reads a JSON config whose keys mirror `ng_run_cross_prediction()` (see
`docs/frontend/contracts/config_schema.json`; unknown keys error) and writes an `ng_run_result.v1`
JSON plus the optional workbook/figures. The exporters below produce the supporting decision/benchmark
payloads.

Backend capability registry:

```r
registry <- ng_backend_capability_registry()
ng_write_backend_capability_registry_json("results/backend_capabilities.json")
```

Head-to-head benchmark:

```powershell
Rscript tools\run_head_to_head_benchmark.R
```

Dashboard JSON:

```powershell
Rscript tools\export_head_to_head_dashboard_json.R
```

Standalone visual report:

```powershell
Rscript tools\render_head_to_head_visual_report.R
```

Use these outputs for frontend integration, evidence review, and method
comparison.

## Validation Runners

Run focused multi-trait validation:

```powershell
Rscript tools\run_multitrait_validation.R
```

Run replicated multi-trait grid:

```powershell
Rscript tools\run_multitrait_validation_grid.R
```

Run crop-aware multi-trait grid:

```powershell
Rscript tools\run_multitrait_crop_validation_grid.R
```

Run formal head-to-head benchmark:

```powershell
Rscript tools\run_head_to_head_benchmark.R
```

Recommended validation defaults:

- start with small parent sizes and one replicate for a smoke test;
- increase `*_REPS` only after the pipeline is stable;
- keep exact external baselines optional for CI-scale runs;
- use head-to-head outputs before claiming superiority over older tools.

## Polyploid Workflows

### Autotetraploid 4x

Use these for potato/cassava-like experimental screens:

- `ng_poly4x_as_dosage_matrix()`;
- `ng_poly4x_parent_relationship()`;
- `ng_poly4x_score_crosses()`;
- `ng_poly4x_policy()`;
- `ng_poly4x_select_topn()`.

Policy modes:

- `gain`;
- `diversity`;
- `ocs`.

Run smoke/benchmark scripts:

```powershell
Rscript tools\run_poly4x_benchmark.R
powershell -ExecutionPolicy Bypass -File tools\run_poly4x_grid.ps1
```

Recommended default:

- use 4x modules for research and validation;
- do not generalize 4x validation to all polyploids;
- validate the dosage model and pairing biology for the target crop.

### Allopolyploid/Subgenome

Use these when subgenomes can be modeled as separate dosage matrices:

- `ng_poly_subgenome_as_dosage_list()`;
- `ng_poly_subgenome_parent_relationship()`;
- `ng_poly_subgenome_score_crosses()`;
- `ng_poly_policy()`.

Recommended default:

- use for wheat-like disomic/subgenome experiments;
- keep subgenome weights explicit;
- treat outputs as model-family-specific, not universal.

## Choosing A Default Method

Use this decision path:

1. Run `ng_preflight_input_tables(putative_duplicate_check = TRUE)`.
2. If blocker issues exist, fix data first.
3. If putative duplicate genotype profiles are present, use
   `putative_duplicate_action = "remove"` unless a breeder/data-manager
   intentionally keeps both records and documents why.
4. Run scoring on `qc$cleaned_tables` when duplicate removal was requested.
5. If the target is diploid DH/RIL and one trait dominates, start with
   `ng_design_crosses()`.
6. If multiple traits matter and weights are unknown, use multi-trait `auto`.
7. If economic weights are credible, use `economic_index`.
8. If desired changes are credible, use `desired_gain`.
9. Add `lambda_group` and parent-use constraints when diversity matters.
10. Use native `var_complex`, SimpleMating-style, and AlphaMate-style
    components as add-ons or benchmarks, not as hard dependencies.
11. Use the Excel workbook for breeder-facing crossing recommendations.
12. Use validation runners before making method superiority claims.

## Troubleshooting

### R CMD check reports genomicMateSelectR NOTE

This is expected when `genomicMateSelectR` is not installed. It is listed as an
optional enhancement, not a required dependency.

### PopVar or SimpleMating is unavailable

Use native engines:

```r
ng_add_popvar_scores(..., engine = "native")
ng_add_simplemating_scores(..., engine = "native")
```

For the main user workflow, prefer
`ng_run_cross_prediction(trait_value_metric = "var_complex")` instead of
requiring PopVar. Use `engine = "auto"` only for comparison workflows that try
exact external calls and fall back natively.

### AlphaMate executable is unavailable

Use:

```r
ng_alphamate_style_select(...)
```

Only use `ng_select_alphamate()` when the executable is installed and the path
is configured.

### No feasible mating plan

Relax one or more constraints:

- increase `max_crosses_per_parent`;
- reduce `min_unique_parents`;
- increase `max_pair_kinship`;
- lower `lambda_group`;
- switch `method = "greedy_local"` for a diagnostic run.

### Strict thresholds select nothing

Use soft thresholds first:

```r
ng_add_multitrait_score(..., strict_thresholds = FALSE)
```

Then inspect `multi_trait_threshold_violation`.

### Workbook cannot be written

Install `openxlsx`:

```r
install.packages("openxlsx")
```

Then rerun `ng_write_crossing_plan_workbook()`.

## Backend Verification Checklist

Before trusting a backend change:

```powershell
Rscript tests\data_preflight.R
Rscript tests\external_compatibility.R
Rscript tests\crossing_plan_workbook.R
Rscript tools\check_r_package.R
```

For the full regression directory:

```powershell
$tests = Get-ChildItem tests -Filter *.R | Where-Object { $_.Name -ne 'helper_load.R' } | Sort-Object Name
foreach ($test in $tests) {
  Rscript $test.FullName
  if ($LASTEXITCODE -ne 0) { throw "Test failed: $($test.Name)" }
}
```

Expected current package-check status:

```text
PASS with one NOTE if genomicMateSelectR is not installed.
```

No backend change should be considered complete until the focused test and
package check pass.
