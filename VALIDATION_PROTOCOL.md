# Diagnostic-First Validation Protocol

This protocol validates `nextgen_cross_design` as a framework before changing
router or score logic. The primary decision target is realized `top10_gv` for
DH/RIL cross design. Secondary checks are `mean_gv`, `max_gv`, parent-use
balance, group coancestry, and runtime.

## Principle

Do not tune `ng_meta_router_ocs*` or any other method until diagnostics show
which component is failing:

1. marker-effect and variance calibration;
2. score ranking value;
3. mate-allocation value;
4. multi-cycle robustness across parent counts.

`var_simple`, PopVar, SimpleMating, PMV-balanced, recombination-GEBV,
meta-portfolio, meta-selector, meta-router, real AlphaMate when executable
dependencies are available, and the frontier policy remain controls. A
promoted method must be Pareto-competitive against these controls across
parent-size grids, not only in one tuned run.

## Runnable Entry Point

Use the wrapper from the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_framework_validation.ps1
```

By default it runs only smoke checks. Heavy phases are selected with
`NG_VALIDATION_PHASE`:

```powershell
$env:NG_VALIDATION_PHASE = "family"    # family-level calibration
$env:NG_VALIDATION_PHASE = "allocator" # score-vs-allocation crosscheck
$env:NG_VALIDATION_PHASE = "grid"      # replicated parent-size grid
$env:NG_VALIDATION_PHASE = "all"       # family, allocator, then grid
```

For full validation, run outside the sandbox with C++ enabled:

```powershell
$env:NG_VALIDATION_USE_CPP = "1"
$env:NG_VALIDATION_REPS = "3"
$env:NG_VALIDATION_PARENT_SIZES = "20,30,40,50,60,70,80"
$env:NG_VALIDATION_EFFECT_TRAINING_N = "400"
$env:NG_VALIDATION_PHASE = "all"
powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_framework_validation.ps1
```

## Phases

### 1. Family Calibration

Purpose: determine whether predicted variance and top-tail scores are
statistically trustworthy before allocator tuning.

Default wrapper outputs use prefix:
`results/diagnostic_family_5k_*`.

Primary reads:

- `*_metric_summary_avg.csv`
- `*_metric_summary_by_rep.csv`
- `*_families.csv`
- `*_status.csv`
- `*_timing_summary.csv`

If PMV or recombination variance has poor correlation, bad slope, or unstable
top-tail enrichment, fix marker-effect/variance calibration before changing the
router.

### 2. Allocator Crosscheck

Purpose: separate score quality from mate-allocation quality by comparing
top-N, adaptive OCS, and SimpleMating-style constrained selection under the
same score families.

Default wrapper outputs use prefix:
`results/diagnostic_allocator_5k_*`.

Primary reads:

- `*_overall_avg.csv`
- `*_comparison_avg.csv`
- `*_selection_avg.csv`

If a score wins under top-N but loses under OCS, tune allocation. If it loses
under every allocator, tune the score/calibration instead.

### 3. Parent-Size Grid

Purpose: validate multi-cycle robustness across breeding-program sizes.

Default wrapper outputs use prefix:
`results/diagnostic_parent_grid_5k_*`.

Primary reads:

- `diagnostic_parent_grid_5k_overall_avg.csv`
- `diagnostic_parent_grid_5k_comparison_avg.csv`
- `diagnostic_parent_grid_5k_winner_summary.csv`
- `diagnostic_parent_grid_5k_selection_avg.csv`

Promotion requires Pareto-competitive `top10_gv` against `var_simple`,
PopVar/SimpleMating-style usefulness, and existing meta controls across parent
sizes. If the result is mixed, report a decision rule or Pareto frontier rather
than claiming a universal winner. `ng_frontier_policy_ocs*` is the current
implementation of that decision-rule approach, but its source tag and policy
bands must be revalidated for each materially different crop, dataset, or
breeding setup.

### 4. Real AlphaMate External Screen

Purpose: separate the repo's AlphaMate-style OCS allocator from the official
AlphaMate executable.

Use `alphamate_opt*` methods after setting `NG_ALPHAMATE_EXE` if the binary is
not under `external/AlphaMate/binaries/AlphaMate.exe`. On Windows, also set
`NG_ALPHAMATE_RUNTIME_PATH` if `libiomp5md.dll` is not on `PATH`.

The standard replicated 5K screen is:

```powershell
powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_alphamate_external_grid_5k.ps1
```

Primary reads:

- `*_overall_avg.csv`
- `*_comparison_avg.csv`
- `*_selection_summary.csv`
- `*_selection_avg.csv`

Require `alphamate_exit_code = 0` in the per-scenario selection summaries.
Treat `ModeMaxCriterion` and `ModeMinCoancestry` outputs as diagnostic modes,
not promoted cross-design baselines, unless they produce valid non-self,
non-repeated crosses for the benchmark candidate table.

### 5. Multi-Trait Head-To-Head Screen

Purpose: compare practical multi-trait breeding objectives under the same crop
scenario, parent set, realized family means, trait directions, and crossing
budget.

Use the CI-friendly style-proxy runner first:

```powershell
Rscript nextgen_cross_design\tools\run_head_to_head_benchmark.R
```

Primary reads:

- `*_summary.csv`
- `*_selections.csv`
- `*_scores.csv`
- `*_comparisons.csv`
- `*_winner_summary.csv`
- `*_config.csv`
- `*_method_registry.csv`

The default registry compares `nextgen_auto_ocs`,
`nextgen_economic_index_ocs`, and `nextgen_desired_gain_ocs` against
`popvar_style_weighted_topn`, `simplemate_style_threshold_topn`, and
`alphamate_style_weighted_ocs`. The style-proxy baseline rows must keep
`implementation = "style_proxy"`, `exact_external_status`, and
`fallback_reason` fields. Do not report these rows as exact PopVar,
SimpleMating, or AlphaMate evidence. Exact external jobs should be added as a
separate optional tier when the required packages or executable are installed.

## Defaults

- Objective: realized `top10_gv`.
- Breeding target: DH/RIL from inbred-line crosses.
- AlphaSimR threads: `1` for deterministic comparisons.
- Effect training: fixed `400` unless overridden.
- C++: enabled by wrapper default for heavy phases, but smoke checks report
  whether C++ is actually available.

## Crop-Genome Portability

Use `tools/run_crop_genome_scenarios.R` after the standard parent-size grid to
test whether the frontier policy remains reasonable under different crop-like
genome architectures. Use `tools/run_crop_genome_alphamate_grid.ps1` when the
screen should include real AlphaMate target-degree controls. The built-in
scenarios vary chromosome number, genetic map length, marker density, QTL
density, heritability, founder count, and effect-training size.

This screen is intentionally scoped. It does not validate universal crop
agnosticism, polyploid inheritance, clonal propagation, or crop-specific
recombination landscapes. A production recommendation for a new crop should
produce a new crop validation grid and either a new `NG_FRONTIER_POLICY_SPEC` or
an explicit `ng_crop_aware_policy_ocs*` dispatch rule rather than relying on the
default 2026-05-03 policy bands.

Current crop templates include diploid DH/RIL scenarios for compact selfing,
maize, barley, field pea, and cassava, plus diploidized stress approximations
for bread wheat, potato, cassava tetraploid, and sugarcane. Do not treat those
stress approximations as true polyploid validation until the benchmark has a
polyploid dosage and progeny-generation path.

The crop-aware policy must report its decision through `crop_policy_*` and
`pred_ng_crop_policy_*` columns. The current smoke-derived rules are limited:
ordinary diploid DH/RIL templates route to `ng_frontier_policy_ocs10_lps2`,
`cassava_tetraploid_stress` routes first to `alphamate_opt60`, and
`sugarcane_polyploid_stress` routes first to `alphamate_opt45`.

## Autotetraploid 4x Validation

Use `tools/run_poly4x_benchmark.R` for true 4x validation. This runner must
not call `makeDH()` and must report `used_dh = 0` in selection summaries.

Minimum smoke validation:

```powershell
Rscript nextgen_cross_design\tests\poly4x_dosage.R
Rscript nextgen_cross_design\tests\poly4x_simulation.R
Rscript nextgen_cross_design\tests\poly4x_scoring.R
Rscript nextgen_cross_design\tests\poly4x_runner_smoke.R
```

Minimum scenario screen:

```powershell
$env:NG_POLY4X_GRID_SCENARIOS = "potato_autotetraploid_4x,cassava_autotetraploid_4x"
$env:NG_POLY4X_GRID_PARENT_SIZES = "20,40"
$env:NG_POLY4X_REPS = "1"
$env:NG_POLY4X_CYCLES = "1"
powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_poly4x_grid.ps1
```

Promotion requires sampled-progeny predictions to rank independent realized 4x
family outcomes better than simple 4x variance top-N and to keep parent-use and
coancestry diagnostics within the configured limits.
