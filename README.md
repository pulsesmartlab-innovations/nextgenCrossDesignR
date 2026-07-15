# nextgen_cross_design

Fresh implementation track for cross prediction and mate allocation.

Current package version: `nextgenCrossDesign` 0.4.0.

Shareable backend source tarball for outside-repo developers:

```text
dist/nextgenCrossDesign_0.4.0.tar.gz
```

The tarball is backend-only. It contains the R package source, compiled-kernel
source, help files, and backend README. For an outside frontend developer, share
the tarball together with:

- [`docs/frontend/FRONTEND_DEVELOPER_TEMPLATE.md`](docs/frontend/FRONTEND_DEVELOPER_TEMPLATE.md)
- [`docs/BACKEND_USER_GUIDE.md`](docs/BACKEND_USER_GUIDE.md)
- [`VALIDATED_STATE.md`](VALIDATED_STATE.md)
- [`docs/V0_3_4_QC_PRIORITY_WORKFLOW.md`](docs/V0_3_4_QC_PRIORITY_WORKFLOW.md)
- [`docs/V0_3_2_DIMINISHING_RETURNS.md`](docs/V0_3_2_DIMINISHING_RETURNS.md)

For the current evidence-backed claims, incorporated scripts, and limits on
PopVar/SimpleMating/AlphaMate/polyploid comparisons, see
[`VALIDATED_STATE.md`](VALIDATED_STATE.md). That file is the short-form status
record; `BENCHMARK_NOTES.md` contains the long benchmark history.

For a guided backend user manual with recommended defaults by breeding
scenario, input schemas, native `var_complex`, SimpleMating-style and
AlphaMate-style components, validation runners, and report exports, see
[`docs/BACKEND_USER_GUIDE.md`](docs/BACKEND_USER_GUIDE.md).

For an R-user vignette in the style of breeding-package tutorials, see
[`vignettes/nextgenCrossDesign.Rmd`](vignettes/nextgenCrossDesign.Rmd). It
now starts with the recommended first run using `ng_run_cross_prediction()`:
files, package QC, duplicate-genotype handling, marker-effect estimation,
cross scoring, optimizer/allocation settings, priority ranking, and outputs. The
same wrapper-first workflow is available as a copy-paste R script at
[`inst/examples/00_user_friendly_cross_prediction.R`](inst/examples/00_user_friendly_cross_prediction.R).
The detailed lower-level workflow is available at
[`inst/examples/01_chronological_cross_prediction_workflow.R`](inst/examples/01_chronological_cross_prediction_workflow.R).
The rest of the vignette walks through `var_complex`, PMV/VPM/`var_simple`
choices, AlphaMate-style and external AlphaMate allocation, cross-number
sweeps, validation, external-style baselines, policy runners, and polyploid
research workflows.

For the frontend developer handoff, use
[`docs/frontend/FRONTEND_DEVELOPER_TEMPLATE.md`](docs/frontend/FRONTEND_DEVELOPER_TEMPLATE.md)
as the main product/UX specification, and
[`docs/frontend/contracts/`](docs/frontend/contracts/) as the machine-readable
run contract. This repository is the backend only; the frontend is the separate
`nextgenCrossWorkbench` Shiny package.

The design is intentionally separated into four layers:

1. phenotype and marker-effect estimation;
2. cross-level mean, variance, usefulness, and uncertainty metrics;
3. candidate screening;
4. constrained mate allocation and parent contribution optimization.

## v0.4.0 Polyploid, Dominance, and COMA Removal

Version 0.4.0 rebuilds the polyploid path correctly (ploidy-general, additive by
default, dominance optional) and removes the COMA integration (benchmarked, never
beat the native allocator; a third-party package the backend no longer depends on).
Highlights: correct allele-frequency polyploid GRM (`ng_polyploid_grm`, VanRaden or
Yang) plus a digenic dominance GRM and ploidy-aware QC; additive+dominance marker
effects and genotypic-value prediction for clonal crops; dominance-aware cross
scoring (heterosis + within-family variance, double reduction, C++-accelerated); and
a one-call `ng_design_crosses_poly(dominance =, gain =, double_reduction =, grm_method =)`.
See `inst/examples/29_polyploid_qc_and_grm.R` through `31_polyploid_dominance_crossing.R`
and the vignette "Polyploid and dominance-aware design" section.

## v0.3.13 User-Workflow Additions

Version 0.3.13 adds the remaining practical examples users need before
independent testing: duplicate QC reporting/removal, exact input matching and
map units, cross-number sweeps, workbook and figure outputs, allocation-method
comparison, posterior robust mate allocation, and full-posterior PMV shortlist
reruns. These are available at `inst/examples/14_qc_duplicate_removal_and_reporting.R`
through `inst/examples/20_full_posterior_pmv_shortlist.R`.

## v0.3.12 User-Workflow Additions

Version 0.3.12 adds sequential, runnable examples for the wrapper choices that
users commonly confuse. `inst/examples/08_prediction_mode_trait_by_trait.R`
and `inst/examples/09_prediction_mode_index_as_trait.R` show the two
`prediction_mode` paths. `inst/examples/10_multitrait_method_auto.R` through
`inst/examples/13_multitrait_method_desired_gain.R` show `auto`, `weighted`,
`economic_index`, and `desired_gain` as separate scripts with method-specific
parameter blocks at the top.

## v0.3.11 User-Workflow Additions

Version 0.3.11 exposes the PMV/posterior workflow knobs directly through
`ng_run_cross_prediction()`: `method_varPMV`, `ril_mode`,
`run_posterior_prediction`, `posterior_method`, `nIter`, `burnIn`, and
`use_parallel`. These are wired into the executed workflow. `method_varPMV =
"full_posterior"` fits and uses the full marker-effect covariance in PMV
scoring, and `run_posterior_prediction = TRUE` returns posterior cross-score
tables in `result$posterior_predictions`. The new
`inst/examples/07_trait_value_metric_parameter_guide.R` script compares the
trait-value metrics and runs a small posterior-prediction example.

## v0.3.10 User-Workflow Additions

Version 0.3.10 makes `lpSolve` a required package import because the OCS MIP
optimizers are part of the user-facing optimizer choices. The optimizer guide at
`inst/examples/06_optimizer_parameter_guide.R` now runs `auto`,
`greedy_local`, `repair_local`, `mip_linear`, and `mip_contribution` as normal
package functionality instead of treating the MIP rows as optional.

## v0.3.9 User-Workflow Additions

Version 0.3.9 adds a dedicated package OCS example at
`inst/examples/05_ocs_user_run.R`. It keeps `allocation_method = "ocs"` separate
from the AlphaMate-style examples and shows the OCS-specific crossing capacity,
parent-use, kinship, optimizer, and penalty settings in one runnable script.

## v0.3.8 User-Workflow Additions

Version 0.3.8 adds a dedicated external AlphaMate example at
`inst/examples/04_alphamate_executable_user_run.R`. It keeps the package-native
`alphamate_style` path separate from `allocation_method = "alphamate_executable"`
and shows all executable-specific settings, including `NG_ALPHAMATE_EXE`,
runtime paths, work directory retention, evolutionary-search controls, and
thread count.

## v0.3.7 User-Workflow Additions

Version 0.3.7 adds a native `trait_value_metric = "var_complex"` option to
`ng_run_cross_prediction()`. This is the package-native PopVar-inspired
usefulness metric; it does not call PopVar and does not make PopVar a required
dependency. The same wrapper now exposes `allocation_method`, with
`"ocs"` for package OCS, `"alphamate_style"` for the package-developed
AlphaMate-style gain-diversity allocator, and `"alphamate_executable"` when the
user wants to call an installed AlphaMate executable directly.

The immediate goal is not to clone `genomicMateSelectR`, `SimpleMating`, or AlphaMate. The goal is to combine their strongest ideas in a faster and more explicit framework:

- `genomicMateSelectR`: phased haplotype plus recombination-aware progeny variance, with VPM/PMV distinction.
- `SimpleMating`: breeder-facing criteria such as MPV, TGV, and usefulness, followed by constrained cross selection.
- AlphaMate: joint selection, diversity maintenance, and mate allocation with valid mating-plan constraints.

The current working direction is not a single tuned method. The active default
is now `ng_frontier_policy_ocs*`, an empirical policy layer that delegates to
the best validated method family for the current parent-count band. It records
the delegated method, source, band, candidates, and fallback status through
`frontier_policy_*` and `pred_ng_frontier_policy_*` columns. This policy is not
crop-agnostic or dataset-agnostic; it is a validated default for the current
DH/RIL AlphaSimR evidence and should be revalidated when crop biology, marker
density, training population size, trait architecture, generation scheme, or
external package availability changes.

For crop portability screens, `ng_crop_aware_policy_ocs*` adds a second,
auditable dispatch layer. It reads `NG_CROP_*` metadata from the crop-genome
harness, records `crop_policy_*` and `pred_ng_crop_policy_*` diagnostics, routes
ordinary diploid DH/RIL screens through the frontier policy, and routes the
current cassava tetraploid and sugarcane polyploid stress templates to the best
real-AlphaMate target seen in the smoke grid. This is a practical stress-screen
policy, not proof of true polyploid generality.

Within the DH/RIL frontier policy, `ng_meta_router_ocs*` remains a core control
method. It builds complete candidate mating plans and chooses among
recombination-GEBV, scaled PMV-balanced usefulness, the blended meta portfolio,
PopVar/SimpleMating-style usefulness when available, and `var_simple`. It uses
prior-cycle method-family performance, current predicted meta rank, candidate
gain rank, diversity, and SNP-effect reliability. Keep `ng_meta_router_ocs*`,
`ng_meta_selector_ocs*`, and `ng_meta_portfolio_ocs*` in validation grids as
controls. Hard minimum unique-parent constraints are opt-in because they can
erase too much selection intensity.

## Multi-Trait Selection

Breeding programs can now build a single selection objective from multiple
score columns while increasing some traits and decreasing others. Use
`ng_multitrait_spec()` to declare trait columns, directions, optional weights,
and optional thresholds; then use `ng_add_multitrait_score()`,
`ng_multitrait_select_topn()`, or `ng_optimize_multitrait_mating_plan()`.

For user-facing workflows, prefer `ng_breeder_selection_objective()`. It wraps
the lower-level multi-trait spec, requires explicit trait directions by default,
auto-selects the strongest practical method supported by the supplied inputs,
and sends the resulting score through the same top-N or OCS allocation path.

## Data QC and Putative Duplicate Removal

The current package keeps putative duplicate genotype handling embedded in package QC. Use
`ng_preflight_input_tables(putative_duplicate_check = TRUE,
putative_duplicate_action = "remove")` before scoring. The package detects
near-identical genotype profiles, keeps one parent per duplicate cluster, and
returns cleaned inputs in `qc$cleaned_tables` plus audit metadata in
`qc$cleaning$putative_duplicates`.

The direct duplicate APIs remain exported for reporting:

- `ng_detect_putative_duplicates()` returns duplicate pairs, clusters, nearest
  neighbors, marker summaries, and optional similarity matrices.
- `ng_plot_putative_duplicates()` writes PNG/PDF duplicate heatmaps.
- `ng_write_data_preflight_json()` records duplicate-removal metadata for UI or
  batch logs; full cleaned genotype/phenotype tables are returned by
  `ng_preflight_input_tables()` in R rather than embedded in JSON.

The recommended reproducible user workflow is:

```powershell
Rscript inst/examples/00_user_friendly_cross_prediction.R
```

That script calls `ng_run_cross_prediction()`, runs package QC, uses cleaned
tables after duplicate removal, estimates marker effects trait-by-trait, selects
crosses, assigns priority tiers, writes figures, and writes the Excel workbook
when `openxlsx` is installed.

The same workflow now requires an explicit input-matching contract when file
column names are not obvious: `phenotype_id_col`, `genotype_id_col`,
`direction_trait_col`, `direction_column_col`, `direction_direction_col`,
`map_marker_col`, `map_chr_col`, `map_pos_bp_col`, `map_position_unit`, and
`bp_per_cm`. The result includes `result$input_match_audit`, which records the
matched parent order, marker order, direction-file columns, marker-map columns,
and base-pair to cM conversion used for the run.

When economic weights are known:

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "lodging", "disease"),
  direction = c("maximize", "minimize", "minimize"),
  weight = c(2, 1, 1)
)
plan <- ng_optimize_multitrait_mating_plan(scores, traits, n_crosses = 100)
```

When weights are unknown, use the default `method = "auto"`. It rank-normalizes
traits, assigns equal normalized weights unless weights are supplied, and keeps
the chosen weights in `attr(scored, "multi_trait")`.

When relative economic weights are known but exact desired gains are not, use
`method = "economic_index"`. This keeps the weighted-index interface but
estimates the covariance among oriented/scaled traits and solves
covariance-aware coefficients from `economic_weight` or `weight`.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "dry_matter"),
  direction = c("maximize", "minimize", "maximize"),
  economic_weight = c(1, 4, 2)
)
scored <- ng_add_multitrait_score(scores, traits, method = "economic_index")
```

For a full desired-gain/economic objective, use `method = "desired_gain"` with
`desired_change` and `economic_weight`. The implementation converts traits to
the beneficial direction, estimates the trait covariance from the candidate
cross table, solves stabilized desired-gain index coefficients, and records the
coefficients plus target and predicted response vectors in the multi-trait
diagnostics.

```r
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "dry_matter"),
  direction = c("maximize", "minimize", "maximize"),
  desired_change = c(5, 35, 3),
  economic_weight = c(1, 4, 2)
)
plan <- ng_optimize_multitrait_mating_plan(
  scores, traits, n_crosses = 100, multitrait_method = "desired_gain"
)
```

Thresholds are soft by default: a cross that misses a minimum or maximum target
is penalized but not automatically discarded. This is safer when no candidate
cross satisfies every trait target. Use `strict_thresholds = TRUE` only for hard
quality, disease, or market cutoffs.

For practical crossing programs, request a realistic list such as 100 crosses
and then rank the selected list with `ng_rank_cross_priority()`. The default
tiers are `highly_priority`, `priority`, `medium_priority`, and `low_priority`;
with 100 crosses these default to 10, 25, 35, and 30 crosses. Users can change
`n_crosses`, priority breaks, and score/kinship/threshold weights for their
own crossing capacity.

Use `ng_write_cross_priority_workbook()` when the result needs to be reviewed
outside R. The workbook gives objective evidence columns and blank editable
fields for `Breeder_Rationale`, `Breeder_Notes`, `Final_Decision`, and
`Crossing_Status`; it does not impose one fixed rationale sentence on every
breeding program.

## Portfolio Size and Diminishing Returns

Version 0.3.2 adds a portfolio-size sweep for deciding how many crosses to make
without guessing. `ng_optimize_mating_plan_curve()` runs
`ng_optimize_mating_plan()` across a K range and returns total gain, marginal
gain, group coancestry, estimated effective population size, unique parent use,
and the recommended K under one of four criteria:

- `elbow_relative`: smallest K where marginal gain drops below a relative
  threshold.
- `elbow_kneedle`: normalized curvature detection for smooth diminishing-return
  curves.
- `ne_target`: smallest K that reaches an effective-population-size floor.
- `coancestry_budget`: largest K that stays below a group-coancestry ceiling.

`ng_plot_diminishing_returns()` renders the curve when `ggplot2` is installed.
These helpers do not change the existing optimizer defaults; they are opt-in
decision support for breeders who need to balance gain against diversity and
parent-use pressure.

```r
curve <- ng_optimize_mating_plan_curve(
  scores,
  K_range = 3:20,
  parent_K = parent_K,
  lambda_group = 0.05,
  max_crosses_per_parent = 4L,
  criterion = "ne_target",
  ne_min = 50
)

attr(curve, "elbow_K")
ng_plot_diminishing_returns(curve)
```

See [`docs/V0_3_2_DIMINISHING_RETURNS.md`](docs/V0_3_2_DIMINISHING_RETURNS.md)
for the release note and validation summary.

Smoke validation runner:

```text
Rscript tools/run_multitrait_validation.R
```

The runner writes `*_summary.csv`, `*_selections.csv`, and `*_scores.csv` to
`results/` by default. It compares `auto`, `weighted`, `economic_index`,
`desired_gain`, and `threshold` on a deterministic synthetic scenario with
mixed positive and negative trait directions.

Replicated validation grid:

```text
Rscript tools/run_multitrait_validation_grid.R
```

Useful grid overrides:

```text
NG_MULTITRAIT_GRID_PARENT_SIZES=20,40,60,80
NG_MULTITRAIT_GRID_REPS=3
NG_MULTITRAIT_GRID_CROSSES=6
NG_MULTITRAIT_GRID_METHODS=auto,weighted,economic_index,desired_gain,threshold
NG_MULTITRAIT_GRID_OCS_LAMBDA_GROUP=0.05
```

The grid writes `*_summary.csv`, `*_selections.csv`, `*_scores.csv`,
`*_winner_summary.csv`, and `*_config.csv`. Winner summaries average replicates
within each parent size, maximize realized index, yield, quality, and unique
parents, and minimize disease and max parent use. Ties are recorded in
`tied_methods` and `tied_method_count`.

AlphaSimR crop multi-trait validation grid:

```text
Rscript tools/run_multitrait_crop_validation_grid.R
```

Useful crop-grid overrides:

```text
NG_MULTITRAIT_CROP_GRID_SCENARIOS=compact_selfing,maize_like,cassava_diploid
NG_MULTITRAIT_CROP_GRID_PARENT_SIZES=12,20
NG_MULTITRAIT_CROP_GRID_REPS=1
NG_MULTITRAIT_CROP_GRID_CROSSES=4
NG_MULTITRAIT_CROP_GRID_REALIZED_PROGENY=8
NG_MULTITRAIT_CROP_GRID_OCS_LAMBDA_GROUP=0.05
```

This runner uses the crop-genome metadata, AlphaSimR multi-trait additive
genetic values with crop-specific covariance profiles, DH base parents, and
realized DH progeny family means. It writes `*_summary.csv`,
`*_selections.csv`, `*_scores.csv`, `*_winner_summary.csv`, and
`*_config.csv`. Winner summaries keep a primary tie-broken `method` column and
also record `tied_methods` plus `tied_method_count`.

The crop runner realizes every candidate pair before selection. This is useful
for audit-scale validation and small/medium screens, but broad parent-size
grids should keep parent counts, replicate counts, and progeny per cross under
control or move to a shortlist realization design.

Formal multi-trait head-to-head benchmark:

```text
Rscript tools/run_head_to_head_benchmark.R
```

Useful overrides:

```text
NG_HEAD_TO_HEAD_SCENARIOS=compact_selfing,cassava_diploid,potato_tetraploid_stress
NG_HEAD_TO_HEAD_PARENT_SIZES=12,20
NG_HEAD_TO_HEAD_REPS=1
NG_HEAD_TO_HEAD_CROSSES=4
NG_HEAD_TO_HEAD_REALIZED_PROGENY=8
NG_HEAD_TO_HEAD_METHODS=nextgen_auto_ocs,nextgen_economic_index_ocs,nextgen_desired_gain_ocs,popvar_style_weighted_topn,simplemate_style_threshold_topn,alphamate_style_weighted_ocs
NG_HEAD_TO_HEAD_OCS_LAMBDA_GROUP=0.05
```

This runner compares the NextGen multi-trait OCS candidates against
CI-friendly PopVar-, SimpleMating-, and AlphaMate-style baselines on the same
crop scenario, parent set, realized family means, trait directions, and crossing
budget. It writes `*_summary.csv`, `*_selections.csv`, `*_scores.csv`,
`*_comparisons.csv`, `*_winner_summary.csv`, `*_config.csv`, and
`*_method_registry.csv`. The registry and summary outputs carry
`implementation`, `exact_external_status`, and `fallback_reason` fields so
style-proxy baselines are never confused with exact package or executable runs.
Use the existing external runners for exact PopVar/SimpleMating/AlphaMate jobs
when those tools are installed. `potato_tetraploid_stress` in this runner is a
diploidized crop stress approximation; use `tools/run_poly4x_benchmark.R` for
true autotetraploid 4x claims.

Breeder-facing visual report:

```text
NG_HEAD_TO_HEAD_PREFIX=head_to_head_multitrait_grid_20260507
Rscript tools/render_head_to_head_visual_report.R
```

The renderer reads the benchmark CSVs and writes
`*_visual_report.html`. The report is the first frontend boundary over the R
backend: it shows multi-trait response, gain-diversity trade-offs, parent
contribution, mate allocation, ranked crossing decisions, and evidence caveats
for style-proxy baselines.

Interactive frontend:

A point-and-click Shiny workbench for this engine is maintained separately as
the `nextgenCrossWorkbench` package:
<https://github.com/pulsesmartlab-innovations/NextGenCrossDesign>

The workbench does not link against this package. It drives it out-of-process
through the headless JSON contract in `docs/frontend/contracts/`
(`ng_run_config.v1` in, `ng_run_result.v1` out) via
`tools/run_cross_prediction_json.R`, so the frontend installs as pure R while
the compiled kernel stays here.

## Autotetraploid 4x Model

`ng_poly4x_*` is a separate model family for true autotetraploid 4x cross
design. It uses AlphaSimR 4x dosage genotypes, direct non-DH `makeCross()`
families, sampled-progeny family scoring, and AlphaMate-style constrained
allocation through the existing optimizer. It is intended first for
potato/cassava-like autotetraploid clonal crops.

This path is not a universal polyploid model. Sugarcane-like mixed or aneuploid
inheritance and wheat-like allopolyploid subgenomes remain separate modeling
problems.

The production-facing 4x policy wrapper is `ng_poly4x_policy()`. Its default
`mode = "gain"` uses usefulness top-N and is the recommended gain default from
the 10-rep diagonal recurrent screen. `mode = "diversity"` uses the stronger
OCS diversity setting and is the recommended diversity-control mode.
`mode = "ocs"` exposes the default OCS setting for explicit experiments, but it
is not auto-selected for potato p20 because the observed gain edge was weak.
`ng_poly4x_policy_select()` returns auditable policy names for these modes
without changing the older diploid/crop stress dispatcher.

Tiny smoke test:

```text
Rscript tests/poly4x_runner_smoke.R
```

Direct benchmark runner:

```text
Rscript tools/run_poly4x_benchmark.R
```

The direct runner keeps `NG_POLY4X_METHODS` for legacy methods and also accepts
`NG_POLY4X_POLICY_MODES=gain,diversity,ocs`. Policy modes are recorded as
`poly4x_policy_<mode>` methods and preserve `poly4x_policy_*` diagnostics in
selection outputs plus `pred_poly4x_policy_*` diagnostics in family outputs.

Controlled OCS comparison runner:

```text
Rscript tools/run_poly4x_controlled_ocs.R
```

This controlled runner compares decisions on a common parent population by
reusing one scored cross table and cached realized family outcomes for
overlapping selected pairs across `usefulness`, `default`, and `strongdiv`
configurations. The shared parent context advances with
`NG_POLY4X_CONTROLLED_ADVANCE_CONFIG`, which defaults to `usefulness`. This is
not an independent long-term trajectory benchmark for each config: in cycle 2
and later, non-advance configs are evaluated as decisions in the shared context
created by the configured advancement strategy.

Controlled configs still accept `usefulness`, `default`, and `strongdiv`, and
now also accept policy aliases `gain`, `diversity`, and `ocs`.
`NG_POLY4X_CONTROLLED_ADVANCE_CONFIG` can use any config listed in
`NG_POLY4X_CONTROLLED_CONFIGS`, including the policy aliases.

Small grid:

```text
powershell -ExecutionPolicy Bypass -File tools/run_poly4x_grid.ps1
```

By default the grid runs the legacy 4x methods plus policy modes
`gain`, `diversity`, and `ocs`. Override `NG_POLY4X_METHODS` and
`NG_POLY4X_POLICY_MODES` to narrow or expand the comparison.

## Key Decision

The default target here is DH/RIL cross design, not generic F1 progeny variance. For inbred-line crossing, the most relevant segregation variance is the recombination variance of the F1 haplotype contrast. The C++ kernel computes the full Haldane version of:

```text
Var(DH_ij) = a_ij' R a_ij
a_ijk = 0.5 * (x_ik - x_jk) * beta_k
R_kl = 1 - 2r_kl for loci on the same chromosome, else 0
```

The PMV extension uses `E(beta_k^2) = beta_k^2 + Var(beta_k)` on diagonal terms. Off-diagonal posterior covariances are deliberately not assumed unless supplied later by an exact shortlist scorer.

Because Haldane decay is exponential, each chromosome quadratic form is computed exactly with a linear recursion rather than a dense marker-by-marker matrix.

## Benchmark Target

Use both focused realistic runs and parent-size grids. Breeding programs may have 20, 30, 40, 50, 60, 70, 80, or more available parents, so no method should be promoted from one parent-count scenario alone.

The standard parent-size screen compares:

- `var_simple`;
- expected top-k `var_simple`;
- recombination-aware DH PMV and calibrated hybrid usefulness;
- AlphaMate-style gain-diversity mate allocation with adaptive parent-use penalties;
- state-of-art external baselines where installable.

The first correctness check is not whether a metric wins immediately. It is whether predicted within-family variance has the right direction and scale against realized simulated family variance.

## Quick Smoke Test

```r
source("R/load.R")
ng_load("nextgen_cross_design")
```

Or from the repository root:

```text
Rscript tests/smoke_test.R
```

The standalone AlphaSimR harness defaults to pure R to avoid accidental compiler failures in locked-down sessions. For realistic 5K-marker runs, use the C++ kernel after verifying Rtools/Rcpp:

```text
NG_USE_CPP=1 Rscript tools/run_alphasimr_benchmark.R
```

For fair reproducible AlphaSimR comparisons, keep `NG_ALPHASIMR_THREADS=1`
unless you are only doing exploratory runtime checks. AlphaSimR's default
thread count can change realized QTL/recombination outcomes under the same
R seed, which invalidates head-to-head method comparisons.

Shared score caching is enabled by default with `NG_SHARED_SCORING=1`. When
multiple methods start from the same parent population in a cycle, the benchmark
now fits marker effects, computes internal cross scores, and runs external
PopVar/SimpleMating shortlist rescoring once, then remaps the shared score table
back to each method's parent IDs. Set `NG_SHARED_SCORING=0` only for debugging.

The realistic benchmark wrapper configures the current OCS frontier:

```text
powershell -ExecutionPolicy Bypass -File tools/run_realistic_benchmark.ps1
```

The general-applicability screen runs multiple parent counts and summarizes the results:

```text
powershell -ExecutionPolicy Bypass -File tools/run_parent_size_grid.ps1
```

The frontier-policy validation grid can be launched directly with:

```text
powershell -ExecutionPolicy Bypass -File tools/run_frontier_policy_validation_grid.ps1
```

Real AlphaMate can also be used as an external mating-plan baseline when the
official binary is available. By default the wrapper looks for
`external/AlphaMate/binaries/AlphaMate.exe`, or use `NG_ALPHAMATE_EXE` to point
to another copy. On this Windows setup the bundled executable also needs
`libiomp5md.dll`, so set `NG_ALPHAMATE_RUNTIME_PATH` to the directory that
contains it, for example:

```text
NG_METHODS=var_simple_topn,alphamate_opt45,ng_frontier_policy_ocs10_lps2
NG_ALPHAMATE_RUNTIME_PATH=C:\Python\Lib\site-packages\torch\lib
Rscript tools/run_alphasimr_benchmark.R
```

Supported AlphaMate benchmark methods are `alphamate_opt`, `alphamate_opt30`,
`alphamate_opt45`, `alphamate_opt60`, and other numeric target-degree variants.
These parse `ModeOptTarget1` output and record `alphamate_*` diagnostics in the
selection summary.

The replicated 5K real-AlphaMate comparison wrapper is:

```text
powershell -ExecutionPolicy Bypass -File tools/run_alphamate_external_grid_5k.ps1
```

Useful environment overrides:

```text
NG_GRID_PARENT_SIZES=20,30,40,50,60,70,80
NG_GRID_EFFECT_TRAINING_MODE=same   # or min80 / min160
NG_GRID_PREFIX=nextgen_parent_grid
NG_GRID_GENOME_LENGTH_M=1.0
NG_BALANCED_DIVERSITY_COL=var_simple_cal
NG_BALANCED_DIVERSITY_WEIGHT=0.30
NG_BALANCED_PAIR_KINSHIP_WEIGHT=0.10
NG_BALANCED_AUTO_MIN_UNIQUE=0
NG_ADAPTIVE_SCORE_COLS=uc_recomb_gebv,etk_dh_recomb_var_gebv_cal,etk_dh_recomb_var_blend_cal,etk_dh_pmv_scaled_var_blend_cal,etk_var_simple_cal
NG_ADAPTIVE_TARGET=realized_top10
NG_ADAPTIVE_FALLBACK_COL=etk_var_simple_cal
NG_ADAPTIVE_FALLBACK_MAX_WEIGHT=0.35
NG_META_METHOD_MIN_HISTORY_N=10
NG_META_METHOD_RECENT_CYCLES=3
NG_META_ROUTER_FAMILIES=recomb_gebv,pmv_balanced,portfolio,popvar_uc,simple_usefa,var_simple
NG_META_ROUTER_MIN_HISTORY_N=10
NG_META_ROUTER_RECENT_CYCLES=3
NG_META_ROUTER_HISTORY_WEIGHT=0.50
NG_META_ROUTER_PRIOR_WEIGHT=0.25
NG_META_ROUTER_PLAN_WEIGHT=0.05
NG_META_ROUTER_GAIN_WEIGHT=0.10
NG_META_ROUTER_DIVERSITY_WEIGHT=0.12
NG_META_ROUTER_RELIABILITY_WEIGHT=0.03
NG_META_ROUTER_REGRET_GUARD=1
NG_META_ROUTER_REGRET_PLAN_WEIGHT=0.60
NG_META_ROUTER_REGRET_GAIN_WEIGHT=0.40
NG_META_ROUTER_REGRET_MIN_ADVANTAGE=0.75
NG_META_ROUTER_REGRET_MAX_ROUTER_PENALTY=0.40
NG_FRONTIER_POLICY_SOURCE=validated_5k_parent_grid_2026_05_03
NG_FRONTIER_POLICY_FALLBACK_METHOD=ng_meta_router_ocs10_lps2
```

The built-in frontier policy currently dispatches:

| Parent band | Primary method |
| ---: | --- |
| 1-25 | `ng_recomb_gebv_ocs10_lps2` |
| 26-35 | `ng_meta_selector_ocs10_lps2` |
| 36-45 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` |
| 46-55 | `ng_meta_router_ocs10_lps2` |
| 56-65 | `popvar_uc_ocs10_lps1`, then `simple_usefa_ocs10_lps1`, then router fallback |
| 66-75 | `ng_meta_selector_ocs10_lps2` |
| 76+ | `ng_meta_portfolio_ocs10_lps2`, then router fallback |

Override it for a crop/site validation with `NG_FRONTIER_POLICY_SPEC`, for
example:

```text
NG_FRONTIER_POLICY_SPEC=20-40=ng_meta_router_ocs10_lps2|var_simple_ocs10_lps2;41+=ng_meta_portfolio_ocs10_lps2
```

To stress-test portability across crop-like genome architectures, run:

```text
Rscript tools/run_crop_genome_scenarios.R
```

The built-in scenarios now include compact selfing, maize, wheat, barley, field
pea, potato, cassava, and sugarcane-like templates. They vary chromosome count,
genetic map length, marker density, QTL density, heritability, founders, and
effect-training size while keeping the DH/RIL diploid assumption. Wheat,
potato, cassava tetraploid, and sugarcane entries are marked as diploidized
stress approximations; this is a portability screen, not a true polyploid or
crop-specific biological simulator. Key overrides:

```text
NG_CROP_GRID_SCENARIOS=compact_selfing,maize_like,bread_wheat_hexaploid_approx,barley_like,field_pea_like,potato_tetraploid_stress,cassava_diploid,cassava_tetraploid_stress,sugarcane_polyploid_stress
NG_CROP_GRID_PARENT_SIZES=20,60,80
NG_CROP_GRID_REPS=1
NG_CROP_GRID_CYCLES=2
NG_CROP_GRID_PREFIX=crop_genome_frontier_grid
```

To run the same crop screen with real AlphaMate target-degree controls:

```text
powershell -ExecutionPolicy Bypass -File tools/run_crop_genome_alphamate_grid.ps1
```

Both crop wrappers include `ng_crop_aware_policy_ocs10_lps2` by default. In the
current smoke policy it delegates `cassava_tetraploid_stress` to
`alphamate_opt60`, `sugarcane_polyploid_stress` to `alphamate_opt45`, and other
templates to `ng_frontier_policy_ocs10_lps2` before trying AlphaMate fallbacks.

External PopVar and SimpleMating baselines can be run with:

```text
powershell -ExecutionPolicy Bypass -File tools/run_external_baseline_benchmark.ps1
```

The default exact external screen uses 40 parents and 5K markers. Full 80-parent all-pair PopVar/SimpleMating usefulness runs can be much slower, so they should be treated as large validation jobs rather than default smoke tests.

The replicated exact external parent-size screen runs 20, 30, and 40 parents by default:

```text
powershell -ExecutionPolicy Bypass -File tools/run_external_parent_size_grid.ps1
```

For 50 to 80+ parents, use the shortlist/exact-rescore screen. It scores all
pairs with the fast internal metrics, then runs PopVar and SimpleMating only on
the strongest shortlist. The default shortlist is the union of top pairs from
calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV so
the external methods are not restricted to one internal score ranking. The
default is a one-replicate exploratory screen; raise `NG_EXTERNAL_GRID_REPS` or
`NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER` for deeper validation:

```text
powershell -ExecutionPolicy Bypass -File tools/run_external_shortlist_parent_size_grid.ps1
```

To diagnose whether a result comes from the score itself or from the mate
allocation algorithm, run the allocator crosscheck grid. It compares top-N,
adaptive OCS, and SimpleMating-style constrained selection across the internal,
PopVar, and SimpleMating score families:

```text
powershell -ExecutionPolicy Bypass -File tools/run_allocator_crosscheck_grid.ps1
```

Before promoting any metric, run the family-level calibration benchmark. It
scores the same sampled crosses with the internal metrics plus PopVar,
SimpleMating, and the genomicMateSelectR-derived F1/DH variance check, then
compares predictions against realized AlphaSimR families:

```text
powershell -ExecutionPolicy Bypass -File tools/run_family_calibration_grid.ps1
```

The diagnostic-first validation protocol is documented in
`VALIDATION_PROTOCOL.md`. The wrapper below defaults to
smoke checks; set `NG_VALIDATION_PHASE` to `family`, `allocator`, `grid`, or
`all` for the heavier validation phases:

```text
powershell -ExecutionPolicy Bypass -File tools/run_framework_validation.ps1
```

After a diagnostic run, generate the markdown summary report with:

```text
Rscript tools/summarize_validation_report.R
```
