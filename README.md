# nextgenCrossDesign

Genomic cross prediction and mate allocation for breeding programs.

Version 0.4.0.

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

## Quick Smoke Test

```r
source("R/load.R")
ng_load()
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

The validation and benchmark harness (parent-size grids, frontier-policy and
crop-portability screens, external PopVar/SimpleMating/AlphaMate baselines, and
the family-calibration diagnostics) is documented in `VALIDATION_PROTOCOL.md`
and `BENCHMARK_NOTES.md`.
