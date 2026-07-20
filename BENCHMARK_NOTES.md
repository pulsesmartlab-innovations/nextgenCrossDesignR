# Benchmark Notes

## COMA vs native mate allocation, outbred+dominance recurrent GS (2026-07-03)

Runner: `tools/run_mate_allocation_study.R` (recurrent mode). Outbred diploid AND autotetraploid,
runMacs founders (realistic LD), addTraitAD dominance. Each method = an independent recurrent
GS lineage from shared founders: allocate on GEBV (augmented additive+dominance ridge;
simulation-based predicted GPMP + within-family variance), evaluate on true gv; common GEBV
parent selection so only the ALLOCATOR differs. Methods: `coma_oma` (COMA exact convex OMA,
installed jendelman/COMA), `ng_useful` (mean + i*SD), `ng_mean`. 200 founders -> 180 pop,
14 parents, 12 crosses, family 18, 6 chr / ~720 SNP / 90 QTL, h2 0.4, dF 0.05.

Final-cycle realized gain (mean +/- SE), by config:
- A (8 cyc, meanDD 0.5, 10 reps): diploid coma 5.48+/-.16 ~ useful 5.47+/-.24 ~ mean 5.39+/-.18;
  4x mean 6.52+/-.30 ~ useful 6.50+/-.35 ~ coma 6.32+/-.46. All within SE. `ng_useful` retains
  most variance (paired +0.037 var vs coma diploid; +0.058 var & +0.18 gain vs coma 4x).
- B (15 cyc, meanDD 0.5, 8 reps): diploid coma 5.99+/-.31 ~ useful 5.82+/-.27 > mean 5.52+/-.38;
  4x useful 8.13+/-.44 ~ coma 7.93+/-.52 ~ mean 7.94+/-.45. `ng_useful` highest var AND lowest
  coancestry both ploidies; variance edge persists all 15 cycles but does NOT become a gain lead.
- C (8 cyc, STRONG dominance meanDD 1.2, 10 reps): diploid mean 3.42+/-.22 > coma 2.97+/-.24 =
  useful 2.97+/-.25 (~2 SE); 4x useful 4.20+/-.18 ~ mean 4.18+/-.28 > coma 3.98+/-.29. COMA does
  NOT gain an edge on its heterosis home turf.

Takeaways: PARITY on gain (native greedy ~ COMA convex, within SE; native leads under strong
dominance diploid); `ng_useful` robustly retains more variance at lower coancestry but the
retained variance does not convert to a gain lead by cycle 15. Caveats: COMA's continuous OMA is
discretized to fixed equal-family n crosses (not its variable-family strength); augmented-ridge
dominance estimation is weak; approximate coancestry matching; single architecture. See
VALIDATED_STATE.md "COMA Mate-Allocation Study Note" for allowed/disallowed claims.

2x2 follow-up (usefulness-aware OMA): added `coma_useful` = our within-family-variance usefulness
fed as the per-mating merit into COMA's exact convex ΔF solve (a criterion COMA's package does not
offer -- it optimizes only the cross mean), giving {mean, usefulness} x {greedy, COMA convex}.
Final-cycle gain, main 10 reps x 10 cyc | long 8 reps x 18 cyc:
  diploid  : coma_oma 5.72|6.24  coma_useful 5.67|6.11  ng_mean 5.51|6.37  ng_useful 5.34|5.70
  tetraploid: coma_oma 7.07|8.28  ng_useful 6.60|8.65  coma_useful 6.51|8.17  ng_mean 6.32|8.62
`coma_useful` does NOT beat `coma_oma` on gain in either ploidy/horizon (behind, while retaining
more variance). Diploid paired coma_useful-coma_oma over cyc 6/12/18: gain -0.33/-0.29/-0.13,
var +0.03/+0.08/+0.12 -- gain gap shrinks as variance banks, a directional hint of a crossover
beyond 18 cycles (unconfirmed). VERDICT: fusing our variance into COMA's OMA is a diversity lever,
NOT a gain improvement over COMA's mean OMA at tested horizons. The mean criterion is at least as
good for gain; do not claim the usefulness-aware OMA as an improvement.

VARIABLE-FAMILY-SIZE follow-up (the realistic setting -- family size is never equal in practice):
NG_MATEALLOC_FAMILY_MODE=variable gives each method a fixed total progeny budget split across
matings -- COMA by its optimal contributions (ng_select_coma total_progeny=, its NATIVE unequal-
family mode; ng_select_coma gained a total_progeny arg + largest-remainder progeny allocation),
native by OCS selection + ng_allocate_family_sizes(score_weighted). 6 reps x 6 cyc, diploid+4x,
16 parents / portfolio<=20 / total_progeny 160. Final-cycle gain | var | coancestry:
  diploid  : ng_useful 4.83|.27|.016  ng_mean 4.72|.26|.017  coma_useful 4.54|.23|.031  coma_oma 4.54|.19|.031
  tetraploid: ng_mean 4.71|.34|.013  ng_useful 4.64|.43|.010  coma_useful 4.51|.35|.048  coma_oma 4.46|.37|.050
Native is STRICTLY BETTER than COMA even here: more gain at LOWER coancestry and more retained
variance, both ploidies (paired COMA-native: diploid -0.18..-0.29 gain, +0.014 coancestry;
tetraploid -0.18..-0.24 gain, +0.038 coancestry). COMA runs near its ΔF ceiling (~0.05), spending
inbreeding for one-generation contribution-optimality; native stays diverse (~0.015) and out-gains
it over cycles. (A 2-rep pilot spuriously showed COMA ahead; 6 reps reversed it.) So across BOTH
equal-family and variable-family settings COMA does not beat the native allocator on gain. Note:
the 24-parent variable run was killed for cost (per-cross progeny simulation scales badly); this
16-parent run is directional (6 reps), single config.

## Metric variants in recurrent selection, 11-method run (2026-07-03)

Runner: `tools/run_ril_breeding_program_benchmark.R` (10 reps x 25 cycles; 100 founders ->
60 parents, 40 crosses, 10 chr, ~2500 SNP / 300 QTL, 500 training lines; seed 2025). All 11
methods completed all cycles (0 skips). This runs the full cross-scoring metric family in the
RECURRENT loop -- the multi-cycle counterpart to the single-generation metric-merit ranking
study (below). Pilot first (3 reps x 8 cycles) confirmed the setup and clustering.

Cycle-25 cumulative gain (mean +/- SE, n=10), ranked:

| method                | gain          | He    | polymorphic | acc(top10) |
|-----------------------|---------------|-------|-------------|------------|
| strategy_balanced     | 11.63 +/- 0.20 | 0.044 | 0.123 | 0.871 |
| var_simple_ocs        | 11.17 +/- 0.25 | 0.065 | 0.161 | 0.867 |
| var_complex_evolution | 10.48 +/- 0.16 | 0.014 | 0.038 | 0.866 |
| package_user_api      | 10.35 +/- 0.28 | 0.023 | 0.060 | 0.858 |
| pmv_ocs               | 10.20 +/- 0.27 | 0.014 | 0.037 | 0.861 |
| var_complex_ocs       | 10.08 +/- 0.25 | 0.014 | 0.039 | 0.859 |
| uc_pmv_ocs            | 10.07 +/- 0.35 | 0.015 | 0.040 | 0.843 |
| uc_vpm_ocs            | 10.04 +/- 0.20 | 0.011 | 0.030 | 0.865 |
| vpm_ocs               |  9.56 +/- 0.21 | 0.012 | 0.033 | 0.865 |
| simple_usefa_select   |  9.45 +/- 0.27 | 0.002 | 0.006 | 0.854 |
| mean_ocs              |  9.19 +/- 0.43 | 0.004 | 0.013 | 0.856 |

HEADLINE -- single-generation ranking != recurrent realized gain: `var_simple`, the WEAKEST
metric in the single-generation study (ranking accuracy 0.443, last), is the SECOND-BEST
recurrent performer (11.17, ~4 SE above the usefulness cluster) and retains by far the most
polymorphism (0.161). It is a relationship-DISTANCE criterion favouring unrelated/complementary
parents -> preserves variance -> sustains long-term gain, despite being a poor single-generation
MERIT predictor. The PMV/VPM/uc usefulness variants cluster tightly (~9.6-10.2, overlapping SE
-- they rank crosses nearly identically) and deplete diversity hard (He ~0.01, poly ~0.03),
plateauing lower. mean/SimpleMating lowest, most depleted. Accuracy ~0.84-0.87 across all --
the metric differences are about variance management over cycles, not prediction quality. So:
short-term single-gen selection -> PMV/VPM usefulness (var_simple weak); long-term recurrent ->
explicit diversity management on a good merit metric (strategy dial 11.63), not var_simple's
incidental diversity effect nor a pure-usefulness metric. Caveat: single config, 10 reps,
directional; absolute gains run-specific (see cross-run note in the 6-method entry).

## Metric-merit study (2026-07-03)

Runner: `tools/run_metric_merit_study.R` (8 conditions x 15 reps = 120 runs, 0 failures).
Design: AlphaSimR `runMacs` GENERIC founders (quickHaplo gives ~0 marker-QTL LD -> unusable),
DH lines, RR-BLUP effects on a 200-SNP/chr chip; candidate crosses among 40 parents. Judge =
unbiased true-effect yardstick (each cross re-scored with the TRUE QTL effects via
`ng_dh_recomb_variance_pairs`, true usefulness = mid-parent true GV + i * true within-family
SD). Measures: Spearman ranking accuracy vs true usefulness (allocation-free) and top-K true
merit. Conditions: training size {120, 400} x h2 {0.25, 0.60} x architecture
{oligogenic 5 QTL/chr, polygenic 40 QTL/chr}.

Overall ranking accuracy (mean +/- SE across all conditions):

| metric              | rank accuracy |
|---------------------|---------------|
| UC-VPM = SimpleMating | 0.471 +/- 0.016 |
| UC-PMV (default)      | 0.469 +/- 0.015 |
| mean (gain-only)      | 0.455 +/- 0.014 |
| UC-var_simple         | 0.443 +/- 0.016 |

Usefulness advantage over mean-only (UC-PMV minus mean, paired within rep, by condition):

| condition                       | UC-PMV - mean      |
|---------------------------------|--------------------|
| oligogenic, train=400, h2=0.60  | +0.058 +/- 0.006   |
| oligogenic, train=400, h2=0.25  | +0.025 +/- 0.005   |
| oligogenic, train=120, h2=0.60  | +0.019 +/- 0.006   |
| polygenic,  train=400, h2=0.60  | +0.019 +/- 0.007   |
| polygenic,  train=400, h2=0.25  | +0.009 +/- 0.005   |
| polygenic,  train=120, h2=0.25  | +0.006 +/- 0.010   |
| polygenic,  train=120, h2=0.60  | -0.003 +/- 0.004   |
| oligogenic, train=120, h2=0.25  | -0.020 +/- 0.013   |

Findings: (1) UC-VPM and SimpleMating additive usefulness are IDENTICAL in ranking
(max |diff| 1e-4 across all 120 runs) -- independent cross-validation. (2) UC-PMV ~= UC-VPM
in ranking (PMV's posterior term shifts magnitudes, not order; its value is calibration).
(3) Usefulness beats mean-only only when the within-family variance is both meaningful and
well-estimated (oligogenic + adequate training + heritability); it is ~= mean for polygenic
traits and can slightly HURT when effects are poorly estimated (oligogenic + small training +
low h2). (4) var_simple is the weakest merit metric -- a diversity/QC proxy, not a cross-merit
criterion. Recommendation: default UC-PMV; use mean-only for highly polygenic or small/low-h2
training; var_simple/mean for diversity/baseline only. Caveats: single-generation ranking on
DH progeny, one demography; not multi-cycle realized gain.

## RIL recurrent study, 6 methods incl. user-API (2026-07-03)

Runner: `tools/run_ril_breeding_program_benchmark.R` (10 reps x 25 cycles; 100 founders ->
60 parents, 40 crosses, 10 chr, ~2500 SNP / 300 QTL, 500 training lines; seed 2025;
RGL_USE_NULL for SimpleMating). Six methods, all completed all cycles (0 skips):
`var_complex_ocs` (internal calls), `var_complex_evolution`, `strategy_balanced`, `mean_ocs`,
`simple_usefa_select` (SimpleMating package API), and `package_user_api` -- the package driven
through its TOP-LEVEL user API (`ng_run_cross_prediction`, data frames in / selected crosses
out), the same black-box way SimpleMating is invoked.

Cycle-25 cumulative gain (mean +/- SE, n=10):

| method                     | gain          | He   | acc(top10) |
|----------------------------|---------------|------|------------|
| strategy_balanced          | 11.60 +/- 0.22 | 0.05 | 0.869 |
| var_complex_ocs (internal) | 10.10 +/- 0.19 | 0.01 | 0.873 |
| package_user_api           | 10.10 +/- 0.12 | 0.02 | 0.847 |
| var_complex_evolution      |  9.89 +/- 0.22 | 0.01 | 0.859 |
| mean_ocs                   |  9.49 +/- 0.24 | 0.01 | 0.865 |
| simple_usefa_select (RIL)  |  9.26 +/- 0.23 | 0.00 | 0.857 |

USER-API VALIDATION: `package_user_api` reproduces the internal `var_complex_ocs` path on
long-term gain -- both land at 10.10 by cycle 25, with mean absolute per-cycle difference 0.57
(~5% of the ~10 scale), concentrated in early cycles where they converge (by-cycle user-API:
2.02 / 6.92 / 9.39 / 10.0 / 10.1 / 10.1 at c1/5/10/15/20/25 vs internal 2.28 / 7.43 / 9.64 /
10.1 / 10.1 / 10.1). The small early gap is honest: the user API trains effects on the
candidate parents' own phenotypes (what a user has) and scores recurrent RIL parents with
`assume_inbred = FALSE`. So the package used as a black-box user API matches the internal
result -- the same package-to-package basis on which SimpleMating is compared.

Other findings: the strategy dial sustains the highest long-term gain by preserving variance
(fixed-lambda methods plateau by ~cycle 15). Evolution and OCS are WITHIN NOISE on realized
gain (this run OCS 10.10 vs evolution 9.89; the earlier 5-method run had evolution 9.92 vs OCS
9.66 -- ordering flips within ~1 SE; do not claim evolution superiority on realized gain).
SimpleMating plateaus lowest under default culling. NOTE: SimpleMating was fixed in commit
aabe2f0 (Type="RIL", unfixed RIL loci coded NA as getUsefA requires 0/2/NA); an earlier
Type="DH" run had it dying at cycle 2 (invalid column). Caveats: single config, 10 reps,
untuned OCS/SimpleMating diversity settings -> directional, not a delta-F-matched comparison
(use `target_coancestry`). Cross-run note: absolute gains shifted ~+0.4 vs the earlier 5-method
run for shared methods, so treat absolute values as run-specific; relative orderings (strategy
top, SimpleMating lowest, user-API == internal) are the takeaways.

## RIL recurrent study, first 5-method run (2026-07-03, superseded)

(Superseded by the 6-method run above; retained for provenance.) 5 methods
(`var_complex_ocs`, `var_complex_evolution`, `strategy_balanced`, `mean_ocs`,
`simple_usefa_select`). Cycle-25 gain: strategy_balanced 11.32 +/- 0.23 > evolution
9.92 +/- 0.23 > OCS 9.66 +/- 0.21 > mean 9.02 +/- 0.22 > SimpleMating 8.79 +/- 0.36.

## v0.1.0 metric break (2026-05-22)

The metric layer was corrected in v0.1.0; see `docs/V0_1_0_MIGRATION.md` for
the full list of fixed bugs and removed columns. Every benchmark run logged
below was produced with one or more of the following statistical bugs in the
pipeline:

- Posterior marker-effect variance `beta_var` understated under LD (T1.3);
  PMV inherited the bias and ranked high-`dh_pmv_var` crosses too favorably.
- `pmv_scale`-rescaled PMV mixed (genetic value)² units with the `var_simple`
  relatedness-distance metric; `uc_dh_scaled`, `uc_gated`, `uc_hybrid` (and
  `_gebv`/`_adj`/`_blend` variants) were dimensionally invalid (T1.4).
- `recomb_model = "kosambi"` was silently treated as Haldane in both the R
  recursion and the C++ kernel (T1.1).
- `target = "RIL"` was silently treated as DH (T1.2).
- Multi-trait `economic_index` and `desired_gain` solved against sample
  Pearson covariance of cross predictions, not Smith-Hazel `P` / Pesek-Baker
  `G`, and applied the ridge penalty twice (T1.6).

Runs below should be repeated on v0.1.0 before being cited externally. Method
names that were removed (e.g. `ng_pmv_scaled_blend_balanced_ocs10_lps2`,
`ng_recomb_rel_*`, `ng_hybrid_*`) have been renamed or deleted. New methods
to use as the in-units replacements:
- `ng_pmv_blend_balanced_ocs10_lps2` instead of `ng_pmv_scaled_blend_balanced_ocs10_lps2`
- `uc_dh_gebv` / `uc_dh_blend` instead of `uc_gated` / `uc_hybrid`
- `etk_dh_pmv_var_*` instead of `etk_dh_pmv_scaled_var_*`

---

## Run: `multitrait_validation`

Date: 2026-05-06.

This is the first deterministic smoke harness for the multi-trait objective
layer. It is not a replicated crop validation grid. The purpose is to make the
multi-trait methods runnable and auditable on a controlled trade-off scenario
before scaling to AlphaSimR or crop-specific validation.

Command:

```text
Rscript nextgen_cross_design/tools/run_multitrait_validation.R
```

Default configuration:

- Parent count: 20.
- Selected crosses: 6.
- Seed: 1.
- Allocator: top-N.
- Methods: `auto`, `weighted`, `economic_index`, `desired_gain`, and
  `threshold`.
- Traits: increase yield, decrease disease, and increase quality.

Observed smoke output:

| Method | Mean realized index | Mean yield | Mean disease | Mean quality | Unique parents | Max parent use |
|---|---:|---:|---:|---:|---:|---:|
| `auto` | 1.643 | 96.677 | 31.289 | 41.185 | 7 | 6 |
| `weighted` | 1.643 | 96.677 | 31.289 | 41.185 | 7 | 6 |
| `economic_index` | 1.256 | 96.283 | 34.713 | 39.783 | 9 | 4 |
| `desired_gain` | 1.116 | 97.675 | 36.361 | 39.133 | 9 | 4 |
| `threshold` | 1.643 | 96.677 | 31.289 | 41.185 | 7 | 6 |

Interpretation:

- The smoke harness is working and produces method-level realized summaries.
- Covariance-aware methods produced different parent-use behavior in this
  scenario, with more unique parents and lower max parent use.
- This is not evidence of multi-trait superiority. It is a deterministic
  contract and a launch point for a replicated multi-trait validation grid.

## Run: `multitrait_grid_smoke_20260506`

Date: 2026-05-06.

This is the first replicated synthetic grid for the multi-trait objective
layer. It repeats the mixed-direction validation scenario across parent sizes
and replicates, then summarizes which method wins each realized metric after
averaging replicates. It is still synthetic validation, not a crop-specific
AlphaSimR multi-trait breeding benchmark.

Command:

```text
NG_MULTITRAIT_GRID_PREFIX=multitrait_grid_smoke_20260506
NG_MULTITRAIT_GRID_PARENT_SIZES=8,10
NG_MULTITRAIT_GRID_REPS=2
NG_MULTITRAIT_GRID_CROSSES=3
Rscript nextgen_cross_design/tools/run_multitrait_validation_grid.R
```

Configuration:

- Parent sizes: 8 and 10.
- Replicates: 2.
- Selected crosses: 3.
- Allocator: top-N.
- Methods: `auto`, `weighted`, `economic_index`, `desired_gain`, and
  `threshold`.
- Output prefix: `multitrait_grid_smoke_20260506`.

Winner summary:

| Parents | Metric | Direction | Winning method | Value |
| ---: | --- | --- | --- | ---: |
| 8 | Mean realized index | Maximize | `economic_index` | 0.864 |
| 8 | Mean realized yield | Maximize | `desired_gain` | 99.808 |
| 8 | Mean realized disease | Minimize | `auto` | 39.834 |
| 8 | Mean realized quality | Maximize | `economic_index` | 37.296 |
| 8 | Unique parents | Maximize | `economic_index` | 5.000 |
| 8 | Max parent use | Minimize | `desired_gain` | 2.000 |
| 10 | Mean realized index | Maximize | `auto` | 1.106 |
| 10 | Mean realized yield | Maximize | `desired_gain` | 95.609 |
| 10 | Mean realized disease | Minimize | `auto` | 32.137 |
| 10 | Mean realized quality | Maximize | `auto` | 38.612 |
| 10 | Unique parents | Maximize | `desired_gain` | 4.500 |
| 10 | Max parent use | Minimize | `auto` | 2.000 |

Interpretation:

- The replicated grid harness is runnable and records per-parent-size winners
  by realized breeding objective metric.
- No one multi-trait method dominated all metrics, which is expected when
  breeding programs balance gain, disease reduction, quality, and parent-use
  constraints.
- This strengthens workflow readiness, not external superiority. The next
  scientific validation should run crop-specific multi-trait AlphaSimR scenarios
  with realistic trait covariance and economic weights.

## Run: `multitrait_crop_grid_smoke_20260506`

Date: 2026-05-06.

This is the first AlphaSimR-backed crop multi-trait smoke grid. It uses
crop-genome scenario metadata, crop-specific additive genetic covariance
profiles, multi-trait selection methods, and realized DH progeny family means.
The run used reduced genome settings so it is a contract and wiring check, not
a production-scale crop benchmark.

Command:

```text
NG_MULTITRAIT_CROP_GRID_PREFIX=multitrait_crop_grid_smoke_20260506
NG_MULTITRAIT_CROP_GRID_SCENARIOS=compact_selfing,cassava_diploid
NG_MULTITRAIT_CROP_GRID_PARENT_SIZES=6,8
NG_MULTITRAIT_CROP_GRID_REPS=1
NG_MULTITRAIT_CROP_GRID_CROSSES=2
NG_MULTITRAIT_CROP_GRID_REALIZED_PROGENY=3
NG_MULTITRAIT_CROP_GRID_N_FOUNDERS=10
NG_MULTITRAIT_CROP_GRID_N_CHR=2
NG_MULTITRAIT_CROP_GRID_SEG_SITES=30
NG_MULTITRAIT_CROP_GRID_SNP_PER_CHR=8
NG_MULTITRAIT_CROP_GRID_QTL_PER_CHR=2
Rscript nextgen_cross_design/tools/run_multitrait_crop_validation_grid.R
```

Configuration:

- Scenarios: `compact_selfing` and `cassava_diploid`.
- Parent sizes: 6 and 8.
- One replicate.
- Selected crosses: 2.
- Realized DH progeny per cross: 3.
- Methods: `auto`, `weighted`, `economic_index`, `desired_gain`, and
  `threshold`.
- Trait directions: increase yield, decrease disease, and increase quality.

Winner summary:

| Scenario | Parents | Metric | Direction | Winning method | Value |
| --- | ---: | --- | --- | --- | ---: |
| `compact_selfing` | 6 | Mean realized index | Maximize | `auto` | 1.100 |
| `compact_selfing` | 6 | Mean realized yield | Maximize | `desired_gain` | 102.157 |
| `compact_selfing` | 6 | Mean realized disease | Minimize | `auto` | 29.685 |
| `compact_selfing` | 6 | Mean realized quality | Maximize | `auto` | 37.025 |
| `compact_selfing` | 8 | Mean realized index | Maximize | `auto` | 1.004 |
| `compact_selfing` | 8 | Mean realized yield | Maximize | `economic_index` | 104.615 |
| `compact_selfing` | 8 | Mean realized disease | Minimize | `auto` | 20.380 |
| `compact_selfing` | 8 | Mean realized quality | Maximize | `economic_index` | 36.909 |
| `cassava_diploid` | 6 | Mean realized index | Maximize | `auto` | 0.230 |
| `cassava_diploid` | 6 | Mean realized yield | Maximize | `auto` | 101.334 |
| `cassava_diploid` | 6 | Mean realized disease | Minimize | `auto` | 27.684 |
| `cassava_diploid` | 6 | Mean realized quality | Maximize | `desired_gain` | 37.088 |
| `cassava_diploid` | 8 | Mean realized index | Maximize | `economic_index` | 0.916 |
| `cassava_diploid` | 8 | Mean realized yield | Maximize | `auto` | 115.361 |
| `cassava_diploid` | 8 | Mean realized disease | Minimize | `economic_index` | 24.540 |
| `cassava_diploid` | 8 | Mean realized quality | Maximize | `economic_index` | 33.484 |

Interpretation:

- The crop multi-trait AlphaSimR harness is now runnable end to end using DH
  base parents. The original smoke prefix was rerun after the DH-parent fix.
- The winning method changes by crop scenario, parent size, and metric. This is
  the desired diagnostic behavior for a multi-objective breeding tool.
- This is smoke evidence only. The next publishable run should use real crop
  genome settings, more parents, more replicates, more progeny per cross, and
  possibly external allocation baselines.

## Run: `multitrait_crop_grid_extensive_20260506`

Date: 2026-05-06.

This extensive pass was run after code review found that the first crop harness
used heterozygous founder parents while documenting DH/RIL behavior. The harness
now makes DH base parents, passes marker-derived `parent_K` into OCS allocation,
uses a configurable OCS group-coancestry penalty, uses scenario-order-stable
seeds, writes candidate score tables, and records tie-aware winner summaries.

Baseline verification:

- Full R test sweep over `nextgen_cross_design/tests/*.R`: passed.
- Focused tests passed for synthetic multi-trait grids, crop multi-trait grids,
  and package-root test invocation.

Reduced-genome expanded grid:

```text
NG_MULTITRAIT_CROP_GRID_PREFIX=multitrait_crop_grid_extensive_20260506
NG_MULTITRAIT_CROP_GRID_SCENARIOS=compact_selfing,maize_like,barley_like,cassava_diploid
NG_MULTITRAIT_CROP_GRID_PARENT_SIZES=8,12,16
NG_MULTITRAIT_CROP_GRID_REPS=2
NG_MULTITRAIT_CROP_GRID_CROSSES=3
NG_MULTITRAIT_CROP_GRID_REALIZED_PROGENY=6
NG_MULTITRAIT_CROP_GRID_N_FOUNDERS=24
NG_MULTITRAIT_CROP_GRID_N_CHR=4
NG_MULTITRAIT_CROP_GRID_SEG_SITES=80
NG_MULTITRAIT_CROP_GRID_SNP_PER_CHR=20
NG_MULTITRAIT_CROP_GRID_QTL_PER_CHR=5
Rscript nextgen_cross_design/tools/run_multitrait_crop_validation_grid.R
```

Output checks:

| Output | Rows | Notes |
| --- | ---: | --- |
| `*_summary.csv` | 120 | 4 scenarios x 3 parent sizes x 2 reps x 5 methods; no non-finite numeric summary values |
| `*_selections.csv` | 360 | 3 selected crosses per summary row |
| `*_scores.csv` | 1712 | all candidate pairs with predicted and realized trait means |
| `*_winner_summary.csv` | 72 | 4 scenarios x 3 parent sizes x 6 metrics; no non-finite numeric values |
| `*_config.csv` | 24 | includes parent generation, methods, `n_crosses`, and AlphaSimR thread count |

Primary winner counts in `winner_summary`:

| Method | Winner rows |
| --- | ---: |
| `auto` | 41 |
| `desired_gain` | 17 |
| `economic_index` | 14 |

Tie-size distribution:

| Tied method count | Winner rows |
| ---: | ---: |
| 1 | 25 |
| 2 | 6 |
| 3 | 30 |
| 4 | 4 |
| 5 | 7 |

Interpretation:

- The extensive reduced-genome grid is internally consistent and auditable:
  selected crosses, full candidate scores, config, and tie-aware winners are all
  written.
- The primary winner changes by metric. `auto` is the most frequent primary
  winner for realized index, disease reduction, quality, and unique parents;
  `economic_index` is most frequent for realized yield; no method dominates all
  objectives.
- Many cells tie exactly because `auto`, `weighted`, and soft-threshold scoring
  can resolve to the same selection when supplied weights and thresholds do not
  change the ranking. The `tied_methods` column should be used when interpreting
  winner counts.

Metadata-scale crop grid:

- Prefix: `multitrait_crop_grid_metadata_20260506`.
- Scenarios: `compact_selfing`, `maize_like`, and `cassava_diploid`.
- Parent sizes: 12 and 20.
- One replicate, 4 selected crosses, 8 realized DH progeny per cross.
- Output rows: 30 summary, 120 selections, 768 full candidate scores, 36 winner
  rows, and 6 config rows.
- Primary winners: `auto` 21 rows, `economic_index` 9 rows, `desired_gain` 6
  rows.

OCS allocator check:

- Prefix: `multitrait_crop_grid_ocs_smoke_20260506`.
- Scenarios: `compact_selfing` and `cassava_diploid`, 8 parents, one replicate.
- The OCS path completed with marker-derived parent relationships and finite
  `group_coancestry` summaries using the default `ocs_lambda_group = 0.05`.
- Output rows: 6 summary, 18 selections, 56 full candidate scores, 12 winner
  rows, and 2 config rows.

Remaining limitations:

- This is still reduced-scale validation, not a production superiority claim.
- The crop runner realizes every candidate pair before selection. That is useful
  for audit-scale testing, but broad realistic grids should use a shortlist
  realization design before raising parent counts, replicate counts, or progeny
  per cross.

## Run: `head_to_head_multitrait_smoke_20260506`

Date: 2026-05-06.

This is the first formal multi-trait head-to-head harness. It compares NextGen
multi-trait OCS candidates with PopVar-, SimpleMating-, and AlphaMate-style
baselines on the same crop scenario, parent set, realized family means, trait
directions, and crossing budget. The external-tool rows in this smoke are
explicit style proxies for CI, not exact package or executable runs.

Command:

```text
$env:NG_HEAD_TO_HEAD_PREFIX = "head_to_head_multitrait_smoke_20260506"
$env:NG_HEAD_TO_HEAD_SCENARIOS = "compact_selfing,cassava_diploid,potato_tetraploid_stress"
$env:NG_HEAD_TO_HEAD_PARENT_SIZES = "6"
$env:NG_HEAD_TO_HEAD_REPS = "1"
$env:NG_HEAD_TO_HEAD_CROSSES = "2"
$env:NG_HEAD_TO_HEAD_REALIZED_PROGENY = "3"
$env:NG_HEAD_TO_HEAD_SEED = "101"
$env:NG_HEAD_TO_HEAD_N_FOUNDERS = "10"
$env:NG_HEAD_TO_HEAD_N_CHR = "2"
$env:NG_HEAD_TO_HEAD_SEG_SITES = "30"
$env:NG_HEAD_TO_HEAD_SNP_PER_CHR = "8"
$env:NG_HEAD_TO_HEAD_QTL_PER_CHR = "2"
Rscript nextgen_cross_design/tools/run_head_to_head_benchmark.R
```

Methods:

- NextGen candidates: `nextgen_auto_ocs`, `nextgen_economic_index_ocs`, and
  `nextgen_desired_gain_ocs`.
- Style-proxy baselines: `popvar_style_weighted_topn`,
  `simplemate_style_threshold_topn`, and `alphamate_style_weighted_ocs`.
- The method registry records `implementation`, `exact_external_status`, and
  `fallback_reason` so style proxies are not mistaken for exact PopVar,
  SimpleMating, or AlphaMate outputs.

Output checks:

| Output | Rows | Notes |
| --- | ---: | --- |
| `*_summary.csv` | 18 | 3 scenarios x 1 parent size x 1 rep x 6 methods |
| `*_selections.csv` | 36 | 2 selected crosses per summary row |
| `*_scores.csv` | 45 | 3 scenarios x all 15 candidate pairs at 6 parents |
| `*_comparisons.csv` | 198 | candidate rows compared against each baseline across metrics |
| `*_winner_summary.csv` | 18 | 3 scenarios x 6 realized metrics |
| `*_config.csv` | 3 | one row per scenario/parent-size/rep case |
| `*_method_registry.csv` | 6 | candidate/baseline labels plus exact-external caveats |

Primary winner counts in `winner_summary`:

| Method | Winner rows |
| --- | ---: |
| `alphamate_style_weighted_ocs` | 14 |
| `nextgen_desired_gain_ocs` | 3 |
| `nextgen_economic_index_ocs` | 1 |

Tie-size distribution:

| Tied method count | Winner rows |
| ---: | ---: |
| 1 | 2 |
| 2 | 2 |
| 4 | 4 |
| 6 | 10 |

Interpretation:

- The formal head-to-head harness is now runnable and produces summary,
  selection, score, comparison, winner, config, and method-registry artifacts.
- The smoke is intentionally tiny and contains many exact ties because some
  scalarized methods select the same two crosses under reduced settings.
- `potato_tetraploid_stress` here is a diploidized crop stress approximation.
  True autotetraploid claims still require the separate `ng_poly4x_*` runner.
- This does not prove superiority over PopVar, SimpleMating, AlphaMate, or
  AlphaMate-style polyploid workflows. It creates the auditable comparison
  harness needed to run larger style-proxy grids and optional exact-external
  jobs when those tools are installed.

## Run: `head_to_head_multitrait_grid_20260507`

Date: 2026-05-07.

This is the first replicated formal multi-trait head-to-head grid using the new
runner. It keeps the same CI-friendly style-proxy baseline design as the smoke
test, but expands to three scenarios, three parent sizes, and two replicates.

Command:

```text
$env:NG_HEAD_TO_HEAD_PREFIX = "head_to_head_multitrait_grid_20260507"
$env:NG_HEAD_TO_HEAD_SCENARIOS = "compact_selfing,cassava_diploid,potato_tetraploid_stress"
$env:NG_HEAD_TO_HEAD_PARENT_SIZES = "8,12,16"
$env:NG_HEAD_TO_HEAD_REPS = "2"
$env:NG_HEAD_TO_HEAD_CROSSES = "3"
$env:NG_HEAD_TO_HEAD_REALIZED_PROGENY = "6"
$env:NG_HEAD_TO_HEAD_SEED = "211"
$env:NG_HEAD_TO_HEAD_N_FOUNDERS = "24"
$env:NG_HEAD_TO_HEAD_N_CHR = "4"
$env:NG_HEAD_TO_HEAD_SEG_SITES = "80"
$env:NG_HEAD_TO_HEAD_SNP_PER_CHR = "20"
$env:NG_HEAD_TO_HEAD_QTL_PER_CHR = "5"
Rscript nextgen_cross_design/tools/run_head_to_head_benchmark.R
```

Output checks:

| Output | Rows | Notes |
| --- | ---: | --- |
| `*_summary.csv` | 108 | 3 scenarios x 3 parent sizes x 2 reps x 6 methods |
| `*_selections.csv` | 324 | 3 selected crosses per summary row |
| `*_scores.csv` | 1284 | all candidate pairs for 8, 12, and 16 parents across scenarios/reps |
| `*_comparisons.csv` | 1188 | candidate rows only, compared against each style-proxy baseline across available metrics |
| `*_winner_summary.csv` | 54 | 3 scenarios x 3 parent sizes x 6 realized metrics |
| `*_config.csv` | 18 | one row per scenario/parent-size/rep case |
| `*_method_registry.csv` | 6 | candidate/baseline labels plus exact-external caveats |

Primary winner counts in `winner_summary`:

| Method | Winner rows |
| --- | ---: |
| `alphamate_style_weighted_ocs` | 30 |
| `nextgen_desired_gain_ocs` | 10 |
| `popvar_style_weighted_topn` | 9 |
| `nextgen_economic_index_ocs` | 5 |

Tie-size distribution:

| Tied method count | Winner rows |
| ---: | ---: |
| 1 | 10 |
| 2 | 29 |
| 3 | 1 |
| 4 | 13 |
| 6 | 1 |

Candidate-vs-style-baseline comparison checks:

- All 1188 comparison rows had `method_role = candidate`.
- Baseline metadata in `*_comparisons.csv` preserved
  `baseline_implementation = style_proxy` and
  `baseline_exact_external_status = not_run_multitrait_ci_proxy` for PopVar,
  SimpleMating, and AlphaMate rows.
- For `mean_realized_index`, `nextgen_auto_ocs` tied
  `alphamate_style_weighted_ocs` in all 18 scenario/parent-size/rep cells,
  beat the PopVar/SimpleMating-style top-N baselines in 5 cells, and tied those
  baselines in 8 cells.
- Across all available comparison metrics, `nextgen_desired_gain_ocs` and
  `nextgen_economic_index_ocs` had positive average deltas for realized yield,
  unique parents, max parent use, and group coancestry, but negative average
  deltas for realized disease, realized index, and realized quality.

Interpretation:

- The larger style-proxy grid is working and produces auditable comparison
  artifacts at a practical replicated scale.
- The results are mixed by scenario, parent size, and metric. The
  AlphaMate-style weighted OCS proxy is the most frequent primary winner, while
  NextGen desired-gain and economic-index OCS are strongest on yield-oriented
  and diversity-use metrics.
- This grid does not prove superiority over PopVar, SimpleMating, AlphaMate, or
  true polyploid AlphaMate workflows. It is evidence that the same-input
  benchmark process is now operational. Exact external package/executable runs
  remain a separate validation tier.

## Implementation: head-to-head visual decision report

Date: 2026-05-07.

The multi-trait head-to-head benchmark now has a breeder-facing HTML report
renderer:

```text
Rscript nextgen_cross_design/tools/render_head_to_head_visual_report.R
```

It reads `*_summary.csv`, `*_selections.csv`, `*_comparisons.csv`,
`*_winner_summary.csv`, and `*_method_registry.csv` for the configured
`NG_HEAD_TO_HEAD_PREFIX`, then writes `*_visual_report.html`. The report is a
frontend layer over the R backend, not a replacement for the benchmark CSVs. It
shows selection-index response, trait-direction response, gain-diversity
frontier, parent contribution, mate allocation matrix, ranked crossing
decisions, winner summaries, candidate-vs-baseline deltas, and explicit
style-proxy evidence caveats.

## Implementation: lab-server frontend boundary

Date: 2026-05-07.

The benchmark now exports a stable dashboard JSON artifact with
`schema_version = ng_head_to_head_dashboard.v1`:

```text
Rscript nextgen_cross_design/tools/export_head_to_head_dashboard_json.R
```

The new `frontend/` Next.js app reads these dashboard artifacts, lists runs,
opens an interactive breeder/QG decision workspace, exposes health/run/method
API routes, and includes a Postgres schema plus a separate R worker scaffold
for lab-server deployment. This keeps the R package as the scientific backend
while giving the UI a typed artifact contract instead of scraping generated
HTML.

## Implementation: `ng_poly4x_*`

Date: 2026-05-04.

The autotetraploid model is intentionally separate from the DH/RIL benchmark.
The implementation target is:

- `ng_poly4x_var_topn`
- `ng_poly4x_usefulness_topn`
- `ng_poly4x_ocs`

The runner validates potato/cassava-like 4x AlphaSimR scenarios by simulating
4x progeny directly with `makeCross()`, pulling marker dosage in `0..4`, scoring
candidate crosses from sampled progeny, and validating selected crosses against
an independent realized progeny sample.

## Run: `poly4x_grid`

Date: 2026-05-06.

This is the first default true-autotetraploid grid after making
`run_poly4x_grid.ps1` include policy modes by default.

Configuration:

- Scenarios: `potato_autotetraploid_4x`, `cassava_autotetraploid_4x`.
- Parent sizes: 20 and 40.
- One replicate, one cycle.
- Default policy modes: `gain`, `diversity`, and `ocs`, alongside legacy
  `ng_poly4x_var_topn`, `ng_poly4x_usefulness_topn`, and `ng_poly4x_ocs`.
- Small direct-4x AlphaSimR validation: sampled scoring families and
  independent realized progeny samples.

Realized selected-family top10 winners:

| Scenario | Parents | Best top10 method | Mean GV | Top10 GV | Max GV |
| --- | ---: | --- | ---: | ---: | ---: |
| `cassava_autotetraploid_4x` | 20 | `poly4x_policy_gain` | 1.658 | 2.681 | 3.192 |
| `cassava_autotetraploid_4x` | 40 | `ng_poly4x_ocs` | 2.153 | 3.216 | 3.841 |
| `potato_autotetraploid_4x` | 20 | `poly4x_policy_diversity` | 1.521 | 2.934 | 3.559 |
| `potato_autotetraploid_4x` | 40 | `poly4x_policy_diversity` | 1.874 | 2.938 | 3.720 |

Interpretation:

This is smoke evidence that the true-4x grid now evaluates the policy modes
end to end. It is not a production recommendation or a superiority claim. The
results are mixed by scenario and parent size, with policy modes winning three
of four top10 cells and legacy `ng_poly4x_ocs` winning the cassava 40-parent
cell. The next useful 4x run should increase replicates and cycles before
promoting a default beyond the current policy-mode labels.

## Run: `poly4x_grid_2rep2cycle_20260506`

Date: 2026-05-06.

This run repeats the default true-4x grid with more evidence than the initial
smoke screen.

Configuration:

- Scenarios: `potato_autotetraploid_4x`, `cassava_autotetraploid_4x`.
- Parent sizes: 20 and 40.
- Two replicates, two cycles.
- Methods: `ng_poly4x_var_topn`, `ng_poly4x_usefulness_topn`,
  `ng_poly4x_ocs`, `poly4x_policy_gain`, `poly4x_policy_diversity`, and
  `poly4x_policy_ocs`.
- Small direct-4x AlphaSimR validation with sampled scoring families and
  independent realized progeny samples.

Average selected-family winners across replicates and cycles:

| Scenario | Parents | Best top10 method | Top10 GV | Best mean method | Mean GV | Best max method | Max GV |
| --- | ---: | --- | ---: | --- | ---: | --- | ---: |
| `cassava_autotetraploid_4x` | 20 | `ng_poly4x_usefulness_topn` | 3.051 | `poly4x_policy_ocs` | 1.991 | `ng_poly4x_usefulness_topn` | 3.862 |
| `cassava_autotetraploid_4x` | 40 | `ng_poly4x_ocs` | 3.677 | `ng_poly4x_ocs` | 2.482 | `ng_poly4x_ocs` | 4.485 |
| `potato_autotetraploid_4x` | 20 | `poly4x_policy_ocs` | 3.239 | `poly4x_policy_ocs` | 2.179 | `poly4x_policy_ocs` | 3.833 |
| `potato_autotetraploid_4x` | 40 | `poly4x_policy_gain` | 3.380 | `poly4x_policy_gain` | 2.158 | `poly4x_policy_gain` | 4.170 |

Final-cycle top10 winners were the same method families: `ng_poly4x_usefulness_topn`
for cassava p20, `ng_poly4x_ocs` for cassava p40, `poly4x_policy_ocs` for
potato p20, and `poly4x_policy_gain` for potato p40.

Interpretation:

This run strengthens the conclusion that true-4x policy modes are functional,
but it does not support a single universal 4x default. Policy modes were best
for potato in this grid, while legacy usefulness/OCS methods were best for
cassava. The next 4x decision layer should therefore be scenario-aware rather
than promoting one global mode. A stronger validation should use more
replicates, more cycles, larger parent sizes, and possibly real AlphaMate
comparison if the 4x scoring-to-AlphaMate interface is enabled for the same
parent population.

Implementation follow-up:

- `ng_poly4x_policy_select(mode = "auto")` now uses this grid as an
  evidence-scoped scenario-aware selector.
- Potato-like scenarios route to `poly4x_policy_ocs` at 20 parents and
  `poly4x_policy_gain` above 20 parents.
- Cassava-like scenarios route to legacy `ng_poly4x_usefulness_topn` at
  20 parents and legacy `ng_poly4x_ocs` above 20 parents.
- Explicit `mode = "gain"`, `"diversity"`, or `"ocs"` remains unchanged.
- Selector diagnostics record `source = "poly4x_grid_2rep2cycle_20260506"`,
  scenario-specific `reason`, primary method/family, fallback status, crop,
  and parent-count band.

## Implementation: External AlphaMate Baseline

Date: 2026-05-03.

The official AlphaGenes AlphaMate repository is cloned under
`external/AlphaMate`, and `ng_select_alphamate()` can now run
`external/AlphaMate/binaries/AlphaMate.exe` as a real external baseline. The
wrapper writes AlphaMate `Criterion`, relationship matrix, and spec files,
executes the binary in its working directory, parses `ModeOptTarget1`, and
maps the selected mating plan back to the benchmark score table. Selection
summaries record `alphamate_mode`, `alphamate_target_degree`,
`alphamate_criterion_col`, `alphamate_exit_code`, executable path, runtime
path, and work directory.

On this machine the AlphaMate binary requires `libiomp5md.dll`; the successful
runs used `NG_ALPHAMATE_RUNTIME_PATH=C:\Python\Lib\site-packages\torch\lib`.
`NG_ALPHAMATE_EXE` can point to another binary. The active R library reported
`PopVar=FALSE` and `SimpleMating=FALSE` during this pass, so this section is
about real AlphaMate only, not a refreshed external PopVar/SimpleMating run.

Implemented benchmark method names:

- `alphamate_opt` and numeric target variants such as `alphamate_opt30`,
  `alphamate_opt45`, and `alphamate_opt60`;
- `alphamate_maxcriterion` and `alphamate_mincoancestry` for diagnostics, but
  these should not be promoted unless their outputs satisfy the benchmark's
  non-self, non-repeated cross constraints.

Smoke checks:

- `Rscript nextgen_cross_design/tests/alphamate_external.R` first failed on a
  missing plan-matching diagnostic, then passed after adding explicit
  self-cross, repeated-pair, and unmatched-pair validation.
- `alphamate_wrapper_smoke`: 20 parents, one replicate, one cycle, two
  chromosomes, 80 SNP/chr, 10 QTL/chr. `alphamate_opt45` ran with
  `alphamate_exit_code = 0` and `ModeOptTarget1`.
- `ModeMaxCriterion` was tested separately and produced self-cross/repeated
  output despite the spec disabling selfing and repeated matings; the wrapper
  now reports those invalid pairs directly instead of a generic mismatch.

Compact real-AlphaMate parent-size screen:

- Run prefix: `alphamate_external_grid_smoke`.
- Parent sizes: 20 and 40.
- One replicate, one breeding cycle, reduced two-chromosome genome, pure-R
  variance kernel.
- Methods: `var_simple_topn`, `alphamate_opt30`, `alphamate_opt45`,
  `alphamate_opt60`, and `ng_frontier_policy_ocs10_lps2`.
- All six AlphaMate runs exited with code 0.

Cycle-1 top10 results:

| Parents | Best Top10 Method | Best Top10 | `var_simple_topn` | Best AlphaMate | Frontier | Note |
| ---: | --- | ---: | ---: | ---: | ---: | --- |
| 20 | `alphamate_opt45` tied with `var_simple_topn` | 3.415 | 3.415 | 3.415 (`opt45`) | 3.411 | Single-rep tie on top10; `var_simple_topn` had higher mean and transgressive rate. |
| 40 | `ng_frontier_policy_ocs10_lps2` | 3.727 | 3.528 | 3.564 (`opt30`) | 3.727 | Frontier beat the best AlphaMate target by +0.163 top10 in this smoke screen. |

This confirms the real AlphaMate integration and gives an initial benchmark
signal, but it is not sufficient to claim general superiority over AlphaMate.
The next publishable comparison should run replicated parent-size grids with
real AlphaMate target-degree sweeps under the same 5K-marker protocol used for
the frontier-policy validation.

Full 5K real-AlphaMate parent-size validation:

- Run prefix: `alphamate_external_grid_5k`.
- Parent sizes: 20, 30, 40, 50, 60, 70, 80.
- 3 AlphaSimR replicates, 3 breeding cycles, 5K SNPs, 400 effect-training
  individuals, single-threaded AlphaSimR, C++ recombination kernel.
- Methods: `var_simple_topn`, `alphamate_opt30`, `alphamate_opt45`,
  `alphamate_opt60`, and `ng_frontier_policy_ocs10_lps2`.
- `tools/run_alphamate_external_grid_5k.ps1` captures the standard run
  settings. The 80-parent scenario required a longer timeout but used the same
  1000-iteration AlphaMate settings as the other parent sizes.
- All 189 AlphaMate selections exited with code 0: 63 runs each for
  `alphamate_opt30`, `alphamate_opt45`, and `alphamate_opt60`. All parsed
  `ModeOptTarget1` with criterion `cross_mean_adjusted_pheno`.

Average cycle top10 across cycles 1 to 3:

| Parents | Frontier Top10 | Var TopN | Best AlphaMate | Best AlphaMate Top10 | Frontier - Var | Frontier - Best AlphaMate |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 20 | 5.850 | 5.130 | `alphamate_opt30` | 5.123 | +0.719 | +0.726 |
| 30 | 5.992 | 5.362 | `alphamate_opt30` | 5.530 | +0.630 | +0.462 |
| 40 | 6.060 | 5.696 | `alphamate_opt30` | 5.666 | +0.364 | +0.393 |
| 50 | 5.961 | 5.610 | `alphamate_opt45` | 5.718 | +0.351 | +0.243 |
| 60 | 6.278 | 6.029 | `alphamate_opt30` | 5.952 | +0.249 | +0.326 |
| 70 | 6.461 | 6.353 | `alphamate_opt30` | 5.989 | +0.109 | +0.473 |
| 80 | 6.548 | 6.327 | `alphamate_opt30` | 6.266 | +0.220 | +0.282 |

The frontier policy won all seven parent sizes for `top10_gv`, and also won
`mean_gv`, `max_gv`, and transgressive rate in the combined winner summary.
Average top10 delta was +0.377 versus `var_simple_topn` and +0.415 versus the
best real AlphaMate target degree at each parent size.

This is now direct evidence against the official AlphaMate executable for the
validated DH/RIL AlphaSimR setup. It is still evidence-scoped, not universal:
new crops, ploidy systems, training designs, trait architectures, and
recombination landscapes should generate a new validation grid before treating
the frontier policy as a production default.

## Implementation: `ng_frontier_policy_ocs*`

Date: 2026-05-03.

The empirical frontier is now implemented as `ng_frontier_policy_ocs*`. This is
a policy layer, not a new biological score. It chooses a delegated method from
the replicated 5K-marker parent-size frontier and then runs that method through
the normal selection path.

Default policy source: `validated_5k_parent_grid_2026_05_03`.

| Parent band | Primary method | Fallback candidates |
| ---: | --- | --- |
| 1-25 | `ng_recomb_gebv_ocs10_lps2` | guarded router |
| 26-35 | `ng_meta_selector_ocs10_lps2` | meta-portfolio, guarded router |
| 36-45 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | guarded router |
| 46-55 | `ng_meta_router_ocs10_lps2` | meta-portfolio |
| 56-65 | `popvar_uc_ocs10_lps1` | SimpleMating `usefa`, guarded router |
| 66-75 | `ng_meta_selector_ocs10_lps2` | meta-portfolio, guarded router |
| 76+ | `ng_meta_portfolio_ocs10_lps2` | guarded router |

The policy records `frontier_policy_*` columns in the selection summary and
`pred_ng_frontier_policy_*` columns in family outputs. `NG_FRONTIER_POLICY_SPEC`
can replace the built-in bands with a crop/site-specific validation rule, and
`NG_FRONTIER_POLICY_FALLBACK_METHOD` controls the last-resort method.

This is deliberately not described as agnostic to every dataset or crop. The
current policy is evidence-scoped to the validated DH/RIL AlphaSimR setup. It
must be revalidated when parent counts, trait architecture, marker density,
training population size, recombination landscape, selection target, or external
PopVar/SimpleMating availability changes.

Smoke checks:

- `Rscript nextgen_cross_design/tests/frontier_policy.R` passed after first
  failing on the missing API and then on missing external-baseline detection.
- `frontier_policy_smoke`: 20 parents, 2 chromosomes, 80 SNP/chr, 1 cycle, 1
  rep. `ng_frontier_policy_ocs10_lps2` delegated to
  `ng_recomb_gebv_ocs10_lps2` and wrote the policy diagnostics.
- `frontier_policy_wrapper_smoke`: the dedicated validation-grid wrapper ran a
  reduced 20-parent grid and wrote frontier-policy diagnostics.

Full parent-size validation:

- Run prefix: `frontier_policy_parent_grid_5k`.
- Parent sizes: 20, 30, 40, 50, 60, 70, 80.
- 3 AlphaSimR replicates, 3 breeding cycles, 5K SNPs, 400 effect-training
  individuals, single-threaded AlphaSimR, C++ recombination kernel.
- Methods included `var_simple_topn`, `var_simple_ocs10_lps2`,
  PopVar usefulness, SimpleMating `usefa`, recombination-GEBV OCS,
  PMV-scaled balanced OCS, meta portfolio, meta selector, guarded meta router,
  and `ng_frontier_policy_ocs10_lps2`.

Average top10 across cycles 1 to 3:

| Parents | Frontier Top10 | Var TopN | PopVar | SimpleMating | Router | Frontier - Var | Frontier - PopVar | Delegated Method |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 20 | 5.957 | 5.130 | 5.906 | 5.906 | 5.876 | +0.826 | +0.050 | `ng_recomb_gebv_ocs10_lps2` |
| 30 | 6.091 | 5.362 | 5.986 | 5.986 | 5.568 | +0.729 | +0.104 | `ng_meta_selector_ocs10_lps2` |
| 40 | 6.030 | 5.696 | 5.644 | 5.644 | 5.963 | +0.334 | +0.386 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` |
| 50 | 6.012 | 5.610 | 5.993 | 5.993 | 6.012 | +0.402 | +0.019 | `ng_meta_router_ocs10_lps2` |
| 60 | 6.278 | 6.029 | 6.278 | 6.278 | 6.238 | +0.249 | +0.000 | `popvar_uc_ocs10_lps1` |
| 70 | 6.500 | 6.353 | 6.310 | 6.310 | 6.495 | +0.147 | +0.190 | `ng_meta_selector_ocs10_lps2` |
| 80 | 6.603 | 6.327 | 6.495 | 6.495 | 6.603 | +0.275 | +0.107 | `ng_meta_portfolio_ocs10_lps2` |

The frontier policy won or tied the best observed top10 method at all seven
parent sizes. Average top10 deltas were +0.423 versus `var_simple_topn`, +0.122
versus PopVar/SimpleMating usefulness, and +0.102 versus the guarded router.
Fallback rate was zero in all parent-size bands.

This is the strongest internal evidence so far. It supports using
`ng_frontier_policy_ocs10_lps2` as the current validated default for the
deterministic DH/RIL AlphaSimR setup. It is not a claim that the framework has
beaten the real AlphaMate software, because this grid compared against
AlphaMate-style gain-diversity allocation implemented in this repo, not an
external AlphaMate run.

## Implementation: Crop-Genome Scenario Screen

Date: 2026-05-03.

Added a crop-genome portability harness:

- `NG_GENOME_LENGTH_M` is now consumed by the AlphaSimR benchmark and passed to
  `quickHaplo(..., genLen=...)`.
- `ng_crop_genome_scenarios()` defines diploid DH/RIL-compatible templates
  plus explicit diploidized stress approximations for crops whose production
  biology is polyploid or clonal.
- `tools/run_crop_genome_scenarios.R` runs `run_parent_size_grid.ps1` across
  those templates and combines `overall_avg`, `comparison_avg`,
  `winner_summary`, and `selection_avg` outputs by scenario.
- `tools/run_crop_genome_alphamate_grid.ps1` runs the crop scenario harness
  with real AlphaMate target-degree controls.

The templates vary chromosome count, genetic map length, marker density, QTL
density, heritability, founder count, and effect-training size. The current
AlphaSimR benchmark path still uses a diploid DH/RIL crossing engine, so
`bread_wheat_hexaploid_approx`, `potato_tetraploid_stress`,
`cassava_tetraploid_stress`, and `sugarcane_polyploid_stress` are scale and
diversity stress tests, not true polyploid dosage-inheritance models.

Smoke checks:

- `Rscript nextgen_cross_design/tests/crop_genome_scenarios.R` first failed on
  the missing scenario API, then failed on a named-environment propagation bug,
  and now passes.
- `crop_genome_smoke`: compact-selfing template, 20 parents, one replicate,
  one cycle, two methods. The run completed and wrote combined scenario output.
- `crop_genome_three_scenario_smoke`: compact-selfing, maize-like, and large
  diploid approximation at 20 parents, one replicate, one cycle, two methods.
  The frontier policy delegated to `ng_recomb_gebv_ocs10_lps2` in all three
  cases with no fallback. This confirms harness wiring only; it is not enough
  to claim cross-crop superiority.

Three-scenario smoke top10 results versus `var_simple_topn`:

| Scenario | n_chr | Map length M/chr | SNP/chr | QTL/chr | h2 | Frontier top10 | Delta top10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `compact_selfing` | 7 | 0.75 | 650 | 24 | 0.55 | 4.476 | +0.006 |
| `maize_like` | 10 | 1.60 | 700 | 36 | 0.45 | 4.784 | +0.696 |
| `large_diploid_approx` | 21 | 1.20 | 350 | 22 | 0.35 | 4.414 | +0.124 |

Expanded real-AlphaMate crop smoke:

- Run prefix: `crop_genome_alphamate_smoke`.
- Scenarios: compact selfing, maize, bread wheat approximation, barley, field
  pea, potato tetraploid stress, cassava diploid, cassava tetraploid stress,
  and sugarcane polyploid stress.
- Parent size: 20.
- One replicate, one cycle, C++ kernel, real AlphaMate `opt30`, `opt45`, and
  `opt60` with 200 evolutionary iterations.
- All 27 AlphaMate selections exited with code 0 and parsed `ModeOptTarget1`.

Smoke top10 versus best real AlphaMate target:

| Scenario | Crop | Harness model | Frontier top10 | Best AlphaMate | Delta |
| --- | --- | --- | ---: | ---: | ---: |
| `barley_like` | barley | diploid DH/RIL | 4.951 | 4.523 (`opt60`) | +0.428 |
| `bread_wheat_hexaploid_approx` | wheat | diploidized approximation | 5.563 | 5.249 (`opt30`) | +0.315 |
| `cassava_diploid` | cassava | diploid DH/RIL | 4.524 | 3.916 (`opt45`) | +0.608 |
| `cassava_tetraploid_stress` | cassava | diploidized approximation | 3.764 | 4.215 (`opt60`) | -0.451 |
| `compact_selfing` | generic selfing | diploid DH/RIL | 4.476 | 4.422 (`opt30`) | +0.054 |
| `field_pea_like` | field pea | diploid DH/RIL | 6.445 | 6.445 (`opt45`) | +0.000 |
| `maize_like` | maize | diploid DH/RIL | 4.784 | 4.189 (`opt30`) | +0.595 |
| `potato_tetraploid_stress` | potato | diploidized approximation | 4.693 | 4.410 (`opt45`) | +0.283 |
| `sugarcane_polyploid_stress` | sugarcane | diploidized approximation | 3.478 | 3.529 (`opt45`) | -0.050 |

This smoke screen strengthens the diploid-crop signal but shows two important
stress failures: cassava tetraploid and sugarcane polyploid approximations.
The next scientific step is therefore not to retune the frontier bands on the
same DH/RIL assumptions. It is to add a true polyploid/clonal simulation path
or, at minimum, run a replicated crop grid and produce a separate
`NG_FRONTIER_POLICY_SPEC` for the diploidized stress scenarios.

## Implementation: Crop-Aware Policy Layer

Date: 2026-05-04.

Added `ng_crop_aware_policy_ocs*` as a separate dispatch family rather than
modifying `ng_frontier_policy_ocs*`. The crop policy reads `NG_CROP_SCENARIO`,
`NG_CROP`, `NG_CROP_HARNESS_MODEL`, and `NG_CROP_VALIDATION_SCOPE`, then records
`crop_policy_*` and `pred_ng_crop_policy_*` diagnostics in the benchmark output.
`ng_crop_genome_env()` now emits those metadata values for every crop-grid
subprocess.

Current smoke-derived dispatch:

- ordinary diploid DH/RIL templates route first to `ng_frontier_policy_ocs10_lps2`;
- `cassava_tetraploid_stress` routes first to `alphamate_opt60`;
- `sugarcane_polyploid_stress` routes first to `alphamate_opt45`;
- all crop-aware candidates retain fallbacks and error logs, so failed
  candidate methods are visible instead of silently changing the plan.

This is intentionally a policy layer over the current diploidized stress
harness. It is not a true polyploid dosage-inheritance engine and should not be
reported as crop-agnostic proof.

Verification:

- `Rscript nextgen_cross_design/tests/crop_aware_policy.R` passes.
- `Rscript nextgen_cross_design/tests/alphamate_external.R` now includes a
  long-parent-ID regression. `ng_select_alphamate()` aliases parent IDs before
  writing AlphaMate files and restores the original IDs after matching, because
  the external executable truncates long IDs such as
  `ng_crop_aware_policy_ocs10_lps2_R1_C0_P1`.
- `crop_aware_policy_smoke`: maize, cassava tetraploid stress, and sugarcane
  polyploid stress at 20 parents, one replicate, one cycle, real AlphaMate
  `opt30`, `opt45`, and `opt60`. The crop-aware policy routed maize to
  `ng_frontier_policy_ocs10_lps2`, cassava tetraploid stress to
  `alphamate_opt60`, and sugarcane polyploid stress to `alphamate_opt45`, all
  with `crop_policy_fallback = 0`.

## Implementation: `ng_meta_router_ocs*`

The router is now implemented as a separate method, not as a replacement for the previous meta selector. It builds complete candidate mating plans, then chooses one plan family using only pre-selection information:

- prior-cycle method-family performance, standardized within cycle to avoid breeding-value drift;
- current predicted `ng_meta_score` rank of the selected plan;
- current candidate-family gain rank;
- plan diversity from parent use, group coancestry, and pair kinship;
- SNP-effect reliability, which downweights PMV/recombination families when marker effects are weak.

Default candidate families are `recomb_gebv`, `pmv_balanced`, `portfolio`, `popvar_uc`, `simple_usefa`, and `var_simple`. PopVar/SimpleMating columns are now computed automatically when the router is requested, if those packages are installed. The router records its choice in the selection and family CSV outputs through `meta_router_*` and `pred_ng_meta_router_*` columns. Prior routed choices are classified by `pred_ng_meta_router_family`, so the router can learn from its own previous cycles as well as from standalone control methods.

The router now has a regret guard enabled by default. The guard keeps the
existing router score as the primary decision, then compares the chosen
candidate against the other candidate plans using a weighted branch-local
plan/meta score and native gain score. If another plan is better by at least
`NG_META_ROUTER_REGRET_MIN_ADVANTAGE` and is not too far behind on the router
score (`NG_META_ROUTER_REGRET_MAX_ROUTER_PENALTY`), the router uses that plan
instead. This is meant to prevent large branch-local misses without hard-coding
the 20/30/40/50/60/70/80 parent-size frontier. Guard diagnostics are recorded
as `meta_router_guard_*` and `pred_ng_meta_router_guard_*` columns.

Smoke checks:

- `Rscript nextgen_cross_design/tests/smoke_test.R` passed.
- `Rscript nextgen_cross_design/tests/meta_router_regret_guard.R` passed.
- `regret_guard_smoke`: 20 parents, 2 chromosomes, 40 SNP/chr, 1 cycle, 1
  rep. The benchmark harness completed and wrote `meta_router_guard_*`
  diagnostics.
- `meta_router_smoke`: 20 parents, 100 SNP/chr, 2 cycles, 1 rep. Router diagnostics were written; with default `NG_META_ROUTER_MIN_HISTORY_N=10`, cycle-2 history remained neutral because each method had only 6 selected families.
- `meta_router_history_smoke`: same setup with `NG_META_ROUTER_MIN_HISTORY_N=3`. Cycle-2 router consumed prior-cycle family history (`meta_router_history_score` nonzero) and selected the portfolio candidate.

This is a functional router, not yet a validated state-of-art win. The next evidence step is a deterministic 5K-marker parent-size grid comparing `ng_meta_router_ocs10_lps2` against `var_simple`, PopVar/SimpleMating usefulness, PMV-balanced, recombination-GEBV, `ng_meta_portfolio`, and `ng_meta_selector`.

First 5K grid result: the router was competitive at 60/80 parents but missed the winner at 20/40 parents because the history family for `var_simple` also included the non-allocator `var_simple_topn`, and the score leaned too strongly on the portfolio-derived current `ng_meta_score`. Router v2 therefore uses allocator-aware history patterns, caps the minimum history requirement by the selected cross count, and shifts default weights toward prior family performance, plan diversity, and native family gain. Router v3 adds a smooth parent-size/reliability prior: small-parent blocks favor robust `var_simple`, middle parent counts favor PMV-balanced, and larger parent counts favor portfolio/recombination. Current router weights are `history=0.50`, `prior=0.25`, `plan=0.05`, `gain=0.10`, `diversity=0.12`, `reliability=0.03`. The same allocator-aware method-history mapping is now used inside the meta portfolio score.

## Run: `diagnostic_framework_validation_5k_3rep_20_80`

Date: 2026-05-03.

The full diagnostic-first validation wrapper was run outside the sandbox with
C++ enabled:

- `NG_VALIDATION_PHASE=all`
- `NG_VALIDATION_USE_CPP=1`
- `NG_VALIDATION_REPS=3`
- `NG_VALIDATION_PARENT_SIZES=20,30,40,50,60,70,80`
- `NG_VALIDATION_EFFECT_TRAINING_N=400`

Outputs are summarized in
`results/diagnostic_parent_grid_5k_validation_report.md`.

Family calibration completed for all 21 parent-size/replicate combinations:
PopVar, SimpleMating, and the genomicMateSelectR-derived variance check all
reported `ok`.

Replicated parent-size top10 frontier:

| Parents | Top10 winner | Top10 | Router top10 | Router gap |
| ---: | --- | ---: | ---: | ---: |
| 20 | `ng_recomb_gebv_ocs10_lps2` | 5.970 | 5.871 | -0.098 |
| 30 | `ng_meta_selector_ocs10_lps2` | 6.091 | 5.606 | -0.485 |
| 40 | `ng_meta_router_ocs10_lps2` | 6.022 | 6.022 | 0.000 |
| 50 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 6.084 | 6.035 | -0.049 |
| 60 | `popvar_uc_ocs10_lps1` | 6.278 | 5.990 | -0.287 |
| 70 | `ng_meta_selector_ocs10_lps2` | 6.500 | 6.495 | -0.005 |
| 80 | `ng_meta_portfolio_ocs10_lps2` | 6.603 | 6.596 | -0.006 |

Allocator crosscheck top10 winners were also mixed: `ng_ocs_mip10_lps1` at 20,
`popvar_musp_ocs10_lps1` at 30, `ng_hybrid_select4` at 40,
`ng_hybrid_cal_topn` at 50, `var_simple_ocs10_lps1` at 60 and 70, and
`popvar_musp_topn` at 80.

Interpretation:

This validates the framework diagnostic protocol, not a universal router win.
The router is competitive at 40, 70, and 80 parents, but the replicated grid
shows material misses at 30 and 60 parents. The next method change should not
be a local patch to the current router weights. It should formalize a
parent-size and live-evidence decision rule, or add a regret guard that prevents
the router from choosing a family whose branch-local plan score is materially
below the strongest recombination, PMV-balanced, PopVar/SimpleMating, or
portfolio candidate. Until that beats the frontier above in replicated
validation, the project should report the frontier instead of claiming a single
default winner.

## Run: `diagnostic_parent_grid_guard_5k_3rep_20_80`

Date: 2026-05-03.

This run repeated the replicated parent-size grid after adding the router
regret guard. It used the same core validation settings as
`diagnostic_framework_validation_5k_3rep_20_80`, with C++ enabled, three
replicates, three cycles, parent sizes 20/30/40/50/60/70/80, 5K SNPs, and
400 effect-training individuals. The outputs use prefix
`results/diagnostic_parent_grid_guard_5k_*`. The markdown report is
`results/diagnostic_parent_grid_guard_5k_validation_report.md`.

Replicated parent-size top10 frontier with the guard enabled:

| Parents | Top10 winner | Top10 | Guarded router top10 | Router gap |
| ---: | --- | ---: | ---: | ---: |
| 20 | `ng_recomb_gebv_ocs10_lps2` | 5.970 | 5.891 | -0.078 |
| 30 | `ng_meta_selector_ocs10_lps2` | 6.091 | 5.568 | -0.523 |
| 40 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.994 | 5.961 | -0.033 |
| 50 | `ng_meta_router_ocs10_lps2` | 6.020 | 6.020 | 0.000 |
| 60 | `popvar_uc_ocs10_lps1` | 6.278 | 6.238 | -0.040 |
| 70 | `ng_meta_selector_ocs10_lps2` | 6.500 | 6.495 | -0.005 |
| 80 | `ng_meta_portfolio_ocs10_lps2` | 6.603 | 6.603 | 0.000 |

Guarded router change versus the previous router grid:

| Parents | Old router top10 | Guarded router top10 | Delta |
| ---: | ---: | ---: | ---: |
| 20 | 5.871 | 5.891 | +0.020 |
| 30 | 5.606 | 5.568 | -0.038 |
| 40 | 6.022 | 5.961 | -0.062 |
| 50 | 6.035 | 6.020 | -0.015 |
| 60 | 5.990 | 6.238 | +0.248 |
| 70 | 6.495 | 6.495 | +0.000 |
| 80 | 6.596 | 6.603 | +0.006 |

Interpretation:

The regret guard is useful but incomplete. It corrected most of the severe
60-parent router miss and made the router top10 winner at 50 parents, while
keeping the router essentially tied with the frontier at 70 and 80 parents.
However, it worsened the 30-parent router result and no longer wins the
40-parent top10 comparison, where PMV-balanced is better. This is evidence
against more blind threshold tuning. The most practical product behavior is now
an empirical frontier/default-policy layer: use recombination-GEBV at 20,
meta-selector at 30, PMV-balanced at 40, guarded router at 50,
PopVar/SimpleMating-style usefulness at 60, meta-selector at 70, and
meta-portfolio/guarded-router at 80, while retaining all controls in validation.

## Runs: `meta_router_v4_prior_5k_1rep_3cyc_20_40` and `meta_router_v4_prior_5k_1rep_3cyc_60_80`

Date: 2026-04-30 and 2026-05-03.

This continued the deterministic 5K router validation across 20, 40, 60, and
80 parents with one AlphaSimR replicate, three cycles, 400 effect-training
individuals, 3 phenotype reps, and shortlist/exact external rescoring with
multiplier 20. The 20/40 results were already present; the 60/80 continuation
was run from the handoff state.

The 60/80 continuation emitted an Rtools/make warning inside the sandbox and
therefore used the R fallback despite `NG_USE_CPP=1`. The dedicated C++/R
kernel consistency test passed when run outside the sandbox, so the ranking
results are usable, but timing from the sandboxed 60/80 continuation should
not be used as C++ performance evidence.

Average cycle 1-3 metrics:

| Parents | Best Mean | Best Top10 | Router Family Choices | Router Mean | Router Top10 |
| ---: | --- | --- | --- | ---: | ---: |
| 20 | `ng_meta_router_ocs10_lps2` / `var_simple_ocs10_lps2` 4.557 | `ng_meta_router_ocs10_lps2` / `var_simple_ocs10_lps2` 5.165 | `var_simple`, `var_simple`, `var_simple` | 4.557 | 5.165 |
| 40 | `var_simple_ocs10_lps2` 4.641 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` 5.372 | `pmv_balanced`, `var_simple`, `pmv_balanced` | 4.607 | 5.266 |
| 60 | `ng_recomb_gebv_ocs10_lps2` 5.201 | `ng_meta_portfolio_ocs10_lps2` 5.820 | `pmv_balanced`, `portfolio`, `portfolio` | 4.836 | 5.437 |
| 80 | `ng_meta_portfolio_ocs10_lps2` / `ng_meta_router_ocs10_lps2` 5.465 | `var_simple_topn` 6.113 | `portfolio`, `portfolio`, `portfolio` | 5.465 | 6.101 |

Interpretation:

- Router v4 fixed the small-parent failure by routing 20-parent plans to the
  allocator-aware `var_simple` family.
- At 40 parents, the router was competitive but still below PMV-balanced on
  top10 and slightly below `var_simple_ocs10_lps2` on mean.
- At 60 parents, the router did not recover the best recombination/portfolio
  behavior. Its chosen family is applied to the router branch's current parent
  population, so it should not be expected to match the standalone portfolio
  branch after prior cycles diverge.
- At 80 parents, the router tracked the portfolio branch and nearly matched
  the `var_simple_topn` top10 winner while improving mean.

This is still not enough evidence to promote `ng_meta_router_ocs10_lps2` as a
default general winner. The next validation should be a replicated,
outside-sandbox C++ run over 20/30/40/50/60/70/80 parents. The next method
change should focus on the 60-parent miss: either make the router compare
branch-local leader/portfolio uncertainty more conservatively, or add a regret
guard that avoids routing to a family whose live plan score is materially below
the best available recombination or portfolio candidate.

## Run: `nextgen_80p_5k_3rep5cycle`

This run used the first implementation with in-sample marker-effect reliability. It was not statistically valid because in-sample reliability was treated as if it were predictive reliability.

Key result:

- `var_simple_topn` still won most cycles.
- The apparent marker-effect reliability was inflated, often above 0.5 and near 1.0 in cycle 1.

Conclusion:

The framework must not gate GEBV/PMV by in-sample R2.

## Run: `nextgen_80p_5k_3rep5cycle_cvrel`

This run used cross-validated marker-effect reliability and made the `var_simple` baseline use adjusted phenotype mean.

Key result:

- Cross-validated reliability was very low, about 0.02 to 0.11.
- In-sample reliability remained high, about 0.50 to 0.99.
- `var_simple_topn` still won the early and mid cycles.

Conclusion:

With 80 individuals and 5K markers, SNP-effect estimates are too weak to drive PMV directly. The earlier concern from the user is confirmed: when SNP effects cannot be estimated precisely, adjusted phenotype/BLUE/BLUP should dominate the mean component.

## Run: `nextgen_80p_5k_3rep5cycle_hybrid`

This run added:

- map-aware DH PMV scaled onto the `var_simple` variance scale;
- hybrid variance: 25% scaled PMV plus 75% `var_simple` when effect reliability is low;
- common phenotype-noise seeds across methods within replicate/cycle;
- AlphaMate-style allocation with max parent use of 4.

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_allocator` | +0.145 | +0.207 | +0.171 |
| `ng_uc_topn` | +0.115 | +0.177 | +0.175 |
| `ng_pmv_scaled_topn` | +0.054 | +0.079 | +0.071 |
| `ng_hybrid_topn` | -0.004 | +0.085 | +0.140 |
| `ng_allocator_hybrid` | -0.062 | +0.052 | +0.075 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_allocator` | -0.023 | -0.017 | -0.096 |
| `ng_uc_topn` | -0.054 | -0.087 | -0.162 |
| `ng_pmv_scaled_topn` | -0.100 | -0.150 | -0.285 |
| `ng_hybrid_topn` | -0.214 | -0.220 | -0.343 |
| `ng_allocator_hybrid` | -0.211 | -0.214 | -0.298 |

Interpretation:

The new framework can beat `var_simple_topn` at cycle 5 in this run, especially the `ng_allocator` and `ng_uc_topn` branches, but it does not yet dominate across all cycles. This is not yet a state-of-the-art win. It is evidence that the allocation layer and calibrated PMV signal can help later-cycle gain, but the early-cycle penalty is still too large.

The next required comparison is `var_simple_allocator`, which separates the allocation effect from the metric effect under the same parent-use cap.

## Run: `nextgen_80p_5k_3rep5cycle_hybrid_alloc`

This run added `var_simple_allocator` as a constrained allocation control.

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV | Max Parent Use |
| --- | ---: | ---: | ---: | ---: |
| `ng_allocator_hybrid` | +0.039 | +0.079 | +0.098 | 4 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 | 9 |
| `var_simple_allocator` | -0.102 | -0.007 | +0.056 | 4 |
| `ng_allocator` | -0.143 | -0.118 | -0.148 | 4 |
| `ng_uc_topn` | -0.166 | -0.151 | -0.106 | 10 |
| `ng_hybrid_topn` | -0.178 | -0.237 | -0.268 | 12 |
| `ng_pmv_scaled_topn` | -0.354 | -0.363 | -0.417 | 9 |

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_hybrid_topn` | -0.046 | -0.074 | -0.146 |
| `ng_uc_topn` | -0.108 | -0.126 | -0.124 |
| `ng_allocator_hybrid` | -0.136 | -0.130 | -0.176 |
| `ng_pmv_scaled_topn` | -0.192 | -0.257 | -0.265 |
| `ng_allocator` | -0.202 | -0.189 | -0.179 |
| `var_simple_allocator` | -0.240 | -0.250 | -0.205 |

Interpretation:

The allocation layer achieved the intended parent-use constraint. `ng_allocator_hybrid` was the only branch to beat unconstrained `var_simple_topn` at cycle 5 in this run, and it did so while capping max parent use at 4. However, no new method beat `var_simple_topn` averaged across cycles 1 to 5. This is not yet a robust win.

The most useful signal remains:

- raw `dh_pmv_var` has the best variance-calibration RMSE, but its scale is too small for direct usefulness scoring;
- scaled PMV alone over-weights variance and loses;
- the allocator can trade a small early-cycle penalty for a later-cycle gain and substantially better parent-use balance.

Next statistical work:

1. Calibrate PMV from historical realized family variance instead of median-scaling it to `var_simple`.
2. Add expected top-family value for actual family size, not only fixed usefulness intensity.
3. Evaluate across more replicates before claiming an advantage.
4. Compare `var_simple_allocator` and `ng_allocator_hybrid` over longer cycles only after PMV calibration is improved.

## Run: `nextgen_80p_5k_3rep5cycle_topk_cal`

This run added branch-specific variance calibration from prior realized families and expected top-k family value. With 40 progeny per cross and `NG_TOPK_PROP=0.10`, the score targets the expected mean of the top 4 DH progeny from each family.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `var_simple_etk_topn` | +0.014 | -0.019 | -0.009 |
| `ng_hybrid_cal_topn` | -0.010 | +0.043 | +0.064 |
| `ng_uc_topn` | -0.068 | -0.045 | -0.097 |
| `ng_pmv_cal_topn` | -0.148 | -0.228 | -0.229 |
| `ng_cal_allocator` | -0.307 | -0.211 | -0.224 |
| `var_simple_allocator` | -0.404 | -0.332 | -0.342 |

Interpretation:

Branch-specific calibration improved the ranking behavior for top-tail criteria but did not give a clear PMV mean-gain win. The selected families from one branch are too few and too biased to calibrate aggressively after only one or two cycles.

## Run: `nextgen_80p_5k_3rep5cycle_topk_globalcal`

This run used pooled historical calibration across all previous-cycle selected families in the replicate. It avoids same-cycle leakage by updating the calibration pool only after all branches finish a cycle.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_hybrid_cal_topn` | +0.082 | +0.090 | +0.076 |
| `var_simple_etk_topn` | +0.014 | -0.019 | -0.009 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_uc_topn` | -0.002 | -0.020 | +0.030 |
| `var_simple_allocator` | -0.096 | -0.024 | +0.002 |
| `ng_cal_allocator` | -0.137 | -0.113 | -0.087 |
| `ng_pmv_cal_topn` | -0.408 | -0.505 | -0.574 |

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV | Max Parent Use |
| --- | ---: | ---: | ---: | ---: |
| `ng_hybrid_cal_topn` | +0.208 | +0.177 | +0.207 | 11 |
| `var_simple_etk_topn` | +0.175 | +0.170 | +0.165 | 19 |
| `var_simple_allocator` | +0.096 | +0.203 | +0.378 | 4 |
| `ng_uc_topn` | +0.053 | +0.085 | +0.251 | 12 |
| `ng_cal_allocator` | +0.046 | +0.147 | +0.233 | 4 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 | 10 |
| `ng_pmv_cal_topn` | -0.662 | -0.732 | -0.717 | 13 |

Calibration:

- `dh_pmv_var` remained the best raw variance predictor by RMSE (`0.372`).
- Calibrated hybrid variance had better calibration RMSE than raw hybrid variance (`0.824` vs `2.033`).
- Pure calibrated PMV selected too much variance and poor means, so it lost badly.

Interpretation:

This is the first benchmark where a PMV-containing method beats `var_simple_topn` averaged across cycles. The winning metric is not pure PMV; it is historical-calibrated hybrid expected top-k value. The best constrained branch, `ng_cal_allocator`, protected parent use at 4 and beat `var_simple_topn` in the final cycle, but it still lost averaged across cycles.

Current conclusion:

The path forward is calibrated hybrid expected-top-k plus a better constrained allocator. Pure PMV should not be used alone in this low-training-size setting. The allocator needs to preserve the `ng_hybrid_cal_topn` gain while enforcing parent-use and diversity constraints.

## Run: `nextgen_80p_5k_3rep5cycle_repair_frontier`

This run tested a fast local repair allocator after unconstrained top-cross ranking. The repair variants enforced maximum parent-use caps of 4, 6, 8, or 10 and then used local upgrades to recover objective value.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_hybrid_cal_topn` | +0.027 | +0.008 | -0.019 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_cal_repair8` | -0.039 | -0.030 | -0.084 |
| `ng_cal_repair4` | -0.135 | -0.107 | -0.104 |
| `var_simple_repair4` | -0.138 | -0.077 | -0.106 |
| `ng_cal_repair6` | -0.241 | -0.236 | -0.306 |
| `ng_cal_repair10` | -0.312 | -0.351 | -0.474 |

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV | Max Parent Use |
| --- | ---: | ---: | ---: | ---: |
| `var_simple_repair4` | +0.091 | +0.005 | -0.079 | 4 |
| `ng_hybrid_cal_topn` | +0.006 | -0.153 | -0.277 | 13 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 | 12 |
| `ng_cal_repair8` | -0.079 | -0.116 | -0.200 | 8 |
| `ng_cal_repair4` | -0.229 | -0.234 | -0.171 | 4 |

Interpretation:

The repair allocator was fast and respected parent-use caps, but it did not recover enough objective value after removing over-used parents. This is a useful engineering baseline, but it is not the allocator we should rely on for final recommendations.

## Run: `nextgen_80p_5k_3rep5cycle_mip_frontier`

This run replaced the repair heuristic with an exact linear MIP allocator through `lpSolve`. The tested frontier used parent-use caps of 4, 6, 8, and 10 for the calibrated hybrid expected top-k score, plus a `var_simple_mip4` constrained control.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `var_simple_etk_topn` | +0.073 | +0.063 | +0.024 |
| `ng_cal_mip6` | +0.001 | -0.022 | -0.037 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_hybrid_cal_topn` | -0.028 | -0.086 | -0.093 |
| `ng_cal_mip10` | -0.029 | -0.051 | -0.111 |
| `ng_cal_mip8` | -0.117 | -0.151 | -0.238 |
| `var_simple_mip4` | -0.152 | -0.183 | -0.242 |
| `ng_cal_mip4` | -0.183 | -0.177 | -0.249 |

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV | Max Parent Use |
| --- | ---: | ---: | ---: | ---: |
| `var_simple_etk_topn` | +0.029 | -0.007 | -0.093 | 17 |
| `ng_hybrid_cal_topn` | +0.029 | -0.092 | -0.152 | 20 |
| `ng_cal_mip6` | +0.011 | -0.049 | -0.055 | 6 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 | 11 |
| `ng_cal_mip10` | -0.122 | -0.177 | -0.222 | 10 |
| `ng_cal_mip8` | -0.187 | -0.272 | -0.198 | 8 |
| `ng_cal_mip4` | -0.193 | -0.292 | -0.337 | 4 |
| `var_simple_mip4` | -0.249 | -0.283 | -0.284 | 4 |

Timing:

Exact MIP allocation took about `0.07` seconds per allocation stage in this benchmark. Optimization speed is no longer the limiting issue for the tested problem size.

Interpretation:

Exact allocation is working and fast, but hard parent-use caps remove too much short-term selection intensity. The best constrained PMV-containing branch was `ng_cal_mip6`: it was essentially neutral for mean gain across cycles, but still negative for top-tail and maximum progeny value. In this run the best average gain came from `var_simple_etk_topn`, which is still unconstrained and used parents heavily.

Updated conclusion:

The next improvement should be objective design, not another local search heuristic. We need an optimum-contribution style formulation that trades gain, calibrated family variance, parent contribution, and family size allocation in one objective. Hard caps can remain as guardrails, but they should not be the main diversity mechanism.

## Run: `nextgen_80p_5k_3rep5cycle_ocs_frontier`

This run added an optimum-contribution-style MIP. The objective uses calibrated hybrid expected top-k value, pair kinship, a soft parent-contribution penalty, and a loose parent-use guardrail. It also tested variable family-size allocation under the same total progeny budget.

Configuration:

- 80 parents, 5K markers, 3 phenotype reps, 5 breeding cycles, 3 simulation replicates.
- Global historical variance calibration.
- `NG_LAMBDA_PARENT_USE=10`, `NG_LAMBDA_GROUP=1`.
- OCS guardrail cap `10`.
- Family-size branches kept total progeny at `800`, with 20 selected crosses and `10` to `80` DH progeny per cross.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10` | +0.280 | +0.258 | +0.386 |
| `ng_hybrid_cal_topn` | +0.141 | +0.104 | +0.039 |
| `var_simple_etk_topn` | +0.122 | +0.063 | -0.070 |
| `ng_ocs_fam10` | +0.031 | -0.057 | -0.134 |
| `var_simple_ocs10` | +0.024 | -0.024 | -0.065 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_cal_mip6` | -0.049 | -0.063 | -0.086 |
| `var_simple_ocs_fam10` | -0.124 | -0.210 | -0.234 |

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_hybrid_cal_topn` | +0.105 | -0.016 | -0.157 |
| `ng_ocs_mip10` | +0.101 | +0.010 | -0.023 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `var_simple_etk_topn` | -0.043 | -0.149 | -0.276 |
| `var_simple_ocs10` | -0.106 | -0.158 | -0.143 |
| `ng_cal_mip6` | -0.115 | -0.089 | -0.208 |
| `ng_ocs_fam10` | -0.259 | -0.414 | -0.525 |
| `var_simple_ocs_fam10` | -0.342 | -0.364 | -0.384 |

Allocation summary:

| Method | Avg Unique Parents | Avg Max Parent Use | Avg Parent-Use Sq | Avg Group Coancestry |
| --- | ---: | ---: | ---: | ---: |
| `ng_ocs_fam10` | 16.1 | 7.8 | 0.106 | 0.201 |
| `ng_cal_mip6` | 13.3 | 6.0 | 0.108 | 0.178 |
| `ng_ocs_mip10` | 15.7 | 8.3 | 0.111 | 0.186 |
| `var_simple_topn` | 15.3 | 8.6 | 0.115 | 0.185 |
| `var_simple_etk_topn` | 14.9 | 10.1 | 0.133 | 0.225 |
| `ng_hybrid_cal_topn` | 14.6 | 10.1 | 0.134 | 0.245 |

Timing:

- C++ scoring was about `1.05` seconds per scoring stage.
- OCS MIP selection was about `0.59` to `0.68` seconds per allocation stage.
- Optimization is still not the bottleneck.

Interpretation:

This is the first clear win for the new design. `ng_ocs_mip10` beat `var_simple_topn`, `var_simple_etk_topn`, the unconstrained calibrated hybrid ranking, and the previous hard-cap PMV allocator on average across cycles. The gain is not coming from pure PMV alone; it comes from calibrated hybrid expected top-k plus a soft contribution-aware allocation objective.

The simple family-size allocator is not validated. It produced very balanced parent use and correctly preserved the total progeny budget, but it lost gain. The current family-size rule is too aggressive because it allocates `10` to `80` progeny from noisy cross scores without modeling marginal expected value. Family sizes should be revised as a marginal-value optimization, not a simple score-proportional allocation.

Conclusion at this stage:

Use `ng_ocs_mip10` as the next method to investigate in the 80-parent track. Keep `var_simple_topn`, `var_simple_etk_topn`, `ng_hybrid_cal_topn`, and `ng_cal_mip6` as benchmark controls. Do not promote `ng_ocs_fam10` until family-size allocation is redesigned.

## Run: `nextgen_80p_5k_3rep5cycle_ocs_lpu_frontier`

This run tuned the OCS parent-use penalty and replaced score-proportional family-size allocation with a marginal top-selection allocator. The marginal allocator assigns each extra DH progeny to the family with the largest expected excess above the current predicted threshold for selecting the top 80 out of 800 progeny.

Configuration:

- 80 parents, 5K markers, 3 phenotype reps, 5 cycles, 3 simulation replicates.
- Global historical variance calibration.
- OCS guardrail cap `10`.
- `lambda_parent_use` tested at `0`, `2`, `5`, `10`, and `20`.
- Family-size branches used `10` to `80` DH progeny per cross under fixed total progeny `800`.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lpu20` | +0.128 | +0.169 | +0.178 |
| `var_simple_ocs10_lpu10` | +0.070 | +0.074 | +0.104 |
| `ng_ocs_mip10_lpu10` | +0.061 | +0.067 | +0.099 |
| `ng_ocs_mip10_lpu0` | +0.042 | +0.035 | +0.038 |
| `ng_ocs_mip10_lpu5` | +0.040 | +0.031 | -0.026 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_ocs_fam10_lpu10` | -0.011 | -0.026 | -0.059 |
| `ng_ocs_mip10_lpu2` | -0.068 | -0.064 | -0.070 |
| `ng_ocs_fam10_lpu20` | -0.083 | -0.124 | -0.159 |
| `ng_cal_mip6` | -0.112 | -0.132 | -0.179 |
| `ng_hybrid_cal_topn` | -0.119 | -0.165 | -0.124 |
| `var_simple_etk_topn` | -0.121 | -0.153 | -0.163 |

Final cycle average over 3 replicates:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lpu20` | +0.219 | +0.122 | +0.075 |
| `ng_ocs_mip10_lpu0` | +0.166 | +0.061 | -0.037 |
| `ng_ocs_mip10_lpu10` | +0.126 | +0.055 | -0.064 |
| `ng_ocs_mip10_lpu5` | +0.023 | -0.127 | -0.238 |
| `var_simple_ocs10_lpu10` | +0.022 | -0.074 | -0.184 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `var_simple_etk_topn` | -0.040 | -0.070 | -0.123 |
| `ng_ocs_fam10_lpu10` | -0.048 | -0.163 | -0.237 |
| `ng_ocs_fam10_lpu20` | -0.073 | -0.163 | -0.221 |
| `ng_ocs_mip10_lpu2` | -0.096 | -0.124 | -0.167 |
| `ng_cal_mip6` | -0.178 | -0.250 | -0.333 |
| `ng_hybrid_cal_topn` | -0.299 | -0.436 | -0.548 |

Allocation summary averaged across cycles:

| Method | Avg Unique Parents | Avg Max Parent Use | Avg Parent-Use Sq | Avg Group Coancestry |
| --- | ---: | ---: | ---: | ---: |
| `ng_ocs_fam10_lpu20` | 16.3 | 6.8 | 0.097 | 0.179 |
| `ng_ocs_mip10_lpu20` | 15.7 | 7.3 | 0.101 | 0.175 |
| `var_simple_ocs10_lpu10` | 15.0 | 7.4 | 0.107 | 0.227 |
| `ng_cal_mip6` | 12.5 | 6.0 | 0.109 | 0.180 |
| `ng_ocs_mip10_lpu10` | 14.2 | 8.3 | 0.113 | 0.180 |
| `var_simple_topn` | 15.8 | 9.0 | 0.121 | 0.250 |
| `var_simple_etk_topn` | 15.0 | 10.9 | 0.146 | 0.268 |
| `ng_hybrid_cal_topn` | 15.7 | 11.7 | 0.156 | 0.308 |

Interpretation:

In this 80-parent run, the best setting was `ng_ocs_mip10_lpu20`. The stronger parent-use penalty did not just improve diversity; it also gave the best gain, top-tail, and max-progeny outcomes in this run. This supports the idea that the allocation objective was previously too permissive of repeated parent use and that soft contribution penalties are more useful than hard caps alone.

The marginal family-size allocator is still not validated. It improved the logic relative to score-proportional allocation, but the `10` to `80` family-size range was too aggressive and both family-size branches lost to equal family sizes.

## Run: `nextgen_80p_5k_3rep5cycle_fam_bounded`

This focused follow-up tested the same marginal family-size logic with narrower family-size bounds: `25` to `55` DH progeny per cross, still under fixed total progeny `800`.

Average across cycles 1 to 5:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_ocs_fam10_lpu20` | -0.097 | -0.170 | -0.178 |
| `ng_ocs_mip10_lpu20` | -0.180 | -0.184 | -0.114 |

Interpretation:

This run should not be compared directly to the full frontier because the global calibration pool is smaller when fewer branches are simulated. It does show that simply narrowing the family-size bounds does not rescue family-size allocation. The family-size problem needs a more conservative and better-calibrated approach, likely with shrinkage toward equal family size and a minimum predicted advantage before reallocating progeny.

Updated conclusion:

The fixed `ng_ocs_mip10_lpu20` result should be treated as an 80-parent diagnostic, not a general recommendation. Keep family-size allocation experimental and disabled for recommendations until it beats equal family sizes under the same calibration pool. The next family-size work should add shrinkage or a no-reallocation threshold rather than forcing all available progeny into variable family sizes.

## Run: `nextgen_parent_grid_screen`

This screen tested general applicability across parent counts instead of optimizing for one 80-parent case.

Configuration:

- 20, 30, 40, 50, 60, 70, and 80 parents.
- 5K markers, 3 phenotype reps, 3 breeding cycles, 2 simulation replicates.
- Effect-training size equal to parent count.
- Selected crosses scaled with parent count: 5, 8, 10, 12, 15, 18, and 20.
- Adaptive `lps` OCS penalties compared against `var_simple_topn`, `var_simple_etk_topn`, `var_simple_ocs10_lps2`, `ng_hybrid_cal_topn`, and `ng_cal_mip6`.

Best method by parent size, averaged over cycles:

| Parents | Crosses | Best Mean-Gain Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| ---: | ---: | --- | ---: | ---: | ---: |
| 20 | 5 | `ng_ocs_mip10_lps1` | +0.464 | +0.612 | +0.612 |
| 30 | 8 | `ng_ocs_mip10_lps2` | +0.075 | +0.118 | +0.039 |
| 40 | 10 | `ng_ocs_mip10_lps1` | +0.054 | +0.137 | +0.247 |
| 50 | 12 | `var_simple_etk_topn` | +0.084 | +0.046 | +0.236 |
| 60 | 15 | `var_simple_ocs10_lps2` | +0.204 | +0.289 | +0.299 |
| 70 | 18 | `var_simple_etk_topn` | +0.141 | +0.106 | +0.020 |
| 80 | 20 | `ng_hybrid_cal_topn` | +0.127 | +0.140 | -0.025 |

Average over all parent sizes:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lps1` | +0.050 | +0.085 | +0.080 |
| `ng_ocs_mip10_lps2` | +0.045 | +0.080 | +0.037 |
| `var_simple_etk_topn` | +0.045 | +0.055 | +0.064 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_hybrid_cal_topn` | -0.004 | +0.028 | +0.030 |
| `ng_ocs_mip10_lps4` | -0.009 | +0.042 | +0.070 |
| `ng_cal_mip6` | -0.011 | +0.001 | -0.020 |
| `var_simple_ocs10_lps2` | -0.023 | -0.005 | +0.027 |

Interpretation:

Adaptive `lps1` and `lps2` are more general than fixed `lpu20`, but no single method wins every parent size. The current screening frontier is `ng_ocs_mip10_lps1`, `ng_ocs_mip10_lps2`, and `var_simple_etk_topn`, with `ng_hybrid_cal_topn` kept as an unconstrained PMV-containing control.

## Run: `nextgen_parent_grid_screen_min80`

This screen repeated the parent-size grid but allowed extra individuals to estimate marker effects when the crossing-parent set was small.

Configuration:

- Same parent sizes, markers, reps, cycles, and selected-cross scaling as `nextgen_parent_grid_screen`.
- Effect-training size was `max(n_parents, 80)`.

Best method by parent size, averaged over cycles:

| Parents | Crosses | Effect Training N | Best Mean-Gain Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| ---: | ---: | ---: | --- | ---: | ---: | ---: |
| 20 | 5 | 80 | `ng_cal_mip6` | +0.083 | +0.065 | +0.003 |
| 30 | 8 | 80 | `var_simple_ocs10_lps2` | +0.519 | +0.810 | +0.859 |
| 40 | 10 | 80 | `ng_cal_mip6` | +0.211 | +0.241 | +0.264 |
| 50 | 12 | 80 | `ng_ocs_mip10_lps2` | +0.279 | +0.435 | +0.452 |
| 60 | 15 | 80 | `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| 70 | 18 | 80 | `ng_hybrid_cal_topn` | +0.222 | +0.277 | +0.464 |
| 80 | 20 | 80 | `ng_ocs_mip10_lps1` | +0.104 | +0.161 | +0.222 |

Average over all parent sizes:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_hybrid_cal_topn` | +0.064 | +0.065 | +0.080 |
| `ng_ocs_mip10_lps4` | 0.000 | +0.016 | +0.039 |
| `var_simple_ocs10_lps2` | 0.000 | +0.041 | +0.040 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_cal_mip6` | -0.003 | +0.008 | -0.033 |
| `ng_ocs_mip10_lps1` | -0.038 | -0.022 | -0.018 |
| `ng_ocs_mip10_lps2` | -0.047 | -0.043 | -0.064 |
| `var_simple_etk_topn` | -0.087 | -0.099 | -0.134 |

Interpretation:

Adding extra marker-effect training individuals changes which method wins, but it does not create a universal PMV winner. This reinforces the main design point: the software should recommend from a calibrated frontier conditional on parent count, training population size, and diversity constraints. It should not expose one hard-coded crossing-block method as state of the art.

## External Baseline Integration

Implemented package-backed baseline branches:

- PopVar: `popvar_mu_topn`, `popvar_var_topn`, `popvar_uc_topn`, `popvar_musp_topn`.
- SimpleMating: `simple_mpv_topn`, `simple_usefa_topn`, `simple_mpv_selectN`, `simple_usefa_selectN`.

The wrappers call the installed packages directly when available. PopVar is used through deterministic `pop_predict2()` with marker effects supplied from the current effect-estimation layer. SimpleMating is used through `getMPV()`, `getUsefA()`, and `selectCrosses()`.

Important runtime finding:

- Exact SimpleMating usefulness and PopVar all-pair scoring are expensive at 5K markers.
- The full 80-parent/5K all-pair exact external run did not complete within the interactive timing budget.
- The practical exact-package screen currently defaults to 40 parents, 5K markers, 1 replicate, and 1 cycle.

## Run: `external_40p_5k_all_exact`

This is a package-backed exact baseline screen, not yet a definitive validation run.

Configuration:

- 40 parents, 5K markers, 3 phenotype reps.
- 1 simulation replicate, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Compared `var_simple_topn`, PopVar superior progeny and usefulness-style scores, SimpleMating MPV/usefulness/`selectCrosses`, and adaptive `ng_ocs_mip10_lps1/lps2`.

Cycle 1 comparison versus `var_simple_topn`:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV | Max Parent Use |
| --- | ---: | ---: | ---: | ---: |
| `popvar_uc_topn` | +0.187 | -0.399 | -0.882 | 5 |
| `popvar_musp_topn` | +0.146 | -0.570 | -0.675 | 5 |
| `ng_ocs_mip10_lps1` | +0.140 | -0.412 | -0.634 | 4 |
| `simple_mpv_topn` | +0.102 | -0.593 | -1.063 | 6 |
| `ng_ocs_mip10_lps2` | +0.054 | -0.287 | -0.658 | 4 |
| `simple_usefa_topn` | +0.032 | -0.539 | -1.003 | 5 |
| `simple_usefa_select4` | +0.027 | -0.663 | -0.917 | 4 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 | 7 |

Interpretation:

This does not validate novelty by metric alone. PopVar's usefulness-style and superior-progeny criteria are competitive and slightly ahead of the current adaptive OCS branch for mean gain in this screen. `var_simple_topn` still produced the best top-tail and maximum progeny value, although with heavier parent reuse. The current novelty claim should therefore be limited to the integrated framework: modular effect reliability, calibrated hybrid PMV, adaptive OCS allocation, parent-size screening, and speed engineering. The next validation must show that this integrated system beats PopVar/SimpleMating controls across repeated parent-size scenarios, not just against `var_simple`.

## Run: `external_parent_grid_exact_5k_2rep`

This is the first replicated exact-package parent-size grid.

Configuration:

- Parent sizes: 20, 30, and 40.
- 5K markers, 3 phenotype reps.
- 2 simulation replicates, 1 breeding cycle.
- Selected crosses scaled by parent size: 5, 8, and 10.
- Exact package-backed branches: PopVar superior-progeny/usefulness-style scores, SimpleMating MPV/usefulness/`selectCrosses`, and adaptive `ng_ocs_mip10_lps1/lps2`.

Best method by parent size, averaged over the two replicates:

| Parents | Crosses | Best Mean-Gain Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| ---: | ---: | --- | ---: | ---: | ---: |
| 20 | 5 | `ng_ocs_mip10_lps1` | +0.720 | +0.915 | +1.018 |
| 30 | 8 | `ng_ocs_mip10_lps2` | +0.335 | +0.604 | +0.926 |
| 40 | 10 | `simple_usefa_select4` | +0.012 | -0.133 | -0.327 |

Average over the 20/30/40-parent screen:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lps2` | +0.248 | +0.284 | +0.343 |
| `simple_usefa_topn` | +0.238 | +0.302 | +0.320 |
| `simple_usefa_select4` | +0.212 | +0.229 | +0.198 |
| `popvar_musp_topn` | +0.208 | +0.327 | +0.402 |
| `ng_ocs_mip10_lps1` | +0.199 | +0.296 | +0.282 |
| `popvar_uc_topn` | +0.183 | +0.246 | +0.154 |
| `simple_mpv_topn` | +0.007 | +0.162 | +0.174 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |

Timing summary:

| Method | Avg Score Sec | Max Score Sec |
| --- | ---: | ---: |
| `simple_usefa_select4` | 78.6 | 160.7 |
| `popvar_uc_topn` | 76.9 | 136.7 |
| `simple_usefa_topn` | 73.8 | 142.7 |
| `popvar_musp_topn` | 69.3 | 124.3 |
| `simple_mpv_topn` | 29.5 | 49.6 |
| `ng_ocs_mip10_lps2` | 28.0 | 58.1 |
| `ng_ocs_mip10_lps1` | 26.2 | 51.6 |
| `var_simple_topn` | 23.7 | 44.3 |

Interpretation:

This run does not support the claim that PopVar or SimpleMating consistently has an accuracy advantage over the new framework. It also does not yet prove that the new framework is superior. The current evidence is head-to-head:

- `ng_ocs_mip10_lps1/lps2` won mean gain, top-10 gain, and max progeny value at 20 and 30 parents.
- SimpleMating's constrained usefulness branch had the best mean gain at 40 parents, but it was only slightly above `var_simple_topn` and lost top-tail/max value.
- PopVar's superior-progeny branch had the best 40-parent max progeny value and the best overall average max value across the 20/30/40 screen.
- Exact PopVar/SimpleMating scoring is 2.6 to 3.3 times slower than the adaptive OCS branches in this screen.

The next methodological task is to make the 50 to 80+ parent comparison feasible. Exact all-pair external usefulness is too slow at that scale. The correct next design is a two-stage benchmark:

1. fast candidate screening across all pairs;
2. exact PopVar/SimpleMating usefulness rescoring only on a large shortlist;
3. same final allocation constraints across all methods.

## Run: `external_parent_grid_shortlist_smoke`

This is the first scalable external-baseline smoke test.

Configuration:

- 50 parents, 5K markers, 3 phenotype reps.
- 1 simulation replicate, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Exact PopVar/SimpleMating scores were computed only on a 5x shortlist ranked by calibrated hybrid expected top-k value. After this smoke test, the default scalable runner was changed to use a union shortlist across calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV so the external methods are not restricted to one internal ranking.

Cycle 1 comparison versus `var_simple_topn`:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `popvar_uc_topn` | +0.437 | +0.548 | +0.399 |
| `simple_usefa_topn` | +0.365 | +0.591 | +0.403 |
| `ng_ocs_mip10_lps2` | +0.343 | +0.466 | +0.271 |
| `popvar_musp_topn` | +0.289 | +0.459 | +0.467 |
| `ng_ocs_mip10_lps1` | +0.208 | +0.515 | +0.623 |
| `simple_mpv_topn` | +0.107 | +0.370 | +0.437 |
| `simple_usefa_select4` | +0.070 | +0.009 | -0.160 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |

Interpretation:

The shortlist/rescore path is functional and makes larger parent-size external comparisons feasible. In this one-seed 50-parent smoke, PopVar UC had the best mean gain, SimpleMating usefulness had the best top-10 gain, and `ng_ocs_mip10_lps1` had the best maximum progeny value. This should not be treated as a final ranking. The value of the smoke test is that external package criteria can now be benchmarked at larger parent counts without exact all-pair scoring.

## Run: `external_parent_grid_shortlist_union_5k_1rep_m2`

This is the first large-parent union-shortlist screen. The 50-parent row comes from `external_parent_grid_shortlist_union_smoke`; the 60/70/80 rows come from `external_parent_grid_shortlist_union_5k_1rep_m2`. Settings were otherwise matched.

Configuration:

- Parent sizes: 50, 60, 70, and 80.
- 5K markers, 3 phenotype reps.
- 1 simulation replicate, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Shortlist: union of the top `2 * n_crosses` pairs from calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV.

Average over the 50/60/70/80-parent screen:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `popvar_musp_topn` | +0.451 | +0.445 | +0.489 |
| `simple_usefa_topn` | +0.443 | +0.334 | +0.325 |
| `ng_ocs_mip10_lps1` | +0.442 | +0.501 | +0.534 |
| `popvar_uc_topn` | +0.434 | +0.540 | +0.524 |
| `simple_mpv_topn` | +0.388 | +0.293 | +0.531 |
| `ng_ocs_mip10_lps2` | +0.384 | +0.440 | +0.588 |
| `simple_usefa_select4` | +0.109 | +0.244 | +0.249 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |

Best method by parent size:

| Parents | Best Mean-Gain Method | Delta Mean GV | Best Top10 Method | Delta Top10 GV | Best Max Method | Delta Max GV |
| ---: | --- | ---: | --- | ---: | --- | ---: |
| 50 | `simple_usefa_topn` | +0.654 | `ng_ocs_mip10_lps1` | +0.495 | `ng_ocs_mip10_lps1` | +0.771 |
| 60 | `ng_ocs_mip10_lps1` | +0.829 | `popvar_uc_topn` | +1.416 | `ng_ocs_mip10_lps2` | +1.698 |
| 70 | `simple_mpv_topn` | +0.166 | `popvar_uc_topn` | +0.260 | `popvar_uc_topn` | +0.711 |
| 80 | `simple_mpv_topn` | +0.590 | `ng_ocs_mip10_lps1` | +0.579 | `simple_mpv_topn` | +0.876 |

Interpretation:

This screen is a real head-to-head comparison against PopVar and SimpleMating package criteria, but it is still exploratory because it uses one replicate and shortlist rescoring. It does not show that PopVar or SimpleMating consistently dominate our framework. It also does not justify claiming that our framework consistently dominates them. The strongest current statement is that every non-`var_simple` method had positive average gain over the 50/60/70/80 screen, but individual parent sizes were mixed. The new OCS branches are competitive on top-tail and maximum progeny value while PopVar/SimpleMating remain very competitive for mean gain.

## Runtime Change: Shared Scoring

The AlphaSimR harness now caches score tables by identical parent population and calibration-history state. In cycle 1, all methods start from the same parents, so the runner now:

1. fits marker effects once;
2. computes internal cross scores once;
3. adds all requested PopVar/SimpleMating package scores once on the shortlist;
4. remaps the shared parent IDs back to each method before method-specific selection.

This keeps the selection and progeny simulation branches separate while removing duplicated all-pair scoring.

Validation:

- `shared_cache_smoke`: internal methods shared one score group.
- `shared_cache_external_smoke`: PopVar/SimpleMating methods shared one score group.
- `external_parent_grid_shared_smoke`: 50 parents, 5K markers, 1 rep, 1 cycle completed in about 3.1 minutes. The comparable pre-cache union-shortlist run took about 7.5 minutes.

## Run: `external_parent_grid_shared_5k_2rep_m2`

This is the first replicated 50/60/70/80 screen after shared scoring.

Configuration:

- Parent sizes: 50, 60, 70, and 80.
- 5K markers, 3 phenotype reps.
- 2 simulation replicates, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Shortlist: union of the top `2 * n_crosses` pairs from calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV.
- Shared scoring enabled.

Average over the 50/60/70/80-parent screen:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lps1` | +0.035 | -0.058 | -0.132 |
| `simple_usefa_topn` | +0.021 | -0.053 | -0.062 |
| `popvar_musp_topn` | +0.009 | -0.041 | -0.169 |
| `simple_mpv_topn` | +0.008 | -0.120 | -0.138 |
| `simple_usefa_select4` | +0.006 | +0.170 | +0.094 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |
| `ng_ocs_mip10_lps2` | -0.016 | -0.086 | -0.226 |
| `popvar_uc_topn` | -0.034 | -0.183 | -0.347 |

Best method by parent size:

| Parents | Best Mean-Gain Method | Delta Mean GV | Best Top10 Method | Delta Top10 GV | Best Max Method | Delta Max GV |
| ---: | --- | ---: | --- | ---: | --- | ---: |
| 50 | `simple_usefa_select4` | +0.165 | `simple_usefa_select4` | +0.397 | `simple_usefa_select4` | +0.197 |
| 60 | `var_simple_topn` | 0.000 | `simple_mpv_topn` | +0.111 | `simple_usefa_topn` | +0.452 |
| 70 | `simple_usefa_topn` | +0.332 | `simple_usefa_select4` | +0.183 | `var_simple_topn` | 0.000 |
| 80 | `simple_usefa_select4` | +0.012 | `simple_usefa_select4` | +0.158 | `simple_usefa_select4` | +0.155 |

Timing:

| Parents | Reps | Total Score Sec | Shared Score Groups |
| ---: | ---: | ---: | ---: |
| 50 | 2 | 270 | 2 |
| 60 | 2 | 382 | 2 |
| 70 | 2 | 558 | 2 |
| 80 | 2 | 524 | 2 |

Interpretation:

The stronger replicated result does not support a simple claim that PMV/PopVar/SimpleMating automatically beats `var_simple`. `var_simple_topn` remained the mean-gain winner at 60 parents and the max-progeny winner at 70 parents. SimpleMating's constrained usefulness branch was the most consistently useful external branch in this run, especially for top-tail and max progeny. The current `ng_ocs_mip10_lps1` branch had the best average mean gain across parent sizes, but the advantage was small and not enough to claim superiority.

The next methodological step is to separate metric quality from allocation quality: run PopVar, SimpleMating, and `var_simple` scores through the same OCS allocator, and run our calibrated hybrid score through the same constrained SimpleMating-style selection. That will show whether the remaining differences come from the metric or from the mate-allocation algorithm.

## Run: `allocator_crosscheck_5k_2rep_m2`

This run separates score quality from mate-allocation quality by crossing the
main score families with the same allocation styles:

- top-N ranking;
- adaptive OCS allocation;
- SimpleMating-style constrained selection with max parent use 4.

Configuration:

- Parent sizes: 50, 60, 70, and 80.
- 5K markers, 3 phenotype reps.
- 2 simulation replicates, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Shortlist: union of the top `2 * n_crosses` pairs from calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV.
- Shared scoring enabled.

Average over the 50/60/70/80-parent screen:

| Method | Delta Mean GV | Delta Top10 GV | Delta Max GV |
| --- | ---: | ---: | ---: |
| `ng_ocs_mip10_lps1` | +0.199 | +0.166 | +0.165 |
| `simple_usefa_topn` | +0.196 | +0.150 | +0.102 |
| `popvar_musp_topn` | +0.135 | +0.120 | +0.192 |
| `popvar_musp_ocs10_lps1` | +0.127 | +0.034 | -0.010 |
| `var_simple_ocs10_lps1` | +0.124 | +0.108 | +0.109 |
| `ng_hybrid_cal_topn` | +0.117 | +0.078 | +0.026 |
| `simple_usefa_ocs10_lps1` | +0.113 | +0.071 | +0.101 |
| `var_simple_select4` | +0.098 | -0.049 | -0.024 |
| `ng_hybrid_select4` | +0.084 | -0.011 | +0.047 |
| `simple_usefa_select4` | +0.042 | +0.088 | +0.162 |
| `popvar_musp_select4` | +0.011 | +0.059 | +0.138 |
| `var_simple_topn` | 0.000 | 0.000 | 0.000 |

Best method by parent size:

| Parents | Best Mean-Gain Method | Delta Mean GV | Best Top10 Method | Delta Top10 GV | Best Max Method | Delta Max GV |
| ---: | --- | ---: | --- | ---: | --- | ---: |
| 50 | `ng_ocs_mip10_lps1` | +0.058 | `var_simple_topn` | 0.000 | `var_simple_select4` | +0.439 |
| 60 | `simple_usefa_topn` | +0.272 | `popvar_musp_topn` | +0.321 | `popvar_musp_topn` | +0.543 |
| 70 | `var_simple_select4` | +0.231 | `var_simple_ocs10_lps1` | +0.328 | `simple_usefa_ocs10_lps1` | +0.488 |
| 80 | `simple_usefa_topn` | +0.382 | `ng_ocs_mip10_lps1` | +0.538 | `simple_usefa_select4` | +0.676 |

Timing:

| Parents | Reps | Total Score Sec | Shared Score Groups |
| ---: | ---: | ---: | ---: |
| 50 | 2 | 204 | 2 |
| 60 | 2 | 273 | 2 |
| 70 | 2 | 349 | 2 |
| 80 | 2 | 462 | 2 |

Interpretation:

The crossed run shows that optimizer choice materially changes the result.
`ng_ocs_mip10_lps1` had the best average mean gain, but the margin over
`simple_usefa_topn` was small. PopVar superior-progeny top-N remained very
competitive for top-tail and maximum progeny value, and SimpleMating usefulness
was the strongest external baseline for several parent sizes.

This is not yet a state-of-art win. It is stronger evidence that the framework
is competitive when the same allocation logic is applied across metrics. The
next ranking claim should come from a larger replicated parent-size grid,
including 20, 30, 40, 50, 60, 70, and 80 parents, and should report a Pareto
frontier across mean gain, top-tail gain, parent use, and coancestry.

## Family-Level Calibration Benchmark

The allocation benchmarks can hide whether a metric is statistically correct.
The family calibration runner samples a common set of candidate crosses, scores
those same crosses with all available criteria, simulates DH families in
AlphaSimR, and compares each score against realized family mean, variance,
top-10 progeny value, and maximum progeny value.

Runner:

```text
powershell -ExecutionPolicy Bypass -File nextgen_cross_design/tools/run_family_calibration_grid.ps1
```

Outputs:

- `_families.csv`: one row per simulated family with all predicted scores and realized outcomes.
- `_metric_summary_by_rep.csv`: per-replicate correlations, slopes, top-decile enrichment, and top-overlap.
- `_metric_summary_avg.csv`: average metric diagnostics across replicates.
- `_status.csv`: package status and effect-reliability diagnostics.

The runner uses cross-validated marker-effect reliability by default
(`NG_EFFECT_KFOLD=5`). This matters because in-sample reliability can be near 1
with 80 individuals and 5K markers even when cross-validated reliability is near
zero.

For genomicMateSelectR, direct inbred-parent `calcCrossLD(sire, dam)` is not
appropriate for fully inbred DH parents because each parent has identical
haplotypes and contributes zero within-parent gametic LD. The benchmark uses
genomicMateSelectR's LD/quadform calculation on the synthetic F1 haplotype
contrast and scales it to DH genotype variance.

## Run: `family_calibration_80p_5k_5rep_cv`

Configuration:

- 80 parents, 5K markers, 3 phenotype reps.
- 5 simulation replicates.
- 80 sampled candidate families per replicate.
- 80 DH progeny per sampled cross.
- Effect training set: the 80 crossing parents only.
- Cross-validated effect reliability.

Status:

- PopVar, SimpleMating, and genomicMateSelectR-derived F1/DH VPM were ok in all 5 reps.
- Mean CV reliability: 0.078.
- Mean in-sample reliability: about 0.996.
- Mean-source gate: all 400 sampled families used adjusted phenotype rather than GEBV.

Primary-target summary:

| Target | Best Score | Pearson | Spearman | Top-Overlap |
| --- | --- | ---: | ---: | ---: |
| realized mean | `cross_mean` / adjusted phenotype | 0.935 | 0.897 | 0.800 |
| realized variance | `dh_recomb_var` / PopVar var / SimpleMating var / GMS F1-DH VPM | 0.289 | 0.236 | 0.150 |
| realized top10 | `uc_dh_scaled` | 0.904 | 0.864 | 0.725 |
| realized max | `uc_dh_scaled` | 0.882 | 0.842 | 0.700 |

Variance-target contrast:

| Variance Score | Pearson | Spearman |
| --- | ---: | ---: |
| `dh_recomb_var` | +0.289 | +0.236 |
| `dh_pmv_scaled_var` | +0.281 | +0.208 |
| `hybrid_var` | +0.201 | +0.135 |
| `gated_var` | +0.023 | +0.011 |
| `var_simple` | -0.014 | -0.034 |

Interpretation:

With only 80 individuals to estimate 5K marker effects, the marker-effect model
is not reliable under cross-validation. The fallback to adjusted phenotype is
justified. In this setting, `var_simple` is not predicting realized family
variance; its previous competitiveness is more likely coming from the mean
component and allocation behavior than from correct segregation-variance
prediction.

The recombination-aware variance terms are directionally better than
`var_simple`, but variance prediction is still weak. Top-tail prediction is
high mainly because family mean is highly predictive; PMV/recombination
variance adds only marginal extra signal in the parent-only training setting.

## Run: `family_calibration_80p_5k_5rep_train400_cv`

Configuration is the same as above except the effect-training set is augmented
to 400 individuals.

Status:

- PopVar, SimpleMating, and genomicMateSelectR-derived F1/DH VPM were ok in all 5 reps.
- Mean CV reliability: 0.381.
- Mean in-sample reliability: about 0.823.
- Mean-source gate: 4 of 5 reps used GEBV; 1 of 5 still fell back to adjusted phenotype.

Primary-target summary:

| Target | Best Score | Pearson | Spearman | Top-Overlap |
| --- | --- | ---: | ---: | ---: |
| realized mean | `mpv` / PopVar mean / SimpleMating usefulness mean | 0.925 | 0.891 | 0.825 |
| realized variance | `dh_recomb_var` / PopVar var / SimpleMating var / GMS F1-DH VPM | 0.325 | 0.279 | 0.275 |
| realized top10 | `etk_dh_pmv_var_cal` | 0.895 | 0.859 | 0.700 |
| realized max | `etk_dh_pmv_var_cal` | 0.869 | 0.828 | 0.675 |

Variance-target contrast:

| Variance Score | Pearson | Spearman |
| --- | ---: | ---: |
| `dh_recomb_var` | +0.325 | +0.279 |
| `dh_pmv_scaled_var` | +0.320 | +0.271 |
| `hybrid_var` | +0.303 | +0.257 |
| `gated_var` | +0.238 | +0.195 |
| `var_simple` | -0.036 | -0.031 |

Interpretation:

Augmenting the effect-training set materially improves reliability and makes
GEBV-based means usable in most replicates. It also improves variance ranking
for recombination-aware and PMV-style metrics, while `var_simple` remains
non-informative for realized variance. This supports warning users that crossing
parents alone may be too few for marker-effect estimation, and that additional
genotyped/phenotyped individuals should be included for effect estimation when
possible.

The PMV family is statistically defensible relative to `var_simple` for the
variance target. The remaining issue is effect size: variance correlations are
only moderate, and top-tail ranking is still dominated by family mean. The next
optimization work should therefore treat PMV/recombination variance as a
secondary top-tail/diversity signal, not as a large standalone replacement for
mean-based selection.

## Follow-Up: Mean Plus Recombination Variance

The scoring code now includes explicit family-mean plus recombination-variance
criteria:

- `uc_recomb = cross_mean + i * sqrt(dh_recomb_var)`;
- `dh_recomb_var_cal`;
- `etk_dh_recomb_var_cal = cross_mean + top-k-intensity * sqrt(dh_recomb_var_cal)`.

This makes the recombination-variance path directly comparable with PMV
usefulness (`uc_dh`) and calibrated PMV expected top-k
(`etk_dh_pmv_var_cal`).

### Run: `family_calibration_mean_recomb_80p_5k_5rep_cv`

Configuration: 80 parents, 5K markers, 5 replicates, 80 sampled families per
replicate, 80 DH progeny per family, parent-only effect training.

Mean CV reliability was 0.086.

| Target | Score | Pearson | Spearman | Top-Overlap |
| --- | --- | ---: | ---: | ---: |
| realized top10 | `etk_dh_pmv_var_cal` | 0.901 | 0.877 | 0.775 |
| realized top10 | `uc_recomb` | 0.898 | 0.875 | 0.750 |
| realized top10 | `etk_dh_recomb_var_cal` | 0.898 | 0.875 | 0.750 |
| realized top10 | `popvar_musp_high` / `simple_usefa` | 0.887 | 0.851 | 0.800 |
| realized max | `etk_dh_recomb_var_cal` | 0.871 | 0.846 | 0.700 |
| realized max | `etk_dh_pmv_var_cal` | 0.873 | 0.846 | 0.725 |
| realized max | `popvar_musp_high` / `simple_usefa` | 0.866 | 0.833 | 0.775 |

Interpretation: adding mean to recombination variance makes it competitive with
mean plus PMV. In the parent-only low-reliability setting, PMV has a very small
edge for realized top10, while recombination expected top-k has the best
Spearman by a negligible margin for realized max. The practical difference is
small.

### Run: `family_calibration_mean_recomb_80p_5k_5rep_train400_cv`

Configuration is the same except effect training was augmented to 400
individuals.

Mean CV reliability was 0.345.

| Target | Score | Pearson | Spearman | Top-Overlap |
| --- | --- | ---: | ---: | ---: |
| realized top10 | `popvar_musp_high` / `simple_usefa` | 0.913 | 0.894 | 0.775 |
| realized top10 | `uc_hybrid` / `etk_hybrid_var_cal` | 0.901 | 0.882 | 0.775 |
| realized top10 | `uc_recomb` | 0.900 | 0.878 | 0.800 |
| realized top10 | `etk_dh_recomb_var_cal` | 0.900 | 0.878 | 0.800 |
| realized max | `popvar_musp_high` / `simple_usefa` | 0.890 | 0.866 | 0.775 |
| realized max | `etk_dh_pmv_scaled_var_cal` | 0.870 | 0.846 | 0.700 |
| realized max | `etk_dh_recomb_var_cal` | 0.869 | 0.845 | 0.725 |
| realized max | `uc_recomb` | 0.869 | 0.844 | 0.725 |

Interpretation: when effect training is augmented, external PopVar/SimpleMating
usefulness is strongest for top-tail ranking in this sampled-family diagnostic.
Mean plus recombination remains useful and improves top-overlap relative to
mean plus PMV in some rows, but it does not create a large advantage. This
supports a hybrid decision rule: use family mean as the dominant term, use
recombination/PMV variance as a secondary tie-breaker and diversity/top-tail
signal, and benchmark the allocator separately.

## Diagnosis: Why PopVar/SimpleMating Usefulness Looked Better

The local package code shows that PopVar and SimpleMating usefulness are both
essentially:

```text
GEBV cross mean + selection_intensity * sqrt(VPM)
```

where VPM is the point-effect progeny variance from recombination. They are not
using our `dh_pmv_var` quantity, and they are not using the adjusted phenotype
fallback for the mean when marker-effect reliability is low.

The previous comparison therefore mixed two differences:

1. Mean source: our `uc_dh`, `uc_recomb`, and calibrated `etk_*` scores used
   `cross_mean`, which can be adjusted phenotype after reliability gating.
   PopVar/SimpleMating used marker-effect GEBV mean.
2. Variance target: PopVar/SimpleMating used VPM/recombination variance. Our
   PMV term added a diagonal marker-effect uncertainty approximation. That is
   not a full posterior covariance PMV and can mostly reward uncertainty or
   segregating-marker count rather than exploitable segregation variance.

An apples-to-apples calculation on the existing 80-parent, 5K-marker,
400-effect-training run shows the mean-source issue clearly:

| Score | Spearman Top10 | Spearman Max |
| --- | ---: | ---: |
| `popvar_uc` | 0.895 | 0.868 |
| `simple_usefa` | 0.881 | 0.852 |
| `uc_recomb` using gated `cross_mean` | 0.873 | 0.847 |
| `mpv + i * sqrt(dh_recomb_var)` | 0.894 | 0.872 |
| `mpv + i * sqrt(dh_pmv_var)` | 0.896 | 0.871 |
| `mpv + i * sqrt(dh_pmv_scaled_var)` | 0.898 | 0.874 |

So the PMV usefulness concept was not truly losing once the same GEBV mean was
used. The main mistake was benchmarking a reliability-gated adjusted-phenotype
mean against external GEBV-mean usefulness and interpreting the difference as a
PMV variance failure.

The code now exposes explicit GEBV-mean usefulness variants:

- `uc_recomb_gebv`
- `uc_dh_gebv`
- `uc_dh_scaled_gebv`
- `uc_gated_gebv`
- `uc_hybrid_gebv`
- `etk_*_gebv_cal`

The next allocator benchmark should include these GEBV-mean variants beside
the reliability-gated versions. The statistical recommendation is to keep the
mean-source decision separate from the variance metric: compare GEBV-mean,
adjusted-phenotype mean, and blended-mean scores explicitly rather than hiding
that choice inside one `cross_mean` column.

## Run: `mean_source_allocator_80p_5k_5rep_train400_cv`

Date: 2026-04-28.

Configuration: 80 parents, 5K SNPs, 5 replicates, 1 breeding cycle, 3 phenotype
replicates, 400 individuals for marker-effect training, 5-fold effect CV, 10
selected crosses, 40 progeny per cross. External PopVar/SimpleMating-style
scores were benchmarked beside local recombination and PMV variants. Local
variants were run with three explicit family-mean sources:

- GEBV mean: `*_gebv_*`
- adjusted phenotype mean: `*_adj_*`
- reliability blend: `*_blend_*`

Mean effect reliability across the five replicates was 0.405.

| Method | Mean GV | Delta Mean | Top10 GV | Delta Top10 | Max GV | Delta Max | Transgressive Rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `var_simple_topn` | 3.635 | 0.000 | 4.511 | 0.000 | 4.977 | 0.000 | 0.7875 |
| `popvar_uc_topn` | 3.761 | 0.126 | 4.682 | 0.170 | 5.134 | 0.157 | 0.8750 |
| `simple_usefa_topn` | 3.765 | 0.130 | 4.659 | 0.147 | 5.227 | 0.250 | 0.8775 |
| `popvar_uc_ocs10_lps1` | 3.795 | 0.159 | 4.595 | 0.083 | 5.044 | 0.067 | 0.8925 |
| `simple_usefa_ocs10_lps1` | 3.726 | 0.090 | 4.624 | 0.113 | 5.048 | 0.071 | 0.8750 |
| `ng_recomb_blend_cal_topn` | 3.797 | 0.161 | 4.698 | 0.186 | 5.158 | 0.181 | 0.9225 |
| `ng_pmv_scaled_blend_cal_topn` | 3.811 | 0.175 | 4.755 | 0.243 | 5.334 | 0.357 | 0.9500 |
| `ng_recomb_blend_ocs10_lps1` | 3.822 | 0.187 | 4.873 | 0.362 | 5.490 | 0.513 | 0.9200 |
| `ng_pmv_scaled_blend_ocs10_lps1` | 3.858 | 0.222 | 4.755 | 0.243 | 5.288 | 0.311 | 0.9675 |

Best observed methods in this run:

- Mean GV: `ng_pmv_scaled_blend_ocs10_lps1`
- Top10 GV: `ng_recomb_blend_ocs10_lps1`
- Max GV: `ng_recomb_blend_ocs10_lps1`
- Transgressive rate: `ng_pmv_scaled_blend_ocs10_lps1`

Interpretation:

This allocator run confirms that the earlier usefulness weakness was primarily
a mean-source problem, not proof that PMV/recombination usefulness is inferior.
Once the family mean is explicit, the reliability-blended mean paired with
scaled PMV or recombination variance outperforms both external PopVar-style and
SimpleMating-style usefulness in this 80-parent scenario.

The adjusted-phenotype-only mean was not enough. The GEBV-only mean improved
several rows, but the best results came from blending GEBV and adjusted
phenotype using marker-effect reliability. That is biologically and
statistically reasonable here: with effect reliability around 0.4, pure marker
effects still carry noise, while pure adjusted phenotype throws away useful
marker information.

The practical recommendation is to keep usefulness in the form:

```text
mean_source + selection_intensity * sqrt(segregation_variance)
```

but expose `mean_source` as a deliberate modeling decision:

- use GEBV mean when marker-effect reliability is high;
- use adjusted phenotype or BLUP/BLUE mean when marker-effect reliability is
  poor;
- use a reliability blend in the intermediate zone;
- use recombination/PMV variance as a top-tail and transgressive-potential
  term, preferably inside a constrained mate-allocation optimizer rather than
  as unconstrained top-N ranking only.

This is still one scenario, not a general proof. The next validation step is a
parent-size grid across 20, 30, 40, 50, 60, 70, and 80 parents, with the same
mean-source variants and external PopVar/SimpleMating-style baselines. The
current evidence supports advancing `ng_pmv_scaled_blend_ocs` and
`ng_recomb_blend_ocs` as the primary candidates for that grid.

## Run: `mean_source_parent_grid_5k_5rep_train400`

Date: 2026-04-28.

Configuration: 20, 30, 40, 50, 60, 70, and 80 parents; 5K SNPs; 5 replicates;
1 breeding cycle; 3 phenotype replicates; 400 individuals for marker-effect
training; 5-fold effect CV; fixed 10 selected crosses; 40 progeny per cross.
The grid compared:

- `var_simple_topn`;
- PopVar/SimpleMating-style usefulness, top-N and OCS;
- local recombination expected top-k, with GEBV, adjusted, and blended means;
- local scaled-PMV expected top-k, with GEBV, adjusted, and blended means;
- OCS versions of the local recombination and scaled-PMV scores.

Mean marker-effect reliability by parent size:

| Parents | Mean Reliability |
| ---: | ---: |
| 20 | 0.542 |
| 30 | 0.439 |
| 40 | 0.434 |
| 50 | 0.449 |
| 60 | 0.422 |
| 70 | 0.343 |
| 80 | 0.379 |

Best local method versus best external method:

| Parents | Metric | Best Local | Local | Best External | External | Local - External | Local - VarSimple |
| ---: | --- | --- | ---: | --- | ---: | ---: | ---: |
| 20 | mean | `ng_pmv_scaled_gebv_cal_topn` | 3.753 | `simple_usefa_topn` | 3.766 | -0.014 | +0.184 |
| 20 | top10 | `ng_recomb_gebv_ocs10_lps1` | 4.498 | `popvar_uc_ocs10_lps1` | 4.501 | -0.003 | +0.231 |
| 20 | max | `ng_recomb_gebv_ocs10_lps1` | 4.667 | `simple_usefa_ocs10_lps1` | 4.641 | +0.025 | +0.326 |
| 30 | mean | `ng_recomb_gebv_cal_topn` | 3.942 | `simple_usefa_ocs10_lps1` | 3.930 | +0.012 | +0.203 |
| 30 | top10 | `ng_pmv_scaled_blend_ocs10_lps1` | 4.770 | `simple_usefa_ocs10_lps1` | 4.715 | +0.055 | +0.355 |
| 30 | max | `ng_pmv_scaled_gebv_ocs10_lps1` | 5.133 | `simple_usefa_ocs10_lps1` | 4.915 | +0.218 | +0.521 |
| 40 | mean | `ng_recomb_gebv_ocs10_lps1` | 4.008 | `popvar_uc_topn` | 3.997 | +0.011 | +0.198 |
| 40 | top10 | `ng_recomb_blend_ocs10_lps1` | 4.836 | `popvar_uc_topn` | 4.848 | -0.012 | +0.107 |
| 40 | max | `ng_recomb_blend_ocs10_lps1` | 5.165 | `popvar_uc_topn` | 5.175 | -0.010 | +0.135 |
| 50 | mean | `ng_recomb_gebv_ocs10_lps1` | 3.844 | `popvar_uc_topn` | 3.823 | +0.022 | +0.148 |
| 50 | top10 | `ng_recomb_gebv_ocs10_lps1` | 4.640 | `simple_usefa_topn` | 4.646 | -0.005 | +0.234 |
| 50 | max | `ng_pmv_scaled_blend_cal_topn` | 5.011 | `popvar_uc_ocs10_lps1` | 4.976 | +0.036 | +0.314 |
| 60 | mean | `ng_pmv_scaled_blend_ocs10_lps1` | 3.880 | `popvar_uc_topn` | 3.844 | +0.036 | +0.137 |
| 60 | top10 | `ng_pmv_scaled_adj_ocs10_lps1` | 4.795 | `popvar_uc_topn` | 4.718 | +0.078 | +0.265 |
| 60 | max | `ng_pmv_scaled_gebv_cal_topn` | 5.252 | `simple_usefa_topn` | 5.240 | +0.012 | +0.395 |
| 70 | mean | `ng_recomb_gebv_cal_topn` | 3.802 | `popvar_uc_ocs10_lps1` | 3.807 | -0.005 | +0.229 |
| 70 | top10 | `ng_recomb_gebv_ocs10_lps1` | 4.702 | `popvar_uc_ocs10_lps1` | 4.765 | -0.063 | +0.296 |
| 70 | max | `ng_recomb_gebv_ocs10_lps1` | 5.047 | `popvar_uc_ocs10_lps1` | 5.234 | -0.187 | +0.187 |
| 80 | mean | `ng_pmv_scaled_blend_ocs10_lps1` | 3.833 | `simple_usefa_ocs10_lps1` | 3.800 | +0.033 | +0.059 |
| 80 | top10 | `ng_pmv_scaled_blend_ocs10_lps1` | 4.915 | `popvar_uc_ocs10_lps1` | 4.802 | +0.113 | +0.112 |
| 80 | max | `ng_recomb_blend_ocs10_lps1` | 5.713 | `popvar_uc_topn` | 5.428 | +0.285 | +0.207 |

Interpretation:

The new local metrics are no longer losing to `var_simple`. The best local
method beat `var_simple` for mean, top10, and max at every parent size in this
grid. The top-tail advantage was material: local minus `var_simple` ranged
from +0.107 to +0.355 for top10 and from +0.135 to +0.521 for max.

Against PopVar/SimpleMating-style usefulness, the result is competitive rather
than dominant. Local methods won 5 of 7 parent sizes for mean GV, 3 of 7 for
top10 GV, and 5 of 7 for max GV. The losses at 20, 40, and 50 parents for
top10 were small; the 70-parent scenario was the main concern, where
`popvar_uc_ocs10_lps1` was clearly stronger for mean, top10, and max. That
scenario also had the lowest marker-effect reliability in the grid.

The practical design implication is not to lock onto one score globally. A
robust allocator should use a candidate portfolio:

- `ng_pmv_scaled_blend_ocs` as the default when reliability is moderate and
  the objective is top-tail or transgressive potential;
- `ng_recomb_gebv_ocs` or `ng_recomb_blend_ocs` when recombination variance
  is giving stronger max-family signal;
- PopVar/SimpleMating-style recombination usefulness as a retained benchmark
  and possible fallback, especially when cross-validation shows the local PMV
  calibration is weak;
- `var_simple` only as a baseline or emergency fallback, not as a primary
  design criterion.

This grid supports the mean-source fix and the blended usefulness direction,
but it also shows that the final production framework should be adaptive. The
next statistical step is to add cross-validated metric stacking or a
reliability-driven selector that learns when to weight PMV, recombination VPM,
and external-style usefulness rather than assuming one metric dominates all
parent sizes.

## Diagnostic: PMV/Usefulness Equivalence and Corrected Evaluation

Date: 2026-04-28.

The direct diagnostic in `results/usefulness_equivalence_comparisons.csv`
shows that the core recombination variance is not the source of the previous
PMV weakness:

- `dh_recomb_var` is the same additive DH recombination variance as
  SimpleMating `getUsefA()` and PopVar `pred_varG` for inbred `0/2`
  genotypes, DH generation 1, and Haldane cM positions.
- `dh_recomb_var` versus `simple_usefa_var` and `dh_recomb_var` versus
  `popvar_varG` both had Spearman `1.0`, top-overlap `1.0`, and max absolute
  variance deltas at machine precision.
- `uc_recomb_gebv`, PopVar UC, and SimpleMating UsefA have the same cross
  ranking; observed mean/UC shifts are constant intercept shifts, not ranking
  differences.
- `dh_pmv_var` is intentionally different from SimpleMating/PopVar because it
  includes marker-effect uncertainty. It should not be interpreted as the same
  formula.

The benchmark evaluator also had a fairness bug: identical selected crosses
under different methods could receive different AlphaSimR DH draws. This made
PopVar/SimpleMating appear different even when they selected the same mating
plan. The evaluator now uses a realization cache keyed by replicate, cycle,
parent population state, parent indices, and progeny count, so identical
crosses in the same population get identical realized family outcomes.

## Run: `recomb_rel_80p_5k_commonrng`

Configuration: 80 DH parents, 5K SNPs, 400 effect-training individuals, 3
phenotype reps, 1 cycle, corrected common-realization evaluator, full external
scoring with no shortlist.

Key result:

- `ng_recomb_gebv_topn`, `ng_recomb_gebv_cal_topn`, `popvar_uc_topn`, and
  `simple_usefa_topn` selected equivalent plans and had identical outcomes:
  mean GV `3.3226`, top10 GV `4.1392`, max GV `4.4468`.
- `var_simple_topn` was clearly behind: mean GV `2.7981`, top10 GV `3.6539`,
  max GV `4.0603`.
- `ng_portfolio_topn` and `ng_recomb_rel_blend_cal_topn` beat
  `var_simple_topn` but did not beat exact GEBV usefulness in this run.

Interpretation:

- The old conclusion that PMV/usefulness was worse than `var_simple` was
  confounded by benchmark realization noise and by comparing different mean
  sources.
- The validated native metric for PopVar/SimpleMating-style usefulness is
  `uc_recomb_gebv` / `ng_recomb_gebv_topn`.
- The reliability-corrected/blended variants are experimental; they are not
  yet better than exact GEBV usefulness.

## Run: 5-Cycle C++ Exact-Usefulness Allocation Checks

Configuration: 80 DH parents, 5K SNPs, 400 effect-training individuals, 3
phenotype reps, 5 cycles, native C++ recombination kernel enabled. External
packages were omitted because `ng_recomb_gebv_topn` is now verified equivalent
to PopVar/SimpleMating usefulness for this setting.

Key result:

- Exact usefulness (`ng_recomb_gebv_topn`) strongly beat `var_simple_topn`
  early and remained competitive.
- OCS on exact usefulness improved long-run behavior versus exact top-N in
  some settings, especially with stronger parent-use penalties/caps.
- `var_simple_ocs10_lps1` was still the best final-cycle method in the
  5-cycle cap sweep. This means the remaining weakness is the mate-allocation
  objective, not the recombination variance formula.

Current recommendation:

- Use `ng_recomb_gebv_topn` as the native PopVar/SimpleMating-equivalent
  usefulness baseline.
- Keep `var_simple` as a diversity/relationship baseline and diagnostic, not
  as the primary usefulness metric.
- Next improvement should target allocation: combine exact usefulness with a
  stronger diversity/parent-contribution objective and validate across parent
  sizes before claiming superiority.

## Run: `balanced_allocator_relaxed_80p_5k_5cycle_cpp`

Date: 2026-04-29.

This run added a balanced-usefulness allocation layer in
`R/04_optimizers.R`. The layer leaves the PMV/recombination metric unchanged
and builds a rank-normalized allocation objective from:

- a primary usefulness gain column, for example `uc_recomb_gebv` or
  `etk_dh_recomb_var_blend_cal`;
- a secondary diversity column, defaulting to `var_simple_cal`;
- a pair-kinship penalty;
- the existing adaptive parent-use and group-coancestry penalties.

The first implementation forced an automatic hard minimum of 24 unique parents
for 20 crosses in the 80-parent run. That was too conservative: it preserved
diversity but lost gain. The default has therefore been changed so hard
minimum unique-parent constraints are opt-in through
`NG_BALANCED_MIN_UNIQUE_PARENTS` or `NG_BALANCED_AUTO_MIN_UNIQUE=1`.

Configuration for the useful follow-up:

- 80 DH parents, 5K SNPs, 400 effect-training individuals, 3 phenotype reps,
  5 cycles, 1 simulation replicate.
- Native C++ recombination kernel enabled.
- 20 selected crosses and 40 DH progeny per cross.
- OCS cap 10, adaptive parent-use penalty suffix `_lps2`.
- Balanced diversity weight `0.30`, pair-kinship weight `0.10`, no hard
  minimum unique-parent constraint.

Average across cycles 1 to 5:

| Method | Mean GV | Top10 GV | Max GV |
| --- | ---: | ---: | ---: |
| `ng_recomb_gebv_ocs10_lps2` | 7.342 | 7.955 | 8.221 |
| `ng_recomb_blend_balanced_ocs10_lps2` | 7.299 | 7.893 | 8.123 |
| `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 7.196 | 7.774 | 8.008 |
| `ng_recomb_gebv_topn` | 7.087 | 7.587 | 7.827 |
| `var_simple_ocs10_lps1` | 7.064 | 7.616 | 7.891 |
| `ng_useful_balanced_ocs10_lps2` | 6.999 | 7.498 | 7.709 |
| `var_simple_topn` | 6.692 | 7.261 | 7.497 |

Final cycle:

| Method | Mean GV | Top10 GV | Max GV |
| --- | ---: | ---: | ---: |
| `ng_recomb_gebv_ocs10_lps2` | 10.011 | 10.254 | 10.322 |
| `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 9.734 | 10.046 | 10.198 |
| `ng_recomb_blend_balanced_ocs10_lps2` | 9.683 | 9.902 | 9.972 |
| `ng_recomb_gebv_topn` | 9.555 | 9.758 | 9.991 |
| `var_simple_ocs10_lps1` | 9.444 | 9.659 | 9.748 |
| `var_simple_topn` | 9.262 | 9.528 | 9.633 |
| `ng_useful_balanced_ocs10_lps2` | 8.986 | 9.151 | 9.287 |

Allocation summary averaged across cycles:

| Method | Unique Parents | Max Parent Use | Parent-Use Sq | Group Coancestry |
| --- | ---: | ---: | ---: | ---: |
| `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 17.2 | 7.8 | 0.099 | 0.182 |
| `var_simple_ocs10_lps1` | 15.6 | 8.4 | 0.112 | 0.264 |
| `ng_recomb_blend_balanced_ocs10_lps2` | 15.6 | 9.0 | 0.115 | 0.253 |
| `ng_recomb_gebv_ocs10_lps2` | 14.4 | 9.2 | 0.128 | 0.308 |
| `ng_useful_balanced_ocs10_lps2` | 15.8 | 9.4 | 0.116 | 0.258 |
| `var_simple_topn` | 16.0 | 9.4 | 0.120 | 0.219 |
| `ng_recomb_gebv_topn` | 14.6 | 10.4 | 0.139 | 0.415 |

Interpretation:

The clearest winner in this replicate was not the first balanced composite; it
was exact recombination usefulness inside the adaptive OCS allocator:
`ng_recomb_gebv_ocs10_lps2`. It beat both `var_simple_topn` and
`var_simple_ocs10_lps1` on average and at cycle 5. This supports the current
diagnosis that the PMV/usefulness formula is not the bottleneck; the critical
piece is the mate-allocation objective and avoiding overly hard diversity
constraints.

The balanced PMV/recombination variants were still useful. They reduced parent
concentration and group coancestry relative to exact GEBV top-N, and
`ng_recomb_blend_balanced_ocs10_lps2` was the second-best method averaged
across cycles. However, adding `var_simple_cal` directly to `uc_recomb_gebv`
as `ng_useful_balanced_ocs10_lps2` underperformed, so the diversity term should
remain a secondary tie-breaker/constraint signal rather than a replacement for
the exact usefulness score.

Next validation:

- repeat this relaxed cap-10 OCS comparison over multiple simulation
  replicates;
- run the parent-size grid for 20 to 80 parents with
  `ng_recomb_gebv_ocs10_lps2`, `ng_recomb_blend_balanced_ocs10_lps2`,
  `ng_pmv_scaled_blend_balanced_ocs10_lps2`, PopVar/SimpleMating-style
  usefulness, and `var_simple` controls;
- keep hard minimum unique-parent constraints off by default unless the user
  explicitly requests a diversity-first mating plan.

## Run: `balanced_external_parent_grid`

Date: 2026-04-29.

This run directly addresses the concern that an 80-parent result is not enough
to call a method better than PopVar or SimpleMating. It compared the current
local candidates against PopVar/SimpleMating-style usefulness over multiple
parent sizes.

Configuration:

- Parent sizes: 20, 30, 40, 50, 60, 70, and 80.
- 5K SNPs, 400 effect-training individuals, 3 phenotype reps.
- 3 AlphaSimR replicates per parent size, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Native C++ recombination kernel enabled.
- For 20, 30, and 40 parents: exact all-pair PopVar/SimpleMating scoring.
- For 50, 60, 70, and 80 parents: shortlist/exact-rescore with multiplier 20,
  using a broad shortlist from local exact usefulness, calibrated recombination,
  calibrated PMV, MPV, and `var_simple`.

Compared methods:

- local: `ng_recomb_gebv_topn`, `ng_recomb_gebv_ocs10_lps2`,
  `ng_recomb_blend_balanced_ocs10_lps2`,
  `ng_pmv_scaled_blend_balanced_ocs10_lps2`;
- external: `popvar_uc_topn`, `popvar_uc_ocs10_lps1`,
  `simple_usefa_topn`, `simple_usefa_ocs10_lps1`,
  `simple_usefa_select10`;
- controls: `var_simple_topn`, `var_simple_ocs10_lps1`.

Best local versus best external:

| Parents | Tier | Metric | Best Local | Local | Best External | External | Local - External |
| ---: | --- | --- | --- | ---: | --- | ---: | ---: |
| 20 | exact | mean | `ng_recomb_gebv_topn` | 3.538 | `popvar_uc_ocs10_lps1` | 3.538 | 0.000 |
| 20 | exact | top10 | `ng_recomb_gebv_topn` | 4.219 | `popvar_uc_ocs10_lps1` | 4.219 | 0.000 |
| 20 | exact | max | `ng_recomb_gebv_ocs10_lps2` | 4.364 | `popvar_uc_ocs10_lps1` | 4.364 | 0.000 |
| 30 | exact | mean | `ng_recomb_gebv_topn` | 3.750 | `popvar_uc_topn` | 3.750 | 0.000 |
| 30 | exact | top10 | `ng_recomb_gebv_ocs10_lps2` | 4.563 | `popvar_uc_ocs10_lps1` | 4.563 | 0.000 |
| 30 | exact | max | `ng_recomb_gebv_ocs10_lps2` | 4.742 | `popvar_uc_ocs10_lps1` | 4.742 | 0.000 |
| 40 | exact | mean | `ng_recomb_gebv_topn` | 3.842 | `popvar_uc_topn` | 3.842 | 0.000 |
| 40 | exact | top10 | `ng_recomb_gebv_ocs10_lps2` | 4.825 | `popvar_uc_ocs10_lps1` | 4.825 | 0.000 |
| 40 | exact | max | `ng_recomb_gebv_ocs10_lps2` | 5.231 | `popvar_uc_ocs10_lps1` | 5.231 | 0.000 |
| 50 | shortlist | mean | `ng_recomb_gebv_topn` | 3.828 | `popvar_uc_topn` | 3.828 | 0.000 |
| 50 | shortlist | top10 | `ng_recomb_gebv_topn` | 4.641 | `popvar_uc_topn` | 4.641 | 0.000 |
| 50 | shortlist | max | `ng_recomb_gebv_topn` | 4.877 | `popvar_uc_topn` | 4.877 | 0.000 |
| 60 | shortlist | mean | `ng_recomb_gebv_topn` | 3.407 | `popvar_uc_topn` | 3.407 | 0.000 |
| 60 | shortlist | top10 | `ng_recomb_blend_balanced_ocs10_lps2` | 4.222 | `popvar_uc_topn` | 4.122 | +0.101 |
| 60 | shortlist | max | `ng_recomb_blend_balanced_ocs10_lps2` | 4.650 | `popvar_uc_ocs10_lps1` | 4.419 | +0.231 |
| 70 | shortlist | mean | `ng_recomb_gebv_ocs10_lps2` | 4.045 | `popvar_uc_ocs10_lps1` | 4.004 | +0.041 |
| 70 | shortlist | top10 | `ng_recomb_gebv_ocs10_lps2` | 4.904 | `popvar_uc_ocs10_lps1` | 4.899 | +0.005 |
| 70 | shortlist | max | `ng_recomb_gebv_ocs10_lps2` | 5.412 | `popvar_uc_ocs10_lps1` | 5.412 | 0.000 |
| 80 | shortlist | mean | `ng_recomb_gebv_topn` | 4.226 | `popvar_uc_ocs10_lps1` | 4.272 | -0.046 |
| 80 | shortlist | top10 | `ng_recomb_gebv_topn` | 5.124 | `popvar_uc_ocs10_lps1` | 5.239 | -0.115 |
| 80 | shortlist | max | `ng_recomb_gebv_topn` | 5.598 | `popvar_uc_ocs10_lps1` | 5.650 | -0.052 |

Summary:

- Local methods beat the best PopVar/SimpleMating comparator in 4 of 21
  parent-size/objective cells.
- Local methods tied the best PopVar/SimpleMating comparator in 14 of 21
  cells.
- Local methods lost to PopVar/SimpleMating in 3 of 21 cells, all in the
  80-parent shortlist tier.
- Local methods beat the best `var_simple` control in 15 of 21 cells.

Interpretation:

This is not a win over PopVar/SimpleMating. The local exact recombination
usefulness is now correctly matching PopVar/SimpleMating usefulness in the
exact all-pair tier, which is good evidence that the metric is implemented
correctly. But matching an external baseline is not superiority.

The balanced allocation variants helped in some larger-parent shortlist cells,
especially 60-parent top-tail and max value, but they did not dominate. The
80-parent shortlist tier was a clear warning: `popvar_uc_ocs10_lps1` beat the
best local method for mean, top10, and max.

Current conclusion:

- Do not call the current framework better than PopVar or SimpleMating.
- The current framework is competitive and faster/modular, but not yet
  statistically superior.
- The next improvement should not be another small hand-tuned weight. We need
  an adaptive allocation/stacking rule that learns from reliability,
  parent-size, and calibration diagnostics when to use exact recombination
  usefulness, PMV-scaled/blended usefulness, and contribution penalties.

## Run: `adaptive_external_parent_grid`

Date: 2026-04-29.

This run added the adaptive stack to the same parent-size comparison. The stack
uses rank-normalized local scores, reliability-aware prior weights, optional
history weights, and an adaptive parent-use input for OCS. It was compared
against PopVar/SimpleMating-style usefulness and `var_simple`.

Configuration:

- Parent sizes: 20, 30, 40, 50, 60, 70, and 80.
- 5K SNPs, 400 effect-training individuals, 3 phenotype reps.
- 3 AlphaSimR replicates per parent size, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- 20/30/40 parents: exact all-pair PopVar/SimpleMating scoring.
- 50/60/70/80 parents: shortlist/exact-rescore with multiplier 20.
- Added methods: `ng_adaptive_stack_topn`, `ng_adaptive_stack_ocs10`,
  `ng_adaptive_stack_ocs10_lps1`, and `ng_adaptive_stack_ocs10_lps2`.

Best local versus best external:

| Parents | Tier | Metric | Best Local | Local | Best External | External | Local - External |
| ---: | --- | --- | --- | ---: | --- | ---: | ---: |
| 20 | exact | mean | `ng_adaptive_stack_ocs10` | 3.813 | `popvar_uc_ocs10_lps1` | 3.758 | +0.055 |
| 20 | exact | top10 | `ng_adaptive_stack_ocs10` | 4.662 | `popvar_uc_ocs10_lps1` | 4.526 | +0.136 |
| 20 | exact | max | `ng_adaptive_stack_ocs10` | 4.743 | `popvar_uc_ocs10_lps1` | 4.729 | +0.014 |
| 30 | exact | mean | `ng_recomb_gebv_topn` | 3.830 | `popvar_uc_topn` | 3.830 | 0.000 |
| 30 | exact | top10 | `ng_adaptive_stack_topn` | 4.592 | `popvar_uc_topn` | 4.545 | +0.047 |
| 30 | exact | max | `ng_adaptive_stack_topn` | 4.773 | `popvar_uc_topn` | 4.759 | +0.014 |
| 40 | exact | mean | `ng_adaptive_stack_topn` | 3.866 | `popvar_uc_topn` | 3.805 | +0.061 |
| 40 | exact | top10 | `ng_adaptive_stack_topn` | 4.486 | `popvar_uc_topn` | 4.480 | +0.006 |
| 40 | exact | max | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 4.637 | `popvar_uc_topn` | 4.608 | +0.029 |
| 50 | shortlist | mean | `ng_adaptive_stack_ocs10` | 3.823 | `popvar_uc_ocs10_lps1` | 3.806 | +0.017 |
| 50 | shortlist | top10 | `ng_adaptive_stack_ocs10` | 4.615 | `popvar_uc_ocs10_lps1` | 4.612 | +0.003 |
| 50 | shortlist | max | `ng_adaptive_stack_ocs10` | 5.088 | `popvar_uc_ocs10_lps1` | 5.088 | 0.000 |
| 60 | shortlist | mean | `ng_adaptive_stack_ocs10_lps1` | 4.168 | `popvar_uc_ocs10_lps1` | 4.106 | +0.062 |
| 60 | shortlist | top10 | `ng_adaptive_stack_ocs10_lps1` | 4.901 | `popvar_uc_ocs10_lps1` | 4.817 | +0.084 |
| 60 | shortlist | max | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.405 | `popvar_uc_ocs10_lps1` | 5.314 | +0.090 |
| 70 | shortlist | mean | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 3.876 | `popvar_uc_topn` | 3.685 | +0.192 |
| 70 | shortlist | top10 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 4.938 | `popvar_uc_topn` | 4.543 | +0.395 |
| 70 | shortlist | max | `ng_adaptive_stack_ocs10_lps2` | 5.547 | `popvar_uc_ocs10_lps1` | 4.971 | +0.576 |
| 80 | shortlist | mean | `ng_recomb_gebv_ocs10_lps2` | 4.140 | `popvar_uc_ocs10_lps1` | 4.128 | +0.012 |
| 80 | shortlist | top10 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.277 | `popvar_uc_ocs10_lps1` | 5.207 | +0.070 |
| 80 | shortlist | max | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.875 | `popvar_uc_ocs10_lps1` | 5.863 | +0.012 |

Summary:

- Local methods beat the best PopVar/SimpleMating comparator in 19 of 21
  parent-size/objective cells.
- Local methods tied PopVar/SimpleMating in the remaining 2 cells.
- Local methods lost to PopVar/SimpleMating in 0 cells.
- Local methods beat the best `var_simple` control in 17 of 21 cells, tied in
  1 cell, and lost in 3 cells.

Interpretation:

This validates the adaptive-stack direction against the external software
baselines in this AlphaSimR design. The win is not coming from copying
PopVar/SimpleMating; the best local method changes with parent size and
objective. The adaptive stack dominates many cells, while PMV-scaled balanced
usefulness becomes important in larger-parent top-tail and max-value cells.

This is still not a universal win over `var_simple`. `var_simple_ocs10_lps1`
remained best for 50-parent mean, 60-parent mean, and 60-parent top10. That
supports keeping the adjusted phenotype/BLUE fallback as an explicit part of
the adaptive logic when SNP effects are not estimated precisely.

## Run: `adaptive_guarded_failure_50_60_3rep`

Date: 2026-04-29.

This run added an explicit reliability-guarded fallback inside
`ng_adaptive_stack_score`. The stack still computes the raw adaptive consensus,
but it now blends in `etk_var_simple_cal` when cross-validated marker-effect
reliability is below 0.55. The fallback weight is bounded by
`NG_ADAPTIVE_FALLBACK_MAX_WEIGHT` and is recorded in the selection and family
outputs.

Configuration:

- Parent sizes: 50 and 60.
- 5K SNPs, 400 effect-training individuals, 3 phenotype reps.
- 3 AlphaSimR replicates, 1 breeding cycle.
- 10 selected crosses, 40 DH progeny per cross.
- Shortlist/exact-rescore with multiplier 20.

Result:

| Parents | Best Mean | Mean | Best Top10 | Top10 | Best Max | Max |
| ---: | --- | ---: | --- | ---: | --- | ---: |
| 50 | `ng_adaptive_stack_ocs10` | 4.092 | `ng_adaptive_stack_ocs10` | 4.999 | `ng_adaptive_stack_ocs10` | 5.396 |
| 60 | `ng_adaptive_stack_ocs10` | 4.012 | `var_simple_ocs10_lps1` | 4.866 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.523 |

The guarded fallback fixed the earlier 50-parent mean and 60-parent mean
losses to `var_simple`. The remaining `var_simple` edge was 60-parent top10,
and it was very small: 4.866 versus 4.864 for `ng_adaptive_stack_ocs10`.

Observed fallback weights were about 0.09 to 0.21 as reliability ranged from
0.48 down to 0.40. That is the intended behavior: enough BLUE/adjusted-pheno
signal to prevent overtrusting marker effects, but not enough to collapse the
stack into `var_simple`.

## Run: `adaptive_guarded_5cycle_20_80_1rep`

Date: 2026-04-29.

This was a 5-cycle stress smoke, not a definitive replicated validation. It
used 20, 40, 60, and 80 parents; 5K SNPs; 400 effect-training individuals; one
AlphaSimR replicate; and shortlist/exact-rescore with multiplier 20.

Average across cycles 1 to 5:

| Parents | Metric | Best Local | Local | Best External | External | Local - External |
| ---: | --- | --- | ---: | --- | ---: | ---: |
| 20 | mean | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 6.965 | `popvar_uc_ocs10_lps1` | 6.785 | +0.180 |
| 20 | top10 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 7.543 | `popvar_uc_ocs10_lps1` | 7.283 | +0.260 |
| 20 | max | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 7.700 | `popvar_uc_ocs10_lps1` | 7.455 | +0.245 |
| 40 | mean | `ng_adaptive_stack_ocs10_lps2` | 6.762 | `popvar_uc_ocs10_lps1` | 6.898 | -0.136 |
| 40 | top10 | `ng_adaptive_stack_ocs10_lps2` | 7.311 | `popvar_uc_ocs10_lps1` | 7.473 | -0.162 |
| 40 | max | `ng_adaptive_stack_ocs10_lps2` | 7.468 | `popvar_uc_ocs10_lps1` | 7.690 | -0.222 |
| 60 | mean | `ng_recomb_gebv_ocs10_lps2` | 7.213 | `popvar_uc_ocs10_lps1` | 7.310 | -0.097 |
| 60 | top10 | `ng_recomb_gebv_ocs10_lps2` | 7.782 | `popvar_uc_ocs10_lps1` | 7.847 | -0.065 |
| 60 | max | `ng_recomb_gebv_ocs10_lps2` | 7.942 | `popvar_uc_ocs10_lps1` | 8.049 | -0.107 |
| 80 | mean | `ng_recomb_gebv_ocs10_lps2` | 7.182 | `popvar_uc_ocs10_lps1` | 7.041 | +0.140 |
| 80 | top10 | `ng_recomb_gebv_ocs10_lps2` | 7.743 | `popvar_uc_ocs10_lps1` | 7.617 | +0.126 |
| 80 | max | `ng_recomb_gebv_ocs10_lps2` | 7.919 | `popvar_uc_ocs10_lps1` | 7.856 | +0.063 |

Interpretation:

The guarded fallback is useful, but it is not sufficient for robust multi-cycle
superiority. In this one-replicate stress run, local methods beat the external
baseline at 20 and 80 parents, but PopVar/SimpleMating-style usefulness beat
the best local method at 40 and 60 parents. The next work should be an
across-cycle meta-selector/portfolio allocator that can choose among
`ng_adaptive_stack`, `ng_recomb_gebv`, `ng_pmv_scaled_balanced`, and
PopVar/SimpleMating-style usefulness based on live calibration, reliability,
and parent-size context.

## Runs: `meta_method_base_5cycle_40p`, `meta_method_base_5cycle_60p`, `meta_selector_5cycle_40p`, `meta_selector_alloc_5cycle_60p`

Date: 2026-04-29.

These runs added method-family history to the meta layer. The key correction was
to exclude `ng_adaptive_score` from direct method-history learning because it is
already a composite score. Letting it compete as a base method created
self-reinforcement. The meta layer now learns over base families only:
recombination usefulness, PMV-scaled usefulness, var-simple, PopVar usefulness,
and SimpleMating usefulness.

Two meta outputs are now available:

- `ng_meta_portfolio_*`: blended rank-consensus score across the available
  base families.
- `ng_meta_selector_*`: chooses the current leading base family from realized
  method history and allocates mates with the matching framework. PMV-scaled
  blend leaders use the balanced-usefulness allocator; other leaders use the
  regular MIP contribution allocator.

Focused one-replicate 5-cycle results with 5K SNPs, 400 effect-training
individuals, 3 phenotype reps, 10 crosses, and 40 DH progeny per cross:

| Parents | Run | Best average top10 | Top10 | Best final top10 | Top10 |
| ---: | --- | --- | ---: | --- | ---: |
| 40 | `meta_method_base_5cycle_40p` | `ng_meta_portfolio_ocs10_lps2` | 7.424 | `popvar_uc_ocs10_lps1` / `simple_usefa_ocs10_lps1` | 9.530 |
| 40 | `meta_selector_5cycle_40p` | `ng_meta_selector_ocs10_lps1` / `popvar_uc_ocs10_lps1` / `simple_usefa_ocs10_lps1` | 7.500 | `ng_recomb_gebv_ocs10_lps2` | 9.540 |
| 60 | `meta_method_base_5cycle_60p` | `popvar_uc_ocs10_lps1` / `simple_usefa_ocs10_lps1` | 7.391 | `popvar_uc_ocs10_lps1` / `simple_usefa_ocs10_lps1` | 9.743 |
| 60 | `meta_selector_alloc_5cycle_60p` | `ng_meta_selector_ocs10_lps2` | 7.984 | `var_simple_ocs10_lps1` | 9.868 |

Interpretation:

The selector is more useful than a blended consensus when one method family is
clearly leading. In the allocator-aware 60-parent run, `ng_meta_selector_ocs10_lps2`
was the best average top10 method, beating PopVar/SimpleMating-style usefulness
and `var_simple_ocs10_lps1`. It did not win the final cycle, where
`var_simple_ocs10_lps1` had the best top10 value. This means the current claim
should be "promising across-cycle selector" rather than "validated universal
winner." The next validation should replicate the selector grid across 20 to 80
parents before promoting it as the default.

## Run: `fair_thread1_meta_selector_5k_1rep_20_40_60_80`

Date: 2026-04-30.

This run fixes a benchmark reproducibility problem found during selector
validation. AlphaSimR's `SimParam$new()` defaulted to 20 threads on this
machine, and same-seed runs were not reproducible. The harness now records
`NG_ALPHASIMR_THREADS` and defaults it to 1. It also seeds effect-training and
ridge fitting from a stable population/cache key instead of the score-group
order, so a method's result does not change just because other methods are
included in the panel.

Sanity checks:

- `rngcheck_thread1_pmvonly_b` and `rngcheck_thread1_pmvonly_c`: identical
  PMV-only cycle 0 and cycle 1 metrics.
- `rngcheck_seedkey_full_20p` versus `rngcheck_seedkey_pmvonly_20p`: identical
  PMV trajectory through cycles 0 to 2 whether PMV was run alone or in the full
  method panel.

Configuration:

- Parent sizes: 20, 40, 60, 80.
- 5 cycles, 1 AlphaSimR replicate.
- 5K SNPs, 400 effect-training individuals, 3 phenotype reps.
- 40 DH progeny per cross; top crosses scaled by parent count.
- `NG_USE_CPP=1`, `NG_ALPHASIMR_THREADS=1`.
- PopVar/SimpleMating shortlist/exact-rescore multiplier 20.

Average top10 across cycles 1 to 5:

| Parents | Best Method | Top10 | PMV Balanced | PopVar/Simple | Meta Portfolio | Meta Selector LPS2 | Var OCS | Var TopN |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | `ng_recomb_gebv_ocs10_lps2` | 6.009 | 5.657 | 5.810 | 5.700 | 5.981 | 5.579 | 5.043 |
| 40 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 6.806 | 6.806 | 6.074 | 6.263 | 5.811 | 6.429 | 6.261 |
| 60 | `ng_recomb_gebv_ocs10_lps2` | 7.249 | 6.563 | 7.047 | 7.041 | 7.029 | 6.464 | 6.318 |
| 80 | `ng_meta_portfolio_ocs10_lps2` | 7.340 | 7.258 | 7.274 | 7.340 | 7.212 | 6.690 | 7.324 |

Interpretation:

The previous non-thread-controlled grids should be treated as exploratory only.
Under the deterministic harness, the result is not a universal PMV or
meta-selector win. PMV-balanced clearly wins at 40 parents, recombination-GEBV
wins at 20 and 60, and the blended meta portfolio has only a small top10 edge
at 80. The selector did learn plausible leaders in several cases, but the
single-leader allocation still leaves value on the table. The next method
change should make the meta layer choose between leader allocation and portfolio
allocation based on live uncertainty, then rerun a replicated deterministic
20/30/40/50/60/70/80 grid before making any superiority claim.

## Run: `meta_gate_thread1_5k_1rep_3cyc_20_40_60_80`

Date: 2026-04-30.

This run added an uncertainty gate to `ng_meta_selector_*`. The selector now
builds both a leader plan and a portfolio plan. If leader confidence is weak
and the portfolio plan has a better meta-consensus score, the selector uses the
portfolio allocator; otherwise it keeps the leader allocator. The selection
summary records `meta_selector_decision`, `meta_selector_confidence`,
`meta_selector_leader_margin`, and `meta_selector_portfolio_score_delta`.

Configuration:

- Parent sizes: 20, 40, 60, 80.
- 3 cycles, 1 AlphaSimR replicate.
- 5K SNPs, 400 effect-training individuals, 3 phenotype reps.
- `NG_USE_CPP=1`, `NG_ALPHASIMR_THREADS=1`.
- PopVar/SimpleMating shortlist/exact-rescore multiplier 20.

Average top10 across cycles 1 to 3:

| Parents | Best Method | Top10 | PMV Balanced | Recomb-GEBV | PopVar/Simple | Meta Portfolio | Gated Selector | Var OCS | Var TopN |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | `ng_recomb_gebv_ocs10_lps2` | 4.782 | 4.576 | 4.782 | 4.570 | 4.412 | 4.402 | 4.655 | 4.125 |
| 40 | `ng_pmv_scaled_blend_balanced_ocs10_lps2` | 5.342 | 5.342 | 4.892 | 4.934 | 5.077 | 5.065 | 5.265 | 5.061 |
| 60 | `ng_recomb_gebv_ocs10_lps2` | 5.825 | 5.380 | 5.825 | 5.750 | 5.820 | 5.781 | 5.386 | 5.235 |
| 80 | `var_simple_topn` | 6.113 | 6.028 | 5.924 | 6.019 | 6.105 | 6.043 | 5.389 | 6.113 |

Selector decisions:

| Parents | Portfolio Decisions | Leader Decisions |
| ---: | ---: | ---: |
| 20 | 2 | 1 |
| 40 | 2 | 1 |
| 60 | 0 | 3 |
| 80 | 0 | 3 |

Interpretation:

The gate is useful diagnostically but is not yet a validated improvement. It
helped move the 40-parent selector toward the stronger portfolio/PMV behavior,
but it hurt the 20-parent selector where recombination-GEBV was the best top10
method. At 80 parents, the portfolio remained better than the gated selector
and was close to `var_simple_topn`. This means the software is not yet
demonstrably better than PopVar, SimpleMating, or AlphaMate-style allocation
as a general product claim. This result motivated the separate
`ng_meta_router_ocs*` implementation documented at the top of this file:
estimate live predictive accuracy/risk for each method family and choose among
complete recombination-GEBV, PMV-balanced, portfolio, external-usefulness, and
`var_simple` plans rather than using a single confidence threshold.
