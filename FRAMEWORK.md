# Metrics and Mate-Allocation Framework

This track is designed around one principle: cross prediction and mate allocation must be evaluated as separate problems.

## What We Learn From Existing Tools

`genomicMateSelectR` is strongest on the biological variance model. Its `PredCrossVar()` formulation is based on parental haplotypes, recombination, and marker effects. This is the correct direction for PMV because cross variance is not just relationship distance. The weakness for our target use is speed and the need to distinguish immediate progeny variance from DH/RIL recombination variance.

`SimpleMating` is strongest on practical breeder-facing criteria. Its MPV, TGV, and usefulness-style scores make clear that variance alone is not the objective. The objective is usually superior progeny or multi-cycle gain while keeping a feasible crossing block.

AlphaMate is strongest on the allocation problem. It treats selection, diversity maintenance, and mating as one constrained decision problem. That is closer to a real breeding program than simply taking the top ranked crosses.

The framework must also be robust to breeding-program scale. Some programs may design a crossing block from 20 to 40 parents, while others may have 80 or more. Parameters that work only in one parent-count scenario are treated as diagnostics, not defaults.

## Why `var_simple` Keeps Winning

`var_simple` is a robust relationship-distance proxy:

```text
0.5 * (K_ii + K_jj) - K_ij
```

It has three advantages that noisy PMV variants often lack:

- it does not depend on poorly estimated SNP effects;
- it is stable in low-training, high-marker settings;
- it indirectly favors complementary/diverse parents.

PMV should beat it only when the variance model, marker effects, and target progeny type are aligned. If PMV is based on noisy effects, wrong variance target, or uncalibrated scale, `var_simple` can correctly beat it.

## Proposed Metric Set

The framework keeps all raw components and derives composite scores only at the end.

| Metric | Purpose | Notes |
| --- | --- | --- |
| `mpv` | mid-parent GEBV | only trusted when marker-effect reliability is sufficient |
| `cross_mean` | BLUE/BLUP/adjusted phenotype or GEBV mean | reliability-gated; avoids using weak GEBV blindly |
| `var_simple` | diversity/segregation proxy | retained as a baseline and fallback |
| `dh_recomb_var` | DH/RIL recombination VPM | exact Haldane chromosome recursion over the F1 haplotype contrast |
| `dh_pmv_var` | PMV with SNP-effect uncertainty | diagonal posterior uncertainty added by default |
| `uc_dh` | usefulness criterion from DH PMV | `mean + i * sqrt(PMV)` |
| `uc_gated` | reliability-gated usefulness | blends PMV and `var_simple` by effect reliability |
| `pair_kinship` | mating inbreeding risk | optimized separately, not hidden inside PMV |
| `etk_*` | expected top-k family value | uses actual family size and calibrated family variance |
| `uc_recomb` / `etk_dh_recomb_var_cal` | family mean plus DH recombination variance | explicit mean+recombination criteria |
| `popvar_*` | external PopVar criteria | deterministic/package PopVar mean, variance, and superior progeny value |
| `simple_*` | external SimpleMating criteria | package MPV, additive usefulness, and `selectCrosses` allocation |

Future exact-shortlist metrics should add:

- full haplotype `genomicMateSelectR`-style PMV;
- off-diagonal marker-effect posterior covariance;
- probability of producing at least one progeny above a threshold;
- expected top-k family value;
- multi-trait restricted selection index.

## Allocation Objective

The recommended allocation objective is:

```text
maximize sum(cross usefulness)
       - lambda_mating * mean(pair kinship)
       - lambda_group  * group coancestry(parent contributions)
       - lambda_parent_use * sum(parent contribution^2)
```

subject to:

- exactly `n_crosses` selected;
- soft parent-contribution penalties, with hard max/min crosses per parent only as guardrails;
- optional minimum unique parents;
- optional maximum pair kinship;
- optional sex/role compatibility;
- optional family-size allocation under a fixed total progeny budget.

This is an AlphaMate-style gain-diversity tradeoff, not just a pair ranking.

The current fast implementation uses a mixed-integer linear model. Parent contribution is handled with binary incremental count variables, so `sum(parent contribution^2)` is represented exactly as a convex piecewise-linear penalty. Group coancestry is handled by iterative linearization around the current contribution vector, which keeps the optimization fast while still penalizing concentrated related contributions.

Implemented allocation branches:

| Branch | Meaning |
| --- | --- |
| `ng_cal_mipN` | hard-cap MIP using calibrated hybrid expected top-k, max parent use `N` |
| `ng_ocs_mipN_lpuX` | soft optimum-contribution MIP using calibrated hybrid expected top-k, guardrail cap `N`, absolute parent-use penalty `X` |
| `ng_ocs_mipN_lpsX` | same OCS MIP, but `X` is an adaptive multiplier scaled by score spread and crossing-block size |
| `ng_ocs_famN_lpuX` / `ng_ocs_famN_lpsX` | same OCS MIP plus experimental variable family sizes under fixed total progeny |
| `var_simple_ocsN_lpuX` / `var_simple_ocsN_lpsX` | OCS controls using calibrated `var_simple` expected top-k |
| `var_simple_ocs_famN_lpuX` / `var_simple_ocs_famN_lpsX` | `var_simple` OCS controls plus experimental variable family sizes |

For general use, prefer the `lps` adaptive scale over fixed `lpu` penalties. The current adaptive form is:

```text
actual lambda_parent_use = lps_multiplier * n_crosses * robust_score_scale
```

where `robust_score_scale` is estimated from the high-value candidate pool. This prevents an 80-parent tuning value from being reused blindly in 20-parent or 50-parent crossing blocks.

Family-size allocation currently supports two modes:

- `score_weighted`: proportional allocation from a cross score;
- `marginal_topk`: greedy allocation by expected excess above the predicted global selection threshold.

The marginal allocator is biologically more defensible than score weighting, but it is still experimental. Current AlphaSimR results show that varying family size from noisy cross predictions can lose gain even when the mating plan itself is strong.

## Benchmark Discipline

Every AlphaSimR benchmark should report three classes of results:

1. Metric calibration:
   correlation and slope between predicted variance and realized family variance.
2. Ranking value:
   top-N family mean, top-10 progeny value, and max progeny value using the same optimizer.
3. Allocation value:
   gain, group coancestry, pair kinship, unique parent count, and parent-use distribution using the same metric.

If PMV loses, these diagnostics should say whether the failure is from effect estimation, variance prediction, candidate screening, or the optimizer.

The dedicated family-level calibration runner is:

```text
powershell -ExecutionPolicy Bypass -File nextgen_cross_design/tools/run_family_calibration_grid.ps1
```

It should be run before allocator tuning. Its main outputs are per-metric
correlations, calibration slopes, top-decile enrichment, and top-overlap against
the same realized families. For inbred DH parents, genomicMateSelectR is used
through a synthetic F1 haplotype contrast rather than direct inbred-parent
`calcCrossLD`, because direct inbred-parent gametic LD is zero.

Every promoted method should also pass a parent-size screen. The default grid uses 20, 30, 40, 50, 60, 70, and 80 parents, scales the number of selected crosses roughly as `round(n_parents / 4)` with practical bounds, and compares both:

- same-size effect training, where the effect-training set is only the crossing parents;
- augmented effect training, where extra individuals can be used only to estimate marker effects.

The expected recommendation is a Pareto frontier or decision rule, not a universal single winner.

## External Baselines

PopVar and SimpleMating are now explicit benchmark branches rather than informal references.

Implemented branches:

| Branch | Meaning |
| --- | --- |
| `popvar_mu_topn` | PopVar predicted cross mean |
| `popvar_var_topn` | PopVar predicted genetic variance |
| `popvar_uc_topn` | PopVar mean plus usefulness-style variance term |
| `popvar_musp_topn` | PopVar superior progeny mean |
| `simple_mpv_topn` | SimpleMating mid-parental value |
| `simple_usefa_topn` | SimpleMating additive usefulness |
| `simple_mpv_selectN` | SimpleMating `selectCrosses` using MPV with max parent use `N` |
| `simple_usefa_selectN` | SimpleMating `selectCrosses` using additive usefulness with max parent use `N` |
| `var_simple_selectN` | SimpleMating-style constrained selection using calibrated `var_simple` expected top-k |
| `ng_hybrid_selectN` | SimpleMating-style constrained selection using calibrated hybrid expected top-k |
| `popvar_musp_selectN` | SimpleMating-style constrained selection using PopVar superior progeny value |
| `popvar_uc_selectN` | SimpleMating-style constrained selection using PopVar usefulness criterion |
| `popvar_musp_ocsN_lpsX` | adaptive OCS allocation using PopVar superior progeny value |
| `popvar_uc_ocsN_lpsX` | adaptive OCS allocation using PopVar usefulness criterion |
| `popvar_mu_ocsN_lpsX` | adaptive OCS allocation using PopVar predicted cross mean |
| `simple_mpv_ocsN_lpsX` | adaptive OCS allocation using SimpleMating MPV |
| `simple_usefa_ocsN_lpsX` | adaptive OCS allocation using SimpleMating additive usefulness |

These branches use the installed packages when available. They are not substitutes for the main framework; they are validation controls. A claim of novelty requires beating or clearly improving on these baselines in at least one defensible dimension:

- realized gain or top-tail value under the same AlphaSimR design;
- parent-use and coancestry at comparable gain;
- runtime/scalability for 5K-marker, 80-parent or larger crossing blocks;
- robustness across parent counts and effect-training sizes.

The exact external parent-size screen should start with 20, 30, and 40 parents. For 50 to 80+ parents, exact all-pair SimpleMating usefulness should be replaced by a shortlist/exact-rescore design or by formula-compatible fast surrogates before running replicated grids. The current shortlist is a union across calibrated hybrid expected top-k, hybrid usefulness, `var_simple`, and MPV, which avoids benchmarking PopVar/SimpleMating inside only one internal score ranking. The current shortlist runner is:

```text
powershell -ExecutionPolicy Bypass -File nextgen_cross_design/tools/run_external_shortlist_parent_size_grid.ps1
```

The benchmark harness uses shared scoring by default. In the first cycle, all
methods have the same parent population, so fitting effects and computing
5K-marker all-pair scores separately for every method wastes time and can
confound runtime comparisons. Shared scoring computes the score table once per
identical parent-population/calibration-history group and remaps parent labels
before method-specific selection.

To separate score quality from allocation quality, use the allocator crosscheck
runner. It compares each score family under top-N ranking, adaptive OCS, and
SimpleMating-style constrained selection where available:

```text
powershell -ExecutionPolicy Bypass -File nextgen_cross_design/tools/run_allocator_crosscheck_grid.ps1
```
