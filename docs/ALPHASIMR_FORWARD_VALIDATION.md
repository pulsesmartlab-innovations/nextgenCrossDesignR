# AlphaSimR Forward Validation

Generated: 2026-08-25 18:09:44 CDT
Installed package: `nextgenCrossDesign 0.22.0`
Installed namespace: `/private/tmp/ngcd_alphasimr_gate.LsRRIC/nextgenCrossDesign`
AlphaSimR: `2.1.0`
Source Git commit: `a1c530dffa49` (working tree dirty: `TRUE`)

## Decision

- Controlled AlphaSimR forward-validation gate: **PASS** (8/8 checks passed).
- This validates the stated additive DH and infinite-RIL models in the simulated scenario; it does not replace crop-, population-, environment-, and cycle-specific field validation.

## Leakage and fairness audit

Predictions were completed first using only parent records and 300 explicitly disjoint training-only lines. AlphaSimR progeny were then generated. No progeny phenotype, progeny genotype, QTL effect, realized family statistic, or later-cycle record was passed to `ng_run_cross_prediction()`. Every candidate cross was graded on the same independently generated DH family sample. True QTL effects were used only in the separate formula-oracle comparison and in post-prediction evaluation.

The public DH and RIL calls used the installed one-call workflow, the graph LD-pruning backend, and `method_varPMV = "fast"`; neither production feature was disabled.

## Simulation design

- runMacs GENERIC founders, 3 chromosomes, 500 segregating sites/chromosome, 20 additive QTL/chromosome, and 100 SNP/chromosome.
- 300 training-only lines, 18 candidate inbred parents, phenotype h2 = 0.7.
- All 153 DH crosses with 600 progeny/cross; 30 crosses with 600 F10 lines/cross as the finite approximation to RIL infinity.
- R seed 20260825; AlphaSimR threads fixed at 1. On AlphaSimR 2.1.0 in this environment, runMacs/MaCS founder panels are stochastic across fresh processes even with the same R seed, so exact numeric output is not claimed to be bitwise reproducible.

## Results

| Domain | Check | Acceptance | Observed | Status |
| --- | --- | --- | --- | --- |
| Forward design | Prediction precedes independent progeny realization | exact training-ID audit; no parent overlap; DH/RIL targets; fast PMV; map-aware LD; constraints TRUE | training_only=300; overlap=0; fast=TRUE; map_aware=TRUE; DH_audit=TRUE; RIL_audit=TRUE | PASS |
| AlphaSimR DH | Causal-locus cross mean | maximum absolute Monte Carlo z <= 6 | max\|z\|= 2.68 | PASS |
| AlphaSimR DH | Causal-locus within-family variance | pooled ratio in [0.94, 1.06] and Pearson >= 0.95 | ratio=0.9741; r=0.9923 | PASS |
| AlphaSimR RIL | Infinite-RIL cross mean approximated by F10 | maximum absolute Monte Carlo z <= 6 | max\|z\|=2.116 | PASS |
| AlphaSimR RIL | Infinite-RIL within-family variance approximated by F10 | pooled ratio in [0.88, 1.12] and Pearson >= 0.9 | ratio=0.9804; r=0.9962 | PASS |
| Public API | Out-of-progeny RIL-infinity ranking against F10 families | mean Spearman >= 0.35 and usefulness Spearman >= 0.1 | mean_rho=0.8843; usefulness_rho=0.774 | PASS |
| Public API | Out-of-progeny cross-mean ranking | Spearman(predicted mean, realized DH family mean) >= 0.35 | rho=0.8954 | PASS |
| Public API | Out-of-progeny usefulness ranking and top-decile enrichment | Spearman >= 0.2 and predicted top-decile enrichment > 0 | rho=0.8424; enrichment=1.057; variance_rho=0.7107 | PASS |

## Independent stochastic repeat audit

Four fresh runMacs founder panels were run through the complete gate; all 32 panel-level checks passed. Across panels, DH causal-variance ratios were 0.9955-1.0070 (family Pearson 0.9860-0.9927), F10/RIL-infinity ratios were 0.9699-1.0200 (Pearson 0.9873-0.9920), public DH mean Spearman was 0.7238-0.8658, and public DH usefulness Spearman was 0.6652-0.8409. The repeat table is `docs/ALPHASIMR_FORWARD_VALIDATION_REPEATS.csv`.

## Interpretation boundary

A pass is positive forward-simulation evidence that the installed implementation predicts the AlphaSimR meiosis model correctly at causal loci and produces useful out-of-progeny rankings in this high-information additive diploid scenario. It is not mathematical proof for all breeding programs and does not establish dominance, epistasis, GxE, low-relatedness transfer, sparse training, crop-specific maps, autotetraploids, or operational field gain. Those require separate pre-registered validations.

The RIL comparison is deliberately described as an approximation: the package target is infinite selfing, whereas AlphaSimR families here are F10 with small residual heterozygosity. The wider RIL tolerance was declared before the run.

## Reproduce

```sh
gate_lib=$(mktemp -d /tmp/ngcd_alphasimr_gate.XXXXXX)
R CMD INSTALL --preclean --no-multiarch --library="$gate_lib" .
NGCD_RELEASE_LIB="$gate_lib" Rscript tools/run_alphasimr_forward_gate.R
```

Machine-readable checks: `docs/ALPHASIMR_FORWARD_VALIDATION_RESULTS.csv`
Family-level evidence: `docs/ALPHASIMR_FORWARD_VALIDATION_FAMILIES.csv`
Independent-panel summary: `docs/ALPHASIMR_FORWARD_VALIDATION_REPEATS.csv`
