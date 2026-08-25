# Statistical Release Gate

Generated: 2026-08-25 15:45:45 CDT
Package: `nextgenCrossDesign 0.21.0`
Installed namespace: `/private/tmp/ngcd_release_gate_lib_20260825/nextgenCrossDesign`
Source Git commit: `f88da827df7e` (working tree dirty: `TRUE`)
R: `R version 4.6.1 (2026-06-24)`
Compiled LD backend available: `TRUE`

## Decision

- Mathematical/software quantitative-genetics gate: **PASS** (29/29 checks passed).
- Historical forward-validation gate: **NOT RUN**. No historical cross-by-progeny outcome table was supplied to this gate.
- Unrestricted worldwide production release: **HOLD** until the historical gate and an independent quantitative-genetics review pass.

A code-gate pass proves that the installed implementation satisfies the identities and invariants below within stated numerical tolerances. It does not prove prediction accuracy in every germplasm, crop, environment, generation, or breeding program.

## Scope and method

This runner loads an installed package and calls its production namespace, including compiled kernels. It does not source package test files. Mathematical identities use independent calculations; randomized graph checks compare exact retained-marker sets; the final check uses the public one-call workflow.

Numerical identity checks use an absolute tolerance of `1e-10`; probability identities use `1e-12`; positive-semidefinite checks allow minimum eigenvalues down to `-1e-10` for floating-point roundoff. Discrete graph, domain, and hard-constraint checks require exact agreement.

Graph LD pruning and fast PMV are retained. The gate explicitly checks them; neither is disabled.

## Results

| Domain | Check | Evidence | Acceptance | Observed | Status |
| --- | --- | --- | --- | --- | --- |
| Diploid variance | DH exact gamete enumeration | All 2^3 F1 gametic paths with exact Haldane transition probabilities | max absolute error <= 1e-10 | 8.882e-16 | PASS |
| Diploid variance | RIL-infinity Haldane-Waddington identity | Independent pairwise RIL recombination transformation | max absolute error <= 1e-10 |     0 | PASS |
| Diploid variance | Fast PMV retained: compiled recursion versus dense form | Installed compiled Haldane-DH kernel versus chromosome-block quadratic form | max absolute error <= 1e-10 | 2.22e-16 | PASS |
| Diploid variance | Diagonal posterior reduction | Full posterior covariance specialized to diagonal marker-effect covariance | max absolute error <= 1e-10 | 2.22e-16 | PASS |
| LD pruning | Graph pruning versus independent connected-component oracle | Randomized map-aware cases with missing calls, MAF filtering, and LD edges | 60/60 exact marker-set matches | 60/60 | PASS |
| LD pruning | Installed C++/auto versus R graph parity | Same randomized cases through both production backends | 60/60 exact marker-set matches | 60/60 | PASS |
| LD pruning | Map-order invariance | Marker columns permuted while chromosome/position map is held fixed | 60/60 invariant marker sets | 60/60 | PASS |
| LD pruning | Compiled graph backend loaded | Installed namespace backend registry | TRUE | TRUE | PASS |
| Marker effects | Fold-local centering translation invariance | Every marker shifted by a different constant; installed fit and CV rerun | prediction/CV difference <= 1e-10 | 1.776e-15 | PASS |
| Marker effects | Predictive diagnostics are not mislabeled reliability | Installed ridge fit metadata | reliability=NA, calibrated flag FALSE, CV R2/correlation separately reported | reliability=NA; calibrated=FALSE | PASS |
| Marker effects | Small-sample CV is not replaced by in-sample fit | Eight training records, below the package CV floor | both CV diagnostics are NA | R2=NA; cor=NA | PASS |
| Marker effects | Reliability gating requires explicit calibration | Same numeric reliability with calibration flag FALSE then TRUE | uncalibrated cannot gate; calibrated can gate | adjusted_pheno / GEBV | PASS |
| Relationships | Additive relationship, kinship, and progeny-F scale | Named relationship matrix with G[A,B]=0.4 | relationship=0.4; kinship=F=0.2 | relationship=  0.4; kinship=  0.2 | PASS |
| Relationships | Group relationship/coancestry/effective-size identity | Four equally contributing unrelated parents with G=I | c'Gc=0.25; group coancestry=0.125; Ne=4 |  0.25 / 0.125 /     4 | PASS |
| Mate allocation | Final hard-constraint audit | Installed greedy/local allocator with size, parent-use, and unique-parent constraints | 6 crosses; cap <=3; >=5 parents; final audit TRUE | n=6; cap=3; parents=6; audit=TRUE | PASS |
| Multi-trait | Smith-Hazel coefficient identity | Direct independent solve of b = P^-1 G a | max coefficient error <= 1e-10 | 5.551e-17 | PASS |
| Multi-trait | Pesek-Baker desired-gain identity | Direct independent solve of b proportional to G^-1 d | max coefficient error <= 1e-10 | 1.11e-16 | PASS |
| Multi-trait | Zero economic/desired target is preserved while thresholds remain active | Production public scorer metadata and per-trait violation diagnostics | both second-trait targets exactly zero and threshold violation detected | economic=    0; desired=    0; threshold_active=TRUE | PASS |
| Multi-trait | Formal indices fail closed without P and G | Public economic-index and desired-gain calls with covariance matrices omitted | both calls error | TRUE | PASS |
| Polyploid | Diploid reduction to Vitezica dominance coding | Independent construction of D=(-2p^2,2pq,-2q^2) basis and GRM | max absolute error <= 1e-10 |     0 | PASS |
| Polyploid | Dominance GRM positive semidefinite | Minimum eigenvalue of installed-kernel GRM | minimum eigenvalue >= -1e-10 | -4.768e-15 | PASS |
| Polyploid | Autotetraploid double-reduction probability and mean identities | All parental dosages at conventional alpha=1/6 | PMF sum, gamete mean, and progeny mean error <= 1e-12 | 8.882e-16 | PASS |
| Polyploid | Double-reduction model domain fails closed | Nonzero alpha at ploidy 6 and alpha above 1/6 at ploidy 4 | both calls error | TRUE | PASS |
| Polyploid | Additive-dominance R/C++ kernel parity and phase label | Same autotetraploid crosses through both production kernels | max absolute error <= 1e-10 and dosage-only phase-prior label present | error=8.882e-16; label=TRUE | PASS |
| Probability | At-least-one superior progeny closed form | Independent 1 - Phi((tau-mu)/sigma)^k calculation | max absolute error <= 1e-12 |     0 | PASS |
| Probability | Zero-SD deterministic boundary | Degenerate family distributions below and above threshold | exactly c(0,1) | 0,1 | PASS |
| Probability | Expected excess above threshold identity | Independent truncated-normal first moment plus deterministic limit | max absolute error <= 1e-12 |     0 | PASS |
| Probability | Invalid probability domains fail closed | Negative SD and fractional progeny count | all invalid calls error | TRUE | PASS |
| End-to-end | Installed production runner | One-call DH workflow with compiled scoring, graph LD pruning, fast PMV, and allocation | 6 finite selected crosses; fast PMV; map-aware LD; hard audit TRUE | n=6; fast=TRUE; map_aware=TRUE; audit=TRUE | PASS |

## Claims permitted by this gate

- The installed diploid DH and infinite-selfed-RIL variance implementation matches the stated quantitative-genetic formulas for the checked models.
- The fast Haldane-DH PMV recursion matches the dense quadratic form, and graph LD pruning matches its connected-component specification.
- Relationship/coancestry scaling, formal Smith-Hazel and Pesek-Baker coefficients, probability metrics, and tested optimizer constraints are internally coherent.
- The tested polyploid single-locus/dosage identities and additive-dominance R/C++ parity hold inside the explicitly reported model domain.

## Claims not permitted by this gate

- No claim of universal prediction accuracy, realized genetic gain, or superiority over external packages follows from algebraic or synthetic checks.
- No empirical recommendation is authorized for a crop, target population, generation interval, training design, or environment not represented in a forward-validation dataset.
- Polyploid additive-plus-dominance fitting remains experimental because additive and dominance components share one ridge penalty.
- Dosage-only autopolyploid within-family variance remains a uniform compatible-phase prior expectation; exact linkage-phase claims require phased homologues.
- UCPC remains experimental/internal until calibrated against independent breeding outcomes.

## Required historical gate

For each breeding population, train only on records available before the target cycle, predict candidate crosses, and compare predictions with subsequently observed progeny. At minimum report cross-mean calibration intercept/slope, RMSE, Pearson and Spearman correlation, within-family variance calibration, top-decile enrichment, selected-vs-control realized response, inbreeding/coancestry, family sizes, uncertainty intervals, and results by cycle/environment. Pre-register thresholds with breeders before examining outcomes; do not invent universal cutoffs after seeing the data.

## Reproduce

```sh
mkdir -p /tmp/ngcd_release_gate_lib
R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib .
NGCD_RELEASE_LIB=/tmp/ngcd_release_gate_lib Rscript tools/run_statistical_release_gate.R
```

Machine-readable results: `docs/STATISTICAL_RELEASE_GATE_RESULTS.csv`
