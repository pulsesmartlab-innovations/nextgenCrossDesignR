# nextgenCrossDesign 0.22.0 Release Candidate

Date: 2026-08-25

## Status

Version 0.22.0 is the frozen research/validation release candidate produced by the direct
installed-package quantitative-genetics audit. It is suitable for independent review and
leakage-free historical validation. It is not yet an unrestricted worldwide production
release.

Graph LD pruning and fast PMV remain enabled. Both are exercised by the release gates; neither
was removed or replaced with a slower fallback.

## Why this release is necessary

The central cross-mean and within-family variance formulas were sound, but the audit found
breeder-facing portfolio-risk inconsistencies around them:

- selected crosses and the candidate pool were normalized separately, allowing the same cross
  to receive different confidence, risk, and portfolio labels in one result;
- tied uncertainty values could be split by input row order;
- documented confidence-method labels did not match the values emitted by the runner;
- posterior sampling changed the global RNG state and could therefore perturb a later unseeded
  stochastic allocator even though portfolio reporting was not in its objective;
- the robust allocator's former 0.25 default was reconstructed by a normal approximation while
  documentation described an empirical lower credible bound.

Version 0.22.0 corrects each issue and adds installed-package checks that fail if it returns.

## Migration table

| 0.21 behavior | 0.22 behavior | Required action |
| --- | --- | --- |
| Workbook `risk_bin_within_plan` | `risk_bin_within_candidate_pool` | Update spreadsheet imports and frontend column mappings. |
| Workbook `cross_confidence_within_plan` | `cross_confidence_within_candidate_pool` | Update spreadsheet imports and frontend column mappings. |
| Selected rows were re-normalized after allocation | Candidate-pool annotations are copied unchanged to selected rows | Do not compare old and new selected risk bins as if their reference sets were identical. |
| Robust default `robustness_quantile = 0.25` | `robustness_quantile = NULL`, resolving to the cached empirical lower credible bound | Omit the argument for the statistically preferred exact default. |
| Any uncached robust quantile was silently normal-approximated | Normal approximation fails closed unless explicitly allowed | Prefer re-running posterior prediction with a matching `ci_level`; otherwise set `allow_normal_approximation = TRUE` and retain the warning/summary flag. |
| Accidental `relative_*` confidence-method labels | `posterior_ci`, `midparent_pev`, or `midparent_pev_partial` (and index equivalents) | Update code that matched the accidental strings. |

To obtain an exact empirical 25th percentile, run posterior prediction with `ci_level = 0.50` and
then request `robustness_quantile = 0.25`. The lower tail of that central interval is the desired
quantile. Using `allow_normal_approximation = TRUE` is explicitly less defensible for nonlinear
or skewed posterior merit distributions.

## Statistical interpretation

- `cross_confidence` is a within-run min-max transformation of uncertainty, not a calibrated
  probability.
- `risk_bin` is a candidate-pool relative label, not an absolute risk threshold and not
  comparable across runs. Ties remain together, so bins may be unequal.
- `prob_top_tier` remains a separate joint merit-and-uncertainty quantity and is never used as
  the confidence axis.
- Single-trait posterior confidence uses the posterior SD of the metric actually ranked.
- Multi-trait production runs retain the disclosed block-diagonal index-PEV approximation from
  independent univariate trait fits. Multi-trait posterior confidence is not claimed.
- The explicit P(top-N) robust objective maximizes the expected number of selected crosses that
  belong to the posterior top-N set. It is not the joint probability that the entire selected
  plan is top-N stable.

## Evidence frozen with this candidate

The release gate installs the package and calls its production namespace; it does not source
package test files.

- Statistical/QG gate: 44/44 PASS. Portfolio-specific checks cover the centered PEV quadratic
  form and marker alignment, monotone and tie-stable confidence, independence from P(top-N),
  exact multi-trait `sqrt(w' S w)`, exact robust lower-quantile and P(top-N) reductions,
  probability-domain guards, posterior RNG isolation, allocation non-interference, and
  selected/candidate consistency in single- and multi-trait one-call runs.
- AlphaSimR forward gate: 8/8 PASS. Predictions used 300 disjoint training-only lines and were
  completed before DH or F10 validation progeny were generated. Four additional stochastic
  founder panels passed 32/32 panel-level checks.
- Package build/check: `Status: OK`.

See `STATISTICAL_RELEASE_GATE.md`, `ALPHASIMR_FORWARD_VALIDATION.md`, and their machine-readable
CSV files in this directory.

## Reproduce the frozen gates

```sh
gate_lib=$(mktemp -d /tmp/ngcd_022_gate.XXXXXX)
R CMD INSTALL --preclean --no-multiarch --library="$gate_lib" .
NGCD_RELEASE_LIB="$gate_lib" Rscript tools/run_statistical_release_gate.R
NGCD_RELEASE_LIB="$gate_lib" Rscript tools/run_alphasimr_forward_gate.R
Rscript tools/check_r_package.R
```

## Remaining release gate

Before removing the production hold, perform a pre-registered, leakage-free validation in every
intended crop/population domain. Train using only information available before the target cycle;
predict crosses before examining progeny; then report cross-mean calibration intercept/slope,
RMSE, Pearson and Spearman correlations, within-family variance calibration, top-decile
enrichment, selected-versus-control realized response, coancestry/inbreeding, family sizes, and
uncertainty intervals by cycle and environment.

An independent quantitative geneticist should review the frozen source, this protocol, and the
field results before production approval.
