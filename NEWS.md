# nextgenCrossDesign 0.22.0 (release candidate)

## Quantitative-genetics correctness

- Hardened the diploid DH/RIL variance kernels against marker-order and malformed-input errors
  while retaining graph LD pruning and fast PMV.
- Made posterior multi-trait and disomic-subgenome scoring invariant to genotype-column order.
- Validated the installed DH and RIL-infinity production workflows against independently
  generated AlphaSimR families with progeny withheld until after prediction.

## Portfolio risk and robust allocation

- Resolve confidence normalization, risk tertiles, and portfolio medians once on the full
  post-filter candidate pool, then copy the annotations unchanged to selected crosses.
- Keep tied uncertainty values together so risk labels no longer depend on row order.
- Enforce marker-name alignment and the centered-genotype quadratic form for mid-parent PEV.
- Correct `confidence_method` provenance labels, including `posterior_ci` and
  `midparent_pev[_partial]`.
- Isolate posterior samplers from the caller's random-number stream so reporting cannot
  indirectly perturb later stochastic allocation.
- Make robust posterior allocation exact by default: `robustness_quantile = NULL` uses the
  empirical lower credible bound cached from posterior draws. An uncached normal approximation
  now requires `allow_normal_approximation = TRUE`, emits a warning, and is recorded in the plan
  summary.
- Validate P(top-N) inputs as probabilities. The corresponding objective maximizes expected
  overlap with the posterior top-N set, not a joint probability for the whole plan.

## User-visible migration

- Workbook columns `risk_bin_within_plan` and `cross_confidence_within_plan` are renamed to
  `risk_bin_within_candidate_pool` and `cross_confidence_within_candidate_pool`.
- Selected-cross confidence, risk, and portfolio labels can change relative to 0.21.0 because
  they now retain the full candidate-pool reference frame instead of being re-normalized on the
  selected subset.
- See `docs/V0_22_0_RELEASE_CANDIDATE.md` for migration and interpretation details.

## Verification status

- Installed-package statistical/quantitative-genetics gate: 44/44 checks passed.
- Controlled AlphaSimR forward gate: 8/8 checks passed, with four additional independent
  founder panels passing all 32 panel-level checks.
- Staged `R CMD build` and `R CMD check --no-manual --no-build-vignettes`: `Status: OK`.

This is a research/validation release candidate. Unrestricted worldwide production use remains
on hold pending leakage-free crop/population-specific field validation and independent
quantitative-genetics review.
