# nextgenCrossDesign 0.23.0

## Breaking

* Check lines are now **references, not filters**. `exclude_threshold_violators` and
  `check_basis` are removed, and checks no longer drop or reorder candidate crosses. A check
  is a benchmark genotype (a released variety) supplied in its own `check_geno` matrix, separate
  from the candidate parents, and is never crossed (check x check would be a self, which the
  package does not enumerate).
* `ng_apply_trait_checks()` is removed; use `ng_attach_check_reference()`, which never changes
  `nrow(scores)` or its order.
* `ng_trait_check_spec()` loses its `basis` argument: the check always follows the run's own
  per-trait `mean_source`, which is what keeps the reference on the same scale as the parent
  means it is compared against.

## New

* `check_geno`, `check_pheno`, `check_records`, `check_progeny_size` on
  `ng_run_cross_prediction()`. `check_geno` is the check's own genotype matrix; `check_pheno` is
  a phenotype-shaped table for phenotypic mean sources; `check_records` is the equivalent
  pre-built structure. Which input is actually consulted is decided by the **run's own resolved
  mean source** per trait, not by the caller, so the check reference is always on the same
  scale as the parent means it is compared against. `check_progeny_size` is **required**
  whenever `trait_checks` is supplied and has no default -- progeny per family scales
  P(beat check) directly, so it must be the breeding program's own figure rather than a number
  the package invents.
* Per-trait reference columns on the candidate table: `<trait>_check_value`, `<trait>_vs_check`
  (direction-aware: positive always means better), `<trait>_check_ok`,
  `<trait>_p_beat_check`, plus `checks_all_ok` and, for multi-trait runs,
  `p_beat_all_checks` -- the probability that a progeny beats every check at once.
* `Checks` sheet in the cross-priority workbook, a horizontal (or vertical, following whichever
  axis carries the mean) check reference line on the score-versus-kinship plot, and
  `ng_plot_check_panels()` for per-trait small multiples on multi-trait runs.
* Exported the check-reference public API: `ng_align_check_geno()`,
  `ng_check_reference_value()`, `ng_attach_check_reference()`, `ng_check_tau_bounds()`,
  `ng_attach_joint_check_probability()`, `ng_check_records_from_pheno()`,
  `ng_plot_check_panels()`, `ng_check_line_value()`.

## Bug fixes

* `ng_p_superior_progeny()` compared the wrong threshold for zero-variance crosses when `tau`
  varied per cross. `mu` was subset to the zero-variance rows while `tau` was not, so each such
  cross was scored against whichever threshold sat at position 1 rather than its own --
  silently returning a wrong probability at exactly the crosses whose answer is unambiguous.
  Scalar `tau`, which every prior in-package caller passed, was unaffected.
* Stage state passed through the runner's internal `ctx` list no longer loses `NULL` fields.
  `ctx$key <- NULL` deletes the key in R, so a stage reading that field by bare name after
  `list2env()` raised an unbound-variable error rather than seeing `NULL`. Two paths were
  affected: any run with `write_outputs = TRUE` and no `trait_checks`, and single-trait runs
  reading the exact cross-trait covariance.

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
