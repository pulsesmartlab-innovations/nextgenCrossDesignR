# nextgenCrossDesign 0.28.0

## Reproducibility: a trait's result no longer depends on its row position

A trait's science depended on WHERE the breeder happened to list it in the
direction file. `ng_run_cross_prediction()`'s per-trait loop derived both of its
seeds from the trait's row position `i`:

* `seed + i - 1L` was passed to `ng_fit_ridge_effects()` and reached
  `ng_choose_ridge_lambda()`, which does `set.seed(seed); sample(...)` to build
  the k-fold CV partition used to select the ridge penalty. A different partition
  can select a different lambda, which changes the fitted marker effects and
  therefore the progeny variance, the usefulness, the index, the mid-parent PEV,
  the robust allocation and P(beat check).
* `seed + 1000L + i - 1L` seeded the per-trait posterior draws.

`ng_posterior_multitrait_cross_predict()` (R/32) had the same defect
(`seed + j` over the traits table), so the index posterior moved with the file
order too.

Measured on a 24-parent / 40-marker panel with `seed = 909`, the trait named
`disease` scored `disease_vpm` = **2.512** when it was listed FIRST and
**3.20e-07** when it was listed SECOND -- seven orders of magnitude -- and the
selected crossing plan was a different set of crosses. The two orders were both
"right" only in the sense that neither was: the split was arbitrary.

Both seeds are now derived from the trait's IDENTITY, and the two sites are
treated differently on purpose:

* **Lambda CV (`R/39`, `R/02`)**: every trait now shares ONE fold partition (the
  run's `seed`, unmodified). Folds are a nuisance parameter of lambda selection,
  not a source of innovation -- given lambda, each trait's `beta` is a
  deterministic function of its own `y` -- so sharing them introduces no
  cross-trait coupling, and it makes per-trait CV comparisons paired.
* **Posterior draws (`R/39`, `R/32`)**: traits keep DISTINCT streams, keyed on
  the trait name via the new internal `ng_trait_rng_seed()` /
  `ng_name_hash32()`. A shared stream would give every trait identical random
  innovations and manufacture cross-trait correlation in exactly the quantity
  `ng_posterior_multitrait_cross_predict()` builds. The hash is an explicit
  byte polynomial over the name's UTF-8 bytes, computed in exact double
  arithmetic, so it is identical across sessions, platforms and R versions and
  does not depend on any R hashing internal.

**This changes numbers.** In a multi-trait run, every trait after the first now
reports different variances, usefulness, index contributions, posterior
intervals and (potentially) a different crossing plan. The previous values were
position-dependent and therefore arbitrary; the new ones are the values that
trait gets in a single-trait run. Specifically:

* A single-trait run is BIT-IDENTICAL to 0.27.0 for every deterministic output
  (variances, usefulness, index, plan). Only its posterior draw stream moves,
  because the stream is now keyed on the trait's name -- a different sample from
  the same posterior, i.e. Monte Carlo noise, not a different model. A one-trait
  exception was rejected deliberately: it would make a trait's draws depend on
  how many other traits share the file, which is the same class of context
  dependence this release removes.
* In a multi-trait run, the FIRST-listed trait's deterministic columns are also
  bit-identical to 0.27.0, because sharing the run seed is exactly what position
  1 already received.

`tests/trait_order_invariance.R` pins this shut: it permutes the direction file
and asserts every per-trait column, the whole per-trait posterior table, the
index posterior and the selected plan are unchanged with `tolerance = 0`.

## Known, not fixed

`ng_estimate_genetic_covariance(method = "two_stage_ridge")` and
`ng_posterior_genetic_covariance()` (R/31) carry the same position-derived
seeding (`seed + j` over the columns of `Y`), plus residual bootstrap draws taken
from one stream in trait order. Neither has a non-test caller inside the package
-- they are opt-in diagnostic heuristics that already warn they are not formal
estimators -- and de-positioning their seeds necessarily moves the numbers in an
estimator whose accuracy guard is currently a coin flip: on unmodified 0.27.0,
`tests/genetic_covariance_estimator.R`'s Frobenius acceptance
(`ratio <= 0.5`) is met for only ~53% of base seeds (120-seed sweep: mean ratio
0.625, sd 0.346), so seed 2026 passes by luck. That estimator needs its own
stability work; fixing the seeding without it would just relabel the luck.

# nextgenCrossDesign 0.27.0

## Safety: user-supplied P and G

The phenotypic (P) and additive-genetic (G) covariance matrices that the
`economic_index` (Smith-Hazel) and `desired_gain` (Pesek-Baker) selection
indices require can now arrive from a GUI over a JSON bridge. Three changes make
that safe.

* **The covariance transform is now label-aware.** `ng_multitrait_cov_to_value_z()`
  applied each trait's direction sign and value scale by row/column ORDER and
  never consulted `dimnames`, and `ng_add_multitrait_score()` fell back to a
  POSITIONAL read whenever the supplied dimnames did not happen to contain every
  trait -- so a matrix carrying a misspelled trait label, or labels on one
  dimension only, was silently accepted in whatever order it arrived and produced
  a plausible but wrong index. (`jsonlite::fromJSON(jsonlite::toJSON(M))` returns
  an unnamed matrix, so this is not hypothetical for a GUI.) A labelled matrix is
  now REORDERED BY NAME and must name the trait set exactly: a missing, extra,
  misspelled or one-sided label is a hard error that names the offending labels.
  An unlabelled matrix is still accepted and is documented as POSITIONAL, in
  `traits$trait` order. `ng_posterior_multitrait_cross_predict()` no longer stamps
  trait names onto a user-supplied P positionally; it aligns by name too.
  Results for a correctly-ordered P/G are bit-identical.
* **`desired_gain` no longer requires P.** The Pesek-Baker coefficient solve is
  `b = G^{-1} d` with `target_matrix = G` and `projection = NULL`: P never enters
  it. `desired_gain` now solves with G alone. P is used solely to standardise the
  REPORTED predicted response by the index SD `sqrt(b' P b)`, so without P the
  coefficients and the emitted index are unchanged and only
  `desired_gain_predicted_response` degrades to `NA`, with
  `desired_gain_unavailable` / `..._unavailable_reason` naming what is missing and
  why. `economic_index` still requires BOTH P and G (`b = P^{-1} G a` genuinely
  needs P), and the refusal to substitute candidate-cross score covariance for
  quantitative-genetic covariance is unchanged for both methods.
* **`ng_estimate_genetic_covariance()` documents and guards its `Y` contract.**
  `Y` must be phenotypes or BLUEs. BLUPs are invalid as supplied
  (`Var(BLUP) = sigma2_g - PEV`, so REML sees too little genetic variance and
  almost no residual: sigma2_g biased LOW, h2 pushed toward 1, genetic
  correlations distorted when reliabilities differ) and must be deregressed first
  (Garrick-Taylor-Dekkers) or corrected by `sqrt(r_t * r_s)`. GEBVs are worse and
  circular: they are a linear function of the same markers the GRM is built from,
  so the residual goes to ~0, h2 -> 1, and G-hat becomes the covariance of
  PREDICTIONS rather than of true breeding values. A new engine-independent guard
  fits a per-trait profile REML heritability against the GRM and WARNS (never
  errors) when the residual variance is at or near zero or the implied h2 is at or
  near 1. The returned G-hat gains provenance attributes: `traits`, `n_markers`,
  `y_input_contract`, `genetic_variance`, `residual_variance`,
  `implied_heritability` (+ `_method` / `_note`), `reml_genetic_variance`,
  `reml_residual_variance`, `shrunk_input_suspected` and `shrunk_input_traits`.

## Tests

* New `tests/user_supplied_pg_contract.R` (20 checks) covers all three changes:
  permuted-but-labelled P/G give a bit-identical index, mislabelled/incomplete/
  one-sided-label matrices error by name, unlabelled matrices still work
  positionally, `desired_gain` succeeds with G alone and errors informatively
  without it, `economic_index` still requires both, and the shrunk-input guard
  fires on GEBV-like `Y` through both estimator engines while staying silent on
  phenotypes.

# nextgenCrossDesign 0.26.0

## New features

* The multi-trait selection index now has a posterior.
  `ng_posterior_multitrait_cross_predict()` was exported, documented, tested and exercised by
  the statistical release gate, but had **zero non-test callers**, so `multi_trait_score` -- the
  merit a multi-trait plan is actually ranked on -- reached the user as a point estimate with no
  uncertainty at all. Anything that wanted an uncertainty on the plan's merit had to fall back
  to a PER-TRAIT posterior, i.e. `posterior_predictions[[1]]`, which is whichever trait happens
  to sit in row 1 of the breeder's direction file. `ng_run_cross_prediction()` now calls it on a
  multi-trait run with posterior prediction on, and returns the table as
  `result$posterior_multitrait`. The new columns
  (`multi_trait_score_post_mean` / `_post_lower` / `_post_upper` / `_post_sd`, the
  `multitrait_posterior_topn_prob_<N>` family, and any requested
  `multi_trait_score_post_q<prob>`) also ride `candidate_crosses`, so they survive candidate
  filtering and allocation. Purely additive: no existing column, name or list element changes.
  Gated so a single-trait run is bit-identical and a multi-trait run with
  `run_posterior_prediction = FALSE` pays nothing.
* `ng_posterior_multitrait_cross_predict()` gains `robustness_quantile`, in the same shape
  0.25.0 gave `ng_posterior_cross_predict()`: the requested quantile (and its mirror `1 - q`) is
  cached as an extra empirical quantile of the SAME draws that produced the credible interval,
  as `multi_trait_score_post_q<prob>`, leaving `ci_level` -- and therefore the reported interval
  -- untouched. It also gains `multi_trait_score_post_sd` and a second, narrower `"posterior"`
  metadata attribute in the shape `ng_optimize_robust_mating_plan()` reads
  (`gain_col = "multi_trait_score"`, `direction = "maximize"`), so a robust allocation can now
  be run on the INDEX the plan ranks on rather than on one trait's posterior.
  **No `direction` argument** is added, deliberately: `multi_trait_score` is
  direction-normalised higher-is-better for every index method (a stated design invariant of
  `ng_add_multitrait_score()`, since both combination paths apply the trait direction upstream
  of the combination), so the orientation is hard-coded at `"maximize"` with a comment rather
  than plumbed as a knob a caller could set wrong.
* `prob_top_tier` and posterior-based `cross_confidence` are now produced on multi-trait runs.
  `ng_annotate_cross_priority_multitrait()` gains `post_sd` and `prob_top_tier`, and
  `ng_run_cross_prediction()` passes the index posterior SD and the index top-N probability the
  same way the single-trait branch has always passed its per-trait equivalents. Previously the
  multi-trait branch passed neither, so `prob_top_tier` was ABSENT from a multi-trait
  `candidate_crosses` and `confidence_method` always fell back to `"midparent_pev_index"` even
  with posterior prediction on. `prob_top_tier` is now always emitted (`NA` when no index
  posterior was run) so the multi-trait and single-trait reporting surfaces carry the same
  columns.

## Bug fixes

* `ng_add_multitrait_score()`: the soft threshold penalty divided a deficit measured on
  `threshold_column` by the robust spread of `column` -- the RANKING column. Those are different
  quantities in different units, so `trait_value_metric`, a knob with nothing to do with
  thresholds, silently rescaled how hard the thresholds bit: a measured mean penalty of 0.35
  under `trait_value_metric = "usefulness"` against 18.0 under `"pmv"` on identical data with
  identical thresholds (a 52x swing; under a variance metric the advisory soft threshold became
  de facto hard). The deficit is now scaled by the dispersion of the column it is measured on.
  `ng_multitrait_value_scale()` already falls back IQR -> MAD -> SD -> range -> 1, so a
  degenerate threshold column cannot divide by zero. When `threshold_column` IS `column` -- the
  default for every direct caller -- the two scales are the same quantity and the result is
  bit-identical to 0.25.0.
* `ng_posterior_multitrait_cross_predict()`, `value_mode = "usefulness"`: the per-trait,
  per-draw value was built as `mean + i * sigma` with **no direction sign**, so a minimize trait
  (disease, lodging) was scored on the unfavourable tail of its family -- within-family variance
  charged as a liability instead of credited as an opportunity -- and `value_z` then applied
  `-1` on top of that. Now `mean + sign * i * sigma` with `sign = -1` for a minimize trait, the
  sign applied only to the `i * sigma` term and never to the mean, matching
  `ng_run_cp_trait_value()`. This was latent in 0.25.0 (the function had no callers); wiring it
  above makes it live, which is why it is fixed first.

## Documentation

* The index posterior's **per-draw re-standardisation** is now documented where the columns are
  produced and recorded in the posterior metadata (`index_rescaling`, `index_rescaling_note`).
  `ng_add_multitrait_score()` is called inside the draw loop -- which is what makes each draw a
  coherent joint draw rather than a recombination of per-trait marginals -- and it re-derives
  its own rank-normalisation / IQR centring from the rows it is handed. Any component of
  posterior uncertainty that shifts or rescales the whole candidate pool together is therefore
  removed before the quantile is taken, so `multi_trait_score_post_lower` / `_upper` are an
  interval on a PER-DRAW RELATIVE index: narrower than a genuine index credible interval, and
  not on the same scale as the point-estimate `multi_trait_score`. Rank stability
  (`multitrait_posterior_topn_prob_*`) and conservative-tail ordering are unaffected. Removing
  the re-standardisation is a design change to the index itself and is deliberately out of scope.

# nextgenCrossDesign 0.25.0

## New features

* `ng_run_cross_prediction()` gains `ci_level` (default `0.95`) and `robustness_quantile`
  (default `NULL`). They are separate controls on purpose. `ci_level` owns the REPORTED
  posterior credible interval (`<gain>_post_lower` / `_post_upper`); `robustness_quantile`
  owns the tail a robust mate allocation is optimised against. Previously only the CI tails
  (0.025 / 0.975 at the fixed `ci_level = 0.95`) were cached, so
  `ng_optimize_robust_mating_plan()` correctly refused every other quantile -- which is every
  quantile a caller like the Shiny frontend's robustness slider actually requests -- and robust
  allocation produced no plan in practice. Serving it by bending `ci_level` to `1 - 2q` would
  have silently relabelled a reported "95% credible interval" as, say, a 50% one, so the
  requested quantile is instead extracted from the SAME posterior draws as an extra empirical
  quantile. `ng_posterior_cross_predict()` gains the matching `robustness_quantile` argument and
  caches `<gain_col>_post_q<prob>` columns (both `q` and `1 - q`, so one run serves either trait
  direction), recorded in the `"posterior"` metadata attribute. No extra sampling, no
  approximation, and the reported interval is bit-for-bit unchanged.
* `ng_optimize_robust_mating_plan()` prefers that exact cached column ahead of its lower/upper
  CI matching, so the normal-approximation path (`allow_normal_approximation = TRUE`) is now a
  genuine last resort rather than the only route to a non-default quantile. The refusal and its
  opt-in are unchanged; the error message now points at `robustness_quantile` instead of
  `ci_level`. The plan summary records `robust_tail_probability` and `robust_quantile_source`.

## Bug fixes

* `ng_optimize_robust_mating_plan()` gains `direction` (`"maximize"` / `"minimize"`, default
  `"maximize"`, accepting the package's usual direction vocabulary) and is now direction-correct.
  The ranked value is not normalised to higher-is-better: `ng_run_cp_trait_value()` returns the
  raw mean for `trait_value_metric = "mean"` and `mean + sign * i * sqrt(var)` for usefulness,
  with the sign applied only to the `i*SD` term, so for a minimize trait (disease, lodging) a
  LOWER value is better. The allocator treated `<gain_col>_post_lower` as the conservative value
  unconditionally, which for a minimize trait is the OPTIMISTIC tail -- it would have selected
  crosses on their best case and labelled the plan robust. The conservative tail is now `q` under
  `maximize` and `1 - q` under `minimize`, and because `ng_optimize_mating_plan()` maximizes the
  sum of the gain column, a minimize objective enters negated so the optimiser prefers LOWER
  robust values. This was latent: robust allocation never ran before this release.
* `ng_posterior_cross_predict()` gains the matching `direction` argument and now orients
  `posterior_topn_prob_<N>` by it: under `minimize` the stable top-N is the N SMALLEST ranked
  values, not the largest. `ng_run_cross_prediction()` derives the orientation per trait with the
  new internal `ng_run_cp_value_orientation()`, which is the trait's direction for the mean and
  usefulness metrics but always `"maximize"` for the pure-variance and parent-distance metrics
  (more within-family variance is more opportunity whichever way the trait points, and
  `ng_run_cp_trait_value()` applies no sign there). For a decrease trait scored on mean or
  usefulness, `posterior_topn_prob_<N>` -- and therefore `<trait>_post_topn` /
  `prob_top_tier` in the cross-priority risk layer -- previously ranked the WORST crosses as the
  most stable; those numbers change for such runs. Maximize traits and pure-variance metrics are
  unaffected. `ng_optimize_robust_mating_plan()` refuses a `posterior_topn_prob` objective whose
  plan direction disagrees with the direction the top-N column was built under.

## Documentation

* The robust objective's aggregation approximation is now stated where it is defined and recorded
  in the plan summary (`robustness_quantile_aggregation`, `..._note`): it sums per-cross
  posterior quantiles, which is not the quantile of the plan's total. Quantiles are additive only
  under comonotonicity; here the crosses' posteriors are positively but imperfectly correlated
  through the shared marker-effect draws, so the summed objective is a conservative bound on the
  plan total's quantile rather than a plan-level coverage statement. Computing a joint
  plan-level quantile remains out of scope.

# nextgenCrossDesign 0.24.1

## Bug fixes

* `ng_posterior_multitrait_cross_predict()`'s multivariate threshold path (`tau_lower_vec`/
  `tau_upper_vec`) no longer consumes or perturbs the caller's ambient random-number stream.
  `mvtnorm::pmvnorm()`, called once per (posterior draw, cross pair) at 3 or more traits,
  switches to a randomised algorithm that both draws from and advances `.Random.seed`; left
  unguarded, this could shift anything a caller does with the RNG after the call (a second,
  unscoped call site of the same bug fixed for the check-reference path in 0.24.0 -- see I1
  above). `p_superior_progeny_mt_post_mean`/`_lower`/`_upper` move by ~1e-5 to ~1e-4 at >= 3
  traits versus the unfixed code (measured on a fixed-seed probe), reflecting genuinely different
  Monte Carlo lattice noise rather than a bug in the new number; the 2-trait case (where pmvnorm
  is a deterministic closed form) is bit-identical before and after.

# nextgenCrossDesign 0.24.0

## Breaking

* Threshold columns are now scored against the family MEAN, never the ranking metric.
  `min_value`/`max_value` are trait-unit thresholds, but `ng_add_multitrait_score()` compared
  them against `traits$column` -- under the default `trait_value_metric = "usefulness"` that is
  `mean + i*sqrt(pmv)`, not the family mean, and under `pmv`/`vpm` it is a variance. A cross
  whose predicted mean genuinely missed a threshold could clear the (unit-mismatched) comparison
  and, under `threshold_policy = "strict"`, still get selected -- a cross that previously
  survived can now be excluded, and one that previously survived on a false pass can now be
  correctly excluded. `ng_run_cross_prediction()` sets the new `threshold_column` to
  `<trait>_mean` (or `selection_index_mean` for `index_as_trait`) automatically; a direct caller
  of `ng_multitrait_score()`/`ng_multitrait_spec()` is unaffected unless it opts in.
* `min`/`minimum` and `max`/`maximum` are now recognised as threshold-column aliases in the
  multi-trait scorer (`ng_multitrait_spec()`), matching the aliases `ng_preflight_input_tables()`
  already accepted. A spec written with an alias previously passed preflight's consistency check
  and then had its threshold silently discarded (never enforced); it is now enforced, so a run
  using an alias can newly exclude crosses that previously slipped through unfiltered.

## Bug fixes

* `<trait>_p_beat_check` and `p_beat_all_checks` now integrate shared posterior marker-effect
  uncertainty (PEV) correctly instead of raising it to the k-th power, which is only valid for
  variance components independent across a family's k progeny (D1-D3/D6). This changes the
  reported numbers materially -- e.g. 0.9997 -> 0.678 on the documented worked case -- because
  the old value systematically overstated confidence for any family with non-trivial PEV. The
  joint `p_beat_all_checks` is also now guaranteed `<=` every per-trait marginal (previously it
  could exceed a single-trait `<trait>_p_beat_check`, which is impossible for a genuine
  joint-vs-marginal pair) and integrates the SAME shared per-trait effect uncertainty the
  marginals do, via a fixed, seeded 150-draw Monte Carlo (D7). Whenever no checked trait's PEV
  is available, both paths reproduce the pre-0.23.0 numbers exactly.
* `p_beat_all_checks` is `NA_real_`, not a false-affirmative `1.0`, whenever any active checked
  trait cannot be evaluated (D3) -- integrating an unbounded `(-Inf, Inf)` region previously
  reported certainty instead of "unknown".
* The check reference line on the main score-versus-kinship plot now refuses to draw under any
  `trait_value_metric` except `"mean"` (D4), instead of ranking the check's mean-scale value
  against the candidates' usefulness-scaled rank axis -- previously wrong on ~22% of crosses per
  the QG review's simulation. The corresponding weighted-index branch also refuses (`NA`) unless
  every weighted trait is represented among the checked traits (D5), rather than silently
  imputing 0 for an unweighted trait.
* The joint check probability (`p_beat_all_checks`) no longer consumes or perturbs the caller's
  ambient random-number stream. `mvtnorm::pmvnorm()`, used at 3 or more checked traits, switches
  to a randomised algorithm that both draws from and advances `.Random.seed`; left unguarded,
  this could shift a with-checks run's later stochastic allocation (e.g.
  `optimizer = "evolution"`) away from a without-checks run's, breaking the "checks are a
  reference, never a filter" invariant. The joint probability is also now reproducible run to
  run at >= 3 traits, which it previously was not.
* `%||%` is now defined at package scope (`R/00_utils.R`). It only entered base R in 4.4.0, but
  `DESCRIPTION` declares R (>= 4.1.0) and the package used it in 15 places outside the one
  function carrying a local copy -- so the package failed to load at all on R 4.1-4.3, invisibly
  to anyone developing on a current R.

## New

* `priority_check_weight` (default `0`) on `ng_run_cross_prediction()`: lets `check_violation`
  (the count of checks a cross fails) influence a cross's priority tier, blended the same way as
  the existing score/kinship/threshold weights and never a veto. The default of `0` leaves every
  existing run's tiers unchanged; `ng_rank_cross_priority()`'s own `check_weight` argument
  (added in 0.23.0) was previously unreachable from the runner.
* `check_violation` itself (a 4th, off-by-default priority-tier component computed from the same
  per-check pass/fail matrix as `checks_all_ok`) was added in 0.23.0's development but is now
  actually wireable end to end via `priority_check_weight`, above.

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
