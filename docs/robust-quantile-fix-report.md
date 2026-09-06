# Robust posterior allocation: arbitrary quantile + direction correctness (0.25.0)

Backend repo (`nextgenCrossDesign`), branch `main`. Version bumped 0.24.1 -> 0.25.0.

## The defect, restated

The Shiny frontend orchestrates robust allocation itself: it runs a prediction, then calls
`ng_optimize_robust_mating_plan()` with the breeder's `robustness_quantile` slider
(0.05-0.50, default 0.25).

`ng_posterior_cross_predict()` cached exactly two empirical tails of the ranked value, at
`alpha = (1 - ci_level)/2` and `1 - alpha`, with `ci_level` fixed at 0.95 because
`ng_run_cross_prediction()` did not expose it. Every quantile the slider can produce therefore
missed both cached tails, and the allocator refused it -- correctly, since it declines to fake a
quantile without `allow_normal_approximation`. Net effect: robust allocation silently yielded no
plan for real users.

## Requirement 1 -- serve any robustness quantile exactly, without collateral damage

### (a) `ci_level` is now an explicit argument of `ng_run_cross_prediction()`

Default `0.95`, validated in `ng_cp__build_ctx()` as one finite probability in `(0, 1)` -- the
same rule `ng_posterior_cross_predict()` already applied. Echoed in `result$settings$ci_level`
and documented in `docs/frontend/contracts/config_schema.json` (group `posterior`).

### (b) The robustness quantile is threaded through and cached from the same draws

Setting `ci_level = 1 - 2q` -- which the old error message suggested -- was deliberately NOT the
fix. `ci_level` has two jobs: it sets the robust tail AND the reported credible interval. Using
it to reach a 0.25 robust quantile would turn every reported "95% credible interval" into a 50%
one: a displayed statistic corrupted to serve an unrelated internal need.

Instead, `ng_posterior_cross_predict()` gains `robustness_quantile`. `uc_mat` already holds the
per-draw ranked values, so one more `row_quantile()` call per requested probability costs no
extra sampling and involves no approximation. The columns are named
`<gain_col>_post_q<prob>` (decimal point rendered as `_`, so `0.25` -> `..._post_q0_25` and
`0.025` -> `..._post_q0_025`; `formatC(..., digits = 10)` absorbs floating-point noise such as
`(1 - 0.95)/2 = 0.02500000000000002`). Both `q` and `1 - q` are cached from a single request, so
one prediction run serves a robust allocation in either trait direction (see Requirement 2). The
probability-to-column map is recorded in the `"posterior"` metadata attribute as
`posterior_quantiles` (a `prob` / `column` data frame), alongside `robustness_quantile` and
`direction`.

`ng_run_cross_prediction()` passes its own `robustness_quantile` (validated as one probability in
`(0, 1)`, `NULL` = off) straight through, per trait.

### (c) The optimiser prefers the exact cached column

`ng_optimize_robust_mating_plan()` resolves the tail in this order:

1. exact cached `<gain_col>_post_q<prob>` column, located via the metadata map (with a
   name-convention fallback for tables written before the metadata existed);
2. `<gain_col>_post_lower` when the needed tail coincides with the CI's lower tail;
3. `<gain_col>_post_upper` when it coincides with the upper one;
4. the normal approximation -- unchanged, still gated on `allow_normal_approximation = TRUE`,
   still warned about, still flagged in the plan summary;
5. otherwise the original hard error.

`allow_normal_approximation` and its refusal are intact. The refusal message now points at
`robustness_quantile` on `ng_posterior_cross_predict()` / `ng_run_cross_prediction()` and says
explicitly that this route leaves `ci_level`, and therefore the reported interval, untouched. The
summary gained `robust_tail_probability` and `robust_quantile_source` so an auditor can see which
tail was used and where it came from.

## Requirement 2 -- direction correctness

### The analysis was correct. Evidence

`ng_run_cp_trait_value()` (`R/39_cross_prediction_runner.R:508-543`) is not normalised to
higher-is-better:

* `trait_value_metric = "mean"` returns `cross_mean_blend` raw, no sign;
* usefulness returns `mean_value + sign * i * sqrt(var)` with
  `sign <- if (identical(direction, "maximize")) 1 else -1` applied only to the `i*SD` term.

So for a minimize trait the ranked value is `mean - i*SD` and a LOWER value is better. The
package's own test asserts exactly this for a decrease trait
(`tests/user_cross_prediction_posterior_parameters.R`):

```r
disease_expected <- result$candidate_crosses$disease_mean -
  i * sqrt(pmax(result$candidate_crosses$disease_pmv_full_posterior, 0))
stopifnot(max(abs(result$candidate_crosses$disease_value - disease_expected)) < 1e-8)
```

`ng_optimize_robust_mating_plan()` had no `direction` argument and took
`<gain_col>_post_lower` as the conservative value unconditionally. For a minimize trait the lower
tail is the OPTIMISTIC case, so the old code would have ranked crosses by their best case and
called the plan robust. It never surfaced because robust allocation never ran.

The main (non-robust) allocation path is unaffected: `ng_multitrait_score()`
(`R/19_multi_trait_selection.R:634`) applies `sign_i <- if (maximize) 1 else -1` when building
the index the allocator actually optimises. The robust path bypassed that, which is why only it
was wrong.

### How the optimiser ranks, and how minimize is expressed

Read before deciding: `ng_optimize_mating_plan()` (`R/04_optimizers.R`) folds every per-cross
penalty into `scores$.linear_gain <- scores[[gain_col]] - lambda_mating * pair_kinship - ...`
(line 277) and every allocator MAXIMIZES the sum of it -- greedy takes
`order(scores$.linear_gain, decreasing = TRUE)`, the MIP passes `objective.in = scores$.linear_gain`
with `max`, and the local search accepts swaps that raise `sum(.linear_gain[selected])`.

So minimize is expressed by negation:

* conservative tail probability = `q` under `maximize`, `1 - q` under `minimize`;
* `.robust_gain <- (if (maximize) 1 else -1) * conservative_value`, so maximizing the objective
  prefers LOWER robust values;
* the un-negated conservative value stays on the table as `.robust_gain_value`, and the summary
  reports `robust_total_value` / `robust_mean_value` on the native scale plus
  `robust_objective_is_negated`, so a reader never has to know about the sign convention.

`direction` defaults to `"maximize"`, so existing callers are unchanged. It accepts the package's
existing vocabulary by reusing `ng_multitrait_direction()` (max/maximize/increase/higher/+ and
min/minimize/decrease/lower/-), not a new normaliser.

### Adjacent defect found and fixed: `posterior_topn_prob` orientation

`posterior_topn_prob_<N>` was computed as `sort(col, decreasing = TRUE)[N]` regardless of
direction. For a minimize trait that counts the fraction of draws in which a cross is among the N
LARGEST -- i.e. the worst crosses -- and it feeds `<trait>_post_topn` and `prob_top_tier` in the
cross-priority risk layer. `ng_posterior_cross_predict()` now takes `direction` and orients the
top-N by it. `ng_optimize_robust_mating_plan()` refuses a `posterior_topn_prob` objective whose
plan direction disagrees with the direction the column was built under, because that column (unlike
a quantile) bakes an orientation in and cannot be reinterpreted after the fact. Quantile columns
are orientation-free, so no such guard is imposed there; both directions are recorded in the
summary (`robust_direction`, `robust_posterior_direction`) for audit.

### Orientation is not always the trait's breeding direction

Care is needed here, and the new internal `ng_run_cp_value_orientation()`
(`R/39_cross_prediction_runner.R`) is what the runner uses. `ng_run_cp_trait_value()` applies no
sign to the pure-variance and parent-distance metrics -- more within-family variance is more
opportunity whichever way the trait points -- so `trait_value_metric` in
`{"pmv", "vpm", "parent_distance", "le"}` is ALWAYS `"maximize"`, even for a decrease trait. Only
`mean` and `usefulness` carry the trait's units and follow its direction. Blindly passing the
trait direction would have flipped the top-N for a minimize trait ranked on family variance.
`direction` on both public functions is documented as the orientation of the ranked value.

## Requirement 3 -- the aggregation caveat, stated

The objective is a SUM of per-cross quantiles, which is not the quantile of the plan's total.
This is now stated in the comment block where the objective is defined (`R/30`, above
`ng_optimize_robust_mating_plan()`) and recorded in the plan summary the way the
normal-approximation fallback already was:

* `robustness_quantile_aggregation = "sum_of_per_cross_quantiles"`
* `robustness_quantile_aggregation_note` -- a plain-language sentence naming the tail probability
  and saying the total does not carry that coverage.

The statement made is: quantiles are additive only when the per-cross posteriors are comonotonic
(perfectly rank-correlated). Here they are positively but imperfectly correlated, through the
shared marker-effect draws. Under the approximately jointly Gaussian posterior this layer
produces (the draws are linear functionals of a Gaussian ridge posterior), the plan total's
q-quantile is `mean_total + z * sd_total` with `sd_total <= sum(sd_i)`, so the summed objective
sits on the pessimistic side of the plan's true quantile in BOTH orientations: it is a
conservative bound, not a plan-level coverage statement. No joint plan-level quantile is
computed -- that is a larger design change and out of scope.

## Verification

New plain-script tests (the repository style; `tests/` is not testthat):

* `tests/robust_quantile_direction.R`
  * an exact non-default quantile (`q = 0.25`) is served with no `allow_normal_approximation`,
    and `robust_quantile_source` names the cached `..._post_q0_25` column;
  * COLLATERAL-DAMAGE GUARD: requesting a robustness quantile leaves
    `<gain>_post_lower/_upper/_mean`, `pmv_post_lower/_upper` and `ranked_value_post_sd`
    bit-for-bit identical (`all.equal(..., tolerance = 0)`), and `ci_level` stays 0.95;
  * the cached quantile is a genuine interior quantile (strictly inside the 95% interval, and not
    a copy of a CI tail);
  * DIRECTION: on a genuinely minimize-oriented ranked value (`mean - i*SD`, built via
    `value_fun` + `direction = "minimize"`), the plan uses tail 0.75, sources the `..._post_q0_75`
    column, and is exactly the N crosses with the SMALLEST upper-tail value; it differs from the
    maximize-direction plan on the same data and is strictly better at the pessimistic tail
    (mean upper-tail value -2.7565 vs -0.8770 for the direction-blind plan);
  * minimize-oriented `posterior_topn_prob` selects the smallest ranked values, and a mismatched
    plan direction on that objective is refused;
  * the refusal still fires for `q = 0.40` (uncached, not a CI tail) and names
    `allow_normal_approximation`; the opt-in still works and is still flagged;
  * the `NULL` default still resolves to the exact cached lower CI tail (unchanged behaviour);
  * `ng_run_cp_value_orientation()` unit checks, including the pure-variance exception.
* `tests/runner_robust_quantile.R`
  * `ci_level` and `robustness_quantile` reach `ng_posterior_cross_predict()` through
    `ng_run_cross_prediction()` and are echoed in `result$settings`;
  * a decrease trait's posterior metadata carries `direction = "minimize"`, an increase trait's
    `"maximize"`;
  * the reported interval is identical between a run with and without `robustness_quantile`, for
    both traits;
  * `ci_level = 0.80` narrows the reported interval and leaves `_post_mean` untouched;
  * end-to-end, the frontend's orchestration (run, then robust-allocate at the breeder's
    quantile) succeeds for BOTH directions without `allow_normal_approximation`;
  * invalid `ci_level` / `robustness_quantile` are rejected with the documented messages.

Existing scripts run (all PASS): `posterior_cross_predict`, `posterior_multitrait_cross_predict`,
`posterior_multitrait_rng_scoping`, `posterior_pmv_full`, `posterior_genetic_covariance`,
`user_cross_prediction_posterior_parameters`, `priority_risk_portfolio`,
`priority_risk_portfolio_multitrait`, `threshold_probability_multitrait`, `staged_pipeline`,
`runner_integration`, `run_cross_prediction_json_messages`, `crossing_plan_workbook`,
`strategy_exports`, `shiny_app`, `user_cross_prediction_workflow`, `ctx_null_fields`,
`contract_schema_drift`, plus `tests/testthat` (80 pass, 1 skip: the opt-in C++ cross-check) and
`tools/run_statistical_release_gate.R` (44/44 PASS), whose "Portfolio risk" block covers the
robust allocator's exact-lower-quantile reduction, its approximation opt-in, and the P(top-N)
objective.

No existing statistical guard was weakened. The only behaviour change to an existing output is
the `posterior_topn_prob_<N>` orientation fix for minimize traits scored on mean or usefulness,
where the previous numbers ranked the worst crosses as most stable; this is recorded in NEWS.

## Known limitation (documented, not fixed)

The plan-level aggregation caveat above. The objective remains a conservative bound on the plan
total's quantile rather than the quantile itself.
