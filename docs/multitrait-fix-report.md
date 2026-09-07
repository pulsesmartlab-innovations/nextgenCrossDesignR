# Multi-trait fix report — nextgenCrossDesign 0.26.0

Scope: four defects from the multi-trait path audit
(`~/NextGenCrossDesign/.superpowers/sdd/2026-09-03-check-reference-lines-frontend/multitrait-audit.md`,
findings 6, 3, 7 and 5), plus wiring the index posterior that finding 5 identified as unreachable.
Extends the 0.25.0 exact-quantile + direction work recorded in `docs/robust-quantile-fix-report.md`.

Backend only. No frontend file was touched (frontends are owned elsewhere and out of scope), so
the frontend defects the audit lists (E1, C2, C3, 2b/2d) are **not** addressed here — but the
backend capability E1 needs (a posterior on the index) now exists and is returned by the runner.

---

## Fix 1 — direction sign in the index posterior's usefulness value

**File:** `R/32_posterior_multitrait.R` (the `use_pmv` branch of the per-draw, per-trait loop).

**Defect.** Under `value_mode = "usefulness"` the per-draw per-trait value was
`mean + i * sigma` with no direction sign. For a minimize trait (disease, lodging) that is the
*unfavourable* tail of the family: more within-family variance was charged as a liability instead
of credited as an opportunity, and `value_z` then applied `-1` on top of it.

**Fix.**

```r
sign_j <- if (identical(traits$direction[[j]], "maximize")) 1 else -1
trait_cross_value[, j] <- mp_js + sign_j * intensity * sigma_js
```

Matches `ng_run_cp_trait_value()` (`R/39_cross_prediction_runner.R`) exactly: the sign multiplies
**only** the `i * sigma` term, never the mean. `traits$direction` is already normalised to
`"maximize"` / `"minimize"` by `ng_multitrait_spec()` (`R/19:172`), which `R/32` calls on entry.

**Evidence.** `tests/posterior_multitrait_usefulness_direction.R`. With one trait and one draw,
`multi_trait_score` is exactly the rank-normalisation of that draw's usefulness value, so the
function's output identifies the expression it used. The test rebuilds that single draw outside
the package and asserts:

* `multi_trait_score_post_mean == rank_normalize(mean - i*sigma, bigger_is_better = FALSE)` for a
  minimize trait (tolerance 1e-10), and is **not** equal to the `mean + i*sigma` version;
* flipping the direction moves the value by exactly `2 * i * sigma` — i.e. the mean is not signed;
* a maximize trait is unchanged (`sign = +1`).

Measured on the test data: `i*sigma` range `[5.714, 13.502]`, so the sign has real magnitude to
bite on; `|fixed - broken|` range `[11.428, 27.004]`.

This defect was latent in 0.25.0 (no callers). Fix 4 makes it live, which is why it went first.

---

## Fix 2 — threshold penalty scaled by the wrong column's dispersion

**File:** `R/19_multi_trait_selection.R`, the per-trait threshold block.

**Defect.** The DEFECT-1 fix had correctly moved the *comparison* onto `threshold_column`
(`<trait>_mean`, trait units) but left the *normalising scale* on `column` — whatever
`trait_value_metric` produced. Under `trait_value_metric` in `{pmv, vpm, var_complex}` the ranking
column is a **variance**, so a trait-unit deficit was divided by the IQR of a variance.

**Fix.** Scale the deficit by the dispersion of the column the deficit is measured on:

```r
thresh_scale <- if (identical(traits$threshold_column[[i]], traits$column[[i]])) scale else
  ng_multitrait_value_scale(x_thresh)
```

Division by zero is impossible: `ng_multitrait_value_scale()` (`R/19:333`) already falls back
IQR → MAD → SD → range → 1 and can never return 0 or a non-finite value. The explicit
`identical()` short-circuit makes the no-op case exact rather than merely equal, so every direct
caller whose `threshold_column` IS `column` (the `ng_multitrait_spec()` default) is bit-identical
to 0.25.0.

**Evidence.** `tests/threshold_penalty_scale_invariance.R`. Same threshold column, same thresholds,
only the ranking metric changed:

| ranking column | 0.25.0 mean violation | 0.26.0 mean violation |
|---|---|---|
| `protein_usefulness` (trait units) | 0.466177 | 0.493942 |
| `protein_pmv` (variance units) | 141.126154 | 0.493942 |
| **swing from an unrelated knob** | **302.7x** | **1.000x (exact)** |

(The audit measured 52x on its own data; this test's variance column is on a smaller scale, so the
same defect shows up larger. The 0.25.0 usefulness figure also differs slightly from the 0.26.0 one
because `IQR(protein_usefulness) != IQR(protein_mean)` — that difference is exactly the defect.)

The test also asserts the closed form
`violation == pmax(min_value - x_thresh, 0) / (IQR(x_thresh)/1.349)`, a non-regression case where
`threshold_column == column`, and a degenerate (zero-spread) threshold column that must stay finite.

**Not claimed:** the *applied* penalty `effective_penalty * violation` is still metric-dependent,
because `threshold_penalty_autoscale` keys off `IQR(weighted_score)` and the index's own scale
legitimately changes with the ranking metric. The invariant fixed here is the one that was wrong:
the deficit is now dimensionless in its own units.

---

## Fix 3 — multi-trait branch never passed top-N / posterior-SD

**Files:** `R/39_cross_prediction_runner.R` (`ng_cp__stage_rank`, multi-trait branch) and
`R/37_cross_portfolio_risk.R` (`ng_annotate_cross_priority_multitrait()`).

**Defect.** The single-trait branch passes `post_sd` and `prob_top_tier`; the multi-trait branch
passed neither, and `ng_annotate_cross_priority_multitrait()` did not even accept them. So
`prob_top_tier` was **absent** from a multi-trait `candidate_crosses`, and `confidence_method`
always fell back to `"midparent_pev_index"` even with posterior prediction on.

**Fix.** `ng_annotate_cross_priority_multitrait()` gains `post_sd` and `prob_top_tier` and mirrors
the single-trait contract exactly: posterior SD (when present) drives `ng_cross_confidence()` under
`method_prefix = "posterior_ci"`, otherwise the mid-parent-PEV index path is used unchanged; and
`prob_top_tier` is always emitted (`NA` when no index posterior was run), never binned into
`risk_bin` and never relabelled as confidence. The runner reads
`multi_trait_score_post_sd` and `multi_trait_score_post_topn` off the cross table and passes them.

**Evidence.** `tests/runner_multitrait_index_posterior.R`, real output:

```
multi-trait prob_top_tier: n_finite = 45, range [0.020, 0.213]
confidence_method: single-trait = posterior_ci, multi-trait = posterior_ci,
                   multi-trait posterior off = midparent_pev_index_partial
```

and `priority_risk_diagnostics$posterior_used == TRUE` on the multi-trait run. With posterior off,
`prob_top_tier` is present and all-`NA` and the PEV fallback is retained.

---

## Fix 4 — wiring the index posterior

**Files:** `R/32_posterior_multitrait.R`, `R/39_cross_prediction_runner.R`.

**Defect.** `ng_posterior_multitrait_cross_predict()` was exported, documented, tested and
release-gated with **zero non-test callers**. `multi_trait_score` therefore reached the user as a
point estimate with no uncertainty, and anything needing an uncertainty on the plan's merit had to
reach for `posterior_predictions[[1]]` — one trait, selected by direction-file row order.

### 4a. Called from the runner

`ng_cp__stage_predict()` now calls it after the per-trait loop, gated on
`run_posterior_prediction && nrow(trait_spec) > 1 && nrow(cross_table) > 0`, with
`value_mode = "mean"` (cost `O(n*m*T*S)` — the same order as the per-trait posteriors already
computed, not the per-draw PMV rescoring `"usefulness"` needs). Result shape is **additive**:

* new top-level `result$posterior_multitrait` (the full table plus both metadata attributes);
* `multi_trait_score_post_mean/_lower/_upper/_post_sd`, `multitrait_posterior_topn_prob_<N>` and
  any `multi_trait_score_post_q<prob>` columns ride `candidate_crosses`, so they survive candidate
  filtering and allocation by row alignment, exactly like `<trait>_post_sd` / `<trait>_post_topn`;
* nothing existing is removed or renamed; `result$posterior_predictions` keeps its trait names and
  its order.

The fit is given the same marker-effect training augmentation the per-trait posteriors use
(parents + training rows, with only parents ever crossed), so the index posterior is fitted on the
same information rather than on the parent panel alone.

### 4b. `robustness_quantile`, in the shape 0.25.0 gave `R/30`

`ng_posterior_multitrait_cross_predict()` gains `robustness_quantile`. Both `q` and `1 - q` are
cached as extra **empirical** quantiles of the same `index_mat` draws that produced the credible
interval, named by the existing `ng_posterior_quantile_col()` convention
(`multi_trait_score_post_q0_25`). `ci_level` is untouched. A second, narrower `"posterior"`
attribute is attached in the shape `ng_optimize_robust_mating_plan()` reads, so a robust allocation
can now run on the index. Verified end to end in the test:

```
robust plan on the INDEX: quantile source = multi_trait_score_post_q0_25
                          (exact, not a normal approximation)
```

`ng_posterior_cross_predict()`'s `ranked_value_post_sd` also gets its index counterpart,
`multi_trait_score_post_sd`, which is what Fix 3 consumes.

### 4c. `direction` for the index: **"maximize"**, hard-coded, verified not assumed

`multi_trait_score` is direction-normalised higher-is-better for **every** index method, so a
`direction` argument here could only ever take one value. It is hard-coded with a comment rather
than plumbed. Three independent lines of evidence:

1. **Stated design invariant in the code being called.** `ng_add_multitrait_score()`
   (`R/19_multi_trait_selection.R`) carries a "Design invariant (deliberate; do not normalize this
   away)" comment: *"the emitted index is oriented higher = better for EVERY method... guaranteed
   upstream of the combination"*.
2. **Both combination paths, read directly.** Rank path: `ng_rank_normalize(x, bigger_is_better =
   identical(direction, "maximize"))` negates a minimize trait before ranking, then combines with
   weights that `ng_multitrait_resolve_weights()` forces `>= 0` and normalises to sum 1. Solved
   path: `oriented <- if (maximize) x else -x`, `value_z <- (oriented - center)/scale`, combined
   with `b = P^-1 G a` / `b = G^-1 d` whose `a, d >= 0` are enforced at `R/19:181-188`. In both
   cases the direction is applied *upstream* of the combination, so the output points up.
3. **Measured on a real 2-trait run** (yield `increase`, disease `decrease`), asserted in
   `tests/runner_multitrait_index_posterior.R`:
   `cor(multi_trait_score, yield_mean) = +0.988`, `cor(multi_trait_score, disease_mean) = -0.988`.

This matches the audit's own conclusion at B1 ("It does **not** need `direction`").

### 4d. The re-standardisation caveat is documented, not "fixed"

`ng_add_multitrait_score()` is called **inside** the draw loop — that is what makes each column of
`index_mat` a coherent joint draw rather than a recombination of per-trait marginal intervals — and
it re-derives its own rank-normalisation / IQR centring from the rows it is handed. So the index is
re-standardised within every draw, and any component of posterior uncertainty that shifts or
rescales the whole candidate pool together is removed before the quantile is taken. Consequences:

* `multi_trait_score_post_lower/_upper` are an interval on a **per-draw relative** index, narrower
  than a genuine index credible interval;
* they are **not** on the same scale as the point-estimate `multi_trait_score` and must not be
  differenced against it;
* rank stability (`multitrait_posterior_topn_prob_*`) and conservative-tail **ordering** are
  unaffected.

This is documented in a block comment where the columns are produced, and recorded in the metadata
as `index_rescaling = "per_draw_restandardized"` plus a prose `index_rescaling_note`, so the
frontend can badge it. Removing the re-standardisation is a design change to the index itself and
was deliberately left out of scope.

### 4e. Scope notes recorded rather than silently absorbed

* `R/32` refits its own per-trait ridge posteriors (same model, same training rows, different
  seed). Those are draws from the same posterior — not a different one — but they are not the
  *same* draws as `posterior_predictions[[trait]]`.
* `R/32` scores the per-draw index with `ng_add_multitrait_score()`'s **default** threshold penalty
  weight, not the run's `threshold_penalty_weight`. It is a posterior on the index, not on the
  run's threshold-penalised objective.
* The index posterior is built in `ng_cp__stage_predict`, i.e. **before** `ng_cp__stage_index`.
  For a multi-trait run using `economic_index` / `desired_gain` with posterior prediction on and
  no `genetic_covariance` supplied, the hard error now comes from `R/32`
  ("posterior economic/desired-gain prediction requires genetic_covariance or matched
  genetic_covariance_draws") instead of from `ng_multitrait_index_covariance()` a stage later.
  No previously-working run breaks: `ng_add_multitrait_score()` already refused that
  configuration (`R/19:340-345`), so such a run failed in 0.25.0 too — only the message moves.
* `ng_optimize_robust_mating_plan(objective = "posterior_topn_prob")` looks for
  `posterior_topn_prob_<N>`; the index table names its columns
  `multitrait_posterior_topn_prob_<N>`, so only the `"posterior_quantile"` objective is served
  directly today.

---

## Single-trait results are bit-identical

Asserted, not assumed. A baseline was produced by running an identical script against a pristine
`git archive HEAD` tree (0.25.0) and compared with `all.equal(..., tolerance = 0)`:

| run | result |
|---|---|
| single trait, posterior ON (+ `robustness_quantile = 0.25`) | **BIT-IDENTICAL** |
| single trait, posterior OFF | **BIT-IDENTICAL** |
| multi trait, posterior OFF | only change is the added all-`NA` `prob_top_tier` column; every shared column, `plan_summary`, `trait_scores`, `posterior_predictions` and `priority_risk_diagnostics` **bit-identical** |

Compared objects: `candidate_crosses`, `selected_crosses`, `plan_summary`, `trait_scores`,
`posterior_predictions`, `priority_risk_diagnostics`.

---

## Test output

New tests:

```
tests/posterior_multitrait_usefulness_direction.R
  i*sigma range [5.71413, 13.50218]; |fixed - broken| range [11.42825, 27.00435]
  posterior_multitrait_usefulness_direction: OK

tests/threshold_penalty_scale_invariance.R
  mean violation, ranking column = protein_usefulness : 0.493942
  mean violation, ranking column = protein_pmv        : 0.493942
  threshold_penalty_scale_invariance: OK

tests/runner_multitrait_index_posterior.R
  multi-trait prob_top_tier: n_finite = 45, range [0.020, 0.213]
  confidence_method: single-trait = posterior_ci, multi-trait = posterior_ci,
                     multi-trait posterior off = midparent_pev_index_partial
  cor(multi_trait_score, yield_mean) = +0.988 ; cor(multi_trait_score, disease_mean) = -0.988
  robust plan on the INDEX: quantile source = multi_trait_score_post_q0_25
                            (exact, not a normal approximation)
  runner_multitrait_index_posterior: OK
```

Existing tests re-run (all PASS): `multi_trait_selection`, `multitrait_index_math`,
`multitrait_exact_cov_integration`, `multitrait_economic_index_external_cov`,
`threshold_column_units`, `threshold_alias_columns`, `threshold_probability_multitrait`,
`posterior_cross_predict`, `posterior_multitrait_cross_predict`,
`posterior_multitrait_rng_scoping`, `posterior_pmv_full`, `priority_risk_portfolio`,
`priority_risk_portfolio_multitrait`, `cross_priority`, `cross_priority_workbook`,
`cross_priority_plots`, `robust_quantile_direction`, `runner_robust_quantile`,
`runner_integration`, `runner_training_set`, `runner_cross_cost`, `check_reference_runner`,
`check_reference_invariant`, `multi_trait_validation_smoke`, `multi_trait_crop_validation`,
`user_cross_prediction_posterior_parameters`, `statistical_invariants`.

Statistical release gate, against the installed 0.26.0:

```
Statistical code gate: PASS (44/44)
Wrote: docs/STATISTICAL_RELEASE_GATE_RESULTS.csv
Wrote: docs/STATISTICAL_RELEASE_GATE.md
```

Install check: `packageVersion("nextgenCrossDesign")` reports **0.26.0**.

---

## Where the brief's analysis differed from the code

* **Fix 2, "handle a zero or non-finite IQR without dividing by it."** Already handled:
  `ng_multitrait_value_scale()` falls back IQR → MAD → SD → range → 1 and never returns 0 or a
  non-finite value. No new guard was needed; a test covers the degenerate case anyway.
* **Fix 3, "pass them, matching the single-trait branch."** The arguments could not simply be
  passed: `ng_annotate_cross_priority_multitrait()` (`R/37`) did not accept `post_sd` or
  `prob_top_tier` at all, so the parameters had to be added there first. And the values the
  multi-trait branch needs did not exist before Fix 4 — there was no index-level posterior SD or
  index top-N probability anywhere in the run. Fix 3 is therefore genuinely downstream of Fix 4,
  not independent of it.
* **Fix 4, `direction`.** The brief asked to verify rather than assume. Verified: the index is
  always `"maximize"`, hard-coded with a comment; no argument plumbed. (The audit reached the same
  conclusion at B1.)
* Everything else in the brief matched the code as written.

## Known limitations

* The frontend defects the audit ranks highest (E1's trait-1 robust plan, C2's PMV-as-per-progeny
  variance, C3's direction vocabulary, 2b's silent independence fallback) are untouched — that tree
  is out of scope here. The backend now *offers* the index posterior E1 needs; consuming it is a
  frontend change.
* `multi_trait_score_post_*` on a multi-trait run adds one extra set of per-trait ridge posterior
  fits plus `S` index recomputations. Cheap relative to the per-trait posterior loop already
  running, but not free.
