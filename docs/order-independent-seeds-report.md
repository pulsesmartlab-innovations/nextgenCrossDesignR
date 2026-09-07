# Order-independent per-trait seeds (0.28.0)

A trait's scientific result depended on **where the breeder happened to list it** in the
trait-direction file. This report records the reproduction, the fix, why the two seed sites were
treated differently, what else was audited, and what was deliberately left alone.

Package root: `cross_prediction/` (the tracked package `nextgenCrossDesign` lives at the repo
root). Baseline: `52c63a5` (0.27.0).

---

## 1. The defect

`R/39_cross_prediction_runner.R`, inside `ng_cp__stage_predict()`'s per-trait job, where `i` is
the trait's **row position** in the direction file:

| line (0.27.0) | expression | consumer |
| --- | --- | --- |
| 970 | `seed = seed + i - 1L` | `ng_fit_ridge_effects()` → `ng_choose_ridge_lambda()` → `set.seed(seed); folds <- sample(rep(seq_len(kfold), length.out = n))` (`R/02_effects.R:135-136`), and `ng_ridge_cv_predict()` (`R/02_effects.R:166-167`) |
| 1008 | `seed = seed + 1000L + i - 1L` | `ng_fit_ridge_effects_posterior()` — the per-trait posterior draw stream |

A different seed gives a different CV fold partition → possibly a different selected ridge
lambda → different fitted marker effects → **a different progeny variance**, and from there a
different usefulness, index contribution, mid-parent PEV, robust allocation, P(beat check) and
crossing plan.

`R/32_posterior_multitrait.R:188` had the same defect one level down: `seed = seed + j` over the
rows of the caller's `traits` table. `ng_posterior_multitrait_cross_predict()` is called by the
runner (`R/39:1261`) to build the **index** posterior, so the file order moved that too.

### The mechanism is not hypothetical

Unit-level probe of `ng_choose_ridge_lambda()` on the user's own demo panel
(`test_data/rich_geno.csv`, 160 lines × 3 191 markers), lambda grid `10^seq(-2, 5, length.out = 20)`,
`kfold = 5`, seeds 1…6:

```
demo:YIELD     lambda(seed 1..6) = 7848 7848 7848 7848 7848 7848    var(beta) max/min = 1.00
demo:FOL_DIS   lambda(seed 1..6) = 3360 3360 7848 3360 3360 3360    var(beta) max/min = 2.80
demo:Extract   lambda(seed 1..6) = 3360 7848 1438 3360 7848 3360    var(beta) max/min = 6.44
demo:STM_BRK   lambda(seed 1..6) =  264  616  616  616  616 1438    var(beta) max/min = 3.21
synth n=40 m=60 (seed 5)  lambda = 264 1438 264 616 616 1e5         var(beta) max/min = 7.8e+04
synth n=40 m=60 (seed 6)  lambda = 616 616 1e5 1e5 113 1438         var(beta) max/min = 2.6e+05
```

So a one-step change in the seed routinely moves the selected lambda by one or more grid steps,
and the resulting marker-effect variance by up to five orders of magnitude on small panels.

---

## 2. Reproduction (my numbers)

### 2.1 On the user's demo data — the order effect is real but hits only `cv_predictive_r2`

`test_data/*` (160 lines, 3 191 markers), two traits `YIELD` (increase) and `FOL_DIS` (decrease),
`seed = 1`, defaults otherwise. Genotype/map aligned to each other and the 9 parents above the
2 % inbred het tolerance dropped, leaving 151 parents. Running the direction file as
`(YIELD, FOL_DIS)` and then as `(FOL_DIS, YIELD)`:

```
YIELD_cv_predictive_r2     (Y,D) 0.10390457   (D,Y) 0.089654211   identical = FALSE
FOL_DIS_cv_predictive_r2   (Y,D) 0.1151876    (D,Y) 0.12066972    identical = FALSE
FOL_DIS_vpm / _pmv / _value / _midparent_pev ......................  identical = TRUE
YIELD_vpm  / _pmv / _value / _midparent_pev  ......................  identical = TRUE
```

At n = 151 both fold splits happened to select the same lambda for both traits, so on this
particular panel only the reported CV predictive R² (itself a per-trait output, and the input to
`min_effect_reliability`) flipped. **The variance columns did not move here.**

### 2.2 On a panel where lambda genuinely flips — seven orders of magnitude

24 parents, 40 markers over 2 chromosomes, data seed 346, run `seed = 909`, `progeny = "DH"`,
`trait_value_metric = "usefulness"`, `uc_variance_source = "pmv"`, `optimizer = "greedy_local"`,
`n_crosses = 5`. Traits `yield` (increase) and `disease` (decrease). Mean over the 276 candidate
crosses, comparing **the same trait** in the two file orders:

**BEFORE (0.27.0), posterior off**

```
column                    disease listed 2nd   disease listed 1st   ratio        identical
disease_vpm               3.197428e-07         2.5124748            1.27e-07     FALSE
disease_pmv               0.00092880902        2.5403264            3.66e-04     FALSE
disease_var_complex       0.00092880902        2.5403264            3.66e-04     FALSE
disease_value             -1.686178            -3.7593466           0.449        FALSE
disease_cv_predictive_r2  -0.1045129           0.056829355         -1.839        FALSE
disease_midparent_pev     0.000850763          1.5039769e-05       56.6          FALSE
yield_vpm                 0.0030084744 (1st)   0.05773146 (2nd)     0.052        FALSE
disease_mean              -1.643669            -1.643669            1.00         TRUE

selected plan (yield, disease):  P10|P11 P11|P19 P08|P11 P10|P19 P08|P10
selected plan (disease, yield):  P02|P11 P10|P11 P02|P10 P11|P19 P10|P19
```

`disease_vpm` moves by a factor of **7.9 × 10⁶** — the same trait, the same data, the same run
seed, one row swapped — and **three of the five selected crosses change**. With the posterior on,
26 of the 35 numeric columns of `posterior_predictions$disease` differed.

The trait sitting in **row 1** always reproduced the single-trait run exactly; row 2 was the one
that drifted.

**AFTER (0.28.0), same workload**

```
disease_vpm               2.5124748            2.5124748            1.0          TRUE
disease_pmv               2.5403264            2.5403264            1.0          TRUE
disease_var_complex       2.5403264            2.5403264            1.0          TRUE
disease_value             -3.7593466           -3.7593466           1.0          TRUE
disease_cv_predictive_r2  0.056829355          0.056829355          1.0          TRUE
disease_midparent_pev     1.5039769e-05        1.5039769e-05        1.0          TRUE
yield_vpm                 0.0030084744         0.0030084744         1.0          TRUE

selected plan (yield, disease):  P02|P11 P10|P11 P02|P10 P11|P19 P10|P19
selected plan (disease, yield):  P02|P11 P10|P11 P02|P10 P11|P19 P10|P19

posterior on: posterior_predictions$disease  35/35 numeric columns identical
              disease_post_sd  0.014378521 == 0.014378521
single-trait baseline:  disease alone == disease in either multi-trait order  (all columns)
```

### 2.3 Note on the brief's figure

The brief quoted `disease_vpm` as "2.1e-09 when disease is listed first and 0.147 when listed
second". I could not reproduce that pairing: on the shipped demo data the variance columns are
order-*invariant* (§2.1), and 2.1e-09 / 0.147 look like values from two different panels rather
than two orders of one panel (my own small synthetic panel gives `disease_vpm` ≈ 2.26e-09 in
*both* orders, and the demo panel gives `FOL_DIS_pmv` ≈ 0.124 in *both*). The **defect** is
exactly as described and its magnitude can exceed the brief's claim (§2.2 measures 7.9e+06×);
only the specific pair of numbers is not reproducible as stated.

---

## 3. The fix

### 3.1 New helper — `R/00_utils.R`

```r
ng_name_hash32(x)                       # UTF-8 byte polynomial: h <- (h * 131 + byte) %% (2^31 - 1)
ng_trait_rng_seed(base_seed, trait, salt = 0L)   # (base + salt + hash) %% (2^31 - 1), as integer
```

Written out in plain R arithmetic on `charToRaw(enc2utf8(s))` **specifically so it depends on no
hashing internal** that R, a platform or a locale is free to change. Every intermediate is at most
`131 * (2^31 - 2) + 255 ≈ 2.8e11`, exactly representable in a double, so the value is bit-identical
everywhere. The return value is in `[0, 2^31 - 2]`, always safe for `set.seed()` (checked in the
regression test with a 1 350-character trait name and a base seed of 2 147 483 000).

### 3.2 Line 970 — lambda CV: **share one partition** (`seed = seed`)

Adopted, and the brief's reasoning holds against the real code:

* The seed's only reach is `ng_choose_ridge_lambda()` and `ng_ridge_cv_predict()`
  (`R/02_effects.R:131-175`), i.e. it *only* selects the fold partition. It is a nuisance
  parameter of lambda selection.
* It cannot couple the traits. Once lambda is fixed, `beta_j = X'(XX' + lambda I)^{-1} y_c,j`
  (`R/02_effects.R:39-41`) is a **deterministic** function of that trait's own `y`. There is no
  random innovation to share, so a common partition induces no correlation between traits'
  fitted effects. (This is exactly the property that does *not* hold at line 1008.)
* It makes the per-trait CV comparisons paired: differences in `cv_predictive_r2` between traits
  now reflect the traits, not the split.
* It preserves bit-identity for a single-trait run *and* for the first-listed trait of any run,
  because position 1 already received `seed + 1 - 1 = seed`.

### 3.3 Line 1008 — posterior draws: **distinct, name-derived** streams

`seed = ng_trait_rng_seed(seed, paste(trait, column, sep = "\r"), salt = 1000L)`.

Distinctness is mandatory here: this seed drives the actual posterior innovations, and
`ng_posterior_multitrait_cross_predict()` combines per-trait draws into the index posterior. Equal
seeds would give the traits identical innovations and manufacture cross-trait correlation in
exactly that quantity. The key is the trait name **and** its phenotype column, because nothing
upstream (`ng_run_cp_trait_spec()`, `R/39:566`) forces trait names to be unique.

### 3.4 `R/32_posterior_multitrait.R:188` — same treatment

`seed = ng_trait_rng_seed(seed, trait_names[[j]])`. Distinct per trait, keyed on the name.
(`R/39:1274`'s `seed = seed + 2000L` for the whole index-posterior call is **not** a defect: it
is one call, and `2000L` is a constant salt, not a position.)

---

## 4. Single-trait bit-identity

Verified against a pristine `git archive HEAD` tree at `/tmp/ngcdseed` (0.27.0), running the same
workload in both trees and comparing every numeric column of the candidate table with
`identical()`:

| workload | result |
| --- | --- |
| single trait, posterior off | **BIT-IDENTICAL** — 26/26 numeric columns |
| single trait, selected plan | **BIT-IDENTICAL** |
| single trait, posterior on | 23/28 identical; the 5 that differ are all posterior-derived: `disease_post_sd`, `disease_post_topn`, `cross_confidence`, `relative_precision`, `prob_top_tier` |
| multi trait, posterior off | trait 1 (`yield_*`) identical; trait 2 (`disease_*`) and everything downstream of it changed — **this is the fix** |
| multi trait, posterior on | as above, plus the index-posterior columns (R/32) |

**The single-trait posterior draw stream does change, deliberately.** It was `seed + 1000L`; it is
now `seed + 1000L + hash(trait)`. Justification: the alternative — special-casing a one-trait run
back onto the old stream — would make a trait's draws depend on **how many other traits share the
file**, reintroducing precisely the context dependence this change removes. Nothing statistical
changes: it is a different sample from the *same* posterior (Monte Carlo noise), and every
deterministic single-trait output — variances, usefulness, index, mid-parent PEV, plan — stays
bit-identical.

---

## 5. Other position-derived randomness — what was audited

Every RNG-consuming call in the 56 tracked `R/*.R` files was enumerated
(`set.seed`, `sample`, `rnorm`, `runif`, `rbinom`, `rmultinom`, `rWishart`, `mvtnorm::`) and every
`seed = <expr>` containing a loop index.

| site | verdict |
| --- | --- |
| `R/39:970`, `R/39:1008` | **FIXED** (the reported defect) |
| `R/32:188` `seed + j` | **FIXED** — same defect, and it is on the runner's production path |
| `R/39:1274` `seed + 2000L` (index posterior) | OK — constant salt on a single call, no index |
| `R/39:936` `set.seed(seed)` before the trait loop | OK — each job re-seeds inside `ng_fit_ridge_effects*`, so the ambient stream is not what the fits consume |
| `R/02:135`, `R/02:166` fold sampling | OK — this is the consumer that was being mis-seeded, not itself position-aware |
| `R/30_posterior_prediction.R:84/147` `rnorm` in the draw loop | OK — the whole sampler runs under `ng_with_rng_seed(seed, …)` with the caller's (now identity-derived) seed; the loop index is a *draw* counter |
| `R/33:189-195` `ng_fixed_seed_normal_matrix()` | OK-with-note. Fixed constant seed, no loop index. It fills an `n_draws × t_traits` z-matrix column-wise, so permuting traits reassigns which block of the antithetic stream each trait gets. This is a Monte-Carlo quadrature of a *joint* multivariate-normal probability, not a per-trait model fit: permutation perturbs only the MC error term (antithetic, fixed seed), never a fitted quantity. Left alone. |
| `R/09_family_calibration.R:79/104` | OK — one seed for one family-sampling call, no per-trait index |
| `R/40_evolutionary_optimizer.R` | OK — driven by `evol_seed`, indices are population/generation counters |
| `R/45_staged_runner.R` | OK — contains no seed at all; it delegates to the runner |
| `R/46`, `R/48`, `R/50_polyploid_staged_runner.R` | OK — a single shared `ridge_seed` regardless of trait, no position term |
| `R/34`, `R/36`, `R/37`, `R/19`, `R/49`, `R/50_ucpc.R` | OK — no RNG |

### `R/31_genetic_covariance.R` — FOUND, deliberately NOT fixed

Four sites carry the identical defect:

* `:108` `ng_fit_ridge_effects(..., seed = seed + j)` in `ng_genetic_cov_two_stage_ridge()`
* `:623` `ng_fit_ridge_effects_posterior(..., seed = seed + j)` (`beta_posterior`)
* `:655` `ng_fit_ridge_effects(..., seed = seed + j)` (`parametric_bootstrap`)
* `:665-668` residual draws taken from **one** `rnorm` stream consumed in trait order inside the
  `b`/`j` loops — the same position dependence one level down.

I implemented the fix, measured it, and then reverted it. Reasons, with evidence:

1. **Not on the runner's path.** `ng_estimate_genetic_covariance()` and
   `ng_posterior_genetic_covariance()` have **no non-test caller** inside `R/`; the runner requires
   user-supplied P and G. Nothing in a `ng_run_cross_prediction()` run reaches these seeds.
2. **They are opt-in heuristics that already say so.** `ng_posterior_genetic_covariance()` refuses
   to run without `allow_heuristic = TRUE`, and both warn they are "not a multivariate
   variance-component estimator".
3. **Their accuracy guard is a coin flip, so any seed change flips it.** Sweeping the base seed
   over 120 values on the exact simulation in `tests/genetic_covariance_estimator.R`
   (n = 100, m = 200, 3 traits), Frobenius ratio `||G_hat − G_true||_F / ||G_true||_F`:

   | scheme | mean | median | sd | P(ratio > 0.5) |
   | --- | --- | --- | --- | --- |
   | `seed + j` (0.27.0, shipped) | 0.625 | 0.408 | 0.346 | **0.47** |
   | shared `seed` | 0.624 | 0.402 | 0.347 | 0.46 |
   | `hash(name)` | 0.506 | 0.352 | 0.315 | 0.28 |

   The test asserts `ratio <= 0.5` at the single base seed 2026, where the shipped scheme lands on
   0.374 — i.e. **it passes by luck; ~47 % of base seeds would fail it on unmodified 0.27.0.**
   Both order-independent schemes fail at 2026 (shared 1.112, hash 0.943) despite being no worse —
   indeed `hash(name)` is *better* on average (Wilcoxon vs shared, p = 0.015).

   So de-positioning R/31's seeds cannot be shipped without also re-tuning that acceptance
   criterion, and re-tuning a science gate as a drive-by inside a reproducibility fix is the wrong
   trade. Recommended follow-up, as its own change: replace the single-seed Frobenius assertion
   with an ensemble statistic that passes on 0.27.0 *and* on the fixed code (so it is a genuine
   invariant, not a rescue), then switch all four sites to `ng_trait_rng_seed(seed, trait_names[[j]])`
   — `hash(name)`, not the shared seed, because R/31 correlates the traits' `Beta` columns
   (`cor(Beta)`) and independent per-trait splits are the structure the code already intends.

---

## 6. Verification

### New regression test — `tests/trait_order_invariance.R`

Permutes the direction file and asserts, with `all.equal(..., tolerance = 0)`:

* every numeric `<trait>_*` column of `candidate_crosses`, posterior off and on;
* the full per-trait `posterior_predictions[[trait]]` table (35 numeric columns each);
* the index posterior `posterior_multitrait` (8 numeric columns);
* `multi_trait_score` and the selected plan;
* that a single-trait run reproduces the multi-trait deterministic columns;
* the hash helper's literal values, uniqueness, salting, overflow safety and locale independence.

```
$ Rscript tests/trait_order_invariance.R
posterior off:
  candidate_crosses         13 disease columns identical (tolerance = 0)
  candidate_crosses         13 yield   columns identical (tolerance = 0)
  disease_vpm mean: yield-first = 2.512474803, disease-first = 2.512474803
posterior on:
  candidate_crosses         15 disease columns identical (tolerance = 0)
  candidate_crosses         15 yield   columns identical (tolerance = 0)
  posterior_predictions     35 yield   columns identical (tolerance = 0)
  posterior_predictions     35 disease columns identical (tolerance = 0)
  posterior_multitrait       8 index-posterior columns identical (tolerance = 0)
  single-trait run reproduces multi-trait disease  12 columns (tolerance = 0)
  single-trait run reproduces multi-trait yield    12 columns (tolerance = 0)
  ng_name_hash32: yield = 1512549153, disease = 794204123 (stable across sessions)
trait_order_invariance: OK
```

The same file run against the pristine 0.27.0 tree fails, as it must:

```
$ cd /tmp/ngcdseed && Rscript tests/trait_order_invariance.R
Error: candidate_crosses: disease columns depend on the trait's row position:
  disease_value, disease_pmv, disease_pmv_fast, disease_pmv_used, disease_vpm,
  disease_mean_gebv, disease_var_complex, disease_cv_predictive_r2, disease_midparent_pev
```

### Existing suites

Fast gate, `tests/testthat/` against the **installed** 0.28.0
(`/tmp/ngcd_lib_028`): `TOTAL failed: 0  errors: 0  passed: 80`
(one skip: the C++ cross-check, which needs `NG_TEST_CPP=1`).

Deep harness, `tests/*.R`, source tree — all PASS:

```
runner_integration                       runner_multitrait_index_posterior
posterior_multitrait_cross_predict       posterior_multitrait_rng_scoping
posterior_cross_predict                  posterior_pmv_full
posterior_genetic_covariance             genetic_covariance_estimator
staged_pipeline                          polyploid_staged_equivalence
multi_trait_selection                    multi_trait_validation_smoke
multi_trait_crop_validation              multitrait_exact_cov_integration
multitrait_economic_index_external_cov   multitrait_index_math
priority_risk_portfolio                  priority_risk_portfolio_multitrait
runner_robust_quantile                   runner_training_set
runner_cross_cost                        runner_marker_ld_exposure
user_cross_prediction_workflow           user_cross_prediction_matching_contract
user_cross_prediction_posterior_parameters
user_cross_prediction_var_complex_alphamate
user_supplied_pg_contract                check_reference_runner
check_reference_outputs                  threshold_probability_multitrait
smoke_test                               statistical_invariants
marker_effect_estimation                 ridge_beta_var
vignette_chronological_examples          poly4x_runner_smoke
trait_order_invariance (new)
```

Selected timings: `trait_order_invariance` 32 s, `polyploid_staged_equivalence` 21 s,
`ridge_beta_var` 20 s, `smoke_test` 21 s.

`tests/pipeline_integration.R` was NOT run to completion: it makes ~10 full
`ng_design_crosses()` calls on a 40 x 300 panel and exceeds a 10-minute budget on this machine,
independently of this change. It is provably unaffected — the only behavioural edits are inside
`ng_cp__stage_predict()` (R/39) and `ng_posterior_multitrait_cross_predict()` (R/32), and
`ng_design_crosses()` (R/05) reaches neither: `ng_cp__stage_predict` has exactly one reference,
the `ng_cp_pipeline` table at `R/39:2111` used by `ng_run_cross_prediction()`, and R/32's entry
point has exactly one caller, `R/39:1294`. The R/00_utils.R edit is purely additive
(two new functions, nothing existing touched).

### Statistical release gate

`R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_lib_028 .` from a clean `/tmp`
copy of the working tree, then
`NGCD_RELEASE_LIB=/tmp/ngcd_lib_028 Rscript tools/run_statistical_release_gate.R`:

```
Statistical code gate: PASS (44/44)
Package: `nextgenCrossDesign 0.28.0`
```

`docs/STATISTICAL_RELEASE_GATE.md` / `.csv` were regenerated against 0.28.0. The only numeric
movement in the whole table is one identity check's rounding noise
(`Fold-local centering translation invariance`, 3.553e-15 -> 1.776e-15, acceptance <= 1e-10).

`packageVersion("nextgenCrossDesign")` from the installed library: **0.28.0**.
