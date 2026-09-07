# Making user-supplied P and G safe to accept from a GUI (0.27.0)

Backend: `nextgenCrossDesign`, branch `main`, version bumped 0.26.0 -> 0.27.0.

Context: the Shiny frontend is about to let breeders supply the phenotypic (P) and
additive-genetic (G) covariance matrices the Smith-Hazel (`economic_index`) and
Pesek-Baker (`desired_gain`) selection indices require. Those matrices arrive over a
JSON bridge. Three changes harden the backend against the two ways a GUI-supplied
matrix goes wrong (wrong order/labels, and a caller who only has one of the two
matrices), plus one documentation/guard gap on the estimator that produces G.

---

## Where the brief matched the code, and where it did not

The brief was right about the transform itself and about `desired_gain`. Two points
need correcting, and one is more serious than the brief assumed.

**1. `ng_multitrait_cov_to_value_z()` was indeed purely positional** -- confirmed at
`R/19_multi_trait_selection.R`. It applied `signs` and `value_scales` by row/column
order and never looked at `dimnames(M)`.

**2. But the main entry point already had a partial by-name reorder.**
`ng_add_multitrait_score()` contained a local `ng_multitrait_check_cov()` that did:

```r
if (!is.null(dimnames(M)) && all(traits$trait %in% rownames(M)) &&
    all(traits$trait %in% colnames(M))) {
  M <- M[traits$trait, traits$trait, drop = FALSE]
} else if (nrow(M) != nrow(traits) || ncol(M) != nrow(traits)) {
  ng_stop(...)
}
```

So a **fully** labelled, permuted matrix reaching `ng_add_multitrait_score()` was
already reordered correctly before it hit the positional transform. The brief's
scenario ("a correctly-labelled matrix that arrives reordered silently gets the wrong
sign and scale") was therefore *not* a live defect on that path.

**The real hole was the `else if` fallthrough**, and it is worse than a permutation
bug because it is silent for input that *looks* labelled:

* a matrix whose labels contain a **misspelling** (`"lodgingg"`) fails the
  `%in%` test, matches on dimension, and was then read **positionally in whatever
  order it arrived**;
* a matrix with labels on **one dimension only** (rownames but no colnames -- a
  realistic partial JSON serialisation) did exactly the same;
* `ng_multitrait_cov_to_value_z()` itself, and any direct caller of the
  `*_fit()` functions, had no protection at all;
* and `ng_posterior_multitrait_cross_predict()` (`R/32_posterior_multitrait.R`)
  actively **stamped** trait names onto a user-supplied P positionally
  (`dimnames(P_hat) <- list(trait_names, trait_names)`), destroying whatever labels
  the caller sent and defeating any by-name check downstream.

**3. `ng_estimate_genetic_covariance()` is NOT internal.** It is exported
(`NAMESPACE:39`) and aliased in `man/nextgenCrossDesign-api.Rd`. The brief described
it as `:::`-internal with no help page. It has a help page (the shared API page); what
it lacked was any statement of what `Y` must contain.

**4. The brief's proposed guard trigger does not work on the `two_stage_ridge`
engine.** See Change 3 below -- the fix was to use an engine-independent REML fit
rather than that engine's own `sigma_e2`.

Confirmed, as the brief stated: `jsonlite::fromJSON(jsonlite::toJSON(M))` returns an
**unnamed** matrix. Asserted in the new test.

---

## Change 1 -- the covariance transform is label-aware

New shared helpers in `R/19_multi_trait_selection.R`:

* `ng_multitrait_align_cov(M, trait_names, name)` -- alignment and the label contract.
* `ng_multitrait_validate_cov(M, trait_names, name)` -- alignment plus the existing
  validity checks (finite, symmetric within `1e-8`, positive diagonal, positive
  semidefinite by `eigen`).

Contract, now enforced everywhere P/G enters:

| input | behaviour |
|---|---|
| dimnames on **both** dimensions | reordered **by name** to `traits$trait`; labels must match the trait set **exactly** |
| missing / extra / misspelled label | **hard error**, naming the missing and the unexpected labels |
| duplicated label | hard error, naming the duplicate |
| dimnames on **one** dimension only | **hard error** (previously a silent positional read) |
| **no** dimnames | accepted, **interpreted positionally** in `traits$trait` order, dimension must match; documented in the code and in the error text |

Applied at four sites:

1. `ng_multitrait_cov_to_value_z()` -- aligns `M` by name before building
   `L = diag(sign / scale)`, and takes `value_scales` by name too so scale and sign
   can never disagree with the matrix.
2. `ng_multitrait_index_covariance()` -- `ensure_pxp()` replaced by
   `ng_multitrait_validate_cov()`.
3. `ng_add_multitrait_score()` -- the fallthrough above replaced by
   `ng_multitrait_validate_cov()`.
4. `R/32_posterior_multitrait.R` -- the positional `dimnames(P_hat) <- ...` stamp
   replaced by `ng_multitrait_align_cov()`.

### The numeric difference this prevents

Three traits (`yield` maximize, `disease` minimize, `lodging` minimize), non-diagonal
P and G, permutation `(lodging, yield, disease)`. Measured in
`tests/user_supplied_pg_contract.R`:

* permuted **labelled** matrix, before and after: reordered by name -> Smith-Hazel
  coefficients and `multi_trait_score` **bit-identical** to the correctly-ordered run
  (`identical()`, not a tolerance).
* the same permuted matrix read **positionally** (i.e. what an unlabelled JSON round
  trip delivers, and what a mislabelled matrix used to get):
  **max |db| = 0.1708** on coefficients normalised to `sum(|b|) = 1`, and
  **Spearman rank correlation of the emitted index = 0.9168** against the correct
  run. That is a different ranking of candidate crosses from the same inputs, with
  no error and no warning.
* misspelled label, missing trait, extra trait, one-sided dimnames: all now error and
  name the offending labels.

---

## Change 2 -- `desired_gain` no longer requires P

Verified against the code, not assumed. In `ng_multitrait_index_covariance()` the
`desired_gain` branch returns `target_matrix = G`, `projection = NULL`. In
`ng_multitrait_solve_index()`:

```r
rhs <- if (!is.null(cov_info$projection)) as.numeric(cov_info$projection %*% target) else target
```

With `projection = NULL` the right-hand side is the desired-gain target itself and the
solve is against `target_matrix = G` alone: **`b = G^{-1} d`. P is genuinely absent
from the Pesek-Baker coefficient solve.**

**What `response_P` actually feeds** -- exactly one thing:

```r
sigma_i   <- sqrt(max(as.numeric(crossprod(coefficients, cov_info$response_P %*% coefficients)), 0))
predicted <- as.numeric(cov_info$response_G %*% coefficients) / sigma_i
```

`response_P` supplies the index standard deviation `sigma_I = sqrt(b' P b)` that
standardises the **reported** `predicted_response`. That value is stored in
`attr(scores, "multi_trait")$desired_gain_predicted_response` and copied into
`summary$multitrait_desired_gain_predicted_response`. A grep of every package file
that mentions `predicted_response` shows those are its **only** consumers -- nothing
downstream computes with it, so it can degrade to `NA` safely.

Now:

* `economic_index` requires **both** P and G (message unchanged: still contains
  `"requires both phenotypic_covariance"`, so the existing regression assertions in
  `tests/multitrait_index_math.R` and
  `tests/multitrait_economic_index_external_cov.R` still hold).
* `desired_gain` requires **only G**. Without G it errors with
  `"desired_gain requires genetic_covariance (G) for the Pesek-Baker solve b = G^{-1} d"`.
* Without P, `desired_gain` still solves. Coefficients and the emitted
  `multi_trait_score` are **bit-identical** to the P+G run (asserted with
  `identical()`); `predicted_response` and `index_sd` come back `NA`, and the result
  reports which quantities are unavailable and why:
  * `attr(...)$desired_gain_unavailable` -> `"predicted_response"`
  * `attr(...)$desired_gain_unavailable_reason` -> a sentence naming P and the
    `sqrt(b' P b)` scaling, and stating the coefficients are unaffected
  * both are carried into the plan summary as
    `multitrait_desired_gain_unavailable` / `..._unavailable_reason`.
* The "candidate-score covariance is not a substitute for quantitative-genetic
  covariance" refusal is intact in **both** error messages.

---

## Change 3 -- what `Y` must be for `ng_estimate_genetic_covariance()`

### Documentation

A prominent block above the (exported) function in `R/31_genetic_covariance.R` states
the contract: **BLUEs** valid (with the note that unbalanced stage-one designs ideally
enter inverse-SE weighted, which this function does not do); **BLUPs** invalid as
supplied (`Var(BLUP) = sigma2_g - PEV` -> too little genetic variance, almost no
residual, sigma2_g biased low, h2 -> 1, genetic correlations distorted when
reliabilities differ), remedy deregress (Garrick-Taylor-Dekkers 2009) or correct
correlations by `sqrt(r_t * r_s)` or supply a validated G; **GEBVs** worst and
circular (linear function of the same markers the GRM is built from, so residual -> 0,
h2 -> 1, and G-hat is the covariance of predictions, understating sigma2_g by roughly
the reliability). `man/nextgenCrossDesign-api.Rd` carries a condensed version.

### The guard, and why the obvious implementation does not work

The brief's suggested trigger was "the fitted residual variance is at or near zero, or
the implied h2 is at or near 1". Taking those from the fitting engine's own outputs
**fails on `two_stage_ridge`**: its `sigma_e2` is an in-sample ridge residual at a
CV-selected lambda, which does not collapse on shrunk input. Measured on noiseless
GEBV-like `Y` (n=80, m=300), the engine reported residual/observed variance ratios of
**0.758 and 0.934** and implied h2 of **0.858 and 0.322** -- the guard would never have
fired. A Haseman-Elston moment estimator was also prototyped and rejected: far too
noisy at breeding-programme n (GEBV input returned h2 0.65 and 0.43).

The implemented guard is therefore **engine-independent**:
`ng_genetic_cov_grm_reml_h2()` runs an EMMA-style profile REML of each trait against
the GRM,

```
y = 1 mu + g + e,   Var(y) = sigma_g2 K + sigma_e2 I
```

profiled over `delta = sigma_e2 / sigma_g2` on a single `eigen(K)` shared across
traits. It is computed **only** for this check and never used to build G-hat.

**Trigger:** any trait with implied `h2 >= 0.99` **or** REML residual variance /
observed variance `<= 0.01`.

**What it says:** names the offending traits, prints their implied h2 and
residual/observed ratio, states that `Var(BLUP) = sigma2_g - PEV` biases sigma2_g low
and pushes h2 to 1, that GEBVs are circular and make G-hat the covariance of
predictions, lists the three remedies (deregress / `sqrt(r_t * r_s)` / supply a
validated G), and ends "Not an error: G-hat is returned, and is biased low."
It **warns and never errors**.

**Skipped, with a stated reason in `implied_heritability_note`**, for `n < 100`
(estimate too erratic to act on -- measured: at n=80 the REML h2 for GEBV input was
0.588, a false negative) and for `n > 2000` (the n x n eigen stops being cheap).

Measured behaviour at n=200, m=500 (in the new test):

| `Y` | guard REML h2 | REML residual / observed | guard |
|---|---|---|---|
| phenotypes (`g + e`) | 0.5650, 0.5151 | 0.440, 0.477 | silent |
| GEBV-like (`g`, no residual) | 1.0000, 0.9861 | 1.04e-09, 1.31e-02 | **fires** |

The sommer REML engine corroborates independently: on the same GEBV input its own
residual covariance comes back **0.00e+00** for both traits.

### Provenance attributes

Existing (`method`, `requested_method`, `formal_variance_component_estimate`,
`genetic_correlation`, `n_used`) were **not** enough -- a caller could not say which
traits the matrix covers, how many markers backed it, what `Y` was assumed to be, or
what the fit implied about heritability. Added:

`traits`, `n_markers`, `y_input_contract`, `genetic_variance`, `residual_variance`
(from the fitting engine), `implied_heritability`, `implied_heritability_method`,
`implied_heritability_note`, `reml_genetic_variance`, `reml_residual_variance` (from
the guard), `shrunk_input_suspected`, `shrunk_input_traits`.

---

## Verification

### Bit-identity against a pristine `git archive HEAD` tree

A harness ran the same fixed-seed workload under `/tmp/ngcdpg` (pristine 0.26.0 HEAD)
and under the patched working tree: a full single-trait `ng_run_cross_prediction()`
(OCS allocation, usefulness metric), a multi-trait `weighted` run with no P/G at all,
and `economic_index` + `desired_gain` runs with correctly-ordered **labelled** P and G.
14 captured quantities (selected crosses, scored table, all three index score vectors,
both coefficient/target/predicted-response/covariance sets) compared with
`identical()` and `all.equal(tolerance = 0)`:

```
single_gain    identical=TRUE  all.equal(tol=0)=TRUE
single_pairs   identical=TRUE  all.equal(tol=0)=TRUE
single_scored  identical=TRUE  all.equal(tol=0)=TRUE
mt_weighted    identical=TRUE  all.equal(tol=0)=TRUE
mt_economic    identical=TRUE  all.equal(tol=0)=TRUE
mt_desired     identical=TRUE  all.equal(tol=0)=TRUE
ei_coef        identical=TRUE  all.equal(tol=0)=TRUE
ei_pred        identical=TRUE  all.equal(tol=0)=TRUE
ei_target      identical=TRUE  all.equal(tol=0)=TRUE
ei_cov         identical=TRUE  all.equal(tol=0)=TRUE
dg_coef        identical=TRUE  all.equal(tol=0)=TRUE
dg_pred        identical=TRUE  all.equal(tol=0)=TRUE
dg_target      identical=TRUE  all.equal(tol=0)=TRUE
dg_cov         identical=TRUE  all.equal(tol=0)=TRUE
OVERALL identical: TRUE
```

### New test

`tests/user_supplied_pg_contract.R` -- **PASS (20 checks)**, auto-discovered by
`tools/run_tests.R` (the deep-harness CI job globs `tests/*.R`, so no registry edit
was needed).

### Existing suites

All PASS: `multitrait_index_math`, `multitrait_economic_index_external_cov`,
`multitrait_exact_cov_integration`, `multi_trait_selection`,
`multi_trait_validation_smoke`, `multi_trait_validation_grid`,
`multi_trait_crop_validation`, `genetic_covariance_estimator`,
`breeder_selection_objective`, `runner_multitrait_index_posterior`,
`threshold_probability_multitrait`, `priority_risk_portfolio_multitrait`,
`cross_trait_within_family_cov`, `backend_capability_registry`,
`contract_schema_drift`, `strategy_exports`, `runner_integration`,
`run_cross_prediction_json_messages`, `user_cross_prediction_workflow`,
`crossing_plan_workbook`, `cross_priority_workbook`,
`vignette_chronological_examples`.

### Release gate

`R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib` from a
`/tmp` copy (the repo lives on OneDrive, where in-place builds hang), then
`NGCD_RELEASE_LIB=... Rscript tools/run_statistical_release_gate.R`:

```
Statistical code gate: PASS (44/44)
```

`packageVersion("nextgenCrossDesign")` from the install library reports **0.27.0**.

---

## Notes for the frontend

* Send P and G **with dimnames on both dimensions**. `jsonlite` drops them on a
  matrix round trip, so serialise the labels explicitly (e.g. as a
  `{traits: [...], matrix: [[...]]}` object) and restore them before calling the
  backend. An unlabelled matrix is accepted but is read **positionally in trait
  order**, and the backend cannot detect a reordering.
* `desired_gain` can now be offered with G only. When P is absent the index is fully
  valid; only the reported predicted response is unavailable, and
  `multitrait_desired_gain_unavailable_reason` is a display-ready sentence explaining
  that.
* `economic_index` must still collect both matrices.
* If the GUI offers "estimate G from my data", label the `Y` upload as phenotypes or
  BLUEs, and surface `shrunk_input_suspected` / `shrunk_input_traits` /
  `implied_heritability` from the returned matrix -- that is the user-visible form of
  the guard.

## Residual concerns

* The guard is a heuristic with a documented blind spot at small n (`n < 100`,
  skipped) and near the tolerance boundary. In the test, a second GEBV trait landed at
  h2 = 0.9861 and was **not** flagged while the first (h2 = 1.0000) was. It reduces
  the chance of a silently biased G-hat; it does not eliminate it.
* `ng_posterior_multitrait_cross_predict()` accepts `genetic_covariance_draws` as a
  bare 3-d array. Slices carry dimnames only if the array does; an unlabelled array is
  still read positionally per the documented contract. Labelling the first two
  dimensions of that array is the caller's responsibility.
* Erroring on an **extra** label is a deliberate tightening: a superset G covering
  more traits than the current index used to be silently subset. No test or internal
  caller relied on that, and the brief asked for exact-set matching, but a user
  passing a program-wide G to a two-trait index will now see an error and must subset
  it themselves.
* `ng_estimate_genetic_covariance()` still fits BLUEs unweighted; inverse-SE weighting
  for unbalanced stage-one designs is documented as a known approximation, not
  implemented.
