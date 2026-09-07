# The estimator refuses its own impossible output, and the dev loader loads only the package (0.30.0)

Two changes.

1. `ng_estimate_genetic_covariance()` now checks the `G_hat` it just produced against the
   phenotypic variance of the data it was fitted to, and BLOCKS when the implied per-trait
   heritability exceeds 1. Both engines. No rescaling, no clamping.
2. `ng_load()` now sources exactly what `R CMD build` would include, instead of globbing
   `R/[0-9]*.R` off the filesystem and picking up eleven untracked legacy scripts that are not
   part of the package.

Standard applied throughout, from the package owner: **prefer what is statistically defensible,
and BLOCK rather than warn**, because an invalid `G` propagates into the index weights, the cross
ranking, the allocation and the joint probabilities alike, and every downstream number then looks
plausible.

---

## 1. Change 1 — the rule

After the engine returns and before anything else is attached, the estimator computes

```
h2_t = diag(G_hat)_t / stats::var(Y[, t])
```

on the SAME complete-case rows the fit used (`Y` has already been intersected with `geno`'s ids
and reduced to `complete.cases(Y)` at that point), and refuses when

```
diag(G_hat)_t  >  var(Y[, t])  +  ng_cov_tol(max |c(vp, vg)|)
```

for any trait. `ng_cov_tol(scale) = 1e-8 * max(1, |scale|)` is the package-wide relative
covariance tolerance introduced in 0.29.0 — deliberately the SAME helper, and the same form of
comparison, that `ng_multitrait_validate_cov_pair()` uses for its per-trait `diag(G) <= diag(P)`
special case, so the estimator's self-check and the downstream pair guard cannot disagree about
where the boundary is.

Implemented in `ng_genetic_cov_require_possible_h2()` (`R/31_genetic_covariance.R`), called from
`ng_genetic_cov_estimate()` immediately after `dimnames(G_hat)` is set and BEFORE the advisory
already-shrunk-input warning, so an impossible `G_hat` is reported as the hard error it is rather
than being preceded by a softer complaint about a different problem.

### Which phenotypic variance, and why

The reference is the per-trait SAMPLE variance of the fitted rows — `stats::var(Y[, t])`, which is
already computed in the function as `observed_variance` and is reused rather than recomputed.
(`diag(stats::var(Y))` is the same quantity mathematically but takes a different code path in
`stats` and can differ in the last bit; the per-column form is what the guard uses and what the
tests assert.)

Three reasons for that rather than `ng_estimate_phenotypic_covariance()`:

1. **It is not itself a modelling choice.** No shrinkage intensity, no PSD projection, no second
   estimator whose own error could either excuse or manufacture a violation. The objection is
   between `G_hat` and the data, with nothing in between.
2. **It is never stricter than the downstream pair guard.** With `shrinkage = "none"` the P
   diagonal IS this variance; with `shrinkage = "auto"` the Ledoit-Wolf target is
   `diag(diag(S))` with `S` the `1/n` sample covariance, so the shrunk P's diagonal is SMALLER
   than this unbiased `1/(n-1)` one. Anything this guard passes,
   `ng_multitrait_validate_cov_pair()` also passes on the same data. Verified directly in
   `tests/estimator_selfguard_bit_identity.R`.
3. **It is the more forgiving of the two conventions**, so a violation it reports is not an
   artefact of picking the tighter denominator. The brief's own figures show the difference: at
   seed 2026 the implied `h2` against `P_hat` is `[2.533, 0.531, 1.158]` and against `var(Y)` it
   is `[2.507, 0.526, 1.147]` — smaller by exactly `n/(n-1) = 100/99`.

**The caveat, stated rather than hidden.** Under `y = mu + g + e` with `Var(g) = G (x) K`, the
individual-level phenotypic variance is `G_tt * K_ii + R_tt`, so the exact model-implied ratio
carries a `mean(diag(K))` factor. For this package's VanRaden GRM that factor is ~1 for an outbred
panel (measured **0.9955, range [0.977, 1.016]** over the 200 simulated panels below) and ~2 for a
fully inbred DH/RIL panel — i.e. `>= 1` in the cases that matter, which makes the plain ratio
CONSERVATIVE: it understates `h2` and therefore under-blocks rather than false-blocking. The guard
deliberately does not apply the factor, because it is a statement about the returned matrix
against the observed data, not a re-derivation of the fitted model.

### Example message

`ng_estimate_genetic_covariance(geno, Y, method = "two_stage_ridge", kfold = 5, seed = 2026)` on
the generating model in `tests/genetic_covariance_estimator.R`:

```
ng_estimate_genetic_covariance(method = "two_stage_ridge") produced a genetic covariance matrix
that is IMPOSSIBLE for the data it was fitted to: genetic variance cannot exceed phenotypic
variance, because h2 = sigma2_g / sigma2_p is a ratio of a part to the whole and must lie in
[0, 1]. trait 'trait_a' was given genetic variance 3.59757 against an observed phenotypic variance
of 1.43476 (the sample variance of that trait over the 100 complete-case rows the fit used),
implying h2 = 2.507; trait 'trait_c' was given genetic variance 1.00416 against an observed
phenotypic variance of 0.87571 (the sample variance of that trait over the 100 complete-case rows
the fit used), implying h2 = 1.147. two_stage_ridge builds the diagonal by GBLUP lambda-inversion
(sigma_g2 = sigma_e2 * denom / lambda) from a lambda selected by cross-validated PREDICTION, not
estimated as a variance-component ratio. That heuristic's variance SCALE is unreliable -- measured
over 200 simulated datasets it puts the worst per-trait genetic variance out by a factor of about
6 on average and up to about 1000, and it implies h2 > 1 on roughly 86% of them. Its genetic
CORRELATIONS are usable (relative Frobenius error 0.226 on average, 0.479 at worst over the same
200 datasets); its covariance scale is not. Nothing is returned, and G-hat is deliberately NOT
rescaled, clamped or projected onto the feasible set: a silently corrected G would produce index
weights that cannot be traced back to anything. Remedies: supply a validated genetic_covariance
directly (from a designed multi-environment trial, or from published variance components for the
crop and trait); or use method = "sommer_remml", the formal multivariate variance-component
estimator, which is what method = "auto" resolves to. If you only need the genetic CORRELATION
structure, ng_estimate_genetic_correlation() returns it without the variance scale this check
refuses.
```

### What fraction it blocks

200 independently simulated datasets from the estimator test's own generating model
(n = 100, m = 200, 3 traits, every true `h2 = 0.5`), each with its own estimator seed:

| Quantity | value |
| --- | --- |
| datasets with ANY implied `h2 > 1` — **REFUSED** | **171 / 200 = 85.5%** |
| trait-level violations | 47% of the 600 (trait, dataset) pairs |
| per-trait implied `h2`, median | 0.613 / 1.022 / 0.594 |
| per-trait implied `h2`, mean | 1.038 / 1.142 / 0.968 |
| per-trait implied `h2`, max | 9.84 / 9.86 / 4.61 |
| `mean(diag(K))` across the 200 panels | 0.9955, range [0.977, 1.016] |

The brief's `~87%` was measured against `P_hat`; against `var(Y)` it is 85.5%, the small
difference being the `n/(n-1)` denominator discussed above. Both numbers say the same thing.

The `tests/` suite now asserts the rate rather than only the seed: over the 8 datasets the
correlation assertion already uses, **7 of 8 are refused**, and every verdict is cross-checked
against `any(diag(G_hat) > var(Y))` computed independently — `identical()`, not "approximately
agrees".

### It blocks for `sommer_remml` too

**Decision: yes, the same rule, unconditionally.** Three reasons, in order of weight.

1. **`desired_gain` never sees a `P`.** Pesek-Baker solves `b = G^{-1} d` from `G` alone, so
   `ng_multitrait_validate_cov_pair()` — which only runs when a caller supplies BOTH matrices —
   never runs on that route. On the `desired_gain` path the estimator's self-check is the ONLY
   place an impossible `G` can be caught before it becomes index coefficients. Exempting
   `sommer_remml` would leave the DEFAULT engine on the DEFAULT resolution of `method = "auto"`
   unguarded on exactly that route.
2. **The rule is about the matrix, not the machinery.** `h2 > 1` is impossible whichever
   estimator wrote it. A guard that trusts an engine by name is not a guard on the output; it is
   a guard on provenance, and provenance is what the `formal_variance_component_estimate`
   attribute already reports.
3. **Measurement supports it rather than contradicting it.** Over 118 successful sommer fits from
   the same `h2 = 0.5` generating model (120 attempted, 2 engine errors), exactly **one** violates
   — **0.85%**, against 85.5% for `two_stage_ridge` — and inspection shows it is a fit that did not
   succeed, not a good fit being refused. The worked example (that one, seed 5034): sommer returned `sigma_g2 = [0.800, 0.306, 1.047]`
   with residuals `[1.068, 4.372, 0.054]`, i.e. it assigned `h2 = 0.95` to `trait_c` and
   `h2 = 0.065` to `trait_b` against a truth of 0.5 for both, and `trait_c`'s implied ratio
   against `var(Y)` came out at 1.048. A multivariate REML fit of 3 traits on n = 100 with
   `tolParInv = 1e-3` and no convergence assertion is not above scrutiny.

The message differs only in the diagnosis clause: for `sommer_remml` it says the fit itself should
be read as having failed — too few records for the number of traits, a poorly conditioned GRM, or
a run that stopped short of convergence — rather than blaming a heuristic's variance scale.

Honest statement of the residual risk: at n = 100 the sample variance in the denominator is itself
noisy, so a REML fit genuinely near `h2 = 1` can be tipped over the line by sampling error in the
denominator. That is a real false-positive mode. It is accepted because (a) the failure is loud,
named and numeric rather than silent, (b) the remedies in the message are the right actions in
that situation anyway, and (c) the owner's standard is to block, and a `G` sitting at `h2 ~ 1`
should not be feeding a Smith-Hazel solve without the breeder looking at it.

### No rescaling, and the new attributes

`G_hat` is never clamped, divided through by the implied `h2`, or projected onto the feasible set.
Nothing is returned at all when the guard fires. Two attributes now report the quantity the guard
judged, whether or not it fired, so a caller can see how close its own estimate came:

* `phenotypic_variance_observed` — the per-trait `stats::var(Y[, t])` used as the denominator;
* `implied_h2_vs_observed_variance` — the ratio itself.

These are NOT the existing `implied_heritability`, which remains an independent profile-REML `h2`
against the GRM computed only for the already-shrunk-input warning. The distinction is stated in
the code comment at the attachment site, because two attributes both called a heritability would
otherwise be a trap.

### `ng_estimate_genetic_correlation()` — implemented

The brief asked whether the correlation matrix should become a first-class result, to be
implemented only if small and clean. It was, so it was.

`two_stage_ridge` splits in two, exactly as `R/31`'s own comments predict and as 0.29.0 measured:
its off-diagonals are the Pearson correlation of the per-trait ridge `beta` vectors, and Pearson
correlation is invariant to per-trait multiplicative shrinkage, so it survives the ridge
attenuation (relative Frobenius error mean 0.226, max 0.479 over 200 datasets, against mean 1.88
and max 19.3 for the covariance). Blocking the covariance route without exposing the correlation
would have removed the only usable output of an estimator the package still ships.

The addition is:

* a new exported `ng_estimate_genetic_correlation(geno, Y, method, ridge_lambda, kfold, seed)`;
* the old body of `ng_estimate_genetic_covariance()` moved verbatim into an INTERNAL
  `ng_genetic_cov_estimate(..., require_possible_h2)`, which the two public functions call.

That structure matters: the bypass is not reachable from any exported signature. There is
deliberately no `check = FALSE` argument on the covariance estimator, because an escape hatch on
the public function is how a hard block becomes a suggestion.

The correlation result carries the same provenance attributes plus `scale = "correlation"`, and
the h2 check is not applied to it — a correlation matrix has a unit diagonal and makes no claim
about genetic variance for the check to object to. The docs say plainly that it is not a route
back to a refused covariance: multiplying it up by any variance reintroduces exactly the scale
that was refused.

### The one live caller, checked

`~/NextGenCrossDesign/inst/app/tools/run_cross_prediction_json.R:1104` calls
`nextgenCrossDesign::ng_estimate_genetic_covariance(geno = ..., Y = Y)` with NO `method`
argument, so it takes `method = "auto"` and resolves to `sommer_remml` — the engine this guard
refuses on under 1% of datasets, not the heuristic it refuses on 86%. It is the SECOND rung of
that file's cross-trait covariance ladder, taken only when the exact within-family `wf_var_*` /
`wf_cov_*` columns are absent, and the call is already wrapped so that a failure produces the
frontend's own refusal:

> `no cross-trait covariance is available for this run: the exact within-family columns
> (wf_var_* / wf_cov_*) are absent and ng_estimate_genetic_covariance() failed, so a joint
> probability could only be computed by assuming the traits are independent. It is not computed
> rather than reported as if the traits were uncorrelated.`

So a refusal there degrades exactly the way that file already intends: loudly, into "not
computed", never into a silent independence assumption. Nothing in the frontend needs changing to
be correct. Worth flagging for whoever owns it (frontends are out of scope for this repository):
that rung only ever uses `G_hat` as a CORRELATION — it is passed to
`ng_add_p_superior_progeny_multitrait(G_hat = ...)` — so `ng_estimate_genetic_correlation()` is
the more honest call there, and would not be subject to this guard at all.

Nothing inside the package itself calls `ng_estimate_genetic_covariance()`: `git grep` over `R/`
finds one comment reference in `R/33_threshold_probability_multitrait.R` and no call. The guard
therefore cannot break any pipeline path.

### Out of scope, noted

`ng_posterior_genetic_covariance()` runs the same lambda-inversion and produces the same
impossible variance scale. It is NOT guarded here, for two reasons: the brief scopes Change 1 to
`ng_estimate_genetic_covariance()`, and that function already refuses to run at all without
`allow_heuristic = TRUE` and warns that it returns heuristic draws from separate univariate ridge
fits. It is a candidate for the same treatment in a later release.

---

## 2. Change 2 — `ng_load()` loads the package, and only the package

### What was wrong, and how it actually manifested here

`ng_load()` globbed `R/[0-9]*.R` off the filesystem. That is a directory listing, not a package
manifest. This working tree carries ELEVEN untracked legacy prototype scripts in `R/` whose names
collide with the package's own numbering, so `sort()` interleaved them with the real sources and
`sys.source()`d them into `.GlobalEnv` — several of them AFTER the package file they collide with,
so a colliding definition would have won.

The brief's diagnosis of the mechanism needs one correction, and it is worth stating precisely
because it changes what a fix has to detect.

**The brief says reading these files returns `Operation timed out`. From R, it does not.** They are
OneDrive cloud-only placeholders, and from R they open successfully, report their real size, and
read as ZERO bytes — with no error and no warning:

```
file.size("R/01_relationships.R")   ->  5558
file.access("R/01_relationships.R", 4) -> 0        (readable)
readBin(file(f, "rb"), "raw", 200)  ->  raw(0)     (0 bytes)
readLines(f, warn = TRUE)           ->  character(0), no warning
```

(`grep` and `head` from the shell DO report `Operation timed out`; R's file layer does not.) So
`ng_load()` did not fail — it silently sourced eleven empty files, taking about 1.2 s each. Before:
`ng_load()` returned `TRUE` in **16.2 s**. That is the worst version of the problem: on THIS
machine the legacy files defined nothing (so no result was ever wrong), and on a machine where
OneDrive has materialised them they WOULD be sourced and could shadow package functions. The file
set a dev/test session exercises was environment-dependent, and nothing anywhere said so.

### The rule chosen: respect `.Rbuildignore`

`ng_load()` now excludes exactly what `.Rbuildignore` excludes, via
`ng_load_apply_rbuildignore()`, which applies the same rule `R CMD build` does: each non-blank,
non-comment line is a Perl regular expression matched UNANCHORED against the path relative to the
package root (`grepl(p, rel, perl = TRUE)`; a missing `.Rbuildignore` excludes nothing).

Justification, against the two alternatives the brief offered:

* **It makes dev-load and the built package the SAME definition, not two lists that drift.**
  `.Rbuildignore` is the single existing declaration of what the package contains; deriving the
  load set from it means a file added to `R/` is loaded exactly when it would be shipped. The
  eleven legacy files are already listed there, with a comment explaining why — that declaration
  existed and `ng_load()` was simply not consulting it.
* **It needs no external tooling.** A git-tracked-files rule is correct on a clone and wrong or
  unavailable everywhere else: an unpacked tarball, a CI image, a copied directory, a git-less
  machine. The brief asked for robustness when git is absent; the cleanest way to be robust when
  git is absent is not to need git.
* **An explicit manifest in `load.R` would be a third list to keep in sync** with `.Rbuildignore`
  and the filesystem, and the failure mode of a stale manifest is silently not loading a real
  source file.

Two files the tarball DOES contain are still not sourced, for reasons unrelated to
`.Rbuildignore`, and the code says so: `load.R` itself (the bootstrap) and `RcppExports.R`, whose
`.Call()` stubs dispatch into a compiled DLL a source-loaded session does not have —
`ng_load_cpp()` supplies those symbols through `Rcpp::sourceCpp()` instead. Both are outside the
`^[0-9]` glob, so the rule never sees them.

**Verified equivalence.** The resulting set is exactly the 54 files a `git archive HEAD` tree
contains, and exactly the git-tracked `R/[0-9]*.R` set — `identical()` on both comparisons.

**Verified git-less.** Run against `/tmp/ngcd30` — the exported 0.30.0 tree, which has no `.git`
directory at all — `ng_load()` returns `TRUE` and loads 54 of 54 numbered files. A
git-tracked-files rule would have had nothing to consult there; `.Rbuildignore` travels with the
tree.

### Files no longer sourced

54 of the 65 files matching `R/[0-9]*.R` are loaded. The 11 now excluded, all untracked, all
already `.Rbuildignore`d, all OneDrive placeholders on this machine:

```
R/01_relationships.R              R/06_build_cross_data.R
R/02_duplicate_detection_legacy.R R/07_optimizers.R
R/02b_duplicate_detection.R       R/08_select_optimal_parents.R
R/03_ld.R                         R/09_pmv_cpp.R
R/04_variance_simple_uc.R         R/10_marker_effects.R
R/05_pmv.R
```

None was deleted, moved or modified. They are the user's files and the working tree is untouched.

### Unreadable is now loud, not silent

`ng_load_require_materialised()` runs over every file that WILL be sourced, before any of them is.
`file.size()` alone is not enough — a placeholder reports its real size — so it opens the file and
checks that at least one byte actually comes back:

```
ng_load(): the package source file R/05_pmv.R reports 60582 bytes but reads as EMPTY. This is the
signature of a cloud-storage placeholder (OneDrive/iCloud 'Files On-Demand') whose contents are
not materialised on this machine. Sourcing it would silently define nothing and leave the loaded
package incomplete. Download/pin the file locally, or remove it from R/, and re-run.
```

A file that cannot be stat'ed at all gets its own named error. A genuinely zero-byte file is
allowed through (there is nothing to read).

On hanging: the check adds no new hang risk, because it only touches files that are about to be
sourced anyway. The eleven that used to cost 1.2 s each are no longer opened at all. `ng_load()`
went from **16.2 s to 0.20 s** on this tree, and all of the difference was placeholder I/O buying
nothing.

The same cause explains the `alphamate_external` failure the 0.29.0 report recorded: the untracked
`external/AlphaMate/binaries/AlphaMate.exe` is a placeholder and `sh` returns 126. That is
environmental and is not addressed here.

---

## 3. Bit-identity

`tests/estimator_selfguard_bit_identity.R` compares the 0.30.0 tree against a **pristine 0.29.0
tree** (`git archive HEAD` of `7ec06cd`, unpacked to `/tmp`, loaded with the same helper), using
`identical()` at `tolerance = 0`. Reference values are IEEE-754 hex float literals (`%a`), which
round-trip exactly where decimal `%.17g` does not.

The fixture is a dataset from the same generating model at a seed (5031) chosen because
`two_stage_ridge`'s implied `h2` comes out BELOW 1 on all three traits — one of the ~15% the guard
lets through. The point is that on such input the guard must be invisible.

```
  implied h2 vs var(Y): [0.5686, 0.5725, 0.6158] -- all below 1, so the guard is silent
  OK: the guard-satisfying fixture also passes ng_multitrait_validate_cov_pair()
  OK: G_hat and its genetic_correlation are bit-identical to 0.29.0 (tolerance = 0)
  OK: provenance attributes and P_hat are bit-identical to 0.29.0 (tolerance = 0)
  OK: implied_h2_vs_observed_variance reports exactly what the guard judged
  OK: ng_score_crosses() under the reduced ng_load() file set is bit-identical to 0.29.0
  ng_load() sources 54 of the 65 files matching R/[0-9]*.R
  build-ignored, no longer sourced: 01_relationships.R, 02_duplicate_detection_legacy.R,
    02b_duplicate_detection.R, 03_ld.R, 04_variance_simple_uc.R, 05_pmv.R, 06_build_cross_data.R,
    07_optimizers.R, 08_select_optimal_parents.R, 09_pmv_cpp.R, 10_marker_effects.R
  OK: ng_load() loads exactly the non-build-ignored R/[0-9]*.R files
  OK: ng_load_require_materialised() refuses, by name, a file it cannot read
estimator_selfguard_bit_identity: PASS (7 checks)
```

Covered at `tolerance = 0`: `G_hat` (all 9 entries), its `genetic_correlation`, its
`genetic_variance` / `residual_variance` / `method` / `n_used` attributes, `P_hat` (all 9 entries),
and — for the loader change — a full `ng_score_crosses()` run (190 crosses; `cross_mean`, `pmv`,
`usefulness_pmv`, `rank_score` and the cross keys asserted).

---

## 4. What the guards broke inside the package's own suite, and why

Four `tests/` scripts called `ng_estimate_genetic_covariance(method = "two_stage_ridge")` on
fixtures the new guard correctly refuses. This is the same pattern 0.29.0 saw when Guard 2 found
two impossible `G` fixtures inside the suite on its first sweep.

| Test | Implied `h2` on its fixture | Resolution |
| --- | --- | --- |
| `genetic_covariance_estimator` | `[2.507, 0.526, 1.147]` | The refusal is now ASSERTED (message, traits and remedies). The covariance-scale CHARACTERISATION runs against the internal `ng_genetic_cov_two_stage_ridge()`, so the diagnosis is not lost because the public route closed. The correlation-recovery assertion runs through the new `ng_estimate_genetic_correlation()`. A new check asserts the refusal RATE (7/8) and that each verdict matches `any(diag(G_hat) > var(Y))` exactly. |
| `genetic_covariance_trait_order_invariance` | `yield 3.13, protein 1.94` | A refusal is a verdict and must also be order-invariant, so that is now what is asserted for the covariance entry point: same traits, same numbers, whichever column order `Y` arrived in (compared clause by clause after sorting, since the message iterates traits in supplied order). The matrix-returning assertion moved to `ng_estimate_genetic_correlation()`. `implied_h2_vs_observed_variance` is bounded at 1e-12 relative rather than asserted at `tolerance = 0`, because it inherits the same LAPACK reassociation noise from the PSD projection that `G_hat` already does (its denominator, and the pre-projection `genetic_variance`, ARE bit-identical). |
| `posterior_genetic_covariance` | `y1 1.498` | The invariant under test (posterior diagonals equal point diagonals, because `sigma_g2 = sigma_e2 * denom / lambda` is hyperparameter-conditional) is a property of the ENGINE, so the point estimate is taken from `ng_genetic_cov_two_stage_ridge()` directly. Measured gap after the change: 1.65e-15, unchanged. |
| `user_supplied_pg_contract` | `yield 1.44, protein 1.46` | This section inspects a RETURNED estimate and its provenance attributes, so it is routed through `ng_estimate_genetic_correlation()` — identical pipeline, identical shrunk-input guard, identical provenance. A new assertion first states the distinction the two guards draw: the already-shrunk-input condition is ADVISORY (warns, returns), an impossible variance ratio is FATAL. `required_attrs` drops `genetic_correlation` (the returned object IS the correlation) and gains the two 0.30.0 attributes. |

No fixture was weakened to make a test pass, and no criterion was retuned to a lucky seed.

---

## 5. Test results

Deep harness, `tests/*.R`, excluding `tests/pipeline_integration.R` (exceeds 10 minutes on this
machine, and neither change touches the pipeline's numerics — the loader change is asserted
bit-identical on `ng_score_crosses()` above, and `git grep` over `R/` shows no pipeline path calls
the genetic-covariance estimator at all — one comment reference, no call):

```
136 scripts run

SWEEP: 135 PASS / 1 FAIL
FAILED: alphamate_external.R
```

The single failure is the one 0.29.0 recorded, unchanged and environmental:

```
Error: AlphaMate failed with exit code 126: sh:
.../external/AlphaMate/binaries/AlphaMate.exe: Operation timed out
```

`external/AlphaMate/binaries/AlphaMate.exe` is an untracked OneDrive placeholder — the same cause
as the eleven `R/` files, seen through `sh` rather than through R. It is NOT the estimator guard
and NOT the loader change. 0.29.0's count was 134 PASS / 1 FAIL on 135 scripts; 0.30.0 is 135 PASS
/ 1 FAIL on 136, the extra pass being the new bit-identity file. **Confirmed: `alphamate_external`
is still the only failure.**

`tests/pipeline_integration.R` was excluded from that sweep for time, but was then run separately
for extra assurance and PASSES:

```
pipeline integration test passed
Rscript tests/pipeline_integration.R  1196.12s user 1.70s system 99% cpu 19:58.70 total
```

That is the only script in the suite that exercises the full `ng_run_cross_prediction()` path under
the reduced `ng_load()` file set, so with it the deep harness is **136 PASS / 1 FAIL over 137
scripts**, the one failure still being `alphamate_external`.

`testthat` fast gate, against the INSTALLED 0.30.0:

```
testthat: 80 passed / 0 failed / 0 error / 1 skipped   (NG_TEST_CPP not set)
```

Statistical release gate, against the installed 0.30.0, regenerated into
`docs/STATISTICAL_RELEASE_GATE.md`:

```
Statistical code gate: PASS (44/44)
Package: nextgenCrossDesign 0.30.0
```

Install, from a `/tmp` copy (in-place `R CMD build` does not complete on this OneDrive tree):

```
R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib .   -> * DONE
packageVersion("nextgenCrossDesign")                                            -> 0.30.0
"ng_estimate_genetic_correlation" %in% getNamespaceExports(...)                  -> TRUE
```

---

## 6. Where the brief and the code disagreed

Three places, all minor, none changing what was implemented.

1. **The placeholder files do not time out from R; they read as ZERO BYTES, silently.** The brief
   says reading them returns `Operation timed out`, which is true from `grep`/`head` and false
   from R (`readLines()` returns `character(0)` with no warning, `file.size()` still reports the
   real size). This matters for the fix: a readability check based on `file.exists()`,
   `file.access()` or `try(readLines())` would have passed all eleven. The check has to compare
   the reported size against bytes actually obtained, which is what
   `ng_load_require_materialised()` does. It also explains why the pre-0.30.0 `ng_load()` did not
   fail: it sourced eleven empty files and returned `TRUE`, at a cost of 16 seconds.
2. **The blocking fraction is 85.5%, not ~87%, when measured the way the guard measures.** The
   brief's ~87% is against `P_hat`; against `var(Y)` — the denominator the guard uses, and the
   larger one — it is 171/200. The same difference shifts the seed-2026 implied `h2` from the
   brief's `[2.533, 0.531, 1.158]` to `[2.507, 0.526, 1.147]`, exactly the factor `n/(n-1)`.
3. **"A REML estimate should not normally violate it" holds, but not universally.** It was worth
   checking rather than assuming: the sommer engine violates on 1 of 118 successful fits (0.85%)
   from a true-`h2 = 0.5` model, and that fit is visibly bad (one trait at `h2 = 0.95`, another at
   `0.065`). The brief's instinct — "if it does, that is worth blocking as well" — is what the
   measurement supports.

---

## 7. Files changed

```
DESCRIPTION                                        0.29.0 -> 0.30.0
NAMESPACE                                          export(ng_estimate_genetic_correlation)
NEWS.md                                            0.30.0 entry
R/31_genetic_covariance.R                          self-guard, internal impl split, new estimator, 2 attrs
R/load.R                                           .Rbuildignore rule + materialisation check
man/nextgenCrossDesign-api.Rd                      alias + guard/correlation contract
tests/estimator_selfguard_bit_identity.R           NEW  7 checks
tests/genetic_covariance_estimator.R               refusal asserted; correlation route; refusal rate
tests/genetic_covariance_trait_order_invariance.R  refusal is order-invariant; correlation route
tests/posterior_genetic_covariance.R               point estimate via the internal engine
tests/user_supplied_pg_contract.R                  advisory vs fatal; correlation route; attrs
vignettes/nextgenCrossDesign.Rmd                   guard + ng_estimate_genetic_correlation()
docs/STATISTICAL_RELEASE_GATE.md                   regenerated against 0.30.0
docs/estimator-selfguard-report.md                 this file
```

No untracked file in the working tree was deleted, moved or modified.
