# Blocking validity guards for user-supplied covariance matrices (0.29.0)

Scope: four hard-error guards on `phenotypic_covariance` (P) and
`genetic_covariance` (G) in `R/19_multi_trait_selection.R`, plus two related
defects — the remaining position-derived per-trait seeds in
`R/31_genetic_covariance.R`, and an accuracy criterion in
`tests/genetic_covariance_estimator.R` that was passing by luck.

Standard applied throughout: **prefer what is statistically defensible, and
BLOCK rather than warn**, because an invalid G propagates into the index weights,
the cross ranking, the robust allocation and the joint probabilities, and every
downstream number then looks plausible.

---

## 1. Where the guards live

| Guard | Function | Called from |
| --- | --- | --- |
| 1 symmetry (relative) | `ng_multitrait_validate_cov()` | every P/G entry point |
| 2 `P - G` PSD | `ng_multitrait_validate_cov_pair()` | `ng_add_multitrait_score()` (raw units) and `ng_multitrait_index_covariance()` (post-`value_z`) |
| 3 positive definite | `ng_multitrait_require_positive_definite()` | `ng_multitrait_index_covariance()`, on the inverted matrix only |
| 4 implied correlations | `ng_multitrait_check_implied_correlations()` | `ng_multitrait_validate_cov()`, before the PSD test |

One shared tolerance helper, `ng_cov_tol(scale) = 1e-8 * max(1, |scale|)`, gives
the package a single convention. `scale` is the largest element magnitude for an
ELEMENTWISE comparison (symmetry, implied correlations) and the spectral radius
for an EIGENVALUE comparison (PSD, `P - G`). The `max(1, .)` floor keeps the rule
absolute for small matrices, so this is a strict relaxation of the pre-0.29.0
behaviour and never a tightening.

---

## 2. Guard 1 — relative symmetry tolerance

**The inconsistency.** The symmetry test was absolute (`asym > 1e-8`) while the
PSD test six lines below was relative (`1e-8 * max(1, max|ev|)`). Covariances
carry the square of the trait's unit, so their magnitude reflects the breeder's
choice of unit and nothing about validity. Yield in kg/ha gives covariances of
order 1e6, where agreement to 1e-8 absolute is unreachable for any hand-entered
or spreadsheet-rounded matrix; the same yield in t/ha gives order 1, where it is
trivial.

**The rule.** `max|M - M'| <= 1e-8 * max(1, max|M|)`.

The reference scale is the largest element magnitude, not the spectral radius,
because symmetry is an elementwise comparison: it asks whether `M_ts` and `M_st`
agree relative to the size of the entries being compared. (For a symmetric
matrix `max|M_ij| <= ||M||_2 <= p * max|M_ij|`, so the two references agree up to
a factor of `p`; the elementwise one is the honest match for an elementwise
test.)

**Measured.** A `G` at kg/ha scale (`1e6`) with one cell rounded to 9 significant
digits: absolute asymmetry `1.000e-03`, relative asymmetry `1.111e-10`. Refused
before, accepted now, and the resulting index is identical to the
exactly-symmetric matrix's (the matrix is symmetrised anyway). The same
*relative* asymmetry (1e-3 relative) at kg/ha scale is still refused, so this is
not a blanket relaxation.

**Example message.**

```
genetic_covariance must be symmetric: cov(protein, yield) = 1.2 but
cov(yield, protein) = 1.7, a discrepancy of 0.5 against a tolerance of 9e-08
(1e-08 x the largest entry, 9)
```

---

## 3. Guard 2 — `P - G` must be positive semidefinite (the headline guard)

**The principle.** `P = G + R`, so `R = P - G` is itself a covariance matrix and
must be PSD. This is the multivariate generalisation of `0 <= h2 <= 1`: a
negative eigenvalue of `P - G` means there is a trait contrast `c` with
`c' G c > c' P c` — a linear combination of the traits whose GENETIC variance
exceeds its PHENOTYPIC variance.

**Where it runs.** Whenever both matrices are supplied, and always AFTER by-name
alignment, so trait order can neither confound the verdict nor launder an invalid
pair. It is applied twice, deliberately:

* at `ng_add_multitrait_score()`, in the breeder's RAW units, so the objection is
  raised against the numbers actually typed, regardless of which multi-trait
  method the run goes on to use;
* at `ng_multitrait_index_covariance()`, after the `value_z` congruence
  `M -> L M L` with `L = diag(sign / scale)`. A congruence preserves the Loewner
  order, so `P - G` is PSD in `value_z` units exactly when it is in raw units;
  the second call is a consistency check, not a second opinion.

**Order of the two checks.** The elementwise special case `diag(G) <= diag(P)`
(per-trait `h2 <= 1`) is tested first because its message is far more actionable.
Only when no single trait is individually at fault is the general eigenvalue
objection reported, with the offending eigenvalue, the eigenvector contrast that
realises it, and that contrast's genetic and phenotypic variance.

### A realistic P/G pair 0.28.0 accepted and 0.29.0 catches

The package's own estimators. `two_stage_ridge` `G_hat` paired with
`ng_estimate_phenotypic_covariance(shrinkage = "auto")` `P_hat`, on the
simulation that ships in `tests/genetic_covariance_estimator.R` (n = 100,
m = 200, 3 traits, true `h2 = 0.5` for every trait), at the file's own seed 2026:

```
diag(G_hat) = [3.598, 1.860, 1.004]     diag(G_true) = [1.000, 2.000, 0.500]
implied per-trait h2 from (G_hat, P_hat) = [2.533, 0.531, 1.158]
```

0.28.0 accepted this pair and solved `b = P^{-1} G a` from it. 0.29.0 refuses:

```
genetic variance exceeds phenotypic variance, which implies a heritability above
1 and is impossible: trait 'trait_a' has genetic_covariance variance 3.59757
exceeding its phenotypic_covariance variance 1.42041, implying h2 = 2.533; trait
'trait_c' has genetic_covariance variance 1.00415 exceeding its
phenotypic_covariance variance 0.867174, implying h2 = 1.158. P = G + R requires
the residual variance R = P - G to be a valid covariance matrix. Check that
phenotypic_covariance and genetic_covariance are on the SAME scale and in the
SAME units.
```

This is not a one-seed artefact. Over **200 independently simulated datasets**
from that generating model, `P_hat - G_hat` is PSD for only **13%** of them
(`P_true - G_hat` for 14.5%) — i.e. the shipped pair implies `h2 > 1` on ~87% of
datasets. The cause is diagnosed in section 6.

### The general case a per-trait check cannot see

Every per-trait `h2` legal, but a contrast with `h2 > 1`
(`tests/covariance_validity_guards.R`, check 2c):

```
G = [[4.0, 3.9, 0], [3.9, 4.0, 0], [0, 0, 9]]     (r_g(yield, protein) = 0.975)
P = [[5.0, 1.0, 0], [1.0, 5.0, 0], [0, 0, 20]]    (r_p(yield, protein) = 0.20)
diag(G) < diag(P) elementwise; both PD; but min eigenvalue of P - G = -1.9
```

```
the residual covariance R = phenotypic_covariance - genetic_covariance is not
positive semidefinite, so some combination of the traits is assigned MORE genetic
than phenotypic variance (h2 > 1 for that combination) and no valid P = G + R
decomposition exists. Smallest eigenvalue of P - G is -1.9 against a tolerance of
-2e-07. The trait contrast that realises it is +0.707*yield +0.707*protein
+0.000*lodging, which this pair gives genetic variance 7.9 and phenotypic
variance 6 (implied h2 = 1.317). Check that phenotypic_covariance and
genetic_covariance are on the SAME scale and in the SAME units.
```

Also verified: the same invalid pair permuted (and labelled) is refused
identically, and the guard is silent on a valid pair and on a `G`-only
`desired_gain` run.

---

## 4. Guard 3 — positive DEFINITE for the matrix the index inverts

`economic_index` (Smith-Hazel) solves `b = P^{-1} G a` and inverts **P**.
`desired_gain` (Pesek-Baker) solves `b = G^{-1} d` and inverts **G**. Only PSD was
checked, which permits a singular matrix, and `ng_multitrait_solve_index()`
carries `ridge = 1e-6` and a Moore-Penrose fallback that would produce weights
from it without comment.

### What the un-guarded solve actually returned

Measured directly on a rank-2 `G` (lodging exactly yield + protein genetically),
bypassing the guard:

```
un-guarded b on the singular G = (+0.3333, +0.3333, -0.3333)
|cos angle to G's NULL eigenvector| = 1.000000
b' G b / (b'b * lambda_max)       = 4.938e-13
```

The returned index lies **entirely in G's null space** and carries **zero**
genetic variance. The mechanism: LAPACK reports the zero eigenvalue of that `G`
as `5.33e-15` while the routine's own keep-tolerance is
`3 * eps * lambda_max = 4.0e-15`, so the null direction is retained, receives a
coefficient of order `1/5.33e-15` (≈ 3.8e13), swamps the two real directions, and
survives the `sum(|b|)` normalisation as the whole index. The answer — sign
included — is fixed by floating-point noise in the eigendecomposition.

### The threshold, and why it is set by the ridge

The solve runs on `M + r I` with `r = ridge * mean(diag(M)) = ridge * mean(ev)`.
Along an eigenvector with eigenvalue `lambda` the ridge supplies a fraction
`r / (lambda + r)`:

* `lambda >> r` — the ridge is a rounding correction. This is exactly the MILD
  ill-conditioning the ridge exists to absorb, and it is **not** blocked.
* `lambda <~ r` — the ridge DOMINATES, and the coefficients along that contrast
  come from a software constant rather than from the supplied matrix. The ridge
  does not rescue this solve; it conceals that there is nothing to solve.

In condition-number terms, using `mean(ev) >= lambda_max / p`:

```
r >= 1e-6 * lambda_max / p = (1e-6 / p) * kappa * lambda_min
```

so `r` reaches `lambda_min` once `kappa >= p * 1e6`, i.e. `kappa ~ 1e8` for a
realistic multi-trait index (`p` up to ~100). **`NG_INDEX_MAX_CONDITION = 1e8`**
is therefore the point at which the *default* ridge stops being a correction and
starts being the answer. It coincides with the classical numerical limit:
`1e8 * .Machine$double.eps ~ 2e-8`, so more than half the significant digits in
`b` are already gone. A separate numerical-singularity floor
(`lambda_min <= p * eps * lambda_max`) catches rank deficiency independently.

**The threshold deliberately does not scale with the caller's `ridge`.** Raising
the ridge is a choice of more shrinkage; it does not make a rank-deficient matrix
more informative, it only substitutes more of the constant for more of the data.
An earlier draft of this guard *did* scale the floor with the caller's ridge, and
it wrongly refused `tests/multitrait_index_math.R`'s deliberate `ridge = 1e6`
check on a perfectly conditioned `P` (`kappa = 7.5`). That is the bug the fixed
rule avoids.

**Semidefinite is not blocked for a method that does not invert the matrix.**
Verified: the same singular `G` runs fine under `economic_index`, where `G`
enters only as the forward projection `G a`; a singular `P` is refused there. A
mildly ill-conditioned `G` (`kappa = 6e+04`) is accepted.

**Example message.**

```
desired_gain (Pesek-Baker, b = G^{-1} d) inverts genetic_covariance, so
genetic_covariance must be positive DEFINITE, not merely positive semidefinite.
Its smallest eigenvalue is 3.32161e-17 and its largest is 0.86443, a condition
number of 2.60244e+16 against the limit 1e+08. The matrix is numerically SINGULAR
(rank deficient): the trait contrast +0.577*yield +0.577*protein -0.577*lodging
carries no variance at all, so genetic_covariance does not determine an index
coefficient along it. The ridge legitimately absorbs MILD ill-conditioning, but it
cannot rescue this solve -- it would only substitute that software constant for
the matrix you supplied, and return coefficients that look like a selection index
but are not one. Raising `ridge` does not lift this limit. Drop one of the
redundant traits, or supply a genetic_covariance that is better estimated (more
records, or a shrinkage estimator such as
ng_estimate_phenotypic_covariance(shrinkage = 'auto')).
```

---

## 5. Guard 4 — implied correlations within `[-1, 1]`

`|M_ts| <= sqrt(M_tt M_ss)` is the 2×2 principal-minor condition and is therefore
implied by PSD; this guard adds no mathematical content. It exists entirely for
the message, and runs **before** the eigenvalue test so the actionable objection
is the one that fires:

```
genetic_covariance implies a correlation outside [-1, 1], which no covariance
matrix can have (|cov(t,s)| <= sqrt(var(t) var(s)) for every pair): yield vs
protein: implied correlation 1.31 (covariance 3.70524, sqrt(4 * 2) = 2.82843).
Fix the offending off-diagonal cell(s) or the variances they are referred to.
```

Verified: that same matrix is also non-PSD, and the PSD message does **not**
fire; a correlation of exactly `-1` (a legitimately singular, PSD matrix) is on
the boundary and is **not** refused.

---

## 6. Fix 5 — `R/31` position-derived seeds

### On the brief's premise

The brief states that the frontend calls `ng_estimate_genetic_covariance()`
directly and that "the defect is live in production". **Half of that is right and
half is not, and it matters for how the fix is justified.**

Confirmed: `~/NextGenCrossDesign/inst/app/tools/run_cross_prediction_json.R:1104`
does call `nextgenCrossDesign::ng_estimate_genetic_covariance(geno = ..., Y = Y)`
— as the second rung of the cross-trait covariance ladder, taken only when the
exact within-family `wf_var_*` / `wf_cov_*` columns are absent.

Not confirmed: that call passes **no `method` argument**, so it takes the default
`method = "auto"`, which resolves to `sommer_remml`
(`R/31_genetic_covariance.R:408-415`). That engine does not enter the per-trait
loop and never touches the `seed + j` sites. The defective sites are reachable
only via an explicit `method = "two_stage_ridge"`, or via
`ng_posterior_genetic_covariance()` under `allow_heuristic = TRUE`.

Both are exported, documented, user-callable paths, and a trait's result must not
depend on the column order of `Y` whichever route reaches it, so the fix is
applied — but on that basis, not on "live in production".

### The four sites and how each was treated

Treated exactly as 0.28.0 treated the equivalent sites in `R/39` and `R/32`:

| Site | Was | Now | Why |
| --- | --- | --- | --- |
| `ng_genetic_cov_two_stage_ridge()` per-trait fit | `seed + j` | `seed` | lambda-CV fold split: a nuisance parameter, shared partition makes per-trait CV paired and cannot couple traits |
| `ng_posterior_genetic_covariance(beta_posterior)` per-trait draws | `seed + j` | `ng_trait_rng_seed(seed, trait_names[[j]])` | posterior innovations: streams must stay DISTINCT or they manufacture cross-trait correlation in `R_beta_b` |
| `parametric_bootstrap` per-trait fit | `seed + j` | `seed` | lambda-CV fold split, as above |
| `parametric_bootstrap` residual draws | one `set.seed(seed)`, shared stream consumed in COLUMN ORDER | `ng_with_rng_seed(ng_trait_rng_seed(seed, trait, salt = b * 1000003), rnorm(...))` | trait `j` was getting whatever the stream held after traits `1..j-1` drew theirs |

The `salt = b * 1000003` stride keeps `(trait, draw)` streams apart: two pairs can
collide only if two name hashes differ by an exact multiple of a large prime.
`ng_with_rng_seed()` restores the caller's global RNG state, which the previous
bare `set.seed(seed)` did not.

### Order-invariance evidence

`tests/genetic_covariance_trait_order_invariance.R` permutes the columns of `Y`
(`yield, protein, disease` → `disease, yield, protein`), reorders the answer back
**by name**, and asserts — the same method `tests/trait_order_invariance.R` uses:

```
OK: two_stage_ridge: per-trait lambda, CV r2, sigma_e2 and sigma_g2 identical at tolerance = 0
OK: two_stage_ridge: marker effects, beta correlation matrix and the raw G identical at tolerance = 0
  G_hat (after the PSD projection) differs by 5.773e-15 absolute, 8.891e-16 relative
OK: two_stage_ridge: G_hat agrees to LAPACK reassociation noise (< 1e-12 relative) after the PSD projection
OK: ng_estimate_genetic_covariance(two_stage_ridge): same verdict through the exported entry point
  beta_posterior: max |draw difference| = 1.066e-14 absolute, 8.928e-16 relative
OK: ng_posterior_genetic_covariance(beta_posterior): every draw is order-invariant
  parametric_bootstrap: max |draw difference| = 1.421e-14 absolute, 1.154e-15 relative
OK: ng_posterior_genetic_covariance(parametric_bootstrap): every draw is order-invariant
  beta_posterior implied r(yield, protein) across draws: sd = 0.0595, range = [-0.046, 0.102]
OK: per-trait posterior streams remain DISTINCT (draw-to-draw correlation still varies)
OK: renaming a trait moves its stream -- the key is the trait's identity, not a constant
genetic_covariance_trait_order_invariance.R: PASS (8 checks)
```

Everything the seeds control is bit-identical (`tolerance = 0`). The only residual
difference is `G_hat` **after** `ng_genetic_cov_project_psd()`, at ~9e-16
relative: that is floating-point reassociation inside LAPACK's `eigen()` /
`nearPD` applied to a permuted matrix, not seed dependence. `G_raw`, the input to
that projection, *is* bit-identical. This is asserted below 1e-12 relative and
documented rather than hidden.

Two complementary properties are asserted so the fix cannot be "fixed" the wrong
way: per-trait posterior streams remain **distinct** (a single shared seed would
have made them identical and manufactured correlation), and **renaming** a trait
does move its stream — the key is an identity, not a constant.

---

## 7. Fix 6 — the guard test that passed by luck

### Reconciling the 0.28.0 measurement

The 0.28.0 agent reported that the criterion
`||G_hat - G_true||_F / ||G_true||_F <= 0.5` is met in 53% of 120 base seeds
(mean 0.625, sd 0.346). That number is reproducible, but only under one reading
of "base seed":

| Tree | What was varied | P(ratio ≤ 0.5) | mean | sd | median |
| --- | --- | --- | --- | --- | --- |
| 0.28.0 pristine | estimator seed only, one fixed dataset | 0.600 | 0.578 | 0.330 | 0.402 |
| 0.29.0 | estimator seed only, one fixed dataset | 0.625 | 0.575 | 0.345 | 0.402 |
| 0.28.0 pristine | **fresh dataset per seed** | 0.100 | 1.728 | 1.928 | 1.067 |
| 0.29.0 | **fresh dataset per seed** | 0.175 | 1.698 | 1.887 | 1.119 |

(120 seeds each.) The 53% figure is the *estimator-seed-only* condition. The
criterion the test actually claims — "the estimator recovers `G_true`" — is a
statement about datasets, and under that reading it is met by **10–18%**, not
53%. A 200-seed rerun on the 0.29.0 tree gives 12% (mean 1.88, sd 2.06, median
1.26, max 19.3).

The table also settles a second question: **the `R/31` seed fix does not degrade
the estimator.** In the like-for-like estimator-seed condition the two trees are
statistically indistinguishable (mean 0.578 vs 0.575). But it *is* a reseeding,
and it moves seed 2026's own ratio from **0.374 (passing)** to **1.112
(failing)** — which is precisely the fragility the 0.28.0 agent flagged.

### Diagnosis

The estimator splits cleanly in two, exactly as `R/31`'s own comments predict
(200 independent datasets):

| Quantity | mean | sd | median | q90 | max |
| --- | --- | --- | --- | --- | --- |
| relative Frobenius error, **covariance** | 1.878 | 2.062 | 1.260 | 3.581 | 19.31 |
| relative Frobenius error, **correlation** | 0.226 | 0.067 | 0.221 | 0.311 | **0.479** |
| worst per-trait `log |G_hat_tt / G_true_tt|` | 1.824 | 1.533 | 1.456 | 3.131 | 6.94 |
| implied `max h2` vs `P_hat` | 1.923 | 1.386 | 1.293 | 3.892 | 9.30 |

* **Off-diagonals** are the Pearson correlation of the per-trait ridge `beta`
  vectors. Pearson correlation is invariant to per-trait multiplicative
  shrinkage, so it survives ridge attenuation: the correlation matrix is
  recovered to ≤ 0.479 relative Frobenius error on **200/200** datasets.
* **Diagonals** are the GBLUP lambda inversion
  `sigma_g2 = sigma_e2 * denom / lambda`. All of the error lives here: the worst
  per-trait genetic variance is out by a factor of `exp(1.82) ≈ 6.2` on average
  and up to `exp(6.94) ≈ 1000`. It overshoots often enough that `G_hat` exceeds
  `P_hat` on 87% of datasets — the Guard 2 finding of section 3, seen from the
  other side. The two findings are the same finding.

Secondary properties measured and **rejected** as candidate assertions: all three
correlation signs correct (61%), `|r|` ranking correct (31%), max per-correlation
absolute error ≤ 0.5 (99%, max 0.621 — too thin a margin).

### What the test now does

The criterion was **not** tightened to a new lucky seed. Instead:

* **ASSERTED** — the correlation-structure property, which genuinely holds:
  `||R_hat - R_true||_F / ||R_true||_F <= 0.6`, checked over **8 independent
  datasets** rather than one. The bound is justified by the 200-dataset
  reference (mean 0.226, sd 0.067, max 0.479): 0.6 sits ~5.6 sd above the mean
  and 25% above the worst case observed, so it is a bound on a property, not a
  threshold fitted to a seed. Observed in the test: `0.137, 0.278, 0.155, 0.269,
  0.231, 0.245, 0.044, 0.274` (mean 0.204, max 0.278).
* **CHARACTERISING, not guarding** — the covariance-scale error is computed,
  printed with the 200-dataset reference distribution alongside it, and bounded
  only by a collapse ceiling (three orders of magnitude). The file says in its
  own header that the shipped `<= 0.5` criterion is withdrawn and why.
* **ASSERTED (new)** — that the 0.29.0 `P - G` guard *refuses* the
  `(G_hat, P_hat)` pair on this simulation. That refusal is now a documented
  property of the estimator, not a surprise.
* The Smith-Hazel integration check is re-pointed at the generating model's own
  `G_true` / `P_true` (every true `h2 = 0.5`), because what it tests is that a
  valid pair reaches the Smith-Hazel solve, not the estimator's accuracy.

The structural guarantees (symmetric, PSD, positive diagonals, unit correlation
diagonal, provenance attributes, Ledoit-Wolf shrinkage behaviour) remain asserted
unconditionally.

**Reader's summary:** `two_stage_ridge` is usable for the SHAPE of the genetic
correlation structure and not for the SCALE of genetic variance. That is
consistent with the warning it already emits and with its status as an opt-in
heuristic; `sommer_remml` remains the default and the only formal
variance-component route.

---

## 8. A valid P/G pair is numerically untouched

`tests/covariance_guards_bit_identity.R` compares the 0.29.0 tree against a
**pristine 0.28.0 tree** (`git archive HEAD` of `a42dda6` unpacked to `/tmp`) on a
pair that satisfies every guard. Reference values are embedded as IEEE-754 hex
float literals (`%a`), which round-trip exactly — decimal `%.17g` does not in R,
and was observed to differ in the last bit.

```
OK: economic_index and desired_gain scores are bit-identical to 0.28.0 (tolerance = 0)
OK: index coefficients and index SDs are bit-identical to 0.28.0 (tolerance = 0)
OK: the selected crossing plan is the same set of crosses, in the same order
covariance_guards_bit_identity.R: PASS (3 checks)
```

`identical()`, not `all.equal()`. Covered: both index methods' 60 candidate
scores, the coefficient vectors, the index SDs, the full `multi_trait` metadata
lists, and a 4-cross `ng_optimize_multitrait_mating_plan()` plan with its summary.

---

## 9. What the guards and the seed fix broke, and why

The full `tests/` sweep (135 scripts, `tests/pipeline_integration.R` excluded as
it exceeds 10 minutes on this machine) surfaced six failures. Each was run
against a pristine 0.28.0 tree to establish whether it was a regression.

| Test | Regression? | Cause | Resolution |
| --- | --- | --- | --- |
| `alphamate_external` | No — environmental | The untracked `external/AlphaMate/binaries/AlphaMate.exe` is a OneDrive placeholder; `sh` returns 126 `Operation timed out`. Fails identically with my changes stashed. The pristine tree passes only because `git archive` omits the untracked binary, so the test takes its graceful "binary unavailable" branch. | Not fixed; environmental |
| `metric_name_normalization` | No — pre-existing | `ng_cp__build_ctx()` validates a COMPLETE config (defaults live in `ng_run_cross_prediction()`'s formals, not the builder). A later release added a `ci_level` validation the file's minimal `cfg` never supplied, so it died on `ci_level must be one finite probability in (0, 1)` before reaching a single metric assertion. Fails identically on 0.28.0. | Fixed (out of scope but one line): `ci_level = 0.95` added to the fixture, with a comment about the complete-config contract |
| `posterior_multitrait_usefulness_direction` | No — stale since 0.28.0 | The file reconstructs R/32's per-trait posterior draw independently, using `seed = SEED + 1L`. 0.28.0 changed R/32's key from `seed + j` to `ng_trait_rng_seed(seed, trait)` and left this reconstruction behind, so it has been comparing against a different draw stream — and failing — ever since, independently of the direction property it exists to test. Fails identically on 0.28.0. | Fixed: reconstruction now uses `ng_trait_rng_seed(SEED, "disease")` |
| `posterior_genetic_covariance` | **YES — mine** | See below | Fixed in `R/31` |
| `posterior_multitrait_rng_scoping` | **YES — mine (Guard 2 firing correctly)** | See below | Fixture corrected |
| `threshold_probability_multitrait` | **YES — mine (Guard 2 firing correctly)** | See below | Fixture corrected |

### `posterior_genetic_covariance` — a real invariant the seed fix broke

The file asserts that `beta_posterior`'s `G_mean` diagonals equal the point
`ng_estimate_genetic_covariance()` diagonals exactly, because
`sigma_g2 = sigma_e2 * denom / lambda` is hyperparameter-conditional and fixed
across draws. Under `seed + j` both routines happened to select lambda on the
same fold split. My first pass gave the point estimator a SHARED split (the base
`seed`) and the posterior a per-trait identity-derived one, so they selected
different lambdas and the diagonals diverged by **1.306**.

The cause is that `ng_fit_ridge_effects_posterior()` uses one `seed` for two
jobs that need opposite treatments. The fix separates them: lambda is now
selected up front with `ng_fit_ridge_effects(seed = seed)` — the shared,
order-invariant partition, identical to the point estimator's — and passed in as
an explicit `lambda`, so the identity-derived seed governs only the posterior
innovations. Diagonal gap after the fix: **1.65e-15**.

A second assertion in the same file, that every off-diagonal is attenuated
toward zero relative to the point estimate at a flat `1e-3` slack, was found to
be another luck-dependent criterion. Attenuation is an EXPECTATION statement
about Pearson correlation under added noise, so a pair whose point correlation is
already ~0 has nothing to attenuate and the posterior mean's Monte Carlo error
moves it either way. Measured over 10 posterior seeds on that fixture (point
`|r|` = 0.0041, 0.4496, 0.0643):

```
        pair12   pair13   pair23        <- excess of posterior |r| over point |r|
 [1,]  0.00998 -0.25805 -0.04405
 [2,]  0.00296 -0.27587 -0.04431
 ...
 [8,]  0.01013 -0.35306 -0.06234
max     0.01013 -0.25805 -0.03395
```

The two materially-correlated pairs attenuate on 10/10 seeds; the ~0 pair wanders
in `[-0.004, +0.010]` with no sign preference. The file now asserts the aggregate
attenuation (mean `|r|` 0.173 -> 0.028-0.075 on every seed measured) plus a
per-pair bound at 0.05, ~5x the largest Monte Carlo excess observed.

### Two test fixtures Guard 2 correctly refused

`ng_posterior_multitrait_cross_predict()` pairs the caller's `genetic_covariance`
with a `P_hat` it estimates from `Y` when no P is supplied. Two fixtures passed
an arbitrary identity `G` against phenotypes of a completely different scale:

* `posterior_multitrait_rng_scoping`: `G = diag(n_traits)` (unit genetic
  variance) against phenotypes of variance ~0.2 — implied `h2` = **4.65** and
  **7.60**.
* `threshold_probability_multitrait`: `G = diag(c(yield = 1, disease = 1))`
  against phenotypes of variance 0.92 and 0.97 — implied `h2` = **1.083** and
  **1.026**.

Neither was ever a possible `P = G + R`, and neither test's property (RNG
isolation; threshold-probability calibration) depends on `G`'s scale. Both
fixtures now use `G = 0.5 * P_hat` on the same rows the run uses — a uniform
`h2 = 0.5` with `R = P - G = 0.5 P` positive definite by construction. The guard
found two impossible fixtures inside the package's own suite on its first sweep.

### Final state

```
tests/ sweep: 129 PASS / 6 FAIL before the follow-up fixes
              134 PASS / 1 FAIL after (alphamate_external, environmental)
testthat:     80 passed / 0 failed / 1 skipped (NG_TEST_CPP not set)
release gate: PASS (44/44) against the installed 0.29.0, regenerated
```

---

## 10. Files changed

```
DESCRIPTION                                        0.28.0 -> 0.29.0
NEWS.md                                            0.29.0 entry
R/19_multi_trait_selection.R                       guards 1-4
R/31_genetic_covariance.R                          four seed sites
man/nextgenCrossDesign-api.Rd                      guard contract documented
tests/covariance_validity_guards.R                 NEW  18 checks
tests/covariance_guards_bit_identity.R             NEW   3 checks
tests/genetic_covariance_trait_order_invariance.R  NEW   8 checks
tests/genetic_covariance_estimator.R               rewritten, 6 checks
tests/posterior_genetic_covariance.R               lambda invariant + attenuation bound
tests/posterior_multitrait_rng_scoping.R           impossible G fixture corrected
tests/threshold_probability_multitrait.R           impossible G fixture corrected
tests/posterior_multitrait_usefulness_direction.R  stale 0.28.0 seed reconstruction
tests/metric_name_normalization.R                  missing ci_level in fixture
docs/STATISTICAL_RELEASE_GATE.md                   regenerated against 0.29.0
docs/covariance-guards-report.md                   this file
```
