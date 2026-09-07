# nextgenCrossDesign 0.30.0

## The genetic-covariance estimator refuses its own impossible output

0.29.0 added a hard `P - G` guard for user-supplied covariance pairs. Turning
that same standard on the package's OWN estimator showed that
`ng_estimate_genetic_covariance(method = "two_stage_ridge")` routinely produces
a `G_hat` that no phenotypic covariance can accommodate: on the generating model
that ships in `tests/genetic_covariance_estimator.R` (n = 100, m = 200, 3
traits, every true `h2 = 0.5`), its implied per-trait heritability exceeds 1 on
**85.5% of 200 independently simulated datasets**, and at the test file's own
seed 2026 it returns `diag(G_hat) = [3.598, 1.860, 1.004]` against a truth of
`[1, 2, 0.5]` -- implied `h2 = [2.507, 0.526, 1.147]`. Backend 0.28.0 solved a
Smith-Hazel index `b = P^{-1} G a` from exactly that.

The estimator already holds `Y`, so it can now check itself.

**The rule.** After the engine returns, the per-trait implied heritability

    h2_t = diag(G_hat)_t / stats::var(Y[, t])

is computed on the SAME complete-case rows the fit used, and `h2_t > 1` beyond
the package's relative covariance tolerance (`ng_cov_tol`, the same convention
`ng_multitrait_validate_cov_pair()` uses, so the two checks cannot disagree
about where the boundary is) is a HARD ERROR naming every offending trait, its
genetic variance, its observed phenotypic variance and its implied `h2`.

**Why that phenotypic variance.** The observed sample variance of the fitted
data is the one quantity in scope that is not itself a modelling choice -- no
shrinkage intensity, no PSD projection, no second estimator whose own error
could excuse or manufacture a violation. It is also never STRICTER than the
downstream pair guard on the same data: with `shrinkage = "none"` it is exactly
`diag(P_hat)`, and with `shrinkage = "auto"` the Ledoit-Wolf diagonal is the
`1/n` sample variance, smaller than this unbiased `1/(n-1)` one. Under
`Var(g) = G (x) K` the exact model ratio carries a `mean(diag(K))` factor; for
this package's VanRaden GRM that is ~1 for an outbred panel (measured 0.9955,
range [0.977, 1.016] over 200 panels) and ~2 for a fully inbred DH/RIL panel, so
omitting it makes the guard conservative rather than trigger-happy.

**It blocks for `sommer_remml` too.** Three reasons. The rule is a statement
about the returned matrix, not about the machinery that produced it: an
impossible variance ratio is impossible whichever engine wrote it.
`desired_gain` (Pesek-Baker) solves `b = G^{-1} d` and never sees a `P`, so
`ng_multitrait_validate_cov_pair()` never runs on that route and the estimator's
self-check is the ONLY place an impossible G is caught before it becomes index
weights -- exempting the default engine would leave the default route unguarded.
And measurement supports it rather than contradicting it: over 118 successful
sommer fits from that same `h2 = 0.5` model (120 attempted, 2 engine errors),
exactly one violates (implied `h2 = 1.048`), and inspecting that fit shows it had
assigned `h2 = 0.95` to one trait and `0.065` to another against a truth of 0.5
for both. It is a fit that did not succeed, not a good fit being refused.
Blocking rate: 85.5% for `two_stage_ridge` (171/200 datasets), 0.85% for
`sommer_remml` (1/118). The message says so, and points at
non-convergence, too few records for the number of traits, or a poorly
conditioned GRM rather than at the heuristic's variance scale.

**Nothing is rescaled.** `G_hat` is never clamped, divided through by the
implied `h2`, or projected onto the feasible set. Silently altering a user's
estimate is how a wrong index becomes untraceable; the caller is given the
numbers and chooses.

**`ng_estimate_genetic_correlation()` is new**, because the correlation
structure is the half of `two_stage_ridge` that is estimated acceptably --
its off-diagonals are the Pearson correlation of the per-trait ridge `beta`
vectors, and Pearson correlation is invariant to per-trait multiplicative
shrinkage, so it survives the ridge attenuation that wrecks the diagonal
(relative Frobenius error mean 0.226, max 0.479 over 200 datasets, against mean
1.88 and max 19.3 for the covariance). It returns a correlation matrix, which
has a unit diagonal and therefore makes no claim about genetic variance at all,
so the `h2 <= 1` check has nothing to test and is not applied. It is not a route
back to a refused covariance: multiplying it up by any variance reintroduces
exactly the scale that was refused.

Two new provenance attributes, `phenotypic_variance_observed` and
`implied_h2_vs_observed_variance`, report what the guard judged whether or not
it fired. They are NOT the existing `implied_heritability`, which remains an
independent profile-REML `h2` against the GRM computed only for the
already-shrunk-input warning.

**The one live caller degrades correctly.** The Shiny frontend's
`run_cross_prediction_json.R` calls this estimator with no `method` argument, so it takes
`sommer_remml`, and only as the second rung of its cross-trait covariance ladder when the exact
within-family columns are absent. That call already turns a failure into the frontend's own
refusal ("...it is not computed rather than reported as if the traits were uncorrelated"), so a
refusal there is loud, not silent. Nothing inside this package calls the estimator, so no pipeline
path is affected.

## `ng_load()` sources only what the package actually contains

`ng_load()` globbed `R/[0-9]*.R` off the filesystem. That is a directory
listing, not a package manifest, so on a working tree that carries untracked
files in `R/` it sourced them too -- and it sourced them AFTER the real sources,
so a colliding definition would have won. This repository's working tree carries
eleven such legacy prototype scripts whose names collide with the package's own
numbering (`R/01_relationships.R`, `R/02_duplicate_detection_legacy.R`,
`R/02b_duplicate_detection.R`, `R/03_ld.R`, `R/04_variance_simple_uc.R`,
`R/05_pmv.R`, `R/06_build_cross_data.R`, `R/07_optimizers.R`,
`R/08_select_optimal_parents.R`, `R/09_pmv_cpp.R`, `R/10_marker_effects.R`), so
every source-loaded dev/test session exercised a different file set from the
installed package.

`ng_load()` now excludes exactly what `.Rbuildignore` excludes. Dev-load should
match the BUILT package, and `.Rbuildignore` is the single existing declaration
of what the built package contains, so deriving the load set from it makes the
two the same definition instead of two lists that drift. It needs no external
tooling, unlike a git-tracked-files rule (git is absent from an unpacked
tarball, a CI image or a copied directory) and unlike a manifest in `load.R`
(a third list to keep in sync). Verified: the resulting set is exactly the 54
files a `git archive HEAD` tree contains, and exactly the git-tracked
`R/[0-9]*.R` set.

A file that SHOULD be loaded and cannot be read is now a loud, named failure.
This matters more than it sounds: on OneDrive/iCloud "Files On-Demand" a
cloud-only placeholder OPENS successfully, reports its real size through
`file.size()`, and then reads as ZERO bytes with no error and no warning, so
`sys.source()` on one is a silent no-op. `ng_load()` compares the reported size
against the bytes actually obtained and stops, naming the file and the cause.
Loading also got faster on such a tree -- 16.2 s to 0.20 s here, all of it
placeholder I/O that was buying nothing.

## Verification

`tests/` deep harness: 136 PASS / 1 FAIL over 137 scripts. The one failure is
`alphamate_external`, unchanged and environmental — the untracked
`external/AlphaMate/binaries/AlphaMate.exe` is a OneDrive placeholder and `sh`
returns 126, the same cause as the eleven `R/` files seen through the shell
instead of through R. `tests/pipeline_integration.R` passes (19m59s).
`testthat`: 80 passed / 0 failed / 1 skipped. Statistical release gate: PASS
(44/44) against the installed 0.30.0.

## Bit-identity

`tests/estimator_selfguard_bit_identity.R` compares the 0.30.0 tree against a
pristine 0.29.0 tree (`git archive HEAD` of `7ec06cd`, unpacked to `/tmp`) with
`identical()` at `tolerance = 0`, on a dataset whose implied `h2` is below 1 on
every trait. `G_hat`, its `genetic_correlation`, its `genetic_variance` and
`residual_variance` attributes, `P_hat`, and a full `ng_score_crosses()` run
under the reduced `ng_load()` file set are all bit-for-bit unchanged.

# nextgenCrossDesign 0.29.0

## Blocking validity guards for user-supplied covariance matrices

`phenotypic_covariance` (P) and `genetic_covariance` (G) can arrive from a GUI
over a JSON bridge, from a spreadsheet, or from another package's output. An
invalid one is not a local problem: it propagates into the index weights, the
cross ranking, the robust allocation and the joint probabilities, and every
number downstream looks plausible. All four checks below are therefore HARD
ERRORS, not warnings, and every message names the offending trait(s), gives the
offending number, and states the principle that was violated so a reviewer can
check the objection independently.

**A pair that satisfies every guard is numerically untouched.** Verified against
a pristine 0.28.0 tree at `tolerance = 0` for the `economic_index` and
`desired_gain` scores, the index coefficients, the index SDs and the selected
crossing plan (`tests/covariance_guards_bit_identity.R`).

### Guard 1 -- the symmetry tolerance is now RELATIVE

The symmetry test was ABSOLUTE (`max|M - M'| > 1e-8`) while the positive
semidefinite test six lines below it was already RELATIVE. Covariances carry the
square of the trait's unit, so magnitude is a property of the breeder's unit
choice and not of validity: yield in kg/ha gives covariances of order 1e6, where
agreement to 1e-8 absolute is unreachable for any hand-entered or
spreadsheet-rounded matrix; the same yield in t/ha gives order 1, where it is
trivial. Same principle, two conventions.

Symmetry is now judged against the largest element magnitude -- the right
reference for an ELEMENTWISE comparison, as the spectral radius is for an
EIGENVALUE comparison -- through one shared helper, so the package applies a
single convention: `tol = 1e-8 * max(1, scale)`. The `max(1, .)` floor keeps the
rule absolute for small matrices, which makes this a strict relaxation and never
a tightening: nothing that passed before can fail now.

The message now names both cells and the discrepancy, e.g.

> `genetic_covariance must be symmetric: cov(yield, protein) = 1.7 but cov(protein, yield) = 1.2, a discrepancy of 0.5 against a tolerance of 9e-08`

### Guard 2 -- `P - G` must be positive semidefinite

Since `P = G + R`, the residual covariance `R = P - G` is itself a covariance
matrix and must be PSD. This is the multivariate generalisation of
`0 <= h2 <= 1`: a negative eigenvalue of `P - G` means some linear combination of
the traits is assigned more GENETIC than PHENOTYPIC variance. It is the signature
of P and G having come from different sources, scales or units -- the mistake a
GUI with two upload boxes makes easy.

The check runs whenever both matrices are supplied, AFTER by-name alignment, so
trait order can never confound it or launder an invalid pair.

The elementwise special case `diag(G) <= diag(P)` (per-trait `h2 <= 1`) is tested
first, because it yields a far clearer message:

> `genetic variance exceeds phenotypic variance, which implies a heritability above 1 and is impossible: trait 'protein' has genetic_covariance variance 6 exceeding its phenotypic_covariance variance 5, implying h2 = 1.2`

The general eigenvalue case catches what no per-trait check can see -- every
individual `h2` legal, but a trait CONTRAST with `h2 > 1` -- and reports the
offending eigenvalue together with the contrast that realises it and that
contrast's genetic and phenotypic variance.

**This guard fires on the package's own estimator.** `two_stage_ridge` `G_hat`
paired with `ng_estimate_phenotypic_covariance()` `P_hat` implies `h2 > 1` on 87%
of datasets from the simulation in `tests/genetic_covariance_estimator.R`
(measured over 200 independent draws), because the GBLUP lambda-inversion
diagonal overshoots. Before 0.29.0 that pair was accepted and solved, producing
Smith-Hazel weights from a `P = G + R` decomposition that does not exist.

### Guard 3 -- positive DEFINITE for the matrix the index inverts

`economic_index` (Smith-Hazel) solves `b = P^{-1} G a` and inverts **P**;
`desired_gain` (Pesek-Baker) solves `b = G^{-1} d` and inverts **G**. Only PSD was
checked, which permits a singular matrix, and `ng_multitrait_solve_index()`
carries `ridge = 1e-6` that would quietly produce weights from it.

Measured on a rank-2 `G`: the un-guarded solve returned an index lying **entirely
in G's null space** (`|cos| = 1.000000` to the null eigenvector, `b' G b = 0`),
because LAPACK reports the zero eigenvalue as ~5e-15 -- just above the routine's
own keep-tolerance -- so that direction receives a coefficient of order 1/5e-15
and survives normalisation as the whole index. The coefficients, sign included,
were fixed by floating-point noise, and the index carried no genetic variance at
all. Silently ridging a singular matrix is worse than refusing it.

The inverted matrix must now be positive definite: not numerically singular, and
condition number at most 1e8. That threshold is set by the ridge, not chosen
freely. The solve runs on `M + r I` with `r = ridge * mean(eigenvalues)`, and
using `mean(eigenvalues) >= lambda_max / p`, the default ridge reaches
`lambda_min` once `kappa >= p * 1e6`, i.e. `kappa ~ 1e8` for a realistic
multi-trait index. Below it the ridge is a rounding correction -- the MILD
ill-conditioning it exists to absorb, which is not blocked. Above it the ridge
supplies the answer. (It is also the classical numerical limit:
`1e8 * .Machine$double.eps ~ 2e-8`, so more than half the significant digits in
`b` are already gone.) The threshold deliberately does not scale with a caller's
`ridge`: more shrinkage does not make a rank-deficient matrix more informative.

A merely semidefinite matrix is NOT blocked for a method that does not invert it
-- a singular `G` still runs under `economic_index`, where `G` enters only as the
forward projection `G a`.

### Guard 4 -- implied correlations inside `[-1, 1]`

`|M_ts| <= sqrt(M_tt M_ss)` is implied by PSD (it is the 2x2 principal-minor
condition), so this adds no mathematical content. It exists entirely for the
message, and runs BEFORE the eigenvalue test so the actionable objection wins:

> `genetic_covariance implies a correlation outside [-1, 1] ...: yield vs protein: implied correlation 1.31`

points at the two cells to fix; "must be positive semidefinite" does not. A
correlation of exactly +/-1 is on the boundary and is not refused.

## `R/31_genetic_covariance.R`: per-trait seeds are identity-derived

0.28.0 fixed position-derived per-trait seeds in the runner (`R/39`) and the
multi-trait posterior (`R/32`) but left the same `seed + j` pattern in `R/31`,
on the premise that the file had no non-test caller. That premise was wrong in
one direction and right in another, and both are recorded here.

The Shiny frontend does call `ng_estimate_genetic_covariance()` directly
(`inst/app/tools/run_cross_prediction_json.R`), as the second rung of the
cross-trait covariance ladder behind the exact within-family columns. However,
it calls it with the default `method = "auto"`, which resolves to `sommer_remml`
and never enters the per-trait loop. The defective sites are reachable only via
an explicit `method = "two_stage_ridge"`, or via `ng_posterior_genetic_covariance()`
under `allow_heuristic = TRUE`. Both are exported, documented, user-callable
paths, so they are fixed regardless.

Four sites, treated the same way as 0.28.0 treated the equivalent sites:

* **Lambda-CV fold splits** (`ng_genetic_cov_two_stage_ridge()`, and the
  per-trait fit inside the parametric bootstrap): every trait now shares ONE
  partition, the base `seed` unmodified. Folds are a nuisance parameter of lambda
  selection, not a source of innovation -- given lambda, `beta_j` is a
  deterministic function of `y_j` alone -- so a shared split cannot couple the
  traits, and it makes per-trait `cv_predictive_r2` comparisons paired.
* **Posterior draws** (`ng_posterior_genetic_covariance(method = "beta_posterior")`):
  streams stay DISTINCT per trait, keyed on the trait NAME via
  `ng_trait_rng_seed()`. A shared stream would give identical innovations and
  manufacture cross-trait correlation in exactly the `R_beta_b` being computed.
* **Parametric-bootstrap residual draws**: these were a single `set.seed(seed)`
  outside both loops, with one shared stream consumed in COLUMN ORDER -- so trait
  `j` got whatever remained after traits `1..j-1` had drawn theirs. Each
  `(trait, draw)` now gets its own identity-derived stream, scoped by
  `ng_with_rng_seed()` so the caller's global RNG state is restored.

`tests/genetic_covariance_trait_order_invariance.R` permutes the columns of `Y`,
reorders the answer back by name, and requires: per-trait lambda, CV r2,
`sigma_e2`, `sigma_g2`, the marker effects, the beta correlation matrix and the
raw `G` all identical at `tolerance = 0`; every posterior draw of both methods
invariant. The one quantity that is not bit-identical is `G_hat` after the PSD
projection, which differs by ~9e-16 RELATIVE -- LAPACK reassociation inside
`eigen()`/`nearPD` on a permuted matrix, not seed dependence. It is asserted
below 1e-12 relative and documented rather than hidden.

## `tests/genetic_covariance_estimator.R`: a criterion that passed by luck is withdrawn

The file asserted `||G_hat - G_true||_F / ||G_true||_F <= 0.5` at one hard-coded
seed. Measured over 200 independently simulated datasets from that exact
generating model, the criterion is met by **12%** of them (mean 1.88, sd 2.06,
median 1.26, max 19.3). Seed 2026 passed by chance -- and the `R/31` seed fix
above is itself a reseeding, which moves that seed's ratio from 0.374 to 1.112.
The criterion was not tightened to a new lucky seed.

Diagnosis splits the estimator cleanly in two, as `R/31`'s own comments predict:

* **Off-diagonals** are the Pearson correlation of the per-trait ridge beta
  vectors, which is invariant to per-trait multiplicative shrinkage and so
  survives the ridge attenuation. Relative Frobenius error of the genetic
  CORRELATION matrix over the same 200 datasets: mean 0.226, sd 0.067, max 0.479.
* **Diagonals** are the GBLUP lambda inversion
  `sigma_g2 = sigma_e2 * denom / lambda`, and that is where all the error lives:
  the worst per-trait genetic variance is out by a factor of ~6 on average and by
  up to ~1000. It also overshoots often enough that `G_hat` exceeds `P_hat` on 87%
  of datasets (Guard 2 above).

The file now ASSERTS the correlation-structure property, which genuinely holds,
with a bound (0.6) justified by that 200-dataset measurement and checked over 8
independent datasets rather than one; CHARACTERISES the covariance-scale error by
printing it under a collapse-only ceiling; and asserts that the 0.29.0 `P - G`
guard refuses the `(G_hat, P_hat)` pair. The Smith-Hazel integration check is
re-pointed at the generating model's own `G_true` / `P_true`, since what it tests
is that a valid pair reaches the Smith-Hazel solve, not the estimator's accuracy.

Users should read this as: `two_stage_ridge` is usable for the SHAPE of the
genetic correlation structure and not for the SCALE of genetic variance, which is
consistent with the warning it already emits. It remains an opt-in heuristic;
`sommer_remml` remains the default and the only formal variance-component route.

## Test-suite consequences

`ng_fit_ridge_effects_posterior()` uses ONE `seed` for two jobs that need
opposite treatment -- the lambda-CV fold split and the posterior innovations.
`ng_posterior_genetic_covariance(method = "beta_posterior")` therefore now
selects lambda up front on the shared partition and passes it in explicitly, so
the identity-derived seed governs only the draws. This preserves that function's
documented invariant that its `G_mean` diagonals equal the point
`ng_estimate_genetic_covariance()` diagonals exactly (measured gap 1.65e-15);
without the separation the two routines chose different lambdas and the gap was
1.31.

Three test fixtures were corrected rather than the guards weakened:

* `tests/posterior_multitrait_rng_scoping.R` and
  `tests/threshold_probability_multitrait.R` passed an arbitrary identity `G`
  (unit genetic variance) against phenotypes of variance ~0.2 and ~0.95, i.e.
  implied `h2` of 4.65 / 7.60 and 1.083 / 1.026.
  `ng_posterior_multitrait_cross_predict()` pairs a caller's `G` with a `P`
  estimated from `Y`, so the new `P - G` guard refuses both -- correctly; neither
  was ever a possible `P = G + R`, and neither test's property depends on `G`'s
  scale. Both now use `G = 0.5 * P_hat` on the rows the run uses.
* `tests/posterior_genetic_covariance.R` asserted per-pair off-diagonal
  attenuation at a flat `1e-3` slack. Attenuation is an expectation statement, so
  a pair whose point correlation is already ~0 has nothing to attenuate and the
  posterior mean's Monte Carlo error moves it either way (measured range
  `[-0.004, +0.010]` over 10 posterior seeds, against `-0.26 .. -0.38` for the
  materially-correlated pair). It now asserts the aggregate attenuation plus a
  per-pair bound of 0.05, ~5x the largest Monte Carlo excess observed.

Two failures found in the same sweep were PRE-EXISTING, verified against a
pristine 0.28.0 tree, and fixed here because both are one-line stale fixtures:

* `tests/posterior_multitrait_usefulness_direction.R` reconstructs R/32's
  per-trait posterior draw with `seed = SEED + 1L`. 0.28.0 changed R/32's key to
  `ng_trait_rng_seed(seed, trait)` and left this behind, so the file has been
  comparing against a different draw stream -- and failing -- since 0.28.0,
  independently of the direction property it exists to test.
* `tests/metric_name_normalization.R` builds a minimal `cfg` for
  `ng_cp__build_ctx()`, which validates a COMPLETE config (defaults live in
  `ng_run_cross_prediction()`'s formals, not the builder). A later release added
  a `ci_level` validation the fixture never supplied, so it died before reaching
  a single metric assertion.

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
