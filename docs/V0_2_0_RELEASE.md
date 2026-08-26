# `nextgenCrossDesign` 0.2.0 — Posterior-aware cross prediction

This release moves the package from correctness parity with
PopVar / SimpleMating / AlphaMate / genomicMateSelectR (the goal of v0.1.0)
to a **posterior-aware** cross-prediction layer that none of those tools
currently ship. The aim is to give breeders the uncertainty information
they actually need to make crossing decisions, not just point estimates.

The four additions in this release are independent and each addresses a
gap in the benchmark landscape:

## D1 — Posterior cross prediction (`R/30_posterior_prediction.R`)

`ng_fit_ridge_effects_posterior(geno, y, n_draws, method)` produces S draws
of the marker effects β from the ridge posterior. Two modes:

- `method = "closed_form"` (default, fast): uses the Bhattacharya–Chakraborty–
  Mallick (2016) sampler, factorising `A = XX' + λI` once and producing each
  draw in `O(n²)` back-substitution plus `O(np)` multiply. Conditional on
  (σ²_e, λ) — i.e., empirical Bayes at the cross-validated λ.
- `method = "mcmc"` (slower, more honest): Gibbs over (β, σ²_e, σ²_β) with
  vague inverse-gamma hyperpriors. Integrates over hyperparameter uncertainty.

`ng_posterior_cross_predict(geno, posterior_effects, ...)` then runs the
existing cross-scoring kernel once per draw and returns, alongside every
column from `ng_score_crosses()`:

- `uc_dh_gebv_post_mean` / `_lower` / `_upper` — credible interval on
  usefulness in genetic-value units.
- `dh_pmv_var_post_mean` / `_lower` / `_upper` — credible interval on the
  predicted within-family variance.
- `p_superior_progeny_post_mean` / `_lower` / `_upper` — see D4.
- `posterior_topn_prob_N` for each `N` in `top_n_targets` — fraction of
  posterior draws in which the cross sits in the top-N by usefulness. This
  is a **rank-stability** metric that no benchmark package reports.

Performance: ~1–2 minutes for 80 parents × 5K markers × 500 draws under
Haldane DH with the C++ kernel; Kosambi / RIL go through the dense path and
are O(M²) per draw — the function emits a warning and recommends reducing
`n_draws` or LD-pruning.

## D2 — Full off-diagonal posterior PMV (`R/03_metrics.R`)

The diagonal-only PMV correction `Σ_k d_k² · Var(β_k)` was correct given a
diagonal posterior `Cov(β)`. Marker effects under LD are not independent in
the posterior, though, and ridge produces a fully populated `Σ_β = σ²_e
(X'X + λI)⁻¹`. v0.2.0 adds the genomicMateSelectR formulation natively at
package speed:

```
PMV = a' R a + d' (R ⊙ Σ_β) d
```

via `ng_fit_ridge_effects(return_beta_cov_full = TRUE)` and
`ng_score_crosses(posterior_cov_full = Σ_β)`. The new column
`dh_pmv_var_full_posterior` lands alongside the legacy diagonal `dh_pmv_var`
so legacy pipelines are not invalidated.

Cost: ~80 s for 5K markers × 80 parents on one core (dense `O(M²)` per pair
plus one m×m Hadamard). Verified against the 2-locus closed form to 1.1e-16
in `tests/posterior_pmv_full.R`.

## D3 — Built-in genetic-covariance estimator (`R/31_genetic_covariance.R`)

v0.1.0 added `phenotypic_covariance` and `genetic_covariance` parameters to
the multi-trait index but required the user to supply them. v0.2.0 closes
that gap:

- `ng_estimate_genetic_covariance(geno, Y, method = c("auto", "two_stage_ridge", "sommer_remml"))`
  - `auto` uses `sommer::mmer()` REML if `sommer` is installed; falls back to
    `two_stage_ridge` otherwise.
  - `two_stage_ridge`: per-trait ridge gives β̂; diagonals via the GBLUP
    lambda-inversion `σ²_g = σ²_e · Σ 2p(1-p) / λ`; off-diagonals via the
    correlation of trait β̂ vectors (invariant to per-trait shrinkage).
    Composed as `G = D R_β D` and projected to the nearest PSD matrix
    (eigen-clip).
- `ng_estimate_phenotypic_covariance(Y, shrinkage = c("none", "auto"))` —
  sample covariance with optional Ledoit–Wolf shrinkage toward a diagonal
  target (analytical shrinkage intensity).

Recovers a 3 × 3 G with Frobenius error 0.385 × ‖G_true‖ on the standard
synthetic fixture (n = 100, m = 200).

These plug directly into the v0.1.0 multi-trait solver:

```r
G_hat <- ng_estimate_genetic_covariance(geno, Y)
P_hat <- ng_estimate_phenotypic_covariance(Y, shrinkage = "auto")
ng_add_multitrait_score(scores, traits, method = "economic_index",
                        phenotypic_covariance = P_hat,
                        genetic_covariance = G_hat)
# meta$economic_index_cov_source == "smith_hazel"
```

This is the first end-to-end Smith–Hazel pipeline in the package and is
the layer that lets the cross-prediction outputs feed a real selection
index instead of a phenotypic-proxy one.

## D4 — P(superior progeny ≥ τ)

The breeder-actionable question is "what's the chance this cross gives me
at least one line above yield X?", not "what's the expected mean of the top
decile?". v0.2.0 adds the closed-form metric:

```
P(max ≥ τ) = 1 − Φ((τ − μ) / σ)^k
```

Available as:
- `ng_p_superior_progeny(mu, sigma, tau, k_progeny)` — vectorised helper.
- `ng_add_p_superior_progeny(scores, tau_superior, k_progeny, ...)` — adds
  `p_superior_progeny` to an existing score table.
- Automatic in `ng_posterior_cross_predict(tau_superior, k_progeny, ...)`
  with full credible interval.

The implementation uses log-space `pnorm` to avoid the underflow that hits
naïve `Φ(z)^k` once k > ~50 and z is moderate.

## D1c — Posterior-aware OCS

`ng_optimize_robust_mating_plan()` is a new entry point that calls the
existing MIP optimiser but with a robustness-aware gain column:

- `objective = "posterior_quantile"` (default) maximises the
  `robustness_quantile` quantile of posterior usefulness. The `NULL` default
  uses the exact empirical lower credible bound already cached from posterior
  draws (2.5th percentile for the default 95% interval). A different quantile
  must be cached by choosing the matching `ci_level` during posterior prediction;
  normal reconstruction is available only via explicit
  `allow_normal_approximation = TRUE` and is flagged in the plan summary.
  This selects crosses that perform well even at the pessimistic end of the posterior.
- `objective = "posterior_topn_prob"` maximises Σ `posterior_topn_prob_N`
  over selected crosses. This maximises the expected number of selected crosses
  in the posterior top-N set; it is not a joint probability that the entire plan
  is top-N stable.

Both reuse the gain-vs-diversity penalty machinery from
`ng_optimize_mating_plan()`, so coancestry, parent-use, kinship, and family-
size constraints all still apply.

## Defaults that changed

None in the original 0.2.0 release. Version 0.22.0 later corrected the robust
allocator default to the exact cached empirical lower credible bound and made
normal approximation explicitly opt-in. See `V0_22_0_RELEASE_CANDIDATE.md` for
the current interface and migration guidance.

## New exports

```
ng_fit_ridge_effects_posterior
ng_sample_ridge_posterior_bcm
ng_sample_ridge_posterior_mcmc
ng_posterior_cross_predict
ng_p_superior_progeny
ng_add_p_superior_progeny
ng_optimize_robust_mating_plan
ng_estimate_genetic_covariance
ng_estimate_phenotypic_covariance
```

## Late-stage v0.2.x additions

After the initial v0.2.0 cut three follow-ups were shipped on the same
branch — each tagged with the headline differentiator they unlock.

### Posterior G (`ng_posterior_genetic_covariance`)

The point G estimator landed in v0.2.0 but D1's posterior framework had
no companion posterior over G to feed back into the Smith-Hazel /
Pesek-Baker solver. The follow-up adds a t × t × n_draws posterior
array with two modes:

- `method = "beta_posterior"` (default; fast, coherent with D1). Per
  trait, fit ridge once and draw S samples of β_t from BCM. For draw
  s, build `G_s = D · R_β^s · D` where `D = sqrt(σ²_g)` is fixed across
  draws (lambda-inversion is hyperparameter-conditional) and `R_β^s`
  recomputed per draw. Off-diagonals are biased toward zero relative to
  the point estimate (Pearson attenuation under added noise) and the
  posterior captures the *spread* around that shrunken center.
- `method = "parametric_bootstrap"` (slower; covers hyperparameter
  uncertainty). Resamples residuals around each per-trait fit and
  reruns the full two-stage estimator with CV-tuned λ.

`attr(G_draws, "G_mean")` plugs straight into
`ng_add_multitrait_score(..., genetic_covariance = ...)` at
`cov_source = "smith_hazel"`.

### Banded Kosambi / RIL kernel

Dense `ng_dh_recomb_variance_pairs_dense()` is O(M²) per pair, which
dominated the per-draw cost in `ng_posterior_cross_predict()` for
non-Haldane-DH targets. The follow-up adds:

- `ng_recomb_decay_banded(marker_map, model, target, window_cm)` —
  base-R COO upper-triangle sparse representation.
- `ng_dh_recomb_variance_pairs_banded()` — consumes the COO band,
  doubles off-diagonal contributions for symmetry, returns the same
  `data.frame(dh_recomb_var, dh_pmv_var)` as the dense path.
- Internal dispatch heuristic in `ng_dh_recomb_variance_pairs()`:
  banded fires when `window_cm` is finite, > 0, and
  `window_cm / mean(per_chr_length_cm) < 0.5`. Override via
  `NGCD_BANDED_RATIO` env var.
- `ng_posterior_cross_predict()` warning is silent when the banded
  path is engaged.

Numerical equivalence to dense at 1.11e-16 (machine precision).
Speedup: ~2.7× at m=1000 / window=10; ~5.6× at m=3000.

### End-to-end multi-trait posterior cross prediction (v0.3.0 preview)

`ng_posterior_multitrait_cross_predict()` joins per-trait β posteriors
(D1) + posterior G (above) + the v0.1.0 multi-trait index in a single
call. For each posterior draw s, the function builds per-cross
mid-parent GEBVs from β_t^s for each trait, plugs them into the index
solver with `G^s` as the genetic_covariance argument, and reads off
the per-cross `multi_trait_score`. Output columns:

```
multi_trait_score_post_mean / _lower / _upper
multitrait_posterior_topn_prob_N
```

Two value modes — `"mean"` (fast; cross-level GEBVs only) and
`"usefulness"` (slow; PMV-aware μ + i·σ per trait per draw). The
default is `"mean"`, which is seconds for n=80, m=5000, t=3, S=200.

This is the layer that completes the cross-prediction × multi-trait ×
posterior cube: per-cross posterior CIs on the multi-trait selection
index plus rank-stability probabilities. No benchmark package today
reports any of these for the multi-trait case.

### Late-stage exports

```
ng_posterior_genetic_covariance
ng_posterior_multitrait_cross_predict
```

## New tests

`tests/posterior_cross_predict.R`, `tests/posterior_pmv_full.R`,
`tests/genetic_covariance_estimator.R`. All 36 tests in the suite pass on
v0.2.0.

## Re-validation status vs v0.1.0

`VALIDATED_STATE.md`'s "re-validation required" banner from v0.1.0 still
applies. v0.2.0 only adds posterior outputs alongside existing columns —
it does not change any v0.1.0 numerical contract. The DH/RIL frontier
policy and the multi-trait recovery claims are still subject to the
v0.1.0 re-validation work.

## Why this is the practical novel step

Each of the four benchmark packages stops at point estimates:

| Tool | Cross mean | Cross variance | Usefulness CI | P(superior progeny) | Posterior rank stability | Full off-diag posterior |
|---|---|---|---|---|---|---|
| PopVar | yes | yes | — | — | — | — |
| SimpleMating | yes | yes | — | — | — | — |
| AlphaMate | (input) | — | — | — | — | — |
| genomicMateSelectR | yes | yes | — | — | — | yes (slow) |
| **nextgenCrossDesign 0.2.0** | yes | yes | **yes** | **yes** | **yes** | **yes (fast)** |

A breeder choosing between two crosses with predicted usefulness 12.3 and
12.1 can now ask: "what's the 95 % credible interval on each, and how
often does each end up in the realised top-10?" — and get an answer
without leaving the package. That is the practical differentiator.
