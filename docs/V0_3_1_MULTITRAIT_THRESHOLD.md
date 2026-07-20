# `nextgenCrossDesign` 0.3.1 — Multivariate threshold probability `Pr(G_c ∈ A)`

This release adds the multi-trait extension of `ng_p_superior_progeny()`.
A breeder can now ask: "what is the probability that at least one of k
progeny falls inside a rectangular target region A = ∏_t [τ_t^lo,
τ_t^hi] across all traits at once?" — without ever leaving the package.

## New exports

```
ng_p_superior_progeny_multitrait
ng_build_cross_trait_covariance
ng_add_p_superior_progeny_multitrait
```

`ng_posterior_multitrait_cross_predict()` gains four optional arguments
(`tau_lower_vec`, `tau_upper_vec`, `threshold_k_progeny`,
`threshold_G_hat`) that add posterior credible intervals on the
multi-trait threshold probability via three new columns:
`p_superior_progeny_mt_post_mean / _lower / _upper`.

## Math

For one progeny, `Pr(G ∈ A) = pmvnorm(tau_lower, tau_upper, mu, Sigma_c)`
via `mvtnorm::pmvnorm` under the standard multivariate-normal progeny
assumption (Allier et al. 2019, Wolfe et al. 2021).

For the order-statistic across k progeny:
`Pr(max ∈ A) = 1 − (1 − Pr(one))^k`

Implemented in log-space (`log1p(-p_one)` → `exp(k · log_complement)`)
to mirror the existing `ng_p_superior_progeny()` underflow guard.

## Σ_c construction

`ng_score_crosses()` produces per-trait `dh_pmv_var` (the diagonals of
Σ_c) but not the trait-trait coupling at the cross level. v0.3.1 adds
`ng_build_cross_trait_covariance(per_trait_var, G_hat = NULL)`:

- `G_hat = NULL` ⇒ diagonal Σ_c (independent traits at the cross level).
- `G_hat = <t × t matrix>` ⇒ Σ_c = D · cov2cor(G_hat) · D with
  eigen-clip PSD projection. Pass the output of
  `ng_estimate_genetic_covariance()` directly.

Negative entries in `per_trait_var` are rejected with a clear error.

## New dependency

`mvtnorm` added to `Suggests:`. Clear error if not installed when the
multi-trait threshold path is requested.

## Posterior path runtime note

Supplying `tau_lower_vec` / `tau_upper_vec` to
`ng_posterior_multitrait_cross_predict()` forces per-trait PMV
computation even when `value_mode = "mean"` (the cheaper mid-parent
path). The function emits a `warning()` so the user is aware the
runtime now matches `value_mode = "usefulness"`. There is no
correctness fallback — Σ_c requires per-trait within-family variance.

## Defaults unchanged

Every existing column, function default, and external-tool baseline is
preserved. All new output columns are added alongside existing ones;
opting in requires calling the new helpers.

## Validation

`tests/threshold_probability_multitrait.R` covers:
- Univariate reductions (trait 2 unconstrained reduces to single-trait
  case; upper = +Inf reduces to existing `ng_p_superior_progeny`).
- Σ_c constructor diagonal and G-scaled modes, including PSD projection
  for a slightly non-PSD G_hat input and negative-variance rejection.
- Score-table augmentation against a 3-cross / 2-trait fixture with
  missing-column and length-mismatch guards.
- Posterior CI integration via `ng_posterior_multitrait_cross_predict()`,
  with non-trivial CI width assertion to catch saturation regressions.
- Empirical agreement vs simulated progeny within 0.02 at k=1 (n=5000
  progeny) and within 0.05 at k=20 (1000 simulated families).

## Spec source

This release implements v0.3.1 of the CPP adoption roadmap
(`docs/superpowers/specs/2026-05-24-cpp-spec-adoption-roadmap-design.md`).
The math `Pr(G_c ∈ A)` is from CPP spec §4 (threshold utility); the
"target-attainment utility" framing is from CPP spec §12.
