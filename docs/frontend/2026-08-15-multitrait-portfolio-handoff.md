# Multi-trait portfolio & risk — frontend handoff

**Date:** 2026-08-15
**Backend:** `nextgenCrossDesign` 0.19.0 (`R/37_cross_portfolio_risk.R`, `R/39_cross_prediction_runner.R`)
**Tests:** `tests/priority_risk_portfolio_multitrait.R`, `tests/priority_risk_portfolio.R`
**Design record:** Rev 7 of `2026-07-25-cross-priority-risk-portfolio-design.md`
**Shareable version:** https://claude.ai/code/artifact/197d1ac7-abab-47b7-a9ad-03a4480039f5

## Verdict on "multi-trait should be the same as single-trait"

Right about the **panel**, wrong about the **data**.

- The reporting surface is identical: same six columns, names, types, factor levels, quadrant
  labels and tertile risk bins. Nothing in the display layer needs the trait count.
- The backend previously emitted none of it. `R/39` annotated only when
  `length(trait_spec$column) == 1L`; a two-trait run returned zero portfolio columns and a
  `NULL` diagnostics block. That gate is removed in 0.19.0.

## Contract

Identical column set on `selected_crosses` and `candidate_crosses`.

| Column | Type | Meaning |
| --- | --- | --- |
| `cross_level` | numeric | X axis. Index of mid-parent GEBVs. Unitless. |
| `cross_upside` | numeric | Y axis. `sqrt(w'Sw)`, index SD within the family. ≥ 0. |
| `cross_confidence` | numeric | 0–1, higher = better estimated. Within-run only. |
| `risk_bin` | ordered factor | `low` < `med` < `high`. Tertiles of the plan. |
| `confidence_method` | character | `midparent_pev_index[_partial]` \| `reliability`. |
| `portfolio_profile` | factor | `breakthrough` / `workhorse` / `long_shot` / `deprioritize`. |
| `portfolio_basis` | character | **New.** `single_trait` \| `linear_index` \| `linearized_rank_index`. |
| `risk_driver_trait` | character | **New, multi-trait.** Trait contributing most of this cross's index PEV. |
| `risk_driver_share` | numeric | **New, multi-trait.** That trait's share, 0–1. |
| `cross_level_rule` | character | Provenance: method, basis, per-trait weights, units. |

`priority_risk_diagnostics` gains, for multi-trait runs: `basis = "multi_trait_index"`,
`n_traits`, `index_method`, `index_weights`, `upside_method`, `portfolio_basis`,
`portfolio_basis_note`, `index_traits` (per trait: direction, weight,
`marker_effect_reliability`, `mean_variance_share`, `mean_pev_share`),
`risk_disproportionate_traits`, `pev_concentration`, `pev_concentration_note`.

Supporting per-trait columns ride both tables: `<trait>_mean_gebv`, `<trait>_vpm`,
`<trait>_midparent_pev`, `wf_var_<t>`, `wf_cov_<t>_<s>`.

## The three frontend changes

1. **Badge `linearized_rank_index` — not optional.** `auto`/`weighted`/`threshold` combine
   rank-normalized traits, so no linear index exists in genetic units; the axes come from
   reinterpreting trait weights as standardized-unit coefficients. The quadrant is internally
   consistent but is **not** a decomposition of `multi_trait_score`. `auto` promotes to
   `weighted` whenever trait weights are present, so this is the common case.
   `portfolio_basis_note` carries breeder-facing copy including the remedy (supply
   `economic_weight` → `linear_index`, no badge).
2. **Axis labels lose their units.** Index level is unitless; label by role and put
   `index_weights` in a drawer. Minimized traits carry **negative** weights — render as
   direction, not a bare minus sign.
3. **Confidence tooltip copy.** `midparent_pev_index[_partial]`; the index PEV is the
   block-diagonal `Σ w²·PEV_k` (traits fitted by independent univariate ridges, so no
   cross-trait estimation-error covariance exists). Ranking input, not a calibrated interval.

## Risk attribution — no weakest-link floor, by design

`cross_confidence` is **not** floored by the worst-estimated component trait.
`Σ w_k²·PEV_k` is already the correct propagation of estimation error into the index: a lightly
weighted trait contributes little because it genuinely moves the index little, and a lightly
weighted trait with a large PEV still dominates the sum unaided (verified — a pure-noise trait at
the *smallest* weight took 47% of the index PEV and was the risk driver for every cross). A floor
would double-count uncertainty the weighting already handles and would make `risk_bin` stop
meaning "how well is this cross's *index* estimated". The related breeder question — "is this
cross unacceptable on a minor trait?" — is a **threshold** question, answered separately and
better by `trait_checks` (R/44) and multi-trait `min_value`/`max_value`.

What was missing was attribution, which decomposes exactly:

- **`risk_driver_trait` / `risk_driver_share`** per cross — `w_k²·PEV_k / Σ`, shares sum to 1.
- **`index_traits`** per trait, averaged over the candidate pool — weight, reliability, share of
  index *spread* (`w_k(Sw)_k / w'Sw`, sums to 1, **may be negative** when a trait is antagonistic
  and genuinely removes spread), share of index *error*.
- **`risk_disproportionate_traits`** — above an equal share of the error *and* more than twice as
  much of the error as of the spread. Traits buying more uncertainty than opportunity.

### The commensurability caveat (surface `pev_concentration_note`)

Per-trait PEVs are only comparable across traits when their **residual variances** are, and
`σ²ₑ` is estimated in-sample. A trait whose ridge λ lands on the grid floor interpolates its
training rows, collapsing `σ̂²ₑ` and hence `Σ_β = σ̂²ₑ(X'X+λI)⁻¹`, and reports a near-zero PEV
regardless of how well it is actually predicted. Observed in a small run: λ = 0.695 for one trait
and 0.01 for two others put **100%** of the index PEV on the first. The arithmetic is exact; the
inputs were not comparable. When one trait carries >90% of the index PEV the run now says so via
`pev_concentration_note`, and `risk_bin` must not be read as an index-wide statement.

## UI rule to enforce

`cross_confidence` and `risk_bin` are **within-run** quantities — a min–max normalization and
tertiles of the crosses on screen. A plan always contains roughly one third "high risk"
regardless of whether it is well or badly estimated. Never compare across runs; never phrase a
high-risk count as an absolute quality statement.

## Verification matrix

| Run | Expect |
| --- | --- |
| 1 trait | `basis = "single_trait"`, `portfolio_basis = "single_trait"`, no badge, trait-unit axes |
| 2 traits, `weight` | `index_method = "weighted"`, `portfolio_basis = "linearized_rank_index"`, badge shown |
| 2 traits, `economic_weight` | `index_method = "economic_index"`, `portfolio_basis = "linear_index"`, note `NA` |
| any `decrease` trait | negative entry in `index_weights`, drawer reads "lowers index" |
| 3+ traits | `n_traits` correct, all pairwise `wf_cov_*` present |

## Notes for whoever maintains this next

- **Shared basis.** `cross_level`/`cross_upside` come from one raw-unit `w`, resolved once on the
  candidate pool (`ng_multitrait_index_reference`). Deriving it per table put the plan and the
  candidate pool on different index axes — numbers looked fine in isolation, views were not
  comparable. Guarded by a regression assertion in the multi-trait test.
- **Exact covariance, not the diagonal.** On an antagonistic trait pair the diagonal shortcut
  overstated index SD by 1.85× at the median, up to 7×. `S` is additionally rescaled so its
  diagonal reproduces the run's reported `vpm`, which preserves `T == 1 → sqrt(vpm)` continuity
  and honours het-parent corrected variances.
- **Cost** is linear in T (measured 14.0 / 27.4 / 40.9 s for T = 1/2/3 at 4,005 pairs × 1,500
  markers). The covariance call is chunked over the pair dimension in `R/39` because
  `ng_cross_trait_within_family_cov` materializes `n_pairs × m` plus one `m × n_pairs` per trait.
