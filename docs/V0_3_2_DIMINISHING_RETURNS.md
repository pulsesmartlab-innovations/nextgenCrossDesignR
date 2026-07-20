# `nextgenCrossDesign` 0.3.2 — Diminishing-returns + Ne-aware K selection

This release helps breeders pick K (the number of crosses to make)
**honestly** and **sustainably**. The existing optimizer always fills
exactly the K the user specifies; there is no built-in detection of
whether K=8 vs K=12 actually adds value, and no Ne-aware sanity check
on whether the chosen K keeps the program genetically healthy.

v0.3.2 adds a pure wrapper that iterates the optimizer over a K range,
surfaces the marginal-utility curve plus effective-population-size, and
picks the recommended K under any of four criteria.

## New exports

```
ng_optimize_mating_plan_curve
ng_plot_diminishing_returns
```

## Usage

```r
scores <- ng_score_crosses(...)  # or built any other way
curve  <- ng_optimize_mating_plan_curve(
  scores = scores, K_range = 3:20,
  parent_K = K, lambda_group = 0.05,
  max_crosses_per_parent = 4L,
  criterion = "ne_target",
  ne_min = 50
)
print(curve)              # K, total_gain, group_coancestry, Ne_estimate, ...
attr(curve, "elbow_K")    # recommended K
ng_plot_diminishing_returns(curve)
```

## Criteria

Four ways to pick K, all sharing the same return shape so they're
interchangeable:

- **`"elbow_relative"`** (default): smallest K where the marginal gain
  drops below `relative_threshold` (default 0.05) of the first marginal
  gain in the sweep. Simple, breeder-friendly, scale-stable across runs.
- **`"elbow_kneedle"`**: Satopää et al. (2011) normalized-difference
  algorithm. Robust on smooth concave curves; returns `NA` when the
  curve is too linear or the maximum lies at the boundary.
- **`"ne_target"`**: smallest K with `Ne_estimate ≥ ne_min` (default
  `ne_min = 30`). Sustainability-aware: keeps the program above an
  effective-population-size floor. `Ne_estimate = 1 / (2 *
  group_coancestry)` (Falconer-Mackay).
- **`"coancestry_budget"`**: largest K with `group_coancestry ≤
  coancestry_max` (default 0.05). Hard ceiling on relatedness pressure.

## Why Ne matters

Elbow detection alone watches `total_gain` flatten and says nothing
about whether the breeding program is still genetically healthy at the
chosen K. A plan that maximises one-cycle gain but collapses effective
population size is a long-term loss for the program. The `ne_target`
and `coancestry_budget` criteria let the breeder set the
sustainability floor directly.

The Ne approximation is the standard Falconer-Mackay
`Ne ≈ 1 / (2 · Δf)` with `Δf = group_coancestry`. It's
contribution-weighted via the existing `parent_K` matrix that the
optimiser already uses, so the calculation is consistent with the
allocator's diversity penalty.

## Pure wrapper

`ng_optimize_mating_plan()` is unchanged. The curve helper iterates it
once per K and extracts `attr(plan, "summary")$total_gain` (and a few
other diagnostics) into a single tibble. All optimizer arguments —
`parent_K`, `max_crosses_per_parent`, `lambda_group`, `lambda_mating`,
`gain_col`, `method`, etc. — pass through unchanged via `...`.

## Curve columns

| Column | Meaning |
|---|---|
| `K` | number of crosses in the plan |
| `total_gain` | sum of `gain_col` across the K selected crosses |
| `mean_gain` | `total_gain / K` |
| `group_coancestry` | from the optimizer's summary; relatedness pressure |
| `unique_parents` | distinct parents used |
| `Ne_estimate` | `1 / (2 * group_coancestry)`; `Inf` when coancestry is 0 |
| `marginal_gain` | `total_gain[K] − total_gain[K-1]` |
| `relative_marginal` | `marginal_gain / marginal_gain[K_range[1] + 1]` |

## Validation

`tests/portfolio_size_sweep.R` covers:

- Curve shape: monotone non-decreasing `total_gain`; non-increasing
  `marginal_gain` on a deterministic fixture.
- Elbow attribute attached for every criterion; correct selection
  semantics (smallest K satisfying for `ne_target`; largest K
  satisfying for `coancestry_budget`).
- Dot args propagate (infeasible K signals an error).
- Kneedle detects the elbow at K=8 (±1) on a synthetic log-then-plateau
  curve; returns NA on a purely linear curve.
- `Ne_estimate` column matches the analytical formula on every finite
  row.
- Positive-path tests for both `elbow_relative` (Test 2b) and
  `ne_target` (Test 8b) — exercising the picker actually returning an
  integer K (not just NA-with-guard).
- Legacy `elbow_method` argument signals an error (rename guard).
- `ng_plot_diminishing_returns()` returns a `ggplot` with the elbow
  vline; toggle drops it; graceful handling of NA elbow_K.
- Monotonicity on a realistic 20-parent fixture with `lambda_group =
  0.5` (validation gate 1): at most one strict increase in the
  marginal_gain across the sweep.

## Deferred to v0.3.6

Breeder review during this release also surfaced the need for
**priority parents with counterfactual merit evaluation** — the
ability to specify "use CB29 and CB47 at least once" with the system
computing the cost of obeying so the breeder can decide whether to
insist. This is deferred to v0.3.6 because it touches the optimizer
core (new MIP constraints) and depends on v0.3.5's objective
decomposition for the cost block. See the roadmap spec at
`docs/superpowers/specs/2026-05-24-cpp-spec-adoption-roadmap-design.md`
for the v0.3.6 design.

## Defaults unchanged

Every existing function default is preserved. The new helpers are pure
additions; opting in requires calling them explicitly.

## Spec source

This release implements v0.3.2 of the CPP adoption roadmap, expanded
mid-release with the Ne-aware criteria after breeder review on
2026-05-26. The original spec's "diminishing-returns detection" novelty
path is from CPP spec §12 (permitted-novelty list); the Ne expansion
addresses CPP spec §3 (diversity lesson) and §11 (variance prediction
risk guardrail).
