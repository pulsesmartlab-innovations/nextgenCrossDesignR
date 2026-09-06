# Check lines as references, not filters

**Date:** 2026-09-03
**Status:** Design approved, not yet implemented
**Backend:** nextgenCrossDesign 0.22.0 -> 0.23.0 (breaking)
**Frontend:** NextGenCrossDesign 0.26.0 -> 0.27.0
**Supersedes:** `docs/design/2026-07-25-trait-check-threshold-design.md` (the veto)

## 1. The breeder's question

A **check** is a standard genotype whose performance the breeder already knows: a released
variety, a commercial check, the entry every trial is benchmarked against. It is a yardstick,
not breeding material for the crossing block.

Per trait, the question is:

> Is the mean of the two parents of this cross on the good side of the check?

`increase` traits want the mid-parent **above** the check; `decrease` traits want it **below**.
The purpose is defensive: stop a breeder from relying on a cross whose progeny population mean
lands on the wrong side of a variety they already have.

Single-trait: one check. Multi-trait: one check per trait, possibly a different genotype for
each.

## 2. Why a check can never be a row in the cross table

To appear as a row, a check would have to be crossed. The only cross representing the check
itself is check x check -- a **self**. `ng_make_pairs()` (`R/00_utils.R:18`) enumerates pairs
with `include_self = FALSE`, and self-crosses are not part of this design.

Two further consequences:

- A check has **no variance**. It is a point, not a distribution. Nothing in this design may
  require a variance for it.
- A check placed in the candidate parent matrix becomes a **mating candidate**: the optimizer
  can allocate crosses to it and it consumes parent-use budget. Neither is what a check is for.

Therefore a check is only ever a **scalar reference value per trait**, carried alongside the
results and drawn on top of them.

## 3. What shipped before, and why it is being replaced

Three mechanisms in this codebase use the word "threshold". Only the third is what breeders
asked for:

| Mechanism | Location | Effect on the run |
|---|---|---|
| `min_value` / `max_value` + `threshold_policy` | `R/19_multi_trait_selection.R:222` | Penalises or zeroes crosses **in the objective** |
| `trait_checks` veto | `R/44_trait_checks.R` | **Filters rows out** before the optimizer |
| This design | -- | **Nothing.** Reference only |

The veto (backend 0.14.0, frontend 0.16.0) computes the right comparison but draws the wrong
conclusion from it, and because it reads the check's value out of the parent value vector it
hard-requires the check to be a candidate parent (`R/39_cross_prediction_runner.R:1215`) -- the
exact thing a check is not.

**Decision: replace it.** `exclude_threshold_violators` and `check_basis` are removed. The
direction-resolution logic in `ng_trait_check_spec()` and the mid-parent comparison are kept
and reused. The veto was never in a tagged release and the deployable image is on an older
backend, so the migration cost is near zero.

## 4. Decisions

| # | Decision | Rationale |
|---|---|---|
| D1 | Checks arrive in a **separate check file**, never mixed with parents | `rownames(check_geno)` cannot reach `ng_make_pairs()`; isolation is structural, not a flag |
| D2 | Reference only -- nothing is filtered, penalised, or reordered | A check informs the breeder; it does not decide for them |
| D3 | Direction defaults from the trait's breeding direction, per-trait override | `increase` -> flag if below; `decrease` -> flag if above |
| D4 | Check value follows the run's **`mean_source`**, per trait | Guarantees apples-to-apples with the plotted axis by construction |
| D5 | `check_basis` is deleted | A user-settable basis is precisely the knob that puts the check on the wrong scale |
| D6 | Report **P(beat check)** alongside the mean comparison | A below-check cross with high variance can still throw a superior tail |
| D7 | Checks never enter QC, LD, GRM, or effect estimation | Otherwise adding a check would move every cross's numbers |

## 5. Where the check's value comes from

The cross mean is built at `R/03_metrics.R:243`:

```r
parent_mean <- 0.5 * (mean_source$value[p1] + mean_source$value[p2])
```

`mean_source$value` is a named vector over ids from `ng_choose_mean_source()`
(`R/02_effects.R:195`), which picks **GEBV** when the trait's effects are calibrated and clear
`min_reliability`, else falls back BLUP -> BLUE -> adjusted phenotype. The decision is made
**per trait**, so one run may be GEBV for yield and BLUE for protein.

The check takes its value from the same vector:

```r
check_value <- mean_source$value[check_id]
```

Same source, same trait, same scale as the plotted axis. The user never chooses the basis; the
run does, and the check inherits it.

**Not-evaluable rule.** If the run's source for a trait is BLUE and the check has no BLUE
record, that trait's check reports `NA` with a diagnostic. It must never silently substitute
the check's GEBV against BLUE-scale cross means. A missing check value is reported, never
approximated.

## 6. Architecture

```
parents ──> QC ──> marker set fixed ──> effects_list[[t]] ──> mean_source (per trait)
                          │                     │                      │
                          │                     ▼                      ▼
checks ───────> align to that marker set ──> predict ──────────> check_value[t]
                (never enters QC, LD, GRM, pairs)                      │
                                                                       ▼
                                            reference line + columns + P(beat check)
```

Prediction is **row-wise**: a line's GEBV depends only on its own dosages and the shared
effects, which is why a check can be scored on the parents' model without perturbing it. QC,
LD, and GRM are **set-wise** -- they depend on who else is in the matrix -- which is why checks
must be kept out of them. The design splits exactly along that line.

### Marker alignment (the silent-failure risk)

`ng_predict_gebv()` is a dot product of dosages with marker effects, so the check matrix must
match the parents' on **marker set, column order, and allele coding**. A flipped reference
allele at a subset of markers yields a check GEBV that is wrong but entirely plausible: a
reference line in the right ballpark, in the wrong place, with no warning.

Requirements:

- Hard-fail on marker-set mismatch. Report the overlap count; never quietly intersect.
- Subset the check matrix to the **final** marker set used for effect estimation, in that
  column order.
- Same dosage coding and same reference allele as the parent matrix, validated explicitly.

## 7. P(beat check)

`ng_p_superior_progeny(mu, sigma, tau, k_progeny)` (`R/30_posterior_prediction.R:267`) computes
P(at least one of k progeny >= tau) in closed form. The check enters as **`tau`, a scalar
constant** -- `mu` and `sigma` are the cross's. The check contributes no variance, as required
by section 2.

For `decrease` traits, "beats the check" means below it. The mirrored probability is the same
function with both arguments negated:

```r
ng_p_superior_progeny(-mu, sigma, -tau, k)   # == P(at least one of k progeny <= tau)
```

Proof: the function returns `1 - Phi((tau' - mu')/sigma)^k`. With `mu' = -mu`, `tau' = -tau`
this is `1 - Phi((mu - tau)/sigma)^k = 1 - [1 - Phi((tau - mu)/sigma)]^k`, which is
P(at least one <= tau). No new math.

**Multi-trait.** `ng_p_superior_progeny_multitrait()` (`R/33:25`) already takes both
`tau_lower` and `tau_upper`. `increase` traits fill `tau_lower`, `decrease` traits fill
`tau_upper`, giving **P(a progeny beats every check at once)** natively.

## 8. Plots

A check has a value on the mean axis and **no coordinate on any other axis**. It is therefore a
**reference line**, never a marker -- giving it a diversity position would invent a value it
does not have.

**Its orientation follows the mean-bearing axis, and is not fixed:**

| Plot shape | Mean axis | Check is drawn as |
|---|---|---|
| score / mean on **y**, diversity on **x** | y | a **horizontal** line spanning the x-range |
| diversity on **y**, mean GEBV / phenotype on **x** | x | a **vertical** line spanning the y-range |
| neither axis carries a mean (rank vs rank, kinship vs kinship) | none | **no line** |

The rule: the check line is **perpendicular to the mean-bearing axis**. A plotting helper must
therefore be told which axis that is rather than assuming `h =` -- hard-coding a horizontal line
silently produces a wrong-orientation reference on any plot that puts the mean on x.

**Not every plot is an appropriate host.** The line belongs only where a reader would naturally
ask "is this above or below my check", which requires the axis to carry the trait mean or a
monotone transform of it. Elsewhere it is visual noise at best, and an invitation to misread at
worst.

Target plots:

| Plot | x | y | Point is |
|---|---|---|---|
| `ng_plot_priority_score_vs_kinship` (`R/38:22`) | Pair kinship | Multi-trait score | one cross |
| `ngcd_frontier_plotly` (frontend `R/helpers.R:420`) | Group coancestry | Mean gain | one plan |

Rendering:

- One line per active check, oriented per the table above and labelled with the check id.
- Points on the wrong side get de-emphasised styling (the flag).
- **P(beat check) encoded as marker opacity or size**, so a below-the-line cross with a
  superior tail is visibly distinct from a genuinely dead one. Without this, every point under
  the line looks equally hopeless -- which is the failure mode a mean-only reference line
  introduces.

### What is on the y-axis

The check's value is per-trait; the scatter's y is an index. Validity depends on the index:

- **Single-trait run** -- there is no shortcut here: the axis is still whatever the run's
  scoring family produced, so a single-trait run follows the same rule as its family above. The
  scatter's y is `multi_trait_score`, which is **transformed**
  by `ng_score_breeder_objective()` even for one trait (`R/19_multi_trait_selection.R:576`:
  `value_z <- (oriented - center) / scale`, with a robust per-trait `center`/`scale` taken
  from the candidate distribution). The check must therefore be mapped through that same
  affine transform before it is drawn -- plotting `check_value` in raw trait units puts the
  line entirely off-axis. The transform is deterministic, so this places the check where it
  genuinely falls; it is not an approximation. When the transform is unavailable, draw no
  line (see the rank-based case).
- **`economic_index` / `desired_gain`** -- these genuinely use `value_z`, so map each trait's
  check through its own `center`/`scale` and combine with the **solved** index coefficients
  (`economic_index_coefficients` / `desired_gain_coefficients`), not the raw input weights --
  those are a different vector in general.
  Applying raw coefficients to raw trait values is wrong for the same reason as above.
- **Rank-based scoring** (`auto`, `threshold`, `weighted` -- everything routed through
  `z %*% weights`, where `z` is `ng_rank_normalize()`): the axis is a rank-normal quantile
  scale, not trait units, and not `value_z`. **Verified by running the scoring path**, not
  inferred from the family label: only `economic_index` and `desired_gain` use `value_z`.

  An earlier draft of this design said no line could be drawn here, on the grounds that "a
  check has no rank because it is not a cross." **That reasoning was too strong.** A check's
  *value* can be positioned within the candidate distribution without the check being a member
  of it -- that is a quantile lookup, and it is precisely the question a breeder asks ("where
  does my check sit among these crosses?"). So place it properly: take the check's fractional
  rank among the candidate means, map it through the same rank-normal transform, and apply the
  same standardisation, without perturbing the candidates.

  Getting this wrong is not a rounding matter. Placing a check with the `value_z` transform on a
  rank-normal axis put a measured example at `-0.327` where its true position was `+0.688` --
  opposite sides of zero, about 20% of the axis span, changing how many crosses appear to clear
  the check.

### Multi-trait panel

One axis cannot carry N check lines honestly. Multi-trait checks get **per-trait small
multiples**: one facet per trait with a check, each with its own line on its own scale, plus
the single `p_beat_all_checks` summary column.

**Implementation note (M1, post-ship):** the rank-based quantile placement above (and the
`economic_index`/`desired_gain` affine mapping) is what was DESIGNED, but what actually shipped
(D4, commit 7d2e432) is more conservative: the main scatter's line is drawn only when
`trait_value_metric == "mean"`, and refuses (no line) under every other metric -- including
`"usefulness"`, the package default -- rather than attempting the quantile/affine placement this
section describes. The measured ~22% misplacement rate the D4 fix cites is for exactly that
attempted placement, which is why it was pulled rather than shipped as designed. Consequently the
per-trait panel above is not merely a "multi-trait" convenience: since its y-axis is always
`<trait>_mean` (never the ranking metric), it is commensurable with the check under every
`trait_value_metric`, and is drawn for ANY number of active checks (including exactly one) -- it
is the fallback that makes a check visible at all on a default (non-"mean"-metric) run, not an
enhancement reserved for multi-trait runs.

## 9. Excel output

**New `Checks` sheet** -- the reference itself, one row per trait:

| trait | check_id | direction | value | source | n_crosses_on_wrong_side |
|---|---|---|---|---|---|
| yield | CHK_A | want above | 3.81 | GEBV | 14 / 96 |
| protein | CHK_A | want above | 12.40 | BLUE | 51 / 96 |
| days_to_flower | CHK_B | want below | 71.2 | GEBV | 8 / 96 |

**Columns on the existing cross sheets** (`Selected_All`, `Candidate_Crosses`, tier sheets --
built by `ng_cpw_make_selected`, `R/37_cross_priority_workbook.R:346`), per trait with a check:

- `<trait>_check_value`
- `<trait>_vs_check` -- signed margin, **direction-aware so positive always means better**
- `<trait>_check_ok`
- `<trait>_p_beat_check`

Plus one summary pair: `checks_all_ok`, `p_beat_all_checks`.

The direction-aware margin matters: for a `decrease` trait the raw difference flips sign, and a
breeder scanning a column should not have to remember which traits invert.

## 10. Backend changes

| File | Change |
|---|---|
| `R/44_trait_checks.R` | Rewrite. Keep `ng_trait_check_spec()` direction resolution. Replace `ng_apply_trait_checks()` with `ng_attach_check_reference()`: attaches columns, never subsets rows. Delete the `exclude` branch. |
| `R/02_effects.R` | `ng_choose_mean_source()` gains an optional `extra_ids` / `extra_geno` so check ids land in `$value` without entering `ids` used for pairs. |
| `R/03_metrics.R` | Carry `check_value` per trait alongside `cross_mean_blend`; no change to `parent_mean`. |
| `R/39_cross_prediction_runner.R` | New formals `check_geno`, `check_pheno`. Delete `check_basis`, `exclude_threshold_violators`, and the parent-only validation at :1215. Marker-alignment gate. `trait_check_reference` block in the JSON envelope. |
| `R/30` / `R/33` | No change -- called as-is with `tau` = check value. |
| `R/37_cross_priority_workbook.R` | `Checks` sheet + per-trait columns in `ng_cpw_make_selected()` / `ng_cpw_candidate_table()`. |
| `R/38_cross_priority_plots.R` | Reference line + flag styling + opacity-by-P(beat check) in `ng_plot_priority_score_vs_kinship()`; per-trait facet plot. |
| `NAMESPACE`, `man/` | Export `ng_attach_check_reference`; retire the veto exports. |

## 11. Frontend changes (NextGenCrossDesign 0.26.0 -> 0.27.0)

| File | Change |
|---|---|
| `R/app.R` | Check-file upload card in the guided per-file import flow (its own preview + column mapping + alignment badge, matching the existing per-file pattern). Replace the check picker (`trait_check_pickers`, :1109): choices come from the **check file**, not the genotype table. Drop the `check_basis` control and the exclude toggle. |
| `R/helpers.R` | `ngcd_build_trait_checks()` (:562) drops `bases`; assembles trait + check id + direction override only. |
| `R/ui_charts.R` | Reference-line layer on the results scatter; opacity by `p_beat_check`; per-trait facet chart. |
| `R/report.R` | Reference line in the static/PDF report path. |
| `R/diagnostics.R` | Replace `ngcd_diag_trait_check` veto note with a reference note: n on the wrong side per trait, and any not-evaluable checks. |
| `inst/app/tools/run_cross_prediction_json.R` | Pass the check file through; keep the `as_rows_df` coercion for `trait_checks` (:548). |

**Input-ID conservation guard** applies. The frontend's static input-ID set is baselined in
`.superpowers/sdd/baseline-input-ids.txt`; this change deliberately adds ids (check upload,
check pickers) and removes others (`check_basis`, the exclude toggle), so the baseline must be
re-extracted and the diff reviewed deliberately -- not bypassed. Config keys must also
invalidate a stage, or `tests/testthat/test-pipeline-state.R:262` reports them as orphaned.

## 12. Validation

**The invariant that proves "reference only":**

> A run **with** checks and the same run **without** checks must produce identical cross
> predictions -- same means, same variances, same ranking, same plan, same allocation. The only
> difference is added reference columns.

This is the primary regression test. Additional tests:

1. Check not in the parent set is accepted (the case the old veto rejected).
2. A check id that collides with a parent id is rejected with a clear error.
3. Marker-set mismatch between check and parent matrices hard-fails.
4. Direction: `increase` flags below-check crosses; `decrease` flags above-check.
5. `<trait>_vs_check` is positive-is-better for both directions.
6. Not-evaluable: run source is BLUE, check has no BLUE record -> `NA` + diagnostic, never a
   GEBV substituted against a BLUE axis.
7. `p_beat_check` for a `decrease` trait equals the hand-computed
   `1 - Phi((mu - tau)/sigma)^k`.
8. **Superseded by the shipped D4 behaviour (see the implementation note in section 8):** the
   main scatter's check line is drawn only when `trait_value_metric == "mean"`; every other
   metric (rank-based or not -- including `"usefulness"`, the default) refuses (`NA`, no line),
   not merely when the quantile-placement inputs happen to be missing. The per-trait small-
   multiples panel is the fallback for every other metric, at any number of active checks
   (including one), since its axis is always the trait mean.
9. Multi-trait: `p_beat_all_checks` matches `ng_p_superior_progeny_multitrait()` with
   `tau_lower`/`tau_upper` filled by direction.
10. Excel `Checks` sheet present and populated; per-trait columns on every cross sheet.
11. Zero checks supplied -> output byte-identical to today (no empty sheet, no empty columns).

## 13. Out of scope

- Multiple checks per trait (v1 is one).
- Raw-number thresholds as an alternative to a check genotype.
- Checks influencing selection, ranking, allocation, or the objective in any way.
- Restoring the veto behind a flag.
- A variance for the check (would require selfing).

## 14. Rollout

1. Backend on a working branch; R CMD check + full test suite.
2. Backend 0.23.0, NEWS entry documenting the breaking removal of
   `exclude_threshold_violators` and `check_basis`.
3. Frontend 0.27.0 against backend 0.23.0; bump the backend floor.
4. Mirror PR / CI / merge in both repos.
5. Rebuild the workbench image once both are merged.
