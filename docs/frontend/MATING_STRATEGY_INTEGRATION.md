# Mating-Strategy Integration Spec (for the frontend project)

This document is for the developer building the decision-workbench UI. It describes new
backend capabilities — gain–diversity balancing, progeny-inbreeding management, breeder
mating constraints, marker steering / lethal-allele guarding, and cost/logistics — and
maps each UI control to the exact backend function, parameter, and JSON field that feeds
it. **The backend ships logic + JSON only; no UI code lives here.**

These capabilities are our own implementation of ideas common to tactical mate-selection
tools (cf. Kinghorn 2011, Genet. Sel. Evol. 43:4). Use our vocabulary in the UI:
"diversity emphasis" (not "degrees of emphasis"), "mating-group permissions" (not
"GroupFix"), "min-use-if-used" (not "MinUse").

All new parameters default to off, so existing runs are unchanged.

## Where to call these

Every capability below is reachable from **both** entry points, not just the low-level
allocator:
- `ng_run_cross_prediction(...)` — the user-friendly runner the UI drives. Accepts
  `strategy` / `diversity_emphasis`, `lambda_progeny_inbreeding`, `min_crosses_per_parent`,
  `committed_crosses`, `parent_group` / `group_permission` / `group_quota`, `cost_col` /
  `budget` / `lambda_cost`, `logistic_col` / `lambda_logistic`, `lethal_spec`, and
  `optimizer` (now including `"evolution"`, the native memetic-GA allocator, plus
  `evol_solutions/iterations/stop/seed`).
- `ng_design_crosses(...)` — the scripting pipeline. Same optimizer options via `...`,
  plus `marker_target_spec` / `lethal_spec` / `lambda_marker` for marker steering.

Note on marker **steering** (driving an allele frequency): it is applied at the scoring
layer via `ng_design_crosses` / `ng_apply_marker_management`. In the multi-trait runner,
marker **guarding** (`lethal_spec`) is wired in; steer a favourable major gene by adding it
as a trait in the index rather than blending a single marker into the economic index.

---

## 1. Gain–diversity balance (the strategy dial)

The single most important control. The breeder picks a strategy or a numeric emphasis and
the backend returns the mating plan at that point on the gain-vs-diversity frontier.

- Backend: `ng_optimize_mating_plan(scores, n_crosses, parent_K, strategy=, diversity_emphasis=)`
  - `strategy`: `"high_gain"` | `"balanced"` | `"diversity"` (presets at emphasis 15 / 45 / 75).
  - `diversity_emphasis`: numeric **0 = maximum genetic gain … 100 = maximum diversity**.
  - Returned `attr(plan, "summary")` carries `diversity_emphasis`, `achieved_emphasis`,
    `emphasis_gap`, `selected_lambda_group`.
  - **Alternative (constrained OCS):** `target_coancestry=` picks the SAME frontier by an
    absolute cap instead of a relative dial — "maximize gain at/under this group coancestry"
    (the classic Meuwissen rate-of-inbreeding constraint; `group_coancestry` is ~2× the
    coancestry/ΔF coefficient). Summary adds `target_coancestry`, `achieved_coancestry`,
    `target_coancestry_status` (`met`/`met_slack`/`infeasible`). It takes precedence over the
    dial. Expose EITHER the slider OR a ΔF-target input, not both at once.
- Frontier for the slider/scatter: `ng_write_frontier_json(scores, n_crosses, parent_K, output_path)`
  → schema `ng_mating_frontier.v1`:
  - `points[]`: `{ lambda_group, diversity_emphasis, mean_gain, group_coancestry,
    mean_progeny_inbreeding, unique_parents, max_parent_use }` — plot `mean_gain` (y) vs
    `group_coancestry` (x); the slider position maps to `diversity_emphasis`.
  - `strategies.{high_gain,balanced,diversity}`: the resolved operating point for each preset
    (drop the three marker dots on the frontier).

**UI:** a High Gain ↔ Balanced ↔ Diversity slider bound to `diversity_emphasis`; the frontier
scatter with the current point highlighted and the three preset dots.

## 2. Progeny-inbreeding management + histogram

- Per-cross metric: `ng_score_crosses(...)` now returns `expected_progeny_inbreeding`
  (parental coancestry; 0..~1 relationship scale).
- Objective emphasis: `ng_optimize_mating_plan(..., lambda_progeny_inbreeding=)`. Summary
  exposes `mean_progeny_inbreeding` / `max_progeny_inbreeding`.
  - **Relatedness axes (avoid double-counting):** `lambda_progeny_inbreeding` and the
    pre-existing `lambda_mating` act on the **same** parent-pair-relatedness axis and their
    penalties **add** — expose only ONE "avoid inbred/related matings" control in the UI
    (default to `lambda_progeny_inbreeding`, interpretable coancestry units). `lambda_group`
    is a **separate** population-level diversity control (bind it to the strategy/diversity
    slider). If both `lambda_mating` and `lambda_progeny_inbreeding` are set, the plan summary
    flags `relatedness_penalty_overlap = TRUE`.
- Histogram: `ng_write_progeny_inbreeding_histogram_json(scores_or_plan, output_path, breaks=)`
  → schema `ng_progeny_inbreeding_histogram.v1`: `{ n, mean_progeny_inbreeding,
  max_progeny_inbreeding, bins:[{bin_lower,bin_upper,bin_mid,count,proportion}] }`.

**UI:** a progeny-inbreeding histogram for the current plan; optionally a
"progeny-inbreeding emphasis" control bound to `lambda_progeny_inbreeding`.

## 3. Breeder mating constraints

All passed to `ng_optimize_mating_plan(...)`; all repair-stage (fast) and honest — if a
constraint makes a full-size plan impossible the backend returns a shorter plan and sets
`summary$constraints_reduced_plan = TRUE` plus an R warning (surface it).

| UI control | Parameter | Notes |
|---|---|---|
| Max matings per parent | `max_crosses_per_parent` | existing |
| Min-use-if-used | `min_crosses_per_parent` | a used parent gets ≥ N matings else 0; summary `min_use_if_used` |
| Committed / fixed matings | `committed_crosses` (data.frame `parent1,parent2`) | locked in; summary `n_committed` |
| Mating-group labels | `parent_group` (named vector parent→group) | e.g. heterotic pools, M/F |
| Group permission matrix | `group_permission` (groups × groups logical) | illegal group pairs filtered out |
| Group quotas | `group_quota` (named list, key `ng_group_pair_key(gA,gB)` → cap) | cap crosses per group pair |

**UI:** a committed-matings table editor; a group-assignment + permission-matrix editor; a
per-group-pair quota input. Always show `constraints_reduced_plan` when true.

## 4. Marker steering & lethal-allele guarding

Scoring-layer (needs the genotype matrix). Helper: `ng_apply_marker_management(scores, geno,
marker_target_spec=, lethal_spec=, ploidy=, lambda_marker=, drop_lethal_carrier_crosses=)`
returns the augmented `scores` (feed to the optimizer on `marker_adjusted_gain`).

- Marker steering: `ng_marker_target_spec(marker, direction="increase"|"decrease", target_freq=, weight=)`
  → per-cross `marker_target_score` and `marker_freq_<m>` (expected progeny ALT frequency).
- Lethal guarding: `ng_lethal_recessive_spec(marker, risk_allele="alt"|"ref")` →
  `ng_lethal_recessive_cross_risk` flags `lethal_carrier_cross` / `lethal_risk_loci`;
  carrier × carrier crosses are dropped when `drop_lethal_carrier_crosses=TRUE`.

**UI:** a marker-target list (marker, direction/target, weight); a lethal-recessive list
(marker, risk allele); a toggle to exclude carrier × carrier matings; show per-locus
expected progeny frequencies.

## 5. Cost / logistics

`ng_optimize_mating_plan(..., cost_col=, budget=, lambda_cost=, logistic_col=, lambda_logistic=)`.
- Hard `budget` cap on `sum(cost_col)`; soft `lambda_cost`; soft `lambda_logistic` on a
  per-cross distance/difficulty column. Summary: `total_cost`, `budget`, `over_budget`.

**UI:** a budget input + cost/logistic emphasis controls; show `total_cost` vs `budget`.

## 6. Compare runs

`ng_write_run_comparison_json(named_list_of_plans, output_path)` → schema
`ng_run_comparison.v1`: `{ fields[], runs:[{run, mean_gain, group_coancestry,
mean_progeny_inbreeding, diversity_emphasis, total_cost, n_committed,
constraints_reduced_plan, …}] }`.

**UI:** a side-by-side table of saved runs (one column per run).

## 7. Optimizer choice

`ng_run_cross_prediction(..., optimizer=)` picks the mate-allocation engine:
`"auto"` (default) | `"greedy_local"` | `"repair_local"` | `"mip_linear"` | `"mip_contribution"` |
**`"evolution"`** (the native memetic genetic-algorithm allocator; aliases `"ga"` / `"de"` /
`"memetic"`). The evolution optimizer is tuned with `evol_solutions`, `evol_iterations`,
`evol_stop`, `evol_seed`. It optimizes the same OCS objective, warm-started from greedy with
elitism (so it is never worse than greedy). Evidence (optimizer bake-off): evolution achieves
an objective `>=` greedy/repair and competitive with MIP; in the recurrent study it is within
noise of OCS on realized gain.

**UI:** an optimizer dropdown (default "auto"); reveal `evol_*` tuning only when "evolution" is
chosen. Most users should leave it on "auto".

## 8. Trait-value metric (which scoring metric)

`ng_run_cross_prediction(..., trait_value_metric=, uc_variance_source=)`:
`trait_value_metric` ∈ `"var_complex"` (default) | `"uc"` | `"pmv"` | `"vpm"` | `"var_simple"` |
`"mean"`; `uc_variance_source` ∈ `"pmv"` | `"vpm"` | `"var_simple"`.

Evidence-based default and guidance (see `docs/BACKEND_USER_GUIDE.md` "Choosing A Trait-Value
Metric"; config-scoped): default **`var_complex`** (PMV usefulness); `pmv`/`vpm` rank crosses
almost identically. Use **`mean`** for highly polygenic or small/low-h2 training. `var_simple`
is a diversity proxy, NOT a merit metric. For long-term recurrent selection, manage diversity
explicitly (sections 1-2) on top of a good merit metric rather than switching metric.

**UI:** a metric dropdown defaulting to `var_complex`, with a short "recommended" hint
(polygenic/small training → mean). Don't offer `var_simple` as a merit option.

## 9. Polyploid (any-ploidy) mate design

`ng_design_crosses_poly()` is the single entry for polyploid crops: it takes an allele-dosage
matrix (0..ploidy) + effects or a phenotype, runs ploidy-aware QC (`ng_polyploid_qc`), scores with
a correct allele-frequency polyploid GRM (`ng_polyploid_grm`), and runs the same native allocator +
controls (sections 1–5) — nothing ploidy-specific downstream. Registry: `method_families$id`
includes `poly4x_policy` / `poly_subgenome_policy`.

**UI:** a ploidy selector (2/4/6/…) on the data-import screen; surface the QC summary from
`attr(plan, "qc")` (samples × markers, monomorphic dropped, missingness); everything else reuses
the diploid mate-allocation panels.

---

## Complete new-parameter reference (checklist)

Every parameter added in this round, all on `ng_run_cross_prediction()` unless noted. Group
them in an "Advanced / mate-selection" panel; all default to off so the basic UI is unchanged.

| Parameter | Type | Default | Purpose | Summary field(s) |
|---|---|---|---|---|
| `strategy` | enum | NULL | gain-diversity preset (high_gain/balanced/diversity) | `strategy`, `achieved_emphasis` |
| `diversity_emphasis` | 0–100 | NULL | gain↔diversity dial | `diversity_emphasis`, `emphasis_gap` |
| `target_coancestry` | number | NULL | constrained-OCS inbreeding cap | `target_coancestry`, `achieved_coancestry`, `target_coancestry_status` |
| `lambda_progeny_inbreeding` | ≥0 | 0 | progeny-inbreeding penalty | `mean_progeny_inbreeding`, `max_progeny_inbreeding` |
| `min_crosses_per_parent` | int | 0 | min-use-if-used | `min_use_if_used` |
| `committed_crosses` | df(parent1,parent2) | NULL | locked matings | `n_committed` |
| `parent_group` | named vec | NULL | group label per parent | — |
| `group_permission` | matrix | NULL | legal group×group pairs | — |
| `group_quota` | named list | NULL | max crosses per group pair | — |
| `cost_col` / `budget` / `lambda_cost` | col / num / ≥0 | NULL/Inf/0 | cost cap + soft penalty | `total_cost`, `budget`, `over_budget` |
| `logistic_col` / `lambda_logistic` | col / ≥0 | NULL/0 | distance/difficulty penalty | — |
| `lethal_spec` | spec | NULL | drop carrier×carrier matings | — |
| `marker_ploidy` / `drop_lethal_carrier_crosses` | int / bool | 2 / TRUE | lethal-guard options | — |
| `optimizer` (+`evol_*`) | enum | "auto" | allocation engine incl. "evolution" | `mip_fallback` |
| `trait_value_metric` / `uc_variance_source` | enum | var_complex / pmv | scoring metric | — |
| `grm_method` | enum | vanraden | relationship matrix (vanraden / yang) | — |
| `marker_target_spec` / `lambda_marker` | spec / ≥0 | NULL / 0 | steer a marker's allele frequency (now ON the runner) | `marker_target_score` col on `candidate_crosses` |
| `training_genotype[_file]` / `training_phenotype[_file]` (+`_id_col`) | table/path | NULL | extra effect-only individuals that sharpen the ridge fit but are NEVER candidate parents | `training_only_count`, `effect_training_n`, `marker_effect_training_n` |
| `constraints_reduced_plan` | (output) | — | TRUE when a constraint shrank the plan | surface as a warning banner |

`marker_target_spec` / `lambda_marker` are now first-class `ng_run_cross_prediction()` arguments
(also available at the scores level via `ng_design_crosses` / `ng_apply_marker_management`).
`window_cm` and `min_unique_parents` are now honored on every scoring/allocation path.

**Exact within-family cross-trait covariance** (for the multi-trait threshold probability) is a
companion, not a runner argument: `ng_cross_trait_within_family_cov()` computes the recombination-aware
`a_t' R a_s` per cross, fed to `ng_add_p_superior_progeny_multitrait(cross_trait_cov = ...)` in place of
the population-correlation proxy (example `32_exact_cross_trait_covariance.R`).

---

## Discovery

`ng_backend_capability_registry()` / `ng_write_backend_capability_registry_json()` advertise
the available methods/modules so the UI can feature-detect. The export schemas above are
versioned (`*.v1`); check `schema` before rendering.

The **full run-parameter menu** (not just the mate-selection subset above) is machine-readable in
`docs/frontend/contracts/config_schema.json`, and the whole build order is in `docs/frontend/README.md`.
Drive runs through the headless `tools/run_cross_prediction_json.R` (JSON config → `ng_run_result.v1`).
