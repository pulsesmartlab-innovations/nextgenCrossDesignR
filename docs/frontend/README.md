# Frontend integration kit — start here

This folder is the **framework-agnostic contract** for building a decision-workbench frontend on the
`nextgenCrossDesign` backend. The backend ships **quantitative-genetics logic + versioned JSON/XLSX
only**; the UI is built separately in whatever framework you like (React/Next/Vue/Shiny). You never
reimplement the science in the frontend — you **spawn R and read JSON**.

## What's here

| File | Purpose |
|---|---|
| **`README.md`** (this) | The step-by-step build order and the whole contract at a glance |
| `contracts/config_schema.json` | Machine-readable **input** contract — every run parameter, grouped, with type/default/allowed/gated-by. Drive your config form from this. |
| `contracts/example_config.json` | A runnable example config (points at `contracts/data/`) |
| `contracts/example_result.json` | The **output** contract — a real `ng_run_result.v1` payload to code your result views against |
| `POINT_AND_CLICK_FRONTEND_DEVELOPER_GUIDE.md` | Screen-by-screen build guide (14 screens), runner template, result contract, example-coverage acceptance tests |
| `MATING_STRATEGY_INTEGRATION.md` | UI-control → backend-function/parameter/JSON-field mapping; the versioned decision-export schemas |
| `FRONTEND_DEVELOPER_TEMPLATE.md` | Product position, information architecture (module → backend function), state model, acceptance checklist |

Everything is **schema-versioned** — branch on the `schema` / `schema_version` string, never assume a
single global version.

---

## The 7-step build order

### 0. Load + feature-detect
Use the **installed** package (`library(nextgenCrossDesign)`) in production. Then fetch the capability
manifest once and gate your UI on it:

```bash
NG_BACKEND_CAPABILITIES_OUT=results/backend_capabilities.json \
  Rscript tools/export_backend_capabilities_json.R
```
→ `backend_capabilities.json` (`schema_version: ng_backend_capabilities.v1`): `data_qc`,
`breeding_systems`, `method_families` (24 rows, incl. `training_set_augmentation`), `workflows`,
`external_integrations`, `navigation`. Render menus/capabilities from this instead of hard-coding.

**Optional-dependency gating** (feature-detect; don't assume): `lpSolve` → MIP optimizers (else
greedy/repair), `openxlsx` → xlsx workbook, `ggplot2` → figures, `PopVar`/`SimpleMating`/AlphaMate →
*exact* external baselines (native "style" proxies always available).

### 1. Data readiness (before any run)
```bash
NG_PREFLIGHT_GENO=geno.csv NG_PREFLIGHT_PHENOTYPE=pheno.csv NG_PREFLIGHT_MARKER_MAP=map.csv \
  NG_PREFLIGHT_TRAIT_SPEC=direction.csv NG_PREFLIGHT_OUT=results/data_preflight.json \
  Rscript tools/run_data_preflight.R
```
→ `data_preflight.json` (`ng_data_preflight.v1`): `status` (pass/warning/blocker), `counts`, `tables`,
`issues` (id/severity/table/field/message), `putative_duplicates`. Show a red/yellow/green gate with
per-issue drill-down; a `blocker` must block the run.

### 2. Build the config form
Generate the form from **`contracts/config_schema.json`** — groups → parameters, each with `type`,
`default`, `allowed` (enum choices), and `gated_by`. Do **not** offer `var_simple` as a *merit* metric
(it is a diversity proxy). The parents are the main genotype/phenotype tables; the optional
**training-set** group adds effect-only individuals that are never crossed.

### 3. Run the analysis (headless)
The single entry point — pass a JSON config, get a JSON result:

```bash
NG_RUN_CONFIG=my_config.json NG_RUN_RESULT_OUT=results/run_result.json \
  Rscript tools/run_cross_prediction_json.R
```
The config keys mirror `ng_run_cross_prediction()`; **unknown keys are a hard error** (no silent
drop). This is what your worker spawns. (Advanced object parameters — `committed_crosses`,
`cross_cost`, `parent_group`, `group_permission`, `group_quota`, `trait_weights` — go inline as nested
JSON.)

### 4. Render the result contract
`run_result.json` (`ng_run_result.v1`) — code your views against `contracts/example_result.json`:
- `status` — **`"ok"` or `"error"`. Always branch on this first.** On `"error"` the run was
  **blocked** (e.g. a DH/inbred parent carrying heterozygosity, or a bad input); `error_message`
  holds the human-readable reason to show the user, and the result carries no plan. The runner still
  writes a valid `result.json` on a blocker (exit code is non-zero) — **read the file even on failure**
  so you can display `error_message` instead of a raw stderr dump.
- `warnings` — an array of advisory message strings the run emitted (e.g. the residual-heterozygous-RIL
  variance advisory, the `assume_inbred` deprecation). Present on a successful run; surface each as a
  non-blocking banner/toast so the breeder sees caveats (e.g. "variance biased low — supply phased
  haplotypes"). Empty array when there are none.
- `selected_crosses` — the crossing plan: `parent1, parent2, multi_trait_score, pair_kinship,
  priority_tier, priority_index` (+ per-trait `<t>_value/_mean/_pmv/_vpm/_var_complex`). The breeder deliverable.
- `candidate_crosses` — all scored pairs (for the priority-vs-kinship scatter).
- `plan_summary` — KPI cards: `mean_gain, group_coancestry, unique_parents, max_parent_use,
  mean_progeny_inbreeding, total_cost/budget/over_budget`, and constraint status
  (`constraints_reduced_plan, mip_fallback, target_coancestry_status`, ...).
- `input_match_audit` — `matched_parent_count` (crossed) vs `training_only_count`/`effect_training_n`
  (effects-only). Display "trained on N, crossing K parents".
- `effect_summary` — per trait: `marker_effect_reliability, marker_effect_training_n, ridge_lambda`.
- `qc` — status + issues (mirrors step 1). `objective` — resolved multi-trait method + reason.
- `priority_risk_diagnostics` + the per-cross portfolio columns — see
  `2026-08-15-multitrait-portfolio-handoff.md` and the drop-in `ngcd_portfolio_multitrait.R`.
  Each cross carries `cross_level, cross_upside, cross_confidence, risk_bin, confidence_method,
  portfolio_profile, portfolio_basis` (+ `risk_driver_trait/_share` on multi-trait runs).
  **You must badge `portfolio_basis == "linearized_rank_index"`** — for the rank-based index
  methods (`auto`/`weighted`/`threshold`, and `auto` promotes to `weighted` whenever trait
  weights are present) the axes are indicative, not a decomposition of `multi_trait_score`.
  `cross_confidence`/`risk_bin` are **within-run only** (min–max and tertiles of the full
  post-filter candidate pool, copied unchanged to selected rows). They are not an absolute
  measure of how well a run is estimated; tied uncertainty values can make bins unequal.

  Note on the shipped example: it is a deliberately small 10-parent / 12-marker dataset, so its
  `disease` trait has almost no GEBV spread. That makes it a **useful** reference — it populates
  the warning fields most runs leave null (`pev_concentration_note`,
  `risk_disproportionate_traits`), so you can code those views against real shapes. It also makes
  `index_weights` look alarming (`disease: -3540.7` vs `yield: 0.138`). That is correct, not a
  bug: weights are raw-unit (`coef * sign / scale`), so the 25667x weight ratio is just the
  inverse of the 25670x scale ratio. **Never render weight magnitude as importance** — use the
  `index_traits` shares for that.

### 5. Decision / strategy exports (interactive controls)
See `MATING_STRATEGY_INTEGRATION.md` for the control→field mapping. Key exporters:
- **Gain–diversity frontier** → `ng_write_frontier_json` (`ng_mating_frontier.v1`): the frontier chart +
  `high_gain`/`balanced`/`diversity` preset operating points.
- **Progeny-inbreeding histogram** → `ng_write_progeny_inbreeding_histogram_json`
  (`ng_progeny_inbreeding_histogram.v1`).
- **Compare runs** → `ng_write_run_comparison_json` (`ng_run_comparison.v1`): side-by-side plan metrics.
- **Breeder workbook** → `ng_write_cross_priority_workbook` (.xlsx) / `ng_cross_priority_workbook_tables`
  (per-sheet JSON: Dashboard, tier sheets, Parent_Use_QC, Duplicate_QC, editable breeder columns).

### 6. Benchmark evidence + claims discipline
- **Head-to-head dashboard** → `ng_write_head_to_head_dashboard_json` (`ng_head_to_head_dashboard.v1`):
  `recommended_method`, `method_summary`, per-metric directions.
- **Claims discipline (mandatory):** every benchmark row carries `implementation`,
  `exact_external_status`, `fallback_reason` — badge each result as **exact-external vs style-proxy**
  and never conflate them. Copy the claim boundaries from **`VALIDATED_STATE.md`** (Allowed vs
  Disallowed): the validated path is diploid DH/RIL; multi-trait superiority is not yet claimable at
  production scale; polyploid is experimental/guarded. Encode as an `EvidenceStatus` badge.

---

## Build on it incrementally

- **MVP:** step 0 (registry) → step 2 (config form) → step 3 (headless run) → step 4 render
  `selected_crosses` + `plan_summary` + the workbook download. That is a usable decision tool.
- **Then add:** step 1 QC gate, the priority-vs-kinship scatter (`candidate_crosses`), the strategy
  dial + frontier (step 5), run comparison, and the benchmark/evidence view (step 6).
- **Advanced:** training-set augmentation (step 2 group), posterior prediction, and the exact
  within-family cross-trait covariance for the multi-trait threshold probability
  (`ng_cross_trait_within_family_cov` → `ng_add_p_superior_progeny_multitrait(cross_trait_cov=)`).

## Do not reimplement in the frontend
Marker-effect estimation, recombination variance, PMV/usefulness, OCS/MIP allocation, kinship, priority
ranking, QC/duplicate detection, multi-trait objectives. All of it lives in the backend behind the
contracts above — the UI orchestrates and renders.
