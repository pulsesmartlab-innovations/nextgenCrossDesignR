# Polyploid staged pipeline — design (Phase 1: autotetraploid)

**Date:** 2026-08-08
**Backend:** `nextgenCrossDesign` · **Frontend:** `nextgenCrossWorkbench` (0.21.0 stepped Run tab already shipped)

## Goal

Let the **autotetraploid** cross-design run as the same stepped, compute-once pipeline the diploid
workflow now uses — **QC → Fit effects & score → Allocate & rank** (single-trait, no index step) —
each step run on its own, cached so a later step never re-runs an earlier one, and each step showing
its own figure in the frontend. Disomic-subgenome is Phase 2 (separate spec).

## Why this is safe

The diploid staged runner is **validated and shipped**; it must stay byte-identical. The polyploid
paths are *standalone single-shot functions* (`ng_polyploid_design_crosses`, `R/46`) with no `ctx`/
`run_dir`/stage decomposition — there is no shared code to perturb. We add a **parallel** poly stage
pipeline that reuses the diploid runner's generic persistence helpers (`R/45_staged_runner.R`:
`ng_stage_save`/`ng_stage_load_ctx`/`ng_stage_write_json`/`ng_stage__update_manifest`) but its OWN
stage functions. The diploid `ng_cp_pipeline` is not touched.

## Non-negotiable invariants

1. **Byte-identical equivalence:** driving the autotetraploid design through the staged pipeline
   (qc→predict→allocate→rank) produces a result identical to the current one-shot
   `ng_polyploid_design_crosses(...)`. A golden test asserts this.
2. **Compute-once:** each poly stage runs exactly once; its accumulated `ctx` is snapshotted to
   `run_dir/artifacts/<stage>.rds` and reloaded by the next stage (same mechanism as diploid).
3. **Diploid path untouched:** `ng_cp_pipeline`, `ng_cp__build_ctx`, `ng_run_stage` diploid behavior
   unchanged; a diploid golden test still passes.
4. **Attributes & masks persist across stage boundaries** (the map's flagged hazard): the QC
   marker-keep mask and the `parent_kinship`/`parent_gebv`/`ploidy` attributes must survive the RDS
   round-trip — they are attributes, not columns, so the poly `ctx` stores them explicitly.

## Backend design

### New file `R/50_polyploid_staged_runner.R`

- **`ng_poly_cp_stage_order() -> c("qc","predict","allocate","rank")`** (no `index` — single-trait).
- **`ng_poly_cp__build_ctx(config)`** — normalizes the poly config surface (dosage/effects/phenotype,
  `ploidy`, `poly_dominance`, `poly_grm_method`, `double_reduction`, `selection_prop`, `gain`,
  allocation controls). Mirrors `ng_cp__build_ctx` but for the poly argument set.
- **Stage functions** (each takes `ctx`, returns the augmented `ctx`):
  - `ng_poly_cp__stage_qc(ctx)` — `ng_polyploid_qc(dosage, ploidy, ...)`; stores `ctx$geno` (cleaned),
    `ctx$marker_keep` (mask), `ctx$qc` (status/issues/marker_report). On a blocker sets
    `ctx$qc$status="blocker"` and returns early (no throw), exactly like the diploid qc stage.
  - `ng_poly_cp__stage_predict(ctx)` — realigns supplied `effects` to `ctx$marker_keep` (or fits:
    `ng_fit_ridge_effects` / `ng_polyploid_fit_effects` when `poly_dominance`); scores crosses
    (`ng_polyploid_score_crosses` or `..._dominance`); stores `ctx$scored` (+ `parent_kinship`,
    `parent_gebv`, `ploidy` as explicit `ctx` fields, NOT relying on the frame's attributes surviving).
  - `ng_poly_cp__stage_allocate(ctx)` — `ng_optimize_mating_plan(ctx$scored, n_crosses, gain_col,
    parent_kinship = ctx$parent_kinship, ...)`; stores `ctx$plan`.
  - `ng_poly_cp__stage_rank(ctx)` — assemble the final selected-cross table (priority annotation if
    cheap; otherwise pass through the plan). On `rank`, attach `attr(out,"result") <-
    ng_poly_cp__assemble_result(ctx)` — the SAME result contract the one-shot poly design returns.
- **`ng_poly_cp_pipeline <- list(qc=, predict=, allocate=, rank=)`** dispatch table.
- **`ng_poly_run_stage(stage, run_dir, config)`** — the poly analog of `ng_run_stage`: qc builds ctx
  from config; later stages `ng_stage_load_ctx(run_dir, <upstream>)`, run
  `ng_poly_cp_pipeline[[stage]](ctx)`, then `ng_stage_save`/`ng_stage_write_json`/manifest (reuse R/45
  helpers). QC-blocker gate identical to diploid.
- **`ng_poly_cp__stage_json(stage, ctx)`** — per-stage JSON summary shapes for the frontend:
  qc→`{status, issues, marker_report, putative_duplicates}`; predict→`{effect_summary?,
  n_candidates, poly_metric}`; allocate→`{plan_summary:{n_crosses,mean_gain,...}, parent_use}`.

### Dispatch: the FRONTEND-bundled runner `inst/app/tools/run_cross_prediction_json.R`

**Correction after investigation:** the runtime JSON runner is shipped by the *frontend*
(`ngcd_res("tools", ...)` → `inst/app/tools/run_cross_prediction_json.R`); the backend `tools/` copy is
a stale dev artifact. The poly orchestration lives in that frontend runner: `run_polyploid_design(raw,
result_path)` marshals the config (dosage/ploidy/phenotype/qc/dominance/gain/selection_prop/
double_reduction/grm_method/seed + allocator `...`) and calls the backend one-shot
`ng_polyploid_design_crosses`.

So the dispatch is added in the FRONTEND runner: refactor `run_polyploid_design`'s config→`design_args`
marshaling into a shared helper `poly_design_args(raw)`; add a staged branch that, for `workflow="stage"`
+ poly, calls `nextgenCrossDesign::ng_poly_run_stage(stage, dirname(result_path), poly_design_args(raw))`
and writes the per-stage JSON. The one-shot path keeps calling `ng_polyploid_design_crosses`. This makes
the staged and one-shot poly paths share the exact same argument marshaling (equivalence by
construction at the frontend boundary, on top of the backend golden equivalence test).

## Frontend design (reuse 0.21.0 engine)

- **Route autotetraploid through staged mode.** `run_mode()` returns `"staged"` for autotetraploid
  (not `"single"`); the adaptive cards already render **3 steps** for single-trait (index hidden), so
  poly shows QC → Fit & score → Allocate. Subgenome stays `"single"` until Phase 2.
- **Un-skip poly in the staleness observer** (`app.R:1321`) for autotetraploid only; keep subgenome
  skipped. The poly stage key-subsets (`ngcd_stage_key_patterns`) gain poly keys so a poly config
  change marks the right stage stale (compute-once for poly).
- **Poly-aware stage figures/summaries:** `ngcd_stage_summary()` + the `fig_*` renderers handle the
  poly stage-JSON shapes (poly QC dupes, poly candidate count, poly frontier/parent-use). Where shapes
  match diploid, reuse; add poly branches where they differ.
- **`ngcd_run_stage()`** already sets `workflow="stage"`; it just needs to pass the poly config keys
  (it already forwards `build_params()`, which includes `ploidy`/`poly_*`).

## Testing

- **Backend golden equivalence:** `tests/polyploid_staged_equivalence.R` — run an autotetraploid design
  both one-shot (`ng_polyploid_design_crosses`) and staged (`ng_poly_run_stage` qc→predict→allocate→
  rank over a temp run_dir); assert the selected-cross tables + key metrics are identical.
- **Compute-once:** assert running `allocate` reloads the cached `predict` ctx (does not refit effects)
  — e.g. stamp a sentinel into the predict artifact and confirm allocate sees it without recompute.
- **QC-blocker gate:** a contaminated dosage (out-of-range) blocks at qc; predict/allocate abort.
- **Attribute/mask persistence:** after a full staged run, `parent_kinship` and the marker mask are
  present and correctly shaped in the reloaded ctx.
- **Diploid regression:** the existing diploid staged/equivalence tests still pass (path untouched).
- **Frontend:** the stepped cards render for autotetraploid (3 cards); compute-once testServer for a
  poly allocation-only change; `R CMD check` green both repos.

## Out of scope (Phase 1)

Disomic-subgenome staging (Phase 2 — the list-of-matrices + parent-order/attribute hazards); a
polyploid multi-trait index step (poly scoring is single-trait); the AlphaSimR simulation poly path
(not JSON/file-drivable). No change to diploid staging or the poly *science* (same QC/effects/scoring/
allocation functions, just orchestrated in stages).
