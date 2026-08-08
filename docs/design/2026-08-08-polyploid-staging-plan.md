# Polyploid staged pipeline — Implementation Plan (Phase 1: autotetraploid)

> REQUIRED SUB-SKILL: superpowers:executing-plans (inline). Steps use `- [ ]`.

**Goal:** Run the autotetraploid design as a stepped, compute-once pipeline (QC → Fit & score →
Allocate & rank), each step run on its own with its own figure — reusing the diploid staging engine's
persistence + the frontend's 0.21.0 stepped cards, without touching the validated diploid path or the
poly science.

**Architecture:** Backend gains a *parallel* poly stage pipeline (`ng_poly_run_stage`) that decomposes
`ng_polyploid_design_crosses` and reuses the generic `R/45` persistence helpers. The frontend-bundled
runner (`inst/app/tools/run_cross_prediction_json.R`) dispatches `workflow="stage"`+poly to it, sharing
the SAME arg-marshaling as the one-shot path. The frontend app routes autotetraploid through the staged
cards.

## Global Constraints

- **Byte-identical:** staged autotetraploid == one-shot `ng_polyploid_design_crosses` (golden test).
- **Compute-once:** each poly stage snapshots its `ctx` to `run_dir/artifacts/<stage>.rds`; a later
  stage reloads it (never recomputes upstream). Reuse `ng_stage_save`/`ng_stage_load_ctx`/
  `ng_stage_write_json` semantics (but a poly JSON writer) / `ng_stage__update_manifest`.
- **Diploid untouched:** no edits to `ng_cp_pipeline`, `ng_run_stage`, `ng_cp__build_ctx`.
- **Attributes/masks persist:** the QC marker-keep mask and `parent_kinship`/`ploidy` live as explicit
  `ctx` fields (attributes do not reliably survive being rebuilt in a later stage).
- **Version coordination:** backend minor bump (ships `ng_poly_run_stage`); frontend requires it.

---

### Task 1 (backend): poly stage pipeline `R/50_polyploid_staged_runner.R`

Decompose `ng_polyploid_design_crosses` (`R/46:73`) into a ctx-threaded pipeline. Reuse R/45 helpers.

**Files:** Create `R/50_polyploid_staged_runner.R`; export the public fns in `NAMESPACE`; Test
`tests/polyploid_staged_equivalence.R`.

**Interfaces (Produces):**
- `ng_poly_cp_stage_order() -> c("qc","predict","allocate","rank")`
- `ng_poly_cp__build_ctx(config) -> ctx` (list with `dosage,n_crosses,ploidy,effects,phenotype,pairs,
  max_crosses_per_parent,ridge_seed,run_qc,qc,dominance,gain,selection_prop,double_reduction,
  grm_method,alloc_dots`)
- `ng_poly_cp_pipeline` (named list qc/predict/allocate/rank of `function(ctx) ctx`)
- `ng_poly_run_stage(stage, run_dir, config)` — poly analog of `ng_run_stage`
- `ng_poly_cp__assemble_result(ctx) -> plan` (the same object the one-shot returns)

- [ ] **Step 1 — failing golden test** `tests/polyploid_staged_equivalence.R`:

```r
# staged autotetraploid == one-shot ng_polyploid_design_crosses, byte-for-byte on the plan.
source(...helper_load...); set.seed(11)
P <- 24; M <- 60; ploidy <- 4L
dose <- matrix(rbinom(P*M, ploidy, 0.5), P, M); rownames(dose) <- sprintf("C%02d", 1:P)
pheno <- setNames(rnorm(P), rownames(dose))
cfg <- list(dosage = dose, n_crosses = 12L, ploidy = ploidy, phenotype = pheno,
            run_qc = TRUE, gain = "usefulness", grm_method = "vanraden")
one <- do.call(ng_polyploid_design_crosses, cfg)
rd <- file.path(tempdir(), "polystage"); unlink(rd, recursive = TRUE)
for (s in ng_poly_cp_stage_order()) out <- ng_poly_run_stage(s, rd, cfg)
staged <- attr(out, "result")
stopifnot(identical(as.data.frame(one), as.data.frame(staged)))
stopifnot(identical(attr(one,"summary")$mean_gain, attr(staged,"summary")$mean_gain))
cat("polyploid_staged_equivalence: staged == one-shot\n")
```

- [ ] **Step 2 — run, verify fail** (`ng_poly_run_stage` not found).
- [ ] **Step 3 — implement `R/50_polyploid_staged_runner.R`.** Stage bodies mirror `R/46:91-133`
  exactly, split at the QC / (fit+score) / allocate / result boundaries:

```r
ng_poly_cp_stage_order <- function() c("qc", "predict", "allocate", "rank")

ng_poly_cp__build_ctx <- function(config) {
  d <- config
  d$gain <- match.arg(d$gain %||% "mean", c("mean","usefulness"))
  d$grm_method <- match.arg(d$grm_method %||% "vanraden", c("vanraden","yang"))
  d$ploidy <- as.integer(d$ploidy %||% 2L)
  # everything after the named formals of ng_polyploid_design_crosses is allocator '...'
  known <- c("dosage","n_crosses","ploidy","effects","phenotype","pairs",
             "max_crosses_per_parent","ridge_seed","run_qc","qc","dominance","gain",
             "selection_prop","double_reduction","grm_method")
  d$alloc_dots <- d[setdiff(names(d), known)]
  d
}

ng_poly_cp__stage_qc <- function(ctx) {
  ctx$geno_raw <- ctx$dosage; ctx$marker_keep <- rep(TRUE, ncol(ctx$dosage))
  if (isTRUE(ctx$run_qc %||% TRUE)) {
    qc_res <- do.call(ng_polyploid_qc, c(list(dosage = ctx$dosage, ploidy = ctx$ploidy), ctx$qc %||% list()))
    ctx$dosage <- qc_res$clean
    ctx$marker_keep <- !qc_res$marker_report$dropped
    ctx$qc <- list(status = if (isTRUE(qc_res$summary$blocked)) "blocker" else "pass",
                   issues = qc_res$summary$issues %||% NULL, marker_report = qc_res$marker_report,
                   putative_duplicates = qc_res$putative_duplicates %||% NULL, raw = qc_res$summary)
  } else ctx$qc <- list(status = "pass")
  ctx
}

ng_poly_cp__stage_predict <- function(ctx) {
  if (identical(ctx$qc$status, "blocker")) return(ctx)
  geno <- ng_polyploid_as_dosage_matrix(ctx$dosage, ploidy = ctx$ploidy, name = "dosage")
  if (isTRUE(ctx$dominance)) {
    if (is.null(ctx$phenotype)) ng_stop("dominance = TRUE needs a phenotype ...")
    y <- suppressWarnings(as.numeric(ctx$phenotype)); names(y) <- names(ctx$phenotype) %||% rownames(geno)
    fit <- ng_polyploid_fit_effects(geno, y[rownames(geno)], ploidy = ctx$ploidy,
                                    model = "additive_dominance", seed = ctx$ridge_seed %||% 1L)
    scores <- ng_polyploid_score_crosses_dominance(fit, geno, pairs = ctx$pairs,
                selection_prop = ctx$selection_prop %||% .1, double_reduction = ctx$double_reduction %||% 0,
                grm_method = ctx$grm_method)
    ctx$gain_col <- if (identical(ctx$gain,"usefulness")) "cross_usefulness" else "cross_mean"
  } else {
    eff <- ctx$effects
    if (!is.null(eff) && length(eff) != ncol(geno)) eff <- as.numeric(eff)[ctx$marker_keep]
    if (is.null(eff)) {
      if (is.null(ctx$phenotype)) ng_stop("supply either effects or phenotype")
      y <- suppressWarnings(as.numeric(ctx$phenotype)); names(y) <- names(ctx$phenotype) %||% rownames(geno)
      eff <- ng_fit_ridge_effects(geno, y[rownames(geno)], seed = ctx$ridge_seed %||% 1L)$beta
    }
    scores <- ng_polyploid_score_crosses(geno, eff, ploidy = ctx$ploidy, pairs = ctx$pairs,
                grm_method = ctx$grm_method, selection_prop = ctx$selection_prop %||% .1,
                double_reduction = ctx$double_reduction %||% 0)
    ctx$gain_col <- if (identical(ctx$gain,"usefulness")) "poly_usefulness" else "poly_mean"
  }
  ctx$scored <- scores
  ctx$parent_kinship <- attr(scores, "parent_kinship")   # explicit ctx field (survives RDS)
  ctx
}

ng_poly_cp__stage_allocate <- function(ctx) {
  if (identical(ctx$qc$status, "blocker")) return(ctx)
  ctx$plan <- do.call(ng_optimize_mating_plan, c(list(ctx$scored, n_crosses = ctx$n_crosses,
    gain_col = ctx$gain_col, parent_kinship = ctx$parent_kinship,
    max_crosses_per_parent = ctx$max_crosses_per_parent %||% 4L), ctx$alloc_dots))
  ctx
}

ng_poly_cp__stage_rank <- function(ctx) ctx    # priority annotation deferred; result assembled below

ng_poly_cp__assemble_result <- function(ctx) {
  plan <- ctx$plan
  s <- attr(plan, "summary"); s$ploidy <- ctx$ploidy; s$poly_gain_col <- ctx$gain_col
  s$dominance <- isTRUE(ctx$dominance); attr(plan, "summary") <- s
  attr(plan, "ploidy") <- ctx$ploidy
  if (isTRUE(ctx$run_qc %||% TRUE)) attr(plan, "qc") <- ctx$qc$raw
  plan
}

ng_poly_cp_pipeline <- list(qc = ng_poly_cp__stage_qc, predict = ng_poly_cp__stage_predict,
                            allocate = ng_poly_cp__stage_allocate, rank = ng_poly_cp__stage_rank)

ng_poly_run_stage <- function(stage, run_dir, config = NULL) {
  order <- ng_poly_cp_stage_order()
  stage <- as.character(stage[[1L]]); if (!stage %in% order) ng_stop("bad poly stage")
  ng_stage_store_init(run_dir)
  if (identical(stage, "qc")) ctx <- ng_poly_cp__build_ctx(config)
  else {
    ctx <- ng_stage_load_ctx(run_dir, order[[match(stage, order) - 1L]])
    if (identical(ctx$qc$status, "blocker")) ng_stop("QC blocker -- resolve before running ", stage)
  }
  ctx <- ng_poly_cp_pipeline[[stage]](ctx)
  ng_stage_save(run_dir, stage, ctx)
  ng_poly_stage_write_json(run_dir, stage, ctx)
  status <- if (identical(stage, "qc")) ctx$qc$status else "done"
  ng_stage__update_manifest(run_dir, stage, status)
  out <- list(status = status, stage = stage,
              files = list.files(file.path(run_dir, "artifacts"), full.names = TRUE))
  if (identical(stage, "rank")) attr(out, "result") <- ng_poly_cp__assemble_result(ctx)
  out
}

ng_poly_stage_write_json <- function(run_dir, stage, ctx) {
  art <- ng_stage_store_init(run_dir)
  payload <- switch(stage,
    qc = list(status = ctx$qc$status, issues = ctx$qc$issues,
              putative_duplicates = list(pairs = ctx$qc$putative_duplicates$pairs %||% NULL),
              marker_report = list(kept = sum(ctx$marker_keep), dropped = sum(!ctx$marker_keep))),
    predict = list(n_candidates = if (is.data.frame(ctx$scored)) nrow(ctx$scored) else NA,
                   poly_metric = ctx$gain_col),
    allocate = { s <- attr(ctx$plan, "summary")
      list(plan_summary = list(n_crosses = nrow(as.data.frame(ctx$plan)),
           mean_gain = s$mean_gain, group_coancestry = s$group_coancestry)) },
    rank = NULL)
  jsonlite::write_json(payload, file.path(art, paste0(stage, ".json")),
    auto_unbox = TRUE, dataframe = "rows", digits = 10, null = "null", na = "null")
}
```

(Verify the exact `ng_polyploid_qc` return shape — `$clean`, `$marker_report$dropped`, `$summary`,
`$putative_duplicates` — with `Rscript -e 'ng_load(); str(ng_polyploid_qc(matrix(...),4))'` and adjust
field names if they differ.)

- [ ] **Step 4 — run golden test → PASS.** Also add a **compute-once** assertion: after running qc+
  predict, stamp `readRDS(predict.rds)`, run allocate, confirm it did not refit (predict artifact
  unchanged mtime / effects identical). And a **QC-blocker** test (out-of-range dosage → qc status
  blocker → predict aborts).
- [ ] **Step 5 — diploid regression:** run the existing diploid staged/equivalence tests → still pass.
- [ ] **Step 6 — `Rscript tools/check_r_package.R`** stays 0 WARNING / 0 NOTE. Bump DESCRIPTION minor.
- [ ] **Step 7 — commit.**

---

### Task 2 (frontend runner): staged-poly dispatch in `inst/app/tools/run_cross_prediction_json.R`

**Files:** Modify the frontend-bundled runner.

- [ ] **Step 1:** extract `run_polyploid_design`'s config→`design_args` block (`~lines 66-146`) into a
  pure helper `poly_design_args(raw)` returning the `design_args` list. `run_polyploid_design` calls
  `do.call(ng_polyploid_design_crosses, poly_design_args(raw))` (behavior unchanged).
- [ ] **Step 2:** add staged dispatch. Where the runner reads `workflow`/`stage`, when
  `workflow=="stage"` AND `raw$prediction_family=="polyploid"` (a new meta key the frontend sets for
  staged autotetraploid), call
  `out <- nextgenCrossDesign::ng_poly_run_stage(raw$stage, dirname(result_path), poly_design_args(raw))`,
  then write the per-stage result JSON (status + the stage JSON payload + on `rank`, the full poly
  result payload via the existing `run_polyploid_design` payload builder fed from `attr(out,"result")`).
- [ ] **Step 3:** verify the installed runner path (`ngcd_res("tools","run_cross_prediction_json.R")`)
  is the one edited; add a runner unit test if the frontend has one, else cover via Task 3 e2e.
- [ ] **Step 4 — commit.**

---

### Task 3 (frontend app): route autotetraploid through the staged cards

**Files:** Modify `R/app.R`, `R/helpers.R` (`ngcd_stage_key_patterns`, `ngcd_stage_summary`).

- [ ] **Step 1:** `run_mode()` — return `"staged"` when `is_poly()` (autotetraploid); keep
  `is_subgenome()` → `"single"`. (One clause change.)
- [ ] **Step 2:** staleness observer (`app.R:1321`) — remove the early `return()` for `is_poly()`
  (keep it for `is_subgenome()`), so poly participates in compute-once. Add poly keys to
  `ngcd_stage_key_patterns` (qc: `ploidy`,`poly_run_qc`,`poly_min_maf`,...; predict:
  `dominance`,`gain`,`grm_method`,`selection_prop`,`double_reduction`,`seed`,`poly_trait_col`;
  allocate: the shared allocation keys).
- [ ] **Step 3:** `ngcd_run_stage()` / `build_poly_params()` — for staged poly set
  `prediction_family="polyploid"` + `workflow="stage"` + `stage`, so the runner dispatches to
  `ng_poly_run_stage`. Ensure the staged poly config carries the poly keys `poly_design_args` needs.
- [ ] **Step 4:** `ngcd_stage_summary()` + `fig_*` renderers — handle the poly stage-JSON shapes
  (qc `marker_report$kept/dropped`; predict `n_candidates`/`poly_metric`; allocate `plan_summary`).
  Reuse diploid branches where shapes match; add poly branches otherwise.
- [ ] **Step 5 — testServer:** autotetraploid renders 3 stepped cards; a poly allocation-only input
  change leaves qc/predict `done` (compute-once). Bump `inst/BACKEND_VERSION` to the new backend.
- [ ] **Step 6 — commit.**

---

### Task 4: end-to-end + ship

- [ ] Backend: build tarball, install; frontend: `R CMD check`/testthat green; run a real staged
  autotetraploid demo through the app (QC → Fit & score → Allocate, each its own run + figure).
- [ ] Golden e2e: staged poly result JSON == one-shot poly result JSON for the demo dataset.
- [ ] PRs on both repos, CI green (repo public); merge. Rebuild `ngcd-workbench` image + tarballs.

## Self-review

- Spec coverage: backend staged runner (T1), frontend dispatch (T2), app routing + compute-once (T3),
  e2e/ship (T4). ✓  Byte-identical + compute-once + diploid-untouched + attribute-persistence all have
  explicit tests. Placeholder scan: the one unknown (`ng_polyploid_qc` field names) has a verify
  command. Subgenome explicitly deferred to Phase 2.
