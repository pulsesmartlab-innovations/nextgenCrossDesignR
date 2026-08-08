# ---------------------------------------------------------------------------
# Staged (resumable) runner for the autotetraploid analytic design.
#
# ng_polyploid_design_crosses (R/46) runs QC -> effects -> scoring -> allocation
# in one call. This file decomposes it into stage functions over a shared `ctx`
# so a caller can run ONE stage per subprocess call, cache each stage's ctx to a
# run-dir, and resume from any completed stage -- the same compute-once model as
# the diploid ng_run_stage (R/45), whose generic persistence helpers
# (ng_stage_save / ng_stage_load_ctx / ng_stage__update_manifest) are reused
# here. Driving all stages in order via ng_poly_run_stage() is byte-identical to
# calling ng_polyploid_design_crosses() directly (see tests/polyploid_staged_equivalence.R).
#
# Single-trait, so the `index` stage is a NO-OP passthrough (kept only so the
# poly stage order matches the diploid qc/predict/index/allocate/rank -- the
# frontend's 5-stage pipeline machinery then drives poly unchanged). Poly QC
# CLEANS the dosage and proceeds -- it never blocks -- exactly as the one-shot
# design does, so equivalence holds. The QC-blocker gate is kept structurally
# (harmless) for parity with the diploid runner.
# ---------------------------------------------------------------------------

ng_poly_cp_stage_order <- function() c("qc", "predict", "index", "allocate", "rank")

# Normalize the poly config; everything beyond the named ng_polyploid_design_crosses
# formals is treated as allocator `...` (method/strategy/target_coancestry/...).
ng_poly_cp__build_ctx <- function(config) {
  d <- as.list(config)
  d$gain <- match.arg(d$gain %||% "mean", c("mean", "usefulness"))
  d$grm_method <- match.arg(d$grm_method %||% "vanraden", c("vanraden", "yang"))
  d$ploidy <- as.integer(d$ploidy %||% 2L)
  known <- c("dosage", "n_crosses", "ploidy", "effects", "phenotype", "pairs",
             "max_crosses_per_parent", "ridge_seed", "run_qc", "qc", "dominance", "gain",
             "selection_prop", "double_reduction", "grm_method")
  d$alloc_dots <- d[setdiff(names(d), c(known, "alloc_dots"))]
  d
}

# --- stage: QC (clean dosage; never blocks -- matches the one-shot) ----------
ng_poly_cp__stage_qc <- function(ctx) {
  ctx$marker_keep <- rep(TRUE, ncol(ctx$dosage))
  if (isTRUE(ctx$run_qc %||% TRUE)) {
    qc_res <- do.call(ng_polyploid_qc, c(list(dosage = ctx$dosage, ploidy = ctx$ploidy),
                                         if (is.list(ctx$qc)) ctx$qc else list()))
    ctx$dosage <- qc_res$clean
    ctx$marker_keep <- !qc_res$marker_report$dropped
    dropped <- (qc_res$summary$n_markers_dropped %||% 0L) + (qc_res$summary$n_samples_dropped %||% 0L)
    ctx$qc <- list(status = if (dropped > 0) "warning" else "pass",
                   summary = qc_res$summary,
                   marker_report = qc_res$marker_report,
                   duplicate_samples = qc_res$duplicate_samples)
  } else {
    ctx$qc <- list(status = "pass", summary = NULL)
  }
  ctx
}

# --- stage: predict = fit effects + score crosses ----------------------------
ng_poly_cp__stage_predict <- function(ctx) {
  if (identical(ctx$qc$status, "blocker")) return(ctx)
  geno <- ng_polyploid_as_dosage_matrix(ctx$dosage, ploidy = ctx$ploidy, name = "dosage")
  seed <- ctx$ridge_seed %||% 1L
  if (isTRUE(ctx$dominance)) {
    if (is.null(ctx$phenotype)) ng_stop("dominance = TRUE needs a phenotype to estimate dominance effects")
    y <- suppressWarnings(as.numeric(ctx$phenotype))
    names(y) <- if (!is.null(names(ctx$phenotype))) names(ctx$phenotype) else rownames(geno)
    fit <- ng_polyploid_fit_effects(geno, y[rownames(geno)], ploidy = ctx$ploidy,
                                    model = "additive_dominance", seed = seed)
    scores <- ng_polyploid_score_crosses_dominance(fit, geno, pairs = ctx$pairs,
                selection_prop = ctx$selection_prop %||% 0.10,
                double_reduction = ctx$double_reduction %||% 0, grm_method = ctx$grm_method)
    ctx$gain_col <- if (identical(ctx$gain, "usefulness")) "cross_usefulness" else "cross_mean"
  } else {
    eff <- ctx$effects
    if (!is.null(eff) && length(eff) != ncol(geno)) eff <- as.numeric(eff)[ctx$marker_keep]
    if (is.null(eff)) {
      if (is.null(ctx$phenotype)) ng_stop("supply either effects (per marker) or phenotype (per parent)")
      y <- suppressWarnings(as.numeric(ctx$phenotype))
      names(y) <- if (!is.null(names(ctx$phenotype))) names(ctx$phenotype) else rownames(geno)
      eff <- ng_fit_ridge_effects(geno, y[rownames(geno)], seed = seed)$beta
    }
    scores <- ng_polyploid_score_crosses(geno, eff, ploidy = ctx$ploidy, pairs = ctx$pairs,
                grm_method = ctx$grm_method, selection_prop = ctx$selection_prop %||% 0.10,
                double_reduction = ctx$double_reduction %||% 0)
    ctx$gain_col <- if (identical(ctx$gain, "usefulness")) "poly_usefulness" else "poly_mean"
  }
  ctx$scored <- scores
  # parent_kinship rides as an attribute on `scores`; persist it as an explicit
  # ctx field so it survives the RDS boundary into the allocate stage.
  ctx$parent_kinship <- attr(scores, "parent_kinship")
  ctx
}

# --- stage: allocate ---------------------------------------------------------
ng_poly_cp__stage_allocate <- function(ctx) {
  if (identical(ctx$qc$status, "blocker")) return(ctx)
  ctx$plan <- do.call(ng_optimize_mating_plan, c(
    list(ctx$scored, n_crosses = ctx$n_crosses, gain_col = ctx$gain_col,
         parent_kinship = ctx$parent_kinship,
         max_crosses_per_parent = ctx$max_crosses_per_parent %||% 4L),
    ctx$alloc_dots))
  ctx
}

# --- stage: index (no-op for single-trait poly; keeps 5-stage parity) --------
ng_poly_cp__stage_index <- function(ctx) ctx

# --- stage: rank (priority annotation deferred; result assembled below) ------
ng_poly_cp__stage_rank <- function(ctx) ctx

# Reassemble the exact object ng_polyploid_design_crosses returns.
ng_poly_cp__assemble_result <- function(ctx) {
  plan <- ctx$plan
  s <- attr(plan, "summary")
  s$ploidy <- as.integer(ctx$ploidy); s$poly_gain_col <- ctx$gain_col; s$dominance <- isTRUE(ctx$dominance)
  attr(plan, "summary") <- s
  attr(plan, "ploidy") <- as.integer(ctx$ploidy)
  if (isTRUE(ctx$run_qc %||% TRUE)) attr(plan, "qc") <- ctx$qc$summary
  plan
}

ng_poly_cp_pipeline <- list(
  qc = ng_poly_cp__stage_qc, predict = ng_poly_cp__stage_predict,
  index = ng_poly_cp__stage_index,
  allocate = ng_poly_cp__stage_allocate, rank = ng_poly_cp__stage_rank)

#' @export
ng_poly_run_stage <- function(stage, run_dir, config = NULL) {
  order <- ng_poly_cp_stage_order()
  stage <- as.character(stage[[1L]])
  if (!(stage %in% order)) ng_stop("stage must be one of: ", paste(order, collapse = ", "))
  ng_stage_store_init(run_dir)
  if (identical(stage, "qc")) {
    ctx <- ng_poly_cp__build_ctx(config)
  } else {
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

# Per-stage JSON summary for the frontend cards (poly shapes).
ng_poly_stage_write_json <- function(run_dir, stage, ctx) {
  art <- ng_stage_store_init(run_dir)
  payload <- switch(stage,
    qc = list(status = ctx$qc$status,
              marker_report = list(kept = sum(ctx$marker_keep), dropped = sum(!ctx$marker_keep)),
              duplicate_samples = ctx$qc$duplicate_samples %||% NULL),
    predict = list(n_candidates = if (is.data.frame(ctx$scored)) nrow(ctx$scored) else NA_integer_,
                   poly_metric = ctx$gain_col),
    allocate = {
      s <- attr(ctx$plan, "summary")
      list(plan_summary = list(n_crosses = nrow(as.data.frame(ctx$plan)),
                               mean_gain = s$mean_gain, group_coancestry = s$group_coancestry))
    },
    rank = NULL)
  jsonlite::write_json(payload, file.path(art, paste0(stage, ".json")),
    auto_unbox = TRUE, dataframe = "rows", digits = 10, null = "null", na = "null")
}
