# ---------------------------------------------------------------------------
# Resumable staged runner + run-dir artifact store.
#
# ng_run_cross_prediction (R/39) is decomposed into stage functions over a
# shared `ctx` list (ng_cp_pipeline / ng_cp_stage_order / ng_cp__build_ctx /
# ng_cp__assemble_result). This file adds a resumable entry point that runs
# ONE stage per call, persisting the accumulated ctx to a run-dir so a caller
# (e.g. a frontend driving one subprocess call per stage) can compute each
# stage once and resume from any completed stage. Driving all stages in order
# via ng_run_stage() is byte-identical to calling ng_run_cross_prediction()
# directly, because both paths call the exact same ng_cp_pipeline functions.
# ---------------------------------------------------------------------------

#' @export
ng_run_stage <- function(stage, run_dir, config = NULL) {
  stage <- as.character(stage[[1L]])
  stage_order <- ng_cp_stage_order()
  if (!(stage %in% stage_order)) {
    ng_stop("stage must be one of: ", paste(stage_order, collapse = ", "))
  }
  ng_stage_store_init(run_dir)

  if (identical(stage, "qc")) {
    ctx <- ng_cp__build_ctx(config)
  } else {
    upstream <- stage_order[[match(stage, stage_order) - 1L]]
    ctx <- ng_stage_load_ctx(run_dir, upstream)
    # Gate: a persisted QC blocker must be resolved before any downstream
    # stage runs. ctx$qc is set by the qc stage and threaded through every
    # later ctx snapshot, so this reads the QC status exactly as it was
    # persisted to disk by the qc stage (via the qc.rds artifact).
    if (identical(ctx$qc$status, "blocker")) {
      ng_stop("QC blocker -- resolve before running ", stage)
    }
  }

  ctx <- ng_cp_pipeline[[stage]](ctx)
  ng_stage_save(run_dir, stage, ctx)
  ng_stage_write_json(run_dir, stage, ctx)

  status <- if (identical(stage, "qc")) ctx$qc$status else "done"
  ng_stage__update_manifest(run_dir, stage, status)

  files <- list.files(file.path(run_dir, "artifacts"), full.names = TRUE)
  out <- list(status = status, stage = stage, files = files)
  if (identical(stage, "rank")) {
    attr(out, "result") <- ng_cp__assemble_result(ctx)
  }
  out
}

ng_stage_store_init <- function(run_dir) {
  run_dir <- as.character(run_dir[[1L]])
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  art_dir <- file.path(run_dir, "artifacts")
  dir.create(art_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(art_dir)
}

ng_stage_save <- function(run_dir, stage, ctx) {
  art_dir <- ng_stage_store_init(run_dir)
  saveRDS(ctx, file.path(art_dir, paste0(stage, ".rds")))
  invisible(TRUE)
}

ng_stage_load_ctx <- function(run_dir, stage) {
  path <- file.path(run_dir, "artifacts", paste0(stage, ".rds"))
  if (!file.exists(path)) {
    ng_stop("staged pipeline: artifact for stage '", stage, "' was not found at ", path,
            " -- run that stage before resuming from it")
  }
  readRDS(path)
}

ng_stage_status <- function(run_dir) {
  path <- file.path(run_dir, "artifacts", "manifest.json")
  if (!file.exists(path)) {
    return(list(schema = "ng_stage_manifest.v1", stages = list()))
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

ng_stage__update_manifest <- function(run_dir, stage, status) {
  art_dir <- ng_stage_store_init(run_dir)
  manifest <- ng_stage_status(run_dir)
  if (is.null(manifest$schema)) manifest$schema <- "ng_stage_manifest.v1"
  if (is.null(manifest$stages)) manifest$stages <- list()
  manifest$stages[[stage]] <- list(
    status = status,
    mtime = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )
  jsonlite::write_json(
    manifest, file.path(art_dir, "manifest.json"),
    auto_unbox = TRUE, dataframe = "rows", digits = 10, null = "null", na = "null"
  )
  invisible(manifest)
}

ng_stage_write_json <- function(run_dir, stage, ctx) {
  art_dir <- ng_stage_store_init(run_dir)
  payload <- switch(
    stage,
    qc = list(
      status = ctx$qc$status,
      issues = ctx$qc$issues,
      putative_duplicates = list(
        pairs = if (!is.null(ctx$qc$putative_duplicates)) ctx$qc$putative_duplicates$pairs else NULL
      )
    ),
    predict = list(
      effect_summary = if (is.data.frame(ctx$effect_summary)) ctx$effect_summary else do.call(rbind, ctx$effect_summary),
      ld_pruning_report = ctx$ld_pruning_report,
      n_candidates = if (!is.null(ctx$cross_table)) as.integer(nrow(ctx$cross_table)) else NA_integer_
    ),
    index = list(
      multi_trait_score = if (!is.null(ctx$scored_crosses)) ctx$scored_crosses$multi_trait_score else NULL,
      objective = list(method = ctx$objective$method)
    ),
    allocate = list(
      plan_summary = attr(ctx$plan, "summary"),
      parent_use = ng_run_cp_parent_use(ctx$plan)
    ),
    rank = NULL
  )
  if (is.null(payload)) return(invisible(NULL))
  out_path <- file.path(art_dir, paste0(stage, ".json"))
  jsonlite::write_json(
    payload, out_path,
    auto_unbox = TRUE, dataframe = "rows", digits = 10, null = "null", na = "null"
  )
  invisible(out_path)
}
