# Run many INDEPENDENT single-trait analyses as one batch.
#
# A breeder with 17 traits who wants 17 separate crossing plans -- not one multi-trait
# index -- had no way to ask for that. They set the trait, launched, waited, and repeated
# 17 times. Every one of those runs rebuilt the same QC, the same cleaned genotypes, the
# same duplicate scan, the same LD pruning and the same GRM, none of which depend on which
# trait is being scored.
#
# WHY A BATCH IS THE SAME ANSWER AS 17 RUNS, not an approximation of it. Two properties
# were built into the runner deliberately, for other reasons, and together they are
# exactly what a fan-out needs:
#
#   * the candidate parent set is fixed at QC (geno n pheno). A trait's missing phenotypes
#     enlarge or shrink only its TRAINING augmentation, never the candidates or the pairs.
#   * the ridge-lambda CV fold partition is keyed on the run seed, not the trait's row
#     position (R/39, 0.28.0), and the posterior draw stream is keyed on the trait's NAME.
#
# So trait k scored alone and trait k scored inside a batch see the same parents, the same
# folds and the same draws. tests/batch_matches_individual_runs.R asserts it rather than
# trusting it.
#
# WHAT IS SHARED AND WHAT IS NOT. Everything up to and including ng_cp__predict_prologue()
# is computed once in the parent and handed to every job: cleaned and aligned genotypes,
# the marker map, LD pruning, the training-set alignment, the pair table and the GRM. From
# the per-trait ridge fit onward each job is on its own. That split is not an optimisation
# guess -- the per-cross recombination kernel is O(m^2) per pair and takes the trait's beta,
# so it is irreducibly per-trait and it is where the hours go.

# The run-result envelope, in the package rather than in tools/, because the batch writes
# one per job and a second copy of this whitelist would drift from the first. The JSON
# runner is now a caller.
ng_run_result_envelope <- function(r, warnings = character(0), generated_at = NULL,
                                   package_version = NULL) {
  audit <- r$input_match_audit
  audit$marker_order <- NULL   # drop the (potentially huge) marker-name vector; count kept
  list(
    schema          = "ng_run_result.v1",
    status          = "ok",
    ok              = TRUE,               # frontend keys success on `ok`
    generated_at    = generated_at,
    package_version = package_version,
    prediction_mode = r$prediction_mode,
    settings          = r$settings,
    input_match_audit = audit,
    qc = list(status = r$qc$status, counts = r$qc$counts,
              tables = r$qc$tables, issues = r$qc$issues),
    effect_summary  = r$effect_summary,
    trait_direction = r$trait_direction,
    objective       = if (!is.null(r$objective)) r$objective$diagnostics else NULL,
    plan_summary    = r$plan_summary,
    constraint_diagnostics = r$constraint_diagnostics,
    priority_risk_diagnostics = r$priority_risk_diagnostics,
    trait_check_reference = r$trait_check_reference,
    candidate_crosses = r$candidate_crosses,
    selected_crosses  = r$selected_crosses,
    ld_pruning_report = r$ld_pruning_report,
    warnings          = warnings,
    output_files      = r$output_files
  )
}

# Settings a job may NOT override, because the batch already spent them: they are consumed
# in stage_qc or in the shared predict prologue, so a job that changed one would be scored
# against inputs it did not ask for. Rejecting is the honest response -- silently ignoring
# the override would hand back numbers that quietly contradict the request, and silently
# honouring it would mean the "shared" artefact was not shared at all.
#
# A job wanting different QC belongs in its own batch. That is a real limitation and it is
# stated, not discovered.
ng_cp__batch_shared_keys <- c(
  # inputs and their column mappings -- read in stage_qc
  "genotype", "genotype_file", "genotype_id_col",
  "phenotype", "phenotype_file", "phenotype_id_col",
  "marker_map", "map_file", "map_marker_col", "map_chr_col", "map_pos_col",
  "map_pos_cm_col", "map_pos_bp_col", "map_position_unit", "map_pos_cm_divisor",
  "bp_per_cm",
  "trait_direction", "direction_file", "direction_trait_col", "direction_column_col",
  "direction_direction_col", "direction_value_kind_col", "id_col",
  # genotype QC
  "duplicate_action", "duplicate_threshold", "duplicate_maf_min",
  "duplicate_max_missing_prop", "duplicate_min_compared_markers",
  # the shared predict prologue
  "ld_pruning", "ld_window", "ld_r2_threshold", "ld_maf_threshold", "ld_ploidy",
  "ld_backend", "training_genotype", "training_phenotype", "training_genotype_file",
  "training_phenotype_file", "training_genotype_id_col", "training_phenotype_id_col",
  "grm_method"
)

# Restrict a shared context to one job's traits and apply its overrides.
ng_cp__batch_apply_job <- function(ctx, job, output_root) {
  overrides <- job$overrides %||% list()
  bad <- intersect(names(overrides), ng_cp__batch_shared_keys)
  if (length(bad)) {
    ng_stop("job '", job$id, "' overrides setting(s) the batch has already computed: ",
            paste(bad, collapse = ", "),
            ". These are consumed by quality control or by the shared marker-effect setup, ",
            "which run once for the whole batch. Run such a job as its own batch.")
  }
  if (length(overrides)) ctx <- do.call(ng_ctx_put, c(list(ctx), overrides))

  traits <- job$traits
  spec <- ctx$trait_spec
  keep <- spec$trait %in% traits | spec$column %in% traits
  if (!any(keep)) {
    ng_stop("job '", job$id, "' names trait(s) absent from the direction table: ",
            paste(traits, collapse = ", "))
  }
  job_dir <- file.path(output_root, job$id)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  ng_ctx_put(ctx,
             trait_spec = spec[keep, , drop = FALSE],
             traits_to_use = traits,
             output_dir = job_dir,
             write_outputs = ctx$write_outputs %||% TRUE,
             write_figures = ctx$write_figures %||% FALSE)
}

# One job, start to finish, inside a daemon.
#
# It RETURNS its conditions rather than signalling them. A warning() raised in a worker
# never reaches the parent's calling handlers, so it would never reach the JSON envelope --
# R/52_advisories.R calls a silent advisory channel "the exact defect class this work
# exists to remove". The parent re-emits what comes back, attributed to the trait.
#
# It also returns a SUMMARY, not the result object: seventeen full results would be
# serialised back through the daemon socket for no reason, when each has already written
# itself to disk.
ng_cp__batch_run_one <- function(job, shared_path, output_root,
                                 generated_at = NULL, package_version = NULL) {
  t0 <- Sys.time()
  seen <- character(0)
  value <- withCallingHandlers(
    tryCatch({
      shared <- get0(".ngcd_batch_shared", envir = globalenv(), inherits = FALSE)
      if (is.null(shared)) shared <- readRDS(shared_path)
      ctx <- ng_cp__batch_apply_job(shared, job, output_root)
      for (s in setdiff(ng_cp_stage_order(), "qc")) ctx <- ng_cp_pipeline[[s]](ctx)
      r <- ng_cp__assemble_result(ctx)
      job_dir <- file.path(output_root, job$id)
      env <- ng_run_result_envelope(r, warnings = seen, generated_at = generated_at,
                                    package_version = package_version)
      result_path <- file.path(job_dir, "result.json")
      ng_write_result_json(env, result_path)
      es <- r$effect_summary
      list(
        status = "ok",
        traits = job$traits,
        cv_predictive_r2 = if (!is.null(es) && "cv_predictive_r2" %in% names(es))
          as.numeric(es$cv_predictive_r2) else NA_real_,
        mean_source = if (!is.null(es) && "mean_source" %in% names(es))
          as.character(es$mean_source) else NA_character_,
        n_selected = if (is.data.frame(r$selected_crosses)) nrow(r$selected_crosses) else NA_integer_,
        result_json = result_path,
        output_files = r$output_files
      )
    }, error = function(e) list(status = "error", traits = job$traits,
                                error_message = conditionMessage(e))),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  value$id <- job$id
  value$warnings <- seen
  value$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  value
}

# Non-finite (Inf/NaN) -> null, so the JSON is valid and consumable everywhere. Moved here
# from tools/run_cross_prediction_json.R alongside the envelope: a batch writes one result
# per job, and a second copy of these rules would drift from the single-run path that the
# frontend's contract was written against.
ng_run_result_sanitize <- function(x) {
  if (is.data.frame(x)) {
    for (j in seq_along(x)) if (is.numeric(x[[j]])) x[[j]][!is.finite(x[[j]])] <- NA
    return(x)
  }
  if (is.list(x)) return(lapply(x, ng_run_result_sanitize))
  if (is.numeric(x)) x[!is.finite(x)] <- NA
  x
}

# Write a result envelope. The options are NOT free parameters -- dataframe = "rows" and
# digits = 8 are the shape the frontend parses and the contract documents, so a batch job's
# result.json must be byte-comparable with a single run's.
#
# Written to a .part file and renamed, because rename is atomic: a frontend polling the
# output directory of a running batch would otherwise eventually parse a truncated file.
ng_write_result_json <- function(env, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".part")
  jsonlite::write_json(ng_run_result_sanitize(env), tmp, auto_unbox = TRUE, na = "null",
                       null = "null", dataframe = "rows", pretty = TRUE, digits = 8)
  file.rename(tmp, path)
  path
}

# Same atomicity for the batch manifest, which a caller may poll while jobs are still
# running. The manifest is not the run contract, so it keeps full precision.
ng_write_json_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".part")
  jsonlite::write_json(ng_run_result_sanitize(x), tmp, auto_unbox = TRUE, na = "null",
                       null = "null", dataframe = "rows", pretty = TRUE)
  file.rename(tmp, path)
  path
}

# mirai::everywhere(expr, ...) injects its named arguments into the expression at daemon
# evaluation time, so codetools cannot see those bindings any more than it can see the ones
# list2env() creates for the ng_cp__* stages. Same convention, declared here beside its use.
utils::globalVariables(c(".root", ".shared", ".use_cpp", ".rng"))

ng_run_cross_prediction_batch <- function(config,
                                          jobs = NULL,
                                          output_root,
                                          batch_workers = NULL,
                                          memory_budget_bytes = NULL,
                                          generated_at = NULL,
                                          package_version = NULL) {
  if (!is.list(config)) ng_stop("config must be a list of ng_run_cross_prediction() arguments")
  if (!length(output_root) || !nzchar(output_root)) ng_stop("output_root is required")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (is.null(generated_at)) generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  # Back-fill every missing formal from the runner's own defaults. ng_cp__build_ctx()
  # match.arg()s every enum, and an ABSENT enum arrives as its whole choices vector rather
  # than as its default -- the same trap tools/run_cross_prediction_json.R documents.
  fm <- formals(ng_run_cross_prediction)
  for (nm in setdiff(names(fm), names(config))) {
    config[nm] <- list(tryCatch(eval(fm[[nm]], envir = environment()), error = function(e) NULL))
  }

  # Shared work, once: build, QC, then the trait-independent predict prologue.
  ctx <- ng_cp__build_ctx(config)
  ctx <- ng_cp__stage_qc(ctx)
  if (identical(ctx$qc$status, "blocker")) {
    ng_stop("QC blocker -- resolve before running a batch: ",
            paste(utils::head(ctx$qc$issues$message, 3L), collapse = "; "))
  }
  ctx <- ng_ctx_put(ctx, predict_prologue = ng_cp__predict_prologue(ctx))

  # Default to one job per trait in the direction table -- the case this exists for.
  if (is.null(jobs)) {
    jobs <- lapply(ctx$trait_spec$trait, function(tr) list(id = tr, traits = tr))
    names(jobs) <- ctx$trait_spec$trait
  }
  jobs <- lapply(seq_along(jobs), function(i) {
    j <- jobs[[i]]
    if (is.null(j$id)) j$id <- names(jobs)[[i]] %||% as.character(i)
    if (is.null(j$traits)) j$traits <- j$id
    j
  })
  ids <- vapply(jobs, function(j) as.character(j$id), character(1))
  if (anyDuplicated(ids)) ng_stop("job ids must be unique; duplicated: ",
                                  paste(unique(ids[duplicated(ids)]), collapse = ", "))
  names(jobs) <- ids

  # Fail before spending hours, not after: validate every job's overrides up front.
  for (j in jobs) invisible(ng_cp__batch_apply_job(ctx, j, output_root))

  shared_dir <- file.path(output_root, "_shared")
  dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
  shared_path <- file.path(shared_dir, "shared.rds")
  saveRDS(ctx, shared_path)   # on disk, so a killed batch resumes without redoing QC

  pro <- ctx$predict_prologue
  per_job <- ng_batch_job_bytes(
    n_parents = nrow(pro$geno), n_markers = ncol(pro$geno),
    dense_beta_cov = identical(ctx$method_varPMV, "full_posterior") || ncol(pro$geno) <= 6000L)
  if (is.null(batch_workers)) {
    batch_workers <- ng_batch_worker_count(length(jobs), per_job, memory_budget_bytes)
  } else {
    batch_workers <- structure(max(1L, as.integer(batch_workers)), basis = "caller-specified")
  }
  workers <- min(as.integer(batch_workers), length(jobs))
  basis <- attr(batch_workers, "basis")

  results <- ng_cp__batch_dispatch(jobs, workers, shared_path, output_root,
                                   generated_at, package_version, use_cpp = ctx$use_cpp)

  # Re-emit what the workers could only return. Attributed to the trait, because an
  # unattributed advisory in a 17-job batch tells the breeder nothing actionable.
  for (r in results) {
    for (w in r$warnings %||% character(0)) {
      warning(sprintf("[%s] %s", r$id, w), call. = FALSE)
    }
  }

  manifest <- list(
    schema = "ng_batch_manifest.v1",
    generated_at = generated_at,
    package_version = package_version,
    output_root = output_root,
    n_jobs = length(jobs),
    workers = workers,
    worker_basis = basis,
    per_job_bytes_estimate = per_job,
    jobs = lapply(results, function(r) {
      r$warnings <- NULL
      r
    })
  )
  manifest_path <- ng_write_json_atomic(manifest, file.path(output_root, "manifest.json"))
  structure(list(manifest = manifest, manifest_path = manifest_path, jobs = results),
            class = c("ng_cross_prediction_batch", "list"))
}

# The fan-out itself. Kept separate so the parallel backend is one swappable seam rather
# than a decision spread through the runner.
#
# mirai rather than the trait-level ng_run_cp_apply(): daemons are ordinary processes on
# every platform, so Windows and Linux take the SAME path (no fork-vs-PSOCK split), one
# job's crash cannot take the batch down with it, and mirai_map returns per-job errors so
# only the failures need re-running. At hours per job its dispatch cost is irrelevant --
# what it buys here is isolation and uniformity, not speed.
ng_cp__batch_dispatch <- function(jobs, workers, shared_path, output_root,
                                  generated_at, package_version, use_cpp = FALSE) {
  if (workers <= 1L || !requireNamespace("mirai", quietly = TRUE)) {
    if (workers > 1L) {
      warning("mirai is not available; running the batch one job at a time.", call. = FALSE)
    }
    return(lapply(jobs, ng_cp__batch_run_one, shared_path = shared_path,
                  output_root = output_root, generated_at = generated_at,
                  package_version = package_version))
  }
  # In a dev tree the package is sourced, not installed, so a daemon must load it the same
  # way. Detecting the namespace is what tells the two apart.
  #
  # The root is resolved here rather than via tools/ng_find_package_root.R, which is not part
  # of the package and so is not on a daemon's path -- and not reachable from an installed
  # library either.
  dev_root <- ""
  if (!isNamespaceLoaded("nextgenCrossDesign")) {
    d <- normalizePath(getwd(), mustWork = FALSE)
    repeat {
      if (file.exists(file.path(d, "DESCRIPTION")) &&
          file.exists(file.path(d, "R", "load.R"))) break
      up <- dirname(d)
      if (identical(up, d)) { d <- ""; break }
      d <- up
    }
    dev_root <- d
    if (!nzchar(dev_root)) {
      warning("the package is neither installed nor locatable as a source tree, so batch ",
              "daemons could not load it; running the batch one job at a time.", call. = FALSE)
      return(lapply(jobs, ng_cp__batch_run_one, shared_path = shared_path,
                    output_root = output_root, generated_at = generated_at,
                    package_version = package_version))
    }
  }

  mirai::daemons(workers)
  on.exit(mirai::daemons(0), add = TRUE)
  # THE RNG KIND MUST MATCH THE PARENT'S, and this is not a detail.
  #
  # mirai sets daemons to L'Ecuyer-CMRG so that concurrent workers draw independent
  # streams. That is the right default for code which relies on ambient randomness -- and
  # the wrong one here, because set.seed(s) does not mean the same thing under two
  # different generators. This package seeds every stochastic step explicitly from a
  # derived seed (the ridge-lambda CV fold partition from the run seed; posterior draws
  # from a hash of the trait name), so it needs no stream independence at all; what it
  # needs is that set.seed(s) produce the same folds everywhere.
  #
  # Left unpinned, a batch silently selects a different ridge lambda from the identical
  # data -- measured at 8.86 against the parent's 1.62 on the test fixture -- and every
  # marker-derived number moves with it. Nothing errors; the plan is simply a different
  # plan. tests/batch_is_order_and_worker_independent.R is the regression guard.
  parent_rng <- RNGkind()
  mirai::everywhere({
    RNGkind(kind = .rng[[1L]], normal.kind = .rng[[2L]], sample.kind = .rng[[3L]])
    if (nzchar(.root)) {
      source(file.path(.root, "R", "load.R")); ng_load(.root, use_cpp = .use_cpp, verbose = FALSE)
    } else {
      library(nextgenCrossDesign)
    }
    # Each job's linear algebra is already the machine's whole workload; letting BLAS also
    # thread inside every daemon oversubscribes the box N-fold.
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
      try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
    }
    # Read the shared context ONCE per daemon rather than once per job: it carries the
    # genotype matrix, and a daemon may serve several jobs in turn.
    assign(".ngcd_batch_shared", readRDS(.shared), envir = globalenv())
  }, .root = dev_root, .use_cpp = isTRUE(use_cpp), .shared = shared_path, .rng = parent_rng)

  m <- mirai::mirai_map(
    jobs, ng_cp__batch_run_one,
    .args = list(shared_path = shared_path, output_root = output_root,
                 generated_at = generated_at, package_version = package_version))
  out <- m[]
  # A daemon that died leaves a miraiError in place of the job's result. Convert it to the
  # same failure shape a caught error produces, so one dead worker cannot make the manifest
  # a different shape from a failed run.
  lapply(seq_along(out), function(i) {
    r <- out[[i]]
    if (inherits(r, "miraiError") || inherits(r, "errorValue")) {
      return(list(id = jobs[[i]]$id, traits = jobs[[i]]$traits, status = "error",
                  error_message = paste("worker failed:", conditionMessage(r)),
                  warnings = character(0), elapsed_sec = NA_real_))
    }
    r
  })
}
