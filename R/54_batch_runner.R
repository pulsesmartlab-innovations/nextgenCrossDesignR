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
    # Same pinning as the manifest's plural fields: one warning is an array of one, not a
    # bare string. The frontend renders this list, and a list that is sometimes a string is a
    # bug waiting for the first single-warning run.
    warnings          = I(as.character(warnings)),
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

# The content-addressed artefact key must cover MORE than a job may not override.
#
# Those look like the same question -- "what has the batch already spent?" -- but they are
# not. ng_cp__batch_shared_keys answers "what may a job not override", and it deliberately
# EXCLUDES traits_to_use: choosing which traits a job scores is the entire reason a job
# exists, and ng_cp__batch_apply_job() subsets ctx$trait_spec per job precisely so a job can
# override it.
#
# ng_shared_artifact_key() (R/55) answers a different question: "what does the persisted
# artefact depend on". ng_cp__stage_qc() builds trait_spec from direction, traits_to_use and
# trait_weights (ng_run_cp_trait_spec(), R/39) -- or, under prediction_mode =
# "index_as_trait", from index_col and index_direction (ng_run_cp_index_spec(), R/39) -- and
# trait_spec is one of the fields ng_shared_artifact_fields (R/55) persists into the
# artefact. So the artefact's content depends on all five even though a job may freely
# override them.
#
# Keying the artefact on ng_cp__batch_shared_keys alone -- the bug this set exists to
# prevent -- would let two batches that differ ONLY in traits_to_use collide on the same
# key: the second batch would silently reuse the first batch's artefact, inheriting the
# FIRST batch's trait_spec (and hence its traits) instead of computing its own. No error, no
# warning -- a crossing plan scored against the wrong traits. Because that failure is
# artefact-key-shaped, not override-guard-shaped, it needs its own set rather than a
# broadened ng_cp__batch_shared_keys: broadening the shared-keys set would also make
# traits_to_use unoverridable, which breaks the one thing a job is for.
#
# So there are two sets, deliberately different, and the artefact set is a superset of the
# override-guard set.
ng_cp__batch_artifact_keys <- c(
  ng_cp__batch_shared_keys,
  "traits_to_use", "trait_weights", "prediction_mode", "index_col", "index_direction"
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
                                 generated_at = NULL, package_version = NULL,
                                 config = NULL) {
  t0 <- Sys.time()
  seen <- character(0)
  value <- withCallingHandlers(
    tryCatch({
      ng_job_trait_status_write(output_root, job$id, "running")
      # Rebuild the context the same way the parent does: THIS batch's config, with the
      # artefact's shared work laid over it. The artefact is never the config -- see the
      # comment on ng_shared_artifact_fields in R/55. The expensive half (a genotype matrix
      # and a GRM) is cached per DAEMON by the everywhere() block below; ng_cp__build_ctx()
      # is match.arg() bookkeeping and is cheap enough to redo per job, which is what keeps
      # it inside this withCallingHandlers -- a warning raised here reaches the breeder
      # attributed to the trait rather than dying unheard in a daemon.
      fields <- get0(".ngcd_batch_artifact", envir = globalenv(), inherits = FALSE)
      if (is.null(fields)) {
        fields <- ng_shared_artifact_read(shared_path)
        if (is.null(fields)) {
          ng_stop("the shared artefact at ", shared_path, " is missing, unreadable, or was ",
                  "written by a different version of this package. Re-run the batch without ",
                  "resume so quality control is recomputed.")
        }
      }
      cfg <- get0(".ngcd_batch_config", envir = globalenv(), inherits = FALSE)
      if (is.null(cfg)) cfg <- config
      if (is.null(cfg)) ng_stop("no batch configuration reached the worker for job '", job$id, "'")
      ctx <- ng_cp__batch_ctx(cfg, fields)
      ctx <- ng_cp__batch_apply_job(ctx, job, output_root)
      for (s in setdiff(ng_cp_stage_order(), "qc")) ctx <- ng_cp_pipeline[[s]](ctx)
      r <- ng_cp__assemble_result(ctx)
      job_dir <- file.path(output_root, job$id)
      env <- ng_run_result_envelope(r, warnings = seen, generated_at = generated_at,
                                    package_version = package_version)
      result_path <- file.path(job_dir, "result.json")
      ng_write_result_json(env, result_path)
      es <- r$effect_summary
      summary_fields <- list(
        cv_predictive_r2 = if (!is.null(es) && "cv_predictive_r2" %in% names(es))
          as.numeric(es$cv_predictive_r2)[[1L]] else NA_real_,
        mean_source = if (!is.null(es) && "mean_source" %in% names(es))
          as.character(es$mean_source)[[1L]] else NA_character_,
        # The gate verdict needs no new computation -- effect_gate already sits on the
        # effect_summary row this function reads for cv_predictive_r2 (R/39:1432) -- but it
        # is a new FIELD, and it is what tells a reader whether a plan rests on markers the
        # engine would otherwise have refused.
        effect_gate = if (!is.null(es) && "effect_gate" %in% names(es))
          as.character(es$effect_gate)[[1L]] else NA_character_,
        n_selected = if (is.data.frame(r$selected_crosses)) nrow(r$selected_crosses) else NA_integer_,
        # Recorded here so a RESUMED job can reproduce it (see the carried-forward branch of
        # ng_run_cross_prediction_batch()). Without this, output_files exists only on a job's
        # in-memory return value -- fine for a run that finishes in one process, but a
        # carried-forward entry is rebuilt from THIS FILE, and a manifest is how a consumer
        # finds the deliverable. Losing it here means resuming a job silently drops the
        # workbook path for every already-finished trait, even though the file is on disk.
        output_files = r$output_files)
      ng_job_trait_status_write(output_root, job$id, "done", summary_fields)
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
    }, error = function(e) {
      ng_job_trait_status_write(output_root, job$id, "error",
                                list(error_message = conditionMessage(e)))
      list(status = "error", traits = job$traits, error_message = conditionMessage(e))
    }),
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
utils::globalVariables(c(".root", ".shared", ".use_cpp", ".rng", ".config"))

ng_run_cross_prediction_batch <- function(config,
                                          jobs = NULL,
                                          output_root,
                                          batch_workers = NULL,
                                          memory_budget_bytes = NULL,
                                          job_dir = NULL,
                                          shared_dir = NULL,
                                          generated_at = NULL,
                                          package_version = NULL,
                                          resume = FALSE) {
  if (!is.list(config)) ng_stop("config must be a list of ng_run_cross_prediction() arguments")
  if (!length(output_root) || !nzchar(output_root)) ng_stop("output_root is required")
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  if (is.null(generated_at)) generated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  # When a job directory is given, this batch IS a job: it says so before the expensive work
  # starts, and says how it ended. Without one the function behaves exactly as in 0.36.0, so
  # a direct caller pays nothing for machinery it is not using.
  if (!is.null(job_dir)) {
    if (!file.exists(file.path(job_dir, "job.json"))) ng_job_create(job_dir, config = NULL)
    # `phase` is what a reader needs to tell "still doing shared setup" from "dead". R is
    # single-threaded, so nothing can beat the heartbeat while quality control, the duplicate
    # scan, LD pruning and the GRM run -- the most expensive phase on a real panel -- and a
    # reader with only a heartbeat mtime has to guess. See ng_job_stale_after_default.
    ng_job_mark(job_dir, "running", list(pid = Sys.getpid(), phase = "shared_setup"))
    ng_job_heartbeat(job_dir)
  }

  # Back-fill every missing formal from the runner's own defaults. ng_cp__build_ctx()
  # match.arg()s every enum, and an ABSENT enum arrives as its whole choices vector rather
  # than as its default -- the same trap tools/run_cross_prediction_json.R documents.
  fm <- formals(ng_run_cross_prediction)
  for (nm in setdiff(names(fm), names(config))) {
    config[nm] <- list(tryCatch(eval(fm[[nm]], envir = environment()), error = function(e) NULL))
  }

  # The key is computed on the RESOLVED config, after back-fill -- not on what the caller
  # happened to type. 14 names in ng_cp__batch_shared_keys carry non-NULL package defaults
  # (duplicate_threshold, duplicate_maf_min, ld_pruning, ld_window among them), and every one
  # of those materially changes the artefact. Keying on the pre-back-fill config would mean a
  # config that OMITS duplicate_threshold gets the same key forever -- so if a future release
  # changes that default, the stale key would still resolve to an artefact built under the
  # OLD default, and the batch would silently reuse it instead of rebuilding. Wrong numbers,
  # no error: the same shape of defect as the mean-source GEBV bug and the 0.36.0 RNG-kind
  # bug, both of which produced different results from what the caller asked for with
  # nothing to signal it. The key must describe the settings the artefact was ACTUALLY built
  # with, and those are the resolved ones.
  if (is.null(shared_dir)) shared_dir <- file.path(output_root, "_shared")
  key <- ng_shared_artifact_key(config)

  # Reuse before recomputing. The key covers every setting the batch spends plus the content
  # of each input file, so a hit means the SHARED WORK was built from exactly this data.
  #
  # A hit reuses that work and nothing else. ng_cp__batch_ctx() builds a fresh context from
  # THIS batch's config and lays only the artefact-owned fields over it, so n_crosses, the
  # scoring metric, the acceptability bar, the seed and every other per-batch setting come
  # from the caller who is running now -- see ng_shared_artifact_fields in R/55 for why that
  # distinction is the difference between a cache and a silent wrong answer.
  reuse_path <- file.path(shared_dir, key, "shared.rds")
  artifact <- ng_shared_artifact_read(reuse_path)
  if (!is.null(artifact)) {
    ctx <- ng_cp__batch_ctx(config, artifact)
  } else {
    ctx <- ng_cp__build_ctx(config)
    ctx <- ng_cp__stage_qc(ctx)
  }
  if (identical(ctx$qc$status, "blocker")) {
    ng_stop("QC blocker -- resolve before running a batch: ",
            paste(utils::head(ctx$qc$issues$message, 3L), collapse = "; "))
  }
  if (is.null(ctx$predict_prologue)) {
    ctx <- ng_ctx_put(ctx, predict_prologue = ng_cp__predict_prologue(ctx))
  }

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

  # On disk, so a killed batch resumes without redoing QC -- and CONTENT-ADDRESSED when a
  # store is given, so a second batch on the same data reuses it rather than recomputing
  # quality control, LD pruning and the GRM it has already paid for.
  # `key` was computed above for the reuse check -- do not recompute it. For an in-memory
  # genotype matrix that means one serialisation per batch rather than two.
  art_dir <- ng_shared_artifact_dir(shared_dir, key)
  shared_path <- file.path(art_dir, "shared.rds")
  # Only on a miss. Rewriting an artefact a concurrent batch may be reading buys nothing --
  # the content is a function of the key -- and a schema-rejected artefact is rebuilt here,
  # which is the one case where the existing file must be replaced.
  if (is.null(artifact)) ng_shared_artifact_write(ctx, shared_path)
  # Record the reference unconditionally, and be exact about what that buys.
  #
  # ng_job_prune() reads meta.json$referenced_by and honours any entry naming a job directory
  # that still exists, so for a batch running AS A JOB this is a genuine second source of
  # protection, independent of the job's own shared_ref file (which prune deliberately treats
  # as referencing nothing when it cannot be read). For a batch invoked DIRECTLY there is no
  # job directory, so what gets recorded here is basename(output_root) and prune cannot
  # honour it -- that caller is outside the job store, gets its own private shared_dir by
  # default, and is not something retention over jobs_dir can reason about. The entry is
  # still written because it is the only provenance record of who built the artefact.
  #
  # When this batch IS a job, also drop the pointer beside it so a killed-and-resumed job can
  # find the artefact it already claimed without recomputing the key.
  ref_id <- if (!is.null(job_dir)) basename(job_dir) else basename(output_root)
  if (!is.null(job_dir)) writeLines(key, file.path(job_dir, "shared_ref"))
  ng_shared_artifact_reference(shared_dir, key, ref_id)

  # The job now knows how many traits it is for. Without this, ng_job_status() falls back to
  # counting the status files that exist, so the DENOMINATOR grows as workers start and
  # "k of n done" is underivable -- the one number the overview exists to show. The parent is
  # the only place that knows length(jobs) before any worker has written anything.
  #
  # Marked together with the phase change, so the transition out of shared setup and the
  # trait count reach job.json in one atomic write.
  if (!is.null(job_dir)) {
    ng_job_mark(job_dir, "running",
                list(n_traits = length(jobs), phase = "scoring"))
  }

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

  # Resuming finishes a job rather than restarting it. A batch may have run for hours, so
  # re-running the traits that already succeeded is its own kind of loss. A trait counts as
  # finished only when its own status.json says "done" -- an "error" (or "running", from a
  # dead worker) must be re-run, because skipping it is exactly the loss resume exists to
  # avoid.
  carried <- list()
  run_jobs <- jobs
  if (isTRUE(resume)) {
    done_of <- function(id) {
      f <- file.path(output_root, id, "status.json")
      if (!file.exists(f)) return(FALSE)
      st <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
      identical(st$state, "done")
    }
    keep <- !vapply(names(jobs), done_of, logical(1))
    carried <- lapply(names(jobs)[!keep], function(id) {
      st <- jsonlite::fromJSON(file.path(output_root, id, "status.json"), simplifyVector = TRUE)
      list(id = id, traits = jobs[[id]]$traits, status = "ok",
           cv_predictive_r2 = as.numeric(st$cv_predictive_r2 %||% NA_real_),
           mean_source = as.character(st$mean_source %||% NA_character_),
           n_selected = as.integer(st$n_selected %||% NA_integer_),
           result_json = file.path(output_root, id, "result.json"),
           # FIX B: output_files is what a reader needs to FIND the deliverable, and a
           # carried entry (built here, not by a fresh run) can only reproduce it because
           # ng_cp__batch_run_one() now records it in status.json on the "done" write.
           output_files = st$output_files,
           warnings = character(0), elapsed_sec = NA_real_, resumed = TRUE)
    })
    # No names(carried) <- ...: an R list-name here was never load-bearing (every element
    # already carries its own $id, which is what the ordering below and every caller key on)
    # and naming it would re-introduce exactly the shape inconsistency ng_cp__batch_dispatch()
    # was just fixed to stop producing. One place decides whether `results` is named --
    # nowhere -- rather than a dispatch-shape fix here being undone a few lines later.
    run_jobs <- jobs[keep]
  }

  results <- if (!length(run_jobs)) list() else
    ng_cp__batch_dispatch(run_jobs, min(workers, length(run_jobs)), shared_path, output_root,
                          generated_at, package_version, use_cpp = ctx$use_cpp,
                          job_dir = job_dir, config = config)
  # The manifest describes the whole job. A resumed batch that listed only the traits it
  # re-ran would misrepresent what the breeder actually has on disk.
  results <- c(results, carried)
  # Order by each result's own $id, never by the list's R names(): ng_cp__batch_dispatch()
  # now always returns an unnamed list and `carried` is never named either (see above), but
  # ordering by $id keeps this correct even if that ever changes -- every result, carried or
  # freshly dispatched, always carries its own trait id.
  results <- results[order(match(
    vapply(results, function(r) as.character(r$id), character(1)), names(jobs)))]

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
      # I(): `traits` is conceptually plural and jsonlite's auto_unbox = TRUE renders a
      # length-1 character vector as a bare JSON string. A one-trait job's manifest entry
      # therefore read {"traits": "YIELD"} and a two-trait job's {"traits": ["A","B"]} -- the
      # same defect this file's `jobs` field was already fixed for, one level down. A
      # consumer would have to handle both shapes and would break on whichever it had not
      # tested, and the default batch is ONE trait per job, so the string form is the common
      # case and the array form the surprise. Applied here, at the serialisation boundary,
      # rather than on `results`: the R return value is indexed by callers and must keep its
      # ordinary character vectors.
      r$traits <- I(as.character(r$traits))
      r
    })
  )
  manifest_path <- ng_write_json_atomic(manifest, file.path(output_root, "manifest.json"))

  if (!is.null(job_dir)) {
    any_failed <- any(vapply(results, function(r) identical(r$status, "error"), logical(1)))
    # "finished" means the BATCH completed, not that every trait succeeded -- a failed trait
    # is recorded per trait, and marking the whole job failed would hide sixteen good plans
    # behind one bad one.
    # phase is left meaningful at rest too: "complete" rather than a stale "scoring", so a
    # reader never has to cross-check it against state to know it is out of date. A FAILED
    # job keeps whichever phase it died in -- that is the diagnostic.
    ng_job_mark(job_dir, "finished", list(n_ok = sum(!vapply(
      results, function(r) identical(r$status, "error"), logical(1))),
      any_failed = any_failed, phase = "complete"))
  }

  # `results` is unnamed by construction (see ng_cp__batch_dispatch() and the carried-forward
  # comment above), so there is nothing to normalise here. Callers index this list by each
  # element's own $id field (see by_id() in the batch test suite), never by R list-name.
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
                                  generated_at, package_version, use_cpp = FALSE,
                                  job_dir = NULL, config = NULL) {
  # unname() every return path here, not just some of them: `jobs` carries R list-names (the
  # trait ids) on every path, but only lapply() over a NAMED list propagates them into the
  # result. Leaving even one path named made the manifest's `jobs` field serialise as a JSON
  # OBJECT on that path and a JSON ARRAY on the others -- a shape a consumer would have to
  # handle twice and would break on whichever it had not tested. Every result already carries
  # its own `id` field, so the R list-name was never load-bearing; only jsonlite's rendering
  # of it was.
  if (workers <= 1L || !requireNamespace("mirai", quietly = TRUE)) {
    if (workers > 1L) {
      warning("mirai is not available; running the batch one job at a time.", call. = FALSE)
    }
    return(unname(lapply(jobs, function(j) {
      # Serial jobs beat the heart between traits. The loop over unresolved() below never
      # runs on this path, and a job that goes quiet is indistinguishable from a dead one.
      if (!is.null(job_dir)) ng_job_heartbeat(job_dir)
      ng_cp__batch_run_one(j, shared_path = shared_path, output_root = output_root,
                           generated_at = generated_at, package_version = package_version,
                           config = config)
    })))
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
      return(unname(lapply(jobs, ng_cp__batch_run_one, shared_path = shared_path,
                    output_root = output_root, generated_at = generated_at,
                    package_version = package_version, config = config)))
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
    # Read the shared ARTEFACT once per daemon rather than once per job: it carries the
    # genotype matrix and the GRM, and a daemon may serve several jobs in turn. The batch's
    # own config rides along beside it -- the two are assembled into a context per job by
    # ng_cp__batch_ctx(), because the artefact must never supply a caller's settings.
    #
    # everywhere() evaluates in the daemon's GLOBAL environment. On the dev path ng_load()
    # has just sourced every internal (ng_shared_artifact_read included) into that same
    # global environment, so a plain call resolves. On the installed path library() exposes
    # only EXPORTED functions -- ng_shared_artifact_read is internal, so it must be fetched
    # from the namespace explicitly instead. Getting this wrong fails in exactly one
    # direction -- the dev tree works, production does not -- which is how it reached a
    # release: every batch test lived under tests/*.R, which runs against the source tree,
    # and none lived under tests/testthat/, which runs against the installed package.
    .read_artifact <- if (nzchar(.root)) get("ng_shared_artifact_read", envir = globalenv())
                       else getFromNamespace("ng_shared_artifact_read", "nextgenCrossDesign")
    assign(".ngcd_batch_artifact", .read_artifact(.shared), envir = globalenv())
    assign(".ngcd_batch_config", .config, envir = globalenv())
  }, .root = dev_root, .use_cpp = isTRUE(use_cpp), .shared = shared_path, .rng = parent_rng,
     .config = config)

  m <- mirai::mirai_map(
    jobs, ng_cp__batch_run_one,
    .args = list(shared_path = shared_path, output_root = output_root,
                 generated_at = generated_at, package_version = package_version))
  # NOTE: `config` is NOT in .args -- mirai serialises .args per call, and a config carrying
  # an in-memory genotype matrix would then cross the socket once per trait. It reaches the
  # daemons once, through the everywhere() block above.
  # Collect by polling rather than with a bare m[], for one reason: m[] blocks until every
  # job resolves, and a job that is working hard for two hours would be indistinguishable
  # from a job whose process died an hour ago. Touching an mtime is all this costs -- no file
  # content is rewritten, and progress itself still comes from the workers' own status files.
  if (!is.null(job_dir)) {
    while (any(mirai::unresolved(m))) {
      ng_job_heartbeat(job_dir)
      Sys.sleep(2)
    }
  }
  out <- m[]
  # A daemon that died leaves a miraiError in place of the job's result. Convert it to the
  # same failure shape a caught error produces, so one dead worker cannot make the manifest
  # a different shape from a failed run.
  #
  # unname(): explicit here even though lapply(seq_along(out), ...) already drops names, so
  # this path's shape does not depend on that incidental fact holding across a future edit --
  # see the comment on the serial return above.
  unname(lapply(seq_along(out), function(i) {
    r <- out[[i]]
    if (inherits(r, "miraiError") || inherits(r, "errorValue")) {
      return(list(id = jobs[[i]]$id, traits = jobs[[i]]$traits, status = "error",
                  error_message = paste("worker failed:", conditionMessage(r)),
                  warnings = character(0), elapsed_sec = NA_real_))
    }
    r
  }))
}
