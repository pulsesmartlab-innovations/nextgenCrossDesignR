# Optimizer bake-off harness: run every mate-allocation optimizer on the SAME fixed,
# reproducible score + kinship inputs and report, per method, the achieved objective
# (the common yardstick), the OCS summary metrics, wall-clock runtime, and a
# gain-vs-coancestry Pareto sweep. This is the measurement instrument used to decide
# which allocator is best and to tune knobs BEFORE the full simulation study.
#
# Package-only fixture (no AlphaSimR needed): a synthetic genotype matrix is scored
# with the package's own ng_fit_ridge_effects + ng_score_crosses, and parent_kinship via
# ng_parent_kinship -- so the candidate table is exactly what the real pipeline feeds
# the optimizer. External comparators (AlphaMate binary, SimpleMating) are included
# only when available and never halt the run.
#
# Run:  Rscript tools/run_optimizer_benchmark.R
# Env:  NG_OPT_* overrides (see optimizer_config_from_env).

find_project_root <- function(start = getwd()) {
  candidates <- unique(c(
    Sys.getenv("NG_REPO_ROOT", unset = NA_character_),
    start, file.path(start, "nextgen_cross_design"),
    file.path(start, ".."), file.path(start, "..", "nextgen_cross_design")
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "load.R")) && file.exists(file.path(cand, "DESCRIPTION"))) {
      return(normalizePath(cand, winslash = "/", mustWork = FALSE))
    }
  }
  stop("Could not locate nextgenCrossDesign package root", call. = FALSE)
}

env_chr <- function(name, default) { v <- Sys.getenv(name, unset = NA_character_); if (is.na(v) || !nzchar(v)) default else v }
env_int <- function(name, default) { v <- suppressWarnings(as.integer(Sys.getenv(name, unset = NA_character_))); if (is.na(v)) as.integer(default) else v }
env_num <- function(name, default) { v <- suppressWarnings(as.numeric(Sys.getenv(name, unset = NA_character_))); if (is.na(v)) as.numeric(default) else v }
env_csv <- function(name, default) { v <- Sys.getenv(name, unset = NA_character_); if (is.na(v) || !nzchar(v)) default else trimws(strsplit(v, ",")[[1]]) }
env_num_vec <- function(name, default) { v <- Sys.getenv(name, unset = NA_character_); if (is.na(v) || !nzchar(v)) default else as.numeric(trimws(strsplit(v, ",")[[1]])) }

optimizer_config_from_env <- function(root = find_project_root()) {
  list(
    root = root,
    seed = env_int("NG_OPT_SEED", 20260701L),
    n_parents = env_int("NG_OPT_N_PARENTS", 60L),
    n_markers = env_int("NG_OPT_N_MARKERS", 600L),
    n_chr = env_int("NG_OPT_N_CHR", 5L),
    training_n = env_int("NG_OPT_TRAINING_N", 200L),
    n_crosses = env_int("NG_OPT_N_CROSSES", 30L),
    max_crosses_per_parent = env_int("NG_OPT_MAX_PER_PARENT", 4L),
    lambda_group = env_num("NG_OPT_LAMBDA_GROUP", 0.5),
    lambda_parent_use = env_num("NG_OPT_LAMBDA_PARENT_USE", 1.0),
    lambda_parent_use_mode = env_chr("NG_OPT_LAMBDA_PARENT_USE_MODE", "adaptive"),
    gain_col = env_chr("NG_OPT_GAIN_COL", "usefulness_pmv_gebv"),
    methods = env_csv("NG_OPT_METHODS", c("mip_contribution", "greedy_local", "repair_local", "evolution", "alphamate_style")),
    lambda_grid = env_num_vec("NG_OPT_LAMBDA_GRID", 10^seq(-2, 2.5, length.out = 10)),
    local_iter = env_int("NG_OPT_LOCAL_ITER", 2000L),
    ocs_iter = env_int("NG_OPT_OCS_ITER", 5L),
    use_cpp = as.logical(env_int("NG_OPT_USE_CPP", 1L)),
    reps = env_int("NG_OPT_REPS", 1L),
    output_dir = env_chr("NG_OPT_OUTPUT_DIR", file.path(root, "results")),
    output_prefix = env_chr("NG_OPT_OUTPUT_PREFIX", "optimizer_bakeoff")
  )
}

# Build a fixed, realistic candidate score table + parent kinship from a synthetic
# genotype matrix using the package's own scoring path. Cached to inst/fixtures.
build_optimizer_fixture <- function(cfg) {
  fixture_dir <- file.path(cfg$root, "inst", "fixtures")
  key <- sprintf("optimizer_bakeoff_p%d_m%d_t%d_s%d", cfg$n_parents, cfg$n_markers, cfg$training_n, cfg$seed)
  path <- file.path(fixture_dir, paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))

  n <- cfg$n_parents; m <- cfg$n_markers; tn <- cfg$training_n
  ids <- paste0("P", seq_len(n)); markers <- paste0("M", seq_len(m))
  fx <- ng_with_rng_seed(cfg$seed, {
    all_ids <- paste0("L", seq_len(n + tn))
    geno_all <- matrix(2L * rbinom((n + tn) * m, 1, 0.45), nrow = n + tn, ncol = m,
                       dimnames = list(all_ids, markers))
    true_beta <- rnorm(m, sd = 0.08)
    y_all <- as.numeric(geno_all %*% true_beta + rnorm(n + tn, sd = 1.0))
    names(y_all) <- all_ids
    list(geno_all = geno_all, y_all = y_all, all_ids = all_ids)
  })
  chr <- rep(seq_len(cfg$n_chr), length.out = m)
  pos <- unlist(lapply(split(seq_len(m), chr), function(idx) seq(0, 150, length.out = length(idx))), use.names = FALSE)
  marker_map <- data.frame(marker = markers, chr = chr, pos_cm = pos, stringsAsFactors = FALSE)

  train_ids <- fx$all_ids[seq_len(tn)]
  parent_ids <- fx$all_ids[tn + seq_len(n)]
  fit <- ng_fit_ridge_effects(fx$geno_all[train_ids, , drop = FALSE], fx$y_all[train_ids],
                              ids = train_ids, kfold = 0L, seed = cfg$seed)
  parent_geno <- fx$geno_all[parent_ids, , drop = FALSE]
  rownames(parent_geno) <- ids
  scores <- ng_score_crosses(parent_geno, fit, marker_map = marker_map, ids = ids,
                             adjusted_pheno = setNames(fx$y_all[parent_ids], ids),
                             selection_prop = 0.10, use_cpp = cfg$use_cpp)
  parent_kinship <- ng_parent_kinship(parent_geno)
  out <- list(scores = scores, parent_kinship = parent_kinship, ids = ids, gain_col = cfg$gain_col,
              fixture_key = key)
  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  tryCatch(saveRDS(out, path), error = function(e) NULL)
  out
}

# Registry: name -> function(scores, n_crosses, parent_kinship, lambda_group, lambda_parent_use, cfg) returning a plan.
optimizer_registry <- function() {
  ocs_call <- function(method) function(scores, n_crosses, parent_kinship, lg, lpu, cfg) {
    ng_optimize_mating_plan(scores = scores, n_crosses = n_crosses, gain_col = cfg$gain_col,
                            parent_kinship = parent_kinship, max_crosses_per_parent = cfg$max_crosses_per_parent,
                            lambda_group = lg, lambda_parent_use = lpu,
                            lambda_parent_use_mode = cfg$lambda_parent_use_mode,
                            method = method, local_iter = cfg$local_iter, ocs_iter = cfg$ocs_iter)
  }
  list(
    mip_contribution = ocs_call("mip_contribution"),
    greedy_local     = ocs_call("greedy_local"),
    repair_local     = ocs_call("repair_local"),
    evolution        = ocs_call("evolution"),
    alphamate_style  = function(scores, n_crosses, parent_kinship, lg, lpu, cfg) {
      ng_alphamate_style_select(scores = scores, criterion_col = cfg$gain_col, n_crosses = n_crosses,
                                parent_kinship = parent_kinship, mode = "ModeOptTarget1", target_degree = 45,
                                max_contributions = cfg$max_crosses_per_parent,
                                method = "greedy_local", local_iter = cfg$local_iter)
    }
  )
}

# Map a returned plan (subset of scores rows) back to scores row indices by pair key.
plan_row_indices <- function(plan, scores) {
  key_all <- ng_pair_key(as.character(scores$parent1), as.character(scores$parent2))
  key_sel <- ng_pair_key(as.character(plan$parent1), as.character(plan$parent2))
  idx <- match(key_sel, key_all)
  idx[!is.na(idx)]
}

# Achieved objective on the SAME objective every method is scored against.
achieved_objective <- function(plan, scores, parent_kinship, lg, lpu, gain_col, lambda_mating = 0) {
  sc <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!("pair_kinship" %in% names(sc))) sc$pair_kinship <- 0
  sc$.linear_gain <- sc[[gain_col]] - lambda_mating * sc$pair_kinship
  idx <- plan_row_indices(plan, sc)
  if (!length(idx)) return(NA_real_)
  ng_plan_objective_contribution(sc, idx, parent_kinship, lambda_group = lg, lambda_parent_use = lpu)
}

run_one_method <- function(name, fn, fixture, cfg, lg, lpu) {
  err <- NULL
  t <- system.time(plan <- tryCatch(
    ng_with_rng_seed(cfg$seed, fn(fixture$scores, cfg$n_crosses, fixture$parent_kinship, lg, lpu, cfg)),
    error = function(e) { err <<- conditionMessage(e); NULL }
  ))["elapsed"]
  if (is.null(plan)) {
    return(data.frame(method = name, status = "error",
                      reason = if (!is.null(err)) err else "returned NULL",
                      achieved_objective = NA_real_, mean_gain = NA_real_, group_coancestry = NA_real_,
                      parent_use_sq = NA_real_, unique_parents = NA_integer_, max_parent_use = NA_integer_,
                      n_selected = 0L, elapsed_sec = as.numeric(t), stringsAsFactors = FALSE))
  }
  s <- attr(plan, "summary"); if (is.null(s)) s <- list()
  data.frame(
    method = name, status = "ok", reason = "",
    achieved_objective = achieved_objective(plan, fixture$scores, fixture$parent_kinship, lg, lpu, cfg$gain_col),
    mean_gain = if (!is.null(s$mean_gain)) s$mean_gain else mean(plan[[cfg$gain_col]], na.rm = TRUE),
    group_coancestry = if (!is.null(s$group_coancestry)) s$group_coancestry else NA_real_,
    parent_use_sq = if (!is.null(s$parent_use_sq)) s$parent_use_sq else NA_real_,
    unique_parents = if (!is.null(s$unique_parents)) s$unique_parents else length(unique(c(plan$parent1, plan$parent2))),
    max_parent_use = if (!is.null(s$max_parent_use)) s$max_parent_use else NA_integer_,
    n_selected = nrow(plan), elapsed_sec = as.numeric(t), stringsAsFactors = FALSE
  )
}

run_optimizer_benchmark <- function(cfg = optimizer_config_from_env()) {
  source(file.path(cfg$root, "R", "load.R"))
  ng_load(cfg$root, use_cpp = cfg$use_cpp, verbose = FALSE)

  fixture <- build_optimizer_fixture(cfg)
  reg <- optimizer_registry()
  methods <- intersect(cfg$methods, names(reg))
  message(sprintf("Optimizer bake-off: %d parents, %d candidate crosses, %d markers; methods: %s",
                  cfg$n_parents, nrow(fixture$scores), cfg$n_markers, paste(methods, collapse = ", ")))

  lg <- cfg$lambda_group; lpu <- cfg$lambda_parent_use
  summary_rows <- lapply(methods, function(nm) run_one_method(nm, reg[[nm]], fixture, cfg, lg, lpu))
  summary_df <- do.call(rbind, summary_rows)
  ok <- summary_df$status == "ok" & is.finite(summary_df$achieved_objective)
  summary_df$objective_rank <- NA_integer_
  if (any(ok)) summary_df$objective_rank[ok] <- rank(-summary_df$achieved_objective[ok], ties.method = "min")

  # Pareto sweep: gain vs group_coancestry across a lambda_group grid, per method.
  frontier_rows <- list()
  for (nm in methods) {
    for (g in cfg$lambda_grid) {
      r <- tryCatch(run_one_method(nm, reg[[nm]], fixture, cfg, g, lpu), error = function(e) NULL)
      if (is.null(r) || r$status != "ok") next
      frontier_rows[[length(frontier_rows) + 1L]] <- data.frame(
        method = nm, lambda_group = g, mean_gain = r$mean_gain,
        group_coancestry = r$group_coancestry, unique_parents = r$unique_parents,
        max_parent_use = r$max_parent_use, elapsed_sec = r$elapsed_sec, stringsAsFactors = FALSE)
    }
  }
  frontier_df <- if (length(frontier_rows)) do.call(rbind, frontier_rows) else data.frame()

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(cfg$output_dir, cfg$output_prefix)
  write.csv(summary_df, paste0(prefix, "_summary.csv"), row.names = FALSE)
  if (nrow(frontier_df)) write.csv(frontier_df, paste0(prefix, "_frontier.csv"), row.names = FALSE)
  meta <- list(fixture_key = fixture$fixture_key, seed = cfg$seed, n_parents = cfg$n_parents,
               n_candidate_crosses = nrow(fixture$scores), n_crosses = cfg$n_crosses,
               lambda_group = lg, lambda_parent_use = lpu, gain_col = cfg$gain_col,
               package_version = tryCatch(as.character(read.dcf(file.path(cfg$root, "DESCRIPTION"))[1, "Version"]),
                                          error = function(e) NA_character_),
               methods = methods)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE), paste0(prefix, "_meta.json"))
  }
  writeLines(capture.output(sessionInfo()), paste0(prefix, "_session_info.txt"))

  message("Wrote optimizer bake-off outputs to: ", normalizePath(cfg$output_dir, winslash = "/", mustWork = FALSE))
  print(summary_df[, c("method", "status", "achieved_objective", "mean_gain", "group_coancestry",
                       "unique_parents", "max_parent_use", "elapsed_sec", "objective_rank")], row.names = FALSE)
  invisible(list(summary = summary_df, frontier = frontier_df, meta = meta, fixture = fixture))
}

is_this_script <- function() any(grepl("run_optimizer_benchmark\\.R$", commandArgs(FALSE)))
if (is_this_script()) invisible(run_optimizer_benchmark())
