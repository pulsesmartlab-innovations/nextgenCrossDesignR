# Two-axis benchmark: cross-prediction metric versus mating optimizer.
#
# Metrics:
# - var_simple
# - SimpleMating getMPV, getTGV, getUsefA, getUsefAD nonphased/phased
#
# Optimizers/selectors:
# - topn: current project behavior, rank by metric score and keep top N crosses.
# - simplemating_select: SimpleMating::selectCrosses() on the same scored table.
#
# This isolates:
# 1. Same optimizer, different metrics.
# 2. Same metric, different optimizers.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
Sys.setenv(CROSSPRED_SKIP_CPP = Sys.getenv("OX_SKIP_CPP", unset = "1"))
suppressPackageStartupMessages(library(AlphaSimR))
source(file.path("R", "load.R"))

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.integer(value)
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.numeric(value)
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

env_count <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value <- trimws(tolower(value))
  if (value %in% c("all", "inf", "infinite")) return(Inf)
  as.integer(value)
}

bind_rows_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    for (m in missing) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  do.call(rbind, x)
}

merge_pop_list <- function(pop_list) {
  pop_list <- Filter(Negate(is.null), pop_list)
  if (!length(pop_list)) stop("No populations to merge.", call. = FALSE)
  Reduce(c, pop_list)
}

make_training_population <- function(parent_pop, target_n, sim_param) {
  n_parent <- nInd(parent_pop)
  if (target_n <= n_parent) return(parent_pop)
  extra_n <- target_n - n_parent
  p1 <- sample.int(n_parent, extra_n, replace = TRUE)
  p2 <- sample.int(n_parent, extra_n, replace = TRUE)
  same <- p1 == p2
  while (any(same)) {
    p2[same] <- sample.int(n_parent, sum(same), replace = TRUE)
    same <- p1 == p2
  }
  f1 <- makeCross(parent_pop, cbind(p1, p2), nProgeny = 1, simParam = sim_param)
  aux <- makeDH(f1, nDH = 1, keepParents = FALSE, simParam = sim_param)
  c(parent_pop, aux)
}

make_dh_base_parents <- function(founder_pop, n_parents, sim_param) {
  n_dh_per_founder <- ceiling(n_parents / nInd(founder_pop))
  dh <- makeDH(founder_pop, nDH = n_dh_per_founder, keepParents = FALSE, simParam = sim_param)
  dh[seq_len(n_parents)]
}

as_snp_haplo <- function(pop, sim_param, marker_names) {
  haplo <- pullSnpHaplo(pop, simParam = sim_param)
  colnames(haplo) <- marker_names
  rownames(haplo) <- as.vector(rbind(paste0(pop@id, "_HapA"), paste0(pop@id, "_HapB")))
  haplo
}

population_metrics <- function(pop, base_best, cycle, rep_id, metric, optimizer, branch) {
  g <- as.numeric(gv(pop)[, 1])
  data.frame(
    rep = rep_id,
    metric = metric,
    optimizer = optimizer,
    branch = branch,
    cycle = cycle,
    n_ind = length(g),
    mean_gv = mean(g),
    max_gv = max(g),
    top10_gv = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
    var_gv = stats::var(g),
    transgressive_rate = mean(g > base_best),
    stringsAsFactors = FALSE
  )
}

make_selected_progeny <- function(parent_pop, selected, progeny_per_cross, sim_param) {
  parent_ids <- parent_pop@id
  progeny <- vector("list", nrow(selected))
  family_metrics <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- match(selected$parent1[i], parent_ids)
    p2 <- match(selected$parent2[i], parent_ids)
    if (is.na(p1) || is.na(p2)) stop("Selected parent ID missing from parent population.", call. = FALSE)
    f1 <- makeCross(parent_pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
    dh <- makeDH(f1, nDH = progeny_per_cross, keepParents = FALSE, simParam = sim_param)
    g <- as.numeric(gv(dh)[, 1])
    family_metrics[[i]] <- data.frame(
      cross_rank = i,
      parent1 = as.character(selected$parent1[i]),
      parent2 = as.character(selected$parent2[i]),
      pred_mean = selected$cross_mean[i],
      pred_var = selected$cross_var[i],
      pred_usefulness = selected$cross_usefulness[i],
      realized_mean_gv = mean(g),
      realized_max_gv = max(g),
      realized_top10_gv = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
      realized_var_gv = stats::var(g),
      stringsAsFactors = FALSE
    )
    progeny[[i]] <- dh
  }
  list(pop = merge_pop_list(progeny), family_metrics = bind_rows_fill(family_metrics))
}

quiet_simplemating <- function(expr) {
  value <- NULL
  utils::capture.output({
    value <- force(expr)
  })
  value
}

is_simple_metric <- function(metric) {
  metric %in% c("simple_mpv", "simple_tgv", "simple_usefa",
                "simple_usefad_nonphased", "simple_usefad_phased")
}

needs_effects <- function(metric) {
  metric %in% c("simple_tgv", "simple_usefa",
                "simple_usefad_nonphased", "simple_usefad_phased")
}

as_simplemating_map <- function(marker_map) {
  data.frame(
    Chromosome = marker_map$chr,
    Position = marker_map$pos * 100,
    Marker = marker_map$marker,
    stringsAsFactors = FALSE
  )
}

simplemating_k <- function(geno) {
  K <- VanRadenKin(geno)
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  K
}

add_pair_relationship <- function(cross_df, geno) {
  K <- simplemating_k(geno)
  i <- match(cross_df$parent1, rownames(K))
  j <- match(cross_df$parent2, rownames(K))
  cross_df$pair_relationship <- as.numeric(K[cbind(i, j)])
  cross_df
}

estimate_effects <- function(training_pop, sim_param, marker_names, y, cfg) {
  training_geno <- pullSnpGeno(training_pop, simParam = sim_param)
  colnames(training_geno) <- marker_names
  rownames(training_geno) <- training_pop@id
  effects <- estimate_marker_effect_draws(
    geno_mat = training_geno,
    y = y,
    method = cfg$effect_method,
    n_draws = cfg$effect_draws,
    lambda = cfg$ridge_lambda,
    prior_genetic_var = stats::var(y) * cfg$prior_genetic_var_scale
  )
  beta <- attr(effects, "posterior_mean")
  if (is.null(beta)) beta <- colMeans(as.matrix(effects), na.rm = TRUE)
  beta <- as.numeric(beta)
  names(beta) <- marker_names
  list(effects = effects, beta = beta)
}

simple_mean_df <- function(sm_df, parent_values, ids, original_from_simple) {
  out <- data.frame(
    parent1 = unname(original_from_simple[as.character(sm_df$Parent1)]),
    parent2 = unname(original_from_simple[as.character(sm_df$Parent2)]),
    cross_mean_package = as.numeric(sm_df$Y),
    cross_var = 0,
    cross_sigma = 0,
    cross_usefulness = as.numeric(sm_df$Y),
    pair_relationship = if ("K" %in% names(sm_df)) as.numeric(sm_df$K) else NA_real_,
    stringsAsFactors = FALSE
  )
  original_ids <- unname(original_from_simple[ids])
  original_values <- stats::setNames(as.numeric(parent_values[ids]), original_ids)
  out$cross_mean <- parent_cross_mean(original_values, out$parent1, out$parent2, original_ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "simplemating_package_score"
  out
}

simple_usefulness_df <- function(sm_detail, parent_values, ids, original_from_simple, var_cols) {
  cross_var <- Reduce(`+`, lapply(var_cols, function(v) {
    if (v %in% names(sm_detail)) as.numeric(sm_detail[[v]]) else rep(0, nrow(sm_detail))
  }))
  out <- data.frame(
    parent1 = unname(original_from_simple[as.character(sm_detail$Parent1)]),
    parent2 = unname(original_from_simple[as.character(sm_detail$Parent2)]),
    cross_mean_package = as.numeric(sm_detail$Mean),
    cross_var = pmax(0, cross_var),
    cross_usefulness = as.numeric(sm_detail$Usefulness),
    stringsAsFactors = FALSE
  )
  original_ids <- unname(original_from_simple[ids])
  original_values <- stats::setNames(as.numeric(parent_values[ids]), original_ids)
  out$cross_sigma <- sqrt(pmax(0, out$cross_var))
  out$cross_mean <- parent_cross_mean(original_values, out$parent1, out$parent2, original_ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "simplemating_package_score"
  out
}

predict_simple_metric <- function(metric, parent_values, geno, haplo, marker_map, beta, pair_df, tau) {
  if (!requireNamespace("SimpleMating", quietly = TRUE)) {
    stop("SimpleMating is required for ", metric, call. = FALSE)
  }
  original_ids <- rownames(geno)
  simple_ids <- paste0("SM", seq_along(original_ids))
  names(simple_ids) <- original_ids
  original_from_simple <- stats::setNames(original_ids, simple_ids)

  geno_simple <- as.matrix(geno)
  rownames(geno_simple) <- simple_ids
  parent_values_simple <- stats::setNames(as.numeric(parent_values[original_ids]), simple_ids)
  pair_simple <- data.frame(
    Parent1 = unname(simple_ids[as.character(pair_df$parent1)]),
    Parent2 = unname(simple_ids[as.character(pair_df$parent2)]),
    stringsAsFactors = FALSE
  )
  K <- simplemating_k(geno_simple)
  sm_map <- as_simplemating_map(marker_map)
  beta <- as.numeric(beta[colnames(geno_simple)])
  names(beta) <- colnames(geno_simple)
  dom <- rep(0, length(beta))
  names(dom) <- names(beta)

  if (metric == "simple_mpv") {
    crit <- data.frame(Id = simple_ids, Criterion = as.numeric(parent_values[original_ids]))
    sm <- quiet_simplemating(SimpleMating::getMPV(pair_simple, crit, K))
    return(simple_mean_df(sm, parent_values_simple, simple_ids, original_from_simple))
  }

  if (metric == "simple_tgv") {
    sm <- quiet_simplemating(SimpleMating::getTGV(pair_simple, geno_simple, beta, dom, K, ploidy = 2))
    return(simple_mean_df(sm, parent_values_simple, simple_ids, original_from_simple))
  }

  if (metric == "simple_usefa") {
    sm <- quiet_simplemating(SimpleMating::getUsefA(
      MatePlan = pair_simple,
      Markers = geno_simple,
      addEff = beta,
      K = K,
      Map.In = sm_map,
      propSel = tau,
      Type = "DH",
      Generation = 1,
      n_threads = 1,
      display_progress = FALSE
    ))
    return(simple_usefulness_df(sm[[1]], parent_values_simple, simple_ids, original_from_simple, "Variance"))
  }

  if (metric %in% c("simple_usefad_nonphased", "simple_usefad_phased")) {
    markers <- geno_simple
    method <- "NonPhased"
    if (metric == "simple_usefad_phased") {
      markers <- as.matrix(haplo)
      rownames(markers) <- as.vector(rbind(paste0(simple_ids, "_HapA"), paste0(simple_ids, "_HapB")))
      method <- "Phased"
    }
    sm <- quiet_simplemating(SimpleMating::getUsefAD(
      MatePlan = pair_simple,
      Markers = markers,
      addEff = beta,
      domEff = dom,
      K = K,
      Map.In = sm_map,
      propSel = tau,
      Method = method,
      ploidy = 2,
      n_threads = 1,
      display_progress = FALSE
    ))
    return(simple_usefulness_df(sm[[1]], parent_values_simple, simple_ids, original_from_simple,
                                c("Var_A", "Var_D")))
  }

  stop("Unsupported SimpleMating metric: ", metric, call. = FALSE)
}

score_crosses <- function(metric, parent_pop, parent_values, geno, marker_map, haplo,
                          pair_df, effects, beta, cfg) {
  if (metric == "var_simple") {
    out <- build_cross_data(
      pheno = parent_values,
      geno_mat = geno,
      candidate_ids = rownames(geno),
      variance_method = "var_simple",
      include_self = FALSE,
      use_parallel = FALSE,
      tau = cfg$tau,
      h2 = cfg$h2,
      trait_direction = "increase"
    )
    return(add_pair_relationship(out, geno))
  }

  if (metric %in% c("simple_usefa", "simple_usefad_nonphased", "simple_usefad_phased") &&
      is.finite(cfg$simple_pair_limit) && cfg$simple_pair_limit < nrow(pair_df)) {
    var_screen <- score_crosses("var_simple", parent_pop, parent_values, geno, marker_map, NULL,
                                pair_df, NULL, NULL, cfg)
    tgv_screen <- predict_simple_metric("simple_tgv", parent_values, geno, NULL, marker_map, beta,
                                        pair_df, cfg$tau)
    n_each <- max(cfg$top_crosses, as.integer(cfg$simple_pair_limit))
    v <- var_screen[order(var_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]
    t <- tgv_screen[order(tgv_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]
    pair_df <- unique(rbind(
      utils::head(v[, c("parent1", "parent2"), drop = FALSE], n_each),
      utils::head(t[, c("parent1", "parent2"), drop = FALSE], n_each)
    ))
  }

  out <- predict_simple_metric(metric, parent_values, geno, haplo, marker_map, beta, pair_df, cfg$tau)
  out <- add_pair_relationship(out, geno)
  out$screened_n_pairs <- nrow(pair_df)
  out$screened_from_all_pairs <- choose(nrow(geno), 2)
  out
}

select_topn <- function(cross_df, top_crosses) {
  cross_df <- cross_df[order(cross_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  out <- utils::head(cross_df, top_crosses)
  out$optimizer_status <- "topn"
  out
}

select_simplemating <- function(cross_df, cfg) {
  if (!requireNamespace("SimpleMating", quietly = TRUE)) {
    stop("SimpleMating is required for simplemating_select.", call. = FALSE)
  }
  sel_input <- data.frame(
    Parent1 = as.character(cross_df$parent1),
    Parent2 = as.character(cross_df$parent2),
    Y = as.numeric(cross_df$cross_usefulness),
    K = as.numeric(cross_df$pair_relationship),
    stringsAsFactors = FALSE
  )
  sel_input <- sel_input[is.finite(sel_input$Y) & is.finite(sel_input$K), , drop = FALSE]
  selected_plan <- tryCatch({
    out <- quiet_simplemating(SimpleMating::selectCrosses(
      data = sel_input,
      n.cross = cfg$top_crosses,
      max.cross = cfg$sm_max_cross,
      min.cross = cfg$sm_min_cross,
      max.cross.to.search = cfg$sm_max_search,
      culling.pairwise.k = cfg$sm_culling_pairwise_k
    ))
    out$plan
  }, error = function(e) {
    warning("SimpleMating::selectCrosses failed; falling back to top-N. Error: ", conditionMessage(e))
    NULL
  })
  if (is.null(selected_plan) || !nrow(selected_plan)) return(select_topn(cross_df, cfg$top_crosses))

  key <- paste(cross_df$parent1, cross_df$parent2, sep = "\r")
  selected_key <- paste(selected_plan$Parent1, selected_plan$Parent2, sep = "\r")
  idx <- match(selected_key, key)
  out <- cross_df[idx[!is.na(idx)], , drop = FALSE]
  out$optimizer_status <- "simplemating_select"
  out$optimizer_score <- selected_plan$Y[!is.na(idx)]
  out
}

parent_use_summary <- function(selected) {
  use <- table(c(as.character(selected$parent1), as.character(selected$parent2)))
  data.frame(
    n_selected_crosses = nrow(selected),
    n_unique_parents = length(use),
    max_parent_use = max(use),
    mean_parent_use = mean(use),
    stringsAsFactors = FALSE
  )
}

cfg <- list(
  seed = env_int("OX_SEED", 20260426L),
  reps = env_int("OX_REPS", 3L),
  cycles = env_int("OX_CYCLES", 1L),
  n_founders = env_int("OX_N_FOUNDERS", 120L),
  n_parents = env_int("OX_N_PARENTS", 80L),
  n_chr = env_int("OX_N_CHR", 5L),
  seg_sites = env_int("OX_SEG_SITES", 1200L),
  qtl_per_chr = env_int("OX_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("OX_SNP_PER_CHR", 1000L),
  top_crosses = env_int("OX_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("OX_PROGENY_PER_CROSS", 40L),
  effect_training_n = env_int("OX_EFFECT_TRAINING_N", env_int("OX_N_PARENTS", 80L)),
  effect_method = env_chr("OX_EFFECT_METHOD", "ridge_posterior"),
  effect_draws = env_int("OX_EFFECT_DRAWS", 20L),
  ridge_lambda = env_num("OX_RIDGE_LAMBDA", 1),
  prior_genetic_var_scale = env_num("OX_PRIOR_GENETIC_VAR_SCALE", 1),
  tau = env_num("OX_TAU", 0.10),
  h2 = env_num("OX_H2", 1.0),
  metrics = strsplit(env_chr(
    "OX_METRICS",
    "var_simple,simple_mpv,simple_tgv,simple_usefa,simple_usefad_nonphased,simple_usefad_phased"
  ), ",", fixed = TRUE)[[1]],
  optimizers = strsplit(env_chr("OX_OPTIMIZERS", "topn,simplemating_select"), ",", fixed = TRUE)[[1]],
  simple_pair_limit = env_count("OX_SIMPLE_PAIR_LIMIT", 500L),
  sm_max_cross = env_int("OX_SM_MAX_CROSS", 4L),
  sm_min_cross = env_int("OX_SM_MIN_CROSS", 1L),
  sm_max_search = env_count("OX_SM_MAX_SEARCH", 100000L),
  sm_culling_pairwise_k = env_num("OX_SM_CULLING_PAIRWISE_K", Inf),
  output_prefix = env_chr("OX_OUTPUT_PREFIX", "metric_optimizer_2x2")
)
cfg$metrics <- unique(trimws(cfg$metrics[nzchar(cfg$metrics)]))
cfg$optimizers <- unique(trimws(cfg$optimizers[nzchar(cfg$optimizers)]))

valid_metrics <- c("var_simple", "simple_mpv", "simple_tgv", "simple_usefa",
                   "simple_usefad_nonphased", "simple_usefad_phased")
valid_optimizers <- c("topn", "simplemating_select")
bad_metrics <- setdiff(cfg$metrics, valid_metrics)
bad_optimizers <- setdiff(cfg$optimizers, valid_optimizers)
if (length(bad_metrics)) stop("Unknown OX_METRICS: ", paste(bad_metrics, collapse = ", "), call. = FALSE)
if (length(bad_optimizers)) stop("Unknown OX_OPTIMIZERS: ", paste(bad_optimizers, collapse = ", "), call. = FALSE)
if (any(vapply(cfg$metrics, is_simple_metric, logical(1))) &&
    !requireNamespace("SimpleMating", quietly = TRUE)) {
  stop("SimpleMating is required for the requested simple_* metrics.", call. = FALSE)
}

branch_grid <- expand.grid(metric = cfg$metrics, optimizer = cfg$optimizers, stringsAsFactors = FALSE)
branch_grid$branch <- paste(branch_grid$metric, branch_grid$optimizer, sep = "__")

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(name = names(cfg), value = vapply(cfg, function(x) paste(x, collapse = ","), character(1))),
          file.path("results", paste0(cfg$output_prefix, "_config.csv")), row.names = FALSE)

message("Metric x optimizer benchmark config:")
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(rep_id) {
  message("Running replicate ", rep_id, " of ", cfg$reps)
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
  SP <- SimParam$new(founder)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)

  base_pop <- newPop(founder, simParam = SP)
  base_parents <- make_dh_base_parents(base_pop, cfg$n_parents, SP)
  base_parents@id <- paste0("P0_", seq_len(nInd(base_parents)))
  base_best <- max(as.numeric(gv(base_parents)[, 1]))

  snp_map <- getSnpMap(simParam = SP)
  marker_names <- paste0("Chr", snp_map$chr, "_", snp_map$id)
  marker_map <- data.frame(marker = marker_names, chr = snp_map$chr, pos = snp_map$pos,
                           stringsAsFactors = FALSE)

  branches <- setNames(vector("list", nrow(branch_grid)), branch_grid$branch)
  for (b in branch_grid$branch) branches[[b]] <- base_parents

  metrics_out <- list()
  selections <- list()
  timings <- list()
  family_diagnostics <- list()
  parent_use <- list()
  effect_diagnostics <- list()

  for (r in seq_len(nrow(branch_grid))) {
    bg <- branch_grid[r, ]
    metrics_out[[paste(bg$branch, 0, sep = "_")]] <- population_metrics(
      branches[[bg$branch]], base_best, 0L, rep_id, bg$metric, bg$optimizer, bg$branch
    )
  }

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    for (r in seq_len(nrow(branch_grid))) {
      bg <- branch_grid[r, ]
      parent_pop <- branches[[bg$branch]]
      ids <- paste0("B", r, "_R", rep_id, "_C", cycle - 1L, "_P", seq_len(nInd(parent_pop)))
      parent_pop@id <- ids
      parent_values <- as.numeric(gv(parent_pop)[, 1])
      names(parent_values) <- ids
      geno <- pullSnpGeno(parent_pop, simParam = SP)
      colnames(geno) <- marker_names
      rownames(geno) <- ids
      pair_df <- build_parent_pairs(ids, include_self = FALSE)

      effects <- NULL
      beta <- rep(0, ncol(geno))
      names(beta) <- colnames(geno)
      if (needs_effects(bg$metric)) {
        elapsed_effect <- system.time({
          training_pop <- make_training_population(parent_pop, cfg$effect_training_n, SP)
          training_ids <- paste0("B", r, "_R", rep_id, "_C", cycle - 1L, "_T", seq_len(nInd(training_pop)))
          training_pop@id <- training_ids
          training_y <- as.numeric(gv(training_pop)[, 1])
          eff <- estimate_effects(training_pop, SP, marker_names, training_y, cfg)
          effects <- eff$effects
          beta <- eff$beta
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, cycle = cycle, metric = bg$metric, optimizer = bg$optimizer,
          branch = bg$branch, stage = "effects", elapsed_sec = elapsed_effect,
          stringsAsFactors = FALSE
        )
        effect_diag <- effect_scale_diagnostics(effects, geno, parent_values, cv_folds = 0L)
        effect_diag$rep <- rep_id
        effect_diag$cycle <- cycle
        effect_diag$metric <- bg$metric
        effect_diag$optimizer <- bg$optimizer
        effect_diag$branch <- bg$branch
        effect_diagnostics[[paste(rep_id, bg$branch, cycle, sep = "_")]] <- effect_diag
      }

      elapsed_score <- system.time({
        haplo <- if (bg$metric == "simple_usefad_phased") as_snp_haplo(parent_pop, SP, marker_names) else NULL
        cross_df <- score_crosses(bg$metric, parent_pop, parent_values, geno, marker_map, haplo,
                                  pair_df, effects, beta, cfg)
      })[["elapsed"]]
      timings[[length(timings) + 1L]] <- data.frame(
        rep = rep_id, cycle = cycle, metric = bg$metric, optimizer = bg$optimizer,
        branch = bg$branch, stage = "score", elapsed_sec = elapsed_score,
        stringsAsFactors = FALSE
      )

      elapsed_select <- system.time({
        selected <- if (bg$optimizer == "topn") {
          select_topn(cross_df, cfg$top_crosses)
        } else {
          select_simplemating(cross_df, cfg)
        }
      })[["elapsed"]]
      selected$rep <- rep_id
      selected$cycle <- cycle
      selected$metric <- bg$metric
      selected$optimizer <- bg$optimizer
      selected$branch <- bg$branch
      selections[[paste(rep_id, bg$branch, cycle, sep = "_")]] <- selected
      timings[[length(timings) + 1L]] <- data.frame(
        rep = rep_id, cycle = cycle, metric = bg$metric, optimizer = bg$optimizer,
        branch = bg$branch, stage = "select", elapsed_sec = elapsed_select,
        stringsAsFactors = FALSE
      )
      pus <- parent_use_summary(selected)
      pus$rep <- rep_id
      pus$cycle <- cycle
      pus$metric <- bg$metric
      pus$optimizer <- bg$optimizer
      pus$branch <- bg$branch
      parent_use[[paste(rep_id, bg$branch, cycle, sep = "_")]] <- pus

      progeny_result <- make_selected_progeny(parent_pop, selected, cfg$progeny_per_cross, SP)
      family_diag <- progeny_result$family_metrics
      family_diag$rep <- rep_id
      family_diag$cycle <- cycle
      family_diag$metric <- bg$metric
      family_diag$optimizer <- bg$optimizer
      family_diag$branch <- bg$branch
      family_diagnostics[[paste(rep_id, bg$branch, cycle, sep = "_")]] <- family_diag

      progeny <- progeny_result$pop
      g <- as.numeric(gv(progeny)[, 1])
      keep <- order(g, decreasing = TRUE)[seq_len(min(cfg$n_parents, length(g)))]
      next_pop <- progeny[keep]
      next_pop@id <- paste0("B", r, "_R", rep_id, "_C", cycle, "_P", seq_len(nInd(next_pop)))
      branches[[bg$branch]] <- next_pop

      metrics_out[[paste(bg$branch, cycle, sep = "_")]] <- population_metrics(
        next_pop, base_best, cycle, rep_id, bg$metric, bg$optimizer, bg$branch
      )
    }
  }

  list(
    metrics = bind_rows_fill(metrics_out),
    selections = bind_rows_fill(selections),
    timings = bind_rows_fill(timings),
    family_diagnostics = bind_rows_fill(family_diagnostics),
    parent_use = bind_rows_fill(parent_use),
    effect_diagnostics = bind_rows_fill(effect_diagnostics)
  )
})

metrics <- bind_rows_fill(lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
timings <- bind_rows_fill(lapply(rep_results, `[[`, "timings"))
family_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "family_diagnostics"))
parent_use <- bind_rows_fill(lapply(rep_results, `[[`, "parent_use"))
effect_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "effect_diagnostics"))

overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ metric + optimizer + cycle,
  metrics,
  mean
)

topn_baseline <- overall[overall$metric == "var_simple" & overall$optimizer == "topn",
                         c("cycle", "mean_gv", "max_gv", "top10_gv", "transgressive_rate"),
                         drop = FALSE]
vs_var_topn <- merge(overall, topn_baseline, by = "cycle", suffixes = c("", "_var_simple_topn"))
vs_var_topn$delta_mean_gv <- vs_var_topn$mean_gv - vs_var_topn$mean_gv_var_simple_topn
vs_var_topn$delta_max_gv <- vs_var_topn$max_gv - vs_var_topn$max_gv_var_simple_topn
vs_var_topn$delta_top10_gv <- vs_var_topn$top10_gv - vs_var_topn$top10_gv_var_simple_topn
vs_var_topn$delta_transgressive_rate <- (
  vs_var_topn$transgressive_rate - vs_var_topn$transgressive_rate_var_simple_topn
)

optimizer_effect <- reshape(
  overall[overall$cycle > 0, c("metric", "optimizer", "cycle", "mean_gv", "max_gv", "top10_gv")],
  idvar = c("metric", "cycle"),
  timevar = "optimizer",
  direction = "wide"
)
if (all(c("mean_gv.simplemating_select", "mean_gv.topn") %in% names(optimizer_effect))) {
  optimizer_effect$delta_mean_gv_smopt_vs_topn <- (
    optimizer_effect$mean_gv.simplemating_select - optimizer_effect$mean_gv.topn
  )
  optimizer_effect$delta_max_gv_smopt_vs_topn <- (
    optimizer_effect$max_gv.simplemating_select - optimizer_effect$max_gv.topn
  )
  optimizer_effect$delta_top10_gv_smopt_vs_topn <- (
    optimizer_effect$top10_gv.simplemating_select - optimizer_effect$top10_gv.topn
  )
}

timing_summary <- if (nrow(timings)) {
  stats::aggregate(elapsed_sec ~ metric + optimizer + stage, timings, mean)
} else data.frame()

family_summary <- if (nrow(family_diagnostics)) {
  stats::aggregate(
    cbind(pred_var, realized_var_gv, realized_max_gv, realized_top10_gv) ~ metric + optimizer + cycle,
    family_diagnostics,
    mean
  )
} else data.frame()

parent_use_summary_df <- if (nrow(parent_use)) {
  stats::aggregate(
    cbind(n_selected_crosses, n_unique_parents, max_parent_use, mean_parent_use) ~ metric + optimizer + cycle,
    parent_use,
    mean
  )
} else data.frame()

prefix <- file.path("results", cfg$output_prefix)
write.csv(metrics, paste0(prefix, "_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, paste0(prefix, "_overall.csv"), row.names = FALSE)
write.csv(vs_var_topn, paste0(prefix, "_comparison_vs_var_simple_topn.csv"), row.names = FALSE)
write.csv(optimizer_effect, paste0(prefix, "_optimizer_effect.csv"), row.names = FALSE)
write.csv(selections, paste0(prefix, "_selected_crosses.csv"), row.names = FALSE)
write.csv(timings, paste0(prefix, "_timing.csv"), row.names = FALSE)
write.csv(timing_summary, paste0(prefix, "_timing_summary.csv"), row.names = FALSE)
write.csv(family_diagnostics, paste0(prefix, "_family_diagnostics.csv"), row.names = FALSE)
write.csv(family_summary, paste0(prefix, "_family_summary.csv"), row.names = FALSE)
write.csv(parent_use, paste0(prefix, "_parent_use_by_rep.csv"), row.names = FALSE)
write.csv(parent_use_summary_df, paste0(prefix, "_parent_use_summary.csv"), row.names = FALSE)
write.csv(effect_diagnostics, paste0(prefix, "_effect_diagnostics.csv"), row.names = FALSE)

message("Metric x optimizer overall:")
print(overall)
message("Comparison versus var_simple + topn:")
print(vs_var_topn[order(vs_var_topn$cycle, vs_var_topn$metric, vs_var_topn$optimizer), , drop = FALSE])
message("Same metric optimizer effect:")
print(optimizer_effect)
message("Parent-use summary:")
print(parent_use_summary_df)
message("Timing summary:")
print(timing_summary)
