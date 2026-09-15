# AlphaSimR simulation diagnostic for var_simple vs PMV cross variance.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
suppressPackageStartupMessages({
  library(AlphaSimR)
})
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

pick_top_crosses <- function(df, method, top_k, rank_by = c("usefulness", "variance")) {
  rank_by <- match.arg(rank_by)
  score_col <- if (rank_by == "usefulness") "cross_usefulness" else "cross_var"
  if (!score_col %in% names(df)) stop("Missing ranking column: ", score_col, call. = FALSE)
  df <- df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_k)
  data.frame(
    parent1 = as.character(df$parent1),
    parent2 = as.character(df$parent2),
    method = method,
    rank_by = rank_by,
    predicted_score = df[[score_col]],
    predicted_var = df$cross_var,
    predicted_mean = df$cross_mean,
    predicted_usefulness = df$cross_usefulness,
    stringsAsFactors = FALSE
  )
}

pair_key <- function(parent1, parent2) {
  paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = " x ")
}

evaluate_crosses <- function(pop, cross_df, parent_gv, progeny_per_cross, sim_param) {
  out <- vector("list", nrow(cross_df))
  parent_ids <- pop@id
  for (i in seq_len(nrow(cross_df))) {
    p1 <- match(cross_df$parent1[i], parent_ids)
    p2 <- match(cross_df$parent2[i], parent_ids)
    f1 <- makeCross(pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
    dh <- makeDH(f1, nDH = progeny_per_cross, keepParents = FALSE, simParam = sim_param)
    g <- as.numeric(gv(dh)[, 1])
    p_gv <- parent_gv[c(cross_df$parent1[i], cross_df$parent2[i])]
    out[[i]] <- data.frame(
      parent1 = cross_df$parent1[i],
      parent2 = cross_df$parent2[i],
      method = cross_df$method[i],
      rank_by = cross_df$rank_by[i],
      predicted_score = cross_df$predicted_score[i],
      predicted_var = cross_df$predicted_var[i],
      predicted_mean = cross_df$predicted_mean[i],
      predicted_usefulness = cross_df$predicted_usefulness[i],
      realized_mean = mean(g),
      realized_var = stats::var(g),
      realized_sd = stats::sd(g),
      realized_best = max(g),
      top10_mean = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
      better_than_best_parent = mean(g > max(p_gv)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

summarize_methods <- function(evaluated) {
  methods <- unique(evaluated$method)
  do.call(rbind, lapply(methods, function(m) {
    x <- evaluated[evaluated$method == m, , drop = FALSE]
    data.frame(
      method = m,
      n_crosses = nrow(x),
      mean_predicted_var = mean(x$predicted_var),
      mean_realized_var = mean(x$realized_var),
      mean_realized_best = mean(x$realized_best),
      mean_top10 = mean(x$top10_mean),
      mean_transgressive_rate = mean(x$better_than_best_parent),
      stringsAsFactors = FALSE
    )
  }))
}

run_rep <- function(rep_id, cfg) {
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(
    nInd = cfg$n_founders,
    nChr = cfg$n_chr,
    segSites = cfg$seg_sites,
    genLen = cfg$gen_len
  )
  SP <- SimParam$new(founder)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = cfg$trait_var)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  pop_all <- newPop(founder, simParam = SP)
  parent_pop <- pop_all[seq_len(cfg$n_parents)]
  parent_ids <- paste0("P", seq_len(nInd(parent_pop)))
  parent_pop@id <- parent_ids
  
  snp_geno <- pullSnpGeno(parent_pop, simParam = SP)
  snp_map <- getSnpMap(simParam = SP)
  colnames(snp_geno) <- snp_map$id
  rownames(snp_geno) <- parent_ids
  snp_map <- data.frame(
    marker = snp_map$id,
    chr = snp_map$chr,
    pos = snp_map$pos,
    stringsAsFactors = FALSE
  )
  
  parent_gv <- as.numeric(gv(parent_pop)[, 1])
  names(parent_gv) <- parent_ids
  y_effects <- parent_gv
  if (cfg$effect_target == "phenotype") {
    parent_pop <- setPheno(parent_pop, h2 = cfg$h2, simParam = SP)
    y_effects <- as.numeric(pheno(parent_pop)[, 1])
  }
  effects <- estimate_marker_effects_ridge(snp_geno, y_effects, lambda = cfg$ridge_lambda)
  if (nrow(snp_geno) < ncol(snp_geno)) {
    warning(
      "Marker-effect training n (", nrow(snp_geno), ") is smaller than SNP marker count (", ncol(snp_geno), "). ",
      "For production use, include additional genotyped/phenotyped individuals for marker-effect estimation; ",
      "the crossing candidate parent set can remain smaller.",
      call. = FALSE
    )
  }
  
  t_var <- system.time({
    var_simple <- build_cross_data(
      pheno = parent_gv,
      geno_mat = snp_geno,
      candidate_ids = parent_ids,
      variance_method = "var_simple",
      include_self = FALSE,
      use_parallel = FALSE,
      trait_direction = "increase"
    )
  })
  t_simple <- system.time({
    pmv_simple <- build_cross_data(
      pheno = parent_gv,
      geno_mat = snp_geno,
      candidate_ids = parent_ids,
      variance_method = "pmv_simple",
      marker_effects_mcmc = effects,
      include_self = FALSE,
      use_parallel = FALSE,
      progeny = "DH",
      trait_direction = "increase",
      include_effect_uncertainty = FALSE
    )
  })
  t_pos <- system.time({
    pmv_pos <- build_cross_data(
      pheno = parent_gv,
      geno_mat = snp_geno,
      candidate_ids = parent_ids,
      variance_method = "pmv_pos",
      marker_effects_mcmc = effects,
      map = snp_map,
      include_self = FALSE,
      use_parallel = FALSE,
      progeny = "DH",
      map_function = "kosambi",
      pos_unit = "M",
      window_cM = cfg$pmv_window_cM,
      trait_direction = "increase",
      include_effect_uncertainty = FALSE
    )
  })
  
  var_simple$method <- "var_simple"
  pmv_simple$method <- "pmv_simple"
  pmv_pos$method <- "pmv_pos"
  
  pred_all <- rbind(
    data.frame(var_simple, rep = rep_id, stringsAsFactors = FALSE),
    data.frame(pmv_simple, rep = rep_id, stringsAsFactors = FALSE),
    data.frame(pmv_pos, rep = rep_id, stringsAsFactors = FALSE)
  )
  
  selected <- rbind(
    pick_top_crosses(var_simple, "var_simple", cfg$top_k, cfg$rank_by),
    pick_top_crosses(pmv_simple, "pmv_simple", cfg$top_k, cfg$rank_by),
    pick_top_crosses(pmv_pos, "pmv_pos", cfg$top_k, cfg$rank_by)
  )
  selected$key <- pair_key(selected$parent1, selected$parent2)
  selected <- selected[!duplicated(paste(selected$method, selected$key)), , drop = FALSE]
  evaluated <- evaluate_crosses(parent_pop, selected, parent_gv, cfg$progeny_per_cross, SP)
  evaluated$rep <- rep_id
  
  calibration_keys <- unique(unlist(lapply(split(selected$key, selected$method), identity)))
  calibration <- merge(
    unique(evaluated[, c("parent1", "parent2", "realized_var")]),
    pred_all[, c("parent1", "parent2", "method", "cross_var")],
    by = c("parent1", "parent2")
  )
  calibration$key <- pair_key(calibration$parent1, calibration$parent2)
  calibration <- calibration[calibration$key %in% calibration_keys, , drop = FALSE]
  calib_summary <- do.call(rbind, lapply(unique(calibration$method), function(m) {
    x <- calibration[calibration$method == m, , drop = FALSE]
    data.frame(
      rep = rep_id,
      method = m,
      cor_pred_realized_var = suppressWarnings(stats::cor(x$cross_var, x$realized_var)),
      stringsAsFactors = FALSE
    )
  }))
  
  list(
    pred_all = pred_all,
    evaluated = evaluated,
    method_summary = data.frame(rep = rep_id, summarize_methods(evaluated)),
    calibration_summary = calib_summary,
    timing = data.frame(
      rep = rep_id,
      method = c("var_simple", "pmv_simple", "pmv_pos"),
      elapsed_sec = c(t_var[["elapsed"]], t_simple[["elapsed"]], t_pos[["elapsed"]]),
      stringsAsFactors = FALSE
    ),
      marker_info = data.frame(
      rep = rep_id,
      marker_source = "snp_chip",
      n_effect_training_individuals = nrow(snp_geno),
      n_markers = ncol(snp_geno),
      n_qtl = cfg$n_chr * cfg$qtl_per_chr,
      effect_target = cfg$effect_target,
      ridge_lambda = cfg$ridge_lambda,
      pmv_window_cM = cfg$pmv_window_cM,
      recommendation = if (nrow(snp_geno) < ncol(snp_geno)) {
        "Effect training n is smaller than marker count; include additional genotyped/phenotyped individuals for marker-effect estimation when possible."
      } else {
        "Effect training n is at least marker count."
      },
      stringsAsFactors = FALSE
    )
  )
}

cfg <- list(
  seed = env_int("ALPHASIMR_SEED", 20260425L),
  reps = env_int("ALPHASIMR_REPS", 3L),
  n_founders = env_int("ALPHASIMR_N_FOUNDERS", 120L),
  n_parents = env_int("ALPHASIMR_N_PARENTS", 80L),
  n_chr = env_int("ALPHASIMR_N_CHR", 5L),
  seg_sites = env_int("ALPHASIMR_SEG_SITES", 1200L),
  qtl_per_chr = env_int("ALPHASIMR_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("ALPHASIMR_SNP_PER_CHR", 1000L),
  gen_len = env_num("ALPHASIMR_GEN_LEN", 1),
  trait_var = env_num("ALPHASIMR_TRAIT_VAR", 1),
  effect_target = env_chr("ALPHASIMR_EFFECT_TARGET", "gv"),
  h2 = env_num("ALPHASIMR_H2", 0.5),
  ridge_lambda = env_num("ALPHASIMR_RIDGE_LAMBDA", 1),
  pmv_window_cM = env_num("ALPHASIMR_PMV_WINDOW_CM", 10),
  rank_by = env_chr("ALPHASIMR_RANK_BY", "usefulness"),
  top_k = env_int("ALPHASIMR_TOP_K", 20L),
  progeny_per_cross = env_int("ALPHASIMR_PROGENY_PER_CROSS", 150L)
)

if (cfg$n_parents > cfg$n_founders) stop("ALPHASIMR_N_PARENTS cannot exceed ALPHASIMR_N_FOUNDERS.")
cfg$effect_target <- match.arg(cfg$effect_target, c("gv", "phenotype"))
cfg$rank_by <- match.arg(cfg$rank_by, c("usefulness", "variance"))
if (!load_pmv_cpp(verbose = TRUE)) warning("PMV C++ acceleration unavailable; PMV methods will use R fallbacks.")

dir.create("results", showWarnings = FALSE)
message("AlphaSimR PMV diagnostic config:")
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(i) {
  message("Running replicate ", i, " of ", cfg$reps)
  run_rep(i, cfg)
})

evaluated <- do.call(rbind, lapply(rep_results, `[[`, "evaluated"))
method_summary <- do.call(rbind, lapply(rep_results, `[[`, "method_summary"))
calibration_summary <- do.call(rbind, lapply(rep_results, `[[`, "calibration_summary"))
timing <- do.call(rbind, lapply(rep_results, `[[`, "timing"))
marker_info <- do.call(rbind, lapply(rep_results, `[[`, "marker_info"))

overall_summary <- do.call(rbind, lapply(unique(method_summary$method), function(m) {
  x <- method_summary[method_summary$method == m, , drop = FALSE]
  data.frame(
    method = m,
    reps = length(unique(x$rep)),
    mean_realized_var = mean(x$mean_realized_var),
    mean_realized_best = mean(x$mean_realized_best),
    mean_top10 = mean(x$mean_top10),
    mean_transgressive_rate = mean(x$mean_transgressive_rate),
    stringsAsFactors = FALSE
  )
}))

utils::write.csv(evaluated, file.path("results", "alphasimr_cross_evaluations.csv"), row.names = FALSE)
utils::write.csv(method_summary, file.path("results", "alphasimr_method_summary_by_rep.csv"), row.names = FALSE)
utils::write.csv(calibration_summary, file.path("results", "alphasimr_calibration_summary.csv"), row.names = FALSE)
utils::write.csv(timing, file.path("results", "alphasimr_timing.csv"), row.names = FALSE)
utils::write.csv(marker_info, file.path("results", "alphasimr_marker_info.csv"), row.names = FALSE)
utils::write.csv(overall_summary, file.path("results", "alphasimr_overall_summary.csv"), row.names = FALSE)

message("Overall method summary:")
print(overall_summary)
message("Timing summary:")
print(stats::aggregate(elapsed_sec ~ method, timing, mean))
