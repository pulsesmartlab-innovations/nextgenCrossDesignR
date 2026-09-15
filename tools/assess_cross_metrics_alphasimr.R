# AlphaSimR assessment of available cross prediction metrics.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
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

env_logical <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  tolower(value) %in% c("1", "true", "yes", "y")
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
      metric_target = cross_df$metric_target[i],
      rank_score = cross_df$rank_score[i],
      predicted_mean = cross_df$cross_mean[i],
      predicted_var = cross_df$cross_var[i],
      predicted_usefulness = cross_df$cross_usefulness[i],
      realized_mean = mean(g),
      realized_var = stats::var(g),
      realized_best = max(g),
      top10_mean = mean(utils::head(sort(g, decreasing = TRUE), max(1L, ceiling(0.10 * length(g))))),
      better_than_best_parent = mean(g > max(p_gv)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

tag_metric <- function(df, method, target) {
  df$method <- method
  df$metric_target <- target
  if (!"cross_usefulness" %in% names(df)) {
    df <- add_usefulness_criterion(df, tau = 0.10, h2 = 1.0, trait_direction = "increase")
  }
  df
}

top_by_usefulness <- function(df, top_k) {
  df <- df[order(df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_k)
  df$rank_score <- df$cross_usefulness
  df
}

top_by_score <- function(df, top_k, score_col) {
  df <- df[order(df[[score_col]], decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_k)
  df$rank_score <- df[[score_col]]
  df
}

metric_cols <- c("parent1", "parent2", "cross_mean", "cross_mean_marker_effects",
                 "cross_var", "cross_sigma", "cross_usefulness", "mean_source",
                 "score", "method", "metric_target", "rank_score", "rep")

standardize_metric_df <- function(df) {
  missing <- setdiff(metric_cols, names(df))
  for (nm in missing) df[[nm]] <- NA
  df[, metric_cols, drop = FALSE]
}

metric_summary <- function(evaluated) {
  do.call(rbind, lapply(unique(evaluated$method), function(m) {
    x <- evaluated[evaluated$method == m, , drop = FALSE]
    data.frame(
      method = m,
      metric_target = x$metric_target[1],
      n_crosses = nrow(x),
      mean_rank_score = mean(x$rank_score),
      mean_predicted_var = mean(x$predicted_var),
      mean_predicted_usefulness = mean(x$predicted_usefulness),
      mean_realized_var = mean(x$realized_var),
      mean_realized_best = mean(x$realized_best),
      mean_top10 = mean(x$top10_mean),
      mean_transgressive_rate = mean(x$better_than_best_parent),
      stringsAsFactors = FALSE
    )
  }))
}

calibration_summary <- function(pred, evaluated) {
  eval_unique <- unique(evaluated[, c("parent1", "parent2", "realized_var", "realized_best", "top10_mean")])
  pred$ranking_score <- ifelse(is.finite(as.numeric(pred$score)),
                               as.numeric(pred$score),
                               as.numeric(pred$cross_usefulness))
  merged <- merge(eval_unique, pred[, c("parent1", "parent2", "method", "cross_var", "ranking_score")],
                  by = c("parent1", "parent2"))
  do.call(rbind, lapply(unique(merged$method), function(m) {
    x <- merged[merged$method == m, , drop = FALSE]
    data.frame(
      method = m,
      cor_var_to_realized_var = suppressWarnings(cor(x$cross_var, x$realized_var)),
      cor_usefulness_to_best = suppressWarnings(cor(x$ranking_score, x$realized_best)),
      cor_usefulness_to_top10 = suppressWarnings(cor(x$ranking_score, x$top10_mean)),
      stringsAsFactors = FALSE
    )
  }))
}

metric_requires_map <- function(method) {
  method %in% c("pmv_pos", "ugv_pos", "pmv", "ohv")
}

build_metric_crosses <- function(base_method, parent_gv, geno, ids, effects, map, cfg, ld_mat = NULL) {
  build_cross_data(parent_gv, geno, ids, base_method,
                   marker_effects_mcmc = effects,
                   map = if (metric_requires_map(base_method)) map else NULL,
                   ld_mat = ld_mat,
                   include_self = FALSE,
                   use_parallel = FALSE,
                   progeny = "DH",
                   tau = cfg$tau,
                   h2 = cfg$h2,
                   trait_direction = "increase",
                   pos_unit = "M",
                   window_cM = cfg$pmv_window_cM,
                   method_varPMV = cfg$method_varPMV,
                   method_ld = cfg$ld_method,
                   mean_source = cfg$mean_source,
                   effect_training_n = nrow(geno),
                   mean_underpowered_ratio = cfg$mean_underpowered_ratio,
                   include_effect_uncertainty = FALSE,
                   window_size = cfg$ld_window_size,
                   block_size = cfg$ld_block_size,
                   min_maf = cfg$ld_min_maf,
                   dense_ld_screen_n = if (base_method %in% c("pmv_nopos", "pmv_advanced") &&
                                             cfg$ld_screen_multiplier > 0) {
                     cfg$top_k * cfg$ld_screen_multiplier
                   } else {
                     NULL
                   },
                   dense_ld_screen_by = cfg$ld_screen_by)
}

add_cross_inbreeding <- function(df, G, ids) {
  idx1 <- match(df$parent1, ids)
  idx2 <- match(df$parent2, ids)
  df$inbreeding <- 0.5 * (diag(G)[idx1] + diag(G)[idx2]) - G[cbind(idx1, idx2)]
  df
}

compute_metric <- function(method, parent_gv, geno, ids, effects, map, cfg, G, ld_mat = NULL) {
  base_method <- if (method %in% c("uc_pmv_pos", "EGSI", "EGG-C", "ETS", "CSSR")) {
    cfg$score_base_variance
  } else {
    method
  }
  out <- build_metric_crosses(base_method, parent_gv, geno, ids, effects, map, cfg,
                              ld_mat = if (base_method %in% c("pmv_nopos", "pmv_advanced")) ld_mat else NULL)
  if (method == "uc_pmv_pos") {
    out$cross_var_source <- "pmv_pos"
  } else if (method == "EGSI") {
    out <- add_cross_inbreeding(out, G, ids)
    out <- calculate_egsi(out, selection_intensity = cfg$egsi_selection_percent)
    out$score <- out$EGSI
  } else if (method == "EGG-C") {
    out <- add_cross_inbreeding(out, G, ids)
    out <- calculate_egg_c(out)
    out$score <- out$EGGC
  } else if (method == "ETS") {
    out <- calculate_ets(out, geno, effects, ids, threshold = 0)
    out$score <- out$ETS
  } else if (method == "CSSR") {
    out <- calculate_cssr(out, G, ids, selection_intensity = cfg$tau, lambda = cfg$cssr_lambda)
    out$score <- out$CSSR
  }
  out
}

cfg <- list(
  seed = env_int("ASSESS_SEED", 20260425L),
  reps = env_int("ASSESS_REPS", 2L),
  n_founders = env_int("ASSESS_N_FOUNDERS", 120L),
  n_parents = env_int("ASSESS_N_PARENTS", 80L),
  n_chr = env_int("ASSESS_N_CHR", 5L),
  seg_sites = env_int("ASSESS_SEG_SITES", 1200L),
  qtl_per_chr = env_int("ASSESS_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("ASSESS_SNP_PER_CHR", 1000L),
  progeny_per_cross = env_int("ASSESS_PROGENY_PER_CROSS", 150L),
  top_k = env_int("ASSESS_TOP_K", 20L),
  ridge_lambda = env_num("ASSESS_RIDGE_LAMBDA", 1),
  pmv_window_cM = env_num("ASSESS_PMV_WINDOW_CM", 10),
  tau = env_num("ASSESS_TAU", 0.10),
  h2 = env_num("ASSESS_H2", 1.0),
  mean_source = env_chr("ASSESS_MEAN_SOURCE", "auto"),
  mean_underpowered_ratio = env_num("ASSESS_MEAN_UNDERPOWERED_RATIO", 1),
  methods = env_chr("ASSESS_METHODS", "ready"),
  include_legacy_uc = env_logical("ASSESS_INCLUDE_LEGACY_UC", FALSE),
  score_base_variance = env_chr("ASSESS_SCORE_BASE_VARIANCE", "pmv_pos"),
  method_varPMV = env_chr("ASSESS_METHOD_VARPMV", "fast"),
  ld_method = env_chr("ASSESS_LD_METHOD", "cor"),
  ld_window_size = env_int("ASSESS_LD_WINDOW_SIZE", 250L),
  ld_block_size = env_int("ASSESS_LD_BLOCK_SIZE", 2000L),
  ld_min_maf = env_num("ASSESS_LD_MIN_MAF", 0),
  ld_screen_multiplier = env_int("ASSESS_LD_SCREEN_MULTIPLIER", 0L),
  ld_screen_by = env_chr("ASSESS_LD_SCREEN_BY", "pmv_simple_usefulness"),
  egsi_selection_percent = env_num("ASSESS_EGSI_SELECTION_PERCENT", 10),
  cssr_lambda = env_num("ASSESS_CSSR_LAMBDA", 1)
)

ready_methods <- c("var_simple", "pmv_simple", "pmv_pos", "ugv", "ugv_pos", "gedv",
                   "ohv", "uc_pmv_pos", "EGSI", "EGG-C", "ETS", "CSSR")
slow_methods <- c("pmv_nopos", "pmv_advanced", "pmv")
legacy_uc_methods <- c("uc", "ucpc")
methods <- switch(
  cfg$methods,
  ready = ready_methods,
  all = c(ready_methods, slow_methods, legacy_uc_methods),
  slow = slow_methods,
  legacy_uc = legacy_uc_methods,
  strsplit(cfg$methods, ",", fixed = TRUE)[[1]]
)
methods <- trimws(methods)
if (cfg$include_legacy_uc) methods <- unique(c(methods, legacy_uc_methods))
targets <- c(var_simple = "relationship_diversity_proxy",
             pmv_simple = "effect_weighted_variance_no_recombination",
             pmv_pos = "effect_weighted_variance_with_recombination",
             ugv = "usefulness_genetic_variance_approximation",
             ugv_pos = "ugv_with_recombination_structure",
             gedv = "genetic_diversity_expected_heterozygosity",
             ohv = "optimal_haploid_value_block_assembly",
             uc_pmv_pos = "usefulness_from_recombination_aware_pmv",
             EGSI = "expected_genetic_superiority_index",
             `EGG-C` = "expected_gain_per_unit_coancestry",
             ETS = "expected_transgressive_segregation",
             CSSR = "cross_selection_score_response",
             pmv_nopos = "pmv_using_ld_without_positions",
             pmv_advanced = "pmv_combining_position_and_ld",
             pmv = "varPMV_original_driver",
             uc = "legacy_relationship_usefulness",
             ucpc = "legacy_usefulness_parental_contribution")

if (!load_pmv_cpp(verbose = TRUE)) warning("PMV C++ acceleration unavailable; PMV methods will use R fallbacks.")
dir.create("results", showWarnings = FALSE)
message("Cross metric assessment config:")
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(rep_id) {
  message("Running replicate ", rep_id, " of ", cfg$reps)
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
  SP <- SimParam$new(founder)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  pop_all <- newPop(founder, simParam = SP)
  parent_pop <- pop_all[seq_len(cfg$n_parents)]
  ids <- paste0("P", seq_len(nInd(parent_pop)))
  parent_pop@id <- ids
  
  parent_gv <- as.numeric(gv(parent_pop)[, 1])
  names(parent_gv) <- ids
  geno <- pullSnpGeno(parent_pop, simParam = SP)
  snp_map <- getSnpMap(simParam = SP)
  colnames(geno) <- snp_map$id
  rownames(geno) <- ids
  map <- data.frame(marker = snp_map$id, chr = snp_map$chr, pos = snp_map$pos,
                    stringsAsFactors = FALSE)
  G <- VanRadenKin(geno)
  effects <- estimate_marker_effects_ridge(geno, parent_gv, cfg$ridge_lambda)
  if (nrow(geno) < ncol(geno)) {
    warning("Marker-effect training n (", nrow(geno), ") is smaller than marker count (", ncol(geno),
            "); include additional training individuals for marker-effect estimation when possible.",
            call. = FALSE)
  }
  
  timings <- list()
  ld_cache <- NULL
  if (any(methods %in% c("pmv_nopos", "pmv_advanced"))) {
    timings[["ld_cache"]] <- system.time({
      ld_cache <- calculate_ld_matrix_efficient(
        geno_mat = geno,
        method = cfg$ld_method,
        window_size = cfg$ld_window_size,
        block_size = cfg$ld_block_size,
        use_parallel = FALSE,
        chr = map$chr,
        min_maf = cfg$ld_min_maf,
        verbose = FALSE
      )
    })[["elapsed"]]
  }
  pred_list <- lapply(methods, function(m) {
    timings[[m]] <<- system.time({
      out <- compute_metric(m, parent_gv, geno, ids, effects, map, cfg, G, ld_mat = ld_cache)
    })[["elapsed"]]
    tag_metric(out, m, unname(targets[m]))
  })
  pred_all <- do.call(rbind, lapply(pred_list, standardize_metric_df))
  pred_all$rep <- rep_id
  
  selected_list <- lapply(seq_along(pred_list), function(i) {
    m <- methods[[i]]
    if ("score" %in% names(pred_list[[i]]) && m %in% c("EGSI", "EGG-C", "ETS", "CSSR")) {
      top_by_score(pred_list[[i]], cfg$top_k, "score")
    } else if (m == "uc" && "uc" %in% names(pred_list[[i]])) {
      top_by_score(pred_list[[i]], cfg$top_k, "uc")
    } else if (m == "ucpc" && "ucpc" %in% names(pred_list[[i]])) {
      top_by_score(pred_list[[i]], cfg$top_k, "ucpc")
    } else {
      top_by_usefulness(pred_list[[i]], cfg$top_k)
    }
  })
  selected <- do.call(rbind, lapply(selected_list, standardize_metric_df))
  selected$key <- pair_key(selected$parent1, selected$parent2)
  selected <- selected[!duplicated(paste(selected$method, selected$key)), , drop = FALSE]
  evaluated <- evaluate_crosses(parent_pop, selected, parent_gv, cfg$progeny_per_cross, SP)
  evaluated$rep <- rep_id
  
  list(
    pred_all = pred_all,
    evaluated = evaluated,
    metric_summary = data.frame(rep = rep_id, metric_summary(evaluated)),
    calibration_summary = data.frame(rep = rep_id, calibration_summary(pred_all, evaluated)),
    timing = data.frame(rep = rep_id, method = names(timings),
                        elapsed_sec = as.numeric(timings), stringsAsFactors = FALSE),
    metadata = data.frame(rep = rep_id, n_effect_training_individuals = nrow(geno),
                          n_markers = ncol(geno), n_qtl = cfg$n_chr * cfg$qtl_per_chr,
                          recommendation = if (nrow(geno) < ncol(geno)) {
                            "Use additional genotyped/phenotyped training individuals for SNP effects; candidate parent set can remain smaller."
                          } else {
                            "Training n is at least marker count."
                          },
                          stringsAsFactors = FALSE)
  )
})

metric_by_rep <- do.call(rbind, lapply(rep_results, `[[`, "metric_summary"))
calibration_by_rep <- do.call(rbind, lapply(rep_results, `[[`, "calibration_summary"))
timing <- do.call(rbind, lapply(rep_results, `[[`, "timing"))
metadata <- do.call(rbind, lapply(rep_results, `[[`, "metadata"))
evaluated <- do.call(rbind, lapply(rep_results, `[[`, "evaluated"))

overall <- do.call(rbind, lapply(unique(metric_by_rep$method), function(m) {
  x <- metric_by_rep[metric_by_rep$method == m, , drop = FALSE]
  data.frame(
    method = m,
    metric_target = x$metric_target[1],
    reps = length(unique(x$rep)),
    mean_realized_var = mean(x$mean_realized_var),
    mean_realized_best = mean(x$mean_realized_best),
    mean_top10 = mean(x$mean_top10),
    mean_transgressive_rate = mean(x$mean_transgressive_rate),
    stringsAsFactors = FALSE
  )
}))

write.csv(overall, file.path("results", "cross_metric_assessment_overall.csv"), row.names = FALSE)
write.csv(metric_by_rep, file.path("results", "cross_metric_assessment_by_rep.csv"), row.names = FALSE)
write.csv(calibration_by_rep, file.path("results", "cross_metric_assessment_calibration.csv"), row.names = FALSE)
write.csv(timing, file.path("results", "cross_metric_assessment_timing.csv"), row.names = FALSE)
write.csv(metadata, file.path("results", "cross_metric_assessment_metadata.csv"), row.names = FALSE)
write.csv(evaluated, file.path("results", "cross_metric_assessment_crosses.csv"), row.names = FALSE)

message("Overall assessment:")
print(overall)
message("Timing:")
print(stats::aggregate(elapsed_sec ~ method, timing, mean))
