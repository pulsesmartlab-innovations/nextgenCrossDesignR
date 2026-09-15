# Multi-cycle AlphaSimR assessment of cross prediction metrics.

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

pair_key <- function(parent1, parent2) {
  paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = " x ")
}

merge_pop_list <- function(pop_list) {
  pop_list <- Filter(Negate(is.null), pop_list)
  if (!length(pop_list)) stop("No populations to merge.", call. = FALSE)
  Reduce(c, pop_list)
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

population_metrics <- function(pop, base_best, cycle, rep_id, method) {
  g <- as.numeric(gv(pop)[, 1])
  data.frame(
    rep = rep_id,
    method = method,
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

as_snp_haplo <- function(pop, sim_param, marker_names) {
  haplo <- pullSnpHaplo(pop, simParam = sim_param)
  colnames(haplo) <- marker_names
  rownames(haplo) <- as.vector(rbind(paste0(pop@id, "_HapA"), paste0(pop@id, "_HapB")))
  haplo
}

build_metric_crosses <- function(method, parent_gv, geno, haplo, ids, effects, map, ld_cache, cfg,
                                 effect_training_n) {
  build_cross_data(
    pheno = parent_gv,
    geno_mat = geno,
    candidate_ids = ids,
    variance_method = method,
    marker_effects_mcmc = effects,
    map = if (method %in% c("pmv", "pmv_pos", "pmv_pos_full", "pmv_advanced",
                            "ugv_pos", "ohv")) map else NULL,
    haplo_mat = if (method %in% c("pmv_pos_full", "pmv_advanced")) haplo else NULL,
    ld_mat = if (method %in% c("pmv_advanced", "pmv_nopos")) ld_cache else NULL,
    include_self = FALSE,
    use_parallel = FALSE,
    progeny = "DH",
    map_function = cfg$map_function,
    tau = cfg$tau,
    h2 = cfg$h2,
    trait_direction = "increase",
    pos_unit = "M",
    window_cM = cfg$pmv_window_cM,
    mean_source = cfg$mean_source,
    effect_training_n = effect_training_n,
    mean_underpowered_ratio = cfg$mean_underpowered_ratio,
    include_effect_uncertainty = cfg$include_effect_uncertainty,
    window_size = cfg$ld_window_size,
    block_size = cfg$ld_block_size,
    shrink_lambda = cfg$shrink_lambda,
    dense_ld_screen_n = if (method %in% c("pmv_pos_full", "pmv_advanced", "pmv_nopos") &&
                              cfg$ld_screen_multiplier > 0) {
      cfg$top_crosses * cfg$ld_screen_multiplier
    } else {
      NULL
    },
    dense_ld_screen_by = cfg$ld_screen_by
  )
}

select_crosses <- function(cross_df, top_crosses) {
  if (!"cross_usefulness" %in% names(cross_df)) {
    cross_df <- add_usefulness_criterion(cross_df, tau = 0.10, h2 = 1.0, trait_direction = "increase")
  }
  cross_df <- cross_df[order(cross_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  utils::head(cross_df, top_crosses)
}

make_selected_progeny <- function(parent_pop, selected, progeny_per_cross, sim_param) {
  parent_ids <- parent_pop@id
  progeny <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- match(selected$parent1[i], parent_ids)
    p2 <- match(selected$parent2[i], parent_ids)
    f1 <- makeCross(parent_pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
    progeny[[i]] <- makeDH(f1, nDH = progeny_per_cross, keepParents = FALSE, simParam = sim_param)
  }
  merge_pop_list(progeny)
}

cfg <- list(
  seed = env_int("MC_SEED", 20260425L),
  reps = env_int("MC_REPS", 5L),
  cycles = env_int("MC_CYCLES", 4L),
  n_founders = env_int("MC_N_FOUNDERS", 120L),
  n_parents = env_int("MC_N_PARENTS", 80L),
  effect_training_n = env_int("MC_EFFECT_TRAINING_N", env_int("MC_N_PARENTS", 80L)),
  n_chr = env_int("MC_N_CHR", 5L),
  seg_sites = env_int("MC_SEG_SITES", 1200L),
  qtl_per_chr = env_int("MC_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("MC_SNP_PER_CHR", 1000L),
  top_crosses = env_int("MC_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("MC_PROGENY_PER_CROSS", 40L),
  effect_method = env_chr("MC_EFFECT_METHOD", "ridge_posterior"),
  effect_draws = env_int("MC_EFFECT_DRAWS", 20L),
  effect_cv_folds = env_int("MC_EFFECT_CV_FOLDS", 5L),
  ridge_lambda = env_num("MC_RIDGE_LAMBDA", 1),
  prior_genetic_var_scale = env_num("MC_PRIOR_GENETIC_VAR_SCALE", 1),
  map_function = env_chr("MC_MAP_FUNCTION", "haldane"),
  pmv_window_cM = env_num("MC_PMV_WINDOW_CM", 10),
  ld_window_size = env_int("MC_LD_WINDOW_SIZE", 250L),
  ld_block_size = env_int("MC_LD_BLOCK_SIZE", 2000L),
  ld_source = env_chr("MC_LD_SOURCE", "training"),
  ld_screen_multiplier = env_int("MC_LD_SCREEN_MULTIPLIER", 50L),
  ld_screen_by = env_chr("MC_LD_SCREEN_BY", "union_usefulness"),
  include_effect_uncertainty = as.logical(env_int("MC_INCLUDE_EFFECT_UNCERTAINTY", 1L)),
  shrink_lambda = env_num("MC_SHRINK_LAMBDA", 0.05),
  tau = env_num("MC_TAU", 0.10),
  h2 = env_num("MC_H2", 1.0),
  mean_source = env_chr("MC_MEAN_SOURCE", "auto"),
  mean_underpowered_ratio = env_num("MC_MEAN_UNDERPOWERED_RATIO", 1),
  methods = strsplit(env_chr("MC_METHODS", "var_simple,pmv_pos_full,pmv_advanced,pmv_pos,pmv_nopos,ugv,ugv_pos,gedv,ohv,uc,ucpc"), ",", fixed = TRUE)[[1]]
)
cfg$methods <- trimws(cfg$methods)

dir.create("results", showWarnings = FALSE)
message("Multi-cycle AlphaSimR config:")
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(rep_id) {
  message("Running replicate ", rep_id, " of ", cfg$reps)
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
  SP <- SimParam$new(founder)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  base_pop <- newPop(founder, simParam = SP)
  base_parents <- base_pop[seq_len(cfg$n_parents)]
  base_parents@id <- paste0("P0_", seq_len(nInd(base_parents)))
  base_best <- max(as.numeric(gv(base_parents)[, 1]))

  branches <- setNames(vector("list", length(cfg$methods)), cfg$methods)
  for (m in cfg$methods) branches[[m]] <- base_parents

  metrics <- list()
  selections <- list()
  timings <- list()
  effect_diagnostics <- list()

  for (method in cfg$methods) {
    metrics[[paste(method, 0, sep = "_")]] <- population_metrics(
      branches[[method]], base_best, cycle = 0L, rep_id = rep_id, method = method
    )
  }

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    for (method in cfg$methods) {
      parent_pop <- branches[[method]]
      ids <- paste0(method, "_C", cycle - 1L, "_P", seq_len(nInd(parent_pop)))
      parent_pop@id <- ids
      parent_gv <- as.numeric(gv(parent_pop)[, 1])
      names(parent_gv) <- ids
      geno <- pullSnpGeno(parent_pop, simParam = SP)
      snp_map <- getSnpMap(simParam = SP)
      colnames(geno) <- snp_map$id
      rownames(geno) <- ids
      haplo <- if (method %in% c("pmv_pos_full", "pmv_advanced")) {
        as_snp_haplo(parent_pop, SP, snp_map$id)
      } else {
        NULL
      }
      map <- data.frame(marker = snp_map$id, chr = snp_map$chr, pos = snp_map$pos,
                        stringsAsFactors = FALSE)

      effect_methods <- c("pmv", "pmv_simple", "pmv_pos", "pmv_pos_full", "pmv_nopos", "pmv_advanced",
                          "ugv", "ugv_pos", "gedv", "ohv")
      needs_effects <- method %in% effect_methods
      effects <- NULL
      training_geno <- geno
      if (needs_effects) {
        training_pop <- make_training_population(parent_pop, cfg$effect_training_n, SP)
        training_ids <- paste0(method, "_C", cycle - 1L, "_T", seq_len(nInd(training_pop)))
        training_pop@id <- training_ids
        training_gv <- as.numeric(gv(training_pop)[, 1])
        training_geno <- pullSnpGeno(training_pop, simParam = SP)
        colnames(training_geno) <- snp_map$id
        rownames(training_geno) <- training_ids
        effects <- estimate_marker_effect_draws(
          geno_mat = training_geno,
          y = training_gv,
          method = cfg$effect_method,
          n_draws = cfg$effect_draws,
          lambda = cfg$ridge_lambda,
          prior_genetic_var = stats::var(training_gv) * cfg$prior_genetic_var_scale
        )
        effect_diag <- effect_scale_diagnostics(
          effects,
          training_geno,
          training_gv,
          cv_folds = cfg$effect_cv_folds,
          cv_lambda = cfg$ridge_lambda
        )
        effect_diag$rep <- rep_id
        effect_diag$method <- method
        effect_diag$cycle <- cycle
        effect_diagnostics[[paste(rep_id, method, cycle, sep = "_")]] <- effect_diag
      }

      ld_cache <- NULL
      if (method == "pmv_nopos" || (method == "pmv_advanced" && is.null(haplo))) {
        elapsed <- system.time({
          ld_geno <- if (identical(tolower(cfg$ld_source), "training")) training_geno else geno
          ld_cache <- calculate_ld_matrix_efficient(
            geno_mat = ld_geno,
            method = "cor",
            window_size = cfg$ld_window_size,
            block_size = cfg$ld_block_size,
            use_parallel = FALSE,
            chr = map$chr,
            min_maf = 0,
            verbose = FALSE
          )
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "ld_cache",
          elapsed_sec = elapsed, stringsAsFactors = FALSE
        )
      }

      elapsed <- system.time({
        cross_df <- build_metric_crosses(
          method, parent_gv, geno, haplo, ids, effects, map, ld_cache, cfg,
          effect_training_n = if (needs_effects) nrow(training_geno) else nrow(geno)
        )
      })[["elapsed"]]
      timings[[length(timings) + 1L]] <- data.frame(
        rep = rep_id, method = method, cycle = cycle, stage = "score",
        elapsed_sec = elapsed, stringsAsFactors = FALSE
      )
      selected <- select_crosses(cross_df, cfg$top_crosses)
      selected$rep <- rep_id
      selected$method <- method
      selected$cycle <- cycle
      selected$key <- pair_key(selected$parent1, selected$parent2)
      selections[[paste(rep_id, method, cycle, sep = "_")]] <- selected

      progeny <- make_selected_progeny(parent_pop, selected, cfg$progeny_per_cross, SP)
      g <- as.numeric(gv(progeny)[, 1])
      keep <- order(g, decreasing = TRUE)[seq_len(min(cfg$n_parents, length(g)))]
      next_pop <- progeny[keep]
      next_pop@id <- paste0(method, "_C", cycle, "_P", seq_len(nInd(next_pop)))
      branches[[method]] <- next_pop

      metrics[[paste(method, cycle, sep = "_")]] <- population_metrics(
        next_pop, base_best, cycle = cycle, rep_id = rep_id, method = method
      )
    }
  }

  timing_df <- do.call(rbind, timings)

  list(
    metrics = do.call(rbind, metrics),
    selections = bind_rows_fill(selections),
    timings = timing_df,
    effect_diagnostics = bind_rows_fill(effect_diagnostics)
  )
})

metrics <- do.call(rbind, lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
timings <- do.call(rbind, lapply(rep_results, `[[`, "timings"))
effect_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "effect_diagnostics"))

overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ method + cycle,
  metrics,
  mean
)
timing_summary <- stats::aggregate(elapsed_sec ~ method + stage, timings, mean)

write.csv(metrics, file.path("results", "multicycle_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, file.path("results", "multicycle_metrics_overall.csv"), row.names = FALSE)
write.csv(selections, file.path("results", "multicycle_selected_crosses.csv"), row.names = FALSE)
write.csv(timings, file.path("results", "multicycle_timing.csv"), row.names = FALSE)
write.csv(timing_summary, file.path("results", "multicycle_timing_summary.csv"), row.names = FALSE)
write.csv(effect_diagnostics, file.path("results", "multicycle_effect_diagnostics.csv"), row.names = FALSE)

message("Multi-cycle overall:")
print(overall)
message("Timing summary:")
print(timing_summary)
if (nrow(effect_diagnostics)) {
  message("Effect diagnostics:")
  print(stats::aggregate(
    cbind(effect_sd, parent_value_sd, marker_score_sd, marker_score_cor,
          marker_score_slope, marker_score_cv_cor, marker_score_cv_slope) ~ method + cycle,
    effect_diagnostics,
    mean,
    na.rm = TRUE
  ))
}
