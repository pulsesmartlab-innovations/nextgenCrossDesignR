# Compare exhaustive dense-LD PMV with two-stage screened dense-LD PMV.

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

pair_key <- function(parent1, parent2) {
  paste(pmin(parent1, parent2), pmax(parent1, parent2), sep = " x ")
}

select_top <- function(df, top_k) {
  df <- df[order(df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  utils::head(df, top_k)
}

evaluate_unique_pairs <- function(pop, selected_pairs, parent_gv, progeny_per_cross, sim_param) {
  pair_df <- unique(selected_pairs[, c("parent1", "parent2"), drop = FALSE])
  pair_df$key <- pair_key(pair_df$parent1, pair_df$parent2)
  out <- vector("list", nrow(pair_df))
  parent_ids <- pop@id
  for (i in seq_len(nrow(pair_df))) {
    p1 <- match(pair_df$parent1[i], parent_ids)
    p2 <- match(pair_df$parent2[i], parent_ids)
    f1 <- makeCross(pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
    dh <- makeDH(f1, nDH = progeny_per_cross, keepParents = FALSE, simParam = sim_param)
    g <- as.numeric(gv(dh)[, 1])
    p_gv <- parent_gv[c(pair_df$parent1[i], pair_df$parent2[i])]
    out[[i]] <- data.frame(
      key = pair_df$key[i],
      parent1 = pair_df$parent1[i],
      parent2 = pair_df$parent2[i],
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

summarize_selected <- function(selected, evaluated, full_keys_by_method) {
  selected$key <- pair_key(selected$parent1, selected$parent2)
  merged <- merge(selected, evaluated, by = c("key", "parent1", "parent2"))
  do.call(rbind, lapply(unique(merged$config), function(cfg) {
    x <- merged[merged$config == cfg, , drop = FALSE]
    base_method <- sub("_screen.*$|_full$", "", cfg)
    full_keys <- full_keys_by_method[[base_method]]
    data.frame(
      config = cfg,
      method = base_method,
      screen_multiplier = x$screen_multiplier[1],
      n_selected = nrow(x),
      overlap_with_full = if (is.null(full_keys)) NA_real_ else mean(x$key %in% full_keys),
      mean_predicted_var = mean(x$cross_var),
      mean_predicted_usefulness = mean(x$cross_usefulness),
      mean_realized_var = mean(x$realized_var),
      mean_realized_best = mean(x$realized_best),
      mean_top10 = mean(x$top10_mean),
      mean_transgressive_rate = mean(x$better_than_best_parent),
      stringsAsFactors = FALSE
    )
  }))
}

cfg <- list(
  seed = env_int("SCREEN_SEED", 20260425L),
  reps = env_int("SCREEN_REPS", 10L),
  n_founders = env_int("SCREEN_N_FOUNDERS", 120L),
  n_parents = env_int("SCREEN_N_PARENTS", 80L),
  n_chr = env_int("SCREEN_N_CHR", 5L),
  seg_sites = env_int("SCREEN_SEG_SITES", 1200L),
  qtl_per_chr = env_int("SCREEN_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("SCREEN_SNP_PER_CHR", 1000L),
  progeny_per_cross = env_int("SCREEN_PROGENY_PER_CROSS", 60L),
  top_k = env_int("SCREEN_TOP_K", 10L),
  ridge_lambda = env_num("SCREEN_RIDGE_LAMBDA", 1),
  ld_window_size = env_int("SCREEN_LD_WINDOW_SIZE", 250L),
  screen_multipliers = as.integer(strsplit(Sys.getenv("SCREEN_MULTIPLIERS", "20,50,100"), ",", fixed = TRUE)[[1]])
)

dir.create("results", showWarnings = FALSE)
message("Dense LD screening comparison config:")
print(cfg)

methods <- c("pmv_nopos", "pmv_advanced")
screen_values <- c(NA_integer_, cfg$screen_multipliers)

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
  effects <- estimate_marker_effects_ridge(geno, parent_gv, cfg$ridge_lambda)

  timings <- list()
  timings[["ld_cache"]] <- system.time({
    ld_cache <- calculate_ld_matrix_efficient(
      geno_mat = geno,
      method = "cor",
      window_size = cfg$ld_window_size,
      block_size = 2000,
      use_parallel = FALSE,
      chr = map$chr,
      min_maf = 0,
      verbose = FALSE
    )
  })[["elapsed"]]

  pred_list <- list()
  selected_list <- list()
  for (method in methods) {
    for (screen_n in screen_values) {
      label <- if (is.na(screen_n)) {
        paste0(method, "_full")
      } else {
        paste0(method, "_screen", screen_n)
      }
      dense_n <- if (is.na(screen_n)) NULL else cfg$top_k * screen_n
      timings[[label]] <- system.time({
        pred <- build_cross_data(
          pheno = parent_gv,
          geno_mat = geno,
          candidate_ids = ids,
          variance_method = method,
          marker_effects_mcmc = effects,
          map = NULL,
          ld_mat = ld_cache,
          include_self = FALSE,
          use_parallel = FALSE,
          progeny = "DH",
          tau = 0.10,
          h2 = 1.0,
          trait_direction = "increase",
          mean_source = "auto",
          effect_training_n = nrow(geno),
          mean_underpowered_ratio = 1,
          dense_ld_screen_n = dense_n,
          dense_ld_screen_by = "pmv_simple_usefulness",
          include_effect_uncertainty = FALSE
        )
      })[["elapsed"]]
      pred$config <- label
      pred$method <- method
      pred$screen_multiplier <- if (is.na(screen_n)) NA_integer_ else screen_n
      pred_list[[label]] <- pred
      selected_list[[label]] <- transform(select_top(pred, cfg$top_k),
                                          config = label,
                                          method = method,
                                          screen_multiplier = if (is.na(screen_n)) NA_integer_ else screen_n)
    }
  }

  selected <- do.call(rbind, selected_list)
  full_keys_by_method <- lapply(methods, function(m) {
    pair_key(selected$parent1[selected$config == paste0(m, "_full")],
             selected$parent2[selected$config == paste0(m, "_full")])
  })
  names(full_keys_by_method) <- methods

  evaluated <- evaluate_unique_pairs(parent_pop, selected, parent_gv, cfg$progeny_per_cross, SP)
  summary <- summarize_selected(selected, evaluated, full_keys_by_method)
  summary$rep <- rep_id

  list(
    summary = summary,
    selected = transform(selected, rep = rep_id),
    timing = data.frame(rep = rep_id, method = names(timings),
                        elapsed_sec = as.numeric(timings), stringsAsFactors = FALSE)
  )
})

summary_by_rep <- do.call(rbind, lapply(rep_results, `[[`, "summary"))
selected <- do.call(rbind, lapply(rep_results, `[[`, "selected"))
timing <- do.call(rbind, lapply(rep_results, `[[`, "timing"))

overall <- do.call(rbind, lapply(unique(summary_by_rep$config), function(cfg_name) {
  x <- summary_by_rep[summary_by_rep$config == cfg_name, , drop = FALSE]
  data.frame(
    config = cfg_name,
    method = x$method[1],
    screen_multiplier = x$screen_multiplier[1],
    reps = length(unique(x$rep)),
    mean_overlap_with_full = mean(x$overlap_with_full, na.rm = TRUE),
    mean_realized_var = mean(x$mean_realized_var),
    mean_realized_best = mean(x$mean_realized_best),
    mean_top10 = mean(x$mean_top10),
    mean_transgressive_rate = mean(x$mean_transgressive_rate),
    mean_predicted_usefulness = mean(x$mean_predicted_usefulness),
    stringsAsFactors = FALSE
  )
}))

write.csv(overall, file.path("results", "dense_ld_screening_overall.csv"), row.names = FALSE)
write.csv(summary_by_rep, file.path("results", "dense_ld_screening_by_rep.csv"), row.names = FALSE)
write.csv(timing, file.path("results", "dense_ld_screening_timing.csv"), row.names = FALSE)
write.csv(selected, file.path("results", "dense_ld_screening_selected.csv"), row.names = FALSE)

message("Overall dense-LD screening comparison:")
print(overall)
message("Timing:")
print(stats::aggregate(elapsed_sec ~ method, timing, mean))
