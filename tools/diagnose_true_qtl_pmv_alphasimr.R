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

as_qtl_haplo <- function(pop, sim_param, marker_names) {
  haplo <- pullQtlHaplo(pop, simParam = sim_param)
  colnames(haplo) <- marker_names
  rownames(haplo) <- as.vector(rbind(paste0(pop@id, "_HapA"), paste0(pop@id, "_HapB")))
  haplo
}

true_qtl_effects <- function(sim_param, marker_names) {
  eff <- as.numeric(sim_param$traits[[1]]@addEff)
  if (length(eff) != length(marker_names)) {
    stop("True QTL effect count does not match QTL genotype columns.", call. = FALSE)
  }
  out <- matrix(eff, nrow = 1)
  colnames(out) <- marker_names
  out
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

select_crosses <- function(cross_df, top_crosses) {
  cross_df <- cross_df[order(cross_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  utils::head(cross_df, top_crosses)
}

cfg <- list(
  seed = env_int("TRUEQTL_SEED", 20260426L),
  reps = env_int("TRUEQTL_REPS", 3L),
  cycles = env_int("TRUEQTL_CYCLES", 5L),
  n_founders = env_int("TRUEQTL_N_FOUNDERS", 120L),
  n_parents = env_int("TRUEQTL_N_PARENTS", 80L),
  n_chr = env_int("TRUEQTL_N_CHR", 5L),
  seg_sites = env_int("TRUEQTL_SEG_SITES", 1200L),
  qtl_per_chr = env_int("TRUEQTL_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("TRUEQTL_SNP_PER_CHR", 1000L),
  top_crosses = env_int("TRUEQTL_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("TRUEQTL_PROGENY_PER_CROSS", 40L),
  tau = env_num("TRUEQTL_TAU", 0.10),
  h2 = env_num("TRUEQTL_H2", 1.0)
)

dir.create("results", showWarnings = FALSE)
print(cfg)

rep_results <- lapply(seq_len(cfg$reps), function(rep_id) {
  set.seed(cfg$seed + rep_id)
  founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
  SP <- SimParam$new(founder)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  base_pop <- newPop(founder, simParam = SP)
  base_parents <- base_pop[seq_len(cfg$n_parents)]
  base_parents@id <- paste0("P0_", seq_len(nInd(base_parents)))
  base_best <- max(as.numeric(gv(base_parents)[, 1]))

  branches <- list(var_simple = base_parents, pmv_true_qtl = base_parents)
  metrics <- list()
  selections <- list()
  timings <- list()

  for (method in names(branches)) {
    metrics[[paste(method, 0, sep = "_")]] <- population_metrics(
      branches[[method]], base_best, cycle = 0L, rep_id = rep_id, method = method
    )
  }

  for (cycle in seq_len(cfg$cycles)) {
    for (method in names(branches)) {
      parent_pop <- branches[[method]]
      ids <- paste0(method, "_C", cycle - 1L, "_P", seq_len(nInd(parent_pop)))
      parent_pop@id <- ids
      parent_gv <- as.numeric(gv(parent_pop)[, 1])
      names(parent_gv) <- ids

      elapsed <- system.time({
        if (method == "var_simple") {
          snp_geno <- pullSnpGeno(parent_pop, simParam = SP)
          snp_map <- getSnpMap(simParam = SP)
          colnames(snp_geno) <- snp_map$id
          rownames(snp_geno) <- ids
          cross_df <- build_cross_data(
            pheno = parent_gv,
            geno_mat = snp_geno,
            candidate_ids = ids,
            variance_method = "var_simple",
            include_self = FALSE,
            use_parallel = FALSE,
            tau = cfg$tau,
            h2 = cfg$h2,
            trait_direction = "increase"
          )
        } else {
          qtl_geno <- pullQtlGeno(parent_pop, simParam = SP)
          qtl_map_raw <- getQtlMap(trait = 1, simParam = SP)
          colnames(qtl_geno) <- qtl_map_raw$id
          rownames(qtl_geno) <- ids
          qtl_map <- data.frame(
            marker = qtl_map_raw$id,
            chr = qtl_map_raw$chr,
            pos = qtl_map_raw$pos,
            stringsAsFactors = FALSE
          )
          qtl_haplo <- as_qtl_haplo(parent_pop, SP, qtl_map$marker)
          qtl_effects <- true_qtl_effects(SP, colnames(qtl_geno))
          cross_df <- build_cross_data(
            pheno = parent_gv,
            geno_mat = qtl_geno,
            candidate_ids = ids,
            variance_method = "pmv_pos_full",
            marker_effects_mcmc = qtl_effects,
            map = qtl_map,
            haplo_mat = qtl_haplo,
            pair_subset = build_parent_pairs(ids, include_self = FALSE),
            include_self = FALSE,
            use_parallel = FALSE,
            progeny = "DH",
            map_function = "haldane",
            pos_unit = "M",
            mean_source = "parent_value",
            effect_training_n = Inf,
            include_effect_uncertainty = FALSE,
            tau = cfg$tau,
            h2 = cfg$h2,
            trait_direction = "increase"
          )
        }
      })[["elapsed"]]

      timings[[length(timings) + 1L]] <- data.frame(
        rep = rep_id, method = method, cycle = cycle, stage = "score",
        elapsed_sec = elapsed, stringsAsFactors = FALSE
      )
      selected <- select_crosses(cross_df, cfg$top_crosses)
      selected$rep <- rep_id
      selected$method <- method
      selected$cycle <- cycle
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
  list(metrics = do.call(rbind, metrics),
       selections = bind_rows_fill(selections),
       timings = bind_rows_fill(timings))
})

metrics <- do.call(rbind, lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
timings <- bind_rows_fill(lapply(rep_results, `[[`, "timings"))
overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ method + cycle,
  metrics,
  mean
)
timing_summary <- stats::aggregate(elapsed_sec ~ method + stage, timings, mean)

write.csv(metrics, file.path("results", "true_qtl_pmv_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, file.path("results", "true_qtl_pmv_metrics_overall.csv"), row.names = FALSE)
write.csv(selections, file.path("results", "true_qtl_pmv_selected_crosses.csv"), row.names = FALSE)
write.csv(timings, file.path("results", "true_qtl_pmv_timing.csv"), row.names = FALSE)
write.csv(timing_summary, file.path("results", "true_qtl_pmv_timing_summary.csv"), row.names = FALSE)

message("True-QTL PMV diagnostic overall:")
print(overall)
message("Timing summary:")
print(timing_summary)
