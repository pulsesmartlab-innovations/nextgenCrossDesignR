# Diagnose why var_simple can appear better than pmv_simple in SNP-chip studies.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
suppressPackageStartupMessages(library(AlphaSimR))
source(file.path("R", "load.R"))

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
      predicted_var = cross_df$predicted_var[i],
      predicted_mean = cross_df$predicted_mean[i],
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

set.seed(20260425)
cfg <- list(
  n_founders = as.integer(Sys.getenv("DIAG_N_FOUNDERS", "120")),
  n_parents = as.integer(Sys.getenv("DIAG_N_PARENTS", "80")),
  n_chr = as.integer(Sys.getenv("DIAG_N_CHR", "5")),
  seg_sites = as.integer(Sys.getenv("DIAG_SEG_SITES", "1200")),
  qtl_per_chr = as.integer(Sys.getenv("DIAG_QTL_PER_CHR", "40")),
  snp_per_chr = as.integer(Sys.getenv("DIAG_SNP_PER_CHR", "1000")),
  top_k = as.integer(Sys.getenv("DIAG_TOP_K", "20")),
  random_k = as.integer(Sys.getenv("DIAG_RANDOM_K", "40")),
  progeny_per_cross = as.integer(Sys.getenv("DIAG_PROGENY_PER_CROSS", "150")),
  ridge_lambda = as.numeric(Sys.getenv("DIAG_RIDGE_LAMBDA", "1")),
  gen_len = 1,
  trait_var = 1
)

founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = cfg$gen_len)
SP <- SimParam$new(founder)
SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = cfg$trait_var)
SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
pop_all <- newPop(founder, simParam = SP)
parent_pop <- pop_all[seq_len(cfg$n_parents)]
parent_ids <- paste0("P", seq_len(nInd(parent_pop)))
parent_pop@id <- parent_ids

parent_gv <- as.numeric(gv(parent_pop)[, 1])
names(parent_gv) <- parent_ids

snp_geno <- pullSnpGeno(parent_pop, simParam = SP)
snp_map <- getSnpMap(simParam = SP)
colnames(snp_geno) <- snp_map$id
rownames(snp_geno) <- parent_ids
snp_map <- data.frame(marker = snp_map$id, chr = snp_map$chr, pos = snp_map$pos)

qtl_geno <- pullQtlGeno(parent_pop, simParam = SP)
qtl_map <- getQtlMap(simParam = SP)
colnames(qtl_geno) <- qtl_map$id
rownames(qtl_geno) <- parent_ids
qtl_map <- data.frame(marker = qtl_map$id, chr = qtl_map$chr, pos = qtl_map$pos)

snp_effects <- estimate_marker_effects_ridge(snp_geno, parent_gv, lambda = cfg$ridge_lambda)
qtl_effects <- matrix(SP$traits[[1]]@addEff, nrow = 1)
colnames(qtl_effects) <- qtl_map$marker

var_simple <- build_cross_data(parent_gv, snp_geno, parent_ids, "var_simple",
                               include_self = FALSE, use_parallel = FALSE)
pmv_snp_simple <- build_cross_data(parent_gv, snp_geno, parent_ids, "pmv_simple",
                                   marker_effects_mcmc = snp_effects,
                                   include_self = FALSE, use_parallel = FALSE,
                                   progeny = "DH")
pmv_qtl_simple <- build_cross_data(parent_gv, qtl_geno, parent_ids, "pmv_simple",
                                   marker_effects_mcmc = qtl_effects,
                                   include_self = FALSE, use_parallel = FALSE,
                                   progeny = "DH")
pmv_qtl_pos <- build_cross_data(parent_gv, qtl_geno, parent_ids, "pmv_pos",
                                marker_effects_mcmc = qtl_effects,
                                map = qtl_map,
                                include_self = FALSE, use_parallel = FALSE,
                                progeny = "DH",
                                pos_unit = "M",
                                window_cM = 100)

names(var_simple)[names(var_simple) == "cross_var"] <- "var_simple"
names(pmv_snp_simple)[names(pmv_snp_simple) == "cross_var"] <- "pmv_snp_simple"
names(pmv_qtl_simple)[names(pmv_qtl_simple) == "cross_var"] <- "pmv_qtl_simple"
names(pmv_qtl_pos)[names(pmv_qtl_pos) == "cross_var"] <- "pmv_qtl_pos"

pred <- Reduce(function(x, y) merge(x, y, by = c("parent1", "parent2")),
               list(var_simple[, c("parent1", "parent2", "cross_mean", "var_simple")],
                    pmv_snp_simple[, c("parent1", "parent2", "pmv_snp_simple")],
                    pmv_qtl_simple[, c("parent1", "parent2", "pmv_qtl_simple")],
                    pmv_qtl_pos[, c("parent1", "parent2", "pmv_qtl_pos")]))
names(pred)[names(pred) == "cross_mean"] <- "parent_mean"

pred$uc_var_simple <- pred$parent_mean + 1.755 * sqrt(pmax(pred$var_simple, 0))
pred$uc_pmv_snp_simple <- pred$parent_mean + 1.755 * sqrt(pmax(pred$pmv_snp_simple, 0))
pred$uc_pmv_qtl_simple <- pred$parent_mean + 1.755 * sqrt(pmax(pred$pmv_qtl_simple, 0))
pred$uc_pmv_qtl_pos <- pred$parent_mean + 1.755 * sqrt(pmax(pred$pmv_qtl_pos, 0))

metric_names <- c("var_simple", "pmv_snp_simple", "pmv_qtl_simple", "pmv_qtl_pos",
                  "uc_var_simple", "uc_pmv_snp_simple", "uc_pmv_qtl_simple", "uc_pmv_qtl_pos")

selected <- do.call(rbind, lapply(metric_names, function(m) {
  x <- pred[order(pred[[m]], decreasing = TRUE), ]
  x <- utils::head(x, cfg$top_k)
  data.frame(parent1 = x$parent1, parent2 = x$parent2, method = m,
             predicted_var = if (startsWith(m, "uc_")) NA_real_ else x[[m]],
             predicted_mean = x$parent_mean)
}))
set.seed(20260426)
random_rows <- pred[sample.int(nrow(pred), cfg$random_k), ]
selected <- rbind(selected, data.frame(parent1 = random_rows$parent1,
                                       parent2 = random_rows$parent2,
                                       method = "random",
                                       predicted_var = NA_real_,
                                       predicted_mean = random_rows$parent_mean))
selected$key <- pair_key(selected$parent1, selected$parent2)
eval_crosses <- selected[!duplicated(selected$key), c("parent1", "parent2", "method",
                                                      "predicted_var", "predicted_mean")]
evaluated <- evaluate_crosses(parent_pop, eval_crosses, parent_gv, cfg$progeny_per_cross, SP)
evaluated$key <- pair_key(evaluated$parent1, evaluated$parent2)

calibration <- merge(unique(evaluated[, c("parent1", "parent2", "realized_var",
                                          "realized_best", "top10_mean")]),
                     pred, by = c("parent1", "parent2"))
calibration_summary <- do.call(rbind, lapply(metric_names, function(m) {
  data.frame(
    metric = m,
    cor_realized_var = suppressWarnings(cor(calibration[[m]], calibration$realized_var)),
    cor_realized_best = suppressWarnings(cor(calibration[[m]], calibration$realized_best)),
    cor_top10_mean = suppressWarnings(cor(calibration[[m]], calibration$top10_mean))
  )
}))

selection_summary <- do.call(rbind, lapply(unique(selected$method), function(m) {
  keys <- selected$key[selected$method == m]
  x <- evaluated[evaluated$key %in% keys, ]
  data.frame(
    method = m,
    n_unique_crosses = nrow(x),
    mean_parent_mean = mean(x$predicted_mean, na.rm = TRUE),
    mean_realized_var = mean(x$realized_var),
    mean_realized_best = mean(x$realized_best),
    mean_top10 = mean(x$top10_mean),
    mean_transgressive_rate = mean(x$better_than_best_parent)
  )
}))

effect_summary <- data.frame(
  n_snp = ncol(snp_geno),
  n_qtl = ncol(qtl_geno),
  snp_effect_sd = sd(as.numeric(snp_effects)),
  qtl_effect_sd = sd(as.numeric(qtl_effects)),
  snp_effect_abs_top1pct_share = {
    a <- sort(abs(as.numeric(snp_effects)), decreasing = TRUE)
    sum(utils::head(a, max(1, ceiling(0.01 * length(a))))) / sum(a)
  },
  qtl_effect_abs_top1pct_share = {
    a <- sort(abs(as.numeric(qtl_effects)), decreasing = TRUE)
    sum(utils::head(a, max(1, ceiling(0.01 * length(a))))) / sum(a)
  }
)

dir.create("results", showWarnings = FALSE)
write.csv(selection_summary, file.path("results", "pmv_simple_diagnostic_selection_summary.csv"), row.names = FALSE)
write.csv(calibration_summary, file.path("results", "pmv_simple_diagnostic_calibration.csv"), row.names = FALSE)
write.csv(effect_summary, file.path("results", "pmv_simple_diagnostic_effects.csv"), row.names = FALSE)

message("Selection summary:")
print(selection_summary)
message("Calibration summary:")
print(calibration_summary)
message("Effect summary:")
print(effect_summary)
