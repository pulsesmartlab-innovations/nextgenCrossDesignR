# Multi-cycle AlphaSimR benchmark for corrected PMV with 80 parents and 5K SNPs.
#
# Methods:
# - var_simple: exhaustive relationship-distance baseline.
# - pmv_pos_full: corrected positional PMV screen, then exact posterior-sample PMV.
# - gms_pmv: same screen, then genomicMateSelectR-style cross-specific PMV.
#
# Note: genomicMateSelectR's public cross-variance functions target phased
# full-sib/gametic LD. They are useful here as an external audit of LD handling,
# but they are not a drop-in DH-from-inbred benchmark for AlphaSimR makeDH.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
suppressPackageStartupMessages(library(AlphaSimR))
suppressPackageStartupMessages(library(genomicMateSelectR))
source(file.path("R", "load.R"))

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.integer(value)
}

env_count <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  if (tolower(value) %in% c("all", "inf", "infinite")) return(Inf)
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

as_gms_haplo <- function(pop, sim_param, marker_names) {
  h <- pullSnpHaplo(pop, simParam = sim_param)
  colnames(h) <- marker_names
  rownames(h) <- as.vector(rbind(paste0(pop@id, "_HapA"), paste0(pop@id, "_HapB")))
  h
}

external_gms_pmv <- function(crosses, haplo_mat, recomb_freq_mat, effects_draws,
                             parent_values, tau = 0.10, h2 = 1.0) {
  effects_draws <- as.matrix(effects_draws)
  effects_centered <- scale(effects_draws, center = TRUE, scale = FALSE)
  beta_mu <- attr(effects_centered, "scaled:center")
  beta_cov <- crossprod(effects_centered) / max(1, nrow(effects_centered) - 1)
  names(beta_mu) <- colnames(effects_draws)
  rownames(beta_cov) <- colnames(beta_cov) <- colnames(effects_draws)

  parent_values <- setNames(as.numeric(parent_values), names(parent_values))
  intensity <- selection_intensity(tau)
  out <- vector("list", nrow(crosses))
  for (i in seq_len(nrow(crosses))) {
    sire <- as.character(crosses$parent1[i])
    dam <- as.character(crosses$parent2[i])
    x <- colSums(rbind(
      haplo_mat[grep(paste0("^", sire, "_"), rownames(haplo_mat)), , drop = FALSE],
      haplo_mat[grep(paste0("^", dam, "_"), rownames(haplo_mat)), , drop = FALSE]
    ))
    seg <- names(x[x > 0 & x < 4])
    if (!length(seg)) {
      cross_var <- 0
    } else {
      D <- genomicMateSelectR::calcCrossLD(
        sireID = sire,
        damID = dam,
        recombFreqMat = recomb_freq_mat[seg, seg, drop = FALSE],
        haploMat = haplo_mat[, seg, drop = FALSE]
      )
      vpm <- genomicMateSelectR::quadform(D = D, x = beta_mu[seg], y = beta_mu[seg])
      pmv_add <- sum(diag(D %*% beta_cov[seg, seg, drop = FALSE]))
      cross_var <- as.numeric(vpm + pmv_add)
    }
    if (!is.finite(cross_var) || cross_var < 0) cross_var <- 0
    cross_mean <- 0.5 * (parent_values[sire] + parent_values[dam])
    cross_sigma <- sqrt(cross_var)
    out[[i]] <- data.frame(
      parent1 = sire,
      parent2 = dam,
      cross_mean = cross_mean,
      cross_var = cross_var,
      cross_sigma = cross_sigma,
      cross_usefulness = cross_mean + intensity * sqrt(h2) * cross_sigma,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

select_crosses <- function(cross_df, top_crosses) {
  cross_df <- cross_df[order(cross_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
  utils::head(cross_df, top_crosses)
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

make_selected_progeny <- function(parent_pop, selected, progeny_per_cross, sim_param) {
  parent_ids <- parent_pop@id
  progeny <- vector("list", nrow(selected))
  family_metrics <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- match(selected$parent1[i], parent_ids)
    p2 <- match(selected$parent2[i], parent_ids)
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

cfg <- list(
  seed = env_int("PMV80_SEED", 20260426L),
  reps = env_int("PMV80_REPS", 3L),
  cycles = env_int("PMV80_CYCLES", 5L),
  n_founders = env_int("PMV80_N_FOUNDERS", 120L),
  n_parents = env_int("PMV80_N_PARENTS", 80L),
  n_chr = env_int("PMV80_N_CHR", 5L),
  seg_sites = env_int("PMV80_SEG_SITES", 1200L),
  qtl_per_chr = env_int("PMV80_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("PMV80_SNP_PER_CHR", 1000L),
  top_crosses = env_int("PMV80_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("PMV80_PROGENY_PER_CROSS", 40L),
  screen_n = env_count("PMV80_SCREEN_N", 500L),
  effect_method = env_chr("PMV80_EFFECT_METHOD", "ridge_posterior"),
  effect_training_n = env_int("PMV80_EFFECT_TRAINING_N", env_int("PMV80_N_PARENTS", 80L)),
  effect_draws = env_int("PMV80_EFFECT_DRAWS", 20L),
  gms_eval_n = env_int("PMV80_GMS_EVAL_N", 3L),
  ridge_lambda = env_num("PMV80_RIDGE_LAMBDA", 1),
  prior_genetic_var_scale = env_num("PMV80_PRIOR_GENETIC_VAR_SCALE", 1),
  bglr_model = env_chr("PMV80_BGLR_MODEL", "BRR"),
  bglr_nIter = env_int("PMV80_BGLR_NITER", NA_integer_),
  bglr_burnIn = env_int("PMV80_BGLR_BURNIN", 500L),
  bglr_thin = env_int("PMV80_BGLR_THIN", 5L),
  tau = env_num("PMV80_TAU", 0.10),
  h2 = env_num("PMV80_H2", 1.0),
  methods = strsplit(env_chr("PMV80_METHODS", "var_simple,pmv_pos_full"), ",", fixed = TRUE)[[1]]
)
cfg$methods <- trimws(cfg$methods)

dir.create("results", showWarnings = FALSE)
message("PMV 80-parent benchmark config:")
print(cfg)
if (is.finite(cfg$screen_n) && cfg$screen_n <= cfg$top_crosses &&
    any(cfg$methods %in% c("pmv_pos_full", "gms_pmv"))) {
  warning("PMV80_SCREEN_N is <= PMV80_TOP_CROSSES; exact PMV cannot materially change the screened set.")
}

methods <- cfg$methods

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

  branches <- setNames(vector("list", length(methods)), methods)
  for (m in methods) branches[[m]] <- base_parents

  metrics <- list()
  selections <- list()
  timings <- list()
  calibration <- list()
  family_diagnostics <- list()
  effect_diagnostics <- list()

  for (method in methods) {
    metrics[[paste(method, 0, sep = "_")]] <- population_metrics(
      branches[[method]], base_best, cycle = 0L, rep_id = rep_id, method = method
    )
  }

  snp_map0 <- getSnpMap(simParam = SP)
  marker_names <- paste(snp_map0$chr, snp_map0$id, sep = "_")
  recomb_pos <- snp_map0$pos * 100
  names(recomb_pos) <- marker_names
  recomb_freq_mat <- genomicMateSelectR::genmap2recombfreq(recomb_pos, nChr = cfg$n_chr)
  gms_linkage_mat <- 1 - 2 * recomb_freq_mat
  diag(gms_linkage_mat) <- 1

  for (cycle in seq_len(cfg$cycles)) {
    message("  cycle ", cycle, " of ", cfg$cycles)
    for (method in methods) {
      parent_pop <- branches[[method]]
      ids <- paste0(method, "_C", cycle - 1L, "_P", seq_len(nInd(parent_pop)))
      parent_pop@id <- ids
      parent_gv <- as.numeric(gv(parent_pop)[, 1])
      names(parent_gv) <- ids

      geno <- pullSnpGeno(parent_pop, simParam = SP)
      colnames(geno) <- marker_names
      rownames(geno) <- ids
      haplo <- as_gms_haplo(parent_pop, SP, marker_names)
      map <- data.frame(marker = marker_names, chr = snp_map0$chr, pos = snp_map0$pos,
                        stringsAsFactors = FALSE)

      if (method == "var_simple") {
        elapsed <- system.time({
          cross_df <- build_cross_data(
            pheno = parent_gv,
            geno_mat = geno,
            candidate_ids = ids,
            variance_method = "var_simple",
            include_self = FALSE,
            use_parallel = FALSE,
            tau = cfg$tau,
            h2 = cfg$h2,
            trait_direction = "increase"
          )
        })[["elapsed"]]
      } else {
        elapsed_effects <- system.time({
          training_pop <- make_training_population(parent_pop, cfg$effect_training_n, SP)
          training_ids <- paste0(method, "_C", cycle - 1L, "_T", seq_len(nInd(training_pop)))
          training_pop@id <- training_ids
          training_gv <- as.numeric(gv(training_pop)[, 1])
          training_geno <- pullSnpGeno(training_pop, simParam = SP)
          colnames(training_geno) <- marker_names
          rownames(training_geno) <- training_ids

          effects <- estimate_marker_effect_draws(
            geno_mat = training_geno,
            y = training_gv,
            method = cfg$effect_method,
            n_draws = cfg$effect_draws,
            lambda = cfg$ridge_lambda,
            prior_genetic_var = stats::var(training_gv) * cfg$prior_genetic_var_scale,
            bglr_model = cfg$bglr_model,
            bglr_nIter = if (is.na(cfg$bglr_nIter)) NULL else cfg$bglr_nIter,
            bglr_burnIn = cfg$bglr_burnIn,
            bglr_thin = cfg$bglr_thin
          )
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "effects",
          elapsed_sec = elapsed_effects, stringsAsFactors = FALSE
        )
        effect_diag <- effect_scale_diagnostics(effects, geno, parent_gv)
        effect_diag$rep <- rep_id
        effect_diag$method <- method
        effect_diag$cycle <- cycle
        effect_diagnostics[[paste(rep_id, method, cycle, sep = "_")]] <- effect_diag

        elapsed_screen <- system.time({
          pmv_screen <- build_cross_data(
            pheno = parent_gv,
            geno_mat = geno,
            candidate_ids = ids,
            variance_method = "pmv_pos",
            marker_effects_mcmc = effects,
            map = map,
            include_self = FALSE,
            use_parallel = FALSE,
            progeny = "DH",
            map_function = "haldane",
            pos_unit = "M",
            mean_source = "parent_value",
            tau = cfg$tau,
            h2 = cfg$h2,
            trait_direction = "increase",
            include_effect_uncertainty = TRUE
          )
          pmv_screen <- pmv_screen[order(pmv_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]

          diversity_screen <- build_cross_data(
            pheno = parent_gv,
            geno_mat = geno,
            candidate_ids = ids,
            variance_method = "var_simple",
            include_self = FALSE,
            use_parallel = FALSE,
            tau = cfg$tau,
            h2 = cfg$h2,
            trait_direction = "increase"
          )
          diversity_screen <- diversity_screen[order(diversity_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]

          n_screen <- min(nrow(pmv_screen), max(cfg$top_crosses, cfg$screen_n))
          screen_df <- unique(rbind(
            utils::head(pmv_screen[, c("parent1", "parent2"), drop = FALSE], n_screen),
            utils::head(diversity_screen[, c("parent1", "parent2"), drop = FALSE], n_screen)
          ))
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "screen",
          elapsed_sec = elapsed_screen, stringsAsFactors = FALSE
        )

        if (method == "pmv_pos_full") {
          elapsed <- system.time({
            cross_df <- build_cross_data(
              pheno = parent_gv,
              geno_mat = geno,
              candidate_ids = ids,
              variance_method = "pmv_pos_full",
              marker_effects_mcmc = effects,
              map = map,
              haplo_mat = haplo,
              include_self = FALSE,
              pair_subset = screen_df[, c("parent1", "parent2"), drop = FALSE],
              use_parallel = FALSE,
              progeny = "DH",
              map_function = "haldane",
              pos_unit = "M",
              mean_source = "parent_value",
              tau = cfg$tau,
              h2 = cfg$h2,
              trait_direction = "increase",
              include_effect_uncertainty = TRUE
            )
          })[["elapsed"]]
        } else {
          elapsed <- system.time({
            cross_df <- external_gms_pmv(
              crosses = screen_df[, c("parent1", "parent2"), drop = FALSE],
              haplo_mat = haplo,
              recomb_freq_mat = gms_linkage_mat,
              effects_draws = effects,
              parent_values = parent_gv,
              tau = cfg$tau,
              h2 = cfg$h2
            )
          })[["elapsed"]]
        }

      }

      timings[[length(timings) + 1L]] <- data.frame(
        rep = rep_id, method = method, cycle = cycle, stage = "score",
        elapsed_sec = elapsed, stringsAsFactors = FALSE
      )

      selected <- select_crosses(cross_df, cfg$top_crosses)
      if (method == "pmv_pos_full" && cfg$gms_eval_n > 0) {
        audit_pairs <- utils::head(selected[, c("parent1", "parent2"), drop = FALSE],
                                   min(nrow(selected), cfg$gms_eval_n))
        elapsed_gms <- system.time({
          gms_cmp <- external_gms_pmv(
            crosses = audit_pairs,
            haplo_mat = haplo,
            recomb_freq_mat = gms_linkage_mat,
            effects_draws = effects,
            parent_values = parent_gv,
            tau = cfg$tau,
            h2 = cfg$h2
          )
        })[["elapsed"]]
        pmv_cmp <- merge(
          selected[, c("parent1", "parent2", "cross_var"), drop = FALSE],
          gms_cmp[, c("parent1", "parent2", "cross_var"), drop = FALSE],
          by = c("parent1", "parent2"),
          suffixes = c("_ours", "_gms")
        )
        calibration[[paste(rep_id, method, cycle, sep = "_")]] <- data.frame(
          rep = rep_id,
          cycle = cycle,
          n_pairs = nrow(pmv_cmp),
          cor_var = suppressWarnings(stats::cor(pmv_cmp$cross_var_ours, pmv_cmp$cross_var_gms)),
          mean_abs_diff = mean(abs(pmv_cmp$cross_var_ours - pmv_cmp$cross_var_gms)),
          gms_elapsed_sec = elapsed_gms,
          stringsAsFactors = FALSE
        )
      }
      selected$rep <- rep_id
      selected$method <- method
      selected$cycle <- cycle
      selections[[paste(rep_id, method, cycle, sep = "_")]] <- selected

      progeny_result <- make_selected_progeny(parent_pop, selected, cfg$progeny_per_cross, SP)
      progeny <- progeny_result$pop
      family_diag <- progeny_result$family_metrics
      family_diag$rep <- rep_id
      family_diag$method <- method
      family_diag$cycle <- cycle
      family_diagnostics[[paste(rep_id, method, cycle, sep = "_")]] <- family_diag

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

  list(
    metrics = bind_rows_fill(metrics),
    selections = bind_rows_fill(selections),
    timings = bind_rows_fill(timings),
    calibration = bind_rows_fill(calibration),
    family_diagnostics = bind_rows_fill(family_diagnostics),
    effect_diagnostics = bind_rows_fill(effect_diagnostics)
  )
})

metrics <- bind_rows_fill(lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
timings <- bind_rows_fill(lapply(rep_results, `[[`, "timings"))
calibration <- bind_rows_fill(lapply(rep_results, `[[`, "calibration"))
family_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "family_diagnostics"))
effect_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "effect_diagnostics"))

overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ method + cycle,
  metrics,
  mean
)
timing_summary <- stats::aggregate(elapsed_sec ~ method + stage, timings, mean)

write.csv(metrics, file.path("results", "pmv80_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, file.path("results", "pmv80_metrics_overall.csv"), row.names = FALSE)
write.csv(selections, file.path("results", "pmv80_selected_crosses.csv"), row.names = FALSE)
write.csv(timings, file.path("results", "pmv80_timing.csv"), row.names = FALSE)
write.csv(timing_summary, file.path("results", "pmv80_timing_summary.csv"), row.names = FALSE)
write.csv(calibration, file.path("results", "pmv80_gms_calibration.csv"), row.names = FALSE)
write.csv(family_diagnostics, file.path("results", "pmv80_family_diagnostics.csv"), row.names = FALSE)
write.csv(effect_diagnostics, file.path("results", "pmv80_effect_diagnostics.csv"), row.names = FALSE)

message("PMV 80-parent benchmark overall:")
print(overall)
message("Timing summary:")
print(timing_summary)
if (nrow(calibration)) {
  message("Our full PMV vs genomicMateSelectR-style PMV calibration:")
  calib_summary <- stats::aggregate(mean_abs_diff ~ cycle, calibration, mean, na.rm = TRUE)
  if (any(is.finite(calibration$cor_var))) {
    cor_summary <- stats::aggregate(cor_var ~ cycle, calibration, mean, na.rm = TRUE)
    calib_summary <- merge(calib_summary, cor_summary, by = "cycle", all = TRUE)
  }
  print(calib_summary)
}
if (nrow(family_diagnostics)) {
  message("Predicted variance vs realized family diagnostics:")
  fam_summary <- stats::aggregate(
    cbind(pred_var, realized_var_gv, realized_max_gv, realized_top10_gv) ~ method + cycle,
    family_diagnostics,
    mean
  )
  print(fam_summary)
}
if (nrow(effect_diagnostics)) {
  message("Marker-effect scale diagnostics:")
  eff_summary <- stats::aggregate(
    cbind(effect_sd, parent_value_sd, marker_score_sd, marker_score_cor, marker_score_slope) ~ method + cycle,
    effect_diagnostics,
    mean,
    na.rm = TRUE
  )
  print(eff_summary)
}
