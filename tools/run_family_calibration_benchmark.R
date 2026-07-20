# Family-level metric calibration benchmark.
#
# This script separates statistical score quality from mate allocation. It
# samples the same candidate crosses once, scores those crosses with internal
# and external metrics, simulates AlphaSimR DH families, and compares predicted
# scores against realized family mean, variance, top-tail, and maximum value.

local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

root <- normalizePath(file.path(getwd(), "nextgen_cross_design"), mustWork = FALSE)
if (!dir.exists(root)) root <- normalizePath(file.path(".."), mustWork = TRUE)
source(file.path(root, "R", "load.R"))
ng_use_cpp <- tolower(trimws(Sys.getenv("NG_USE_CPP", unset = "0"))) %in% c("1", "true", "yes", "y")
ng_load(root, use_cpp = ng_use_cpp)

suppressPackageStartupMessages(library(AlphaSimR))

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

env_bool <- function(name, default) {
  value <- tolower(trimws(Sys.getenv(name, unset = NA_character_)))
  if (is.na(value) || !nzchar(value)) return(default)
  value %in% c("1", "true", "yes", "y")
}

env_int_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

attr_chr <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  paste(as.character(x), collapse = "; ")
}

with_rng_seed <- function(seed, expr) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  if (!is.null(seed)) set.seed(seed)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

make_dh_base_parents <- function(founder_pop, n_parents, sim_param) {
  n_dh_per_founder <- ceiling(n_parents / nInd(founder_pop))
  dh <- makeDH(founder_pop, nDH = n_dh_per_founder, keepParents = FALSE, simParam = sim_param)
  dh[seq_len(n_parents)]
}

make_training_population <- function(parent_pop, target_n, sim_param) {
  if (target_n <= nInd(parent_pop)) return(parent_pop)
  extra_n <- target_n - nInd(parent_pop)
  p1 <- sample.int(nInd(parent_pop), extra_n, replace = TRUE)
  p2 <- sample.int(nInd(parent_pop), extra_n, replace = TRUE)
  same <- p1 == p2
  while (any(same)) {
    p2[same] <- sample.int(nInd(parent_pop), sum(same), replace = TRUE)
    same <- p1 == p2
  }
  f1 <- makeCross(parent_pop, cbind(p1, p2), nProgeny = 1, simParam = sim_param)
  aux <- makeDH(f1, nDH = 1, keepParents = FALSE, simParam = sim_param)
  c(parent_pop, aux)
}

phenotype <- function(pop, h2 = 0.5, reps = 3, seed = NULL) {
  g <- as.numeric(gv(pop)[, 1])
  if (h2 >= 0.999) return(g)
  vg <- stats::var(g)
  ve <- if (is.finite(vg) && vg > 0) vg * (1 - h2) / h2 / max(1, reps) else 1
  with_rng_seed(seed, g + stats::rnorm(length(g), sd = sqrt(ve)))
}

simulate_family_realizations <- function(parent_pop,
                                         scores,
                                         n_progeny,
                                         top_prop,
                                         sim_param) {
  parent_ids <- parent_pop@id
  out <- vector("list", nrow(scores))
  for (i in seq_len(nrow(scores))) {
    p1 <- match(scores$parent1[i], parent_ids)
    p2 <- match(scores$parent2[i], parent_ids)
    if (is.na(p1) || is.na(p2)) stop("Selected parent ID missing from population.", call. = FALSE)
    f1 <- makeCross(parent_pop, matrix(c(p1, p2), ncol = 2), nProgeny = 1, simParam = sim_param)
    dh <- makeDH(f1, nDH = n_progeny, keepParents = FALSE, simParam = sim_param)
    g <- as.numeric(gv(dh)[, 1])
    top_n <- max(1L, ceiling(top_prop * length(g)))
    realized <- data.frame(
      family_index = i,
      n_progeny = length(g),
      realized_mean = mean(g),
      realized_var = stats::var(g),
      realized_top10 = mean(utils::head(sort(g, decreasing = TRUE), top_n)),
      realized_max = max(g),
      stringsAsFactors = FALSE
    )
    out[[i]] <- cbind(scores[i, , drop = FALSE], realized)
  }
  ng_bind_rows_fill(out)
}

cfg <- list(
  seed = env_int("NG_SEED", 9201L),
  reps = env_int("NG_CAL_REPS", env_int("NG_REPS", 5L)),
  parent_sizes = env_int_vec("NG_CAL_PARENT_SIZES", env_int("NG_N_PARENTS", 80L)),
  n_founders = env_int("NG_N_FOUNDERS", 120L),
  n_chr = env_int("NG_N_CHR", 5L),
  seg_sites = env_int("NG_SEG_SITES", 1200L),
  snp_per_chr = env_int("NG_SNP_PER_CHR", 1000L),
  qtl_per_chr = env_int("NG_QTL_PER_CHR", 40L),
  n_families = env_int("NG_CAL_N_FAMILIES", 80L),
  random_families = env_int("NG_CAL_RANDOM_FAMILIES", NA_integer_),
  top_per_metric = env_int("NG_CAL_TOP_PER_METRIC", 5L),
  bottom_per_metric = env_int("NG_CAL_BOTTOM_PER_METRIC", 2L),
  progeny_per_cross = env_int("NG_PROGENY_PER_CROSS", 80L),
  phenotype_h2 = env_num("NG_PHENO_H2", 0.5),
  phenotype_reps = env_int("NG_PHENO_REPS", 3L),
  effect_training_n = env_int("NG_EFFECT_TRAINING_N", NA_integer_),
  effect_h2_prior = env_num("NG_EFFECT_H2_PRIOR", 0.5),
  effect_kfold = env_int("NG_EFFECT_KFOLD", 5L),
  selection_prop = env_num("NG_SELECTION_PROP", 0.10),
  min_effect_reliability = env_num("NG_MIN_EFFECT_RELIABILITY", 0.35),
  topk_prop = env_num("NG_TOPK_PROP", 0.10),
  include_external = env_bool("NG_CAL_INCLUDE_EXTERNAL", TRUE),
  include_gms = env_bool("NG_CAL_INCLUDE_GMS", TRUE),
  gms_max_markers = env_int("NG_CAL_GMS_MAX_MARKERS", 5000L),
  gms_ncores = env_int("NG_CAL_GMS_NCORES", 1L),
  alphasimr_threads = env_int("NG_ALPHASIMR_THREADS", 1L),
  simplemating_threads = env_int("NG_SIMPLEMATING_THREADS", 1L),
  use_cpp = ng_use_cpp,
  output_prefix = env_chr("NG_OUTPUT_PREFIX", "family_calibration")
)
if (!is.finite(cfg$random_families)) cfg$random_families <- ceiling(cfg$n_families / 2)
if (!is.finite(cfg$alphasimr_threads) || cfg$alphasimr_threads < 1L) cfg$alphasimr_threads <- 1L

dir.create("results", showWarnings = FALSE)
write.csv(data.frame(name = names(cfg), value = vapply(cfg, function(x) paste(x, collapse = ","), character(1))),
          file.path("results", paste0(cfg$output_prefix, "_config.csv")), row.names = FALSE)

message("nextgen_cross_design family calibration config:")
print(cfg)

rep_results <- list()
run_idx <- 1L
for (n_parents in cfg$parent_sizes) {
  for (rep_id in seq_len(cfg$reps)) {
    message("Parent size ", n_parents, ", replicate ", rep_id, " of ", cfg$reps)
    set.seed(cfg$seed + n_parents * 1000L + rep_id)
    founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
    SP <- SimParam$new(founder)
    SP$nThreads <- max(1L, as.integer(cfg$alphasimr_threads))
    SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
    SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)

    base_pop <- newPop(founder, simParam = SP)
    parents <- make_dh_base_parents(base_pop, n_parents, SP)
    parent_ids <- paste0("R", rep_id, "_N", n_parents, "_P", seq_len(nInd(parents)))
    parents@id <- parent_ids

    snp_map <- getSnpMap(simParam = SP)
    marker_names <- paste0("Chr", snp_map$chr, "_", snp_map$id)
    marker_map <- data.frame(
      marker = marker_names,
      chr = snp_map$chr,
      pos_cm = snp_map$pos * 100,
      stringsAsFactors = FALSE
    )
    geno <- pullSnpGeno(parents, simParam = SP)
    colnames(geno) <- marker_names
    rownames(geno) <- parent_ids

    parent_y <- phenotype(
      parents,
      h2 = cfg$phenotype_h2,
      reps = cfg$phenotype_reps,
      seed = cfg$seed + n_parents * 100000L + rep_id * 100L + 11L
    )
    names(parent_y) <- parent_ids

    train_n <- if (is.finite(cfg$effect_training_n)) cfg$effect_training_n else n_parents
    t_effect <- system.time({
      training_pop <- with_rng_seed(
        cfg$seed + n_parents * 100000L + rep_id * 100L + 19L,
        make_training_population(parents, train_n, SP)
      )
      training_ids <- paste0("R", rep_id, "_N", n_parents, "_T", seq_len(nInd(training_pop)))
      training_pop@id <- training_ids
      training_geno <- pullSnpGeno(training_pop, simParam = SP)
      colnames(training_geno) <- marker_names
      rownames(training_geno) <- training_ids
      training_y <- phenotype(
        training_pop,
        h2 = cfg$phenotype_h2,
        reps = cfg$phenotype_reps,
        seed = cfg$seed + n_parents * 100000L + rep_id * 100L + 17L
      )
      names(training_y) <- training_ids
      effects <- ng_fit_ridge_effects(
        training_geno,
        training_y,
        ids = training_ids,
        h2_prior = cfg$effect_h2_prior,
        kfold = cfg$effect_kfold,
        seed = cfg$seed + rep_id + n_parents
      )
    })[["elapsed"]]

    t_screen <- system.time({
      cheap_scores <- ng_cheap_cross_screen(
        geno = geno,
        effects = effects,
        ids = parent_ids,
        adjusted_pheno = parent_y,
        selection_prop = cfg$selection_prop,
        min_effect_reliability = cfg$min_effect_reliability
      )
      pair_rows <- ng_select_calibration_pair_rows(
        cheap_scores,
        n_families = cfg$n_families,
        random_n = cfg$random_families,
        top_per_metric = cfg$top_per_metric,
        bottom_per_metric = cfg$bottom_per_metric,
        seed = cfg$seed + n_parents * 100000L + rep_id * 100L + 23L
      )
      sampled_pairs <- cheap_scores[pair_rows, c("parent1", "parent2"), drop = FALSE]
    })[["elapsed"]]

    t_score <- system.time({
      scores <- ng_score_crosses(
        geno = geno,
        effects = effects,
        marker_map = marker_map,
        ids = parent_ids,
        pairs = sampled_pairs,
        adjusted_pheno = parent_y,
        selection_prop = cfg$selection_prop,
        min_effect_reliability = cfg$min_effect_reliability,
        use_cpp = cfg$use_cpp
      )
      calibrators <- ng_fit_family_variance_calibrators(NULL)
      scores <- ng_apply_family_variance_calibrators(
        scores,
        calibrators = calibrators,
        n_progeny = cfg$progeny_per_cross,
        top_prop = cfg$topk_prop
      )
      if (isTRUE(cfg$include_external)) {
        row_index <- seq_len(nrow(scores))
        scores <- ng_add_popvar_scores(
          scores = scores,
          geno = geno,
          effects = effects,
          marker_map = marker_map,
          adjusted_pheno = parent_y,
          tail_p = cfg$selection_prop,
          row_index = row_index
        )
        scores <- ng_add_simplemating_scores(
          scores = scores,
          geno = geno,
          effects = effects,
          marker_map = marker_map,
          adjusted_pheno = parent_y,
          prop_sel = cfg$selection_prop,
          n_threads = cfg$simplemating_threads,
          row_index = row_index
        )
      }
      if (isTRUE(cfg$include_gms)) {
        scores <- ng_add_gms_vpm_scores(
          scores = scores,
          geno = geno,
          effects = effects,
          marker_map = marker_map,
          row_index = seq_len(nrow(scores)),
          max_markers = cfg$gms_max_markers,
          n_threads = cfg$gms_ncores
        )
      }
    })[["elapsed"]]

    t_family <- system.time({
      families <- simulate_family_realizations(
        parent_pop = parents,
        scores = scores,
        n_progeny = cfg$progeny_per_cross,
        top_prop = cfg$topk_prop,
        sim_param = SP
      )
    })[["elapsed"]]
    families$rep <- rep_id
    families$n_parents <- n_parents
    families$effect_training_n <- train_n
    families$effect_reliability <- effects$reliability
    families$in_sample_reliability <- effects$in_sample_reliability

    registry <- ng_family_metric_registry(families)
    summary <- ng_family_metric_summary(
      families,
      registry = registry,
      top_prop = cfg$topk_prop
    )
    summary$rep <- rep_id
    summary$n_parents <- n_parents
    summary$effect_training_n <- train_n
    summary$effect_reliability <- effects$reliability

    timing <- data.frame(
      rep = rep_id,
      n_parents = n_parents,
      stage = c("effects", "cheap_screen", "score", "simulate_families"),
      elapsed_sec = c(t_effect, t_screen, t_score, t_family),
      stringsAsFactors = FALSE
    )
    status <- data.frame(
      rep = rep_id,
      n_parents = n_parents,
      n_scores = nrow(scores),
      popvar_status = paste(unique(scores$popvar_status), collapse = ","),
      simple_status = paste(unique(scores$simple_status), collapse = ","),
      gms_status = paste(unique(scores$gms_status), collapse = ","),
      popvar_error = attr_chr(attr(scores, "popvar_error")),
      simple_mpv_error = attr_chr(attr(scores, "simple_mpv_error")),
      simple_usefa_error = attr_chr(attr(scores, "simple_usefa_error")),
      gms_error = attr_chr(attr(scores, "gms_error")),
      effect_reliability = effects$reliability,
      in_sample_reliability = effects$in_sample_reliability,
      stringsAsFactors = FALSE
    )

    rep_results[[run_idx]] <- list(
      families = families,
      summary = summary,
      timings = timing,
      status = status
    )
    run_idx <- run_idx + 1L
  }
}

families <- ng_bind_rows_fill(lapply(rep_results, `[[`, "families"))
summary_by_rep <- ng_bind_rows_fill(lapply(rep_results, `[[`, "summary"))
timings <- ng_bind_rows_fill(lapply(rep_results, `[[`, "timings"))
status <- ng_bind_rows_fill(lapply(rep_results, `[[`, "status"))

summary_avg <- if (nrow(summary_by_rep)) {
  stats::aggregate(
    cbind(pearson, spearman, slope, intercept, linear_rmse, raw_rmse,
          top_pred_mean, all_mean, top_pred_delta, top_overlap,
          effect_reliability) ~ n_parents + effect_training_n + score_col + metric_group + target + primary_target,
    summary_by_rep,
    mean,
    na.rm = TRUE,
    na.action = stats::na.pass
  )
} else data.frame()

timing_summary <- if (nrow(timings)) {
  stats::aggregate(elapsed_sec ~ n_parents + stage, timings, mean)
} else data.frame()

prefix <- file.path("results", cfg$output_prefix)
write.csv(families, paste0(prefix, "_families.csv"), row.names = FALSE)
write.csv(summary_by_rep, paste0(prefix, "_metric_summary_by_rep.csv"), row.names = FALSE)
write.csv(summary_avg, paste0(prefix, "_metric_summary_avg.csv"), row.names = FALSE)
write.csv(timings, paste0(prefix, "_timings.csv"), row.names = FALSE)
write.csv(timing_summary, paste0(prefix, "_timing_summary.csv"), row.names = FALSE)
write.csv(status, paste0(prefix, "_status.csv"), row.names = FALSE)

message("Status:")
print(status)
message("Timing summary:")
print(timing_summary)
message("Primary target summary, ordered by target and Spearman:")
primary <- summary_avg[summary_avg$primary_target %in% TRUE, , drop = FALSE]
if (nrow(primary)) {
  primary <- primary[order(primary$target, -primary$spearman, -primary$top_pred_delta), , drop = FALSE]
  print(primary[, c("n_parents", "score_col", "metric_group", "target", "pearson", "spearman",
                    "slope", "top_pred_delta", "top_overlap", "effect_reliability"), drop = FALSE])
}
