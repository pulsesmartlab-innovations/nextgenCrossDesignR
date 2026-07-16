local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

source("R/load.R")
ng_load(use_cpp = Sys.getenv("NG_USE_CPP", "0") != "0")

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  stop("AlphaSimR is required for this diagnostic.", call. = FALSE)
}

suppressPackageStartupMessages(library(AlphaSimR))

env_int <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.integer(x)
}

env_num <- function(name, default) {
  x <- Sys.getenv(name, unset = NA_character_)
  if (is.na(x) || !nzchar(x)) return(default)
  as.numeric(x)
}

make_dh_base_parents <- function(founder_pop, n_parents, sim_param) {
  n_dh_per_founder <- ceiling(n_parents / nInd(founder_pop))
  dh <- makeDH(founder_pop, nDH = n_dh_per_founder, keepParents = FALSE, simParam = sim_param)
  dh[seq_len(n_parents)]
}

make_training_population <- function(parent_pop, target_n, sim_param) {
  if (nInd(parent_pop) >= target_n) return(parent_pop)
  need <- target_n - nInd(parent_pop)
  p1 <- sample(seq_len(nInd(parent_pop)), need, replace = TRUE)
  p2 <- sample(seq_len(nInd(parent_pop)), need, replace = TRUE)
  same <- p1 == p2
  while (any(same)) {
    p2[same] <- sample(seq_len(nInd(parent_pop)), sum(same), replace = TRUE)
    same <- p1 == p2
  }
  f1 <- makeCross(parent_pop, cbind(p1, p2), nProgeny = 1, simParam = sim_param)
  aux <- makeDH(f1, nDH = 1, keepParents = FALSE, simParam = sim_param)
  c(parent_pop, aux)
}

with_rng_seed <- function(seed, expr) {
  if (is.null(seed) || !is.finite(seed)) return(force(expr))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  set.seed(seed)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

phenotype <- function(pop, h2 = 0.5, reps = 3, seed = NULL) {
  g <- as.numeric(gv(pop)[, 1])
  if (h2 >= 0.999) return(g)
  vg <- stats::var(g)
  ve <- if (is.finite(vg) && vg > 0) vg * (1 - h2) / h2 / max(1, reps) else 1
  with_rng_seed(seed, g + stats::rnorm(length(g), sd = sqrt(ve)))
}

top_set <- function(scores, col, n = 10L) {
  ok <- is.finite(scores[[col]])
  idx <- which(ok)[order(scores[[col]][ok], decreasing = TRUE)]
  ng_pair_key(scores$parent1[idx[seq_len(min(n, length(idx)))]],
              scores$parent2[idx[seq_len(min(n, length(idx)))]])
}

compare_pair <- function(scores, col_a, col_b, n_top = 10L) {
  ok <- is.finite(scores[[col_a]]) & is.finite(scores[[col_b]])
  if (sum(ok) < 3L) {
    return(data.frame(metric_a = col_a, metric_b = col_b, n = sum(ok),
                      max_abs_delta = NA_real_, sd_delta = NA_real_,
                      spearman = NA_real_, top_overlap = NA_real_))
  }
  a <- scores[[col_a]][ok]
  b <- scores[[col_b]][ok]
  top_a <- top_set(scores[ok, , drop = FALSE], col_a, n = n_top)
  top_b <- top_set(scores[ok, , drop = FALSE], col_b, n = n_top)
  data.frame(
    metric_a = col_a,
    metric_b = col_b,
    n = sum(ok),
    max_abs_delta = max(abs(a - b), na.rm = TRUE),
    sd_delta = stats::sd(a - b, na.rm = TRUE),
    spearman = suppressWarnings(stats::cor(a, b, method = "spearman")),
    top_overlap = length(intersect(top_a, top_b)) / max(1L, min(length(top_a), length(top_b))),
    stringsAsFactors = FALSE
  )
}

cfg <- list(
  seed = env_int("NG_DIAG_SEED", 7301L),
  n_founders = env_int("NG_DIAG_N_FOUNDERS", 60L),
  n_parents = env_int("NG_DIAG_N_PARENTS", 30L),
  n_chr = env_int("NG_DIAG_N_CHR", 3L),
  seg_sites = env_int("NG_DIAG_SEG_SITES", 300L),
  snp_per_chr = env_int("NG_DIAG_SNP_PER_CHR", 200L),
  qtl_per_chr = env_int("NG_DIAG_QTL_PER_CHR", 20L),
  effect_training_n = env_int("NG_DIAG_EFFECT_TRAINING_N", 120L),
  pheno_h2 = env_num("NG_DIAG_PHENO_H2", 0.5),
  pheno_reps = env_int("NG_DIAG_PHENO_REPS", 3L),
  effect_kfold = env_int("NG_DIAG_EFFECT_KFOLD", 3L),
  selection_prop = env_num("NG_DIAG_SELECTION_PROP", 0.10),
  top_n = env_int("NG_DIAG_TOP_N", 10L),
  output_prefix = Sys.getenv("NG_DIAG_OUTPUT_PREFIX", "results/usefulness_equivalence")
)

set.seed(cfg$seed)
founder <- quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = 1)
SP <- SimParam$new(founder)
SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)

base_pop <- newPop(founder, simParam = SP)
parents <- make_dh_base_parents(base_pop, cfg$n_parents, SP)
ids <- paste0("P", seq_len(nInd(parents)))
parents@id <- ids

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
rownames(geno) <- ids
parent_y <- phenotype(parents, h2 = cfg$pheno_h2, reps = cfg$pheno_reps, seed = cfg$seed + 2L)
names(parent_y) <- ids

training_pop <- make_training_population(parents, cfg$effect_training_n, SP)
training_ids <- paste0("T", seq_len(nInd(training_pop)))
training_pop@id <- training_ids
training_geno <- pullSnpGeno(training_pop, simParam = SP)
colnames(training_geno) <- marker_names
rownames(training_geno) <- training_ids
training_y <- phenotype(training_pop, h2 = cfg$pheno_h2, reps = cfg$pheno_reps, seed = cfg$seed + 3L)
names(training_y) <- training_ids

effects <- ng_fit_ridge_effects(
  training_geno,
  training_y,
  ids = training_ids,
  h2_prior = cfg$pheno_h2,
  kfold = cfg$effect_kfold,
  seed = cfg$seed + 1L
)

scores <- ng_score_crosses(
  geno = geno,
  effects = effects,
  marker_map = marker_map,
  ids = ids,
  adjusted_pheno = parent_y,
  selection_prop = cfg$selection_prop
)
scores <- ng_add_external_baseline_scores(
  scores = scores,
  geno = geno,
  effects = effects,
  marker_map = marker_map,
  adjusted_pheno = parent_y,
  methods = c("popvar_uc_topn", "simple_usefa_topn"),
  selection_prop = cfg$selection_prop,
  n_crosses = cfg$top_n
)

comparisons <- do.call(rbind, list(
  compare_pair(scores, "dh_recomb_var", "simple_usefa_var", cfg$top_n),
  compare_pair(scores, "uc_recomb_gebv", "simple_usefa", cfg$top_n),
  compare_pair(scores, "dh_recomb_var", "popvar_varG", cfg$top_n),
  compare_pair(scores, "uc_recomb_gebv", "popvar_uc", cfg$top_n),
  compare_pair(scores, "cross_mean_gebv", "simple_usefa_mean", cfg$top_n),
  compare_pair(scores, "cross_mean_gebv", "popvar_mu", cfg$top_n)
))

mean_delta <- scores$simple_usefa_mean - scores$cross_mean_gebv
uc_delta <- scores$simple_usefa - scores$uc_recomb_gebv
summary <- data.frame(
  n_parents = cfg$n_parents,
  n_markers = ncol(geno),
  n_pairs = nrow(scores),
  effect_reliability = effects$reliability,
  ridge_lambda = effects$lambda,
  simple_mean_delta_sd = stats::sd(mean_delta, na.rm = TRUE),
  simple_uc_delta_sd = stats::sd(uc_delta, na.rm = TRUE),
  simple_mean_delta_median = stats::median(mean_delta, na.rm = TRUE),
  simple_uc_delta_median = stats::median(uc_delta, na.rm = TRUE),
  stringsAsFactors = FALSE
)

dir.create(dirname(cfg$output_prefix), showWarnings = FALSE, recursive = TRUE)
write.csv(comparisons, paste0(cfg$output_prefix, "_comparisons.csv"), row.names = FALSE)
write.csv(summary, paste0(cfg$output_prefix, "_summary.csv"), row.names = FALSE)

message("Usefulness equivalence summary:")
print(summary)
message("Pairwise comparisons:")
print(comparisons)
