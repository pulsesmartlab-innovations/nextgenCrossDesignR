# AlphaSimR benchmark for external cross-prediction packages versus var_simple.
#
# Default benchmark:
# - DH/inbred parents, because PopVar deterministic cross prediction expects
#   parent marker genotypes coded as fully inbred alleles.
# - 80 parents, 5K SNPs, 3 replicates, 5 cycles.
# - var_simple and PopVar are run as selection branches.
# - SimpleMating methods can be enabled with EXT_METHODS, for example:
#   EXT_METHODS=var_simple,simple_mpv,simple_tgv,simple_usefa,simple_usefad_nonphased
#
# Optional:
# - Add genomicMateSelectR with EXT_METHODS=var_simple,popvar,gms_pmv.
#   The public predCrossVars() path is expensive for all 80-parent/5K-SNP
#   pairs; use EXT_GMS_PAIR_LIMIT to control whether it is screened or all.

.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
Sys.setenv(CROSSPRED_SKIP_CPP = Sys.getenv("EXT_SKIP_CPP", unset = "1"))
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
  family_metrics <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    p1 <- match(selected$parent1[i], parent_ids)
    p2 <- match(selected$parent2[i], parent_ids)
    if (is.na(p1) || is.na(p2)) {
      stop("Selected cross contains parent IDs not present in the parent population.", call. = FALSE)
    }
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

popvar_marker_matrix <- function(geno) {
  geno <- as.matrix(geno)
  if (any(!geno %in% c(0, 2))) {
    bad <- sort(unique(as.numeric(geno[!geno %in% c(0, 2)])))
    stop(
      "PopVar deterministic prediction needs inbred/DH parent genotypes coded 0/2. Found: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
  M <- geno
  M[M == 0] <- -1
  M[M == 2] <- 1
  storage.mode(M) <- "numeric"
  M
}

estimate_external_effects <- function(training_pop, sim_param, marker_names, y,
                                      method, draws, lambda, prior_scale) {
  training_geno <- pullSnpGeno(training_pop, simParam = sim_param)
  colnames(training_geno) <- marker_names
  rownames(training_geno) <- training_pop@id
  effects <- estimate_marker_effect_draws(
    geno_mat = training_geno,
    y = y,
    method = method,
    n_draws = draws,
    lambda = lambda,
    prior_genetic_var = stats::var(y) * prior_scale
  )
  beta <- attr(effects, "posterior_mean")
  if (is.null(beta)) beta <- colMeans(as.matrix(effects), na.rm = TRUE)
  beta <- as.numeric(beta)
  names(beta) <- marker_names
  list(effects = effects, beta = beta, training_geno = training_geno)
}

predict_popvar_crosses <- function(parent_values, geno, marker_map, beta, pair_df,
                                   tau, h2) {
  if (!requireNamespace("PopVar", quietly = TRUE)) {
    stop("Package 'PopVar' is not installed.", call. = FALSE)
  }
  ids <- rownames(geno)
  M <- popvar_marker_matrix(geno)
  map_in <- data.frame(
    marker = colnames(geno),
    chr = marker_map$chr,
    pos = marker_map$pos * 100,
    stringsAsFactors = FALSE
  )
  y_in <- data.frame(
    Entry = ids,
    Trait = as.numeric(parent_values[ids]),
    stringsAsFactors = FALSE
  )
  marker_effects <- data.frame(
    marker = colnames(geno),
    Trait = as.numeric(beta[colnames(geno)]),
    stringsAsFactors = FALSE
  )
  crossing_table <- data.frame(
    Par1 = as.character(pair_df$parent1),
    Par2 = as.character(pair_df$parent2),
    stringsAsFactors = FALSE
  )
  pv <- suppressWarnings(PopVar::pop_predict2(
    M = M,
    y.in = y_in,
    marker.effects = marker_effects,
    map.in = map_in,
    crossing.table = crossing_table,
    parents = ids,
    tail.p = tau,
    self.gen = 0,
    DH = TRUE,
    models = "rrBLUP"
  ))

  out <- data.frame(
    parent1 = as.character(pv$parent1),
    parent2 = as.character(pv$parent2),
    cross_mean_package = as.numeric(pv$pred_mu),
    cross_var = pmax(0, as.numeric(pv$pred_varG)),
    popvar_musp_high = as.numeric(pv$pred_musp_high),
    stringsAsFactors = FALSE
  )
  out$cross_mean <- parent_cross_mean(parent_values, out$parent1, out$parent2, ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "parent_value"
  add_usefulness_criterion(out, tau = tau, h2 = h2, trait_direction = "increase")
}

extract_gms_predvar <- function(pred_vars) {
  vapply(pred_vars, function(x) {
    xx <- as.data.frame(x)
    as.numeric(xx$predVar[1])
  }, numeric(1))
}

predict_gms_crosses <- function(parent_values, geno, haplo, marker_map, beta, pair_df,
                                tau, h2, ncores = 1L) {
  if (!requireNamespace("genomicMateSelectR", quietly = TRUE)) {
    stop("Package 'genomicMateSelectR' is not installed.", call. = FALSE)
  }
  ids <- rownames(geno)
  effect_mat <- matrix(
    as.numeric(beta[colnames(geno)]),
    nrow = 1,
    dimnames = list("Trait", colnames(geno))
  )
  dose_mat <- as.matrix(geno)
  storage.mode(dose_mat) <- "numeric"
  crosses <- data.frame(
    sireID = as.character(pair_df$parent1),
    damID = as.character(pair_df$parent2),
    stringsAsFactors = FALSE
  )
  means <- genomicMateSelectR::predCrossMeans(
    CrossesToPredict = crosses,
    predType = "BV",
    AddEffectList = list(Trait = effect_mat),
    doseMat = dose_mat,
    ncores = ncores
  )
  recomb_pos <- marker_map$pos * 100
  names(recomb_pos) <- colnames(geno)
  recomb_freq <- genomicMateSelectR::genmap2recombfreq(recomb_pos, nChr = length(unique(marker_map$chr)))
  linkage_mat <- 1 - 2 * recomb_freq
  diag(linkage_mat) <- 1
  vars <- genomicMateSelectR::predCrossVars(
    CrossesToPredict = crosses,
    modelType = "A",
    AddEffectList = list(Trait = effect_mat),
    predType = "VPM",
    haploMat = haplo,
    recombFreqMat = linkage_mat,
    ncores = ncores
  )
  means <- as.data.frame(means)
  vars <- as.data.frame(vars)
  out <- data.frame(
    parent1 = as.character(vars$sireID),
    parent2 = as.character(vars$damID),
    cross_mean_package = as.numeric(means$predMean),
    cross_var = pmax(0, extract_gms_predvar(vars$predVars)),
    n_seg_snps = as.integer(vars$Nsegsnps),
    stringsAsFactors = FALSE
  )
  out$cross_mean <- parent_cross_mean(parent_values, out$parent1, out$parent2, ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "parent_value"
  add_usefulness_criterion(out, tau = tau, h2 = h2, trait_direction = "increase")
}

screen_gms_pairs <- function(var_df, popvar_df, all_pairs, limit) {
  if (!is.finite(limit) || limit >= nrow(all_pairs)) return(all_pairs)
  limit <- max(1L, as.integer(limit))
  screens <- list()
  if (!is.null(var_df) && nrow(var_df)) {
    v <- var_df[order(var_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
    screens[[length(screens) + 1L]] <- utils::head(v[, c("parent1", "parent2"), drop = FALSE], limit)
  }
  if (!is.null(popvar_df) && nrow(popvar_df)) {
    p <- popvar_df[order(popvar_df$cross_usefulness, decreasing = TRUE), , drop = FALSE]
    screens[[length(screens) + 1L]] <- utils::head(p[, c("parent1", "parent2"), drop = FALSE], limit)
  }
  if (!length(screens)) return(utils::head(all_pairs, limit))
  unique(do.call(rbind, screens))
}

is_simplemating_method <- function(method) {
  method %in% c(
    "simple_mpv",
    "simple_tgv",
    "simple_usefa",
    "simple_usefad_nonphased",
    "simple_usefad_phased"
  )
}

as_simplemating_plan <- function(pair_df) {
  data.frame(
    Parent1 = as.character(pair_df$parent1),
    Parent2 = as.character(pair_df$parent2),
    stringsAsFactors = FALSE
  )
}

as_simplemating_map <- function(marker_map) {
  data.frame(
    Chromosome = marker_map$chr,
    Position = marker_map$pos * 100,
    Marker = marker_map$marker,
    stringsAsFactors = FALSE
  )
}

as_simplemating_k <- function(geno) {
  K <- VanRadenKin(geno)
  K <- as.matrix(K)
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  K
}

quiet_simplemating <- function(expr) {
  value <- NULL
  utils::capture.output({
    value <- force(expr)
  })
  value
}

simplemating_mean_df <- function(sm_df, parent_values, ids, score_col = "Y",
                                 package_mean_col = "Y") {
  out <- data.frame(
    parent1 = as.character(sm_df$Parent1),
    parent2 = as.character(sm_df$Parent2),
    cross_mean_package = as.numeric(sm_df[[package_mean_col]]),
    cross_var = 0,
    cross_sigma = 0,
    cross_usefulness = as.numeric(sm_df[[score_col]]),
    simplemating_k = if ("K" %in% names(sm_df)) as.numeric(sm_df$K) else NA_real_,
    stringsAsFactors = FALSE
  )
  out$cross_mean <- parent_cross_mean(parent_values, out$parent1, out$parent2, ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "simplemating_package_score"
  out
}

simplemating_usefulness_df <- function(sm_detail, parent_values, ids,
                                       var_cols, score_col = "Usefulness") {
  cross_var <- Reduce(`+`, lapply(var_cols, function(v) {
    if (v %in% names(sm_detail)) as.numeric(sm_detail[[v]]) else rep(0, nrow(sm_detail))
  }))
  out <- data.frame(
    parent1 = as.character(sm_detail$Parent1),
    parent2 = as.character(sm_detail$Parent2),
    cross_mean_package = as.numeric(sm_detail$Mean),
    cross_var = pmax(0, cross_var),
    cross_usefulness = as.numeric(sm_detail[[score_col]]),
    stringsAsFactors = FALSE
  )
  out$cross_sigma <- sqrt(pmax(0, out$cross_var))
  out$cross_mean <- parent_cross_mean(parent_values, out$parent1, out$parent2, ids)
  out$cross_mean_marker_effects <- out$cross_mean_package
  out$mean_source <- "simplemating_package_score"
  out
}

predict_simplemating_crosses <- function(method, parent_values, geno, haplo, marker_map,
                                         beta, pair_df, tau) {
  if (!requireNamespace("SimpleMating", quietly = TRUE)) {
    stop("Package 'SimpleMating' is not installed.", call. = FALSE)
  }
  original_ids <- rownames(geno)
  simple_ids <- paste0("SM", seq_along(original_ids))
  names(simple_ids) <- original_ids
  original_from_simple <- stats::setNames(original_ids, simple_ids)

  geno <- as.matrix(geno)
  rownames(geno) <- simple_ids
  if (!is.null(haplo)) {
    haplo <- as.matrix(haplo)
    rownames(haplo) <- as.vector(rbind(paste0(simple_ids, "_HapA"), paste0(simple_ids, "_HapB")))
  }
  parent_values <- stats::setNames(as.numeric(parent_values[original_ids]), simple_ids)
  pair_df <- data.frame(
    parent1 = unname(simple_ids[as.character(pair_df$parent1)]),
    parent2 = unname(simple_ids[as.character(pair_df$parent2)]),
    stringsAsFactors = FALSE
  )
  ids <- rownames(geno)
  plan <- as_simplemating_plan(pair_df)
  K <- as_simplemating_k(geno)
  sm_map <- as_simplemating_map(marker_map)
  beta <- as.numeric(beta[colnames(geno)])
  names(beta) <- colnames(geno)
  dom <- rep(0, length(beta))
  names(dom) <- names(beta)

  if (method == "simple_mpv") {
    crit <- data.frame(
      Id = ids,
      Criterion = as.numeric(parent_values[ids]),
      stringsAsFactors = FALSE
    )
    sm <- quiet_simplemating(SimpleMating::getMPV(
      MatePlan = plan,
      Criterion = crit,
      K = K
    ))
    out <- simplemating_mean_df(sm, parent_values, ids)
    out$parent1 <- unname(original_from_simple[out$parent1])
    out$parent2 <- unname(original_from_simple[out$parent2])
    return(out)
  }

  if (method == "simple_tgv") {
    sm <- quiet_simplemating(SimpleMating::getTGV(
      MatePlan = plan,
      Markers = geno,
      addEff = beta,
      domEff = dom,
      K = K,
      ploidy = 2
    ))
    out <- simplemating_mean_df(sm, parent_values, ids)
    out$parent1 <- unname(original_from_simple[out$parent1])
    out$parent2 <- unname(original_from_simple[out$parent2])
    return(out)
  }

  if (method == "simple_usefa") {
    sm <- quiet_simplemating(SimpleMating::getUsefA(
      MatePlan = plan,
      Markers = geno,
      addEff = beta,
      K = K,
      Map.In = sm_map,
      propSel = tau,
      Type = "DH",
      Generation = 1,
      n_threads = 1,
      display_progress = FALSE
    ))
    out <- simplemating_usefulness_df(sm[[1]], parent_values, ids, var_cols = "Variance")
    out$parent1 <- unname(original_from_simple[out$parent1])
    out$parent2 <- unname(original_from_simple[out$parent2])
    return(out)
  }

  if (method == "simple_usefad_nonphased") {
    sm <- quiet_simplemating(SimpleMating::getUsefAD(
      MatePlan = plan,
      Markers = geno,
      addEff = beta,
      domEff = dom,
      K = K,
      Map.In = sm_map,
      propSel = tau,
      Method = "NonPhased",
      ploidy = 2,
      n_threads = 1,
      display_progress = FALSE
    ))
    out <- simplemating_usefulness_df(sm[[1]], parent_values, ids, var_cols = c("Var_A", "Var_D"))
    out$parent1 <- unname(original_from_simple[out$parent1])
    out$parent2 <- unname(original_from_simple[out$parent2])
    return(out)
  }

  if (method == "simple_usefad_phased") {
    sm <- quiet_simplemating(SimpleMating::getUsefAD(
      MatePlan = plan,
      Markers = haplo,
      addEff = beta,
      domEff = dom,
      K = K,
      Map.In = sm_map,
      propSel = tau,
      Method = "Phased",
      ploidy = 2,
      n_threads = 1,
      display_progress = FALSE
    ))
    out <- simplemating_usefulness_df(sm[[1]], parent_values, ids, var_cols = c("Var_A", "Var_D"))
    out$parent1 <- unname(original_from_simple[out$parent1])
    out$parent2 <- unname(original_from_simple[out$parent2])
    return(out)
  }

  stop("Unsupported SimpleMating method: ", method, call. = FALSE)
}

package_status_table <- function(methods, gms_pair_limit) {
  pkgs <- c("PopVar", "genomicMateSelectR", "GenomicMating", "SimpleMating")
  installed <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  has_simplemating <- any(vapply(methods, is_simplemating_method, logical(1)))
  notes <- c(
    "Benchmarked when method 'popvar' is enabled; deterministic pop_predict2() with external ridge effects.",
    if ("gms_pmv" %in% methods) {
      if (is.finite(gms_pair_limit)) {
        paste0("Benchmarked as a screened branch with EXT_GMS_PAIR_LIMIT=", gms_pair_limit, ".")
      } else {
        "Benchmarked on all parent pairs."
      }
    } else {
      "Available but not enabled by default because all-pair predCrossVars() is very slow at 80 parents and 5K SNPs."
    },
    "Not loadable in this R installation from the attempted standard install path.",
    if (installed[["SimpleMating"]]) {
      "Benchmarked when any simple_* method is enabled; local clone patched with C++14 Makevars for this RcppArmadillo toolchain."
    } else {
      "Not loadable here; GitHub install failed against current RcppArmadillo/C++ standard requirements."
    }
  )
  data.frame(
    package = pkgs,
    installed = as.logical(installed),
    benchmarked = c("popvar" %in% methods && installed[["PopVar"]],
                    "gms_pmv" %in% methods && installed[["genomicMateSelectR"]],
                    FALSE,
                    has_simplemating && installed[["SimpleMating"]]),
    note = notes,
    stringsAsFactors = FALSE
  )
}

cfg <- list(
  seed = env_int("EXT_SEED", 20260426L),
  reps = env_int("EXT_REPS", 3L),
  cycles = env_int("EXT_CYCLES", 5L),
  n_founders = env_int("EXT_N_FOUNDERS", 120L),
  n_parents = env_int("EXT_N_PARENTS", 80L),
  n_chr = env_int("EXT_N_CHR", 5L),
  seg_sites = env_int("EXT_SEG_SITES", 1200L),
  qtl_per_chr = env_int("EXT_QTL_PER_CHR", 40L),
  snp_per_chr = env_int("EXT_SNP_PER_CHR", 1000L),
  top_crosses = env_int("EXT_TOP_CROSSES", 20L),
  progeny_per_cross = env_int("EXT_PROGENY_PER_CROSS", 40L),
  effect_training_n = env_int("EXT_EFFECT_TRAINING_N", env_int("EXT_N_PARENTS", 80L)),
  effect_method = env_chr("EXT_EFFECT_METHOD", "ridge_posterior"),
  effect_draws = env_int("EXT_EFFECT_DRAWS", 20L),
  ridge_lambda = env_num("EXT_RIDGE_LAMBDA", 1),
  prior_genetic_var_scale = env_num("EXT_PRIOR_GENETIC_VAR_SCALE", 1),
  tau = env_num("EXT_TAU", 0.10),
  h2 = env_num("EXT_H2", 1.0),
  methods = strsplit(env_chr("EXT_METHODS", "var_simple,popvar"), ",", fixed = TRUE)[[1]],
  gms_pair_limit = env_count("EXT_GMS_PAIR_LIMIT", 300L),
  simple_pair_limit = env_count("EXT_SIMPLE_PAIR_LIMIT", Inf),
  gms_ncores = env_int("EXT_GMS_NCORES", 1L),
  output_prefix = env_chr("EXT_OUTPUT_PREFIX", "external_cross_package")
)
cfg$methods <- trimws(cfg$methods)
cfg$methods <- unique(cfg$methods[nzchar(cfg$methods)])
if (!"var_simple" %in% cfg$methods) cfg$methods <- c("var_simple", cfg$methods)

if ("popvar" %in% cfg$methods && !requireNamespace("PopVar", quietly = TRUE)) {
  warning("Removing method 'popvar' because package PopVar is not installed.")
  cfg$methods <- setdiff(cfg$methods, "popvar")
}
if ("gms_pmv" %in% cfg$methods && !requireNamespace("genomicMateSelectR", quietly = TRUE)) {
  warning("Removing method 'gms_pmv' because package genomicMateSelectR is not installed.")
  cfg$methods <- setdiff(cfg$methods, "gms_pmv")
}
if (any(vapply(cfg$methods, is_simplemating_method, logical(1))) &&
    !requireNamespace("SimpleMating", quietly = TRUE)) {
  warning("Removing SimpleMating methods because package SimpleMating is not installed.")
  cfg$methods <- cfg$methods[!vapply(cfg$methods, is_simplemating_method, logical(1))]
}
valid_methods <- c(
  "var_simple",
  "popvar",
  "gms_pmv",
  "simple_mpv",
  "simple_tgv",
  "simple_usefa",
  "simple_usefad_nonphased",
  "simple_usefad_phased"
)
unknown_methods <- setdiff(cfg$methods, valid_methods)
if (length(unknown_methods)) {
  stop("Unknown EXT_METHODS values: ", paste(unknown_methods, collapse = ", "), call. = FALSE)
}

dir.create("results", showWarnings = FALSE)
status <- package_status_table(cfg$methods, cfg$gms_pair_limit)
write.csv(status, file.path("results", paste0(cfg$output_prefix, "_package_status.csv")), row.names = FALSE)

config_df <- data.frame(name = names(cfg), value = vapply(cfg, function(x) paste(x, collapse = ","), character(1)))
write.csv(config_df, file.path("results", paste0(cfg$output_prefix, "_config.csv")), row.names = FALSE)

message("External package benchmark config:")
print(cfg)
message("External package status:")
print(status)

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
  marker_map <- data.frame(
    marker = marker_names,
    chr = snp_map$chr,
    pos = snp_map$pos,
    stringsAsFactors = FALSE
  )

  branches <- setNames(vector("list", length(cfg$methods)), cfg$methods)
  for (m in cfg$methods) branches[[m]] <- base_parents

  metrics <- list()
  selections <- list()
  timings <- list()
  family_diagnostics <- list()
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
      ids <- paste0(method, "_R", rep_id, "_C", cycle - 1L, "_P", seq_len(nInd(parent_pop)))
      parent_pop@id <- ids
      parent_gv <- as.numeric(gv(parent_pop)[, 1])
      names(parent_gv) <- ids

      geno <- pullSnpGeno(parent_pop, simParam = SP)
      colnames(geno) <- marker_names
      rownames(geno) <- ids

      all_pairs <- build_parent_pairs(ids, include_self = FALSE)
      var_screen <- NULL
      popvar_screen <- NULL
      effects <- NULL
      beta <- NULL

      if (method == "var_simple" || method == "gms_pmv") {
        elapsed_var <- system.time({
          var_screen <- build_cross_data(
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
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle,
          stage = if (method == "var_simple") "score" else "var_simple_screen",
          elapsed_sec = elapsed_var,
          stringsAsFactors = FALSE
        )
      }

      if (method %in% c("popvar", "gms_pmv", "simple_tgv", "simple_usefa",
                        "simple_usefad_nonphased", "simple_usefad_phased")) {
        elapsed_effect <- system.time({
          training_pop <- make_training_population(parent_pop, cfg$effect_training_n, SP)
          training_ids <- paste0(method, "_R", rep_id, "_C", cycle - 1L, "_T", seq_len(nInd(training_pop)))
          training_pop@id <- training_ids
          training_y <- as.numeric(gv(training_pop)[, 1])
          external_effects <- estimate_external_effects(
            training_pop = training_pop,
            sim_param = SP,
            marker_names = marker_names,
            y = training_y,
            method = cfg$effect_method,
            draws = cfg$effect_draws,
            lambda = cfg$ridge_lambda,
            prior_scale = cfg$prior_genetic_var_scale
          )
          effects <- external_effects$effects
          beta <- external_effects$beta
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "effects",
          elapsed_sec = elapsed_effect, stringsAsFactors = FALSE
        )

        effect_diag <- effect_scale_diagnostics(effects, geno, parent_gv, cv_folds = 0L)
        effect_diag$rep <- rep_id
        effect_diag$method <- method
        effect_diag$cycle <- cycle
        effect_diagnostics[[paste(rep_id, method, cycle, sep = "_")]] <- effect_diag
      }

      if (method == "var_simple") {
        cross_df <- var_screen
      } else if (method == "popvar") {
        elapsed <- system.time({
          cross_df <- predict_popvar_crosses(
            parent_values = parent_gv,
            geno = geno,
            marker_map = marker_map,
            beta = beta,
            pair_df = all_pairs,
            tau = cfg$tau,
            h2 = cfg$h2
          )
        })[["elapsed"]]
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "score",
          elapsed_sec = elapsed, stringsAsFactors = FALSE
        )
      } else if (method == "gms_pmv") {
        if ("popvar" %in% cfg$methods) {
          popvar_screen <- tryCatch(
            predict_popvar_crosses(
              parent_values = parent_gv,
              geno = geno,
              marker_map = marker_map,
              beta = beta,
              pair_df = all_pairs,
              tau = cfg$tau,
              h2 = cfg$h2
            ),
            error = function(e) NULL
          )
        }
        gms_pairs <- screen_gms_pairs(var_screen, popvar_screen, all_pairs, cfg$gms_pair_limit)
        elapsed <- system.time({
          cross_df <- predict_gms_crosses(
            parent_values = parent_gv,
            geno = geno,
            haplo = as_snp_haplo(parent_pop, SP, marker_names),
            marker_map = marker_map,
            beta = beta,
            pair_df = gms_pairs,
            tau = cfg$tau,
            h2 = cfg$h2,
            ncores = cfg$gms_ncores
          )
        })[["elapsed"]]
        cross_df$screened_n_pairs <- nrow(gms_pairs)
        cross_df$screened_from_all_pairs <- nrow(all_pairs)
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "score",
          elapsed_sec = elapsed, stringsAsFactors = FALSE
        )
      } else if (is_simplemating_method(method)) {
        simple_pairs <- all_pairs
        if (method %in% c("simple_usefa", "simple_usefad_nonphased", "simple_usefad_phased") &&
            is.finite(cfg$simple_pair_limit) && cfg$simple_pair_limit < nrow(all_pairs)) {
          elapsed_screen <- system.time({
            if (is.null(var_screen)) {
              var_screen <- build_cross_data(
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
            }
            tgv_screen <- predict_simplemating_crosses(
              method = "simple_tgv",
              parent_values = parent_gv,
              geno = geno,
              haplo = NULL,
              marker_map = marker_map,
              beta = beta,
              pair_df = all_pairs,
              tau = cfg$tau
            )
            n_each <- max(cfg$top_crosses, as.integer(cfg$simple_pair_limit))
            v <- var_screen[order(var_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]
            t <- tgv_screen[order(tgv_screen$cross_usefulness, decreasing = TRUE), , drop = FALSE]
            simple_pairs <- unique(rbind(
              utils::head(v[, c("parent1", "parent2"), drop = FALSE], n_each),
              utils::head(t[, c("parent1", "parent2"), drop = FALSE], n_each)
            ))
          })[["elapsed"]]
          timings[[length(timings) + 1L]] <- data.frame(
            rep = rep_id, method = method, cycle = cycle, stage = "screen",
            elapsed_sec = elapsed_screen, stringsAsFactors = FALSE
          )
        }
        elapsed <- system.time({
          cross_df <- predict_simplemating_crosses(
            method = method,
            parent_values = parent_gv,
            geno = geno,
            haplo = if (method == "simple_usefad_phased") as_snp_haplo(parent_pop, SP, marker_names) else NULL,
            marker_map = marker_map,
            beta = if (is.null(beta)) rep(0, ncol(geno)) else beta,
            pair_df = simple_pairs,
            tau = cfg$tau
          )
        })[["elapsed"]]
        if (nrow(simple_pairs) < nrow(all_pairs)) {
          cross_df$screened_n_pairs <- nrow(simple_pairs)
          cross_df$screened_from_all_pairs <- nrow(all_pairs)
        }
        timings[[length(timings) + 1L]] <- data.frame(
          rep = rep_id, method = method, cycle = cycle, stage = "score",
          elapsed_sec = elapsed, stringsAsFactors = FALSE
        )
      }

      selected <- select_crosses(cross_df, cfg$top_crosses)
      selected$rep <- rep_id
      selected$method <- method
      selected$cycle <- cycle
      selections[[paste(rep_id, method, cycle, sep = "_")]] <- selected

      progeny_result <- make_selected_progeny(parent_pop, selected, cfg$progeny_per_cross, SP)
      family_diag <- progeny_result$family_metrics
      family_diag$rep <- rep_id
      family_diag$method <- method
      family_diag$cycle <- cycle
      family_diagnostics[[paste(rep_id, method, cycle, sep = "_")]] <- family_diag

      progeny <- progeny_result$pop
      g <- as.numeric(gv(progeny)[, 1])
      keep <- order(g, decreasing = TRUE)[seq_len(min(cfg$n_parents, length(g)))]
      next_pop <- progeny[keep]
      next_pop@id <- paste0(method, "_R", rep_id, "_C", cycle, "_P", seq_len(nInd(next_pop)))
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
    family_diagnostics = bind_rows_fill(family_diagnostics),
    effect_diagnostics = bind_rows_fill(effect_diagnostics)
  )
})

metrics <- bind_rows_fill(lapply(rep_results, `[[`, "metrics"))
selections <- bind_rows_fill(lapply(rep_results, `[[`, "selections"))
timings <- bind_rows_fill(lapply(rep_results, `[[`, "timings"))
family_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "family_diagnostics"))
effect_diagnostics <- bind_rows_fill(lapply(rep_results, `[[`, "effect_diagnostics"))

overall <- stats::aggregate(
  cbind(mean_gv, max_gv, top10_gv, var_gv, transgressive_rate) ~ method + cycle,
  metrics,
  mean
)
baseline <- overall[overall$method == "var_simple",
                    c("cycle", "mean_gv", "max_gv", "top10_gv", "transgressive_rate"),
                    drop = FALSE]
comparison <- merge(overall, baseline, by = "cycle", suffixes = c("", "_var_simple"))
comparison$delta_mean_gv <- comparison$mean_gv - comparison$mean_gv_var_simple
comparison$delta_max_gv <- comparison$max_gv - comparison$max_gv_var_simple
comparison$delta_top10_gv <- comparison$top10_gv - comparison$top10_gv_var_simple
comparison$delta_transgressive_rate <- comparison$transgressive_rate - comparison$transgressive_rate_var_simple

timing_summary <- if (nrow(timings)) {
  stats::aggregate(elapsed_sec ~ method + stage, timings, mean)
} else {
  data.frame()
}
family_summary <- if (nrow(family_diagnostics)) {
  stats::aggregate(
    cbind(pred_var, realized_var_gv, realized_max_gv, realized_top10_gv) ~ method + cycle,
    family_diagnostics,
    mean
  )
} else {
  data.frame()
}

prefix <- file.path("results", cfg$output_prefix)
write.csv(metrics, paste0(prefix, "_metrics_by_rep.csv"), row.names = FALSE)
write.csv(overall, paste0(prefix, "_overall.csv"), row.names = FALSE)
write.csv(comparison, paste0(prefix, "_comparison_vs_var_simple.csv"), row.names = FALSE)
write.csv(selections, paste0(prefix, "_selected_crosses.csv"), row.names = FALSE)
write.csv(timings, paste0(prefix, "_timing.csv"), row.names = FALSE)
write.csv(timing_summary, paste0(prefix, "_timing_summary.csv"), row.names = FALSE)
write.csv(family_diagnostics, paste0(prefix, "_family_diagnostics.csv"), row.names = FALSE)
write.csv(family_summary, paste0(prefix, "_family_summary.csv"), row.names = FALSE)
write.csv(effect_diagnostics, paste0(prefix, "_effect_diagnostics.csv"), row.names = FALSE)

message("External package benchmark overall:")
print(overall)
message("Comparison versus var_simple:")
print(comparison[order(comparison$cycle, comparison$method), , drop = FALSE])
message("Timing summary:")
print(timing_summary)
if (nrow(family_summary)) {
  message("Selected-family diagnostics:")
  print(family_summary)
}
