#!/usr/bin/env Rscript
# Native mate-allocation study -- recurrent genomic-selection benchmark of the package's own
# allocators. Builds OUTBRED populations WITH DOMINANCE for arbitrary ploidy, scores candidate
# crosses, and races:
#   ng_useful    - usefulness allocation (mean + i*within-family SD)   [variance-aware]
#   ng_mean      - mean-only allocation (control)
# reporting realized gain / usefulness / variance retained / group coancestry.
#
# THIS IS A PILOT-FIRST RESEARCH RUNNER (not part of the installed package). In PILOT mode
# (default) it does ONE outbred generation per ploidy on true genetic values -- validating the
# scenario (dominance signal) and the allocation plumbing before any long recurrent run.
# It uses TRUE values throughout (allocate-on-true, evaluate-on-true), so it grades ALLOCATION,
# not prediction; the recurrent GEBV-based study (Phase 4) adds the estimation layer.
#
# Env config (all optional):
#   NG_MATEALLOC_PLOIDY("2,4")  NG_MATEALLOC_SEED(20260703)  NG_MATEALLOC_N_PARENTS(14)
#   NG_MATEALLOC_N_CROSSES(12)  NG_MATEALLOC_MAX_USE(4)      NG_MATEALLOC_DF(0.05)
#   NG_MATEALLOC_N_SCORE_PROGENY(40)  NG_MATEALLOC_SELECTION_PROP(0.10)
#   NG_MATEALLOC_MEAN_DD(0.5)   NG_MATEALLOC_VAR_DD(0.2)     NG_MATEALLOC_N_CHR(5)
#   NG_MATEALLOC_SEG_SITES(400) NG_MATEALLOC_QTL_PER_CHR(18) NG_MATEALLOC_SNP_PER_CHR(80)
#   NG_MATEALLOC_N_FOUNDERS(60) NG_MATEALLOC_OUTPUT_DIR(results/) NG_MATEALLOC_USE_CPP(0)

find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    start, file.path(start, "nextgen_cross_design"), dirname(start),
    file.path(dirname(start), "nextgen_cross_design"),
    file.path(dirname(dirname(start)), "nextgen_cross_design")
  ), winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    desc <- file.path(candidate, "DESCRIPTION")
    if (file.exists(file.path(candidate, "R", "load.R")) && file.exists(desc) &&
        any(grepl("^Package:\\s*nextgenCrossDesign\\s*$", readLines(desc, n = 20L, warn = FALSE)))) {
      return(candidate)
    }
  }
  stop("Could not locate nextgenCrossDesign root", call. = FALSE)
}

env_chr <- function(name, default) { v <- Sys.getenv(name, unset = ""); if (!nzchar(v)) default else v }
env_int <- function(name, default) { v <- Sys.getenv(name, unset=""); o <- suppressWarnings(as.integer(v)); if (!nzchar(v)||!is.finite(o)) as.integer(default) else o }
env_num <- function(name, default) { v <- Sys.getenv(name, unset=""); o <- suppressWarnings(as.numeric(v)); if (!nzchar(v)||!is.finite(o)) as.numeric(default) else o }
env_int_vec <- function(name, default) { v <- Sys.getenv(name, unset=""); if (!nzchar(v)) return(default); o <- suppressWarnings(as.integer(strsplit(v,"[,;[:space:]]+")[[1]])); o[is.finite(o)] }

root <- find_project_root()
suppressWarnings(suppressMessages({
  source(file.path(root, "tools", "ng_project_libpath.R"))
  ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
  source(file.path(root, "R", "load.R"))
  ng_load(root, use_cpp = identical(env_chr("NG_MATEALLOC_USE_CPP","0"),"1"), verbose = FALSE)
}))
if (!requireNamespace("AlphaSimR", quietly = TRUE)) stop("AlphaSimR is required for the mate-allocation study")

# ------- simulation: outbred founders WITH dominance, arbitrary ploidy -------
# Founder engine defaults to runMacs (coalescent) for REALISTIC marker-QTL LD -- essential for
# genomic prediction (quickHaplo gives ~0 LD and useless GEBV; the metric-merit study hit this).
# quickHaplo is available via NG_MATEALLOC_FOUNDER=quickhaplo for fast plumbing checks.
sim_setup <- function(ploidy, cfg, seed) {
  set.seed(seed)
  founder <- if (identical(cfg$founder_engine, "quickhaplo")) {
    AlphaSimR::quickHaplo(nInd = cfg$n_founders, nChr = cfg$n_chr, segSites = cfg$seg_sites,
                          genLen = 1.0, ploidy = as.integer(ploidy), inbred = FALSE)
  } else {
    AlphaSimR::runMacs(nInd = cfg$n_founders, nChr = cfg$n_chr, segSites = cfg$seg_sites,
                       ploidy = as.integer(ploidy), species = "GENERIC")
  }
  SP <- AlphaSimR::SimParam$new(founder)
  if (ploidy == 4L) SP$quadProb <- 0.1
  # addTraitAD adds directional (meanDD>0) + variable (varDD) dominance -> mid-parent heterosis.
  SP$addTraitAD(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1,
                meanDD = cfg$mean_dd, varDD = cfg$var_dd)
  SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  list(founder = founder, SP = SP)
}

make_parents <- function(setup, n_parents) {
  pop <- AlphaSimR::newPop(setup$founder, simParam = setup$SP)
  if (AlphaSimR::nInd(pop) < n_parents) stop("not enough founders for requested parents")
  pop <- pop[seq_len(n_parents)]
  pop@id <- sprintf("P%03d", seq_len(n_parents))
  pop
}

# VanRaden-style genomic relationship (diagonal ~1) from allele dosage, any ploidy.
genomic_G <- function(dosage) {
  X <- sweep(dosage, 2L, colMeans(dosage), "-")
  denom <- sum(apply(X, 2L, stats::var))
  if (!is.finite(denom) || denom <= 0) denom <- ncol(X)
  K <- tcrossprod(X) / denom
  dimnames(K) <- list(rownames(dosage), rownames(dosage))
  K
}

# Score every candidate cross on TRUE genetic values (mean, within-family var, usefulness),
# plus mid-parent value and heterosis (progeny mean - mid-parent). Returns scores + attrs.
score_crosses <- function(parent_pop, SP, ploidy, cfg, seed) {
  ids <- as.character(parent_pop@id)
  dosage <- AlphaSimR::pullSnpGeno(parent_pop, simParam = SP)
  rownames(dosage) <- ids
  G <- genomic_G(dosage)
  K_coanc <- G / ploidy                       # coancestry scale (diagonal ~1/ploidy)
  parent_gv <- as.numeric(AlphaSimR::gv(parent_pop)[, 1]); names(parent_gv) <- ids
  intensity <- ng_selection_intensity(cfg$selection_prop)
  pairs <- t(utils::combn(length(ids), 2L))
  out <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    a <- ids[pairs[i, 1]]; b <- ids[pairs[i, 2]]
    set.seed(as.integer(((seed %% 30011L) * 7919L + i) %% 2147483647L))
    fam <- AlphaSimR::makeCross(parent_pop, matrix(c(pairs[i, 1], pairs[i, 2]), nrow = 1L),
                                nProgeny = cfg$n_score_progeny, simParam = SP)
    gv <- as.numeric(AlphaSimR::gv(fam)[, 1])
    m <- mean(gv); v <- if (length(gv) > 1L) stats::var(gv) else 0
    mid <- (parent_gv[[a]] + parent_gv[[b]]) / 2
    out[[i]] <- data.frame(parent1 = a, parent2 = b, cross_mean = m, cross_var = v,
                           cross_usefulness = m + intensity * sqrt(max(v, 0)),
                           midparent = mid, heterosis = m - mid,
                           pair_kinship = K_coanc[a, b], stringsAsFactors = FALSE)
  }
  scores <- do.call(rbind, out)
  attr(scores, "parent_kinship") <- K_coanc
  attr(scores, "parent_G") <- G
  scores
}

# realized group coancestry c'Kc for a selected plan (same summary our allocator reports).
group_coancestry <- function(plan, K) {
  cnt <- table(factor(c(plan$parent1, plan$parent2), levels = rownames(K)))
  c_vec <- as.numeric(cnt) / sum(cnt)
  as.numeric(t(c_vec) %*% K %*% c_vec)
}
plan_key <- function(plan) sort(ng_group_pair_key(plan$parent1, plan$parent2))
jaccard <- function(a, b) { u <- length(union(a, b)); if (!u) return(NA_real_); length(intersect(a, b)) / u }

summarize_plan <- function(scores, plan, K, label, method) {
  key <- ng_group_pair_key(scores$parent1, scores$parent2)
  idx <- match(ng_group_pair_key(plan$parent1, plan$parent2), key)
  idx <- idx[!is.na(idx)]
  sel <- scores[idx, , drop = FALSE]
  data.frame(method = method, label = label, n = nrow(sel),
             mean = mean(sel$cross_mean), usefulness = mean(sel$cross_usefulness),
             var = mean(sel$cross_var), heterosis = mean(sel$heterosis),
             group_coancestry = group_coancestry(sel, K), stringsAsFactors = FALSE)
}

run_ploidy <- function(ploidy, cfg, seed) {
  cat(sprintf("\n===== PLOIDY %d (outbred, dominance meanDD=%.2f varDD=%.2f) =====\n",
              ploidy, cfg$mean_dd, cfg$var_dd))
  setup <- sim_setup(ploidy, cfg, seed)
  parents <- make_parents(setup, cfg$n_parents)
  scores <- score_crosses(parents, setup$SP, ploidy, cfg, seed)
  K <- attr(scores, "parent_kinship")
  cat(sprintf("candidate crosses: %d | cross mean range [%.2f, %.2f] | heterosis mean %.3f (range [%.3f, %.3f])\n",
      nrow(scores), min(scores$cross_mean), max(scores$cross_mean),
      mean(scores$heterosis), min(scores$heterosis), max(scores$heterosis)))

  plans <- list()
  # our usefulness allocation (variance-aware) at a coancestry target = dF
  plans$ng_useful <- tryCatch(ng_optimize_mating_plan(scores, cfg$n_crosses, gain_col = "cross_usefulness",
      parent_kinship = K, max_crosses_per_parent = cfg$max_use, target_coancestry = cfg$dF, method = "greedy_local"),
      error = function(e) { cat("ng_useful ERR:", conditionMessage(e), "\n"); NULL })
  # our mean-only allocation (control) at the same coancestry target
  plans$ng_mean <- tryCatch(ng_optimize_mating_plan(scores, cfg$n_crosses, gain_col = "cross_mean",
      parent_kinship = K, max_crosses_per_parent = cfg$max_use, target_coancestry = cfg$dF, method = "greedy_local"),
      error = function(e) { cat("ng_mean ERR:", conditionMessage(e), "\n"); NULL })

  rows <- list()
  for (nm in names(plans)) if (!is.null(plans[[nm]])) rows[[nm]] <- summarize_plan(scores, plans[[nm]], K, sprintf("ploidy%d", ploidy), nm)
  tab <- do.call(rbind, rows)
  cat("\n-- realized (true-value) allocation comparison --\n"); print(tab, row.names = FALSE, digits = 4)
  # pairwise set overlap
  keys <- lapply(plans[!vapply(plans, is.null, logical(1))], plan_key)
  if (length(keys) >= 2L) {
    cat("\n-- selected-set Jaccard overlap --\n")
    nmk <- names(keys)
    for (i in 1:(length(nmk)-1)) for (j in (i+1):length(nmk))
      cat(sprintf("  %-10s vs %-10s : %.2f\n", nmk[i], nmk[j], jaccard(keys[[i]], keys[[j]])))
  }
  tab
}

# ============================ RECURRENT MODE (Phase 4) ============================
# Full recurrent genomic-selection study: allocate on ESTIMATED breeding values, evaluate on
# TRUE genetic values, over multiple cycles and reps. Each method runs an INDEPENDENT breeding
# lineage from the SAME founders (same rep seed) -> methods differ ONLY in mate allocation, so
# realized-gain differences are attributable to the allocator, not to parent selection.
#
# Estimation: package ng_fit_ridge_effects has additive only + no dominance, so we fit an
# AUGMENTED ridge on [additive dosage | heterozygosity basis d*(ploidy-d)] to capture dominance
# (hence mid-parent heterosis / GPMP). Prediction is simulation-based: simulate progeny per
# candidate cross and score them under the ESTIMATED model -> predicted mean (GPMP incl.
# heterosis) + within-family variance -> usefulness (the variance-aware criterion).
# (its GPMP contribution); our usefulness allocates on mean + i*SD; ng_mean is the mean control.

aug_design <- function(dosage, ploidy) cbind(dosage, dosage * (ploidy - dosage))

fit_ghat <- function(dosage, pheno, ploidy, seed) {
  X <- aug_design(dosage, ploidy)
  colnames(X) <- paste0("f", seq_len(ncol(X))); rownames(X) <- rownames(dosage)
  names(pheno) <- rownames(dosage)
  fit <- ng_fit_ridge_effects(X, pheno, seed = as.integer(seed))
  fit$ploidy <- ploidy; fit
}
predict_val <- function(fit, dosage) as.numeric(fit$intercept + aug_design(dosage, fit$ploidy) %*% fit$beta)

add_pheno_noise <- function(gv, h2, seed) {
  set.seed(as.integer(seed)); vg <- stats::var(gv)
  ve <- if (is.finite(vg) && vg > 0) vg * (1 - h2) / h2 else 1
  gv + stats::rnorm(length(gv), 0, sqrt(ve))
}

# Predicted cross scores -- ANALYTIC (no per-cross progeny simulation), any ploidy. The augmented
# ridge fit carries additive effects (first m) on dosage and dominance effects (next m) on
# H = d*(ploidy-d); the package progeny-moment table + C++ kernel give the predicted cross mean and
# within-family variance directly (E[X], Var(X), E[H], Var(H), Cov(X,H) per parental dosage pair),
# so scoring is O(crosses x markers) instead of simulating cfg$n_pred progeny per cross.
score_crosses_pred <- function(parent_pop, SP, ploidy, fit, cfg, seed) {
  ids <- as.character(parent_pop@id)
  pdos <- AlphaSimR::pullSnpGeno(parent_pop, simParam = SP); rownames(pdos) <- ids
  storage.mode(pdos) <- "integer"
  K <- genomic_G(pdos) / ploidy
  intensity <- ng_selection_intensity(cfg$selection_prop)
  m <- ncol(pdos); ba <- fit$beta[seq_len(m)]; bd <- fit$beta[(m + 1L):(2L * m)]
  mt <- ng_polyploid_progeny_moment_table(ploidy)
  prs <- t(utils::combn(length(ids), 2L))
  i1 <- as.integer(prs[, 1] - 1L); i2 <- as.integer(prs[, 2] - 1L); zero <- numeric(m)
  if (exists("ng_poly_dominance_scores_cpp", mode = "function", inherits = TRUE)) {
    res <- ng_poly_dominance_scores_cpp(pdos, i1, i2, mt$mu, mt$varX, mt$EH, mt$varH, mt$covXH,
                                        ba, bd, zero, zero, fit$intercept, TRUE)  # zero centering: intercept absorbs it
    pred_mean <- res[, 1] + res[, 2]
    pred_var <- pmax(res[, 3] + res[, 4] + res[, 5], 0)
  } else {
    np <- nrow(prs); pred_mean <- pred_var <- numeric(np)
    for (i in seq_len(np)) {
      ij <- cbind(pdos[prs[i, 1], ] + 1L, pdos[prs[i, 2], ] + 1L)
      pred_mean[i] <- fit$intercept + sum(ba * mt$mu[ij]) + sum(bd * mt$EH[ij])
      pred_var[i] <- max(sum(ba^2 * mt$varX[ij]) + sum(bd^2 * mt$varH[ij]) + 2 * sum(ba * bd * mt$covXH[ij]), 0)
    }
  }
  a <- ids[prs[, 1]]; b <- ids[prs[, 2]]
  scores <- data.frame(parent1 = a, parent2 = b, pred_mean = pred_mean, pred_var = pred_var,
                       pred_usefulness = pred_mean + intensity * sqrt(pred_var),
                       pair_kinship = K[cbind(a, b)], stringsAsFactors = FALSE)
  attr(scores, "parent_kinship") <- K; scores
}

# Native allocators (ng_mean / ng_useful = greedy allocation on predicted mean / usefulness).
# In VARIABLE-family mode (cfg$family_mode == "variable") the plan carries an n_progeny column
# summing to cfg$total_progeny -- the realistic breeding decision (family size is never equal):
# select crosses, then allocate family sizes by score (ng_allocate_family_sizes).
allocate_plan <- function(method, scores, K, cfg, ploidy) {
  var_mode <- identical(cfg$family_mode, "variable")
  gain_col <- if (method == "ng_useful") "pred_usefulness" else "pred_mean"
  plan <- ng_optimize_mating_plan(scores, cfg$n_crosses, gain_col = gain_col, parent_kinship = K,
                                  max_crosses_per_parent = cfg$max_use, target_coancestry = cfg$dF,
                                  method = "greedy_local")
  if (var_mode) {
    plan <- ng_allocate_family_sizes(plan, total_progeny = cfg$total_progeny,
                                     method = "score_weighted", value_col = gain_col)
  }
  plan
}

# One recurrent breeding lineage for one method, from a shared founder setup.
recurrent_lineage <- function(method, setup, ploidy, cfg, rep_id, rep_seed) {
  SP <- setup$SP
  set.seed(as.integer(rep_seed))
  pop <- AlphaSimR::newPop(setup$founder, simParam = SP)
  if (AlphaSimR::nInd(pop) > cfg$pop_size) pop <- pop[seq_len(cfg$pop_size)]
  pop@id <- sprintf("I%04d", seq_len(AlphaSimR::nInd(pop)))
  base_gv <- mean(AlphaSimR::gv(pop)[, 1])
  rows <- list()
  for (cyc in seq_len(cfg$n_cycles)) {
    dosage <- AlphaSimR::pullSnpGeno(pop, simParam = SP); rownames(dosage) <- pop@id
    true_gv <- AlphaSimR::gv(pop)[, 1]
    pheno <- add_pheno_noise(true_gv, cfg$h2, rep_seed + cyc)
    fit <- fit_ghat(dosage, pheno, ploidy, rep_seed + cyc)
    gebv <- predict_val(fit, dosage)
    # common GEBV parent selection (same rule every method) -> isolate the allocation effect
    parents <- pop[order(gebv, decreasing = TRUE)[seq_len(cfg$n_parents)]]
    parents@id <- sprintf("P%03d", seq_len(cfg$n_parents))
    scores <- score_crosses_pred(parents, SP, ploidy, fit, cfg, as.integer(rep_seed %% 20011L) * 13L + cyc)
    K <- attr(scores, "parent_kinship")
    plan <- tryCatch(allocate_plan(method, scores, K, cfg, ploidy),
                     error = function(e) { cat(sprintf("  [%s p%d rep%d cyc%d] alloc ERR: %s\n", method, ploidy, rep_id, cyc, conditionMessage(e))); NULL })
    if (is.null(plan) || nrow(plan) < 1L) return(if (length(rows)) do.call(rbind, rows) else NULL)
    idx1 <- match(plan$parent1, parents@id); idx2 <- match(plan$parent2, parents@id)
    ok <- is.finite(idx1) & is.finite(idx2)
    idx1 <- idx1[ok]; idx2 <- idx2[ok]
    if (identical(cfg$family_mode, "variable") && "n_progeny" %in% names(plan)) {
      # VARIABLE family sizes: repeat each mating n_progeny times (nProgeny = 1) so the realized
      # families are unequal, as in a real breeding program. Total progeny = cfg$total_progeny.
      np <- pmax(as.integer(plan$n_progeny[ok]), 0L)
      reps <- rep(seq_along(idx1), np)
      cross_plan <- cbind(idx1[reps], idx2[reps])
      newpop <- AlphaSimR::makeCross(parents, cross_plan, nProgeny = 1L, simParam = SP)
    } else {
      newpop <- AlphaSimR::makeCross(parents, cbind(idx1, idx2),
                                     nProgeny = cfg$family_size, simParam = SP)
    }
    newpop@id <- sprintf("C%d_%04d", cyc, seq_len(AlphaSimR::nInd(newpop)))
    tg <- AlphaSimR::gv(newpop)[, 1]
    sel_ids <- unique(c(plan$parent1, plan$parent2))
    Ksel <- K[sel_ids, sel_ids, drop = FALSE]
    # selection-relevant accuracy = cor(GEBV, TRUE gv) on the current population (what drives
    # both parent selection and cross scoring); NOT fit$reliability, which is CV-r2 vs the
    # noisy phenotype and is capped near h2.
    pred_acc <- suppressWarnings(stats::cor(gebv, true_gv))
    rows[[cyc]] <- data.frame(
      rep = rep_id, method = method, ploidy = ploidy, cycle = cyc,
      gain = mean(tg) - base_gv, pop_mean_true = mean(tg), pop_var_true = stats::var(tg),
      sel_group_coancestry = group_coancestry(plan, K),
      mean_parent_rel = if (length(sel_ids) > 1L) mean(Ksel[upper.tri(Ksel)]) else NA_real_,
      pred_accuracy = if (is.finite(pred_acc)) pred_acc else NA_real_, stringsAsFactors = FALSE)
    if (AlphaSimR::nInd(newpop) > cfg$pop_size) newpop <- newpop[seq_len(cfg$pop_size)]
    pop <- newpop
  }
  do.call(rbind, rows)
}

run_recurrent <- function(cfg, seed, ploidies, methods) {
  cat(sprintf("Native mate-allocation RECURRENT study | reps=%d cycles=%d pop=%d parents=%d portfolio<=%d h2=%.2f | family_mode=%s (total_progeny=%d)\n",
              cfg$n_reps, cfg$n_cycles, cfg$pop_size, cfg$n_parents, cfg$n_crosses, cfg$h2,
              cfg$family_mode, cfg$total_progeny))
  all <- list()
  for (ploidy in ploidies) for (rep in seq_len(cfg$n_reps)) {
    setup <- sim_setup(ploidy, cfg, seed + rep)          # shared founders per (ploidy, rep)
    rep_seed <- as.integer((seed + rep * 1000L) %% 2147483647L)
    for (m in methods) {
      lin <- recurrent_lineage(m, setup, ploidy, cfg, rep, rep_seed)
      if (!is.null(lin)) all[[length(all) + 1L]] <- lin
    }
  }
  res <- do.call(rbind, all)
  # summary: final-cycle gain by method x ploidy (mean +/- se over reps)
  fin <- res[res$cycle == cfg$n_cycles, , drop = FALSE]
  agg <- aggregate(cbind(gain, pop_var_true, sel_group_coancestry, pred_accuracy) ~ method + ploidy,
                   data = fin, FUN = function(x) c(mean = mean(x), se = stats::sd(x) / sqrt(length(x))))
  cat("\n==== FINAL-CYCLE realized genetic gain (true values), mean over reps ====\n")
  fmt <- data.frame(ploidy = agg$ploidy, method = agg$method,
                    gain = sprintf("%.3f +/- %.3f", agg$gain[, "mean"], agg$gain[, "se"]),
                    var_true = sprintf("%.3f", agg$pop_var_true[, "mean"]),
                    coancestry = sprintf("%.4f", agg$sel_group_coancestry[, "mean"]),
                    gebv_accuracy = sprintf("%.2f", agg$pred_accuracy[, "mean"]),
                    stringsAsFactors = FALSE)
  fmt <- fmt[order(fmt$ploidy, fmt$method), ]
  print(fmt, row.names = FALSE)
  res
}

cfg <- list(
  founder_engine = tolower(env_chr("NG_MATEALLOC_FOUNDER", "runmacs")),
  n_founders = env_int("NG_MATEALLOC_N_FOUNDERS", 240L), n_chr = env_int("NG_MATEALLOC_N_CHR", 6L),
  seg_sites = env_int("NG_MATEALLOC_SEG_SITES", 400L), qtl_per_chr = env_int("NG_MATEALLOC_QTL_PER_CHR", 18L),
  snp_per_chr = env_int("NG_MATEALLOC_SNP_PER_CHR", 80L), n_parents = env_int("NG_MATEALLOC_N_PARENTS", 14L),
  n_crosses = env_int("NG_MATEALLOC_N_CROSSES", 12L), max_use = env_int("NG_MATEALLOC_MAX_USE", 4L),
  n_score_progeny = env_int("NG_MATEALLOC_N_SCORE_PROGENY", 40L),
  selection_prop = env_num("NG_MATEALLOC_SELECTION_PROP", 0.10), dF = env_num("NG_MATEALLOC_DF", 0.05),
  mean_dd = env_num("NG_MATEALLOC_MEAN_DD", 0.5), var_dd = env_num("NG_MATEALLOC_VAR_DD", 0.2),
  # recurrent-mode fields
  n_reps = env_int("NG_MATEALLOC_N_REPS", 3L), n_cycles = env_int("NG_MATEALLOC_N_CYCLES", 5L),
  pop_size = env_int("NG_MATEALLOC_POP_SIZE", 200L), family_size = env_int("NG_MATEALLOC_FAMILY_SIZE", 20L),
  n_pred = env_int("NG_MATEALLOC_N_PRED", 25L), h2 = env_num("NG_MATEALLOC_H2", 0.4),
  # VARIABLE family size (default) is the biological reality; "equal" is the old simulation
  # convenience. In variable mode n_crosses is a portfolio-size CAP and total_progeny is the seed
  # budget split across matings by score-weighted family sizes.
  family_mode = tolower(env_chr("NG_MATEALLOC_FAMILY_MODE", "variable")),
  total_progeny = env_int("NG_MATEALLOC_TOTAL_PROGENY", 0L)
)
if (cfg$total_progeny <= 0L) cfg$total_progeny <- cfg$pop_size
seed <- env_int("NG_MATEALLOC_SEED", 20260703L)
ploidies <- env_int_vec("NG_MATEALLOC_PLOIDY", c(2L, 4L))
mode <- tolower(env_chr("NG_MATEALLOC_MODE", "recurrent"))
methods <- strsplit(env_chr("NG_MATEALLOC_METHODS", "ng_mean,ng_useful"), "[,;[:space:]]+")[[1]]
methods <- methods[nzchar(methods)]
Sys.setenv(NG_ALPHASIMR_THREADS = "1")
out_dir <- env_chr("NG_MATEALLOC_OUTPUT_DIR", file.path(root, "results"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (identical(mode, "pilot")) {
  cat("Mate-allocation PILOT (single generation, true values).\n")
  res <- do.call(rbind, lapply(ploidies, function(p) run_ploidy(as.integer(p), cfg, seed)))
  out_file <- file.path(out_dir, "mate_allocation_pilot.csv")
} else {
  res <- run_recurrent(cfg, seed, ploidies, methods)
  out_file <- file.path(out_dir, "mate_allocation_recurrent.csv")
}
utils::write.csv(res, out_file, row.names = FALSE)
cat("\nWrote", out_file, "\n")
