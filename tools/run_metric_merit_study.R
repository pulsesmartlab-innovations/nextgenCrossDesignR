#!/usr/bin/env Rscript
# Metric-merit study: which cross-scoring metric earns its place, and under WHICH
# conditions. Races the package's scoring metrics on a controlled AlphaSimR grid and
# grades them on an UNBIASED true-effect yardstick, so any recommendation to end-users is
# evidence-based rather than assumed. This is a research/validation runner (like
# tools/run_ril_breeding_program_benchmark.R), not part of the installed package.
#
# Metrics raced (as cross-selection criteria):
#   mean           - mid-parent GEBV (gain only)
#   uc_vpm         - usefulness with the recombination variance      (mean + i*sqrt(VPM))
#   uc_pmv         - usefulness with the posterior-mean variance      (mean + i*sqrt(PMV)) [default]
#   usefulness_le  - usefulness with the naive relationship-distance "variance"
#   simplemating   - SimpleMating additive usefulness (external; skipped if unavailable)
#
# Conditions swept (the ones that change the answer): training-set size (effect-estimation
# quality), heritability, and trait architecture (oligogenic vs polygenic).
#
# Yardstick: every candidate cross is re-scored with the TRUE QTL effects (true usefulness
# = mid-parent true GV + i * true within-family SD), so a metric cannot flatter itself.
# Two allocation-free measures (metric quality, not allocator quality):
#   (A) ranking accuracy   = Spearman cor(metric, true usefulness) over ALL candidate crosses
#   (B) top-K realized merit = mean true usefulness of the metric's top-K crosses
# Realistic marker-QTL LD comes from runMacs (quickHaplo gives ~0 LD and is unusable here).
#
# Env config (all optional):
#   NG_MERIT_REPS(15) NG_MERIT_WORKERS(0=auto) NG_MERIT_SEED(20260703) NG_MERIT_USE_CPP(1)
#   NG_MERIT_N_PARENTS(40) NG_MERIT_TOPK(20) NG_MERIT_SELECTION_PROP(0.10)
#   NG_MERIT_TRAIN_N("120,400") NG_MERIT_H2("0.25,0.60") NG_MERIT_QTL_PER_CHR("5,40")
#   NG_MERIT_N_CHR(5) NG_MERIT_SEG_SITES(1200) NG_MERIT_SNP_PER_CHR(200)
#   NG_MERIT_INCLUDE_SIMPLEMATING(1) NG_MERIT_OUTPUT_DIR(results/) NG_MERIT_OUTPUT_PREFIX(metric_merit)

find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    start,
    file.path(start, "nextgen_cross_design"),
    dirname(start),
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
  stop("Could not locate nextgenCrossDesign root containing DESCRIPTION and R/load.R", call. = FALSE)
}

env_chr <- function(name, default) { v <- Sys.getenv(name, unset = ""); if (!nzchar(v)) default else v }
env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = ""); if (!nzchar(v)) return(as.integer(default))
  o <- suppressWarnings(as.integer(v)); if (!is.finite(o)) as.integer(default) else o
}
env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = ""); if (!nzchar(v)) return(as.numeric(default))
  o <- suppressWarnings(as.numeric(v)); if (!is.finite(o)) as.numeric(default) else o
}
env_bool <- function(name, default = FALSE) {
  v <- tolower(trimws(Sys.getenv(name, unset = ""))); if (!nzchar(v)) return(isTRUE(default))
  v %in% c("1", "true", "yes", "y")
}
env_num_vec <- function(name, default) {
  v <- Sys.getenv(name, unset = ""); if (!nzchar(v)) return(default)
  o <- suppressWarnings(as.numeric(strsplit(v, "[,;[:space:]]+")[[1]])); o[is.finite(o)]
}

resolve_workers <- function(req, n) {
  if (.Platform$OS.type == "windows") return(1L)   # no forking; run serial
  if (req > 0L) return(min(req, n))
  avail <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  max(1L, min(n, avail - 2L))
}
parallel_lapply <- function(x, fun, workers) {
  if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = FALSE)
  } else lapply(x, fun)
}

metric_merit_config <- function(root = find_project_root()) {
  list(
    root = root,
    reps = env_int("NG_MERIT_REPS", 15L),
    workers = env_int("NG_MERIT_WORKERS", 0L),
    seed = env_int("NG_MERIT_SEED", 20260703L),
    use_cpp = env_bool("NG_MERIT_USE_CPP", TRUE),
    n_parents = env_int("NG_MERIT_N_PARENTS", 40L),
    topk = env_int("NG_MERIT_TOPK", 20L),
    selection_prop = env_num("NG_MERIT_SELECTION_PROP", 0.10),
    train_n = as.integer(env_num_vec("NG_MERIT_TRAIN_N", c(120, 400))),
    h2 = env_num_vec("NG_MERIT_H2", c(0.25, 0.60)),
    qtl_per_chr = as.integer(env_num_vec("NG_MERIT_QTL_PER_CHR", c(5, 40))),
    n_chr = env_int("NG_MERIT_N_CHR", 5L),
    seg_sites = env_int("NG_MERIT_SEG_SITES", 1200L),
    snp_per_chr = env_int("NG_MERIT_SNP_PER_CHR", 200L),
    include_simplemating = env_bool("NG_MERIT_INCLUDE_SIMPLEMATING", TRUE),
    output_dir = env_chr("NG_MERIT_OUTPUT_DIR", file.path(root, "results")),
    output_prefix = env_chr("NG_MERIT_OUTPUT_PREFIX", "metric_merit")
  )
}

arch_label <- function(q) if (q <= 10L) "oligogenic" else "polygenic"

# Score every metric for one population and grade against the true-effect yardstick.
score_metrics <- function(cfg, g_snp, mm_snp, g_qtl, mm_qtl, true_eff, train_g, train_y, cand_y, want_sm) {
  i_int <- ng_selection_intensity(cfg$selection_prop)
  names(true_eff) <- colnames(g_qtl)
  fit <- ng_fit_ridge_effects(train_g, train_y, ids = rownames(train_g), seed = 1L)
  sc  <- ng_score_crosses(g_snp, fit, marker_map = mm_snp, ids = rownames(g_snp),
                          adjusted_pheno = cand_y, target = "DH",
                          selection_prop = cfg$selection_prop, use_cpp = cfg$use_cpp)
  pairs <- sc[, c("parent1", "parent2")]
  # TRUE-effect yardstick (exact, no shrinkage): mid-parent true GV + i * true within-family SD.
  gv_true <- as.numeric(g_qtl %*% true_eff); names(gv_true) <- rownames(g_qtl)
  tmean <- 0.5 * (gv_true[pairs$parent1] + gv_true[pairs$parent2])
  tvar <- ng_dh_recomb_variance_pairs(g_qtl, beta = as.numeric(true_eff),
            beta_var = rep(0, length(true_eff)), marker_map = mm_qtl,
            ids = rownames(g_qtl), pairs = pairs, target = "DH", use_cpp = cfg$use_cpp)$vpm
  true_uc <- as.numeric(tmean) + i_int * sqrt(pmax(tvar, 0))

  metrics <- list(
    mean          = sc$cross_mean_gebv,
    uc_vpm        = sc$usefulness_vpm_gebv,
    uc_pmv        = sc$usefulness_pmv_gebv,
    usefulness_le = sc$cross_mean_gebv + i_int * sqrt(pmax(sc$parent_distance, 0))
  )
  if (isTRUE(want_sm)) {
    sm <- tryCatch(ng_add_simplemating_scores(sc, geno = g_snp, effects = fit, marker_map = mm_snp,
                     adjusted_pheno = cand_y, type = "DH", generation = 1L, engine = "external",
                     fallback = "none")$simple_usefa, error = function(e) rep(NA_real_, nrow(sc)))
    if (any(is.finite(sm))) metrics$simplemating <- sm
  }
  do.call(rbind, lapply(names(metrics), function(nm) {
    m <- metrics[[nm]]; ok <- is.finite(m) & is.finite(true_uc)
    if (sum(ok) < cfg$topk + 5L) return(NULL)
    topk <- order(m, decreasing = TRUE)[seq_len(cfg$topk)]
    data.frame(metric = nm,
               rank_acc = suppressWarnings(cor(m[ok], true_uc[ok], method = "spearman")),
               topk_merit = mean(true_uc[topk], na.rm = TRUE), stringsAsFactors = FALSE)
  }))
}

# One (condition x rep) draw: simulate, estimate effects, score, grade.
run_one <- function(cfg, conditions, idx) {
  ci <- ((idx - 1L) %% nrow(conditions)) + 1L
  rep_i <- ((idx - 1L) %/% nrow(conditions)) + 1L
  cond <- conditions[ci, ]
  tryCatch({
    set.seed(cfg$seed + idx)   # deterministic per idx regardless of worker count / OS
    fp <- AlphaSimR::runMacs(nInd = cond$train_n + cfg$n_parents + 20L, nChr = cfg$n_chr,
                             segSites = cfg$seg_sites, species = "GENERIC")
    SP <- AlphaSimR::SimParam$new(fp)
    SP$addTraitA(nQtlPerChr = cond$qtl_chr); SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
    SP$setVarE(h2 = cond$h2)
    dh <- AlphaSimR::makeDH(AlphaSimR::newPop(fp, simParam = SP), simParam = SP)
    ids <- paste0("L", seq_len(AlphaSimR::nInd(dh))); dh@id <- ids
    ph <- as.numeric(AlphaSimR::setPheno(dh, h2 = cond$h2, simParam = SP)@pheno[, 1]); names(ph) <- ids
    train_ids <- ids[seq_len(cond$train_n)]
    cand_ids  <- ids[cond$train_n + seq_len(cfg$n_parents)]
    snp <- AlphaSimR::pullSnpGeno(dh, simParam = SP); rownames(snp) <- ids; colnames(snp) <- paste0("S", seq_len(ncol(snp)))
    qtl <- AlphaSimR::pullQtlGeno(dh, simParam = SP); rownames(qtl) <- ids; colnames(qtl) <- paste0("Q", seq_len(ncol(qtl)))
    sm_map <- AlphaSimR::getSnpMap(simParam = SP); qm_map <- AlphaSimR::getQtlMap(trait = 1L, simParam = SP)
    mm_snp <- data.frame(marker = colnames(snp), chr = sm_map$chr, pos_cm = sm_map$pos * 100,
                         chr_index = as.integer(as.factor(sm_map$chr)))
    mm_qtl <- data.frame(marker = colnames(qtl), chr = qm_map$chr, pos_cm = qm_map$pos * 100,
                         chr_index = as.integer(as.factor(qm_map$chr)))
    true_eff <- as.numeric(SP$traits[[1]]@addEff)
    res <- score_metrics(cfg, snp[cand_ids, , drop = FALSE], mm_snp, qtl[cand_ids, , drop = FALSE],
                         mm_qtl, true_eff, snp[train_ids, , drop = FALSE], ph[train_ids], ph[cand_ids],
                         want_sm = cfg$include_simplemating)
    if (is.null(res)) return(NULL)
    res$train_n <- cond$train_n; res$h2 <- cond$h2; res$arch <- arch_label(cond$qtl_chr)
    res$rep <- rep_i; res
  }, error = function(e) { message("merit run ", idx, " failed: ", conditionMessage(e)); NULL })
}

run_metric_merit_study <- function(cfg = metric_merit_config()) {
  if (!nzchar(Sys.getenv("RGL_USE_NULL"))) Sys.setenv(RGL_USE_NULL = "TRUE")  # SimpleMating->rgl headless
  source(file.path(cfg$root, "R", "load.R"))
  ng_load(cfg$root, use_cpp = cfg$use_cpp, verbose = FALSE)
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("AlphaSimR is required for the metric-merit study; install it or run the RIL benchmark's dependency setup.", call. = FALSE)
  }
  want_sm <- cfg$include_simplemating && requireNamespace("SimpleMating", quietly = TRUE)
  if (cfg$include_simplemating && !want_sm) {
    message("SimpleMating not available; racing package metrics only.")
  }
  cfg$include_simplemating <- want_sm

  conditions <- expand.grid(train_n = cfg$train_n, h2 = cfg$h2, qtl_chr = cfg$qtl_per_chr,
                            stringsAsFactors = FALSE)
  n_runs <- nrow(conditions) * cfg$reps
  workers <- resolve_workers(cfg$workers, n_runs)
  message(sprintf("Metric-merit study: %d conditions x %d reps = %d runs on %d worker(s)",
                  nrow(conditions), cfg$reps, n_runs, workers))

  res_list <- parallel_lapply(seq_len(n_runs), function(i) run_one(cfg, conditions, i), workers)
  res <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
  if (is.null(res) || !nrow(res)) stop("All metric-merit runs failed; check AlphaSimR install and logs.", call. = FALSE)

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  raw_csv <- file.path(cfg$output_dir, paste0(cfg$output_prefix, "_raw.csv"))
  utils::write.csv(res, raw_csv, row.names = FALSE)
  saveRDS(res, file.path(cfg$output_dir, paste0(cfg$output_prefix, "_raw.rds")))

  se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else stats::sd(x) / sqrt(length(x)) }
  agg <- aggregate(cbind(rank_acc, topk_merit) ~ metric + train_n + h2 + arch, res,
                   function(x) c(mean = mean(x), se = se(x)))
  summ <- data.frame(
    metric = agg$metric, train_n = agg$train_n, h2 = agg$h2, arch = agg$arch,
    rank_acc_mean = agg$rank_acc[, "mean"], rank_acc_se = agg$rank_acc[, "se"],
    topk_merit_mean = agg$topk_merit[, "mean"], topk_merit_se = agg$topk_merit[, "se"],
    stringsAsFactors = FALSE)
  summ <- summ[order(summ$arch, summ$train_n, summ$h2, -summ$rank_acc_mean), ]
  utils::write.csv(summ, file.path(cfg$output_dir, paste0(cfg$output_prefix, "_summary.csv")), row.names = FALSE)

  cat("\n=== RANKING ACCURACY (Spearman cor with true usefulness; mean +/- SE) by condition ===\n")
  for (a in unique(summ$arch)) for (tn in sort(unique(summ$train_n))) for (hh in sort(unique(summ$h2))) {
    sub <- summ[summ$arch == a & summ$train_n == tn & summ$h2 == hh, ]
    if (!nrow(sub)) next
    cat(sprintf("\n[ arch=%s | train_n=%d | h2=%.2f ]\n", a, tn, hh))
    print(data.frame(metric = sub$metric,
                     rank_acc = sprintf("%.3f+/-%.3f", sub$rank_acc_mean, sub$rank_acc_se),
                     topk_true_merit = sprintf("%.3f+/-%.3f", sub$topk_merit_mean, sub$topk_merit_se),
                     stringsAsFactors = FALSE), row.names = FALSE)
  }
  ov <- aggregate(rank_acc ~ metric, res, function(x) c(mean = mean(x), se = se(x)))
  ov <- ov[order(-ov$rank_acc[, "mean"]), ]
  cat("\n=== overall mean ranking accuracy across all conditions ===\n")
  for (i in seq_len(nrow(ov))) cat(sprintf("  %-14s %.3f +/- %.3f\n", ov$metric[i], ov$rank_acc[i, "mean"], ov$rank_acc[i, "se"]))

  message("\nWrote metric-merit outputs to: ", normalizePath(cfg$output_dir, winslash = "/", mustWork = FALSE))
  message("Output prefix: ", cfg$output_prefix)
  invisible(list(raw = res, summary = summ))
}

is_this_script <- function() any(grepl("run_metric_merit_study\\.R$", commandArgs(FALSE)))
if (is_this_script()) invisible(run_metric_merit_study())
