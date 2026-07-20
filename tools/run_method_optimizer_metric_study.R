#!/usr/bin/env Rscript
# =============================================================================
# Method x Optimizer x Metric recurrent RIL breeding-program study vs SimpleMating
# =============================================================================
# A breeder-practical, multi-cycle simulation (AlphaSimR = truth) that races
# METHOD x OPTIMIZER x METRIC combinations -- driving nextgenCrossDesign as a
# USER (ng_run_cross_prediction) -- against the external SimpleMating package.
#
#   * 25 cycles x 10 reps (env-configurable), run in PARALLEL over rep x arm.
#   * RIL pipeline: crosses -> selfing (single-seed descent) F2..F6, NOT DH.
#   * Founders are generated ONCE per rep and SHARED across all arms so every
#     method competes on identical genetics (paired comparison).
#   * SimpleMating: try the real external package; fall back to the package's
#     UC-VPM arm (evidence: metric-merit study shows UC-VPM == SimpleMating),
#     labelled honestly via implementation / exact_external_status columns.
#
# Per cycle we record the 6 breeder tracking metrics (gain + diversity depletion)
# added last session (cf. ril_metrics_row in run_ril_breeding_program_benchmark.R):
#   genetic_gain_index, var_index_gv, expected_heterozygosity,
#   prop_polymorphic_markers, mean_maf, mean_parent_relationship.
#
# >>> A NEW USER/BREEDER honest-performance comparison. EVOLUTION is the standard OCS
#     optimizer on every arm; the only non-evolution arms are the explicit OPTIMIZER
#     axis. The ARM_REGISTRY below is the single place to edit the arms:
#       - METRIC axis: the 6 real scoring metrics (mean/var_simple/uc-vpm/uc-pmv/pmv/
#         vpm/var_complex) on OCS + evolution
#       - STRATEGY dial (strategy_balanced), ALLOCATION axis (alphamate_style vs OCS)
#       - OPTIMIZER axis: uc/pmv under greedy / mip (vs uc_pmv_ocs = evolution)
#       - comparators: SimpleMating (simple_usefa_select) + random control
#
# GS TRAINING (realistic breeder setup): each cycle the GS model is trained on the
# PHENOTYPED candidate parents themselves via ng_run_cross_prediction -- the honest
# user-API path. There is deliberately NO synthetic 500/120-line training population
# (unrealistic for a real program); the breeder trains on the lines they actually
# phenotype and predicts crosses among the elite parents.
#
# Parameters (aligned to last session's recurrent RIL study; tractable defaults --
# scale up via env). All optional:
#   NG_MOM_REPS(10) NG_MOM_CYCLES(25) NG_MOM_WORKERS(0=auto) NG_MOM_SEED(20260703) NG_MOM_USE_CPP(1)
#   NG_MOM_N_FOUNDER(80) NG_MOM_N_PARENTS(30) NG_MOM_N_CROSSES(40)
#   NG_MOM_N_CHR(5) NG_MOM_QTL_PER_CHR(40) NG_MOM_SNP_PER_CHR(200) NG_MOM_SEG_SITES(500) NG_MOM_FOUNDER_H2(0.5)
#   NG_MOM_F2_PER_CROSS(80) NG_MOM_BULK_SIZE(5) NG_MOM_STAGE_H2(0.03,0.15,0.40,0.60,0.70)
#   NG_MOM_WITHIN_FAMILY_PROP(0.10) NG_MOM_BETWEEN_FAMILY_PROP(0.20) NG_MOM_AYT_PROP(0.70)
#   NG_MOM_SELECTION_PROP(0.10) NG_MOM_METHOD_VARPMV(fast) NG_MOM_MIN_EFFECT_RELIABILITY(0.20)
#   NG_MOM_RECOMBINATION_MODEL(haldane) NG_MOM_MAX_USES(6)
#   NG_MOM_LAMBDA_GROUP(0.05) NG_MOM_LAMBDA_MATING(0.02)
#   NG_MOM_EVOL_SOLUTIONS(100) NG_MOM_EVOL_ITERATIONS(150) NG_MOM_EVOL_STOP(30)
#   NG_MOM_INCLUDE_SIMPLEMATING(1) NG_MOM_SIMPLEMATING_GITHUB(1) NG_MOM_SIMPLEMATING_CULLING_K(0.5)
#   NG_MOM_OUTPUT_DIR(results/) NG_MOM_OUTPUT_PREFIX(method_optimizer_metric)
# =============================================================================

suppressWarnings(suppressMessages({
  Sys.setenv(RGL_USE_NULL = "TRUE")   # SimpleMating -> rgl headless
}))

# ---- project root + package load (as a user) --------------------------------
find_project_root <- function(start = getwd()) {
  cands <- unique(normalizePath(c(
    start, file.path(start, "nextgen_cross_design"),
    dirname(start), file.path(dirname(start), "nextgen_cross_design")
  ), winslash = "/", mustWork = FALSE))
  for (c in cands) {
    if (file.exists(file.path(c, "R", "load.R")) && file.exists(file.path(c, "DESCRIPTION")) &&
        any(grepl("^Package:\\s*nextgenCrossDesign", readLines(file.path(c, "DESCRIPTION"), n = 5L, warn = FALSE))))
      return(c)
  }
  stop("Could not locate nextgenCrossDesign root", call. = FALSE)
}

env_chr  <- function(n, d) { v <- Sys.getenv(n, ""); if (!nzchar(v)) d else v }
env_int  <- function(n, d) { v <- Sys.getenv(n, ""); o <- suppressWarnings(as.integer(v)); if (!nzchar(v) || !is.finite(o)) as.integer(d) else o }
env_num  <- function(n, d) { v <- Sys.getenv(n, ""); o <- suppressWarnings(as.numeric(v)); if (!nzchar(v) || !is.finite(o)) as.numeric(d) else o }
env_bool <- function(n, d) { v <- tolower(trimws(Sys.getenv(n, ""))); if (!nzchar(v)) isTRUE(d) else v %in% c("1","true","yes","y") }
env_num_vec <- function(n, d) { v <- Sys.getenv(n, ""); if (!nzchar(v)) return(d); o <- suppressWarnings(as.numeric(strsplit(v, "[,;[:space:]]+")[[1]])); o[is.finite(o)] }

mom_config <- function(root = find_project_root()) {
  list(
    root = root,
    reps = env_int("NG_MOM_REPS", 10L),
    cycles = env_int("NG_MOM_CYCLES", 25L),
    workers = env_int("NG_MOM_WORKERS", 0L),
    seed = env_int("NG_MOM_SEED", 20260703L),
    use_cpp = env_bool("NG_MOM_USE_CPP", TRUE),
    # --- genetic architecture / population (aligned to last session's RIL study;
    #     tractable defaults -- scale to the "realistic" run with the env vars) ---
    n_founder = env_int("NG_MOM_N_FOUNDER", 80L),
    n_parents = env_int("NG_MOM_N_PARENTS", 30L),
    n_crosses = env_int("NG_MOM_N_CROSSES", 40L),
    n_chr = env_int("NG_MOM_N_CHR", 5L),
    qtl_per_chr = env_int("NG_MOM_QTL_PER_CHR", 40L),
    snp_per_chr = env_int("NG_MOM_SNP_PER_CHR", 200L),
    seg_sites = env_int("NG_MOM_SEG_SITES", 500L),
    founder_h2 = env_num("NG_MOM_FOUNDER_H2", 0.5),
    # --- breeding-pipeline structure (F2->F6 selfing/SSD, staged selection) ------
    f2_per_cross = env_int("NG_MOM_F2_PER_CROSS", 80L),
    bulk_size = env_int("NG_MOM_BULK_SIZE", 5L),
    within_family_prop = env_num("NG_MOM_WITHIN_FAMILY_PROP", 0.10),
    between_family_prop = env_num("NG_MOM_BETWEEN_FAMILY_PROP", 0.20),
    ayt_prop = env_num("NG_MOM_AYT_PROP", 0.70),
    h2_stages = env_num_vec("NG_MOM_STAGE_H2", c(0.03, 0.15, 0.40, 0.60, 0.70)),  # F2..F6
    # --- scoring / allocation parameters used in the study ----------------------
    selection_prop = env_num("NG_MOM_SELECTION_PROP", 0.10),
    method_varpmv = env_chr("NG_MOM_METHOD_VARPMV", "fast"),
    min_effect_reliability = env_num("NG_MOM_MIN_EFFECT_RELIABILITY", 0.20),
    recomb_model = env_chr("NG_MOM_RECOMBINATION_MODEL", "haldane"),
    max_uses = env_int("NG_MOM_MAX_USES", 6L),
    lambda_group = env_num("NG_MOM_LAMBDA_GROUP", 0.05),
    lambda_mating = env_num("NG_MOM_LAMBDA_MATING", 0.02),
    evol_solutions = env_int("NG_MOM_EVOL_SOLUTIONS", 100L),
    evol_iterations = env_int("NG_MOM_EVOL_ITERATIONS", 150L),
    evol_stop = env_int("NG_MOM_EVOL_STOP", 30L),
    # --- comparators / output ---------------------------------------------------
    include_simplemating = env_bool("NG_MOM_INCLUDE_SIMPLEMATING", TRUE),
    simplemating_github = env_bool("NG_MOM_SIMPLEMATING_GITHUB", TRUE),
    simplemating_culling_k = env_num("NG_MOM_SIMPLEMATING_CULLING_K", 0.5),
    output_dir = env_chr("NG_MOM_OUTPUT_DIR", file.path(root, "results")),
    output_prefix = env_chr("NG_MOM_OUTPUT_PREFIX", "method_optimizer_metric")
  )
}

# ---- ARM REGISTRY -----------------------------------------------------------
# One row = one competing arm. engine: "package" (ng_run_cross_prediction),
# "random", or "simplemating". For package arms the (trait_value_metric,
# uc_variance_source, allocation_method, optimizer, strategy) tuple defines it.
# Edit this single table to change which method/optimizer/metric combos race.
arm <- function(arm, engine, metric = NA, varsrc = NA, alloc = NA, opt = NA, strategy = NA) {
  data.frame(arm = arm, engine = engine, trait_value_metric = metric,
             uc_variance_source = varsrc, allocation_method = alloc,
             optimizer = opt, strategy = strategy, stringsAsFactors = FALSE)
}

arm_registry <- function() {
  rbind(
    # --- METRIC axis: the 6 real scoring metrics, all on OCS + the EVOLUTION optimizer ---
    arm("mean_ocs",             "package", "mean",        "pmv", "ocs", "evolution"),
    arm("var_simple_ocs",       "package", "var_simple",  "pmv", "ocs", "evolution"),
    arm("uc_vpm_ocs",           "package", "usefulness",          "vpm", "ocs", "evolution"),
    arm("uc_pmv_ocs",           "package", "usefulness",          "pmv", "ocs", "evolution"),
    arm("pmv_ocs",              "package", "pmv",         "pmv", "ocs", "evolution"),
    arm("vpm_ocs",              "package", "vpm",         "pmv", "ocs", "evolution"),
    arm("var_complex_ocs",      "package", "var_complex", "pmv", "ocs", "evolution"),
    # --- STRATEGY dial: explicit diversity management on a good merit metric (evolution) --
    arm("strategy_balanced",    "package", "var_complex", "pmv", "ocs", "evolution", strategy = "balanced"),
    # --- ALLOCATION axis: native AlphaMate-style vs OCS (metric + optimizer held fixed) ---
    arm("var_complex_alphamate","package", "var_complex", "pmv", "alphamate_style", "evolution"),
    # --- OPTIMIZER axis: same uc/pmv metric under greedy / mip (vs uc_pmv_ocs = evolution) -
    arm("uc_pmv_greedy",        "package", "usefulness",          "pmv", "ocs", "greedy_local"),
    arm("uc_pmv_mip",           "package", "usefulness",          "pmv", "ocs", "mip"),
    # --- comparators ----------------------------------------------------------
    arm("simple_usefa_select",  "simplemating"),
    arm("random",               "random")
  )
}

# ---- founders (mirror AlphaSim_2.r, per rep) --------------------------------
build_founders <- function(cfg) {
  fp <- AlphaSimR::runMacs(nInd = cfg$n_founder, nChr = cfg$n_chr,
                           segSites = cfg$seg_sites, inbred = TRUE, split = 100)
  SP <- AlphaSimR::SimParam$new(fp)
  SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1)
  SP$setVarE(h2 = cfg$founder_h2)
  SP$addSnpChip(cfg$snp_per_chr)
  founders <- AlphaSimR::newPop(fp, simParam = SP)
  founders <- AlphaSimR::setPheno(founders, h2 = cfg$founder_h2, simParam = SP)
  founders@id <- paste0("F", seq_len(AlphaSimR::nInd(founders)))
  list(founders = founders, SP = SP)
}

# ---- the 6 breeder tracking metrics (+ raw mean) ----------------------------
breeder_metrics <- function(pop, SP, base_mean) {
  g <- AlphaSimR::gv(pop)[, 1]
  geno <- AlphaSimR::pullSnpGeno(pop, simParam = SP)
  rownames(geno) <- pop@id
  p <- colMeans(geno, na.rm = TRUE) / 2
  poly <- p > 0 & p < 1
  K <- ng_parent_kinship(geno)                     # package (VanRaden G), as a user
  offd <- K[upper.tri(K)]
  data.frame(
    mean_gv = mean(g),
    genetic_gain_index = mean(g) - base_mean,
    var_index_gv = stats::var(g),
    expected_heterozygosity = mean(2 * p * (1 - p), na.rm = TRUE),
    prop_polymorphic_markers = mean(poly),
    mean_maf = mean(pmin(p, 1 - p), na.rm = TRUE),
    mean_parent_relationship = if (length(offd)) mean(offd, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# RIL parents carry a few unfixed (heterozygous) loci, but both the inbred cross-variance
# kernel AND the runner QC require fully homozygous integer 0/2 calls (NA is rejected as a
# non-finite dosage). Impute each residual-het call to that marker's MAJOR homozygote -> the
# parent reads as a fixed RIL. Deterministic, integer, low approximation (~few % of loci).
impute_het_to_hom <- function(geno) {
  het <- which(geno == 1L, arr.ind = TRUE)
  if (nrow(het)) {
    major <- ifelse(colMeans(geno == 2L, na.rm = TRUE) >= colMeans(geno == 0L, na.rm = TRUE), 2L, 0L)
    geno[het] <- major[het[, "col"]]
  }
  storage.mode(geno) <- "integer"
  geno
}

# ---- cross selection dispatch (package as a user / SimpleMating / random) ----
select_crosses_package <- function(arm, parents, SP, cfg, seed) {
  ids <- parents@id
  geno <- impute_het_to_hom(AlphaSimR::pullSnpGeno(parents, simParam = SP)); rownames(geno) <- ids
  ph <- as.numeric(parents@pheno[, 1]); names(ph) <- ids
  smap <- AlphaSimR::getSnpMap(simParam = SP)         # id, chr, site, pos(Morgan)
  genotype  <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
  phenotype <- data.frame(NAME = ids, yield = ph, stringsAsFactors = FALSE)
  marker_map <- data.frame(SNP = colnames(geno), Chr = smap$chr,
                           PosBP = smap$pos * 100 * 1e6,  # Morgan->cM->"bp" (÷1e6 recovers cM)
                           stringsAsFactors = FALSE)
  direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                          stringsAsFactors = FALSE)
  call_args <- list(
    phenotype = phenotype, genotype = genotype, marker_map = marker_map,
    trait_direction = direction, id_col = "NAME",
    map_marker_col = "SNP", map_chr_col = "Chr", map_pos_col = "PosBP",
    map_pos_cm_divisor = 1e6,
    prediction_mode = "trait_by_trait",
    trait_value_metric = arm$trait_value_metric,
    uc_variance_source = arm$uc_variance_source,
    method_varPMV = cfg$method_varpmv,
    progeny = "RIL", ril_mode = "infinite",
    selection_prop = cfg$selection_prop,
    min_effect_reliability = cfg$min_effect_reliability,
    recomb_model = cfg$recomb_model,
    duplicate_action = "none",
    n_crosses = cfg$n_crosses, max_crosses_per_parent = cfg$max_uses,
    allocation_method = arm$allocation_method, optimizer = arm$optimizer,
    use_ocs = TRUE,
    lambda_group = cfg$lambda_group, lambda_mating = cfg$lambda_mating,
    evol_solutions = cfg$evol_solutions, evol_iterations = cfg$evol_iterations,
    evol_stop = cfg$evol_stop,
    write_outputs = FALSE, write_figures = FALSE, seed = seed
  )
  if (!is.na(arm$strategy)) call_args$strategy <- arm$strategy   # strategy dial (e.g. "balanced")
  res <- do.call(ng_run_cross_prediction, call_args)
  sel <- res$selected_crosses
  as.matrix(sel[, c("parent1", "parent2")])
}

select_crosses_random <- function(parents, cfg) {
  ids <- parents@id
  grid <- expand.grid(ids, ids, stringsAsFactors = FALSE)
  grid <- grid[grid$Var1 != grid$Var2, ]
  # de-duplicate reciprocals
  key <- apply(grid, 1L, function(r) paste(sort(r), collapse = "_"))
  grid <- grid[!duplicated(key), ]
  n <- min(cfg$n_crosses, nrow(grid))
  as.matrix(grid[sample.int(nrow(grid), n), ])
}

select_crosses_simplemating <- function(parents, SP, cfg, seed, sm_external) {
  ids <- parents@id
  if (!sm_external) {
    # native evidence-based proxy: UC-VPM (== SimpleMating per metric-merit study)
    proxy <- data.frame(arm = "simplemating_native", engine = "package",
                        trait_value_metric = "usefulness", uc_variance_source = "vpm",
                        allocation_method = "ocs", optimizer = "evolution",
                        strategy = NA, stringsAsFactors = FALSE)
    return(select_crosses_package(proxy, parents, SP, cfg, seed))
  }
  # External SimpleMating via the package's validated wrappers (getUsefA usefulness with
  # residual-het handling; selectCrosses with relatedness culling). RIL type, generation 6.
  geno <- impute_het_to_hom(AlphaSimR::pullSnpGeno(parents, simParam = SP)); rownames(geno) <- ids
  ph <- as.numeric(parents@pheno[, 1]); names(ph) <- ids
  smap <- AlphaSimR::getSnpMap(simParam = SP)
  marker_map <- data.frame(marker = colnames(geno), chr = smap$chr, pos_cm = smap$pos * 100,
                           stringsAsFactors = FALSE)
  effects <- ng_fit_ridge_effects(geno, ph, ids = ids, seed = seed)
  sc <- ng_score_crosses(geno, effects, marker_map = marker_map, ids = ids,
                         adjusted_pheno = ph, target = "RIL",
                         selection_prop = cfg$selection_prop, use_cpp = cfg$use_cpp)
  sc <- ng_add_simplemating_scores(sc, geno = geno, effects = effects, marker_map = marker_map,
                                   adjusted_pheno = ph, type = "RIL", generation = 6L,
                                   prop_sel = cfg$selection_prop, engine = "external", fallback = "none")
  plan <- ng_select_simplemating(sc, score_col = "simple_usefa", n_crosses = cfg$n_crosses,
                                 parent_kinship = ng_parent_kinship(geno),
                                 max_crosses_per_parent = cfg$max_uses, min_crosses_per_parent = 1L,
                                 culling_pairwise_k = cfg$simplemating_culling_k)
  as.matrix(plan[, c("parent1", "parent2")])
}

# ---- RIL advancement: selfing (SSD) F2..F6 with staged selection ------------
# Mirrors the breeder-practical scheme of AlphaSim_2.r / run_ril_breeding_program:
# F2 DH-free family -> within-family select -> F3/F4 bulks -> PYT/AYT/EYT, then
# recycle the best n_parents. All advancement is by selfing (RIL/SSD), never makeDH.
advance_ril <- function(F1, SP, cfg) {
  h2 <- cfg$h2_stages
  npar <- cfg$n_parents
  F2 <- AlphaSimR::setPheno(AlphaSimR::self(F1, nProgeny = cfg$f2_per_cross, simParam = SP),
                            h2 = h2[1], simParam = SP)
  k2 <- max(1L, round(cfg$f2_per_cross * cfg$within_family_prop))
  F2s <- AlphaSimR::selectWithinFam(F2, nInd = k2, use = "pheno", simParam = SP)
  F3 <- AlphaSimR::setPheno(AlphaSimR::self(F2s, nProgeny = cfg$bulk_size, simParam = SP),
                            h2 = h2[2], simParam = SP)
  F3s <- AlphaSimR::selectWithinFam(F3, nInd = 1, use = "pheno", simParam = SP)   # SSD: 1/family
  F4 <- AlphaSimR::setPheno(AlphaSimR::self(F3s, nProgeny = cfg$bulk_size, simParam = SP),
                            h2 = h2[3], simParam = SP)
  n4 <- max(npar, round(AlphaSimR::nInd(F4) * cfg$between_family_prop))
  F4s <- AlphaSimR::selectInd(F4, nInd = min(n4, AlphaSimR::nInd(F4)), use = "pheno", simParam = SP)  # PYT
  F5 <- AlphaSimR::setPheno(AlphaSimR::self(F4s, nProgeny = 1, simParam = SP), h2 = h2[4], simParam = SP)
  n5 <- max(npar, round(AlphaSimR::nInd(F5) * cfg$ayt_prop))
  F5s <- AlphaSimR::selectInd(F5, nInd = min(n5, AlphaSimR::nInd(F5)), use = "pheno", simParam = SP)  # AYT
  F6 <- AlphaSimR::setPheno(AlphaSimR::self(F5s, nProgeny = 1, simParam = SP), h2 = h2[5], simParam = SP)  # EYT
  next_parents <- AlphaSimR::selectInd(F6, nInd = min(npar, AlphaSimR::nInd(F6)),
                                       use = "pheno", simParam = SP)
  list(cohort = F6, next_parents = next_parents)
}

# ---- one arm x one rep: 25-cycle recurrent RIL program ----------------------
run_arm_rep <- function(cfg, arm, rep_i, arm_idx, founder_obj, sm_external) {
  set.seed(cfg$seed + rep_i * 1000L + arm_idx)
  SP <- founder_obj$SP
  parents <- AlphaSimR::selectInd(founder_obj$founders, nInd = cfg$n_parents,
                                  trait = 1, use = "pheno", simParam = SP)
  base_mean <- mean(AlphaSimR::gv(parents)[, 1])
  rows <- vector("list", cfg$cycles)
  for (cyc in seq_len(cfg$cycles)) {
    crossPlan <- tryCatch({
      if (arm$engine == "random")            select_crosses_random(parents, cfg)
      else if (arm$engine == "simplemating") select_crosses_simplemating(parents, SP, cfg, cfg$seed + cyc, sm_external)
      else                                   select_crosses_package(arm, parents, SP, cfg, cfg$seed + cyc)
    }, error = function(e) { message("arm ", arm$arm, " rep ", rep_i, " cyc ", cyc, ": ", conditionMessage(e)); NULL })
    if (is.null(crossPlan) || nrow(crossPlan) < 1L) break
    F1 <- AlphaSimR::makeCross(parents, crossPlan = crossPlan, nProgeny = 1, simParam = SP)
    adv <- advance_ril(F1, SP, cfg)
    m <- breeder_metrics(adv$next_parents, SP, base_mean)
    rows[[cyc]] <- cbind(
      data.frame(rep = rep_i, arm = arm$arm, engine = arm$engine,
                 metric = arm$trait_value_metric %||% NA_character_,
                 allocation = arm$allocation_method %||% NA_character_,
                 optimizer = arm$optimizer %||% NA_character_,
                 strategy = arm$strategy %||% NA_character_,
                 cycle = cyc, n_crosses_used = nrow(crossPlan),
                 n_unique_parents = length(unique(as.vector(crossPlan))),
                 implementation = if (arm$engine == "simplemating")
                   (if (sm_external) "simplemating_external" else "simplemating_native_proxy_uc_vpm") else "package",
                 exact_external_status = if (arm$engine == "simplemating")
                   (if (sm_external) "exact" else "style_proxy") else NA_character_,
                 stringsAsFactors = FALSE),
      m)
    parents <- adv$next_parents
  }
  do.call(rbind, rows)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || is.na(a)) b else a

# ---- SimpleMating availability (build external, native fallback) ------------
resolve_simplemating <- function(cfg) {
  if (!cfg$include_simplemating) return(list(external = FALSE, reason = "disabled"))
  if (requireNamespace("SimpleMating", quietly = TRUE)) return(list(external = TRUE, reason = "installed"))
  if (cfg$simplemating_github && requireNamespace("remotes", quietly = TRUE)) {
    ok <- tryCatch({ remotes::install_github("Resende-Lab/SimpleMating", upgrade = "never", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok && requireNamespace("SimpleMating", quietly = TRUE)) return(list(external = TRUE, reason = "github_build"))
  }
  list(external = FALSE, reason = "unavailable_native_fallback")
}

resolve_workers <- function(req, n) {
  if (.Platform$OS.type == "windows") return(1L)
  if (req > 0L) return(min(req, n))
  max(1L, min(n, tryCatch(parallel::detectCores(), error = function(e) 2L) - 2L))
}

# ---- main -------------------------------------------------------------------
run_method_optimizer_metric_study <- function(cfg = mom_config()) {
  source(file.path(cfg$root, "R", "load.R"))
  ng_load(cfg$root, use_cpp = cfg$use_cpp, verbose = FALSE)
  if (!requireNamespace("AlphaSimR", quietly = TRUE))
    stop("AlphaSimR is required for this study.", call. = FALSE)
  sm <- resolve_simplemating(cfg)
  message(sprintf("SimpleMating: %s (%s)", if (sm$external) "external" else "native fallback", sm$reason))

  arms <- arm_registry()
  if (!sm$external) message("Note: 'simplemating' arm uses the UC-VPM native proxy (labelled style_proxy).")

  # founders once per rep (shared across arms -> paired comparison)
  message(sprintf("Generating %d founder set(s)...", cfg$reps))
  founders_by_rep <- lapply(seq_len(cfg$reps), function(r) {
    set.seed(cfg$seed + r); build_founders(cfg)
  })

  grid <- expand.grid(rep_i = seq_len(cfg$reps), arm_idx = seq_len(nrow(arms)),
                      KEEP.OUT.ATTRS = FALSE)
  n_jobs <- nrow(grid)
  workers <- resolve_workers(cfg$workers, n_jobs)
  message(sprintf("Racing %d arms x %d reps = %d jobs on %d worker(s); %d cycles each.",
                  nrow(arms), cfg$reps, n_jobs, workers, cfg$cycles))

  run_job <- function(i) {
    r <- grid$rep_i[i]; a <- grid$arm_idx[i]
    run_arm_rep(cfg, arms[a, , drop = FALSE], r, a, founders_by_rep[[r]], sm$external)
  }
  res_list <- if (workers > 1L)
    parallel::mclapply(seq_len(n_jobs), run_job, mc.cores = workers, mc.preschedule = FALSE)
  else lapply(seq_len(n_jobs), run_job)

  res <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
  if (is.null(res) || !nrow(res)) stop("All jobs failed; check AlphaSimR / package load.", call. = FALSE)

  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  raw_csv <- file.path(cfg$output_dir, paste0(cfg$output_prefix, "_raw.csv"))
  utils::write.csv(res, raw_csv, row.names = FALSE)
  saveRDS(res, file.path(cfg$output_dir, paste0(cfg$output_prefix, "_raw.rds")))

  # final-cycle summary per arm across reps
  fin <- res[res$cycle == max(res$cycle), ]
  agg <- aggregate(cbind(genetic_gain_index, var_index_gv, expected_heterozygosity,
                         mean_parent_relationship) ~ arm, fin,
                   function(x) mean(x, na.rm = TRUE))
  agg <- agg[order(-agg$genetic_gain_index), ]
  cat(sprintf("\n=== Final cycle (%d) mean over %d reps, ranked by genetic gain ===\n",
              max(res$cycle), cfg$reps))
  cat(sprintf("  %-24s %8s %10s %8s %10s\n", "arm", "gain", "add_var", "He", "meanRel"))
  for (i in seq_len(nrow(agg)))
    cat(sprintf("  %-24s %8.3f %10.3f %8.4f %10.4f\n", agg$arm[i],
                agg$genetic_gain_index[i], agg$var_index_gv[i],
                agg$expected_heterozygosity[i], agg$mean_parent_relationship[i]))
  message("\nWrote: ", raw_csv)
  invisible(list(raw = res, summary = agg))
}

is_this_script <- function() any(grepl("run_method_optimizer_metric_study\\.R$", commandArgs(FALSE)))
if (is_this_script()) invisible(run_method_optimizer_metric_study())
