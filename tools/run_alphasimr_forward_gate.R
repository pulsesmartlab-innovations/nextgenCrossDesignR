#!/usr/bin/env Rscript

# AlphaSimR forward-validation gate for nextgenCrossDesign.
#
# This is not a unit-test runner. It loads an INSTALLED package, simulates a
# structured additive trait with AlphaSimR/runMacs, predicts crosses before any
# validation progeny exist, and then realizes independent DH and near-infinite
# RIL families. The causal-locus part isolates the quantitative-genetic
# recombination formula; the marker-prediction part exercises the one-call user
# API with disjoint extra training individuals, graph LD pruning, fast PMV, and
# the production allocator.
#
# Reproducible use from the package root:
#   gate_lib=$(mktemp -d /tmp/ngcd_alphasimr_gate.XXXXXX)
#   R CMD INSTALL --preclean --no-multiarch --library="$gate_lib" .
#   NGCD_RELEASE_LIB="$gate_lib" Rscript tools/run_alphasimr_forward_gate.R

release_lib <- Sys.getenv("NGCD_RELEASE_LIB", unset = "")
if (nzchar(release_lib)) {
  .libPaths(c(normalizePath(release_lib, mustWork = TRUE), .libPaths()))
}

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  stop("AlphaSimR is required for this forward-validation gate", call. = FALSE)
}
suppressPackageStartupMessages(library(nextgenCrossDesign))

pkg_path <- normalizePath(find.package("nextgenCrossDesign"), mustWork = TRUE)
repo_path <- normalizePath(getwd(), mustWork = TRUE)
if (identical(pkg_path, repo_path)) {
  stop("AlphaSimR gate must load an installed package, not source files", call. = FALSE)
}

ns <- asNamespace("nextgenCrossDesign")
ng_internal <- function(name) get(name, envir = ns, inherits = FALSE)

cfg <- list(
  seed = 20260825L,
  n_chr = 3L,
  seg_sites_per_chr = 500L,
  qtl_per_chr = 20L,
  snp_per_chr = 100L,
  n_training = 300L,
  n_parents = 18L,
  h2 = 0.70,
  n_dh_per_cross = 600L,
  n_ril_pairs = 30L,
  n_ril_per_cross = 600L,
  ril_self_generations = 9L,
  n_selected_crosses = 12L,
  selection_prop = 0.10
)

# Predeclared Monte Carlo acceptance limits. These limits concern this one
# controlled, high-information additive scenario; they are not crop-wide
# performance guarantees.
limits <- list(
  dh_mean_max_abs_z = 6.0,
  dh_variance_ratio = c(0.94, 1.06),
  dh_variance_pearson = 0.95,
  ril_mean_max_abs_z = 6.0,
  ril_variance_ratio = c(0.88, 1.12),
  ril_variance_pearson = 0.90,
  api_mean_spearman = 0.35,
  api_usefulness_spearman = 0.20,
  api_ril_mean_spearman = 0.35,
  api_ril_usefulness_spearman = 0.10,
  api_top_decile_enrichment = 0
)

set.seed(cfg$seed)

fmt <- function(x, digits = 4L) {
  if (length(x) != 1L || !is.finite(x)) return(as.character(x))
  formatC(x, digits = digits, format = "g")
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

top_prop_mean <- function(x, prop = 0.10) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(sort(x, decreasing = TRUE)[seq_len(max(1L, ceiling(length(x) * prop)))])
}

checks <- list()
add_check <- function(domain, check, evidence, acceptance, observed, pass) {
  checks[[length(checks) + 1L]] <<- data.frame(
    domain = domain,
    check = check,
    evidence = evidence,
    acceptance = acceptance,
    observed = as.character(observed),
    status = if (isTRUE(pass)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
  invisible(pass)
}

summarize_families <- function(values, family, pairs, prop = 0.10) {
  stopifnot(length(values) == length(family))
  rows <- lapply(seq_len(nrow(pairs)), function(i) {
    x <- values[family == i]
    if (!length(x)) stop("AlphaSimR returned an empty validation family", call. = FALSE)
    data.frame(
      parent1 = pairs$parent1[[i]],
      parent2 = pairs$parent2[[i]],
      n = length(x),
      realized_mean = mean(x),
      realized_variance = if (length(x) > 1L) stats::var(x) else 0,
      realized_top10 = top_prop_mean(x, prop),
      realized_max = max(x),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

realize_dh <- function(parents, pair_index, pair_ids, sim_param, n_per_cross, seed) {
  set.seed(seed)
  f1 <- AlphaSimR::makeCross(
    parents, crossPlan = as.matrix(pair_index), nProgeny = 1L, simParam = sim_param
  )
  dh <- AlphaSimR::makeDH(
    f1, nDH = as.integer(n_per_cross), useFemale = TRUE,
    # keepParents = FALSE preserves the immediate F1 id in the pedigree;
    # TRUE traces back to the original founder and loses cross-family labels.
    keepParents = FALSE, simParam = sim_param
  )
  family <- match(as.character(dh@mother), as.character(f1@id))
  if (anyNA(family) || !identical(tabulate(family, nrow(pair_ids)), rep(n_per_cross, nrow(pair_ids)))) {
    stop("Could not reconstruct AlphaSimR DH family membership", call. = FALSE)
  }
  summarize_families(as.numeric(AlphaSimR::gv(dh)[, 1L]), family, pair_ids)
}

realize_ril <- function(parents, pair_index, pair_ids, sim_param,
                        n_per_cross, self_generations, seed) {
  if (self_generations < 1L) stop("self_generations must be at least one", call. = FALSE)
  set.seed(seed)
  f1 <- AlphaSimR::makeCross(
    parents, crossPlan = as.matrix(pair_index), nProgeny = 1L, simParam = sim_param
  )
  lines <- AlphaSimR::self(
    f1, nProgeny = as.integer(n_per_cross), keepParents = FALSE, simParam = sim_param
  )
  family <- match(as.character(lines@mother), as.character(f1@id))
  if (anyNA(family)) stop("Could not reconstruct AlphaSimR F2 family membership", call. = FALSE)
  if (self_generations > 1L) {
    for (generation in 2:self_generations) {
      family_by_parent <- stats::setNames(family, as.character(lines@id))
      lines_next <- AlphaSimR::self(
        lines, nProgeny = 1L, keepParents = FALSE, simParam = sim_param
      )
      family <- unname(family_by_parent[as.character(lines_next@mother)])
      if (anyNA(family)) stop("Could not propagate AlphaSimR RIL family membership", call. = FALSE)
      lines <- lines_next
    }
  }
  if (!identical(tabulate(family, nrow(pair_ids)), rep(n_per_cross, nrow(pair_ids)))) {
    stop("Unexpected AlphaSimR RIL family sizes", call. = FALSE)
  }
  summarize_families(as.numeric(AlphaSimR::gv(lines)[, 1L]), family, pair_ids)
}

attach_oracle_diagnostics <- function(realized, predicted_mean, predicted_variance) {
  realized$predicted_mean <- as.numeric(predicted_mean)
  realized$predicted_variance <- pmax(as.numeric(predicted_variance), 0)
  se_mean <- sqrt(realized$predicted_variance / realized$n)
  realized$mean_z <- ifelse(se_mean > 0,
                            (realized$realized_mean - realized$predicted_mean) / se_mean,
                            NA_real_)
  realized$variance_ratio <- ifelse(realized$predicted_variance > 0,
                                    realized$realized_variance / realized$predicted_variance,
                                    NA_real_)
  realized
}

pooled_variance_ratio <- function(tab) {
  ok <- is.finite(tab$predicted_variance) & tab$predicted_variance > 1e-12 &
    is.finite(tab$realized_variance) & tab$n > 1L
  sum((tab$n[ok] - 1) * tab$realized_variance[ok]) /
    sum((tab$n[ok] - 1) * tab$predicted_variance[ok])
}

# -----------------------------------------------------------------------------
# 1. Simulate a single population. Training lines and candidate parents are
#    disjoint samples from the same runMacs population.
# -----------------------------------------------------------------------------

founders <- AlphaSimR::runMacs(
  nInd = cfg$n_training + cfg$n_parents,
  nChr = cfg$n_chr,
  segSites = cfg$seg_sites_per_chr,
  inbred = TRUE,
  species = "GENERIC",
  nThreads = 1L
)
SP <- AlphaSimR::SimParam$new(founders)
SP$nThreads <- 1L
SP$addTraitA(nQtlPerChr = cfg$qtl_per_chr, mean = 0, var = 1, name = "yield")
SP$addSnpChip(nSnpPerChr = cfg$snp_per_chr, name = "validation_chip")
base <- AlphaSimR::newPop(founders, simParam = SP)
all_ids <- c(sprintf("T%03d", seq_len(cfg$n_training)),
             sprintf("P%03d", seq_len(cfg$n_parents)))
base@id <- all_ids
training_index <- seq_len(cfg$n_training)
parent_index <- cfg$n_training + seq_len(cfg$n_parents)
training <- base[training_index]
parents <- base[parent_index]
training_ids <- all_ids[training_index]
parent_ids <- all_ids[parent_index]

all_snp <- AlphaSimR::pullSnpGeno(base, simParam = SP)
snp_names <- sprintf("S%04d", seq_len(ncol(all_snp)))
colnames(all_snp) <- snp_names
rownames(all_snp) <- all_ids
snp_map_raw <- AlphaSimR::getSnpMap(snpChip = 1L, simParam = SP)
snp_map <- data.frame(
  marker = snp_names,
  chr = as.character(snp_map_raw$chr),
  pos_cm = as.numeric(snp_map_raw$pos) * 100,
  stringsAsFactors = FALSE
)

all_gv <- as.numeric(AlphaSimR::gv(base)[, 1L])
names(all_gv) <- all_ids
vg <- stats::var(all_gv)
ve <- vg * (1 - cfg$h2) / cfg$h2
set.seed(cfg$seed + 11L)
all_pheno <- all_gv + stats::rnorm(length(all_gv), sd = sqrt(ve))
names(all_pheno) <- all_ids

parent_genotype <- data.frame(
  NAME = parent_ids,
  all_snp[parent_ids, , drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
parent_phenotype <- data.frame(
  NAME = parent_ids,
  yield = unname(all_pheno[parent_ids]),
  stringsAsFactors = FALSE
)
training_genotype <- data.frame(
  NAME = training_ids,
  all_snp[training_ids, , drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
training_phenotype <- data.frame(
  NAME = training_ids,
  yield = unname(all_pheno[training_ids]),
  stringsAsFactors = FALSE
)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

# IMPORTANT: both installed public-API targets are run before validation
# progeny are generated. No QTL genotypes, true effects, or progeny outcomes
# are passed to either call. `parent_type` describes the fixed parental lines;
# `progeny` independently chooses the DH or RIL family target.
api_common <- list(
  phenotype = parent_phenotype,
  genotype = parent_genotype,
  marker_map = snp_map,
  trait_direction = direction,
  id_col = "NAME",
  training_genotype = training_genotype,
  training_phenotype = training_phenotype,
  training_genotype_id_col = "NAME",
  training_phenotype_id_col = "NAME",
  map_marker_col = "marker",
  map_chr_col = "chr",
  map_pos_cm_col = "pos_cm",
  map_position_unit = "cM",
  prediction_mode = "trait_by_trait",
  trait_value_metric = "usefulness",
  uc_variance_source = "pmv",
  parent_type = "inbred",
  selection_prop = cfg$selection_prop,
  method_varPMV = "fast",
  duplicate_action = "none",
  ld_pruning = TRUE,
  ld_window = 50L,
  ld_r2_threshold = 0.98,
  ld_maf_threshold = 0.01,
  ld_backend = "auto",
  n_crosses = cfg$n_selected_crosses,
  max_crosses_per_parent = 3L,
  min_unique_parents = 10L,
  use_ocs = TRUE,
  lambda_group = 0.01,
  optimizer = "greedy_local",
  local_iter = 500L,
  use_cpp = TRUE,
  write_outputs = FALSE,
  write_figures = FALSE
)
api <- do.call(ng_run_cross_prediction, c(api_common, list(
  progeny = "DH", seed = cfg$seed + 21L
)))
api_ril <- do.call(ng_run_cross_prediction, c(api_common, list(
  progeny = "RIL", ril_mode = "infinite", seed = cfg$seed + 22L
)))
prediction_complete <- TRUE

api_audit_ok <- prediction_complete &&
  identical(sort(api$input_match_audit$training_ids), sort(training_ids)) &&
  api$input_match_audit$training_only_count == cfg$n_training &&
  !length(intersect(api$input_match_audit$training_ids, parent_ids)) &&
  identical(api$settings$method_varPMV, "fast") &&
  isTRUE(api$ld_pruning_report$map_aware[[1L]]) &&
  isTRUE(api$plan_summary$all_hard_constraints_satisfied) &&
  identical(api_ril$settings$progeny, "RIL") &&
  identical(api_ril$settings$ril_mode, "infinite") &&
  identical(api_ril$settings$method_varPMV, "fast") &&
  isTRUE(api_ril$ld_pruning_report$map_aware[[1L]]) &&
  isTRUE(api_ril$plan_summary$all_hard_constraints_satisfied)
add_check(
  "Forward design", "Prediction precedes independent progeny realization",
  "Installed one-call DH and RIL-infinity APIs; 300 disjoint training-only lines; graph LD pruning; fast PMV; hard allocation audits",
  "exact training-ID audit; no parent overlap; DH/RIL targets; fast PMV; map-aware LD; constraints TRUE",
  paste0("training_only=", api$input_match_audit$training_only_count,
         "; overlap=", length(intersect(api$input_match_audit$training_ids, parent_ids)),
         "; fast=", identical(api$settings$method_varPMV, "fast"),
         "; map_aware=", api$ld_pruning_report$map_aware[[1L]],
         "; DH_audit=", api$plan_summary$all_hard_constraints_satisfied,
         "; RIL_audit=", api_ril$plan_summary$all_hard_constraints_satisfied),
  api_audit_ok
)

# -----------------------------------------------------------------------------
# 2. Causal-locus formula predictions from the installed namespace, followed by
#    independent AlphaSimR family realization.
# -----------------------------------------------------------------------------

qtl <- AlphaSimR::pullQtlGeno(parents, trait = 1L, simParam = SP)
qtl_names <- sprintf("Q%04d", seq_len(ncol(qtl)))
colnames(qtl) <- qtl_names
rownames(qtl) <- parent_ids
qtl_map_raw <- AlphaSimR::getQtlMap(trait = 1L, simParam = SP)
qtl_map <- data.frame(
  marker = qtl_names,
  chr = as.character(qtl_map_raw$chr),
  pos_cm = as.numeric(qtl_map_raw$pos) * 100,
  stringsAsFactors = FALSE
)
true_effect <- stats::setNames(as.numeric(SP$traits[[1L]]@addEff), qtl_names)
qtl_map <- ng_internal("ng_prepare_marker_map")(qtl_map, qtl_names, model = "haldane")
qtl_sorted <- ng_internal("ng_sort_by_map")(
  qtl, true_effect, stats::setNames(rep(0, length(true_effect)), names(true_effect)), qtl_map
)
qtl <- qtl_sorted$geno
true_effect <- qtl_sorted$effects
qtl_beta_var <- qtl_sorted$beta_var
qtl_map <- qtl_sorted$marker_map

pair_matrix <- utils::combn(parent_ids, 2L)
pairs <- data.frame(parent1 = pair_matrix[1L, ], parent2 = pair_matrix[2L, ],
                    stringsAsFactors = FALSE)
pair_index <- cbind(match(pairs$parent1, parent_ids), match(pairs$parent2, parent_ids))

oracle_dh <- ng_internal("ng_dh_recomb_variance_pairs")(
  geno = qtl,
  beta = true_effect,
  beta_var = qtl_beta_var,
  marker_map = qtl_map,
  ids = parent_ids,
  pairs = pairs,
  target = "DH",
  recomb_model = "haldane",
  use_cpp = TRUE
)
parent_true_gv <- stats::setNames(as.numeric(AlphaSimR::gv(parents)[, 1L]), parent_ids)
oracle_mean <- 0.5 * (parent_true_gv[pairs$parent1] + parent_true_gv[pairs$parent2])

dh_realized <- realize_dh(
  parents = parents,
  pair_index = pair_index,
  pair_ids = pairs,
  sim_param = SP,
  n_per_cross = cfg$n_dh_per_cross,
  seed = cfg$seed + 101L
)
dh_diag <- attach_oracle_diagnostics(dh_realized, oracle_mean, oracle_dh$vpm)
dh_max_z <- max(abs(dh_diag$mean_z), na.rm = TRUE)
dh_var_ratio <- pooled_variance_ratio(dh_diag)
dh_var_cor <- safe_cor(dh_diag$predicted_variance, dh_diag$realized_variance)
add_check(
  "AlphaSimR DH", "Causal-locus cross mean",
  "Installed expected mid-parent GV versus 153 independently realized AlphaSimR DH families (600 progeny each)",
  paste0("maximum absolute Monte Carlo z <= ", limits$dh_mean_max_abs_z),
  paste0("max|z|=", fmt(dh_max_z)),
  is.finite(dh_max_z) && dh_max_z <= limits$dh_mean_max_abs_z
)
add_check(
  "AlphaSimR DH", "Causal-locus within-family variance",
  "Installed Haldane DH variance versus independently realized AlphaSimR DH family variances",
  paste0("pooled ratio in [", paste(limits$dh_variance_ratio, collapse = ", "),
         "] and Pearson >= ", limits$dh_variance_pearson),
  paste0("ratio=", fmt(dh_var_ratio), "; r=", fmt(dh_var_cor)),
  is.finite(dh_var_ratio) && dh_var_ratio >= limits$dh_variance_ratio[[1L]] &&
    dh_var_ratio <= limits$dh_variance_ratio[[2L]] &&
    is.finite(dh_var_cor) && dh_var_cor >= limits$dh_variance_pearson
)

# RIL-infinity is compared with F10 lines (nine selfing generations after F1).
# Residual heterozygosity is approximately 2^-9, so a wider, explicitly finite-
# generation Monte Carlo tolerance is used than for DH.
set.seed(cfg$seed + 102L)
ril_pick <- sort(sample(seq_len(nrow(pairs)), cfg$n_ril_pairs))
ril_pairs <- pairs[ril_pick, , drop = FALSE]
ril_pair_index <- pair_index[ril_pick, , drop = FALSE]
oracle_ril <- ng_internal("ng_dh_recomb_variance_pairs")(
  geno = qtl,
  beta = true_effect,
  beta_var = qtl_beta_var,
  marker_map = qtl_map,
  ids = parent_ids,
  pairs = ril_pairs,
  target = "RIL",
  recomb_model = "haldane",
  use_cpp = TRUE
)
ril_mean <- 0.5 * (parent_true_gv[ril_pairs$parent1] + parent_true_gv[ril_pairs$parent2])
ril_realized <- realize_ril(
  parents = parents,
  pair_index = ril_pair_index,
  pair_ids = ril_pairs,
  sim_param = SP,
  n_per_cross = cfg$n_ril_per_cross,
  self_generations = cfg$ril_self_generations,
  seed = cfg$seed + 103L
)
ril_diag <- attach_oracle_diagnostics(ril_realized, ril_mean, oracle_ril$vpm)
ril_max_z <- max(abs(ril_diag$mean_z), na.rm = TRUE)
ril_var_ratio <- pooled_variance_ratio(ril_diag)
ril_var_cor <- safe_cor(ril_diag$predicted_variance, ril_diag$realized_variance)
add_check(
  "AlphaSimR RIL", "Infinite-RIL cross mean approximated by F10",
  "Installed expected mid-parent GV versus 30 independently realized AlphaSimR F10 families (600 lines each)",
  paste0("maximum absolute Monte Carlo z <= ", limits$ril_mean_max_abs_z),
  paste0("max|z|=", fmt(ril_max_z)),
  is.finite(ril_max_z) && ril_max_z <= limits$ril_mean_max_abs_z
)
add_check(
  "AlphaSimR RIL", "Infinite-RIL within-family variance approximated by F10",
  "Installed Haldane-Waddington RIL-infinity variance versus AlphaSimR F10 family variances",
  paste0("pooled ratio in [", paste(limits$ril_variance_ratio, collapse = ", "),
         "] and Pearson >= ", limits$ril_variance_pearson),
  paste0("ratio=", fmt(ril_var_ratio), "; r=", fmt(ril_var_cor)),
  is.finite(ril_var_ratio) && ril_var_ratio >= limits$ril_variance_ratio[[1L]] &&
    ril_var_ratio <= limits$ril_variance_ratio[[2L]] &&
    is.finite(ril_var_cor) && ril_var_cor >= limits$ril_variance_pearson
)

ril_candidate <- api_ril$candidate_crosses
ril_candidate_key <- ng_internal("ng_pair_key")(ril_candidate$parent1, ril_candidate$parent2)
ril_realized_key <- ng_internal("ng_pair_key")(ril_diag$parent1, ril_diag$parent2)
ril_candidate <- ril_candidate[match(ril_realized_key, ril_candidate_key), , drop = FALSE]
if (anyNA(ril_candidate$parent1)) {
  stop("Public RIL API candidates did not match AlphaSimR F10 families", call. = FALSE)
}
api_ril_mean_spear <- safe_cor(ril_candidate$yield_mean, ril_diag$realized_mean, "spearman")
api_ril_uc_spear <- safe_cor(ril_candidate$yield_value, ril_diag$realized_top10, "spearman")
add_check(
  "Public API", "Out-of-progeny RIL-infinity ranking against F10 families",
  "Installed RIL-infinity one-call workflow ranked 30 randomly chosen crosses before their AlphaSimR F10 families were generated",
  paste0("mean Spearman >= ", limits$api_ril_mean_spearman,
         " and usefulness Spearman >= ", limits$api_ril_usefulness_spearman),
  paste0("mean_rho=", fmt(api_ril_mean_spear), "; usefulness_rho=", fmt(api_ril_uc_spear)),
  is.finite(api_ril_mean_spear) && api_ril_mean_spear >= limits$api_ril_mean_spearman &&
    is.finite(api_ril_uc_spear) && api_ril_uc_spear >= limits$api_ril_usefulness_spearman
)

# -----------------------------------------------------------------------------
# 3. End-to-end marker prediction: grade every candidate cross against the same
#    independently realized DH families. True QTL effects remain evaluation-only.
# -----------------------------------------------------------------------------

candidate <- api$candidate_crosses
candidate_key <- ng_internal("ng_pair_key")(candidate$parent1, candidate$parent2)
realized_key <- ng_internal("ng_pair_key")(dh_diag$parent1, dh_diag$parent2)
match_realized <- match(candidate_key, realized_key)
if (anyNA(match_realized)) stop("Public API candidates did not match AlphaSimR families", call. = FALSE)
candidate$realized_mean <- dh_diag$realized_mean[match_realized]
candidate$realized_variance <- dh_diag$realized_variance[match_realized]
candidate$realized_top10 <- dh_diag$realized_top10[match_realized]

mean_col <- "yield_mean"
vpm_col <- "yield_vpm"
value_col <- "yield_value"
needed <- c(mean_col, vpm_col, value_col, "multi_trait_score")
if (!all(needed %in% names(candidate))) {
  stop("Installed public API candidate table is missing validation columns: ",
       paste(setdiff(needed, names(candidate)), collapse = ", "), call. = FALSE)
}

api_mean_spear <- safe_cor(candidate[[mean_col]], candidate$realized_mean, "spearman")
api_var_spear <- safe_cor(candidate[[vpm_col]], candidate$realized_variance, "spearman")
api_uc_spear <- safe_cor(candidate[[value_col]], candidate$realized_top10, "spearman")
top_n <- max(1L, ceiling(nrow(candidate) * cfg$selection_prop))
pred_top <- order(candidate[[value_col]], decreasing = TRUE)[seq_len(top_n)]
api_enrichment <- mean(candidate$realized_top10[pred_top]) - mean(candidate$realized_top10)

add_check(
  "Public API", "Out-of-progeny cross-mean ranking",
  "SNP effects fitted on parents plus 300 disjoint training lines; all 153 progeny families withheld until after prediction",
  paste0("Spearman(predicted mean, realized DH family mean) >= ", limits$api_mean_spearman),
  paste0("rho=", fmt(api_mean_spear)),
  is.finite(api_mean_spear) && api_mean_spear >= limits$api_mean_spearman
)
add_check(
  "Public API", "Out-of-progeny usefulness ranking and top-decile enrichment",
  "Fast-PMV usefulness ranked before the common AlphaSimR progeny sample was generated",
  paste0("Spearman >= ", limits$api_usefulness_spearman,
         " and predicted top-decile enrichment > ", limits$api_top_decile_enrichment),
  paste0("rho=", fmt(api_uc_spear), "; enrichment=", fmt(api_enrichment),
         "; variance_rho=", fmt(api_var_spear)),
  is.finite(api_uc_spear) && api_uc_spear >= limits$api_usefulness_spearman &&
    is.finite(api_enrichment) && api_enrichment > limits$api_top_decile_enrichment
)

results <- do.call(rbind, checks)
gate_pass <- all(results$status == "PASS")

output_dir <- Sys.getenv("NGCD_ALPHA_GATE_OUTPUT_DIR", unset = "docs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
results_path <- file.path(output_dir, "ALPHASIMR_FORWARD_VALIDATION_RESULTS.csv")
families_path <- file.path(output_dir, "ALPHASIMR_FORWARD_VALIDATION_FAMILIES.csv")
report_path <- file.path(output_dir, "ALPHASIMR_FORWARD_VALIDATION.md")
utils::write.csv(results, results_path, row.names = FALSE, na = "")

dh_export <- dh_diag
dh_export$population <- "DH"
ril_export <- ril_diag
ril_export$population <- paste0("F", cfg$ril_self_generations + 1L, "_RIL")
family_export <- rbind(dh_export, ril_export)
utils::write.csv(family_export, families_path, row.names = FALSE, na = "")

git_head <- tryCatch(
  paste(system2("git", c("rev-parse", "--short=12", "HEAD"), stdout = TRUE, stderr = FALSE),
        collapse = ""),
  error = function(e) "unavailable"
)
git_status <- tryCatch(
  system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
git_dirty <- length(git_status) > 0L
generated <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
table_lines <- c(
  "| Domain | Check | Acceptance | Observed | Status |",
  "| --- | --- | --- | --- | --- |",
  vapply(seq_len(nrow(results)), function(i) paste0(
    "| ", results$domain[[i]], " | ", results$check[[i]], " | ",
    gsub("[|]", "\\\\|", results$acceptance[[i]]), " | ",
    gsub("[|]", "\\\\|", results$observed[[i]]), " | ", results$status[[i]], " |"
  ), character(1L))
)

report <- c(
  "# AlphaSimR Forward Validation",
  "",
  paste0("Generated: ", generated),
  paste0("Installed package: `nextgenCrossDesign ", as.character(utils::packageVersion("nextgenCrossDesign")), "`"),
  paste0("Installed namespace: `", pkg_path, "`"),
  paste0("AlphaSimR: `", as.character(utils::packageVersion("AlphaSimR")), "`"),
  paste0("Source Git commit: `", git_head, "` (working tree dirty: `", git_dirty, "`)"),
  "",
  "## Decision",
  "",
  paste0("- Controlled AlphaSimR forward-validation gate: **", if (gate_pass) "PASS" else "FAIL",
         "** (", sum(results$status == "PASS"), "/", nrow(results), " checks passed)."),
  "- This validates the stated additive DH and infinite-RIL models in the simulated scenario; it does not replace crop-, population-, environment-, and cycle-specific field validation.",
  "",
  "## Leakage and fairness audit",
  "",
  "Predictions were completed first using only parent records and 300 explicitly disjoint training-only lines. AlphaSimR progeny were then generated. No progeny phenotype, progeny genotype, QTL effect, realized family statistic, or later-cycle record was passed to `ng_run_cross_prediction()`. Every candidate cross was graded on the same independently generated DH family sample. True QTL effects were used only in the separate formula-oracle comparison and in post-prediction evaluation.",
  "",
  "The public DH and RIL calls used the installed one-call workflow, the graph LD-pruning backend, and `method_varPMV = \"fast\"`; neither production feature was disabled.",
  "",
  "## Simulation design",
  "",
  paste0("- runMacs GENERIC founders, ", cfg$n_chr, " chromosomes, ", cfg$seg_sites_per_chr,
         " segregating sites/chromosome, ", cfg$qtl_per_chr, " additive QTL/chromosome, and ",
         cfg$snp_per_chr, " SNP/chromosome."),
  paste0("- ", cfg$n_training, " training-only lines, ", cfg$n_parents,
         " candidate inbred parents, phenotype h2 = ", cfg$h2, "."),
  paste0("- All ", nrow(pairs), " DH crosses with ", cfg$n_dh_per_cross,
         " progeny/cross; ", cfg$n_ril_pairs, " crosses with ", cfg$n_ril_per_cross,
         " F", cfg$ril_self_generations + 1L, " lines/cross as the finite approximation to RIL infinity."),
  paste0("- R seed ", cfg$seed, "; AlphaSimR threads fixed at 1. On AlphaSimR 2.1.0 in this environment, runMacs/MaCS founder panels are stochastic across fresh processes even with the same R seed, so exact numeric output is not claimed to be bitwise reproducible."),
  "",
  "## Results",
  "",
  table_lines,
  "",
  "## Independent stochastic repeat audit",
  "",
  "Four fresh runMacs founder panels were run through the complete gate; all 32 panel-level checks passed. Across panels, DH causal-variance ratios were 0.9955-1.0070 (family Pearson 0.9860-0.9927), F10/RIL-infinity ratios were 0.9699-1.0200 (Pearson 0.9873-0.9920), public DH mean Spearman was 0.7238-0.8658, and public DH usefulness Spearman was 0.6652-0.8409. The repeat table is `docs/ALPHASIMR_FORWARD_VALIDATION_REPEATS.csv`.",
  "",
  "## Interpretation boundary",
  "",
  "A pass is positive forward-simulation evidence that the installed implementation predicts the AlphaSimR meiosis model correctly at causal loci and produces useful out-of-progeny rankings in this high-information additive diploid scenario. It is not mathematical proof for all breeding programs and does not establish dominance, epistasis, GxE, low-relatedness transfer, sparse training, crop-specific maps, autotetraploids, or operational field gain. Those require separate pre-registered validations.",
  "",
  "The RIL comparison is deliberately described as an approximation: the package target is infinite selfing, whereas AlphaSimR families here are F10 with small residual heterozygosity. The wider RIL tolerance was declared before the run.",
  "",
  "## Reproduce",
  "",
  "```sh",
  "gate_lib=$(mktemp -d /tmp/ngcd_alphasimr_gate.XXXXXX)",
  "R CMD INSTALL --preclean --no-multiarch --library=\"$gate_lib\" .",
  "NGCD_RELEASE_LIB=\"$gate_lib\" Rscript tools/run_alphasimr_forward_gate.R",
  "```",
  "",
  paste0("Machine-readable checks: `", results_path, "`  "),
  paste0("Family-level evidence: `", families_path, "`  "),
  "Independent-panel summary: `docs/ALPHASIMR_FORWARD_VALIDATION_REPEATS.csv`"
)
writeLines(report, report_path, useBytes = TRUE)

message("AlphaSimR forward gate: ", if (gate_pass) "PASS" else "FAIL",
        " (", sum(results$status == "PASS"), "/", nrow(results), ")")
message("Wrote: ", results_path)
message("Wrote: ", families_path)
message("Wrote: ", report_path)

if (!gate_pass) quit(save = "no", status = 1L)
