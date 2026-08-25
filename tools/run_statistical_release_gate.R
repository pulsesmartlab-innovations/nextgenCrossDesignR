#!/usr/bin/env Rscript

# Statistical release gate for nextgenCrossDesign.
#
# This is deliberately not a unit-test runner. It loads an INSTALLED package
# and exercises the production namespace (including compiled kernels) against
# independently computed quantitative-genetic identities and hard decision
# invariants. It writes a machine-readable table and a breeder-facing report.
#
# Reproducible use from the package root:
#   R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib .
#   NGCD_RELEASE_LIB=/tmp/ngcd_release_gate_lib \
#     Rscript tools/run_statistical_release_gate.R

release_lib <- Sys.getenv("NGCD_RELEASE_LIB", unset = "")
if (nzchar(release_lib)) .libPaths(c(normalizePath(release_lib, mustWork = TRUE), .libPaths()))

suppressPackageStartupMessages(library(nextgenCrossDesign))
pkg_path <- normalizePath(find.package("nextgenCrossDesign"), mustWork = TRUE)
repo_path <- normalizePath(getwd(), mustWork = TRUE)
if (identical(pkg_path, repo_path)) {
  stop("Release gate must load an installed package, not the source directory", call. = FALSE)
}

ns <- asNamespace("nextgenCrossDesign")
ng_internal <- function(name) get(name, envir = ns, inherits = FALSE)
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

tol_identity <- 1e-10
tol_probability <- 1e-12
tol_psd <- 1e-10
set.seed(20260825)

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

max_abs <- function(x) {
  x <- abs(as.numeric(x))
  if (!length(x) || all(!is.finite(x))) return(Inf)
  max(x[is.finite(x)])
}

fmt <- function(x, digits = 4L) {
  if (length(x) != 1L || !is.finite(x)) return(as.character(x))
  formatC(x, digits = digits, format = "g")
}

expect_error <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

# -----------------------------------------------------------------------------
# 1. Diploid DH/RIL within-family variance and PMV
# -----------------------------------------------------------------------------

markers <- paste0("m", 1:3)
geno_small <- rbind(P1 = c(0, 2, 0), P2 = c(2, 0, 2))
colnames(geno_small) <- markers
beta_small <- stats::setNames(c(0.7, -1.1, 0.35), markers)
beta_var_small <- stats::setNames(c(0.04, 0.02, 0.01), markers)
map_small <- data.frame(marker = markers, chr = "1", pos_cm = c(0, 7, 24),
                        stringsAsFactors = FALSE)
map_small <- ng_internal("ng_prepare_marker_map")(map_small, markers, model = "haldane")
pairs_small <- data.frame(parent1 = "P1", parent2 = "P2", stringsAsFactors = FALSE)

dh_prod <- ng_internal("ng_dh_recomb_variance_pairs")(
  geno = geno_small, beta = beta_small, beta_var = beta_var_small,
  marker_map = map_small, ids = rownames(geno_small), pairs = pairs_small,
  use_cpp = TRUE, recomb_model = "haldane", target = "DH"
)

# Exact enumeration of all two-parent F1 gametic paths. A path starts on either
# parental homologue with probability 1/2 and switches with the Haldane
# recombination fraction in each interval.
states <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(markers))))
interval_r <- 0.5 * (1 - exp(-2 * diff(map_small$pos_cm) / 100))
path_probability <- apply(states, 1L, function(s) {
  0.5 * prod(ifelse(s[-1L] == s[-length(s)], 1 - interval_r, interval_r))
})
d_small <- 0.5 * (geno_small["P1", ] - geno_small["P2", ])
a_small <- d_small * beta_small
dh_deviation <- as.numeric(states %*% a_small)
dh_vpm_oracle <- sum(path_probability * dh_deviation^2)
dh_pmv_oracle <- dh_vpm_oracle + sum(d_small^2 * beta_var_small)
dh_error <- max_abs(c(dh_prod$vpm - dh_vpm_oracle, dh_prod$pmv - dh_pmv_oracle))
add_check(
  "Diploid variance", "DH exact gamete enumeration",
  "All 2^3 F1 gametic paths with exact Haldane transition probabilities",
  paste0("max absolute error <= ", tol_identity), fmt(dh_error), dh_error <= tol_identity
)

# RIL-infinity oracle from the Haldane-Waddington transformation R = 2r/(1+2r).
distance <- abs(outer(map_small$pos_cm, map_small$pos_cm, "-"))
r_meiosis <- 0.5 * (1 - exp(-2 * distance / 100))
R_ril <- 2 * r_meiosis / (1 + 2 * r_meiosis)
decay_ril <- 1 - 2 * R_ril
ril_vpm_oracle <- as.numeric(crossprod(a_small, decay_ril %*% a_small))
ril_pmv_oracle <- ril_vpm_oracle + sum(d_small^2 * beta_var_small)
ril_prod <- ng_internal("ng_dh_recomb_variance_pairs")(
  geno = geno_small, beta = beta_small, beta_var = beta_var_small,
  marker_map = map_small, ids = rownames(geno_small), pairs = pairs_small,
  use_cpp = TRUE, recomb_model = "haldane", target = "RIL"
)
ril_error <- max_abs(c(ril_prod$vpm - ril_vpm_oracle, ril_prod$pmv - ril_pmv_oracle))
add_check(
  "Diploid variance", "RIL-infinity Haldane-Waddington identity",
  "Independent pairwise RIL recombination transformation",
  paste0("max absolute error <= ", tol_identity), fmt(ril_error), ril_error <= tol_identity
)

# Random installed-kernel comparison: the fast Haldane-DH recursion must equal
# the dense quadratic form. The full-posterior path with a diagonal covariance
# must in turn reduce to fast PMV.
n_parent <- 12L
n_marker <- 42L
ids_fast <- sprintf("F%02d", seq_len(n_parent))
markers_fast <- sprintf("q%03d", seq_len(n_marker))
geno_fast <- matrix(2 * stats::rbinom(n_parent * n_marker, 1, 0.5), n_parent,
                    dimnames = list(ids_fast, markers_fast))
map_fast <- data.frame(
  marker = markers_fast,
  chr = rep(1:3, each = n_marker / 3),
  pos_cm = rep(seq(0, 90, length.out = n_marker / 3), 3),
  stringsAsFactors = FALSE
)
map_fast <- ng_internal("ng_prepare_marker_map")(map_fast, markers_fast, model = "haldane")
beta_fast <- stats::setNames(stats::rnorm(n_marker, 0, 0.2), markers_fast)
beta_var_fast <- stats::setNames(stats::runif(n_marker, 0, 0.03), markers_fast)
all_pairs_fast <- utils::combn(ids_fast, 2L)
pair_pick <- sort(sample(seq_len(ncol(all_pairs_fast)), 20L))
pairs_fast <- data.frame(parent1 = all_pairs_fast[1L, pair_pick],
                         parent2 = all_pairs_fast[2L, pair_pick], stringsAsFactors = FALSE)
fast_prod <- ng_internal("ng_dh_recomb_variance_pairs")(
  geno_fast, beta_fast, beta_var_fast, map_fast, ids_fast, pairs_fast,
  use_cpp = TRUE, recomb_model = "haldane", target = "DH"
)
dense_prod <- ng_internal("ng_dh_recomb_variance_pairs_dense")(
  geno_fast, beta_fast, beta_var_fast, map_fast, ids_fast, pairs_fast,
  recomb_model = "haldane", target = "DH"
)
fast_dense_error <- max_abs(as.matrix(fast_prod[c("vpm", "pmv")]) -
                              as.matrix(dense_prod[c("vpm", "pmv")]))
add_check(
  "Diploid variance", "Fast PMV retained: compiled recursion versus dense form",
  "Installed compiled Haldane-DH kernel versus chromosome-block quadratic form",
  paste0("max absolute error <= ", tol_identity), fmt(fast_dense_error),
  fast_dense_error <= tol_identity
)

Sigma_diag <- diag(beta_var_fast, nrow = n_marker)
dimnames(Sigma_diag) <- list(markers_fast, markers_fast)
full_diag <- ng_internal("ng_dh_recomb_variance_pairs_full_posterior")(
  geno_fast, beta_fast, Sigma_diag, map_fast, ids_fast, pairs_fast,
  target = "DH", recomb_model = "haldane", use_cpp = TRUE
)
diag_posterior_error <- max_abs(full_diag$pmv_full_posterior - fast_prod$pmv)
add_check(
  "Diploid variance", "Diagonal posterior reduction",
  "Full posterior covariance specialized to diagonal marker-effect covariance",
  paste0("max absolute error <= ", tol_identity), fmt(diag_posterior_error),
  diag_posterior_error <= tol_identity
)

# -----------------------------------------------------------------------------
# 2. Graph LD pruning: installed R/C++ production paths versus graph oracle
# -----------------------------------------------------------------------------

ld_oracle <- function(geno, marker_map, window, r2_threshold, maf_threshold, ploidy = 2) {
  marker_map <- marker_map[match(colnames(geno), marker_map$marker), , drop = FALSE]
  keep_names <- character()
  for (cc in unique(as.character(marker_map$chr))) {
    idx <- which(as.character(marker_map$chr) == cc)
    idx <- idx[order(as.numeric(marker_map$pos_cm[idx]), as.character(marker_map$marker[idx]))]
    g <- geno[, idx, drop = FALSE]
    maf <- pmin(colMeans(g, na.rm = TRUE) / ploidy,
                1 - colMeans(g, na.rm = TRUE) / ploidy)
    active <- is.finite(maf) & maf >= maf_threshold
    adjacency <- vector("list", ncol(g))
    if (ncol(g) > 1L) {
      for (i in seq_len(ncol(g) - 1L)) {
        if (!active[[i]]) next
        j_end <- min(ncol(g), i + as.integer(window))
        for (j in seq.int(i + 1L, j_end)) {
          if (!active[[j]]) next
          ok <- is.finite(g[, i]) & is.finite(g[, j])
          r2 <- if (sum(ok) >= 2L && stats::sd(g[ok, i]) > 0 && stats::sd(g[ok, j]) > 0) {
            stats::cor(g[ok, i], g[ok, j])^2
          } else NA_real_
          if (is.finite(r2) && r2 > r2_threshold) {
            adjacency[[i]] <- c(adjacency[[i]], j)
            adjacency[[j]] <- c(adjacency[[j]], i)
          }
        }
      }
    }
    unseen <- which(active)
    while (length(unseen)) {
      component <- integer()
      queue <- unseen[[1L]]
      while (length(queue)) {
        node <- queue[[1L]]
        queue <- queue[-1L]
        if (node %in% component) next
        component <- c(component, node)
        queue <- c(queue, setdiff(adjacency[[node]], component))
      }
      best_maf <- max(maf[component])
      representative <- component[which(maf[component] == best_maf)[[1L]]]
      keep_names <- c(keep_names, colnames(g)[representative])
      unseen <- setdiff(unseen, component)
    }
  }
  sort(keep_names)
}

ld_cases <- 60L
ld_pass <- logical(ld_cases)
ld_r_cpp_pass <- logical(ld_cases)
ld_order_pass <- logical(ld_cases)
for (case in seq_len(ld_cases)) {
  n <- sample(28:48, 1L)
  m <- sample(9:19, 1L)
  marker <- sprintf("ld%02d_%03d", case, seq_len(m))
  p <- stats::runif(m, 0.04, 0.5)
  g <- vapply(p, function(pk) stats::rbinom(n, 2L, pk), numeric(n))
  g <- matrix(g, nrow = n, dimnames = list(sprintf("i%03d", seq_len(n)), marker))
  missing <- matrix(stats::runif(n * m) < 0.04, nrow = n)
  g[missing] <- NA_real_
  chr <- sample(seq_len(sample(1:3, 1L)), m, replace = TRUE)
  mm <- data.frame(marker = marker, chr = chr, pos_cm = stats::runif(m, 0, 120),
                   stringsAsFactors = FALSE)
  window <- sample(2:6, 1L)
  threshold <- stats::runif(1L, 0.15, 0.85)
  maf_threshold <- stats::runif(1L, 0, 0.12)
  oracle <- ld_oracle(g, mm, window, threshold, maf_threshold)
  prod_auto <- ng_ld_prune_markers(
    g, window = window, r2_threshold = threshold, maf_threshold = maf_threshold,
    backend = "auto", marker_map = mm
  )$keep_markers
  prod_r <- ng_ld_prune_markers(
    g, window = window, r2_threshold = threshold, maf_threshold = maf_threshold,
    backend = "r", marker_map = mm
  )$keep_markers
  perm <- sample(seq_len(m))
  prod_permuted <- ng_ld_prune_markers(
    g[, perm, drop = FALSE], window = window, r2_threshold = threshold,
    maf_threshold = maf_threshold, backend = "auto", marker_map = mm
  )$keep_markers
  ld_pass[[case]] <- identical(sort(prod_auto), oracle)
  ld_r_cpp_pass[[case]] <- identical(sort(prod_auto), sort(prod_r))
  ld_order_pass[[case]] <- identical(sort(prod_auto), sort(prod_permuted))
}
add_check(
  "LD pruning", "Graph pruning versus independent connected-component oracle",
  "Randomized map-aware cases with missing calls, MAF filtering, and LD edges",
  paste0(ld_cases, "/", ld_cases, " exact marker-set matches"),
  paste0(sum(ld_pass), "/", ld_cases), all(ld_pass)
)
add_check(
  "LD pruning", "Installed C++/auto versus R graph parity",
  "Same randomized cases through both production backends",
  paste0(ld_cases, "/", ld_cases, " exact marker-set matches"),
  paste0(sum(ld_r_cpp_pass), "/", ld_cases), all(ld_r_cpp_pass)
)
add_check(
  "LD pruning", "Map-order invariance",
  "Marker columns permuted while chromosome/position map is held fixed",
  paste0(ld_cases, "/", ld_cases, " invariant marker sets"),
  paste0(sum(ld_order_pass), "/", ld_cases), all(ld_order_pass)
)
compiled_ld_available <- isTRUE(ng_internal("ng_ld_backend_available")())
add_check(
  "LD pruning", "Compiled graph backend loaded",
  "Installed namespace backend registry",
  "TRUE", compiled_ld_available, compiled_ld_available
)

# -----------------------------------------------------------------------------
# 3. Ridge-effect fitting and mean-source semantics
# -----------------------------------------------------------------------------

n_ridge <- 36L
m_ridge <- 15L
ridge_ids <- sprintf("R%02d", seq_len(n_ridge))
ridge_markers <- sprintf("r%02d", seq_len(m_ridge))
X <- matrix(stats::rbinom(n_ridge * m_ridge, 2L, 0.45), n_ridge,
            dimnames = list(ridge_ids, ridge_markers))
y <- as.numeric(12 + X %*% stats::rnorm(m_ridge, 0, 0.3) + stats::rnorm(n_ridge, 0, 0.5))
names(y) <- ridge_ids
shift <- stats::setNames(stats::runif(m_ridge, -5, 5), ridge_markers)
X_shift <- sweep(X, 2L, shift, "+")
fit_ridge <- ng_fit_ridge_effects(X, y, lambda = 5, kfold = 6, seed = 77)
fit_shift <- ng_fit_ridge_effects(X_shift, y, lambda = 5, kfold = 6, seed = 77)
translation_error <- max_abs(c(
  ng_predict_gebv(X, fit_ridge) - ng_predict_gebv(X_shift, fit_shift),
  fit_ridge$cv_predictive_r2 - fit_shift$cv_predictive_r2,
  fit_ridge$cv_predictive_correlation - fit_shift$cv_predictive_correlation
))
add_check(
  "Marker effects", "Fold-local centering translation invariance",
  "Every marker shifted by a different constant; installed fit and CV rerun",
  paste0("prediction/CV difference <= ", tol_identity), fmt(translation_error),
  translation_error <= tol_identity
)

reliability_semantics <- is.na(fit_ridge$reliability) &&
  identical(fit_ridge$reliability_is_calibrated, FALSE) &&
  is.finite(fit_ridge$cv_predictive_r2) && is.finite(fit_ridge$cv_predictive_correlation)
add_check(
  "Marker effects", "Predictive diagnostics are not mislabeled reliability",
  "Installed ridge fit metadata",
  "reliability=NA, calibrated flag FALSE, CV R2/correlation separately reported",
  paste0("reliability=", fit_ridge$reliability,
         "; calibrated=", fit_ridge$reliability_is_calibrated), reliability_semantics
)

fit_small <- ng_fit_ridge_effects(X[1:8, , drop = FALSE], y[1:8], lambda = 5, kfold = 5, seed = 77)
small_n_semantics <- is.na(fit_small$cv_predictive_r2) && is.na(fit_small$cv_predictive_correlation)
add_check(
  "Marker effects", "Small-sample CV is not replaced by in-sample fit",
  "Eight training records, below the package CV floor",
  "both CV diagnostics are NA", paste0("R2=", fit_small$cv_predictive_r2,
                                        "; cor=", fit_small$cv_predictive_correlation),
  small_n_semantics
)

uncalibrated <- fit_ridge
uncalibrated$reliability <- 0.99
uncalibrated$reliability_is_calibrated <- FALSE
chosen_uncalibrated <- ng_internal("ng_choose_mean_source")(
  X, uncalibrated, adjusted_pheno = y, ids = ridge_ids, min_reliability = 0.35
)
calibrated <- uncalibrated
calibrated$reliability_is_calibrated <- TRUE
chosen_calibrated <- ng_internal("ng_choose_mean_source")(
  X, calibrated, adjusted_pheno = y, ids = ridge_ids, min_reliability = 0.35
)
mean_gate_ok <- identical(chosen_uncalibrated$source, "adjusted_pheno") &&
  identical(chosen_calibrated$source, "GEBV")
add_check(
  "Marker effects", "Reliability gating requires explicit calibration",
  "Same numeric reliability with calibration flag FALSE then TRUE",
  "uncalibrated cannot gate; calibrated can gate",
  paste0(chosen_uncalibrated$source, " / ", chosen_calibrated$source), mean_gate_ok
)

# -----------------------------------------------------------------------------
# 4. Relationship scale and mate-allocation hard constraints
# -----------------------------------------------------------------------------

K_pair <- matrix(c(1, 0.4, 0.4, 1), 2L, 2L,
                 dimnames = list(c("A", "B"), c("A", "B")))
pair_scale <- ng_internal("ng_pair_relationship_variance")(
  data.frame(parent1 = "A", parent2 = "B"), K_pair
)
pair_scale_ok <- max_abs(c(pair_scale$pair_relationship - 0.4,
                           pair_scale$pair_kinship - 0.2)) <= tol_identity
add_check(
  "Relationships", "Additive relationship, kinship, and progeny-F scale",
  "Named relationship matrix with G[A,B]=0.4",
  "relationship=0.4; kinship=F=0.2",
  paste0("relationship=", fmt(pair_scale$pair_relationship),
         "; kinship=", fmt(pair_scale$pair_kinship)), pair_scale_ok
)

K_identity <- diag(4L)
dimnames(K_identity) <- list(LETTERS[1:4], LETTERS[1:4])
counts_identity <- stats::setNames(rep(1L, 4L), LETTERS[1:4])
group_relationship <- ng_internal("ng_group_relationship")(counts_identity, K_identity)
group_coancestry <- ng_internal("ng_group_coancestry")(counts_identity, K_identity)
group_scale_ok <- max_abs(c(group_relationship - 0.25, group_coancestry - 0.125,
                            1 / group_relationship - 4)) <= tol_identity
add_check(
  "Relationships", "Group relationship/coancestry/effective-size identity",
  "Four equally contributing unrelated parents with G=I",
  "c'Gc=0.25; group coancestry=0.125; Ne=4",
  paste0(fmt(group_relationship), " / ", fmt(group_coancestry), " / ",
         fmt(1 / group_relationship)), group_scale_ok
)

alloc_parents <- LETTERS[1:6]
K_alloc <- matrix(0.1, 6L, 6L, dimnames = list(alloc_parents, alloc_parents))
diag(K_alloc) <- 1
alloc_pairs <- t(utils::combn(alloc_parents, 2L))
alloc_scores <- data.frame(
  parent1 = alloc_pairs[, 1L], parent2 = alloc_pairs[, 2L],
  merit = seq(nrow(alloc_pairs), 1L), pair_relationship = 0.1,
  pair_kinship = 0.05, expected_progeny_inbreeding = 0.05,
  stringsAsFactors = FALSE
)
alloc_plan <- ng_optimize_mating_plan(
  alloc_scores, n_crosses = 6L, gain_col = "merit", parent_kinship = K_alloc,
  max_crosses_per_parent = 3L, min_unique_parents = 5L,
  method = "greedy_local", local_iter = 200L
)
alloc_summary <- attr(alloc_plan, "summary")
allocation_ok <- nrow(alloc_plan) == 6L &&
  isTRUE(alloc_summary$all_hard_constraints_satisfied) &&
  alloc_summary$max_parent_use <= 3L && alloc_summary$unique_parents >= 5L
add_check(
  "Mate allocation", "Final hard-constraint audit",
  "Installed greedy/local allocator with size, parent-use, and unique-parent constraints",
  "6 crosses; cap <=3; >=5 parents; final audit TRUE",
  paste0("n=", nrow(alloc_plan), "; cap=", alloc_summary$max_parent_use,
         "; parents=", alloc_summary$unique_parents,
         "; audit=", alloc_summary$all_hard_constraints_satisfied), allocation_ok
)

# -----------------------------------------------------------------------------
# 5. Formal multi-trait selection indices
# -----------------------------------------------------------------------------

P_index <- matrix(c(2.0, 0.35, 0.35, 1.4), 2L, 2L)
G_index <- matrix(c(1.2, 0.20, 0.20, 0.8), 2L, 2L)
a_index <- c(0.7, 0.3)
smith_info <- list(target_matrix = P_index, projection = G_index,
                   response_G = G_index, response_P = P_index)
smith_prod <- ng_internal("ng_multitrait_solve_index")(a_index, smith_info, ridge = 0)
smith_oracle <- as.numeric(solve(P_index, G_index %*% a_index))
smith_oracle <- smith_oracle / sum(abs(smith_oracle))
smith_error <- max_abs(smith_prod$coefficients - smith_oracle)
add_check(
  "Multi-trait", "Smith-Hazel coefficient identity",
  "Direct independent solve of b = P^-1 G a",
  paste0("max coefficient error <= ", tol_identity), fmt(smith_error),
  smith_error <= tol_identity
)

d_index <- c(0.8, 0.2)
pesek_info <- list(target_matrix = G_index, projection = NULL,
                   response_G = G_index, response_P = P_index)
pesek_prod <- ng_internal("ng_multitrait_solve_index")(d_index, pesek_info, ridge = 0)
pesek_oracle <- as.numeric(solve(G_index, d_index))
pesek_oracle <- pesek_oracle / sum(abs(pesek_oracle))
pesek_error <- max_abs(pesek_prod$coefficients - pesek_oracle)
add_check(
  "Multi-trait", "Pesek-Baker desired-gain identity",
  "Direct independent solve of b proportional to G^-1 d",
  paste0("max coefficient error <= ", tol_identity), fmt(pesek_error),
  pesek_error <= tol_identity
)

mt_scores <- data.frame(t1 = c(1, 2, 4, 7, 11), t2 = c(8, 6, 5, 3, 2))
economic_spec <- ng_multitrait_spec(
  trait = c("t1", "t2"), column = c("t1", "t2"), direction = c("maximize", "maximize"),
  economic_weight = c(1, 0)
)
economic_scored <- ng_add_multitrait_score(
  mt_scores, economic_spec, method = "economic_index",
  phenotypic_covariance = P_index, genetic_covariance = G_index
)
desired_spec <- ng_multitrait_spec(
  trait = c("t1", "t2"), column = c("t1", "t2"), direction = c("maximize", "maximize"),
  desired_change = c(1, 0), min_value = c(NA, 4)
)
desired_scored <- ng_add_multitrait_score(
  mt_scores, desired_spec, method = "desired_gain",
  phenotypic_covariance = P_index, genetic_covariance = G_index
)
economic_meta <- attr(economic_scored, "multi_trait")
desired_meta <- attr(desired_scored, "multi_trait")
zero_targets_ok <- economic_meta$economic_index_target[[2L]] == 0 &&
  desired_meta$desired_gain_target[[2L]] == 0 &&
  any(desired_scored$multi_trait_t2_violation > 0)
add_check(
  "Multi-trait", "Zero economic/desired target is preserved while thresholds remain active",
  "Production public scorer metadata and per-trait violation diagnostics",
  "both second-trait targets exactly zero and threshold violation detected",
  paste0("economic=", fmt(economic_meta$economic_index_target[[2L]]),
         "; desired=", fmt(desired_meta$desired_gain_target[[2L]]),
         "; threshold_active=", any(desired_scored$multi_trait_t2_violation > 0)),
  zero_targets_ok
)

formal_covariance_required <- expect_error(
  ng_add_multitrait_score(mt_scores, economic_spec, method = "economic_index")
) && expect_error(
  ng_add_multitrait_score(mt_scores, desired_spec, method = "desired_gain")
)
add_check(
  "Multi-trait", "Formal indices fail closed without P and G",
  "Public economic-index and desired-gain calls with covariance matrices omitted",
  "both calls error", formal_covariance_required, formal_covariance_required
)

# -----------------------------------------------------------------------------
# 6. Polyploid inheritance, GRM, and dominance kernel
# -----------------------------------------------------------------------------

n_poly <- 80L
m_poly <- 14L
poly_ids <- sprintf("D%03d", seq_len(n_poly))
poly_markers <- sprintf("d%02d", seq_len(m_poly))
M2 <- matrix(stats::rbinom(n_poly * m_poly, 2L, stats::runif(m_poly, 0.15, 0.85)),
             n_poly, dimnames = list(poly_ids, poly_markers))
Gd_prod <- ng_polyploid_dominance_grm(M2, ploidy = 2L, method = "vanraden",
                                      return_freq = TRUE, use_cpp = TRUE)
p2 <- attr(Gd_prod, "allele_freq")
W2 <- sweep(M2, 2L, 2 * p2, "-")
H2 <- M2 * (2 - M2)
hbar2 <- 2 * p2 * (1 - p2)
b2 <- 1 - 2 * p2
D2 <- sweep(H2, 2L, hbar2, "-") - sweep(W2, 2L, b2, "*")
support2 <- 0:2
varD2 <- vapply(seq_along(p2), function(j) {
  dj <- support2 * (2 - support2) - hbar2[[j]] - b2[[j]] * (support2 - 2 * p2[[j]])
  sum(stats::dbinom(support2, size = 2, prob = p2[[j]]) * dj^2)
}, numeric(1L))
Gd_oracle <- tcrossprod(D2) / sum(varD2)
dominance_grm_error <- max_abs(Gd_prod - Gd_oracle)
dominance_min_eigen <- min(eigen(Gd_prod, symmetric = TRUE, only.values = TRUE)$values)
add_check(
  "Polyploid", "Diploid reduction to Vitezica dominance coding",
  "Independent construction of D=(-2p^2,2pq,-2q^2) basis and GRM",
  paste0("max absolute error <= ", tol_identity), fmt(dominance_grm_error),
  dominance_grm_error <= tol_identity
)
add_check(
  "Polyploid", "Dominance GRM positive semidefinite",
  "Minimum eigenvalue of installed-kernel GRM",
  paste0("minimum eigenvalue >= -", tol_psd), fmt(dominance_min_eigen),
  dominance_min_eigen >= -tol_psd
)

pmf_errors <- numeric()
mean_errors <- numeric()
for (dosage in 0:4) {
  pmf <- ng_internal("ng_poly_gamete_pmf")(dosage, ploidy = 4L, dr = 1 / 6)
  pmf_errors <- c(pmf_errors, sum(pmf) - 1)
  mean_errors <- c(mean_errors, sum((0:2) * pmf) - dosage / 2)
}
moment4 <- ng_polyploid_progeny_moment_table(4L, double_reduction = 1 / 6)
midparent4 <- outer(0:4, 0:4, "+") / 2
poly_moment_error <- max_abs(c(pmf_errors, mean_errors, moment4$mu - midparent4))
add_check(
  "Polyploid", "Autotetraploid double-reduction probability and mean identities",
  "All parental dosages at conventional alpha=1/6",
  paste0("PMF sum, gamete mean, and progeny mean error <= ", tol_probability),
  fmt(poly_moment_error), poly_moment_error <= tol_probability
)

dr_domain_ok <- expect_error(ng_polyploid_progeny_moment_table(6L, double_reduction = 0.01)) &&
  expect_error(ng_polyploid_progeny_moment_table(4L, double_reduction = 0.20))
add_check(
  "Polyploid", "Double-reduction model domain fails closed",
  "Nonzero alpha at ploidy 6 and alpha above 1/6 at ploidy 4",
  "both calls error", dr_domain_ok, dr_domain_ok
)

n4 <- 16L
m4 <- 10L
ids4 <- sprintf("T%02d", seq_len(n4))
markers4 <- sprintf("t%02d", seq_len(m4))
M4 <- matrix(stats::rbinom(n4 * m4, 4L, stats::runif(m4, 0.2, 0.8)), n4,
             dimnames = list(ids4, markers4))
p4 <- colMeans(M4) / 4
W4 <- sweep(M4, 2L, 4 * p4, "-")
H4 <- M4 * (4 - M4)
hbar4 <- colMeans(H4)
Hc4 <- sweep(H4, 2L, hbar4, "-")
ss4 <- colSums(W4^2)
b4 <- ifelse(ss4 > 0, colSums(W4 * Hc4) / ss4, 0)
fit4 <- structure(list(
  beta_add = stats::rnorm(m4, 0, 0.15), beta_dom = stats::rnorm(m4, 0, 0.10),
  intercept = 3, allele_freq = p4, hbar = hbar4, b_orth = b4,
  markers = markers4, ploidy = 4L, model = "additive_dominance"
), class = "ng_polyploid_effects")
pairs4 <- data.frame(parent1 = ids4[1:8], parent2 = ids4[9:16], stringsAsFactors = FALSE)
poly_r <- ng_polyploid_score_crosses_dominance(
  fit4, M4, pairs = pairs4, double_reduction = 1 / 12, use_cpp = FALSE
)
poly_cpp <- ng_polyploid_score_crosses_dominance(
  fit4, M4, pairs = pairs4, double_reduction = 1 / 12, use_cpp = TRUE
)
poly_numeric <- c("cross_mean", "mid_parent_bv", "heterosis", "add_var", "dom_var", "cross_var")
poly_cpp_error <- max_abs(as.matrix(poly_r[poly_numeric]) - as.matrix(poly_cpp[poly_numeric]))
poly_label_ok <- all(poly_cpp$variance_model == "uniform_phase_prior_expectation")
add_check(
  "Polyploid", "Additive-dominance R/C++ kernel parity and phase label",
  "Same autotetraploid crosses through both production kernels",
  paste0("max absolute error <= ", tol_identity, " and dosage-only phase-prior label present"),
  paste0("error=", fmt(poly_cpp_error), "; label=", poly_label_ok),
  poly_cpp_error <= tol_identity && poly_label_ok
)

# -----------------------------------------------------------------------------
# 7. Probability metrics
# -----------------------------------------------------------------------------

mu_prob <- c(-1, 0.5, 2)
sd_prob <- c(0.3, 1.2, 2.5)
tau_prob <- c(0, 1, -0.5)
k_prob <- c(1L, 10L, 100L)
p_prod <- ng_p_superior_progeny(mu_prob, sd_prob, tau_prob, k_prob)
p_oracle <- 1 - stats::pnorm((tau_prob - mu_prob) / sd_prob)^k_prob
p_error <- max_abs(p_prod - p_oracle)
add_check(
  "Probability", "At-least-one superior progeny closed form",
  "Independent 1 - Phi((tau-mu)/sigma)^k calculation",
  paste0("max absolute error <= ", tol_probability), fmt(p_error),
  p_error <= tol_probability
)

zero_sd <- ng_p_superior_progeny(c(0, 2), c(0, 0), tau = 1, k_progeny = 20L)
zero_sd_ok <- identical(as.numeric(zero_sd), c(0, 1))
add_check(
  "Probability", "Zero-SD deterministic boundary",
  "Degenerate family distributions below and above threshold",
  "exactly c(0,1)", paste(zero_sd, collapse = ","), zero_sd_ok
)

ex_mu <- c(-0.2, 1.1, 3)
ex_sd <- c(0.4, 1.5, 0)
ex_tau <- c(0, 0.7, 2)
ex_prod <- ng_expected_excess_above_threshold(ex_mu, ex_sd, ex_tau)
z_ex <- (ex_tau[1:2] - ex_mu[1:2]) / ex_sd[1:2]
ex_oracle <- c(ex_sd[1:2] * stats::dnorm(z_ex) +
                 (ex_mu[1:2] - ex_tau[1:2]) * stats::pnorm(z_ex, lower.tail = FALSE),
               max(ex_mu[[3L]] - ex_tau[[3L]], 0))
ex_error <- max_abs(ex_prod - ex_oracle)
add_check(
  "Probability", "Expected excess above threshold identity",
  "Independent truncated-normal first moment plus deterministic limit",
  paste0("max absolute error <= ", tol_probability), fmt(ex_error),
  ex_error <= tol_probability
)

prob_domain_ok <- expect_error(ng_p_superior_progeny(0, -1, 0, 10L)) &&
  expect_error(ng_p_superior_progeny(0, 1, 0, 2.5)) &&
  expect_error(ng_expected_excess_above_threshold(0, -1, 0))
add_check(
  "Probability", "Invalid probability domains fail closed",
  "Negative SD and fractional progeny count",
  "all invalid calls error", prob_domain_ok, prob_domain_ok
)

# -----------------------------------------------------------------------------
# 8. Installed one-call production workflow with graph LD pruning + fast PMV
# -----------------------------------------------------------------------------

n_e2e <- 24L
m_e2e <- 60L
ids_e2e <- sprintf("E%03d", seq_len(n_e2e))
markers_e2e <- sprintf("e%03d", seq_len(m_e2e))
G_e2e <- matrix(2 * stats::rbinom(n_e2e * m_e2e, 1, 0.5), n_e2e,
                dimnames = list(ids_e2e, markers_e2e))
b_e2e <- stats::rnorm(m_e2e, 0, 0.12)
y_e2e <- as.numeric(50 + G_e2e %*% b_e2e + stats::rnorm(n_e2e, 0, 0.4))
geno_e2e <- data.frame(parent = ids_e2e, G_e2e, check.names = FALSE,
                       stringsAsFactors = FALSE)
pheno_e2e <- data.frame(parent = ids_e2e, yield = y_e2e, stringsAsFactors = FALSE)
map_e2e <- data.frame(
  marker = markers_e2e, chr = rep(1:3, each = 20L),
  pos_cm = rep(seq(0, 95, length.out = 20L), 3L), stringsAsFactors = FALSE
)
direction_e2e <- data.frame(trait = "yield", direction = "increase",
                            stringsAsFactors = FALSE)
e2e <- ng_run_cross_prediction(
  phenotype = pheno_e2e, genotype = geno_e2e, marker_map = map_e2e,
  trait_direction = direction_e2e,
  phenotype_id_col = "parent", genotype_id_col = "parent",
  map_marker_col = "marker", map_chr_col = "chr", map_pos_cm_col = "pos_cm",
  map_position_unit = "cM", progeny = "DH", parent_type = "inbred",
  trait_value_metric = "usefulness", uc_variance_source = "pmv",
  method_varPMV = "fast", ld_pruning = TRUE, ld_backend = "auto",
  ld_window = 8L, ld_r2_threshold = 0.95, ld_maf_threshold = 0.01,
  n_crosses = 6L, max_crosses_per_parent = 3L, min_unique_parents = 5L,
  optimizer = "greedy_local", lambda_group = 0.02, local_iter = 200L,
  duplicate_action = "none", use_cpp = TRUE, seed = 20260825
)
e2e_ok <- nrow(e2e$selected_crosses) == 6L &&
  all(is.finite(e2e$selected_crosses$multi_trait_score)) &&
  identical(e2e$settings$method_varPMV, "fast") &&
  isTRUE(e2e$ld_pruning_report$map_aware[[1L]]) &&
  isTRUE(e2e$plan_summary$all_hard_constraints_satisfied)
add_check(
  "End-to-end", "Installed production runner",
  "One-call DH workflow with compiled scoring, graph LD pruning, fast PMV, and allocation",
  "6 finite selected crosses; fast PMV; map-aware LD; hard audit TRUE",
  paste0("n=", nrow(e2e$selected_crosses),
         "; fast=", identical(e2e$settings$method_varPMV, "fast"),
         "; map_aware=", e2e$ld_pruning_report$map_aware[[1L]],
         "; audit=", e2e$plan_summary$all_hard_constraints_satisfied), e2e_ok
)

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------

results <- do.call(rbind, checks)
code_gate_pass <- all(results$status == "PASS")

csv_path <- Sys.getenv(
  "NGCD_STAT_GATE_CSV",
  unset = file.path("docs", "STATISTICAL_RELEASE_GATE_RESULTS.csv")
)
report_path <- Sys.getenv(
  "NGCD_STAT_GATE_REPORT",
  unset = file.path("docs", "STATISTICAL_RELEASE_GATE.md")
)

dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, csv_path, row.names = FALSE, na = "")

md_escape <- function(x) {
  x <- gsub("\\|", "\\\\|", as.character(x))
  gsub("\n", " ", x, fixed = TRUE)
}
table_lines <- c(
  "| Domain | Check | Evidence | Acceptance | Observed | Status |",
  "| --- | --- | --- | --- | --- | --- |",
  apply(results, 1L, function(row) paste0("| ", paste(md_escape(row), collapse = " | "), " |"))
)

lines <- c(
  "# Statistical Release Gate",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Package: `nextgenCrossDesign ", as.character(packageVersion("nextgenCrossDesign")), "`"),
  paste0("Installed namespace: `", pkg_path, "`"),
  paste0("Source Git commit: `", git_head, "` (working tree dirty: `", git_dirty, "`)"),
  paste0("R: `", R.version.string, "`"),
  paste0("Compiled LD backend available: `", compiled_ld_available, "`"),
  "",
  "## Decision",
  "",
  paste0("- Mathematical/software quantitative-genetics gate: **",
         if (code_gate_pass) "PASS" else "FAIL", "** (", sum(results$status == "PASS"),
         "/", nrow(results), " checks passed)."),
  "- Historical forward-validation gate: **NOT RUN**. No historical cross-by-progeny outcome table was supplied to this gate.",
  "- Unrestricted worldwide production release: **HOLD** until the historical gate and an independent quantitative-genetics review pass.",
  "",
  "A code-gate pass proves that the installed implementation satisfies the identities and invariants below within stated numerical tolerances. It does not prove prediction accuracy in every germplasm, crop, environment, generation, or breeding program.",
  "",
  "## Scope and method",
  "",
  "This runner loads an installed package and calls its production namespace, including compiled kernels. It does not source package test files. Mathematical identities use independent calculations; randomized graph checks compare exact retained-marker sets; the final check uses the public one-call workflow.",
  "",
  "Numerical identity checks use an absolute tolerance of `1e-10`; probability identities use `1e-12`; positive-semidefinite checks allow minimum eigenvalues down to `-1e-10` for floating-point roundoff. Discrete graph, domain, and hard-constraint checks require exact agreement.",
  "",
  "Graph LD pruning and fast PMV are retained. The gate explicitly checks them; neither is disabled.",
  "",
  "## Results",
  "",
  table_lines,
  "",
  "## Claims permitted by this gate",
  "",
  "- The installed diploid DH and infinite-selfed-RIL variance implementation matches the stated quantitative-genetic formulas for the checked models.",
  "- The fast Haldane-DH PMV recursion matches the dense quadratic form, and graph LD pruning matches its connected-component specification.",
  "- Relationship/coancestry scaling, formal Smith-Hazel and Pesek-Baker coefficients, probability metrics, and tested optimizer constraints are internally coherent.",
  "- The tested polyploid single-locus/dosage identities and additive-dominance R/C++ parity hold inside the explicitly reported model domain.",
  "",
  "## Claims not permitted by this gate",
  "",
  "- No claim of universal prediction accuracy, realized genetic gain, or superiority over external packages follows from algebraic or synthetic checks.",
  "- No empirical recommendation is authorized for a crop, target population, generation interval, training design, or environment not represented in a forward-validation dataset.",
  "- Polyploid additive-plus-dominance fitting remains experimental because additive and dominance components share one ridge penalty.",
  "- Dosage-only autopolyploid within-family variance remains a uniform compatible-phase prior expectation; exact linkage-phase claims require phased homologues.",
  "- UCPC remains experimental/internal until calibrated against independent breeding outcomes.",
  "",
  "## Required historical gate",
  "",
  "For each breeding population, train only on records available before the target cycle, predict candidate crosses, and compare predictions with subsequently observed progeny. At minimum report cross-mean calibration intercept/slope, RMSE, Pearson and Spearman correlation, within-family variance calibration, top-decile enrichment, selected-vs-control realized response, inbreeding/coancestry, family sizes, uncertainty intervals, and results by cycle/environment. Pre-register thresholds with breeders before examining outcomes; do not invent universal cutoffs after seeing the data.",
  "",
  "## Reproduce",
  "",
  "```sh",
  "mkdir -p /tmp/ngcd_release_gate_lib",
  "R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib .",
  "NGCD_RELEASE_LIB=/tmp/ngcd_release_gate_lib Rscript tools/run_statistical_release_gate.R",
  "```",
  "",
  paste0("Machine-readable results: `", csv_path, "`")
)
writeLines(lines, report_path, useBytes = TRUE)

message("Statistical code gate: ", if (code_gate_pass) "PASS" else "FAIL",
        " (", sum(results$status == "PASS"), "/", nrow(results), ")")
message("Wrote: ", csv_path)
message("Wrote: ", report_path)

if (!code_gate_pass) quit(save = "no", status = 1L)
