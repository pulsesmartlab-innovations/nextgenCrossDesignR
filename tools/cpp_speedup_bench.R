# Microbenchmark C++ vs R for the three v0.3.0 ports on a realistic fixture.
ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"),
            file.path("nextgen_cross_design", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

bench_results <- list()
bench <- function(label, expr_cpp, expr_r) {
  t_cpp <- system.time(replicate(3L, force(expr_cpp())))["elapsed"] / 3
  t_r   <- system.time(replicate(3L, force(expr_r())))["elapsed"] / 3
  cat(sprintf("%-32s  cpp=%.3fs  r=%.3fs  speedup=%.1fx\n",
              label, t_cpp, t_r, t_r / max(t_cpp, 1e-6)))
  bench_results[[label]] <<- list(cpp = unname(t_cpp), r = unname(t_r))
  invisible(list(cpp = t_cpp, r = t_r))
}

set.seed(20260524L)
n <- 60L; m <- 1500L
ids <- paste0("P", seq_len(n))
markers <- paste0("M", seq_len(m))
geno <- 2L * matrix(rbinom(n * m, 1L, 0.45), n, m, dimnames = list(ids, markers))
y <- as.numeric(geno %*% rnorm(m, sd = 0.08) + rnorm(n, sd = 0.5))
names(y) <- ids
mk <- data.frame(marker = markers, chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 100, length.out = m / 5L), 5L))
fit <- ng_fit_ridge_effects(geno, y, ids = ids, kfold = 5L, seed = 1L,
                            return_beta_cov_full = TRUE)
pairs <- ng_make_pairs(ids)
mm <- ng_prepare_marker_map(mk, marker_ids = markers, model = "haldane")
sorted <- ng_sort_by_map(geno, fit$beta, fit$beta_var, mm)
mm_k <- ng_prepare_marker_map(mk, marker_ids = markers, model = "kosambi")
sorted_k <- ng_sort_by_map(geno, fit$beta, fit$beta_var, mm_k)

cat(sprintf("# Microbenchmark: n=%d, m=%d, n_pairs=%d\n\n",
            n, m, nrow(pairs)))

# ---- 1. Full off-diagonal posterior PMV -----------------------------------
bench("full-posterior PMV",
      function() ng_dh_recomb_variance_pairs_full_posterior(
        geno = sorted$geno, beta = sorted$effects,
        beta_cov_full = fit$beta_cov_full[sorted$marker_map$marker, sorted$marker_map$marker],
        marker_map = sorted$marker_map, ids = ids, pairs = pairs,
        target = "DH", recomb_model = "haldane", use_cpp = TRUE),
      function() ng_dh_recomb_variance_pairs_full_posterior(
        geno = sorted$geno, beta = sorted$effects,
        beta_cov_full = fit$beta_cov_full[sorted$marker_map$marker, sorted$marker_map$marker],
        marker_map = sorted$marker_map, ids = ids, pairs = pairs,
        target = "DH", recomb_model = "haldane", use_cpp = FALSE))

# ---- 2. Banded Kosambi/RIL kernel ------------------------------------------
bench("banded Kosambi (window=20cM)",
      function() ng_dh_recomb_variance_pairs_banded(
        geno = sorted_k$geno, beta = sorted_k$effects, beta_var = sorted_k$beta_var,
        marker_map = sorted_k$marker_map, ids = ids, pairs = pairs,
        target = "DH", recomb_model = "kosambi", window_cm = 20,
        use_cpp = TRUE),
      function() ng_dh_recomb_variance_pairs_banded(
        geno = sorted_k$geno, beta = sorted_k$effects, beta_var = sorted_k$beta_var,
        marker_map = sorted_k$marker_map, ids = ids, pairs = pairs,
        target = "DH", recomb_model = "kosambi", window_cm = 20,
        use_cpp = FALSE))

# ---- 3. BCM posterior sampler ---------------------------------------------
Xc <- sweep(geno, 2L, fit$marker_mean, "-")
storage.mode(Xc) <- "double"
yc <- y - mean(y)
bench("BCM sampler (200 draws)",
      function() ng_sample_ridge_posterior_bcm(
        X = Xc, yc = yc, sigma_e2 = fit$sigma_e2, lambda = fit$lambda,
        n_draws = 200L, seed = 1L, use_cpp = TRUE),
      function() ng_sample_ridge_posterior_bcm(
        X = Xc, yc = yc, sigma_e2 = fit$sigma_e2, lambda = fit$lambda,
        n_draws = 200L, seed = 1L, use_cpp = FALSE))

# ---- 4. ng_local_swap greedy OCS optimizer -------------------------------
scores_b <- ng_score_crosses(geno = geno, effects = fit, marker_map = mk, ids = ids,
                             adjusted_pheno = setNames(rnorm(n), ids),
                             selection_prop = 0.10, recomb_model = "haldane",
                             target = "DH", use_cpp = TRUE)
scores_b$.linear_gain <- scores_b$uc_dh_gebv
parent_K <- ng_parent_kinship(geno)
ord_x <- order(scores_b$.linear_gain, decreasing = TRUE)
sel_init <- integer(0); counts0 <- setNames(integer(n), ids)
for (idx in ord_x) {
  pp <- c(scores_b$parent1[idx], scores_b$parent2[idx])
  if (any(counts0[pp] >= 4L)) next
  sel_init <- c(sel_init, idx); counts0[pp] <- counts0[pp] + 1L
  if (length(sel_init) == 20L) break
}
bench("local_swap (n_crosses=20, iter=2000)",
      function() ng_local_swap(scores = scores_b, selected = sel_init,
                               parents = ids, parent_K = parent_K,
                               max_crosses_per_parent = 4L, lambda_group = 1,
                               local_iter = 2000L, use_cpp = TRUE),
      function() ng_local_swap(scores = scores_b, selected = sel_init,
                               parents = ids, parent_K = parent_K,
                               max_crosses_per_parent = 4L, lambda_group = 1,
                               local_iter = 2000L, use_cpp = FALSE))

# ---- Emit machine-readable CSV for downstream plotting --------------------
# Find the package root (the dir containing R/load.R) and write the CSV to
# its sibling `results/` directory, matching head-to-head/multitrait runners.
.bench_pkg_root <- (function() {
  start <- normalizePath(".", winslash = "/", mustWork = TRUE)
  for (cand in c(start, dirname(start), file.path(start, "nextgen_cross_design"),
                 file.path(dirname(start), "nextgen_cross_design"))) {
    if (file.exists(file.path(cand, "R", "load.R")))
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate nextgenCrossDesign package root for bench CSV output.")
})()
results_dir <- normalizePath(file.path(.bench_pkg_root, "..", "results"),
                             winslash = "/", mustWork = FALSE)
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
bench_df <- do.call(rbind, lapply(names(bench_results), function(label) {
  data.frame(kernel = label,
             r_seconds = bench_results[[label]]$r,
             cpp_seconds = bench_results[[label]]$cpp,
             speedup_x = bench_results[[label]]$r / max(bench_results[[label]]$cpp, 1e-9),
             n = n, m = m, n_pairs = nrow(pairs),
             stringsAsFactors = FALSE)
}))
csv_path <- file.path(results_dir, "cpp_speedup_bench.csv")
write.csv(bench_df, csv_path, row.names = FALSE)
cat(sprintf("\nWrote %s (%d kernels)\n", csv_path, nrow(bench_df)))
