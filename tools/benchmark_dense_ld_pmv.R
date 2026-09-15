.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))

set.seed(11)
n <- as.integer(Sys.getenv("BENCH_N_PARENTS", "80"))
m <- as.integer(Sys.getenv("BENCH_N_MARKERS", "5000"))
window <- as.integer(Sys.getenv("BENCH_LD_WINDOW", "250"))
n_chr <- as.integer(Sys.getenv("BENCH_N_CHR", "5"))

ids <- paste0("P", seq_len(n))
geno <- matrix(sample(c(0, 1, 2), n * m, replace = TRUE), n, m)
rownames(geno) <- ids
colnames(geno) <- paste0("m", seq_len(m))
chr <- rep(seq_len(n_chr), length.out = m)

effects <- matrix(rnorm(m, 0, 0.05), 1, m)
colnames(effects) <- colnames(geno)
mom <- effects_moments(effects)

ld_time <- system.time({
  ld <- calculate_ld_matrix_efficient(
    geno,
    method = "cor",
    window_size = window,
    chr = chr,
    verbose = FALSE
  )
})

edge_time <- system.time({
  edge <- calc_pmv_cross_ld_edges_cpp(
    geno_mat = geno,
    effect_mu = mom$mu,
    effect_var = rep(0, m),
    ld_mat = ld,
    candidate_ids = ids,
    progeny = "DH"
  )
})

batch_time <- system.time({
  batch <- calc_progeny_variance(
    method = "pmv_nopos",
    marker_effects_mcmc = effects,
    geno_mat = geno,
    ld_mat = ld,
    candidate_ids = ids,
    use_parallel = FALSE,
    progeny = "DH",
    include_effect_uncertainty = FALSE
  )
})

cat("ld_cache\n")
print(ld_time)
cat("edge_cpp\n")
print(edge_time)
cat("batch_current\n")
print(batch_time)
cat("maxdiff ", max(abs(edge$cross_var - batch$cross_var)), "\n")
