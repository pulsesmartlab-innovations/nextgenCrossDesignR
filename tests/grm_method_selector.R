# The GRM-method selector (vanraden default / yang) must flow through the additive/diploid
# entry points ng_design_crosses() and ng_run_cross_prediction(), mirroring the polyploid API
# (ng_design_crosses_poly(grm_method=)). Switching the GRM must change the parent kinship the
# allocator sees (pair_kinship), and the default must stay VanRaden (backward compatible).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(909)
n <- 30L; m <- 120L
ids <- sprintf("P%02d", seq_len(n))
# segregating inbred lines (0/2) with a spread of allele frequencies so VanRaden's single
# overall Sum(2p(1-p)) scaling and Yang's per-marker unit-variance weighting genuinely diverge.
freq <- runif(m, 0.1, 0.9)
gm <- vapply(freq, function(p) 2L * rbinom(n, 1, p), integer(n))
rownames(gm) <- ids; colnames(gm) <- sprintf("M%03d", seq_len(m))
beta <- rnorm(m, 0, 0.1)
y <- as.numeric(gm %*% beta) + rnorm(n)
names(y) <- ids
mm <- data.frame(marker = colnames(gm), chr = rep(1:4, length.out = m),
                 pos_cm = rep(seq(0, 100, length.out = ceiling(m / 4)), 4)[seq_len(m)],
                 stringsAsFactors = FALSE)

# --- unit level: the two GRMs differ, and the default is VanRaden -----------------------
Kv <- ng_parent_kinship(gm, method = "vanraden")
Ky <- ng_parent_kinship(gm, method = "yang")
stopifnot(!isTRUE(all.equal(Kv, Ky)))
stopifnot(isTRUE(all.equal(ng_parent_kinship(gm), Kv)))

# --- 1. grm_method flows through the documented pipeline ng_design_crosses() ------------
nx <- 12L
base <- list(geno = gm, y = y, marker_map = mm, ids = ids, n_crosses = nx,
             max_crosses_per_parent = 6L, use_cpp = FALSE)
d_def <- do.call(ng_design_crosses, base)
d_v   <- do.call(ng_design_crosses, c(base, list(grm_method = "vanraden")))
d_y   <- do.call(ng_design_crosses, c(base, list(grm_method = "yang")))
stopifnot(nrow(d_y$plan) == nx)
stopifnot(isTRUE(all.equal(d_def$scores$pair_kinship, d_v$scores$pair_kinship)))   # default == vanraden
stopifnot(!isTRUE(all.equal(d_v$scores$pair_kinship, d_y$scores$pair_kinship)))    # yang changes kinship

# --- 2. grm_method flows through the user-facing runner ng_run_cross_prediction() -------
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE)
runmm     <- data.frame(SNP_code = colnames(gm), Chromosome = mm$chr,
                        Position_BP = mm$pos_cm * 1e6, stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)
run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP_code", map_chr_col = "Chromosome",
  map_pos_col = "Position_BP", map_pos_cm_divisor = 1e6,
  prediction_mode = "trait_by_trait", trait_value_metric = "uc",
  duplicate_action = "none", n_crosses = 8L, max_uses_per_parent = 4L,
  use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L, ...)
r_v <- run(grm_method = "vanraden")
r_y <- run(grm_method = "yang")
stopifnot(inherits(r_y, "ng_cross_prediction_result"))
stopifnot(!isTRUE(all.equal(r_v$candidate_crosses$pair_kinship,
                            r_y$candidate_crosses$pair_kinship)))

cat("grm_method selector test passed\n")
