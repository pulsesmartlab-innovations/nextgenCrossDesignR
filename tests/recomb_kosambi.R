helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# T1.1: recomb_model = "kosambi" must produce a different DH variance than "haldane".
# Before the fix the R recursion and the C++ kernel both hardcoded exp(-2 d / 100)
# (Haldane) regardless of the requested mapping function.

# ---- Map function sanity ------------------------------------------------------------------
stopifnot(abs(ng_haldane_decay(0)   - 1) < 1e-12)
stopifnot(abs(ng_kosambi_decay(0)   - 1) < 1e-12)
stopifnot(abs(ng_haldane_decay(50)  - exp(-1)) < 1e-12)
stopifnot(abs(ng_kosambi_decay(50)  - (1 - tanh(1))) < 1e-12)
stopifnot(ng_haldane_decay(50) != ng_kosambi_decay(50))

# ---- Two-locus closed-form check ---------------------------------------------------------
# Two inbred parents AA-BB vs aa-bb, single chromosome, 50 cM apart, equal positive effects.
# DH variance from F1 of two pure inbreds:
#   Var = a_1^2 + a_2^2 + 2 a_1 a_2 (1 - 2 r)
# with a_k = 0.5 * (x_1k - x_2k) * beta_k = (+/-1) * beta_k for differing inbreds.
geno <- rbind(
  P1 = c(2, 2),  # AABB
  P2 = c(0, 0)   # aabb
)
colnames(geno) <- c("M1", "M2")
beta <- c(1.0, 1.0)
names(beta) <- c("M1", "M2")
beta_var <- c(0.0, 0.0)
names(beta_var) <- c("M1", "M2")
marker_map_h <- ng_prepare_marker_map(
  data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, 50)),
  marker_ids = c("M1", "M2"),
  model = "haldane"
)
marker_map_k <- ng_prepare_marker_map(
  data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, 50)),
  marker_ids = c("M1", "M2"),
  model = "kosambi"
)
pairs <- data.frame(parent1 = "P1", parent2 = "P2", stringsAsFactors = FALSE)

dh_h <- ng_dh_recomb_variance_pairs(
  geno = geno, beta = beta, beta_var = beta_var,
  marker_map = marker_map_h, ids = rownames(geno), pairs = pairs,
  use_cpp = FALSE, recomb_model = "haldane"
)
dh_k <- ng_dh_recomb_variance_pairs(
  geno = geno, beta = beta, beta_var = beta_var,
  marker_map = marker_map_k, ids = rownames(geno), pairs = pairs,
  use_cpp = FALSE, recomb_model = "kosambi"
)

# Analytic: a_1 = a_2 = -1 (since 0.5*(2-0)=1 then times beta=1 ... wait sign:
# d = 0.5*(geno[P1] - geno[P2]) = 0.5*(2 - 0) = +1, so a = (1, 1).
# Var = 1 + 1 + 2*1*1*(1 - 2 r) = 2 + 2*(1 - 2r) = 4 - 4r.
# Equivalently using the package's 1 - 2r decay: 2 + 2 * decay(50).
expect_h <- 2 + 2 * ng_haldane_decay(50)
expect_k <- 2 + 2 * ng_kosambi_decay(50)
stopifnot(abs(dh_h$vpm - expect_h) < 1e-10)
stopifnot(abs(dh_k$vpm - expect_k) < 1e-10)
stopifnot(abs(dh_h$vpm - dh_k$vpm) > 0.05)

# ---- Dense Haldane path must agree with the chromosome recursion (sanity) ----------------
# Force the dense path for a Haldane map and compare to the recursion.
set.seed(42L)
n  <- 8L
m  <- 20L
geno2 <- matrix(2L * (rbinom(n * m, 1, 0.5)), nrow = n)
ids2 <- paste0("P", seq_len(n))
markers2 <- paste0("M", seq_len(m))
rownames(geno2) <- ids2; colnames(geno2) <- markers2
beta2 <- rnorm(m, sd = 0.1); names(beta2) <- markers2
beta_var2 <- rep(0.01, m); names(beta_var2) <- markers2
mm <- ng_prepare_marker_map(
  data.frame(marker = markers2,
             chr = rep(1:2, length.out = m),
             pos_cm = rep(seq(0, 100, length.out = m / 2), 2)),
  marker_ids = markers2, model = "haldane"
)
# Sort by chr/pos as ng_score_crosses does before calling the recursion path; the
# recursion assumes markers are in genomic order along each chromosome.
sorted <- ng_sort_by_map(geno2, beta2, beta_var2, mm)
pairs2 <- ng_make_pairs(ids2)

recursion <- ng_dh_recomb_variance_pairs_r(
  geno = sorted$geno, beta = sorted$effects, beta_var = sorted$beta_var,
  marker_map = sorted$marker_map, ids = ids2, pairs = pairs2
)
dense_h <- ng_dh_recomb_variance_pairs_dense(
  geno = sorted$geno, beta = sorted$effects, beta_var = sorted$beta_var,
  marker_map = sorted$marker_map, ids = ids2, pairs = pairs2,
  recomb_model = "haldane"
)
max_var_err <- max(abs(recursion$vpm - dense_h$vpm))
max_pmv_err <- max(abs(recursion$pmv   - dense_h$pmv))
if (max_var_err > 1e-9 || max_pmv_err > 1e-9) {
  stop(sprintf("dense Haldane path disagrees with recursion: var err=%g, pmv err=%g",
               max_var_err, max_pmv_err))
}

# ---- ng_score_crosses must honor the requested mapping function ---------------------------
sc_h <- ng_score_crosses(
  geno = geno2, effects = list(beta = beta2, beta_var = beta_var2,
                               reliability = 0.5, intercept = 0,
                               marker_mean = colMeans(geno2)),
  marker_map = data.frame(marker = markers2,
                          chr = rep(1:2, length.out = m),
                          pos_cm = rep(seq(0, 100, length.out = m / 2), 2)),
  ids = ids2, adjusted_pheno = setNames(rnorm(n), ids2),
  selection_prop = 0.1, recomb_model = "haldane",
  use_cpp = FALSE
)
sc_k <- ng_score_crosses(
  geno = geno2, effects = list(beta = beta2, beta_var = beta_var2,
                               reliability = 0.5, intercept = 0,
                               marker_mean = colMeans(geno2)),
  marker_map = data.frame(marker = markers2,
                          chr = rep(1:2, length.out = m),
                          pos_cm = rep(seq(0, 100, length.out = m / 2), 2)),
  ids = ids2, adjusted_pheno = setNames(rnorm(n), ids2),
  selection_prop = 0.1, recomb_model = "kosambi",
  use_cpp = FALSE
)
if (max(abs(sc_h$vpm - sc_k$vpm)) < 1e-6) {
  stop("ng_score_crosses produced identical Haldane and Kosambi DH variance")
}

cat("recomb_kosambi: 4/4 checks passed\n")
cat(sprintf("  decay at 50 cM   haldane=%.4f  kosambi=%.4f\n",
            ng_haldane_decay(50), ng_kosambi_decay(50)))
cat(sprintf("  P1xP2 var        haldane=%.4f  kosambi=%.4f\n",
            dh_h$vpm, dh_k$vpm))
