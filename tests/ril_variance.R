helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# T1.2: target = "RIL" must use the Haldane-Waddington (1931) equilibrium
# recombination fraction R = 2r / (1 + 2r), so the covariance kernel is
# (1 - 2 R) = (1 - 2 r) / (1 + 2 r). Before the fix, RIL was treated as DH.

# ---- Map-function sanity: RIL equilibrium decay --------------------------------------------
d <- 50
r_h <- ng_meiosis_r(d, model = "haldane")
expect_dh_h <- 1 - 2 * r_h
expect_ril_h <- (1 - 2 * r_h) / (1 + 2 * r_h)
stopifnot(abs(ng_progeny_decay(d, "haldane", "DH")  - expect_dh_h)  < 1e-12)
stopifnot(abs(ng_progeny_decay(d, "haldane", "RIL") - expect_ril_h) < 1e-12)
stopifnot(expect_dh_h > expect_ril_h)  # RIL kernel decays faster (more meioses)

# Kosambi RIL closed-form: tanh identity gives (1 - tanh(x)) / (1 + tanh(x)) = e^{-2x}.
# With x = 2 d / 100 that means decay_kosambi_RIL = exp(-4 d / 100).
expect_kos_ril <- exp(-4 * d / 100)
stopifnot(abs(ng_progeny_decay(d, "kosambi", "RIL") - expect_kos_ril) < 1e-12)

# ---- Two-locus closed form: variance of additive value under RIL --------------------------
# Two pure inbred parents AABB x aabb, single chromosome, 50 cM, beta_1 = beta_2 = 1.
# d_k = 0.5 * (2 - 0) = 1 for both markers, so a = (1, 1).
# Var(RIL_ij) = a_1^2 + a_2^2 + 2 a_1 a_2 (1 - 2 R) = 2 + 2 * decay_RIL.
geno <- rbind(P1 = c(2, 2), P2 = c(0, 0))
colnames(geno) <- c("M1", "M2")
beta <- setNames(c(1.0, 1.0), c("M1", "M2"))
beta_var <- setNames(c(0.0, 0.0), c("M1", "M2"))
marker_map <- ng_prepare_marker_map(
  data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, d)),
  marker_ids = c("M1", "M2"),
  model = "haldane"
)
pairs <- data.frame(parent1 = "P1", parent2 = "P2", stringsAsFactors = FALSE)

ril_var <- ng_dh_recomb_variance_pairs(
  geno = geno, beta = beta, beta_var = beta_var,
  marker_map = marker_map, ids = rownames(geno), pairs = pairs,
  use_cpp = FALSE, recomb_model = "haldane", target = "RIL"
)
dh_var <- ng_dh_recomb_variance_pairs(
  geno = geno, beta = beta, beta_var = beta_var,
  marker_map = marker_map, ids = rownames(geno), pairs = pairs,
  use_cpp = FALSE, recomb_model = "haldane", target = "DH"
)
expect_ril_var <- 2 + 2 * expect_ril_h
expect_dh_var  <- 2 + 2 * expect_dh_h
stopifnot(abs(ril_var$dh_recomb_var - expect_ril_var) < 1e-10)
stopifnot(abs(dh_var$dh_recomb_var  - expect_dh_var)  < 1e-10)
stopifnot(ril_var$dh_recomb_var < dh_var$dh_recomb_var)  # RIL <= DH at same locus

# ---- ng_score_crosses with target = "RIL" must differ from target = "DH" -----------------
set.seed(7L)
n <- 8L; m <- 30L
g <- matrix(2L * rbinom(n * m, 1, 0.5), nrow = n)
ids2 <- paste0("P", seq_len(n))
markers2 <- paste0("M", seq_len(m))
rownames(g) <- ids2; colnames(g) <- markers2
beta2 <- setNames(rnorm(m, sd = 0.15), markers2)
beta_var2 <- setNames(rep(0.005, m), markers2)
mm2 <- data.frame(
  marker = markers2,
  chr = rep(1:3, length.out = m),
  pos_cm = rep(seq(0, 80, length.out = m / 3), 3)
)
adj <- setNames(rnorm(n), ids2)
effects <- list(beta = beta2, beta_var = beta_var2, reliability = 0.5,
                intercept = 0, marker_mean = colMeans(g))

sc_dh  <- ng_score_crosses(geno = g, effects = effects, marker_map = mm2,
                           ids = ids2, adjusted_pheno = adj,
                           selection_prop = 0.1, target = "DH",
                           recomb_model = "haldane", use_cpp = FALSE)
sc_ril <- ng_score_crosses(geno = g, effects = effects, marker_map = mm2,
                           ids = ids2, adjusted_pheno = adj,
                           selection_prop = 0.1, target = "RIL",
                           recomb_model = "haldane", use_cpp = FALSE)
delta <- abs(sc_dh$dh_recomb_var - sc_ril$dh_recomb_var)
if (max(delta) < 1e-6) {
  stop("ng_score_crosses produced identical DH and RIL recombination variance")
}
# NOTE: we do NOT assert RIL <= DH cross-by-cross. With pairs of effects in
# repulsion (a_k a_l < 0), higher recombination INCREASES the within-family
# variance, so RIL can exceed DH at the individual-cross level. The mean
# direction across many crosses depends on the empirical distribution of
# coupling vs repulsion.

cat("ril_variance: 4/4 checks passed\n")
cat(sprintf("  decay at 50 cM   haldane-DH=%.4f  haldane-RIL=%.4f  kosambi-RIL=%.4f\n",
            expect_dh_h, expect_ril_h, expect_kos_ril))
cat(sprintf("  P1xP2 var        DH=%.4f  RIL=%.4f\n",
            dh_var$dh_recomb_var, ril_var$dh_recomb_var))
cat(sprintf("  mean RIL/DH ratio across %d crosses: %.3f\n",
            nrow(sc_dh), mean(sc_ril$dh_recomb_var / pmax(sc_dh$dh_recomb_var, 1e-12))))
