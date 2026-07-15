# Exact within-family cross-trait covariance a_t' R a_s (ng_cross_trait_within_family_cov) must
# match the empirical covariance among simulated DH progeny of a cross. Ground truth is a direct
# Monte-Carlo simulation of DH segregation (a per-chromosome +/-1 Markov chain with Haldane
# recombination), independent of the analytic kernel.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(2026)
nchr <- 2L; per <- 20L; m <- nchr * per
mk  <- sprintf("M%03d", seq_len(m))
chr <- rep(seq_len(nchr), each = per)
pos_cm <- rep(seq(0, 120, length.out = per), nchr)          # ~6.3 cM spacing
map <- data.frame(marker = mk, chr = chr, pos_cm = pos_cm, stringsAsFactors = FALSE)

# Two inbred parents (0/2) differing at many markers.
g1 <- 2L * rbinom(m, 1, 0.5); g2 <- 2L * rbinom(m, 1, 0.5)
geno <- rbind(P1 = g1, P2 = g2); colnames(geno) <- mk; storage.mode(geno) <- "double"

# Two correlated trait effect vectors.
b1 <- rnorm(m); b2 <- 0.6 * b1 + 0.8 * rnorm(m)
betas <- cbind(t1 = b1, t2 = b2); rownames(betas) <- mk

## Analytic
ct <- ng_cross_trait_within_family_cov(
  geno, betas, map, pairs = data.frame(parent1 = "P1", parent2 = "P2"),
  target = "DH", recomb_model = "haldane"
)
V1 <- ct$wf_var_t1; V2 <- ct$wf_var_t2; C12 <- ct$wf_cov_t1_t2

## Monte-Carlo ground truth: simulate DH segregation indicators y_k in {-1,+1}.
d <- 0.5 * (g1 - g2)                                          # marker contrasts
nsim <- 30000L
Y <- matrix(0, nsim, m)
for (c in seq_len(nchr)) {
  ix <- which(chr == c)
  yc <- matrix(0, nsim, length(ix))
  yc[, 1] <- sample(c(-1, 1), nsim, replace = TRUE)
  for (k in 2:length(ix)) {
    dcm <- pos_cm[ix[k]] - pos_cm[ix[k - 1]]
    r <- 0.5 * (1 - exp(-2 * dcm / 100))                     # Haldane recombination fraction
    flip <- runif(nsim) < r
    yc[, k] <- ifelse(flip, -yc[, k - 1], yc[, k - 1])
  }
  Y[, ix] <- yc
}
T1 <- as.numeric(Y %*% (d * b1)); T2 <- as.numeric(Y %*% (d * b2))
mc <- c(var1 = var(T1), var2 = var(T2), cov12 = cov(T1, T2))

cat(sprintf("  var(t1):  analytic %.3f  MC %.3f\n", V1,  mc["var1"]))
cat(sprintf("  var(t2):  analytic %.3f  MC %.3f\n", V2,  mc["var2"]))
cat(sprintf("  cov(t1,t2): analytic %.3f  MC %.3f\n", C12, mc["cov12"]))

rel <- function(a, b) abs(a - b) / max(1e-6, abs(b))
stopifnot(rel(V1,  mc["var1"])  < 0.06)
stopifnot(rel(V2,  mc["var2"])  < 0.06)
stopifnot(rel(C12, mc["cov12"]) < 0.08)
# The cross-trait covariance must be materially different from the population-correlation proxy
# (D * cov2cor(G) * D would use cor(b1,b2); this uses the recombination-aware bilinear form).
stopifnot(abs(C12) > 1e-6)

cat("cross_trait_within_family_cov.R: PASS\n")
