# Exact within-family cross-trait covariance for the multi-trait threshold probability.
#
# The multi-trait "probability a progeny clears the thresholds" needs the WITHIN-FAMILY covariance
# of the traits among a cross's progeny. By default it is approximated with per-trait within-family
# variances and a POPULATION genetic correlation. `ng_cross_trait_within_family_cov()` instead
# computes the EXACT recombination-aware covariance a_t' R a_s per cross (the two-trait
# generalization of the single-trait a'Ra; its diagonal reproduces vpm). Pass it to
# `ng_add_p_superior_progeny_multitrait(cross_trait_cov = ...)` to build an exact Sigma_c.
#
# This script is self-contained (small simulated inbred lines) and runs against the installed
# package.

library(nextgenCrossDesign)

set.seed(20260706)

# --- Small inbred (DH) parent panel: dosages 0/2 ---
n_par <- 12L; n_chr <- 3L; per_chr <- 30L; m <- n_chr * per_chr
ids <- sprintf("L%02d", seq_len(n_par))
mk  <- sprintf("SNP%03d", seq_len(m))
geno <- matrix(2L * rbinom(n_par * m, 1, 0.5), n_par, m, dimnames = list(ids, mk))
storage.mode(geno) <- "double"

marker_map <- data.frame(
  marker = mk,
  chr    = rep(seq_len(n_chr), each = per_chr),
  pos_cm = rep(seq(0, 120, length.out = per_chr), n_chr)   # ~4 cM spacing
)

# --- Two traits with correlated genetic architectures (yield up, disease down) ---
b_yield   <- rnorm(m)
b_disease <- 0.5 * b_yield + rnorm(m)                       # genetically correlated with yield
pheno <- data.frame(
  yield   = 60 + 5 * as.numeric(scale(geno %*% b_yield))   + rnorm(n_par),
  disease =  3 + 1 * as.numeric(scale(geno %*% b_disease)) + rnorm(n_par),
  row.names = ids
)

# --- Per-trait marker effects, one column per trait (aligned to geno markers) ---
betas <- sapply(c("yield", "disease"), function(tr)
  ng_fit_ridge_effects(geno, setNames(pheno[[tr]], ids), ids)$beta)
stopifnot(nrow(betas) == m)

# --- All candidate crosses among the parents ---
pairs <- ng_make_pairs(ids)

# --- EXACT within-family cross-trait covariance per cross ---
cross_cov <- ng_cross_trait_within_family_cov(
  geno, betas, marker_map, ids = ids, pairs = pairs, target = "DH", recomb_model = "haldane"
)
# diagonal = within-family variance; off-diagonal = recombination-aware cross-trait covariance
head(cross_cov)

# --- Assemble a candidate-cross score table (mid-parent means + within-family variances) ---
gebv_yield   <- as.numeric(geno %*% b_yield)
gebv_disease <- as.numeric(geno %*% b_disease)
names(gebv_yield) <- names(gebv_disease) <- ids
p1 <- pairs$parent1; p2 <- pairs$parent2
scores <- data.frame(
  parent1 = p1, parent2 = p2,
  pred_yield   = 0.5 * (gebv_yield[p1]   + gebv_yield[p2]),
  pred_disease = 0.5 * (gebv_disease[p1] + gebv_disease[p2]),
  var_yield    = cross_cov$wf_var_yield,
  var_disease  = cross_cov$wf_var_disease,
  stringsAsFactors = FALSE, row.names = NULL
)
trait_specs <- data.frame(
  trait    = c("yield", "disease"),
  mean_col = c("pred_yield", "pred_disease"),
  var_col  = c("var_yield", "var_disease"),
  stringsAsFactors = FALSE
)

# Acceptance region: progeny should be high yield AND low disease.
tau_lower <- c(stats::median(scores$pred_yield), -Inf)   # yield >= median
tau_upper <- c(Inf, stats::median(scores$pred_disease))  # disease <= median

# --- P(superior progeny): EXACT within-family covariance vs the population-correlation proxy ---
scored_exact <- ng_add_p_superior_progeny_multitrait(
  scores, trait_specs, tau_lower, tau_upper, k_progeny = 50L,
  cross_trait_cov = cross_cov, out_col = "p_exact"
)
scored_proxy <- ng_add_p_superior_progeny_multitrait(
  scores, trait_specs, tau_lower, tau_upper, k_progeny = 50L,
  G_hat = diag(2), out_col = "p_proxy"                    # zero cross-trait correlation proxy
)

cmp <- data.frame(
  scores[, c("parent1", "parent2")],
  p_exact = round(scored_exact$p_exact, 4),
  p_proxy = round(scored_proxy$p_proxy, 4)
)
cmp <- cmp[order(-cmp$p_exact), ]
cat("Top crosses by P(superior progeny) with the EXACT within-family covariance:\n")
print(utils::head(cmp, 8), row.names = FALSE)

stopifnot(all(scored_exact$p_exact >= 0 & scored_exact$p_exact <= 1))
stopifnot(max(abs(scored_exact$p_exact - scored_proxy$p_proxy)) > 1e-4)  # exact != proxy
message("Exact cross-trait covariance example completed.")
