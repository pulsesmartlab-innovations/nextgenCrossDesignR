# D7: p_beat_all_checks (the joint) now integrates the SAME shared posterior effect uncertainty
# (delta ~ N(0, PEV), diagonal across checked traits -- no cross-trait PEV covariance is
# estimated anywhere in this package) that <trait>_p_beat_check (D1, ng_p_superior_progeny_pev())
# integrates via Gauss-Hermite. The joint integrates by a fixed, seeded Monte Carlo instead
# (ng_p_superior_progeny_multitrait_pev(), R/33), because multivariate Gauss-Hermite is
# impractical beyond a couple of traits, and clips the result to a PROVABLE ceiling -- the
# row-wise minimum of the (independently recomputed) per-trait marginals -- so the invariant
# p_beat_all_checks <= min(marginals) holds EXACTLY, not merely "usually, within Monte Carlo
# noise" (see R/51_check_reference.R::ng_attach_joint_check_probability()'s docstring for the
# subset-and-monotonicity argument that makes the ceiling provable rather than a guess).
#
# BEFORE this fix, a real barley validation run (147 parents x 3189 markers, two correlated
# `decrease` traits, most crosses on the wrong side of the check) failed exactly this assertion:
# "P(beat all checks) never exceeds the B_glucan marginal". Section A below reproduces that shape
# synthetically. The RED evidence (this exact reproduction run against the pre-fix code, showing
# 120 of 200 rows violating the invariant by up to 0.28) is recorded in the fix report,
# .superpowers/sdd/2026-09-03-check-reference-lines-backend/joint-pev-fix-report.md -- re-deriving
# it here would require literally reverting the fix, which is not something a regression test
# should do. This file asserts the fix holds GREEN, going forward.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- Section 1: PEV = 0 (or unavailable) must be EXACT, not approximate --------------------------
# The cheapest and most important correctness check: with no posterior effect uncertainty (or an
# absent _pmv_used column), ng_p_superior_progeny_multitrait_pev() must return EXACTLY what
# ng_p_superior_progeny_multitrait() returns -- no Monte Carlo cost, no approximation error.
mu2 <- c(-1.2, 0.6); Sigma2 <- matrix(c(2, 0.7, 0.7, 1.5), 2)
tl2 <- c(0, -Inf); tu2 <- c(Inf, 1)
plain <- ng_p_superior_progeny_multitrait(mu2, Sigma2, tl2, tu2, k_progeny = 20L)
for (pev0 in list(c(0, 0), c(NA_real_, 0), c(-1, NA_real_))) {
  got <- ng_p_superior_progeny_multitrait_pev(mu2, Sigma2, tl2, tu2, k_progeny = 20L, pev = pev0)
  stopifnot(isTRUE(all.equal(got, plain, tolerance = 1e-14)))
}
cat("Section 1 ok: PEV = 0 / unavailable collapses EXACTLY (no MC draws performed)\n")

# The full attach-level path: yield_vpm == yield_pmv_used (PEV = 0) must reproduce the pre-D7
# joint value exactly, for a real cross_trait_cov table (not just the diagonal fallback).
set.seed(11)
n1 <- 6L; m1 <- 24L
ids1 <- sprintf("R%02d", seq_len(n1)); mk1 <- sprintf("N%03d", seq_len(m1))
geno1 <- matrix(2L * rbinom(n1 * m1, 1, 0.5), n1, m1, dimnames = list(ids1, mk1))
storage.mode(geno1) <- "double"
map1 <- data.frame(marker = mk1, chr = rep(1:2, each = m1 / 2),
                   pos_cm = rep(seq(0, 90, length.out = m1 / 2), 2), stringsAsFactors = FALSE)
bA1 <- rnorm(m1); bB1 <- 0.5 * bA1 + rnorm(m1)
betas1 <- cbind(traitA = bA1, traitB = bB1); rownames(betas1) <- mk1
pairs1 <- as.data.frame(t(utils::combn(ids1, 2)), stringsAsFactors = FALSE)
names(pairs1) <- c("parent1", "parent2")
ctc1 <- ng_cross_trait_within_family_cov(geno1, betas1, map1, pairs = pairs1, target = "DH")
gA1 <- as.numeric(geno1 %*% bA1); names(gA1) <- ids1
gB1 <- as.numeric(geno1 %*% bB1); names(gB1) <- ids1
p1_1 <- pairs1$parent1; p2_1 <- pairs1$parent2
scores1 <- data.frame(
  parent1 = p1_1, parent2 = p2_1,
  traitA_mean = 0.5 * (gA1[p1_1] + gA1[p2_1]), traitA_vpm = ctc1$wf_var_traitA,
  traitA_pmv_used = ctc1$wf_var_traitA,                       # PEV = 0
  traitB_mean = 0.5 * (gB1[p1_1] + gB1[p2_1]), traitB_vpm = ctc1$wf_var_traitB,
  traitB_pmv_used = ctc1$wf_var_traitB,                        # PEV = 0
  stringsAsFactors = FALSE, row.names = NULL)
spec1 <- ng_trait_check_spec(trait = c("traitA", "traitB"), check = c("CHK_A", "CHK_B"),
                             trait_direction = c(traitA = "increase", traitB = "decrease"))
cv1 <- list(traitA = c(CHK_A = unname(stats::median(scores1$traitA_mean))),
           traitB = c(CHK_B = unname(stats::median(scores1$traitB_mean))))
joint1_pev0 <- ng_attach_joint_check_probability(scores1, spec1, cv1, k_progeny = 12L,
                                                 cross_trait_cov = ctc1)
# The pre-D7 code path: pass pev_mat = NULL directly to the lower-level function to recompute
# "what the old code would have given" without any PEV machinery at all.
trait_specs1 <- data.frame(trait = spec1$trait, mean_col = paste0(spec1$trait, "_mean"),
                           var_col = paste0(spec1$trait, "_vpm"), stringsAsFactors = FALSE)
b1 <- ng_check_tau_bounds(spec1, cv1)
pre_d7 <- ng_add_p_superior_progeny_multitrait(
  scores1, trait_specs1, tau_lower = as.numeric(b1$tau_lower[spec1$trait]),
  tau_upper = as.numeric(b1$tau_upper[spec1$trait]), k_progeny = 12L,
  cross_trait_cov = ctc1, pev_mat = NULL, out_col = "p_beat_all_checks")
stopifnot(isTRUE(all.equal(joint1_pev0$p_beat_all_checks, pre_d7$p_beat_all_checks,
                           tolerance = 1e-12)))
cat("Section 1b ok: PEV = 0 on the REAL cross_trait_cov path reproduces the pre-D7 value exactly\n")

# --- Section 2: Monte Carlo estimate matches a LITERAL, independent simulation of the full model -
# Simulate delta ~ N(0, diag(pev)), then k progeny ~ MVN(mu + delta, Sigma_c), then check whether
# ANY of the k progeny clears every threshold simultaneously. This is a direct simulation of the
# MODEL, not a restatement of the closed form under test.
mc_truth_joint <- function(mu, Sigma, pev, tau_lower, tau_upper, k, nsim = 20000L, seed = 1L) {
  set.seed(seed)
  t_n <- length(mu)
  sd_d <- sqrt(pmax(pev, 0))
  L <- t(chol(Sigma))
  hit <- logical(nsim)
  for (i in seq_len(nsim)) {
    delta <- if (any(sd_d > 0)) stats::rnorm(t_n, 0, sd_d) else rep(0, t_n)
    Z <- matrix(stats::rnorm(k * t_n), nrow = t_n, ncol = k)
    X <- (mu + delta) + L %*% Z                     # t_n x k progeny values
    inside <- apply(X, 2L, function(x) all(x >= tau_lower & x <= tau_upper))
    hit[[i]] <- any(inside)
  }
  p_hat <- mean(hit)
  list(p_hat = p_hat, se = sqrt(p_hat * (1 - p_hat) / nsim))
}

cases <- list(
  list(mu = c(-1, -1.5), Sigma = matrix(c(2, 0.9, 0.9, 2.5), 2), pev = c(0.6, 0.8),
       tau_lower = c(0, 0), tau_upper = c(Inf, Inf), k = 15L),
  list(mu = c(0.5, 0.3), Sigma = matrix(c(1.5, -0.6, -0.6, 1.8), 2), pev = c(0.4, 0.3),
       tau_lower = c(-Inf, -Inf), tau_upper = c(0, 0), k = 10L),
  list(mu = c(-2, 1), Sigma = matrix(c(3, 1.2, 1.2, 2), 2), pev = c(1.0, 0.5),
       tau_lower = c(0, -Inf), tau_upper = c(Inf, 0.5), k = 20L)
)
cat(sprintf("\n%-6s %10s %10s %10s %10s %8s\n", "case", "impl", "MC_truth", "MC_se", "|diff|", "within"))
for (i in seq_along(cases)) {
  cs <- cases[[i]]
  impl <- ng_p_superior_progeny_multitrait_pev(cs$mu, cs$Sigma, cs$tau_lower, cs$tau_upper,
                                               cs$k, cs$pev)
  truth <- mc_truth_joint(cs$mu, cs$Sigma, cs$pev, cs$tau_lower, cs$tau_upper, cs$k,
                          nsim = 20000L, seed = 100L + i)
  tol <- max(5 * truth$se, 1e-3)
  d <- abs(impl - truth$p_hat)
  cat(sprintf("%-6d %10.5f %10.5f %10.5f %10.5f %8s\n", i, impl, truth$p_hat, truth$se, d,
              ifelse(d < tol, "yes", "NO")))
  stopifnot(d < tol)
}
cat("Section 2 ok: Monte Carlo joint matches a literal simulation of the full model\n")

# --- Section 3: the invariant on the SHAPE that broke it (real data reproduction) ----------------
# Two correlated DECREASE traits, means mostly on the WRONG side of the check (the check is at 0;
# decrease wants mu <= 0), genuine PEV (~30-40% of VPM, plausible for a moderately-sized training
# set), non-trivial cross-trait correlation. This is the exact configuration that broke the
# invariant pre-fix (see the fix report for the RED transcript: 120 of 200 rows violating it,
# some by as much as 0.28).
mk_stress <- function(seed, n, rho, mu_sd, pev_frac) {
  set.seed(seed)
  vpm1 <- stats::runif(n, 2, 6); vpm2 <- stats::runif(n, 2, 6)
  cov12 <- rho * sqrt(vpm1 * vpm2)
  pev1 <- pev_frac * vpm1; pev2 <- pev_frac * vpm2
  mu1 <- stats::rnorm(n, mean = 2, sd = mu_sd); mu2 <- stats::rnorm(n, mean = 1.5, sd = mu_sd)
  scores <- data.frame(
    parent1 = paste0("P", seq_len(n)), parent2 = paste0("Q", seq_len(n)),
    traitA_mean = mu1, traitA_vpm = vpm1, traitA_pmv_used = vpm1 + pev1,
    traitB_mean = mu2, traitB_vpm = vpm2, traitB_pmv_used = vpm2 + pev2,
    stringsAsFactors = FALSE)
  ctc <- data.frame(parent1 = scores$parent1, parent2 = scores$parent2,
                    wf_var_traitA = vpm1, wf_var_traitB = vpm2,
                    wf_cov_traitA_traitB = cov12, stringsAsFactors = FALSE)
  list(scores = scores, ctc = ctc)
}
spec3 <- ng_trait_check_spec(trait = c("traitA", "traitB"), check = c("CHK_A", "CHK_B"),
                             trait_direction = c(traitA = "decrease", traitB = "decrease"))
cv3 <- list(traitA = c(CHK_A = 0), traitB = c(CHK_B = 0))
k3 <- 30L
d3 <- mk_stress(seed = 1L, n = 200L, rho = 0.7, mu_sd = 2, pev_frac = 0.35)
single3 <- ng_attach_check_reference(d3$scores, spec3, NULL, cv3, k_progeny = k3)
joint3  <- ng_attach_joint_check_probability(d3$scores, spec3, cv3, k_progeny = k3,
                                             cross_trait_cov = d3$ctc)
diff3 <- joint3$p_beat_all_checks - pmin(single3$traitA_p_beat_check, single3$traitB_p_beat_check)
cat(sprintf("Section 3: n=%d crosses, max(joint - min(marginals)) = %.10f, violations = %d\n",
            nrow(d3$scores), max(diff3, na.rm = TRUE), sum(diff3 > 1e-9, na.rm = TRUE)))
stopifnot(all(diff3 <= 1e-9))
# Sanity: this configuration DOES exercise genuine PEV-driven separation (the marginal, once PEV
# is integrated, differs materially from the naive VPM-only closed form) -- otherwise this
# wouldn't be testing anything the pre-fix code got right by accident. traitA is a DECREASE trait
# (sgn = -1: margin = tau - mu), so the naive VPM-only closed form on the same sgn convention is
# ng_p_superior_progeny(-mu, sqrt(vpm), -tau, k).
naive3 <- ng_p_superior_progeny(-d3$scores$traitA_mean, sqrt(d3$scores$traitA_vpm), 0, k3)
stopifnot(max(abs(naive3 - single3$traitA_p_beat_check)) > 0.05)
cat("Section 3 ok: invariant holds on the real-data-shaped stress case (correlated decrease",
    "traits, below-check means, genuine PEV)\n")

# --- Section 4: correlation sweep -- one fixture is not evidence -------------------------------
cat(sprintf("\n%-6s %8s %8s\n", "rho", "max_gap", "n_viol"))
for (rho in c(-0.8, 0, 0.8)) {
  d <- mk_stress(seed = 7L, n = 150L, rho = rho, mu_sd = 1.6, pev_frac = 0.3)
  single <- ng_attach_check_reference(d$scores, spec3, NULL, cv3, k_progeny = 18L)
  joint  <- ng_attach_joint_check_probability(d$scores, spec3, cv3, k_progeny = 18L,
                                              cross_trait_cov = d$ctc)
  diff <- joint$p_beat_all_checks - pmin(single$traitA_p_beat_check, single$traitB_p_beat_check)
  cat(sprintf("%-6.1f %8.2e %8d\n", rho, max(diff, na.rm = TRUE), sum(diff > 1e-9, na.rm = TRUE)))
  stopifnot(all(diff <= 1e-9, na.rm = TRUE))
}
cat("Section 4 ok: invariant holds at rho = -0.8, 0, +0.8\n")

# --- Section 5: wall-clock cost on a realistic candidate count ---------------------------------
n_bench <- 1000L
d5 <- mk_stress(seed = 99L, n = n_bench, rho = 0.6, mu_sd = 1.5, pev_frac = 0.3)
t0 <- Sys.time()
invisible(ng_attach_joint_check_probability(d5$scores, spec3, cv3, k_progeny = 20L,
                                            cross_trait_cov = d5$ctc))
dt <- as.numeric(Sys.time() - t0, units = "secs")
cat(sprintf("Section 5: %d crosses, n_mc_draws = %d -> %.2fs (%.4fs/cross, ~%.0fs per 10k crosses)\n",
            n_bench, NG_JOINT_PEV_MC_DRAWS, dt, dt / n_bench, dt / n_bench * 10000))

cat("\ncheck_reference_joint_pev.R: PASS\n")
