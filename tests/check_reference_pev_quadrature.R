# D1: ng_p_superior_progeny_pev() (R/30_posterior_prediction.R) replaces the invalid
# "raise the WHOLE variance (VPM + shared posterior marker-effect uncertainty) to the k-th power"
# closed form with the correct model:
#   mu_true = mu_hat + delta,  delta ~ N(0, PEV)  shared across the k progeny of the family
#   P = E_delta[ 1 - Phi((tau - (mu_hat + delta)) / sqrt(VPM))^k ]
# via Gauss-Hermite quadrature. This file validates the quadrature against a LITERAL Monte Carlo
# simulation of that model (simulate delta once per replicate, then k progeny given delta, then
# check whether any progeny clears tau) -- not a restatement of the same closed form.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

mc_p_beat <- function(mu, vpm, pev, tau, k_progeny, n_mc = 200000L, seed = 1L) {
  set.seed(seed)
  delta <- if (pev > 0) rnorm(n_mc, 0, sqrt(pev)) else rep(0, n_mc)
  mu_true <- mu + delta
  # k progeny per replicate; TRUE iff at least one of the k clears tau.
  hit <- vapply(mu_true, function(m) {
    if (vpm <= 0) return(m >= tau)
    any(rnorm(k_progeny, m, sqrt(vpm)) >= tau)
  }, logical(1L))
  p_hat <- mean(hit)
  se <- sqrt(p_hat * (1 - p_hat) / n_mc)
  list(p_hat = p_hat, se = se)
}

# --- (mu - tau, k, PEV) grid, vpm = 1, tau = 0 throughout ------------------------------------
points <- data.frame(
  mu_minus_tau = c(-2,   -2,   -1,    0,    1,   -2,   -0.5),
  k            = c(100,  10,   50,   50,   50,   1,    20),
  pev          = c(1,    1,    0.25, 0.5,  2,    1,    0)
)

cat(sprintf("%-14s %6s %6s %10s %10s %10s %8s\n",
            "mu-tau", "k", "PEV", "GH", "MC", "MC_se", "within"))
results <- vector("list", nrow(points))
for (i in seq_len(nrow(points))) {
  mu <- points$mu_minus_tau[[i]]; k <- points$k[[i]]; pev <- points$pev[[i]]
  vpm <- 1; tau <- 0
  gh <- ng_p_superior_progeny_pev(mu, vpm, pev, tau, k)
  mc <- mc_p_beat(mu, vpm, pev, tau, k, n_mc = 200000L, seed = 100L + i)
  # 5 MC standard errors + a small floor for near-0/near-1 saturation.
  tol <- max(5 * mc$se, 1e-3)
  ok <- abs(gh - mc$p_hat) < tol
  cat(sprintf("%-14.4g %6d %6.3g %10.6f %10.6f %10.6f %8s\n",
              mu, k, pev, gh, mc$p_hat, mc$se, ifelse(ok, "yes", "NO")))
  results[[i]] <- ok
  stopifnot(ok)
}

# --- The headline case from the QG review: mu - tau = -2, k = 100, PEV = 1, vpm = 1 ----------
# Old (buggy) closed form on PMV = VPM + PEV = 2 gives ~0.9997; truth is ~0.68.
gh_headline <- ng_p_superior_progeny_pev(-2, 1, 1, 0, 100)
old_buggy <- 1 - stats::pnorm((0 - (-2)) / sqrt(1 + 1))^100
stopifnot(abs(gh_headline - 0.6769) < 0.01)          # matches the QG review's stated truth
stopifnot(old_buggy > 0.999)                          # confirms the bug's magnitude
stopifnot(gh_headline < old_buggy - 0.3)               # GH is materially different (and correct)
cat(sprintf("\nheadline: GH = %.4f  old-buggy-closed-form = %.4f  (QG review truth ~= 0.6769)\n",
            gh_headline, old_buggy))

# --- PEV = 0 (or unavailable) must be EXACT, not approximate, vs the plain closed form --------
for (pev0 in c(0, NA_real_, -1)) {
  # ng_p_superior_progeny_pev treats non-finite / <=0 pev identically: PEV collapses to 0.
  got <- ng_p_superior_progeny_pev(mu = -1.5, vpm = 2, pev = pev0, tau = 0.3, k_progeny = 30)
  want <- ng_p_superior_progeny(-1.5, sqrt(2), 0.3, 30)
  stopifnot(isTRUE(all.equal(got, want, tolerance = 1e-12)))
}
cat("PEV = 0 / unavailable collapses EXACTLY to the plain closed form on sqrt(VPM)\n")

# --- Decrease-trait mirror: sgn negation of both mu and tau must still answer "P(at least one of
# k progeny <= tau)". Cross-validated against an INDEPENDENT direct quadrature of that exact
# question (not a call back into the package's sgn convention), derived from
# P(X <= tau) = Phi((tau - X)/sigma), so P(>=1 of k <= tau) = 1 - Phi((X - tau)/sigma)^k, and
# against a literal Monte Carlo of the same model. -------------------------------------------
direct_low_side <- function(mu, vpm, pev, tau, k, n_nodes = 32L) {
  sigma <- sqrt(vpm)
  if (!is.finite(pev) || pev <= 0) {
    return(1 - stats::pnorm((mu - tau) / sigma)^k)
  }
  rule <- ng_gauss_hermite_rule(n_nodes)
  sigma_d <- sqrt(pev)
  delta <- sqrt(2) * sigma_d * rule$nodes
  f <- 1 - stats::pnorm(((mu + delta) - tau) / sigma)^k
  (1 / sqrt(pi)) * sum(rule$weights * f)
}
mu_r <- 3; vpm_r <- 1.4; pev_r <- 0.6; tau_r <- 3.8; k_r <- 12
p_mirrored <- ng_p_superior_progeny_pev(-mu_r, vpm_r, pev_r, -tau_r, k_r)     # P(<= tau) via sgn=-1
p_direct <- direct_low_side(mu_r, vpm_r, pev_r, tau_r, k_r)                  # same question, independent derivation
stopifnot(isTRUE(all.equal(p_mirrored, p_direct, tolerance = 1e-10)))
mc_dec <- mc_p_beat(-mu_r, vpm_r, pev_r, -tau_r, k_r, n_mc = 200000L, seed = 999L)
stopifnot(abs(p_mirrored - mc_dec$p_hat) < max(5 * mc_dec$se, 1e-3))
cat(sprintf("mirror check: sgn-negated = %.6f, independent direct quadrature = %.6f, MC = %.4f +/- %.4f\n",
            p_mirrored, p_direct, mc_dec$p_hat, mc_dec$se))

# --- Vectorized over rows, and bounds respected -----------------------------------------------
mu_v <- c(-2, 0, 1, 0.3); vpm_v <- c(1, 0.5, 2, 0); pev_v <- c(1, 0, 0.3, 0.2); tau_v <- rep(0, 4)
out_v <- ng_p_superior_progeny_pev(mu_v, vpm_v, pev_v, tau_v, k_progeny = 25)
stopifnot(length(out_v) == 4L, all(out_v >= 0 & out_v <= 1))
# vpm = 0 (row 4): every progeny is deterministic at mu + delta; with pev > 0 this is no longer a
# pure indicator (delta still varies), so it must lie strictly between 0 and 1, not the sigma<=0
# hard indicator ng_p_superior_progeny() alone would give.
stopifnot(out_v[[4L]] > 0, out_v[[4L]] < 1)

cat("check_reference_pev_quadrature.R: PASS\n")
