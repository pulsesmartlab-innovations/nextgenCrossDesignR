# Dominance-aware polyploid cross scoring (R/49): the analytic cross mean + within-family variance
# (from the progeny-moment table) must match a brute-force gamete simulation, for any ploidy.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# sample n gametes from a parent of dosage d, with double-reduction coefficient dr (matches
# ng_poly_gamete_pmf's model: 2 duplicated copies of a random allele + rest hypergeometric).
sim_gamete <- function(d, ploidy, n, dr = 0) {
  g <- ploidy %/% 2L
  out <- integer(n)
  is_dr <- if (dr > 0 && g >= 2L) stats::runif(n) < dr else logical(n)
  nn <- sum(!is_dr)
  if (nn) out[!is_dr] <- stats::rhyper(nn, m = d, n = ploidy - d, k = g)
  for (t in which(is_dr)) {
    alt <- stats::rbinom(1, 1, d / ploidy)             # duplicated allele alt?
    rest <- if (alt == 1L) stats::rhyper(1, m = d - 1L, n = ploidy - d, k = g - 2L)
            else stats::rhyper(1, m = d, n = ploidy - 1L - d, k = g - 2L)
    out[t] <- 2L * alt + rest
  }
  out
}
# simulate n progeny of a cross (gamete_i + gamete_j), with optional double reduction
sim_progeny_dosage <- function(di, dj, ploidy, n, dr = 0) {
  m <- length(di); X <- matrix(0L, n, m)
  for (col in seq_len(m)) {
    X[, col] <- sim_gamete(di[col], ploidy, n, dr) + sim_gamete(dj[col], ploidy, n, dr)
  }
  X
}

set.seed(49)
for (ploidy in c(2L, 4L)) {
  m <- 25L
  ids <- sprintf("P%02d", 1:6)
  dosage <- matrix(rbinom(length(ids) * m, ploidy, 0.5), length(ids), m,
                   dimnames = list(ids, sprintf("m%02d", 1:m)))
  p <- colMeans(dosage) / ploidy
  hbar <- colMeans(dosage * (ploidy - dosage))
  ba <- rnorm(m, 0, 1); bd <- rnorm(m, 0, 1.5)
  fit <- structure(list(beta_add = ba, beta_dom = bd, intercept = 3, allele_freq = p, hbar = hbar,
                        markers = colnames(dosage), ploidy = ploidy, model = "additive_dominance"),
                   class = "ng_polyploid_effects")

  sc <- ng_polyploid_score_crosses_dominance(fit, dosage, selection_prop = 0.10)
  stopifnot(all(c("cross_mean", "heterosis", "add_var", "dom_var", "cross_var") %in% names(sc)))
  stopifnot(isTRUE(attr(sc, "has_dominance")))

  # validate the first few crosses against simulation
  for (r in 1:3) {
    a <- sc$parent1[r]; b <- sc$parent2[r]
    X <- sim_progeny_dosage(dosage[a, ], dosage[b, ], ploidy, n = 40000L)
    H <- X * (ploidy - X)
    gv <- fit$intercept + as.numeric(sweep(X, 2, ploidy * p) %*% ba) +
                          as.numeric(sweep(H, 2, hbar) %*% bd)
    stopifnot(abs(mean(gv) - sc$cross_mean[r]) < 0.15 * (abs(sc$cross_mean[r]) + 1))
    stopifnot(abs(stats::var(gv) - sc$cross_var[r]) < 0.10 * (sc$cross_var[r] + 1))
  }

  # heterosis is nonzero with dominance; additive-only fit gives zero dominance terms
  stopifnot(any(abs(sc$heterosis) > 1e-8), all(sc$dom_var >= 0))
  fit_a <- fit; fit_a$beta_dom <- NULL; fit_a$model <- "additive"
  sc_a <- ng_polyploid_score_crosses_dominance(fit_a, dosage)
  stopifnot(all(sc_a$dom_var == 0), all(abs(sc_a$heterosis) < 1e-12),
            max(abs(sc_a$cross_mean - sc_a$mid_parent_bv)) < 1e-12)
  cat(sprintf("ploidy %d: analytic cross mean/var match simulation; heterosis captured\n", ploidy))
}

# --- DOUBLE REDUCTION (tetraploid): analytic DR scores match a DR-aware simulation, and DR
# reduces expected heterozygosity (more homozygosity) vs no DR ---
set.seed(444)
m <- 25L; ids <- sprintf("P%02d", 1:6)
dosage <- matrix(rbinom(length(ids) * m, 4L, 0.5), length(ids), m, dimnames = list(ids, sprintf("m%02d", 1:m)))
p <- colMeans(dosage) / 4; hbar <- colMeans(dosage * (4 - dosage))
ba <- rnorm(m); bd <- rnorm(m, 0, 1.5)
fit4 <- structure(list(beta_add = ba, beta_dom = bd, intercept = 0, allele_freq = p, hbar = hbar,
                       markers = colnames(dosage), ploidy = 4L, model = "additive_dominance"),
                  class = "ng_polyploid_effects")
dr <- 0.12
sc_dr <- ng_polyploid_score_crosses_dominance(fit4, dosage, double_reduction = dr)
a <- sc_dr$parent1[1]; b <- sc_dr$parent2[1]
Xd <- sim_progeny_dosage(dosage[a, ], dosage[b, ], 4L, 40000L, dr = dr)
Hd <- Xd * (4 - Xd)
gvd <- as.numeric(sweep(Xd, 2, 4 * p) %*% ba) + as.numeric(sweep(Hd, 2, hbar) %*% bd)
stopifnot(abs(mean(gvd) - sc_dr$cross_mean[1]) < 0.15 * (abs(sc_dr$cross_mean[1]) + 1))
stopifnot(abs(stats::var(gvd) - sc_dr$cross_var[1]) < 0.12 * (sc_dr$cross_var[1] + 1))
# DR lowers total expected progeny heterozygosity across crosses
mt0 <- ng_polyploid_progeny_moment_table(4L, double_reduction = 0)
mtd <- ng_polyploid_progeny_moment_table(4L, double_reduction = dr)
stopifnot(sum(mtd$EH) < sum(mt0$EH))
cat(sprintf("double reduction (dr=%.2f): analytic matches DR simulation; heterozygosity reduced\n", dr))

# --- C++ kernel (when compiled) must match the R path to machine precision ---
if (exists("ng_poly_dominance_scores_cpp", mode = "function", inherits = TRUE)) {
  sc_cpp <- ng_polyploid_score_crosses_dominance(fit4, dosage, double_reduction = dr, use_cpp = TRUE)
  sc_r   <- ng_polyploid_score_crosses_dominance(fit4, dosage, double_reduction = dr, use_cpp = FALSE)
  stopifnot(max(abs(sc_cpp$cross_mean - sc_r$cross_mean)) < 1e-8,
            max(abs(sc_cpp$cross_var - sc_r$cross_var)) < 1e-8,
            max(abs(sc_cpp$heterosis - sc_r$heterosis)) < 1e-8)
  cat("C++ dominance-scoring kernel matches R path\n")
}

cat("polyploid dominance cross scoring test passed\n")

# --- analytic cross moments must reproduce SIMULATED progeny (0.19.x) -----------------------
# The strongest available check on the whole chain: gamete pmf -> progeny moment table ->
# orthogonal dominance basis -> cross_mean / cross_var. Simulate real meiosis (including double
# reduction), score the resulting progeny through the fitted model, and compare.
set.seed(41)
Ps <- 4L; ns <- 24L; ms <- 40L; drs <- 0.15
idss <- sprintf("S%02d", seq_len(ns))
fs <- runif(ms, 0.1, 0.9)
Ms <- matrix(rbinom(ns * ms, Ps, rep(fs, each = ns)), ns, ms,
             dimnames = list(idss, sprintf("V%03d", seq_len(ms))))
ys <- as.numeric(Ms %*% rnorm(ms, 0, .25) + (Ms * (Ps - Ms)) %*% rnorm(ms, 0, .12)) + rnorm(ns, 0, .5)
fits <- ng_polyploid_fit_effects(Ms, ys, ploidy = Ps, model = "additive_dominance", seed = 3L)
scs <- ng_polyploid_score_crosses_dominance(fits, Ms, double_reduction = drs, use_cpp = FALSE)

sim_gam <- function(dvec, P, dr) vapply(dvec, function(d) {
  g <- P %/% 2L; al <- c(rep(1L, d), rep(0L, P - d))
  if (g >= 2L && stats::runif(1) < dr) {
    j <- sample.int(P, 1L); rest <- al[-j]
    2L * al[j] + if (g > 2L) sum(sample(rest, g - 2L)) else 0L
  } else sum(sample(al, g))
}, integer(1))

nprog <- 3000
for (i in c(1L, 7L)) {
  a <- scs$parent1[i]; b <- scs$parent2[i]
  prog <- t(vapply(seq_len(nprog), function(z)
    sim_gam(Ms[a, ], Ps, drs) + sim_gam(Ms[b, ], Ps, drs), integer(ms)))
  colnames(prog) <- colnames(Ms); rownames(prog) <- sprintf("p%04d", seq_len(nprog))
  gv <- ng_polyploid_predict_value(fits, prog, type = "genotypic")
  # mean is exact in expectation; sd carries ~1/sqrt(2n) Monte Carlo error (~1.3% at n=3000)
  stopifnot(abs(scs$cross_mean[i] - mean(gv)) / max(abs(scs$cross_mean[i]), 1e-8) < 0.02)
  stopifnot(abs(sqrt(scs$cross_var[i]) - stats::sd(gv)) / sqrt(scs$cross_var[i]) < 0.12)
}

# The decomposition must be exact and non-negative: cross_var = add + dom + 2cov(X,D)
stopifnot(all(scs$add_var >= 0), all(scs$dom_var >= 0), all(scs$cross_var >= 0))
stopifnot(all(abs(scs$cross_mean - (scs$mid_parent_bv + scs$heterosis)) < 1e-10))
# ...and the C++ path must agree on the orthogonal basis (b_orth threaded through the kernel)
scc <- ng_polyploid_score_crosses_dominance(fits, Ms, double_reduction = drs, use_cpp = TRUE)
nm <- vapply(scc, is.numeric, logical(1))
stopifnot(max(abs(as.matrix(scc[, nm]) - as.matrix(scs[, nm]))) < 1e-9)
cat("polyploid analytic-vs-simulated progeny test passed\n")
