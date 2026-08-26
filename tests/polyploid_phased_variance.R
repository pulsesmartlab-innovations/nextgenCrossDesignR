ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# EXACT phased autopolyploid within-family variance (R/49).
#
#   Var(gamete) = [ P * sum_i (b h_i)' R (b h_i) - (b d)' R (b d) ] / (4 (P-1))
#
# Three properties are load-bearing and each is checked below:
#   1. under R = I it collapses EXACTLY to the dosage-only moment table (phase may change the
#      answer only through LINKAGE, never on its own);
#   2. it matches directly simulated meiosis (random bivalent pairing + crossovers);
#   3. with real linkage it REORDERS crosses relative to the dosage path -- which is why it
#      matters: the autopolyploid variance feeds usefulness, which feeds allocation.

set.seed(21)
P <- 4L; n <- 8L; m <- 18L
ids <- sprintf("C%d", seq_len(n))
rn <- as.vector(t(outer(ids, seq_len(P), function(a, b) paste0(a, "_Hap", b))))
H <- matrix(rbinom(n * P * m, 1, 0.5), n * P, m,
            dimnames = list(rn, sprintf("M%02d", seq_len(m))))
D <- t(vapply(ids, function(p) colSums(H[paste0(p, "_Hap", seq_len(P)), , drop = FALSE]),
              numeric(m)))
b <- rnorm(m, 0, 1)
pairs <- ng_make_pairs(ids, include_self = FALSE)
mm_link <- data.frame(marker = colnames(H), chr = 1L, pos_cm = cumsum(c(0, rep(1.5, m - 1))))
mm_far  <- data.frame(marker = colnames(H), chr = 1L, pos_cm = cumsum(c(0, rep(1e6, m - 1))))

# --- 1. unlinked limit == dosage-only moment table, EXACTLY -------------------------------
v_far <- ng_poly_phased_within_family_var(H, pairs, b, mm_far, ploidy = P)
mt <- ng_polyploid_progeny_moment_table(P, double_reduction = 0)
v_dose <- vapply(seq_len(nrow(pairs)), function(i) {
  ij <- cbind(D[pairs$parent1[i], ] + 1L, D[pairs$parent2[i], ] + 1L)
  sum(b^2 * mt$varX[ij])
}, numeric(1))
stopifnot(max(abs(v_far - v_dose)) < 1e-8)
stopifnot(all(v_far > 0))
cat("phased variance collapses to the dosage-only result when unlinked\n")

# single-locus sanity: must be the hypergeometric d(P-d)/(4(P-1)) per parent
for (Pk in c(2L, 4L, 6L)) for (d in 1:(Pk - 1)) {
  h1 <- matrix(c(rep(1, d), rep(0, Pk - d)), Pk, 1,
               dimnames = list(paste0("A_Hap", seq_len(Pk)), "M1"))
  h2 <- matrix(0, Pk, 1, dimnames = list(paste0("B_Hap", seq_len(Pk)), "M1"))
  hh <- rbind(h1, h2)
  vv <- ng_poly_phased_within_family_var(
    hh, data.frame(parent1 = "A", parent2 = "B", stringsAsFactors = FALSE),
    stats::setNames(1, "M1"), data.frame(marker = "M1", chr = 1L, pos_cm = 0), ploidy = Pk)
  stopifnot(abs(vv - d * (Pk - d) / (4 * (Pk - 1))) < 1e-12)
}
cat("single-locus phased variance equals the hypergeometric result for ploidy 2/4/6\n")

# --- 2. agreement with simulated meiosis --------------------------------------------------
v_link <- ng_poly_phased_within_family_var(H, pairs, b, mm_link, ploidy = P)
r_adj <- 0.5 * (1 - exp(-2 * diff(mm_link$pos_cm) / 100))          # Haldane, adjacent
sim_gamete <- function(Hp) {
  perm <- sample(P); g <- numeric(m)
  for (bv in seq_len(P / 2)) {
    a <- perm[2 * bv - 1]; q <- perm[2 * bv]
    cur <- if (stats::runif(1) < 0.5) a else q
    for (k in seq_len(m)) {
      if (k > 1 && stats::runif(1) < r_adj[k - 1]) cur <- if (cur == a) q else a
      g[k] <- g[k] + Hp[cur, k]
    }
  }
  g
}
for (i in c(1L, 9L)) {
  p1 <- pairs$parent1[i]; p2 <- pairs$parent2[i]
  H1 <- H[paste0(p1, "_Hap", seq_len(P)), , drop = FALSE]
  H2 <- H[paste0(p2, "_Hap", seq_len(P)), , drop = FALSE]
  G <- replicate(15000, sum(b * (sim_gamete(H1) + sim_gamete(H2))))
  # Monte-Carlo SE of a variance at n = 15000 is ~1.2%; allow 6%.
  stopifnot(abs(v_link[i] - stats::var(G)) / stats::var(G) < 0.06)
}
cat("phased variance matches simulated bivalent meiosis with crossovers\n")

# --- 3. linkage changes the DECISION, not just the number ----------------------------------
s_dose <- ng_polyploid_score_crosses(D, b, ploidy = P)
s_link <- ng_polyploid_score_crosses(D, b, ploidy = P,
                                     phased_haplotypes = H, marker_map = mm_link)
s_far  <- ng_polyploid_score_crosses(D, b, ploidy = P,
                                     phased_haplotypes = H, marker_map = mm_far)
stopifnot(identical(unique(s_dose$variance_model), "uniform_phase_prior_expectation"))
stopifnot(identical(unique(s_link$variance_model), "phased_exact"))
stopifnot(max(abs(s_far$poly_var - s_dose$poly_var)) < 1e-8)      # consistency through the scorer
stopifnot(max(abs(s_link$poly_var - s_dose$poly_var)) > 1e-6)     # linkage is not a no-op
# the ranking really moves -- this is why the phased path is worth having
stopifnot(stats::cor(s_dose$poly_usefulness, s_link$poly_usefulness, method = "spearman") < 0.999)
stopifnot(all(is.finite(s_link$poly_var)), all(s_link$poly_var >= 0))
cat("phased scoring reorders crosses relative to the dosage path\n")

# --- 4. guards: silence is the enemy here --------------------------------------------------
err <- function(expr) inherits(try(expr, silent = TRUE), "try-error")
stopifnot(err(ng_polyploid_score_crosses(D, b, ploidy = P, phased_haplotypes = H)))   # no map
stopifnot(err(ng_polyploid_score_crosses(D, b, ploidy = P, phased_haplotypes = H,
                                         marker_map = mm_link, double_reduction = 0.1)))
stopifnot(err(ng_poly_phased_within_family_var(H[-1, , drop = FALSE], pairs, b, mm_link, ploidy = P)))
stopifnot(err(ng_poly_phased_within_family_var(H, pairs, b[-1], mm_link, ploidy = P)))
# an unnamed beta of the right length is positional and must NOT silently become zero --
# this exact case returned an all-zero variance during development
v_unnamed <- ng_poly_phased_within_family_var(H, pairs, unname(b), mm_link, ploidy = P)
stopifnot(max(abs(v_unnamed - v_link)) < 1e-12, all(v_unnamed > 0))
# a beta named for markers the haplotypes do not carry is an error, not a silent zero
bad_b <- stats::setNames(b, paste0("Z", seq_len(m)))
stopifnot(err(ng_poly_phased_within_family_var(H, pairs, bad_b, mm_link, ploidy = P)))
stopifnot(err(ng_poly_phased_within_family_var(H, pairs, b, mm_link, ploidy = 3L)))   # odd ploidy
cat("phased variance guards test passed\n")
