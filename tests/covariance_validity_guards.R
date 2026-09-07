helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 0.29.0 -- BLOCKING validity guards for USER-SUPPLIED covariance matrices.
#
#   Guard 1  symmetry tolerance is RELATIVE to the matrix scale (was absolute,
#            while the PSD test six lines below was already relative).
#   Guard 2  P - G must be positive semidefinite, because P = G + R makes
#            R = P - G a covariance matrix. Multivariate 0 <= h2 <= 1.
#   Guard 3  the matrix the index INVERTS must be positive DEFINITE:
#            economic_index inverts P, desired_gain inverts G.
#   Guard 4  implied correlations must lie in [-1, 1], reported as correlations.
#
# Every guard is a hard error. Every guard is also shown NOT to fire on a valid
# matrix, and a valid P/G pair is shown to be numerically untouched.

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }
err <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))
has <- function(msg, ...) {
  stopifnot(!is.na(msg))
  for (frag in c(...)) {
    if (!grepl(frag, msg, fixed = TRUE)) {
      stop("message did not mention '", frag, "'. Message was:\n", msg, call. = FALSE)
    }
  }
  invisible(TRUE)
}

set.seed(20290101)

# ---- Candidate-cross trait table (same shape as user_supplied_pg_contract) ---
n <- 60L
ids <- paste0("P", seq_len(20L))
pairs <- t(utils::combn(ids, 2L))[seq_len(n), , drop = FALSE]
yield   <- stats::rnorm(n, mean = 10, sd = 2.0)
protein <- stats::rnorm(n, mean = 13, sd = 1.5) - 0.3 * yield
lodging <- stats::rnorm(n, mean = 20, sd = 4.0) + 0.4 * yield
scores <- data.frame(parent1 = pairs[, 1L], parent2 = pairs[, 2L],
                     yield = yield, protein = protein, lodging = lodging,
                     stringsAsFactors = FALSE)
traits <- ng_multitrait_spec(
  trait = c("yield", "protein", "lodging"),
  direction = c("maximize", "maximize", "minimize"),
  economic_weight = c(2, 1, 1),
  desired_change = c(4, 2, 3)
)
tn <- traits$trait

# A VALID pair: G and P from the same source and scale, P - G comfortably PD.
G_ok <- matrix(c( 4.0, 1.2, -1.5,
                  1.2, 2.0, -0.8,
                 -1.5, -0.8, 9.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))
P_ok <- matrix(c(10.0, 2.0, -2.5,
                  2.0, 5.0, -1.2,
                 -2.5, -1.2, 20.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))
stopifnot(min(eigen(G_ok, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(P_ok, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(P_ok - G_ok, symmetric = TRUE, only.values = TRUE)$values) > 0)

fit <- function(Pm, Gm, method = "economic_index", ...) {
  ng_add_multitrait_score(scores, traits, method = method,
                          phenotypic_covariance = Pm, genetic_covariance = Gm, ...)
}

# ============================================================================
# The valid pair goes straight through both index methods.
# ============================================================================
ref_ei <- fit(P_ok, G_ok, method = "economic_index")
ref_dg <- fit(P_ok, G_ok, method = "desired_gain")
stopifnot(all(is.finite(ref_ei$multi_trait_score)),
          all(is.finite(ref_dg$multi_trait_score)))
stopifnot(identical(attr(ref_ei, "multi_trait")$economic_index_cov_source, "smith_hazel"))
stopifnot(identical(attr(ref_dg, "multi_trait")$desired_gain_cov_source, "pesek_baker"))
ok("a P/G pair satisfying every guard is accepted by both index methods")

# ============================================================================
# GUARD 1 -- symmetry, judged RELATIVE to the matrix scale
# ============================================================================

# 1a. THE DEFECT THE RELATIVE FORM FIXES. Yield in kg/ha instead of t/ha scales
# every covariance by 1e6. A matrix rounded to 9 significant digits -- what a
# spreadsheet or a JSON bridge delivers -- then disagrees with its transpose by
# ~1e-3 ABSOLUTE while agreeing to ~1e-9 RELATIVE. The old absolute 1e-8 test
# refused it; the relative test accepts it, and (because the matrix is
# symmetrised anyway) gives the same answer as the exactly-symmetric version.
big <- 1e6
G_kg <- G_ok * big
P_kg <- P_ok * big
G_round <- G_kg
G_round[1L, 2L] <- signif(G_round[1L, 2L], 9L) + 1e-3   # ~1e-9 relative
asym_abs <- max(abs(G_round - t(G_round)))
asym_rel <- asym_abs / max(abs(G_round))
cat(sprintf("  kg/ha-scale rounding: absolute asymmetry = %.3e, relative = %.3e\n",
            asym_abs, asym_rel))
stopifnot(asym_abs > 1e-8)          # the old absolute rule would have refused it
stopifnot(asym_rel < 1e-8)          # the relative rule accepts it
rounded <- fit(P_kg, G_round)
exact  <- fit(P_kg, G_kg)
stopifnot(max(abs(rounded$multi_trait_score - exact$multi_trait_score)) < 1e-9)
ok("relative symmetry tolerance accepts a rounded kg/ha matrix the absolute rule refused")

# 1b. Genuine asymmetry is still a hard error, and the message names the cells.
G_asym <- G_ok
G_asym[1L, 2L] <- G_asym[1L, 2L] + 0.5
msg <- err(fit(P_ok, G_asym))
has(msg, "genetic_covariance must be symmetric", "cov(yield, protein)", "cov(protein, yield)",
    "a discrepancy of")
cat("  guard 1 message:", substr(msg, 1L, 150L), "...\n")
ok("genuine asymmetry is still a hard error naming both cells and the discrepancy")

# 1c. The relative rule scales: the SAME relative asymmetry is refused at kg/ha too.
G_asym_kg <- G_kg
G_asym_kg[1L, 2L] <- G_asym_kg[1L, 2L] * 1.001     # 1e-3 relative
msg <- err(fit(P_kg, G_asym_kg))
has(msg, "must be symmetric")
ok("a relatively-large asymmetry is refused at kg/ha scale too (the rule is not a blanket relaxation)")

# ============================================================================
# GUARD 2 -- P - G must be positive semidefinite (P = G + R)
# ============================================================================

# 2a. Per-trait h2 > 1, the elementwise special case, gets the clearest message.
G_h2 <- G_ok
G_h2[2L, 2L] <- 6.0                                # P_ok[2,2] = 5.0
msg <- err(fit(P_ok, G_h2))
has(msg, "heritability above 1", "trait 'protein'", "genetic_covariance variance 6",
    "phenotypic_covariance variance 5", "implying h2 = 1.2")
cat("  guard 2a message:", substr(msg, 1L, 190L), "...\n")
ok("per-trait genetic > phenotypic variance is refused, naming trait, both numbers and h2")

# 2b. REALISTIC SCALE MISMATCH: G supplied in t/ha, P supplied in kg/ha (or vice
# versa) -- the mistake a GUI with two upload boxes makes easy. Here G is 1000x
# too large relative to P.
G_wrong_scale <- G_ok * 1000
msg <- err(fit(P_ok, G_wrong_scale))
has(msg, "heritability above 1", "SAME scale and in the SAME units")
cat("  guard 2b (unit mismatch) message:", substr(msg, 1L, 170L), "...\n")
ok("a G/P unit mismatch (G 1000x too large) is refused with a units diagnosis")

# 2c. THE GENERAL CASE: every per-trait h2 is legal, but a trait CONTRAST has
# h2 > 1. This is the case a per-trait check cannot see. G is built with a much
# stronger genetic correlation than P has phenotypic correlation, so the
# contrast yield - protein has more genetic than phenotypic variance.
G_contrast <- matrix(c(4.0, 3.9, 0.0,
                       3.9, 4.0, 0.0,
                       0.0, 0.0, 9.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))
P_contrast <- matrix(c(5.0, 1.0, 0.0,
                       1.0, 5.0, 0.0,
                       0.0, 0.0, 20.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))
stopifnot(all(diag(G_contrast) < diag(P_contrast)))            # every per-trait h2 < 1
stopifnot(min(eigen(G_contrast, symmetric = TRUE, only.values = TRUE)$values) > 0)
stopifnot(min(eigen(P_contrast, symmetric = TRUE, only.values = TRUE)$values) > 0)
stopifnot(min(eigen(P_contrast - G_contrast, symmetric = TRUE,
                    only.values = TRUE)$values) < 0)           # but P - G is not PSD
msg <- err(fit(P_contrast, G_contrast))
has(msg, "residual covariance R = phenotypic_covariance - genetic_covariance",
    "not positive semidefinite", "Smallest eigenvalue of P - G",
    "trait contrast that realises it", "implied h2 =")
cat("  guard 2c message:", substr(msg, 1L, 260L), "...\n")
ok("a contrast-level h2 > 1 that no per-trait check can see is refused, with the contrast named")

# 2d. The guard is order-blind: the same invalid pair permuted (and labelled) is
# refused identically, so trait order cannot launder an invalid pair.
perm <- c(3L, 1L, 2L)
msg_perm <- err(fit(P_contrast[perm, perm], G_contrast[perm, perm]))
has(msg_perm, "not positive semidefinite")
ok("the P - G guard runs after alignment, so a permuted invalid pair is refused too")

# 2e. It does NOT fire on the valid pair, nor when only one matrix is supplied.
stopifnot(is.na(err(fit(P_ok, G_ok))))
stopifnot(is.na(err(fit(NULL, G_ok, method = "desired_gain"))))
ok("the P - G guard is silent on a valid pair and on a G-only desired_gain run")

# ============================================================================
# GUARD 3 -- positive DEFINITE for the matrix the index inverts
# ============================================================================

# A singular G: lodging is exactly yield + protein in genetic terms.
A <- matrix(c(1, 0, 1,
              0, 1, 1), nrow = 2L, byrow = TRUE)
G_sing <- crossprod(A) * 2
dimnames(G_sing) <- list(tn, tn)
stopifnot(abs(min(eigen(G_sing, symmetric = TRUE, only.values = TRUE)$values)) < 1e-12)
P_for_sing <- G_sing + diag(c(4, 4, 8))
dimnames(P_for_sing) <- list(tn, tn)

# 3a. desired_gain INVERTS G -> singular G is refused.
msg <- err(fit(P_for_sing, G_sing, method = "desired_gain"))
has(msg, "desired_gain (Pesek-Baker, b = G^{-1} d) inverts genetic_covariance",
    "must be positive DEFINITE", "numerically SINGULAR", "condition number",
    "Drop one of the redundant traits")
cat("  guard 3a message:", substr(msg, 1L, 230L), "...\n")
ok("desired_gain refuses a singular G (the matrix it inverts) instead of ridging it")

# 3b. economic_index INVERTS P, and only P: the SAME singular G is accepted there,
# because G enters as the forward projection G a and is never inverted.
scored <- fit(P_for_sing, G_sing, method = "economic_index")
stopifnot(all(is.finite(scored$multi_trait_score)))
ok("economic_index accepts the same singular G -- it projects through G, it does not invert it")

# 3c. ... and refuses a singular P.
P_sing <- G_sing
G_small <- G_sing * 0.25
msg <- err(fit(P_sing, G_small, method = "economic_index"))
has(msg, "economic_index (Smith-Hazel, b = P^{-1} G a) inverts phenotypic_covariance",
    "must be positive DEFINITE")
ok("economic_index refuses a singular P (the matrix it inverts)")

# 3d. WHAT THE UN-GUARDED SOLVE ACTUALLY RETURNED on this matrix, i.e. what the
# guard is preventing. ng_multitrait_solve_index() is called here directly,
# bypassing the guard, on exactly the singular G refused above.
#
# The answer is not merely imprecise; it is entirely the direction G says carries
# NO genetic variance. G_sing's null eigenvector is v0 ~ (1, 1, -1)/sqrt(3), and
# LAPACK reports that eigenvalue as ~5e-15 rather than 0 -- just above the
# routine's own keep-tolerance -- so the null direction is retained with a
# coefficient of order 1/5e-15, swamps the two real directions, and survives the
# sum(|b|) normalisation as the whole index. The returned b is therefore fixed by
# floating-point noise in the eigendecomposition, down to its sign, and the index
# it defines has b' G b = 0: zero predicted genetic response.
b_sing <- ng_multitrait_solve_index(
  stats::setNames(c(0.5, 0.3, 0.2), tn),
  list(target_matrix = G_sing, projection = NULL, response_G = G_sing,
       response_P = P_for_sing), ridge = 1e-6)$coefficients
es_sing <- eigen(G_sing, symmetric = TRUE)
v0 <- es_sing$vectors[, which.min(es_sing$values)]
align <- abs(sum(b_sing * v0)) / sqrt(sum(b_sing^2))
var_frac <- as.numeric(crossprod(b_sing, G_sing %*% b_sing)) /
  (sum(b_sing^2) * max(es_sing$values))
cat(sprintf("  un-guarded b on the singular G = (%s)\n",
            paste(sprintf("%+.4f", b_sing), collapse = ", ")))
cat(sprintf("  |cos angle to G's NULL eigenvector| = %.6f; b' G b / (b'b * lambda_max) = %.3e\n",
            align, var_frac))
stopifnot(align > 0.999)          # the index IS the null direction
stopifnot(abs(var_frac) < 1e-12)  # it carries no genetic variance at all
ok("un-guarded, the singular G returned an index lying entirely in G's null space (b' G b = 0)")

# 3e. MILD ill-conditioning is still accepted -- the ridge legitimately handles it.
G_mild <- G_sing + diag(1e-4, 3L)     # the singular G, lifted just off the boundary
dimnames(G_mild) <- list(tn, tn)
kappa_mild <- {
  ev <- eigen(G_mild, symmetric = TRUE, only.values = TRUE)$values
  max(ev) / min(ev)
}
P_mild <- G_mild + diag(c(6, 3, 8))
dimnames(P_mild) <- list(tn, tn)
stopifnot(kappa_mild > 1e4, kappa_mild < 1e8)
scored_mild <- fit(P_mild, G_mild, method = "desired_gain")
stopifnot(all(is.finite(scored_mild$multi_trait_score)))
cat(sprintf("  mildly ill-conditioned G accepted: condition number = %.4g (limit 1e+08)\n",
            kappa_mild))
ok("a mildly ill-conditioned G (kappa ~ 1e4-1e5) is accepted, as the ridge is designed to allow")

# ============================================================================
# GUARD 4 -- implied correlations inside [-1, 1]
# ============================================================================

# 4a. Impossible genetic correlation of 1.31 between yield and protein.
G_r <- G_ok
G_r[1L, 2L] <- G_r[2L, 1L] <- 1.31 * sqrt(G_ok[1L, 1L] * G_ok[2L, 2L])
msg <- err(fit(P_ok, G_r))
has(msg, "implies a correlation outside [-1, 1]", "yield vs protein",
    "implied correlation 1.31")
cat("  guard 4 message:", substr(msg, 1L, 200L), "...\n")
ok("an impossible implied correlation is refused, reported AS a correlation")

# 4b. Guard 4 wins over the eigenvalue test: the same matrix is also non-PSD, but
# the clearer objection is the one the breeder sees.
stopifnot(min(eigen(G_r, symmetric = TRUE, only.values = TRUE)$values) < 0)
stopifnot(!grepl("positive semidefinite", msg, fixed = TRUE))
ok("guard 4 runs BEFORE the PSD test, so the actionable message wins")

# 4c. A correlation of exactly +/- 1 (a legitimately singular, PSD matrix) passes
# guard 4 -- it is on the boundary, not outside it.
G_unit <- G_ok
G_unit[1L, 3L] <- G_unit[3L, 1L] <- -sqrt(G_ok[1L, 1L] * G_ok[3L, 3L])
stopifnot(is.na(err(ng_multitrait_check_implied_correlations(G_unit, "genetic_covariance"))))
ok("a correlation of exactly -1 is on the boundary and is not refused by guard 4")

# ============================================================================
# BIT-IDENTITY: a valid pair is numerically untouched by all four guards
# ============================================================================
# The reference values are produced by a PRISTINE 0.28.0 tree in
# tests/covariance_guards_bit_identity.R; here we assert internal consistency of
# the guarded path so a later refactor cannot drift it.
stopifnot(identical(fit(P_ok, G_ok)$multi_trait_score, ref_ei$multi_trait_score))
stopifnot(identical(fit(P_ok, G_ok, method = "desired_gain")$multi_trait_score,
                    ref_dg$multi_trait_score))
ok("repeated guarded runs on the valid pair are bit-identical")

cat(sprintf("covariance_validity_guards.R: PASS (%d checks)\n", checks))
