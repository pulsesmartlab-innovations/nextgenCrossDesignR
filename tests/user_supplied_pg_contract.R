helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 0.27.0 -- contract for USER-SUPPLIED P and G arriving from a GUI/JSON bridge.
#
#   Change 1: the covariance -> value_z transform is LABEL-AWARE. A labelled
#             matrix is reordered by name; a mislabelled one is a hard error; an
#             unlabelled one is documented as positional.
#   Change 2: desired_gain (Pesek-Baker, b = G^{-1} d) needs ONLY G. P never
#             enters the coefficient solve, so without P the index is unchanged
#             and only the P-scaled reported response degrades to NA.
#   Change 3: ng_estimate_genetic_covariance() warns when Y looks like
#             already-shrunk predictions (BLUP/GEBV) rather than phenotypes.

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }

set.seed(20270101)

# ---- Candidate-cross trait table ------------------------------------------
n <- 60L
ids <- paste0("P", seq_len(20L))
pairs <- t(utils::combn(ids, 2L))[seq_len(n), , drop = FALSE]
yield   <- stats::rnorm(n, mean = 10, sd = 2.0)
disease <- stats::rnorm(n, mean = 30, sd = 5.0) - 0.3 * yield
lodging <- stats::rnorm(n, mean = 20, sd = 4.0) + 0.4 * yield
scores <- data.frame(parent1 = pairs[, 1L], parent2 = pairs[, 2L],
                     yield = yield, disease = disease, lodging = lodging,
                     stringsAsFactors = FALSE)
traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "lodging"),
  direction = c("maximize", "minimize", "minimize"),
  economic_weight = c(2, 1, 1),
  desired_change = c(4, 6, 3)
)
tn <- traits$trait

# Non-diagonal P and G so a permutation genuinely changes the answer.
G <- matrix(c( 4.0, -3.0,  1.5,
              -3.0, 25.0, -4.0,
               1.5, -4.0, 16.0), 3L, 3L, byrow = TRUE,
            dimnames = list(tn, tn))
P <- matrix(c( 9.0, -4.0,  2.0,
              -4.0, 60.0, -6.0,
               2.0, -6.0, 35.0), 3L, 3L, byrow = TRUE,
            dimnames = list(tn, tn))
stopifnot(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) > 0)

perm <- c(3L, 1L, 2L)   # lodging, yield, disease

fit_index <- function(Pm, Gm, method = "economic_index") {
  ng_add_multitrait_score(scores, traits, method = method,
                          phenotypic_covariance = Pm, genetic_covariance = Gm)
}

# ============================================================================
# Change 1a: a permuted-but-LABELLED matrix gives a bit-identical index
# ============================================================================
ref <- fit_index(P, G)
perm_labelled <- fit_index(P[perm, perm], G[perm, perm])
stopifnot(identical(ref$multi_trait_score, perm_labelled$multi_trait_score))
stopifnot(identical(attr(ref, "multi_trait")$economic_index_coefficients,
                    attr(perm_labelled, "multi_trait")$economic_index_coefficients))
ok("permuted labelled P/G reorder by name -> bit-identical index")

# The numeric difference this prevents: the SAME permuted matrix read
# positionally (which is exactly what an unlabelled JSON round trip delivers).
perm_positional <- fit_index(unname(P[perm, perm]), unname(G[perm, perm]))
b_ref <- attr(ref, "multi_trait")$economic_index_coefficients
b_bad <- attr(perm_positional, "multi_trait")$economic_index_coefficients
coef_gap <- max(abs(b_ref - b_bad))
rank_rho <- stats::cor(ref$multi_trait_score, perm_positional$multi_trait_score,
                       method = "spearman")
cat(sprintf("  permuted-read-positionally vs correct: max |db| = %.4f, spearman(rank) = %.4f\n",
            coef_gap, rank_rho))
stopifnot(coef_gap > 0.05)
ok("positional read of the same permuted matrix is measurably wrong")

# jsonlite really does drop dimnames on a matrix round trip (the GUI risk).
if (requireNamespace("jsonlite", quietly = TRUE)) {
  round_tripped <- jsonlite::fromJSON(jsonlite::toJSON(G))
  stopifnot(is.null(dimnames(round_tripped)))
  ok("jsonlite matrix round trip drops dimnames (unlabelled -> positional)")
}

# ============================================================================
# Change 1b: mislabelled / incomplete labels are hard errors
# ============================================================================
err <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))

G_typo <- G[perm, perm]
rownames(G_typo)[1L] <- colnames(G_typo)[1L] <- "lodgingg"
msg <- err(fit_index(P, G_typo))
stopifnot(!is.na(msg), grepl("do not match the trait set", msg, fixed = TRUE),
          grepl("lodgingg", msg, fixed = TRUE), grepl("lodging", msg, fixed = TRUE))
ok("misspelled trait label errors and names the bad label")

G_missing <- G[c(1L, 2L), c(1L, 2L)]
msg <- err(fit_index(P, G_missing))
stopifnot(!is.na(msg), grepl("Missing: lodging", msg, fixed = TRUE))
ok("incomplete (missing-trait) matrix errors and names the missing trait")

G_extra <- rbind(cbind(G, extra = c(0, 0, 0)), extra = c(0, 0, 0, 1))
msg <- err(fit_index(P, G_extra))
stopifnot(!is.na(msg), grepl("unexpected: extra", msg, fixed = TRUE))
ok("extra-trait matrix errors and names the unexpected label")

# The silent-positional hole this closes: labels on ONE dimension only. Before
# 0.27.0 this fell through to a positional read in whatever order it arrived.
G_rows_only <- G[perm, perm]
colnames(G_rows_only) <- NULL
msg <- err(fit_index(P, G_rows_only))
stopifnot(!is.na(msg), grepl("dimnames on only one dimension", msg, fixed = TRUE))
ok("one-sided dimnames error instead of a silent positional read")

# ============================================================================
# Change 1c: unlabelled input still works, positionally
# ============================================================================
unlabelled <- fit_index(unname(P), unname(G))
stopifnot(identical(ref$multi_trait_score, unlabelled$multi_trait_score))
ok("unlabelled P/G in trait order -> identical to the labelled run (positional)")

msg <- err(fit_index(unname(P), unname(G[c(1L, 2L), c(1L, 2L)])))
stopifnot(!is.na(msg), grepl("interpreted POSITIONALLY", msg, fixed = TRUE))
ok("unlabelled matrix of the wrong dimension errors with the positional contract")

# The low-level transform itself is label-aware, not only the entry point.
vs <- stats::setNames(c(1.5, 3.0, 2.5), tn)
z_ref <- ng_multitrait_cov_to_value_z(G, traits, vs)
z_perm <- ng_multitrait_cov_to_value_z(G[perm, perm], traits, vs)
stopifnot(identical(z_ref, z_perm))
ok("ng_multitrait_cov_to_value_z() reorders labelled input by name")

# ============================================================================
# Change 2: desired_gain needs only G; economic_index still needs both
# ============================================================================
dg_both <- fit_index(P, G, method = "desired_gain")
dg_G_only <- fit_index(NULL, G, method = "desired_gain")
m_both <- attr(dg_both, "multi_trait")
m_gonly <- attr(dg_G_only, "multi_trait")

# The Pesek-Baker solve is b = G^{-1} d: P is genuinely absent from it, so the
# coefficients and the emitted index must be bit-identical with and without P.
stopifnot(identical(m_both$desired_gain_coefficients, m_gonly$desired_gain_coefficients))
stopifnot(identical(dg_both$multi_trait_score, dg_G_only$multi_trait_score))
ok("desired_gain succeeds with G alone; coefficients/index bit-identical to the P+G run")

# P only ever fed the index SD sqrt(b' P b) used to standardise the REPORTED
# response, so that (and only that) degrades to NA.
stopifnot(all(is.finite(m_both$desired_gain_predicted_response)))
stopifnot(all(is.na(m_gonly$desired_gain_predicted_response)))
stopifnot(identical(m_gonly$desired_gain_unavailable, "predicted_response"))
stopifnot(grepl("requires phenotypic_covariance (P)", m_gonly$desired_gain_unavailable_reason,
                fixed = TRUE))
stopifnot(is.na(m_gonly$desired_gain_index_sd), is.finite(m_both$desired_gain_index_sd))
ok("without P the only casualty is predicted_response, reported as NA + a named reason")

msg <- err(fit_index(P, NULL, method = "desired_gain"))
stopifnot(!is.na(msg), grepl("requires genetic_covariance (G)", msg, fixed = TRUE),
          grepl("b = G^{-1} d", msg, fixed = TRUE),
          grepl("candidate-score covariance is not a substitute", msg, fixed = TRUE))
ok("desired_gain without G errors informatively (and keeps the candidate-score refusal)")

msg <- err(fit_index(NULL, G, method = "economic_index"))
stopifnot(!is.na(msg), grepl("requires both phenotypic_covariance", msg, fixed = TRUE))
msg <- err(fit_index(P, NULL, method = "economic_index"))
stopifnot(!is.na(msg), grepl("requires both phenotypic_covariance", msg, fixed = TRUE))
ok("economic_index (Smith-Hazel b = P^{-1} G a) still requires BOTH P and G")

# The summary carries the unavailability forward for a caller/GUI to display.
plan_scores <- scores
plan_scores$pair_kinship <- 0
K <- diag(length(ids)); dimnames(K) <- list(ids, ids)
plan <- ng_optimize_multitrait_mating_plan(
  scores = plan_scores, traits = traits, n_crosses = 2L, parent_kinship = K,
  multitrait_method = "desired_gain", optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L, genetic_covariance = G)
psum <- attr(plan, "summary")
stopifnot(identical(psum$multitrait_desired_gain_unavailable, "predicted_response"))
stopifnot(all(is.na(psum$multitrait_desired_gain_predicted_response)))
ok("mating-plan summary reports which quantities are unavailable without P")

# ============================================================================
# Change 3: already-shrunk-input (BLUP/GEBV) guard
# ============================================================================
set.seed(4242)
# n must clear the guard's min_n = 100; the REML heritability is deliberately
# not attempted below that because it is too erratic to act on.
n_ind <- 200L; n_mrk <- 500L
geno <- matrix(stats::rbinom(n_ind * n_mrk, 2L, 0.3), n_ind, n_mrk,
               dimnames = list(paste0("L", seq_len(n_ind)), paste0("M", seq_len(n_mrk))))
b1 <- stats::rnorm(n_mrk, sd = 0.30)
b2 <- 0.6 * b1 + stats::rnorm(n_mrk, sd = 0.24)
g1 <- as.numeric(scale(geno %*% b1))
g2 <- as.numeric(scale(geno %*% b2))
Y_pheno <- cbind(yield = g1 + stats::rnorm(n_ind, sd = 1.0),
                 protein = g2 + stats::rnorm(n_ind, sd = 1.0))
rownames(Y_pheno) <- rownames(geno)
# GEBV-like input: a pure linear function of the same markers the GRM is built
# from -- no residual left at all.
Y_gebv <- cbind(yield = g1, protein = g2)
rownames(Y_gebv) <- rownames(geno)

catch_shrunk <- function(expr) {
  msgs <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  any(grepl("ALREADY-SHRUNK", msgs, fixed = TRUE))
}

G_pheno <- NULL
fired_pheno <- catch_shrunk(
  G_pheno <- ng_estimate_genetic_covariance(geno, Y_pheno, method = "two_stage_ridge"))
G_gebv <- NULL
fired_gebv <- catch_shrunk(
  G_gebv <- ng_estimate_genetic_covariance(geno, Y_gebv, method = "two_stage_ridge"))
cat(sprintf("  guard REML h2: phenotypes = [%s], GEBVs = [%s]\n",
            paste(sprintf("%.4f", attr(G_pheno, "implied_heritability")), collapse = ", "),
            paste(sprintf("%.4f", attr(G_gebv, "implied_heritability")), collapse = ", ")))
cat(sprintf("  guard REML residual/observed variance: phenotypes = [%s], GEBVs = [%s]\n",
            paste(sprintf("%.3f", attr(G_pheno, "reml_residual_variance") /
                            apply(Y_pheno, 2, stats::var)), collapse = ", "),
            paste(sprintf("%.2e", attr(G_gebv, "reml_residual_variance") /
                            apply(Y_gebv, 2, stats::var)), collapse = ", ")))
stopifnot(isTRUE(fired_gebv))
stopifnot(isTRUE(attr(G_gebv, "shrunk_input_suspected")))
# `protein` here sits just under the tolerance (h2 ~ 0.986), which is exactly the
# point: this is a warning-level heuristic, so it names the traits it can see.
stopifnot("yield" %in% attr(G_gebv, "shrunk_input_traits"))
stopifnot(isFALSE(fired_pheno))
stopifnot(isFALSE(attr(G_pheno, "shrunk_input_suspected")))
ok("guard fires on GEBV-like Y and stays silent on phenotypes")

# The guard is engine-independent: same verdict through the sommer REML engine.
if (requireNamespace("sommer", quietly = TRUE)) {
  G_s <- NULL
  fired_sommer <- catch_shrunk(
    G_s <- ng_estimate_genetic_covariance(geno, Y_gebv, method = "sommer_remml"))
  cat(sprintf("  sommer_remml on GEBVs: guard h2 = [%s], engine residual variance = [%s]\n",
              paste(sprintf("%.4f", attr(G_s, "implied_heritability")), collapse = ", "),
              paste(sprintf("%.2e", attr(G_s, "residual_variance")), collapse = ", ")))
  stopifnot(isTRUE(fired_sommer))
  ok("guard fires through the sommer REML engine too")
}

# The guard is a warning, not an error: G-hat is still returned.
stopifnot(is.matrix(G_gebv), all(dim(G_gebv) == c(2L, 2L)), all(is.finite(G_gebv)))
ok("guard warns without erroring -- G-hat is still returned")

# Guard is skipped, with a stated reason, when n is too small to act on.
small_ids <- rownames(geno)[seq_len(60L)]
G_small <- NULL
fired_small <- catch_shrunk(
  G_small <- ng_estimate_genetic_covariance(geno[small_ids, ], Y_gebv[small_ids, ],
                                            method = "two_stage_ridge"))
stopifnot(isFALSE(fired_small))
stopifnot(is.na(attr(G_small, "implied_heritability_method")))
stopifnot(grepl("skipped: n = 60", attr(G_small, "implied_heritability_note"), fixed = TRUE))
ok("guard is skipped with a stated reason when n is too small")

# Provenance attributes a caller needs to report what this G-hat is.
required_attrs <- c("method", "requested_method", "formal_variance_component_estimate",
                    "genetic_correlation", "n_used", "n_markers", "traits",
                    "y_input_contract", "genetic_variance", "residual_variance",
                    "implied_heritability", "implied_heritability_method",
                    "implied_heritability_note", "reml_genetic_variance",
                    "reml_residual_variance", "shrunk_input_suspected",
                    "shrunk_input_traits")
stopifnot(all(required_attrs %in% names(attributes(G_pheno))))
stopifnot(identical(attr(G_pheno, "traits"), colnames(Y_pheno)))
stopifnot(identical(attr(G_pheno, "n_markers"), ncol(geno)))
stopifnot(grepl("BLUEs", attr(G_pheno, "y_input_contract"), fixed = TRUE))
ok("provenance attributes complete (engine, n, markers, traits, variances, Y contract)")

cat(sprintf("user_supplied_pg_contract.R: PASS (%d checks)\n", checks))
