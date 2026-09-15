# One mean-source policy, every entry point.
#
# ng_choose_mean_source() decides GEBV-vs-phenotype on min_cv_predictive_r2. Three
# callers reach it, and before this fix only ONE of them passed the threshold
# through:
#
#   ng_score_crosses()             R/03_metrics.R  -- passes it
#   ng_posterior_cross_predict()   R/30:536        -- did not even have the formal
#   ng_cheap_cross_screen()        R/09:26         -- has a formal, never passed it
#
# So a run configured with min_cv_predictive_r2 = 0.20 got 0.20 on the main scoring
# path and the hard-coded default of 0.35 everywhere else. Observed in production on
# the 2026-09-11 rerun: for six traits `<TRAIT>_mean` IS the GEBV mid-parent to 1e-9,
# while posterior_predictions[[t]]$mean_source read "adjusted_pheno" and
# effect_summary$mean_source read "GEBV" -- one result.json, two fields, opposite
# answers to the same question.
#
# That is worse than a wrong label. A checker written against the label reported a
# confident FAIL on a run that had succeeded, which is how a silent policy split
# turns into lost confidence in correct results.
#
# The invariant: for identical inputs and an identical threshold, every entry point
# must reach the same decision. Asserted on the decision, not on the plumbing.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(23)
n <- 40L; m <- 60L
geno <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
ids <- sprintf("P%02d", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- sprintf("M%02d", seq_len(m))
beta <- stats::setNames(c(rnorm(10, sd = 1.5), rep(0, m - 10)), colnames(geno))
g_true <- as.numeric(geno %*% beta)
# Noise tuned so cv_predictive_r2 lands BETWEEN the two thresholds the test uses
# (~0.31). A fixture above both, or below both, cannot tell a threaded threshold
# from an ignored one -- which is exactly how this defect stayed hidden.
y <- stats::setNames(g_true + rnorm(n, sd = 1.2 * stats::sd(g_true)), ids)
marker_map <- data.frame(marker = colnames(geno), chr = rep(1:3, each = m / 3),
                         pos_cm = rep(seq(0, 80, length.out = m / 3), 3),
                         stringsAsFactors = FALSE)
fit <- ng_fit_ridge_effects(geno, y, ids = ids)
pairs <- ng_make_pairs(ids)[1:15, , drop = FALSE]

# The fixture must sit BETWEEN the two thresholds, or the test cannot tell a
# threaded threshold from an ignored one.
r2 <- fit$cv_predictive_r2
stopifnot(is.finite(r2), r2 > 0.05, r2 < 0.35)

scored_at <- function(thresh) {
  ng_score_crosses(geno = geno, effects = fit, marker_map = marker_map, ids = ids,
                   pairs = pairs, adjusted_pheno = y, target = "DH",
                   parent_type = "inbred", use_cpp = FALSE,
                   min_cv_predictive_r2 = thresh)$mean_source[[1L]]
}
post_at <- function(thresh) {
  pe <- ng_fit_ridge_effects_posterior(geno, y, ids = ids, n_draws = 4L, seed = 5L)
  ng_posterior_cross_predict(geno = geno, posterior_effects = pe, marker_map = marker_map,
                             ids = ids, pairs = pairs, adjusted_pheno = y,
                             target = "DH", parent_type = "inbred", use_cpp = FALSE,
                             min_cv_predictive_r2 = thresh)$mean_source[[1L]]
}
cheap_at <- function(thresh) {
  ng_cheap_cross_screen(geno = geno, effects = fit, ids = ids, adjusted_pheno = y,
                        min_cv_predictive_r2 = thresh)$mean_source[[1L]]
}

# ---- below the bar: every entry point must refuse the markers ---------------
lo <- c(scored = scored_at(0.99), posterior = post_at(0.99), cheap = cheap_at(0.99))
stopifnot(length(unique(lo)) == 1L)
stopifnot(identical(unname(lo[["scored"]]), "adjusted_pheno"))

# ---- above the bar: every entry point must accept them ----------------------
hi <- c(scored = scored_at(0.01), posterior = post_at(0.01), cheap = cheap_at(0.01))
stopifnot(length(unique(hi)) == 1L)
stopifnot(identical(unname(hi[["scored"]]), "GEBV"))

# ---- and the threshold must actually be doing the work ---------------------
# If any entry point ignored it, that entry point would return the same answer at
# both thresholds -- which is precisely how the defect hid.
stopifnot(!identical(lo[["scored"]],    hi[["scored"]]))
stopifnot(!identical(lo[["posterior"]], hi[["posterior"]]))
stopifnot(!identical(lo[["cheap"]],     hi[["cheap"]]))

# ---- the label must match the values it describes --------------------------
# The failure that surfaced this was a mean_source saying "adjusted_pheno" on a
# table whose cross_mean WAS the GEBV mid-parent. Pin the agreement directly.
s_hi <- ng_score_crosses(geno = geno, effects = fit, marker_map = marker_map, ids = ids,
                         pairs = pairs, adjusted_pheno = y, target = "DH",
                         parent_type = "inbred", use_cpp = FALSE,
                         min_cv_predictive_r2 = 0.01)
stopifnot(identical(s_hi$mean_source[[1L]], "GEBV"))
stopifnot(isTRUE(all.equal(s_hi$cross_mean, s_hi$cross_mean_gebv)))

s_lo <- ng_score_crosses(geno = geno, effects = fit, marker_map = marker_map, ids = ids,
                         pairs = pairs, adjusted_pheno = y, target = "DH",
                         parent_type = "inbred", use_cpp = FALSE,
                         min_cv_predictive_r2 = 0.99)
stopifnot(identical(s_lo$mean_source[[1L]], "adjusted_pheno"))
stopifnot(!isTRUE(all.equal(s_lo$cross_mean, s_lo$cross_mean_gebv)))

cat("mean_source_policy_parity: PASS\n")
