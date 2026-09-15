# One question, one column name, on every scoring path.
#
# 0.32.0 put per-row estimator provenance on the diploid table as
# `variance_estimator`, because ng_apply_het_parent_correction() substitutes a
# different kernel for a subset of rows and two rows in one table could otherwise
# come from two estimators indistinguishably.
#
# The polyploid paths had already answered the same question -- but under a
# different name, `variance_model` (R/46:112, R/49:210). Two names for one fact is
# the defect class this release exists to remove: a consumer that learned
# `variance_estimator` on the diploid path reads NULL on the polyploid one and
# concludes the provenance is absent, when it is merely spelled differently.
#
# So `variance_estimator` must be present and populated wherever a within-family
# variance is scored. `variance_model` stays as-is: it is an existing output
# column and dropping it would break readers for no gain.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(11)
n <- 10L; m <- 24L
ids <- sprintf("P%02d", seq_len(n))
dos <- matrix(sample(0:4, n * m, TRUE), nrow = n,
              dimnames = list(ids, sprintf("M%02d", seq_len(m))))
eff <- stats::setNames(rnorm(m, sd = 0.4), colnames(dos))

# ---- autopolyploid additive path (R/46) -----------------------------------
a <- ng_polyploid_score_crosses(dosage = dos, effects = eff, ploidy = 4L)
stopifnot("variance_estimator" %in% names(a))
stopifnot(all(!is.na(a$variance_estimator)), all(nzchar(a$variance_estimator)))
# and it must agree with the column that already carried the answer
stopifnot(identical(a$variance_estimator, a$variance_model))

# ---- additive + dominance path (R/49) -------------------------------------
fit <- list(dosage = dos, ploidy = 4L,
            beta_additive = eff, beta_dominance = NULL)
d <- tryCatch(ng_polyploid_score_crosses_dominance(fit), error = function(e) e)
if (!inherits(d, "error")) {
  stopifnot("variance_estimator" %in% names(d))
  stopifnot(identical(d$variance_estimator, d$variance_model))
}

# ---- the vocabulary is shared, not merely the column name -----------------
# A reader must be able to compare the two paths' answers. Both use the
# phase-resolution vocabulary; neither invents a private label.
known <- c("phased_exact", "uniform_phase_prior_expectation",
           "dh_recomb_aRa", "dh_recomb_full_posterior", "phased_het_general")
stopifnot(all(a$variance_estimator %in% known))

cat("variance_estimator_column_is_universal: PASS\n")
