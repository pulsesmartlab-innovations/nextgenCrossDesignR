# Every posterior draw must use the cross-mean basis the RUN resolved.
#
# ng_posterior_cross_predict() re-scores each draw through ng_score_crosses(), which
# re-runs ng_choose_mean_source() per draw. The base scoring received the run's
# min_cv_predictive_r2 (R/30:549) but the PER-DRAW scoring did not, so the draws fell
# back to the 0.35 default and could resolve to GEBV while the run itself had resolved
# to the phenotypic mid-parent.
#
# The consequence is not a label mismatch, it is a wrong interval. A phenotypic
# mid-parent does not depend on the marker effects, so it cannot vary across draws and
# its posterior SD is exactly zero -- which is why ng_cross_confidence() reports
# `relative_precision_unavailable` for it rather than inventing a spread. With the
# draws on the GEBV basis instead, the mean varied, and a `mean` run was handed a
# posterior confidence interval describing a quantity that was NOT the number in its
# output. That is the same defect as posterior_predictions$mean_source, one layer down.
#
# This is only observable when the two bases genuinely disagree: good markers (so GEBV
# is reachable) plus a threshold the run cannot meet (so the run falls back). A weak
# fixture hides it, which is why it survived.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(31L)
n <- 60L; m <- 60L; ids <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
             dimnames = list(ids, sprintf("M%03d", seq_len(m))))
gv <- as.numeric(gm %*% c(rnorm(10L, 0, 1), rep(0, m - 10L)))
y  <- gv + rnorm(n, 0, 0.2 * stats::sd(gv))
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE)
mm <- data.frame(SNP = colnames(gm), chr = rep(1:2, length.out = m),
                 bp = rep(seq(0, 100, length.out = 30), 2)[seq_len(m)] * 1e6)
dir1 <- data.frame(trait = "yield", column = "yield", direction = "increase")

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = mm, trait_direction = dir1,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, n_crosses = 8L, max_crosses_per_parent = 3L,
  write_outputs = FALSE, write_figures = FALSE, seed = 5L,
  run_posterior_prediction = TRUE, n_iter = 80L, burn_in = 20L,
  posterior_method = "closed_form", ...)

# ---- the markers are genuinely good, so GEBV is reachable ------------------
r_gebv <- run(trait_value_metric = "mean")
stopifnot(identical(r_gebv$effect_summary$mean_source[[1L]], "GEBV"))
# a GEBV mid-parent IS a function of beta, so it really does vary across draws
stopifnot(identical(r_gebv$selected_crosses$confidence_method[[1L]], "posterior_ci"))
stopifnot(any(is.finite(r_gebv$selected_crosses$cross_confidence)))

# ---- now force the run onto the phenotypic basis ---------------------------
# min_cv_predictive_r2 = 1.1 is unreachable, so the documented fallback tier applies.
# `mean` is exempt from refusal, so this is the ordinary warn-and-fall-back path.
r_ph <- suppressWarnings(run(trait_value_metric = "mean", min_cv_predictive_r2 = 1.1))
stopifnot(identical(r_ph$effect_summary$mean_source[[1L]], "adjusted_pheno"))

# THE INVARIANT: the delivered value is phenotypic, so it cannot have a posterior
# spread, so no interval may be reported for it.
stopifnot(identical(r_ph$selected_crosses$confidence_method[[1L]],
                    "relative_precision_unavailable"))
stopifnot(all(is.na(r_ph$selected_crosses$cross_confidence)))
stopifnot(all(is.na(r_ph$selected_crosses$relative_precision)))
stopifnot(isFALSE(r_ph$priority_risk_diagnostics$posterior_used))

# ...and the spread itself must be degenerate, not merely reported as unavailable.
psd <- r_ph$candidate_crosses$yield_post_sd
if (!is.null(psd)) {
  fin <- psd[is.finite(psd)]
  stopifnot(!length(fin) || max(fin) < 1e-8)
}

# ---- the two runs must not disagree about which basis they used -----------
# One field answers this question. The posterior must not contradict effect_summary.
pp <- r_ph$posterior_predictions[["yield"]]
if (!is.null(pp) && !is.null(pp$mean_source)) {
  stopifnot(all(as.character(pp$mean_source) == "adjusted_pheno"))
}

cat("posterior_draws_follow_the_run_basis: PASS\n")
