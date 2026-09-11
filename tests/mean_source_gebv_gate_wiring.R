# `min_cv_predictive_r2` must be reachable from the scoring entry point.
#
# Companion to mean_source_gebv_gate.R. That test pins the decision inside
# ng_choose_mean_source(); this one pins that a caller can actually REACH the
# decision, which is the half that made the original defect unfixable from
# outside: `min_effect_reliability` was exposed all the way up to
# ng_run_cross_prediction() but gated on a flag the estimator never sets, so no
# value of it -- 0 included -- could restore GEBV means.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(verbose = FALSE)

set.seed(42)
n <- 150L
m <- 80L
# Inbred parents are homozygous by construction: 0/2 only, no hets.
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n))
rownames(geno) <- ids
colnames(geno) <- paste0("M", seq_len(m))

beta_true <- rep(0, m)
beta_true[seq_len(10L)] <- rnorm(10L, sd = 1.5)
g <- as.numeric(scale(geno %*% beta_true))
y <- setNames(g + rnorm(n, sd = 0.2), ids)      # cv_predictive_r2 ~ 0.87

fit <- ng_fit_ridge_effects(geno, y, ids = ids)
stopifnot(fit$cv_predictive_r2 > 0.5)

# Small marker map so ng_score_crosses can compute within-family variance.
marker_map <- data.frame(marker = colnames(geno),
                         chr = rep(1:4, length.out = m),
                         pos_cm = rep(seq(0, 150, length.out = m / 4), times = 4),
                         stringsAsFactors = FALSE)
pairs <- data.frame(parent1 = ids[1:6], parent2 = ids[7:12], stringsAsFactors = FALSE)

score <- function(thresh) {
  ng_score_crosses(geno = geno, effects = fit, marker_map = marker_map, ids = ids,
                   pairs = pairs, adjusted_pheno = y, target = "RIL",
                   parent_type = "inbred", use_cpp = FALSE,
                   min_cv_predictive_r2 = thresh)
}

# THE ASSERTION: the threshold must be reachable from the scoring entry point,
# and must change the mean that comes out of it.
hi <- score(0.99)   # bar the fit cannot clear -> phenotype mid-parent
lo <- score(0.35)   # bar the fit clears        -> GEBV mid-parent

stopifnot("mean_source" %in% names(hi))
stopifnot(identical(unique(hi$mean_source), "adjusted_pheno"))
stopifnot(identical(unique(lo$mean_source), "GEBV"))
stopifnot(!isTRUE(all.equal(hi$cross_mean, lo$cross_mean)))

message("mean_source_gebv_gate_wiring passed")
