# The GEBV branch of ng_choose_mean_source must be reachable.
#
# Regression test for the 2026-09-11 finding: ng_fit_ridge_effects deliberately
# stopped reporting `reliability` (a phenotype CV statistic is not accuracy^2
# against true breeding value, so calling it "reliability" would be wrong), but
# ng_choose_mean_source still gated the GEBV branch on
# `reliability_is_calibrated`. That flag is now always FALSE, so the branch
# became unreachable: every run silently fell back to the adjusted-phenotype
# mid-parent, no matter how well the markers predicted, and
# `min_effect_reliability` became a knob with no effect at any value.
#
# This test pins the behaviour that matters to a breeder: when the fit
# demonstrably predicts out of sample, the cross mean must come from GEBVs.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(verbose = FALSE)

set.seed(42)
n <- 150L
m <- 80L
geno <- matrix(sample(0:2, n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n))
rownames(geno) <- ids
colnames(geno) <- paste0("M", seq_len(m))

# Strong, genuinely genetic signal: 10 causal markers, high heritability, so the
# ridge fit predicts well out of sample and GEBVs are the better mean.
beta_true <- rep(0, m)
beta_true[seq_len(10L)] <- rnorm(10L, sd = 1.5)
g <- as.numeric(scale(geno %*% beta_true))
y <- setNames(g + rnorm(n, sd = 0.2), ids)   # cv_predictive_r2 ~ 0.87

fit <- ng_fit_ridge_effects(geno, y, ids = ids)

# Premise of the test: this fit really does predict out of sample.
stopifnot(is.finite(fit$cv_predictive_r2))
stopifnot(fit$cv_predictive_r2 > 0.35)

chosen <- ng_choose_mean_source(
  geno = geno,
  effects = fit,
  adjusted_pheno = y,
  ids = ids,
  min_reliability = 0.35
)

# THE ASSERTION: a well-predicting fit must yield GEBV-based means.
stopifnot(identical(chosen$source, "GEBV"))

# And the chosen values must actually be the GEBVs, not the phenotypes.
gebv <- ng_predict_gebv(geno[ids, , drop = FALSE], fit)
stopifnot(isTRUE(all.equal(unname(chosen$value), unname(gebv))))
stopifnot(!isTRUE(all.equal(unname(chosen$value), unname(y))))

message("mean_source_gebv_gate passed")
