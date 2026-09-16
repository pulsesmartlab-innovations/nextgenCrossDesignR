# The calibrated-reliability machinery is gone, not merely unreachable.
#
# ng_choose_mean_source() carried a branch selecting GEBV means on a CALIBRATED
# reliability -- r^2 against true breeding value -- gated on a flag no producer in this
# package ever set. It was documented as a reserved hook for a PEV-based reliability
# arriving later.
#
# It is not arriving, because it answers a question this package does not ask. PEV
# reliability exists to say how far to trust a GEBV for an UNPHENOTYPED selection
# candidate. Every parent here already has phenotypic records; marker effects are
# estimated to give those parents GEBVs and to supply beta-hat to the downstream
# quantities -- the a'Ra within-family variance and the cross-trait covariance. Neither
# is "predict this individual".
#
# So the hook was holding space for the wrong quantity, and a reader meeting the branch
# would reasonably conclude the package supports a mode it does not. One statistic
# decides whether a trait's marker effects are good enough, it is cv_predictive_r2, and
# the bar is the breeder's to set.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---- 1. the decision function no longer takes a reliability threshold ----
fm <- names(formals(ng_choose_mean_source))
stopifnot(!("min_reliability" %in% fm))
stopifnot("min_cv_predictive_r2" %in% fm)

# ---- 2. the fit reports no reliability fields ----------------------------
set.seed(9); n <- 30L; m <- 30L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
gv <- as.numeric(g %*% c(rnorm(6, 0, 1), rep(0, m - 6)))
y  <- stats::setNames(gv + rnorm(n, 0, 0.2 * stats::sd(gv)), ids)
fit <- ng_fit_ridge_effects(g, y, seed = 1L)
for (f in c("reliability", "in_sample_reliability", "reliability_is_calibrated")) {
  if (f %in% names(fit)) {
    stop("ng_fit_ridge_effects still reports '", f, "' -- a field that always carried ",
         "NA or FALSE and described a quantity this package does not compute",
         call. = FALSE)
  }
}
# the statistic that DOES decide is still there
stopifnot("cv_predictive_r2" %in% names(fit))

# ---- 3. the scored table carries no calibrated-reliability flag ----------
map <- ng_prepare_marker_map(
  data.frame(marker = colnames(g), chr = rep(1:2, each = m / 2),
             pos_cm = rep(seq(0, 60, length.out = m / 2), 2)), marker_ids = colnames(g))
sc <- ng_score_crosses(geno = g, effects = fit, marker_map = map, ids = ids,
                       pairs = data.frame(parent1 = ids[1:5], parent2 = ids[6:10],
                                          stringsAsFactors = FALSE),
                       adjusted_pheno = y, parent_type = "inbred", use_cpp = FALSE)
stopifnot(!("effect_reliability_is_calibrated" %in% names(sc)))
# and the decision it DOES record is intact
stopifnot("mean_source" %in% names(sc), "mean_source_criterion" %in% names(sc))

# ---- 4. the decision still works, both ways -----------------------------
# Removing the dead branch must not disturb the live one. cv_predictive_r2 alone
# decides, and it still decides differently at different thresholds.
lo <- ng_choose_mean_source(geno = g, effects = fit, ids = ids, adjusted_pheno = y,
                            min_cv_predictive_r2 = 0.01)$source
hi <- ng_choose_mean_source(geno = g, effects = fit, ids = ids, adjusted_pheno = y,
                            min_cv_predictive_r2 = 0.99)$source
stopifnot(identical(lo, "GEBV"), identical(hi, "adjusted_pheno"))

# ---- 5. and the criterion never claims a calibrated reliability ----------
# "calibrated_reliability" was one of the values mean_source_criterion could report.
# It could never actually occur; a reader finding it in the vocabulary would infer a
# mode that does not exist.
for (thr in c(0.01, 0.99)) {
  crit <- ng_choose_mean_source(geno = g, effects = fit, ids = ids, adjusted_pheno = y,
                                min_cv_predictive_r2 = thr)$mean_source_criterion
  stopifnot(!identical(as.character(crit), "calibrated_reliability"))
}

cat("calibrated_reliability_removed: PASS\n")
