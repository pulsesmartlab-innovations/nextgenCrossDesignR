# Every scored row must say which estimator produced its variance, and whether
# that variance is a real PMV.
#
# Two defects, both of the family this release is about: a number that is right
# for its own method, presented with nothing to say which method that was.
#
# DEFECT A -- ng_apply_het_parent_correction() (R/03_metrics.R:263-268) OVERWRITES
# vpm / pmv / pmv_full_posterior for exactly the crosses whose RIL parents carry
# residual heterozygosity, using a different kernel (ng_gms_additive_var_general)
# that is dense and to which window_cm does not apply. It writes no column, no
# flag and no count. Two rows in one table can therefore come from two different
# estimators, under two different truncation regimes, with nothing whatsoever to
# distinguish them. A breeder comparing those rows is comparing incomparable
# numbers and cannot know it.
#
# DEFECT B -- when effects$beta_var is NULL, R/03_metrics.R zero-fills it, which
# makes PMV algebraically identical to VPM: the marker-effect uncertainty term is
# all zeros. The run still reports the column as `pmv`. A degenerate PMV is a
# legitimate quantity, but it must not be indistinguishable from one that carries
# real effect uncertainty.
#
# Note the deliberate asymmetry in the warning: a SUPPLIED beta_var of all zeros
# must NOT warn. R/30_posterior_prediction.R passes explicit zeros on purpose --
# there the marker-effect uncertainty is carried by the posterior draws, so adding
# it to each draw would double-count it. Only an ABSENT beta_var is a surprise.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(19)
n <- 12L; m <- 40L

# Inbred-by-construction dosages (0/2), then make ONE parent residually het so the
# table is forced to mix estimators. That mixture is the condition that used to be
# invisible, so the test must create it rather than avoid it.
geno <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
ids <- sprintf("P%02d", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- sprintf("M%02d", seq_len(m))
het_loci <- 3:9
geno["P01", het_loci] <- 1L

marker_map <- data.frame(marker = colnames(geno), chr = rep(1:2, each = m / 2),
                         pos_cm = rep(seq(0, 90, length.out = m / 2), 2),
                         stringsAsFactors = FALSE)
beta <- stats::setNames(c(rnorm(8, sd = 1.2), rep(0, m - 8)), colnames(geno))
y <- stats::setNames(as.numeric(geno %*% beta) + rnorm(n, sd = 0.25), ids)

# Phased haplotypes consistent with the dosages: het loci split 1/0, others 0/0 or 1/1.
phase <- array(0L, dim = c(n, 2L, m), dimnames = list(ids, NULL, colnames(geno)))
for (i in seq_len(n)) for (k in seq_len(m)) {
  d <- geno[i, k]
  if (d == 2L) phase[i, , k] <- c(1L, 1L)
  else if (d == 1L) phase[i, , k] <- c(1L, 0L)
}

fit <- ng_fit_ridge_effects(geno, y, ids = ids)
pairs <- ng_make_pairs(ids)

score <- function(effects, phased) {
  ng_score_crosses(geno = geno, effects = effects, marker_map = marker_map,
                   ids = ids, pairs = pairs, adjusted_pheno = y,
                   target = "RIL", parent_type = "ril",
                   phased_haplotypes = phased, use_cpp = FALSE)
}

# ---- 1. per-row estimator, and a table that genuinely mixes two -------------
out <- score(fit, phase)
stopifnot("variance_estimator" %in% names(out))
stopifnot("variance_window_cm" %in% names(out))

is_het_cross <- out$parent1 == "P01" | out$parent2 == "P01"
stopifnot(any(is_het_cross), !all(is_het_cross))   # the table must MIX

stopifnot(all(out$variance_estimator[is_het_cross] == "phased_het_general"))
stopifnot(all(out$variance_estimator[!is_het_cross] == "dh_recomb_aRa"))
stopifnot(length(unique(out$variance_estimator)) == 2L)

# window_cm does not apply to the dense phased kernel; the row must say so rather
# than inherit the run's nominal window.
stopifnot(all(is.infinite(out$variance_window_cm[is_het_cross])))

# ---- 2. no phase supplied -> one estimator everywhere ----------------------
# Same genotypes, but inbred parent_type so the het correction cannot engage.
geno_hom <- geno; geno_hom["P01", het_loci] <- 2L
fit_hom <- ng_fit_ridge_effects(geno_hom, y, ids = ids)
out_hom <- ng_score_crosses(geno = geno_hom, effects = fit_hom, marker_map = marker_map,
                            ids = ids, pairs = pairs, adjusted_pheno = y,
                            target = "RIL", parent_type = "inbred", use_cpp = FALSE)
stopifnot(all(out_hom$variance_estimator == "dh_recomb_aRa"))
stopifnot(length(unique(out_hom$variance_estimator)) == 1L)

# ---- 3. degenerate PMV is detectable, and says so --------------------------
stopifnot("pmv_is_degenerate" %in% names(out), "beta_var_source" %in% names(out))
stopifnot(all(out$beta_var_source == "supplied"))
stopifnot(!any(out$pmv_is_degenerate))
stopifnot(any(out$pmv > out$vpm))        # real effect uncertainty is present

no_bv <- fit; no_bv$beta_var <- NULL
w <- character(0)
out_deg <- withCallingHandlers(score(no_bv, phase),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

stopifnot(all(out_deg$pmv_is_degenerate))
stopifnot(all(out_deg$beta_var_source == "absent_zero_filled"))
stopifnot(isTRUE(all.equal(out_deg$pmv, out_deg$vpm)))   # PMV collapsed onto VPM
stopifnot(any(grepl("degenerate", w, fixed = TRUE)))

# ---- 4. a SUPPLIED all-zero beta_var is silent -----------------------------
# The posterior path does this deliberately; warning there would train users to
# ignore the warning that matters.
zero_bv <- fit
zero_bv$beta_var <- stats::setNames(rep(0, m), colnames(geno))
w2 <- character(0)
out_zero <- withCallingHandlers(score(zero_bv, phase),
  warning = function(x) { w2 <<- c(w2, conditionMessage(x)); invokeRestart("muffleWarning") })
stopifnot(all(out_zero$pmv_is_degenerate))               # still degenerate, correctly
stopifnot(all(out_zero$beta_var_source == "supplied"))   # but not a surprise
stopifnot(!any(grepl("degenerate", w2, fixed = TRUE)))   # so no warning

cat("variance_row_provenance: PASS\n")
