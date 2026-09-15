# The run must report what it DID, not what it was asked for.
#
# Two defects, both the same shape: a fact the engine computes and then discards.
#
# DEFECT A -- settings echoes the RAW USER STRING. `trait_value_metric =
# "var_complex"` is rewritten to usefulness + pmv before anything runs, yet
# settings$trait_value_metric still reports "var_complex". Several friendly aliases
# behave the same way (family_variance -> vpm, reliable_family_variance -> pmv,
# le -> parent_distance). A reader cannot tell from the artifact which metric
# actually produced the numbers, and `uc_variance_source` is absent from
# effect_summary entirely -- so per trait you could tell fast-vs-full PMV but not
# PMV-vs-VPM-vs-distance.
#
# DEFECT B -- mean_source_criterion has ZERO consumers. ng_choose_mean_source()
# constructs it in all four return branches and nothing reads it, so the package
# records WHICH basis was chosen while discarding WHY. Those are different
# questions with different remedies:
#
#   "cv_predictive_r2"      the markers cleared the bar
#   "below_threshold"       they were measured and lost
#   "no_phenotype_fallback" GEBV used because nothing else was available
#
# Asserted on the delivered artifact, and -- for the metric -- on the stronger
# invariant that the column name the run REPORTS must actually resolve in the
# table it delivered. A name that does not resolve is worse than no name.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(31)
n <- 30L; m <- 40L
geno <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
ids <- sprintf("P%02d", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- sprintf("M%02d", seq_len(m))
beta <- stats::setNames(c(rnorm(8, sd = 1.5), rep(0, m - 8)), colnames(geno))
g_true <- as.numeric(geno %*% beta)
y <- g_true + rnorm(n, sd = 0.3 * stats::sd(g_true))

gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
ph  <- data.frame(NAME = ids, YIELD = y, stringsAsFactors = FALSE)
map <- data.frame(SNP = colnames(geno), chr = rep(1:2, each = m / 2),
                  cm = rep(seq(0, 60, length.out = m / 2), 2), stringsAsFactors = FALSE)
dir <- data.frame(Trait = "YIELD", Selection_direction = "increase",
                  stringsAsFactors = FALSE)

run <- function(metric = NULL, ucvs = NULL, thresh = 0.01) {
  args <- list(genotype = gdf, genotype_id_col = "NAME",
               phenotype = ph, phenotype_id_col = "NAME", traits_to_use = "YIELD",
               trait_direction = dir, direction_trait_col = "Trait",
               direction_column_col = "Trait",
               direction_direction_col = "Selection_direction",
               marker_map = map, map_marker_col = "SNP", map_chr_col = "chr",
               map_pos_cm_col = "cm", map_position_unit = "cM",
               progeny = "DH", parent_type = "inbred", n_crosses = 8L,
               write_outputs = FALSE, write_figures = FALSE,
               run_posterior_prediction = FALSE,
               min_cv_predictive_r2 = thresh, seed = 20260907L)
  if (!is.null(metric)) args$trait_value_metric <- metric
  if (!is.null(ucvs))   args$uc_variance_source <- ucvs
  do.call(ng_run_cross_prediction, args)
}

# ---- A. the deprecated alias resolves, and the artifact says so ------------
r <- run(metric = "var_complex")
stopifnot(identical(r$settings$trait_value_metric, "usefulness"))
stopifnot(identical(r$settings$uc_variance_source, "pmv"))
stopifnot(identical(r$settings$trait_value_metric_input, "var_complex"))

# The reported column must RESOLVE in the delivered table. A reported name that
# does not exist in the artifact is the defect this whole release is about.
vc <- r$effect_summary$variance_column_used
stopifnot(is.character(vc), nzchar(vc))
stopifnot(paste0("YIELD_", vc) %in% names(r$candidate_crosses))
stopifnot(identical(vc, "pmv"))

# ---- the canonical spelling reports itself unchanged -----------------------
r2 <- run(metric = "usefulness", ucvs = "pmv")
stopifnot(identical(r2$settings$trait_value_metric, "usefulness"))
stopifnot(identical(r2$settings$trait_value_metric_input, "usefulness"))

# ---- a different variance choice must be visible per trait -----------------
# uc_variance_source was absent from effect_summary entirely, so PMV-vs-VPM was
# unreportable even though fast-vs-full PMV was reported.
r3 <- run(metric = "usefulness", ucvs = "vpm")
stopifnot(identical(r3$settings$uc_variance_source, "vpm"))
stopifnot(identical(r3$effect_summary$variance_column_used, "vpm"))
stopifnot(paste0("YIELD_", r3$effect_summary$variance_column_used) %in% names(r3$candidate_crosses))
# ... and it must differ from the PMV run, or the field is decorative.
stopifnot(!identical(r$effect_summary$variance_column_used,
                     r3$effect_summary$variance_column_used))

# ---- B. the REASON for the mean-source decision survives -------------------
hi <- run(thresh = 0.01)     # markers clear the bar
lo <- run(thresh = 0.99)     # markers measured and lost

stopifnot("mean_source_criterion" %in% names(hi$effect_summary))
stopifnot(identical(hi$effect_summary$mean_source[[1L]], "GEBV"))
stopifnot(identical(hi$effect_summary$mean_source_criterion[[1L]], "cv_predictive_r2"))

stopifnot(identical(lo$effect_summary$mean_source[[1L]], "adjusted_pheno"))
stopifnot(identical(lo$effect_summary$mean_source_criterion[[1L]], "below_threshold"))

# The criterion must always be populated -- all four branches of
# ng_choose_mean_source() set it, so an NA here means it was lost in transport.
stopifnot(!is.na(hi$effect_summary$mean_source_criterion[[1L]]))
stopifnot(!is.na(lo$effect_summary$mean_source_criterion[[1L]]))

# It belongs on effect_summary and NOT on the cross table. The mean basis is a
# per-trait decision, so repeating it across every candidate cross would be 5,050
# identical values pretending to be per-row information -- the cross table
# deliberately carries only parent identifiers, kinship, and per-trait values.
stopifnot(!("mean_source_criterion" %in% names(hi$candidate_crosses)))

# The scored table (pre-assembly) is where it is produced, and it must be there
# for effect_summary to have been able to read it.
sc <- ng_score_crosses(geno = geno, effects = ng_fit_ridge_effects(geno, y, ids = ids),
                       marker_map = data.frame(marker = colnames(geno), chr = map$chr,
                                               pos_cm = map$cm, stringsAsFactors = FALSE),
                       ids = ids, pairs = ng_make_pairs(ids)[1:5, , drop = FALSE],
                       adjusted_pheno = y, target = "DH", parent_type = "inbred",
                       use_cpp = FALSE, min_cv_predictive_r2 = 0.01)
stopifnot("mean_source_criterion" %in% names(sc))
stopifnot(identical(sc$mean_source_criterion[[1L]], "cv_predictive_r2"))

cat("settings_report_resolved_method: PASS\n")
