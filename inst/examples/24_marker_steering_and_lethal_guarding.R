# Marker steering and lethal-allele guarding.
#
# Two capabilities:
#   * Lethal-allele guarding (user API): pass lethal_spec to ng_run_cross_prediction() to
#     drop carrier x carrier matings at known deleterious recessive loci.
#   * Marker steering (scores-level API): drive the frequency of a favourable major-gene
#     allele up/down (or toward a target) with ng_marker_target_spec() +
#     ng_apply_marker_management(), then allocate with ng_optimize_mating_plan().

library(nextgenCrossDesign)

if (!("lethal_spec" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build with lethal_spec / marker management. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ids <- sprintf("P%02d", 1:16)
mk <- sprintf("M%02d", 1:12)
gm <- matrix(2L * rbinom(16 * 12, 1, 0.5), 16, 12, dimnames = list(ids, mk))
# lethal locus M12: P01 and P02 carry the risk allele (homozygous, dosage 2), rest safe (0)
gm[, "M12"] <- 0L; gm[c("P01", "P02"), "M12"] <- 2L
# steer locus M01: a MINORITY favourable allele (only 5 carriers) we want to increase, so
# there is genuine room for steering to raise its frequency vs a gain-only plan
gm[, "M01"] <- 0L; gm[1:5, "M01"] <- 2L
genotype <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = rnorm(16, 60, 6), stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = mk, Chr = rep(1:3, each = 4),
                         PosCM = rep(c(0, 20, 40, 60), 3), stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

# --- (a) lethal-allele guarding through the user API ---
lspec <- ng_lethal_recessive_spec("M12", risk_allele = "alt")
res <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM", map_position_unit = "cM",
  prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  n_crosses = 8L, max_crosses_per_parent = 4L, use_ocs = TRUE, lethal_spec = lspec,
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L)
key <- paste(pmin(res$selected_crosses$parent1, res$selected_crosses$parent2),
             pmax(res$selected_crosses$parent1, res$selected_crosses$parent2))
stopifnot(!(paste("P01", "P02") %in% key))   # carrier x carrier excluded
cat("lethal guarding: carrier x carrier (P01xP02) excluded =",
    !(paste("P01", "P02") %in% key), "\n")

# --- (b) marker steering at the scores level (drive M01's ALT allele up) ---
fit <- ng_fit_ridge_effects(gm, phenotype$yield, ids = ids, seed = 1L)
mm <- data.frame(marker = mk, chr = marker_map$Chr, pos_cm = marker_map$PosCM)
scores <- ng_score_crosses(gm, fit, marker_map = mm, ids = ids,
                           adjusted_pheno = setNames(phenotype$yield, ids), target = "DH")
mspec <- ng_marker_target_spec("M01", direction = "increase", weight = 1)
# attach the steering score + expected progeny allele frequency, and blend into a gain col
aug <- ng_apply_marker_management(scores, geno = gm, marker_target_spec = mspec,
                                  gain_col = "usefulness_pmv_gebv", lambda_marker = 0.5)
stopifnot("marker_target_score" %in% names(aug), "marker_adjusted_gain" %in% names(aug))
K <- ng_parent_kinship(gm)
plan_steered <- ng_optimize_mating_plan(aug, n_crosses = 8L, parent_kinship = K,
                                        gain_col = "marker_adjusted_gain")
# the steered plan should favour crosses with a higher marker-target score (and higher
# expected M01 frequency) than a gain-only plan
plan_gain <- ng_optimize_mating_plan(aug, n_crosses = 8L, parent_kinship = K, gain_col = "usefulness_pmv_gebv")
cat(sprintf("mean marker_target_score:      steered=%.3f  gain-only=%.3f\n",
            mean(plan_steered$marker_target_score), mean(plan_gain$marker_target_score)))
cat(sprintf("mean expected M01 progeny freq: steered=%.3f  gain-only=%.3f\n",
            mean(plan_steered$marker_freq_M01), mean(plan_gain$marker_freq_M01)))
stopifnot(mean(plan_steered$marker_target_score) >= mean(plan_gain$marker_target_score) - 1e-8)

message("Marker steering + lethal-allele guarding example completed.")
