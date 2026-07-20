# Gain-diversity balancing: pick the operating point on the mate-allocation frontier.
#
# Instead of hand-tuning the raw coancestry penalty (lambda_group), you can steer the
# gain-vs-diversity trade-off two ways, both through ng_run_cross_prediction():
#   * strategy / diversity_emphasis -- a relative dial ("high_gain" | "balanced" |
#     "diversity", or a number 0 = all gain .. 100 = all diversity).
#   * target_coancestry -- an ABSOLUTE cap (constrained OCS / Meuwissen): maximize gain
#     at or under a target group coancestry (~ a rate-of-inbreeding target).
# This script also exports the frontier as JSON for a decision UI.
#
# WHY it matters (evidence, config-scoped): in a recurrent RIL program pure usefulness metrics
# exhaust genetic variance and gain plateaus, so managing diversity EXPLICITLY is the recommended
# true-use pattern. In the 10-rep x 25-cycle study (tools/run_method_optimizer_metric_study.R) the
# balanced strategy dial (strategy = "balanced") was the 2nd-best arm for cycle-25 realized gain --
# just behind le and clearly ahead of the pure-usefulness cluster -- because it keeps
# variance in the tank. Use the dial or target_coancestry on top of a good merit metric (example 27
# combines it with cost/management), rather than chasing pure usefulness.

library(nextgenCrossDesign)

if (!("target_coancestry" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a nextgenCrossDesign build with the gain-diversity controls ",
       "(strategy / diversity_emphasis / target_coancestry). Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ids <- sprintf("P%02d", 1:16)
# inbred (0/2) candidate lines, 12 markers on 3 chromosomes
geno <- as.data.frame(matrix(2L * rbinom(16 * 12, 1, 0.5), 16, 12,
                             dimnames = list(NULL, sprintf("M%02d", 1:12))))
genotype <- cbind(NAME = ids, geno)
phenotype <- data.frame(NAME = ids,
                        yield = rnorm(16, 60, 6), disease = rnorm(16, 4, 1),
                        stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = colnames(geno), Chr = rep(1:3, each = 4),
                         PosCM = rep(c(0, 20, 40, 60), 3), stringsAsFactors = FALSE)
direction <- data.frame(trait = c("yield", "disease"), column = c("yield", "disease"),
                        direction = c("increase", "decrease"), stringsAsFactors = FALSE)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM", map_position_unit = "cM",
  prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  multi_trait_method = "weighted", trait_weights = c(1, 1),
  n_crosses = 8L, max_crosses_per_parent = 4L, use_ocs = TRUE,
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L, ...)

# --- (a) strategy dial: high_gain vs balanced vs diversity ---
hi <- run(strategy = "high_gain")
ba <- run(strategy = "balanced")
di <- run(strategy = "diversity")
cat("strategy       gain   group_coancestry\n")
for (r in list(hi, ba, di)) {
  s <- r$plan_summary
  cat(sprintf("  %-11s  %6.3f   %6.4f\n", s$strategy, s$mean_gain, s$group_coancestry))
}
# high_gain should give >= gain and >= coancestry than diversity
stopifnot(hi$plan_summary$mean_gain >= di$plan_summary$mean_gain - 1e-8)
stopifnot(di$plan_summary$group_coancestry <= hi$plan_summary$group_coancestry + 1e-8)

# --- (b) numeric diversity_emphasis (0 = all gain .. 100 = all diversity) ---
em <- run(diversity_emphasis = 70)
stopifnot(is.finite(em$plan_summary$achieved_emphasis))

# --- (c) constrained OCS: maximize gain at/under a target group coancestry ---
# choose a target inside the achievable range from the strategy runs above
tc <- mean(c(hi$plan_summary$group_coancestry, di$plan_summary$group_coancestry))
con <- run(target_coancestry = tc)
sc <- con$plan_summary
cat(sprintf("\nconstrained OCS: target=%.4f achieved=%.4f status=%s\n",
            sc$target_coancestry, sc$group_coancestry, sc$target_coancestry_status))
stopifnot(sc$group_coancestry <= tc + 1e-9 || identical(sc$target_coancestry_status, "met_slack"))

# --- (d) export the gain-vs-diversity frontier as JSON for a decision UI ---
# ng_write_frontier_json works on a scored candidate table + parent kinship, which you
# can get from the lower-level scoring API:
geno_mat <- as.matrix(geno); rownames(geno_mat) <- ids
fit <- ng_fit_ridge_effects(geno_mat, phenotype$yield, ids = ids, seed = 1L)
mm <- data.frame(marker = colnames(geno_mat), chr = marker_map$Chr, pos_cm = marker_map$PosCM)
scores <- ng_score_crosses(geno_mat, fit, marker_map = mm, ids = ids,
                           adjusted_pheno = setNames(phenotype$yield, ids), target = "DH")
K <- ng_parent_kinship(geno_mat)
frontier_json <- file.path(tempdir(), "mating_frontier.json")
ng_write_frontier_json(scores, n_crosses = 8L, parent_kinship = K, output_path = frontier_json)
stopifnot(file.exists(frontier_json))

message("Gain-diversity balance + constrained OCS example completed.")
message("Frontier JSON: ", frontier_json)
