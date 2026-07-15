# Manage and report progeny inbreeding.
#
# Every scored cross carries expected_progeny_inbreeding (the coancestry of the two
# parents = the expected inbreeding of their progeny). You can:
#   * penalize it in allocation with lambda_progeny_inbreeding (through
#     ng_run_cross_prediction), and read mean/max progeny F back from plan_summary;
#   * export a progeny-inbreeding histogram (JSON) for a report/UI.
#
# NOTE: lambda_progeny_inbreeding and the existing lambda_mating act on the SAME axis
# (parent-pair relatedness) and their penalties ADD -- prefer setting ONE of them.
# lambda_group is a separate, population-level diversity control.

library(nextgenCrossDesign)

if (!("lambda_progeny_inbreeding" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build with lambda_progeny_inbreeding. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ids <- sprintf("P%02d", 1:16)
geno <- as.data.frame(matrix(2L * rbinom(16 * 12, 1, 0.5), 16, 12,
                             dimnames = list(NULL, sprintf("M%02d", 1:12))))
genotype <- cbind(NAME = ids, geno)
phenotype <- data.frame(NAME = ids, yield = rnorm(16, 60, 6), stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = colnames(geno), Chr = rep(1:3, each = 4),
                         PosCM = rep(c(0, 20, 40, 60), 3), stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM", map_position_unit = "cM",
  prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  n_crosses = 8L, max_uses_per_parent = 4L, use_ocs = TRUE,
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L, ...)

# --- (a) penalize progeny inbreeding: higher lambda -> lower mean progeny F ---
base <- run(lambda_progeny_inbreeding = 0, lambda_mating = 0, lambda_group = 0)
pen  <- run(lambda_progeny_inbreeding = 40, lambda_mating = 0, lambda_group = 0)
cat(sprintf("mean progeny inbreeding: baseline=%.4f  penalized=%.4f\n",
            base$plan_summary$mean_progeny_inbreeding, pen$plan_summary$mean_progeny_inbreeding))
stopifnot(pen$plan_summary$mean_progeny_inbreeding <=
            base$plan_summary$mean_progeny_inbreeding + 1e-8)

# --- (b) the per-cross column is available in the candidate/selected tables ---
stopifnot("expected_progeny_inbreeding" %in% names(base$candidate_crosses) ||
            "pair_kinship" %in% names(base$candidate_crosses))

# --- (c) export a progeny-inbreeding histogram (JSON) for a report/UI ---
geno_mat <- as.matrix(geno); rownames(geno_mat) <- ids
fit <- ng_fit_ridge_effects(geno_mat, phenotype$yield, ids = ids, seed = 1L)
mm <- data.frame(marker = colnames(geno_mat), chr = marker_map$Chr, pos_cm = marker_map$PosCM)
scores <- ng_score_crosses(geno_mat, fit, marker_map = mm, ids = ids,
                           adjusted_pheno = setNames(phenotype$yield, ids), target = "DH")
hist_df <- ng_progeny_inbreeding_histogram(scores, breaks = 10L)
stopifnot(sum(hist_df$count) == nrow(scores))
hist_json <- file.path(tempdir(), "progeny_inbreeding_histogram.json")
ng_write_progeny_inbreeding_histogram_json(scores, output_path = hist_json, breaks = 10L)
stopifnot(file.exists(hist_json))

message("Progeny-inbreeding management example completed.")
message("Histogram JSON: ", hist_json)
