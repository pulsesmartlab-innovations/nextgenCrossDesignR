# Native evolutionary (memetic genetic-algorithm) mate-allocation optimizer.
#
# In addition to greedy_local / repair_local / mip_linear / mip_contribution, the package
# ships a native memetic GA optimizer. Select it with optimizer = "evolution" (aliases
# "ga" / "de" / "memetic") in ng_run_cross_prediction(). It optimizes the same OCS
# objective as the MIP/greedy paths, warm-started from greedy with elitism (so it is never
# worse than greedy) and a C++ local-swap hill-climb, and is tuned with
# evol_solutions / evol_iterations / evol_stop / evol_seed.
#
# WHY it is the recommended default for real programs: on a single-decision objective bake-off
# MIP wins and evolution is a close second, but over a recurrent RIL program (10 reps x 25
# cycles, tools/run_method_optimizer_metric_study.R) evolution MATCHES exact MIP on realized
# gain (cycle-25 7.20 +/- 0.30 vs 7.07 +/- 0.44, not significant) while scaling to large
# candidate sets where MIP is size/time-guarded; both clearly beat greedy (6.78 +/- 0.30).

library(nextgenCrossDesign)

if (!("evol_solutions" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build with the native evolution optimizer. Reinstall the current tarball.",
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

run <- function(optimizer, ...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM", map_position_unit = "cM",
  prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  n_crosses = 8L, max_crosses_per_parent = 4L, use_ocs = TRUE,
  lambda_group = 0.05, lambda_mating = 0.02, optimizer = optimizer,
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L, ...)

# --- native evolution optimizer ---
evo <- run("evolution", evol_solutions = 40L, evol_iterations = 100L, evol_stop = 30L)
stopifnot(inherits(evo, "ng_cross_prediction_result"))
stopifnot(nrow(evo$selected_crosses) == 8L)

# compare its achieved OCS objective to greedy (evolution is warm-started from greedy with
# elitism, so it should be >= greedy)
grd <- run("greedy_local")
cat(sprintf("achieved OCS objective: evolution=%.4f  greedy=%.4f\n",
            evo$plan_summary$total_gain, grd$plan_summary$total_gain))

message("Native evolution optimizer example completed.")
