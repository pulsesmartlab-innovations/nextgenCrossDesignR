# Any-ploidy polyploid mate design in one call: ng_polyploid_design_crosses().
#
# The package separates scoring from allocation, and the allocator is ploidy-agnostic, so a
# polyploid breeder gets the FULL native control suite (strategy dial, target_coancestry,
# committed matings, group quotas, cost/budget, the evolution optimizer) -- for ANY ploidy --
# from a single entry point that also runs ploidy-aware QC and uses a correct allele-frequency
# polyploid GRM (ng_polyploid_grm). Scoring is MEAN-based (mid-parent GEBV); per VALIDATED_STATE.md
# the mean criterion is at least as good as usefulness for realized gain, and it needs no
# simulation and no heavy dependency.
#
# Inputs: an allele-dosage matrix (0..ploidy, parents in rows, markers in columns) plus either
# per-marker additive effects OR a per-parent phenotype (effects are then estimated by ridge).

library(nextgenCrossDesign)

if (!"ng_polyploid_design_crosses" %in% getNamespaceExports("nextgenCrossDesign")) {
  stop("This example needs a build with ng_polyploid_design_crosses(). Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ploidy <- 4L                                  # autotetraploid; set 2 for diploid, 6 for hexaploid
np <- 24L; nm <- 40L
ids <- sprintf("clone%02d", seq_len(np))
dosage <- matrix(rbinom(np * nm, ploidy, 0.5), np, nm,
                 dimnames = list(ids, sprintf("snp%02d", seq_len(nm))))
# A per-parent phenotype (the breeder's trait record); effects are estimated inside the call.
true_effects <- rnorm(nm)
phenotype <- as.numeric(dosage %*% true_effects) + rnorm(np, 0, 2)
names(phenotype) <- ids

# One call: QC -> estimate effects -> score (mean + ploidy GRM) -> native OCS allocation with a
# diversity-aware strategy and a committed mating.
plan <- ng_polyploid_design_crosses(
  dosage = dosage, n_crosses = 15L, ploidy = ploidy,
  phenotype = phenotype,
  max_crosses_per_parent = 4L,
  run_qc = TRUE,                                          # ploidy-aware QC (0..ploidy, MAF, missing)
  strategy = "balanced",                                  # gain-vs-diversity dial (ploidy-general)
  committed_crosses = data.frame(parent1 = "clone01", parent2 = "clone02")
)

stopifnot(nrow(plan) == 15L, identical(attr(plan, "ploidy"), ploidy))
cat("Polyploid (", ploidy, "x) mate design: ", nrow(plan), " crosses.\n", sep = "")
plan_pair_key <- paste(pmin(plan$parent1, plan$parent2), pmax(plan$parent1, plan$parent2), sep = "_")
cat("Committed clone01xclone02 present: ",
    "clone01_clone02" %in% plan_pair_key, "\n", sep = "")

sm <- attr(plan, "summary")
cat(sprintf("Mean cross value = %.3f | group coancestry = %.4f\n", sm$mean_gain, sm$group_coancestry))
print(utils::head(plan[, c("parent1", "parent2", "poly_mean")], 4))

qc <- attr(plan, "qc")
cat(sprintf("QC: %d samples x %d markers, %d monomorphic dropped, %d markers kept\n",
            qc$n_samples, qc$n_markers, qc$n_monomorphic, qc$n_markers_clean))

message("Polyploid mate design example completed.")
