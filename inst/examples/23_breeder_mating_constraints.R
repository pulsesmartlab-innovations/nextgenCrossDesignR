# Breeder mating constraints: min-use-if-used, committed matings, mating-group rules.
#
# All are passed to ng_run_cross_prediction() and applied in the fast repair stage:
#   * min_crosses_per_parent -- a parent that is used at all must be used >= N times
#     (else it is dropped to 0); anti-singleton batch economics.
#   * committed_crosses -- matings you have already decided are locked into the plan and
#     the rest is optimized around them.
#   * parent_group + group_permission -- only allow crosses between permitted groups
#     (e.g. heterotic pools, or male x female); + group_quota caps crosses per group pair.
#
# If a constraint makes a full-size plan impossible the package returns a smaller feasible
# plan and warns (it never silently violates the constraint).

library(nextgenCrossDesign)

if (!all(c("min_crosses_per_parent", "committed_crosses", "group_permission") %in%
         names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build with the breeder mating constraints. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ids <- sprintf("P%02d", 1:20)
geno <- as.data.frame(matrix(2L * rbinom(20 * 12, 1, 0.5), 20, 12,
                             dimnames = list(NULL, sprintf("M%02d", 1:12))))
genotype <- cbind(NAME = ids, geno)
phenotype <- data.frame(NAME = ids, yield = rnorm(20, 60, 6), stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = colnames(geno), Chr = rep(1:3, each = 4),
                         PosCM = rep(c(0, 20, 40, 60), 3), stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM", map_position_unit = "cM",
  prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  n_crosses = 12L, max_crosses_per_parent = 6L, use_ocs = TRUE,
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L, ...)

# --- (a) min-use-if-used: every used parent appears >= 3 times ---
r_minuse <- run(min_crosses_per_parent = 3L)
cnt <- table(c(r_minuse$selected_crosses$parent1, r_minuse$selected_crosses$parent2))
stopifnot(all(as.integer(cnt) >= 3L))
cat("min-use-if-used: min per-parent use =", min(as.integer(cnt)), "\n")

# --- (b) committed matings: a decided cross is always in the plan ---
committed <- data.frame(parent1 = "P01", parent2 = "P02", stringsAsFactors = FALSE)
r_commit <- run(committed_crosses = committed)
key <- paste(pmin(r_commit$selected_crosses$parent1, r_commit$selected_crosses$parent2),
             pmax(r_commit$selected_crosses$parent1, r_commit$selected_crosses$parent2))
stopifnot(paste("P01", "P02") %in% key)
cat("committed matings: P01xP02 present =", paste("P01", "P02") %in% key, "\n")

# --- (c) mating-group permissions: only cross group M with group F ---
grp <- setNames(rep(c("M", "F"), length.out = 20), ids)
perm <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2, 2, dimnames = list(c("M", "F"), c("M", "F")))
r_group <- run(parent_group = grp, group_permission = perm)
g1 <- grp[r_group$selected_crosses$parent1]; g2 <- grp[r_group$selected_crosses$parent2]
stopifnot(all(g1 != g2))   # no within-group matings survived the legality filter
cat("group permissions: all crosses are M x F =", all(g1 != g2), "\n")

message("Breeder mating constraints example completed.")
