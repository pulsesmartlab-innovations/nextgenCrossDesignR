# Combined breeder controls in one run: gain + cost + management + strategy.
#
# A real mate-design run rarely tunes one lever at a time. This example drives the ONE
# user-facing call ng_run_cross_prediction() with several recently-added controls AT ONCE,
# the way a breeder actually would:
#   * GAIN metric  -- trait_value_metric = "var_complex" (usefulness).
#   * OPTIMIZER    -- "evolution" (recommended default; see example 06).
#   * STRATEGY dial -- strategy = "balanced": manage the gain-vs-diversity trade-off explicitly
#                     (the recommended way to sustain gain over recurrent cycles; see example 21).
#   * COST / budget -- cross_cost + cost_col + budget: keep the plan under a crossing budget.
#   * MANAGEMENT    -- min_crosses_per_parent (batch economics); parent_group + group_permission
#                     (only cross between permitted heterotic pools); lambda_progeny_inbreeding
#                     (limit expected progeny inbreeding).
# plan_summary reports what each control achieved.

library(nextgenCrossDesign)

if (!("cross_cost" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build where ng_run_cross_prediction() accepts cross_cost. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260704)
ids <- sprintf("P%02d", 1:20)
gm <- matrix(2L * rbinom(20 * 20, 1, 0.5), 20, 20, dimnames = list(ids, sprintf("M%02d", 1:20)))
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = rnorm(20, 60, 6), stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = colnames(gm), Chr = rep(1:4, each = 5),
                         PosCM = rep(seq(0, 80, length.out = 5), 4), stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase", stringsAsFactors = FALSE)

# two heterotic pools (A, B); only cross ACROSS pools
grp  <- stats::setNames(rep(c("A", "B"), each = 10), ids)
perm <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))

# per-cross cost for every candidate pair
pairs <- t(utils::combn(ids, 2))
cross_cost <- data.frame(parent1 = pairs[, 1], parent2 = pairs[, 2],
                         cost = runif(nrow(pairs), 1, 4), stringsAsFactors = FALSE)

result <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM",
  map_position_unit = "cM", prediction_mode = "trait_by_trait",
  trait_value_metric = "var_complex",                 # GAIN metric
  progeny = "DH", optimizer = "evolution",            # recommended optimizer
  n_crosses = 12L, max_uses_per_parent = 4L,
  strategy = "balanced",                              # STRATEGY dial
  cross_cost = cross_cost, cost_col = "cost", budget = 30,   # COST + budget
  min_crosses_per_parent = 2L,                        # MANAGEMENT: batch economics
  parent_group = grp, group_permission = perm,        # MANAGEMENT: only A x B
  lambda_progeny_inbreeding = 20,                     # MANAGEMENT: limit progeny inbreeding
  duplicate_action = "none", write_outputs = FALSE, write_figures = FALSE, seed = 1L)

s <- result$selected_crosses
pk <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
sel_cost <- sum(result$candidate_crosses$cost[
  match(pk(s$parent1, s$parent2), pk(result$candidate_crosses$parent1, result$candidate_crosses$parent2))],
  na.rm = TRUE)
cnt <- table(c(s$parent1, s$parent2))

# Hard constraints the package always enforces (never silently violated):
stopifnot(all(grp[s$parent1] != grp[s$parent2]))    # legality: only across-pool (A x B) crosses
stopifnot(sel_cost <= 30 + 1e-6)                     # budget: total plan cost stays under cap

# min-use-if-used is a SOFT target here: strategy = "balanced" spreads parents for diversity while
# min_crosses_per_parent concentrates them, and the budget caps the set -- so when these pull
# against each other the package returns the best feasible plan and warns, rather than forcing an
# infeasible one. Report what was achieved (that trade-off IS the point of tuning several levers).
ps <- result$plan_summary
cat(sprintf("selected %d crosses | all A x B: %s | total_cost=%.2f (<=30) | min per-parent use=%d\n",
            nrow(s), all(grp[s$parent1] != grp[s$parent2]), sel_cost, min(as.integer(cnt))))
cat(sprintf("strategy=%s | mean_progeny_inbreeding=%.4f | group_coancestry=%.4f\n",
            ps$strategy, ps$mean_progeny_inbreeding, ps$group_coancestry))

message("Combined breeder controls example completed.")
