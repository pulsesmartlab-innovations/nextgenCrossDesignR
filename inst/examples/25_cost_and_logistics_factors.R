# Cost / budget and logistic factors in mate allocation (breeder path).
#
# When each candidate cross has a cost (or a logistic/geographic penalty), a breeder folds
# those into the allocation THROUGH ng_run_cross_prediction() by supplying cross_cost -- a
# per-cross data frame (parent1/parent2 + one numeric column per factor, e.g. cost, distance).
# The runner generates the candidate crosses internally, so cross_cost is joined onto them by
# UNORDERED parent pair; then:
#   * cost_col + budget  -- a HARD cap on total plan cost (drops worst value-for-money crosses
#     to stay under budget; returns a smaller plan + warning if the budget cannot fund n_crosses).
#   * lambda_cost        -- a SOFT cost penalty.
#   * logistic_col + lambda_logistic -- a soft per-cross penalty (distance, reproductive difficulty).
# Requesting a cost_col you did not supply via cross_cost is an ERROR, never a silent no-op.
# (The lower-level ng_optimize_mating_plan() takes the same knobs if you already hold a scored table.)

library(nextgenCrossDesign)

if (!("cross_cost" %in% names(formals(ng_run_cross_prediction)))) {
  stop("This example needs a build where ng_run_cross_prediction() accepts cross_cost. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ids <- sprintf("P%02d", 1:16)
gm <- matrix(2L * rbinom(16 * 12, 1, 0.5), 16, 12, dimnames = list(ids, sprintf("M%02d", 1:12)))
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = rnorm(16, 60, 6), stringsAsFactors = FALSE)
marker_map <- data.frame(SNP = colnames(gm), Chr = rep(1:3, each = 4),
                         PosCM = rep(c(0, 20, 40, 60), 3), stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase", stringsAsFactors = FALSE)

# per-cross cost + logistic (distance) for every candidate pair
pairs <- t(utils::combn(ids, 2))
cross_cost <- data.frame(parent1 = pairs[, 1], parent2 = pairs[, 2],
                         cost = runif(nrow(pairs), 1, 5), distance = runif(nrow(pairs), 0, 10),
                         stringsAsFactors = FALSE)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = marker_map, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM",
  map_position_unit = "cM", prediction_mode = "trait_by_trait", trait_value_metric = "var_complex",
  progeny = "DH", duplicate_action = "none", n_crosses = 8L, max_crosses_per_parent = 6L,
  optimizer = "evolution", cross_cost = cross_cost,
  write_outputs = FALSE, write_figures = FALSE, seed = 1L, ...)

# helper: total of a per-cross column over the SELECTED plan (unordered-pair join)
pk  <- function(a, b) ifelse(a < b, paste(a, b), paste(b, a))
tot <- function(res, col) {
  s <- res$selected_crosses; cc <- res$candidate_crosses
  sum(cc[[col]][match(pk(s$parent1, s$parent2), pk(cc$parent1, cc$parent2))], na.rm = TRUE)
}

# --- (a) hard budget cap: total plan cost stays under budget ---
budget <- 30
r_budget <- run(cost_col = "cost", budget = budget)
cat(sprintf("budget=%.1f -> total_cost=%.2f (crosses=%d)\n", budget, tot(r_budget, "cost"),
            nrow(r_budget$selected_crosses)))
stopifnot(tot(r_budget, "cost") <= budget + 1e-6)

# --- (b) soft cost penalty: higher lambda_cost lowers total cost ---
c_lo <- tot(run(cost_col = "cost", lambda_cost = 0), "cost")
c_hi <- tot(run(cost_col = "cost", lambda_cost = 5), "cost")
cat(sprintf("soft cost: total cost lambda_cost=0 -> %.2f, lambda_cost=5 -> %.2f\n", c_lo, c_hi))
stopifnot(c_hi <= c_lo + 1e-6)

# --- (c) logistic penalty steers away from high-distance crosses ---
d_lo <- tot(run(logistic_col = "distance", lambda_logistic = 0), "distance")
d_hi <- tot(run(logistic_col = "distance", lambda_logistic = 3), "distance")
cat(sprintf("logistic: total distance lambda_logistic=0 -> %.2f, =3 -> %.2f\n", d_lo, d_hi))
stopifnot(d_hi <= d_lo + 1e-6)

message("Cost / logistics allocation (breeder path) example completed.")
