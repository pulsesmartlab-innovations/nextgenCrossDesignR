# ng_run_cross_prediction() must let a breeder attach per-cross COST / LOGISTIC data so
# cost_col / budget / lambda_cost / logistic_col actually work through the user path. The
# candidate crosses are generated internally, so cost is supplied via cross_cost (a data
# frame with parent1/parent2 + value columns) joined onto them by unordered parent pair.
# Requesting a cost_col that was never supplied must ERROR (not silently do nothing).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(717)
n <- 14L; m <- 60L; ids <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * m, 1, 0.5), n, m, dimnames = list(ids, sprintf("M%03d", seq_len(m))))
y <- as.numeric(gm %*% rnorm(m, 0, 0.1)) + rnorm(n)
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE)
mm  <- data.frame(SNP = colnames(gm), Chr = rep(1:4, length.out = m),
                  PosBP = rep(seq(0, 3e6, length.out = 15), 4)[seq_len(m)], stringsAsFactors = FALSE)
dir <- data.frame(trait = "yield", column = "yield", direction = "increase", stringsAsFactors = FALSE)

pr <- t(utils::combn(ids, 2))
cc <- data.frame(parent1 = pr[, 1], parent2 = pr[, 2],
                 cost = runif(nrow(pr), 5, 10), stringsAsFactors = FALSE)
set.seed(9); cc$cost[sample(nrow(cc), 10)] <- 1     # a pool of cheap crosses

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = mm, trait_direction = dir,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "Chr", map_pos_col = "PosBP",
  map_pos_cm_divisor = 1e6, prediction_mode = "trait_by_trait",
  trait_value_metric = "var_complex", duplicate_action = "none",
  n_crosses = 8L, max_uses_per_parent = 6L, optimizer = "evolution",
  write_outputs = FALSE, write_figures = FALSE, seed = 3L, ...)

sel_cost <- function(res) {
  s <- res$selected_crosses
  k <- ng_group_pair_key(s$parent1, s$parent2)
  kc <- ng_group_pair_key(res$candidate_crosses$parent1, res$candidate_crosses$parent2)
  sum(res$candidate_crosses$cost[match(k, kc)], na.rm = TRUE)
}

# 1. cost_col requested but never supplied -> ERROR (no silent no-op)
err <- tryCatch({ run(cost_col = "cost", budget = 25); "ok" }, error = function(e) "errored")
stopifnot(identical(err, "errored"))

# 2. cross_cost attaches the column; the plan respects the budget
r <- run(cross_cost = cc, cost_col = "cost", budget = 25)
stopifnot("cost" %in% names(r$candidate_crosses))
stopifnot(sel_cost(r) <= 25 + 1e-6)

# 3. a cost penalty lowers total plan cost vs none
r0 <- run(cross_cost = cc, cost_col = "cost", lambda_cost = 0)
r1 <- run(cross_cost = cc, cost_col = "cost", lambda_cost = 5)
stopifnot(sel_cost(r1) <= sel_cost(r0) + 1e-6)

cat("runner cross_cost test passed\n")
