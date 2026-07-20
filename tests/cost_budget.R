# Cost / budget and logistic penalty factors.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(505)
np <- 24L
pp <- sprintf("P%02d", seq_len(np))
cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]],
                     stringsAsFactors = FALSE)
scores$usefulness_pmv_gebv <- rnorm(nrow(scores), 10, 2)
scores$pair_kinship <- 0
scores$cost <- runif(nrow(scores), 1, 5)
K <- diag(np); dimnames(K) <- list(pp, pp)
n_crosses <- 15L

# --- 1. hard budget: total cost <= budget ---
budget <- 40
plan <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                max_crosses_per_parent = 6L,
                                cost_col = "cost", budget = budget,
                                method = "greedy_local")
s <- attr(plan, "summary")
stopifnot(s$total_cost <= budget + 1e-9)
stopifnot(isFALSE(s$over_budget))

# --- 2. soft cost penalty: higher lambda_cost lowers mean plan cost ---
p_lo <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                max_crosses_per_parent = 6L,
                                cost_col = "cost", lambda_cost = 0, method = "greedy_local")
p_hi <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                max_crosses_per_parent = 6L,
                                cost_col = "cost", lambda_cost = 5, method = "greedy_local")
c_lo <- sum(p_lo$cost); c_hi <- sum(p_hi$cost)
stopifnot(c_hi <= c_lo + 1e-9)

# --- 3. logistic penalty steers away from high-penalty crosses ---
scores$distance <- runif(nrow(scores), 0, 10)
pl0 <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                               max_crosses_per_parent = 6L,
                               logistic_col = "distance", lambda_logistic = 0,
                               method = "greedy_local")
pl1 <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                               max_crosses_per_parent = 6L,
                               logistic_col = "distance", lambda_logistic = 3,
                               method = "greedy_local")
stopifnot(sum(pl1$distance) <= sum(pl0$distance) + 1e-9)

# --- 4. infeasible budget returns a shorter plan with a warning ---
tight <- 12  # far below the cost of 15 crosses (min ~15)
w <- NULL
plan_t <- withCallingHandlers(
  ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K, max_crosses_per_parent = 6L,
                          cost_col = "cost", budget = tight, method = "greedy_local"),
  warning = function(cond) { w <<- conditionMessage(cond); invokeRestart("muffleWarning") })
st <- attr(plan_t, "summary")
stopifnot(nrow(plan_t) < n_crosses, st$total_cost <= tight + 1e-9)
stopifnot(!is.null(w), grepl("constraint", w, ignore.case = TRUE))

cat("cost/budget test passed\n")
