# Constrained-OCS operating point: maximize gain subject to a target group coancestry.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(414)
np <- 24L
pp <- sprintf("P%02d", seq_len(np))
cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]], stringsAsFactors = FALSE)
scores$uc_dh_gebv <- rnorm(nrow(scores), 10, 2)
L <- matrix(rnorm(np * np, 0, 0.3), np, np)
K <- crossprod(L) / np; diag(K) <- diag(K) + 1; dimnames(K) <- list(pp, pp)
scores$pair_kinship <- K[cbind(match(scores$parent1, pp), match(scores$parent2, pp))]
n_crosses <- 15L

# Frontier gives the achievable coancestry range.
fr <- ng_pareto_mate_allocation(scores, n_crosses, parent_K = K)$frontier
co_rng <- range(fr$group_coancestry, na.rm = TRUE)
lo_t <- co_rng[1] + 0.35 * diff(co_rng)   # a tight-ish, achievable target
hi_t <- co_rng[1] + 0.80 * diff(co_rng)   # a looser target

# --- 1. constraint is respected: achieved coancestry <= target ---
p_lo <- ng_optimize_mating_plan(scores, n_crosses, parent_K = K, target_coancestry = lo_t)
p_hi <- ng_optimize_mating_plan(scores, n_crosses, parent_K = K, target_coancestry = hi_t)
s_lo <- attr(p_lo, "summary"); s_hi <- attr(p_hi, "summary")
stopifnot(nrow(p_lo) == n_crosses, nrow(p_hi) == n_crosses)
stopifnot(s_lo$group_coancestry <= lo_t + 1e-9)
stopifnot(s_hi$group_coancestry <= hi_t + 1e-9)
stopifnot(identical(s_lo$target_coancestry_status, "met") ||
            identical(s_lo$target_coancestry_status, "met_slack"))

# --- 2. monotone: a looser coancestry cap yields at least as much gain ---
stopifnot(s_hi$mean_gain >= s_lo$mean_gain - 1e-8)

# --- 3. summary records the constraint machinery ---
stopifnot(isTRUE(all.equal(s_lo$target_coancestry, lo_t)))
stopifnot(is.finite(s_lo$achieved_coancestry), is.finite(s_lo$selected_lambda_group))

# --- 4. infeasible target (below the minimum achievable) -> warning + min-coancestry plan ---
too_tight <- co_rng[1] - 0.5 * diff(co_rng)
w <- NULL
p_inf <- withCallingHandlers(
  ng_optimize_mating_plan(scores, n_crosses, parent_K = K, target_coancestry = too_tight),
  warning = function(c) { w <<- conditionMessage(c); invokeRestart("muffleWarning") })
s_inf <- attr(p_inf, "summary")
stopifnot(identical(s_inf$target_coancestry_status, "infeasible"))
stopifnot(!is.null(w), grepl("below the minimum achievable", w))
# it returns the lowest-coancestry plan available
stopifnot(s_inf$group_coancestry <= s_lo$group_coancestry + 1e-9)

# --- 5. additivity: NOT setting target_coancestry leaves the plain path unchanged ---
p_plain <- ng_optimize_mating_plan(scores, n_crosses, parent_K = K, lambda_group = 1)
stopifnot(is.null(attr(p_plain, "summary")$target_coancestry))
# and target_coancestry takes precedence over the emphasis dial (absolute constraint wins)
p_both <- ng_optimize_mating_plan(scores, n_crosses, parent_K = K,
                                  target_coancestry = lo_t, strategy = "high_gain")
stopifnot(!is.null(attr(p_both, "summary")$target_coancestry_status))
stopifnot(attr(p_both, "summary")$group_coancestry <= lo_t + 1e-9)

cat("target coancestry (constrained OCS) test passed\n")
