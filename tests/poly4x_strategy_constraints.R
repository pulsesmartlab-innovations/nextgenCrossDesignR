# Phase 0 parity: the polyploid policy wrappers (ng_poly4x_ocs / ng_poly4x_policy in
# R/16, ng_polyploid_policy in R/18) must forward the shared mate-selection controls
# (strategy / diversity_emphasis / target_coancestry / committed_crosses /
# method = "evolution") straight through to the generic allocator. This test drives
# them on a synthetic poly4x-style score table (no AlphaSimR needed) so it is fast and
# deterministic, and asserts each forwarded control actually takes effect.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(4242)
np <- 18L
pp <- sprintf("P%02d", seq_len(np))

# Relationship matrix with real structure so the gain-vs-coancestry frontier is non-flat.
k_dim <- 5L
loadings <- matrix(rnorm(np * k_dim), np, k_dim)
K <- tcrossprod(loadings) / k_dim
K <- K / mean(diag(K))              # scale so diagonal ~ 1 (coancestry-like units)
dimnames(K) <- list(pp, pp)

cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]],
                     stringsAsFactors = FALSE)
scores$poly4x_usefulness <- rnorm(nrow(scores), 12, 4)
scores$pair_kinship <- K[cbind(scores$parent1, scores$parent2)]
scores$poly4x_pair_coancestry <- scores$pair_kinship
attr(scores, "parent_kinship") <- K

n_crosses <- 8L
eps <- 1e-8

# --- 1. method = "evolution" is forwarded and runs ---
plan_evo <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                          method = "evolution", evol_iterations = 40L,
                          evol_solutions = 30L, evol_stop = 15L)
stopifnot(nrow(plan_evo) == n_crosses)
stopifnot(identical(plan_evo$poly4x_method[[1]], "ng_poly4x_ocs"))

# --- 2. committed_crosses is forwarded (locked matings appear in the plan) ---
committed <- data.frame(parent1 = c("P01", "P17"), parent2 = c("P02", "P18"),
                        stringsAsFactors = FALSE)
plan_com <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                          committed_crosses = committed, method = "greedy_local")
key_plan <- ng_group_pair_key(plan_com$parent1, plan_com$parent2)
key_com <- ng_group_pair_key(committed$parent1, committed$parent2)
stopifnot(nrow(plan_com) == n_crosses, all(key_com %in% key_plan))

# --- 3. target_coancestry is forwarded (constrained OCS lowers group coancestry) ---
base <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                      method = "greedy_local", lambda_group = 0)
base_gc <- attr(base, "summary")$group_coancestry
stopifnot(is.finite(base_gc))
tc <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                    target_coancestry = base_gc * 0.5, method = "greedy_local")
tc_sum <- attr(tc, "summary")
stopifnot(!is.null(tc_sum$target_coancestry_status))
stopifnot(tc_sum$target_coancestry_status %in% c("met", "met_slack", "infeasible"))
# A tighter coancestry target must not INCREASE realized group coancestry.
stopifnot(tc_sum$group_coancestry <= base_gc + 1e-6)

# --- 4. strategy dial is forwarded (diversity <= high_gain on group coancestry) ---
plan_gain <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                           strategy = "high_gain", method = "greedy_local")
plan_div <- ng_poly4x_ocs(scores, n_crosses = n_crosses, parent_kinship = K,
                          strategy = "diversity", method = "greedy_local")
gc_gain <- attr(plan_gain, "summary")$group_coancestry
gc_div <- attr(plan_div, "summary")$group_coancestry
stopifnot(is.finite(gc_gain), is.finite(gc_div))
stopifnot(!is.null(attr(plan_div, "summary")$achieved_emphasis))
stopifnot(gc_div <= gc_gain + 1e-6)

# --- 5. ng_poly4x_policy(mode = "ocs") forwards controls via ... ---
plan_pol <- ng_poly4x_policy(scores, n_crosses = n_crosses, mode = "ocs",
                             parent_kinship = K, method = "greedy_local",
                             committed_crosses = committed)
key_pol <- ng_group_pair_key(plan_pol$parent1, plan_pol$parent2)
stopifnot(nrow(plan_pol) == n_crosses, all(key_com %in% key_pol))
stopifnot(identical(plan_pol$poly4x_policy_mode[[1]], "ocs"))

# --- 6. ng_polyploid_policy (subgenome, R/18) forwards controls via ... ---
sub_scores <- data.frame(parent1 = scores$parent1, parent2 = scores$parent2,
                         poly_usefulness = scores$poly4x_usefulness,
                         pair_kinship = scores$pair_kinship,
                         polyploid_model_family = "allopolyploid_subgenome",
                         supported_scope = "wheat_like_allopolyploid_disomic",
                         polyploid_validation_source = "phase2a_deterministic_core",
                         stringsAsFactors = FALSE)
attr(sub_scores, "parent_kinship") <- K
plan_sub <- ng_polyploid_policy(sub_scores, n_crosses = n_crosses, mode = "ocs",
                           parent_kinship = K, method = "greedy_local",
                           committed_crosses = committed)
key_sub <- ng_group_pair_key(plan_sub$parent1, plan_sub$parent2)
stopifnot(nrow(plan_sub) == n_crosses, all(key_com %in% key_sub))
stopifnot(identical(plan_sub$poly_policy_mode[[1]], "ocs"))

cat("poly4x strategy + constraints forwarding test passed\n")
