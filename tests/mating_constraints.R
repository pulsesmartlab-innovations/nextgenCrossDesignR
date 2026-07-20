# Breeder mating constraints: min-use-if-used, committed matings,
# group mating-permission matrix + quotas.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(202)
np <- 30L
pp <- sprintf("P%02d", seq_len(np))
cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]],
                     stringsAsFactors = FALSE)
scores$usefulness_pmv_gebv <- rnorm(nrow(scores), 10, 3)
scores$pair_kinship <- 0
K <- diag(np); dimnames(K) <- list(pp, pp)

# --- 1. min-use-if-used (via min_crosses_per_parent) ---
n_crosses <- 15L
plan <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                max_crosses_per_parent = 10L,
                                min_crosses_per_parent = 3L, method = "greedy_local")
cnt <- table(c(plan$parent1, plan$parent2))
stopifnot(nrow(plan) == n_crosses)
# every USED parent appears at least 3 times (used-or-zero property)
stopifnot(all(as.integer(cnt) >= 3L))

# --- 2. committed / fixed matings ---
committed <- data.frame(parent1 = c("P01", "P29"), parent2 = c("P02", "P30"),
                        stringsAsFactors = FALSE)
plan_c <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                  max_crosses_per_parent = 6L,
                                  committed_crosses = committed, method = "greedy_local")
key_plan <- ng_group_pair_key(plan_c$parent1, plan_c$parent2)
key_com <- ng_group_pair_key(committed$parent1, committed$parent2)
stopifnot(nrow(plan_c) == n_crosses, all(key_com %in% key_plan))

# --- 3. group mating-permission matrix (only M x F allowed) ---
grp2 <- setNames(rep(c("M", "F"), length.out = np), pp)
perm2 <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2, 2,
                dimnames = list(c("M", "F"), c("M", "F")))
plan_g <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                  max_crosses_per_parent = 6L,
                                  parent_group = grp2, group_permission = perm2,
                                  method = "greedy_local")
g1 <- grp2[plan_g$parent1]; g2 <- grp2[plan_g$parent2]
# no same-group matings survived the legality filter (assert on GROUPS, not ids)
stopifnot(nrow(plan_g) == n_crosses, all(g1 != g2))

# --- 4. group quota with a FEASIBLE setup: 3 groups, all pairings allowed, cap the
# A x B pair so the plan can still reach n_crosses via the other pairings. Assert on
# the group-pair key (grp[ids]), which is what the quota constrains.
grp3 <- setNames(rep(c("A", "B", "C"), length.out = np), pp)
perm3 <- matrix(TRUE, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
ab <- ng_group_pair_key("A", "B")
quota <- setNames(list(3L), ab)
plan_q <- ng_optimize_mating_plan(scores, n_crosses, parent_kinship = K,
                                  max_crosses_per_parent = 6L,
                                  parent_group = grp3, group_permission = perm3,
                                  group_quota = quota, method = "greedy_local")
key_q <- ng_group_pair_key(grp3[plan_q$parent1], grp3[plan_q$parent2])
stopifnot(nrow(plan_q) == n_crosses, sum(key_q == ab) <= 3L)

cat("mating constraints test passed\n")
