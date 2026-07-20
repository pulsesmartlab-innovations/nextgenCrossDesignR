# End-to-end integration: the mate-selection modules must flow through the documented
# pipeline entry point ng_design_crosses(), not only ng_optimize_mating_plan().
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(707)
n <- 40L; m <- 300L
ids <- sprintf("P%02d", seq_len(n))
geno <- matrix(2L * rbinom(n * m, 1, 0.45), n, m,
               dimnames = list(ids, sprintf("M%03d", seq_len(m))))
beta <- rnorm(m, 0, 0.1)
y <- as.numeric(geno %*% beta) + rnorm(n)
names(y) <- ids
mm <- data.frame(marker = colnames(geno), chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 100, length.out = ceiling(m / 5)), 5)[seq_len(m)],
                 stringsAsFactors = FALSE)

nx <- 15L
base_args <- list(geno = geno, y = y, marker_map = mm, ids = ids,
                  n_crosses = nx, max_crosses_per_parent = 6L, use_cpp = FALSE)

# --- 1. metric integration: ng_score_crosses feeds expected_progeny_inbreeding through ---
d0 <- do.call(ng_design_crosses, base_args)
stopifnot("expected_progeny_inbreeding" %in% names(d0$scores))
stopifnot(nrow(d0$plan) == nx)

# --- 2. strategy dial flows through the pipeline ---
dg <- do.call(ng_design_crosses, c(base_args, list(strategy = "high_gain")))
dv <- do.call(ng_design_crosses, c(base_args, list(strategy = "diversity")))
stopifnot(identical(dg$plan_summary$strategy, "high_gain"))
stopifnot(dg$plan_summary$mean_gain >= dv$plan_summary$mean_gain - 1e-8)
stopifnot(dv$plan_summary$group_coancestry <= dg$plan_summary$group_coancestry + 1e-8)

# --- 3. the NEW evolution optimizer is reachable through the pipeline ---
de <- do.call(ng_design_crosses, c(base_args, list(method = "evolution", evol_iterations = 40L,
                                                   evol_solutions = 30L)))
stopifnot(nrow(de$plan) == nx)

# --- 4. progeny-inbreeding emphasis lowers mean progeny F end-to-end ---
dp <- do.call(ng_design_crosses, c(base_args, list(lambda_progeny_inbreeding = 50)))
stopifnot(dp$plan_summary$mean_progeny_inbreeding <=
            d0$plan_summary$mean_progeny_inbreeding + 1e-8)

# --- 5a. committed matings flow through the pipeline ---
committed <- data.frame(parent1 = "P01", parent2 = "P02", stringsAsFactors = FALSE)
dc <- do.call(ng_design_crosses, c(base_args, list(committed_crosses = committed)))
key <- ng_group_pair_key(dc$plan$parent1, dc$plan$parent2)
stopifnot(nrow(dc$plan) == nx, ng_group_pair_key("P01", "P02") %in% key)
stopifnot(identical(dc$plan_summary$n_committed, 1L))

# --- 5b. min-use-if-used flows through (every used parent >= 3, or 0) ---
dmu <- do.call(ng_design_crosses, c(base_args, list(min_crosses_per_parent = 3L)))
cnt <- table(c(dmu$plan$parent1, dmu$plan$parent2))
stopifnot(all(as.integer(cnt) >= 3L))

# --- 6. marker steering + lethal guarding flow through ---
# make M001 a steer target and M300 a lethal locus with two carriers
geno2 <- geno
geno2[, "M300"] <- 0L; geno2[c("P05", "P06"), "M300"] <- 1L
args2 <- base_args; args2$geno <- geno2
mspec <- ng_marker_target_spec("M001", direction = "increase", weight = 1)
lspec <- ng_lethal_recessive_spec("M300", risk_allele = "alt")
dm <- do.call(ng_design_crosses, c(args2, list(marker_target_spec = mspec,
                                               lethal_spec = lspec, lambda_marker = 0.5)))
stopifnot("marker_target_score" %in% names(dm$scores))
# no carrier x carrier (P05 x P06) cross in the plan
stopifnot(!any(ng_group_pair_key(dm$plan$parent1, dm$plan$parent2) ==
                 ng_group_pair_key("P05", "P06")))

# --- 7. constrained-OCS (target_coancestry) flows through ng_design_crosses ---
sweepfr <- ng_pareto_mate_allocation(d0$scores, nx, gain_col = "usefulness_pmv_gebv",
                                     parent_kinship = ng_parent_kinship(geno))$frontier
crng <- range(sweepfr$group_coancestry, na.rm = TRUE)
t_tight <- crng[1] + 0.30 * diff(crng)
t_loose <- crng[1] + 0.75 * diff(crng)
dt <- do.call(ng_design_crosses, c(base_args, list(target_coancestry = t_tight)))
dl <- do.call(ng_design_crosses, c(base_args, list(target_coancestry = t_loose)))
st <- dt$plan_summary; sl <- dl$plan_summary
stopifnot(!is.null(st$target_coancestry_status))                 # reachable end-to-end
stopifnot(st$group_coancestry <= t_tight + 1e-9)                 # constraint respected
stopifnot(sl$mean_gain >= st$mean_gain - 1e-8)                   # looser cap -> >= gain

cat("pipeline integration test passed\n")
