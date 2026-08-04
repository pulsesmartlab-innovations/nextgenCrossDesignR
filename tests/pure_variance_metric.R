helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

st <- data.frame(
  cross_mean_blend = c(10, 2),   # cross A high mean, low var; cross B low mean, high var
  vpm = c(1, 9), pmv = c(1, 9), parent_distance = c(0, 0)
)
dir <- "increase"
# pure variance: cross B (var 9) outranks cross A (var 1), regardless of mean
pv <- ng_run_cp_trait_value(st, dir, trait_value_metric = "vpm", selection_prop = 0.1)
stopifnot(pv[2] > pv[1])                 # ranked by variance alone
stopifnot(abs(pv[1] - 1) < 1e-9)         # returns the variance value itself (vpm=1)
stopifnot(abs(pv[2] - 9) < 1e-9)
# usefulness with the same variance: mean dominates -> cross A wins
uc <- ng_run_cp_trait_value(st, dir, trait_value_metric = "usefulness",
                            uc_variance_source = "vpm", selection_prop = 0.1)
stopifnot(uc[1] > uc[2])
# pmv pure variance also returns the pmv column
pv2 <- ng_run_cp_trait_value(st, dir, trait_value_metric = "pmv", selection_prop = 0.1, method_varPMV = "fast")
stopifnot(abs(pv2[2] - 9) < 1e-9)
cat("pure-variance metric OK\n")
