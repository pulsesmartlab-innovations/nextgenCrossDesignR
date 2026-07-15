helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- ng_test_find_root()
source(file.path(root, "tools", "run_optimizer_benchmark.R"))

out_dir <- file.path(tempdir(), paste0("optbench_", as.integer(Sys.time())))
cfg <- optimizer_config_from_env(root)
cfg$n_parents <- 24L
cfg$n_markers <- 250L
cfg$training_n <- 80L
cfg$n_crosses <- 12L
cfg$lambda_grid <- c(0.1, 1, 10)
cfg$use_cpp <- isTRUE(get0("ng_test_use_cpp", ifnotfound = FALSE, inherits = FALSE))
cfg$output_dir <- out_dir
cfg$output_prefix <- "optbench_smoke"

res <- run_optimizer_benchmark(cfg)
s <- res$summary

stopifnot(is.data.frame(s), nrow(s) >= 1L)
stopifnot(all(c("mip_contribution", "greedy_local", "repair_local") %in% s$method))

ok <- s[s$status == "ok", , drop = FALSE]
stopifnot(nrow(ok) >= 3L)
# every successful optimizer must return exactly n_crosses and a finite achieved objective
stopifnot(all(ok$n_selected == cfg$n_crosses))
stopifnot(all(is.finite(ok$achieved_objective)))
stopifnot(all(ok$max_parent_use <= cfg$max_crosses_per_parent, na.rm = TRUE))

# the evolutionary optimizer, once available, must be feasible and no worse than the
# greedy warm start it is seeded from (elitism guarantees this); skip if not yet built.
ev <- s[s$method == "evolution", , drop = FALSE]
gr <- s[s$method == "greedy_local", , drop = FALSE]
if (nrow(ev) && identical(ev$status[[1]], "ok") && nrow(gr) && identical(gr$status[[1]], "ok")) {
  stopifnot(ev$n_selected[[1]] == cfg$n_crosses)
  stopifnot(ev$achieved_objective[[1]] >= gr$achieved_objective[[1]] - 1e-6)
  cat("evolution feasible and >= greedy_local achieved objective\n")
} else {
  cat("evolution optimizer not available yet; skipped its assertions\n")
}

# outputs were written
stopifnot(file.exists(file.path(out_dir, "optbench_smoke_summary.csv")))

cat("optimizer benchmark smoke passed\n")
