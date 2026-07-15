script <- file.path("tools", "run_optimizer_benchmark.R")
if (!file.exists(script)) stop("Missing optimizer bake-off script: ", script, call. = FALSE)
txt <- paste(readLines(script, warn = FALSE), collapse = "\n")

required_tokens <- c(
  "run_optimizer_benchmark",
  "build_optimizer_fixture",
  "optimizer_registry",
  "achieved_objective",
  "ng_plan_objective_contribution",
  "ng_optimize_mating_plan",
  "ng_alphamate_style_select",
  "ng_score_crosses",
  "ng_parent_kinship",
  "ng_with_rng_seed",
  # every optimizer the bake-off compares must be represented
  "mip_contribution",
  "greedy_local",
  "repair_local",
  "evolution",
  "alphamate_style",
  # outputs and provenance
  "achieved_objective",
  "group_coancestry",
  "objective_rank",
  "elapsed_sec",
  "_summary.csv",
  "_frontier.csv",
  "lambda_grid",
  "fixture_key"
)
for (token in required_tokens) {
  if (!grepl(token, txt, fixed = TRUE)) {
    stop("optimizer bake-off script is missing required token: ", token, call. = FALSE)
  }
}
cat("optimizer benchmark contract passed\n")
