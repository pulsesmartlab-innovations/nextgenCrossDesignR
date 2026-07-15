script <- file.path("tools", "run_metric_merit_study.R")
if (!file.exists(script)) stop("Missing metric-merit study script: ", script, call. = FALSE)
txt <- paste(readLines(script, warn = FALSE), collapse = "\n")

required_tokens <- c(
  # runner conventions shared with the other tools/ runners
  "find_project_root", "run_metric_merit_study", "metric_merit_config",
  "is_this_script", "runMacs",                       # realistic LD (NOT quickHaplo)
  # the metrics being raced
  "uc_vpm", "uc_pmv", "uc_var_simple", "simplemating", "cross_mean_gebv",
  "uc_recomb_gebv", "uc_dh_gebv",
  # unbiased true-effect yardstick + measures
  "true_uc", "ng_dh_recomb_variance_pairs", "rank_acc", "topk_merit",
  "ng_selection_intensity",
  # swept conditions
  "train_n", "qtl_chr", "oligogenic", "polygenic",
  # env config surface
  "NG_MERIT_REPS", "NG_MERIT_TRAIN_N", "NG_MERIT_H2", "NG_MERIT_QTL_PER_CHR",
  "NG_MERIT_INCLUDE_SIMPLEMATING", "NG_MERIT_OUTPUT_DIR",
  # graceful optional-dependency handling + OS-agnostic parallelism
  "AlphaSimR", "requireNamespace(\"SimpleMating\"", "resolve_workers", "RGL_USE_NULL"
)
for (token in required_tokens) {
  if (!grepl(token, txt, fixed = TRUE)) {
    stop("metric-merit study script is missing required token: ", token, call. = FALSE)
  }
}

# quickHaplo produces ~0 marker-QTL LD and must NOT be the generator here.
if (grepl("quickHaplo(", txt, fixed = TRUE)) {
  stop("metric-merit study must not use quickHaplo (near-zero LD); use runMacs", call. = FALSE)
}

cat("metric-merit study contract passed\n")
