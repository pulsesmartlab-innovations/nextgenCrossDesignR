helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ng_run_cross_prediction() now exposes ci_level (the REPORTED credible interval) and
# robustness_quantile (the tail a robust mate allocation is optimised against) as separate
# controls, and threads the ranked value's orientation into the posterior layer so a
# decrease trait gets a minimize-oriented top-N and a minimize-correct robust tail.

tmp <- tempfile("ng_runner_robust_quantile_")
dir.create(tmp, recursive = TRUE)

ids <- paste0("P", sprintf("%02d", 1:10))
set.seed(4242L)
geno <- data.frame(NAME = ids, matrix(2L * rbinom(10L * 12L, 1L, 0.5), nrow = 10L,
                                      dimnames = list(NULL, sprintf("M%02d", 1:12))),
                   check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 60, 53, 68, 55, 64, 59, 63, 61, 57),
  disease = c(4.0, 3.4, 5.1, 2.5, 4.7, 2.9, 3.6, 3.0, 3.9, 4.4),
  stringsAsFactors = FALSE
)
marker_map <- data.frame(
  SNP_code = sprintf("M%02d", 1:12),
  Chromosome = rep(1:2, each = 6L),
  Position_BP = rep(c(0, 2, 5, 9, 14, 20), 2) * 1e6,
  stringsAsFactors = FALSE
)
direction <- data.frame(
  Trait = c("yield", "disease"),
  PhenotypeColumn = c("yield", "disease"),
  Selection_direction = c("increase", "decrease"),
  stringsAsFactors = FALSE
)
files <- list(phenotype = file.path(tmp, "phenotype.csv"),
              genotype = file.path(tmp, "genotype.csv"),
              map = file.path(tmp, "map.csv"),
              direction = file.path(tmp, "trait_direction.csv"))
write.csv(phenotype, files$phenotype, row.names = FALSE, quote = FALSE)
write.csv(geno, files$genotype, row.names = FALSE, quote = FALSE)
write.csv(marker_map, files$map, row.names = FALSE, quote = FALSE)
write.csv(direction, files$direction, row.names = FALSE, quote = FALSE)

run <- function(...) {
  ng_run_cross_prediction(
    phenotype_file = files$phenotype, genotype_file = files$genotype,
    map_file = files$map, direction_file = files$direction,
    phenotype_id_col = "NAME", genotype_id_col = "NAME",
    direction_trait_col = "Trait", direction_column_col = "PhenotypeColumn",
    direction_direction_col = "Selection_direction",
    map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_bp_col = "Position_BP", map_position_unit = "bp", bp_per_cm = 1e6,
    prediction_mode = "trait_by_trait", trait_value_metric = "usefulness",
    uc_variance_source = "pmv", progeny = "DH", recomb_model = "haldane",
    selection_prop = 0.20, duplicate_action = "none",
    run_posterior_prediction = TRUE, posterior_method = "closed_form",
    n_iter = 160L, burn_in = 10L, n_crosses = 5, max_crosses_per_parent = 3,
    optimizer = "greedy_local", write_outputs = FALSE, write_figures = FALSE,
    parent_type = "inbred", use_cpp = FALSE, seed = 909L, ...
  )
}

rq <- 0.25
res_plain <- run()
res_rq <- run(robustness_quantile = rq)

stopifnot(isTRUE(all.equal(res_plain$settings$ci_level, 0.95)))
stopifnot(is.na(res_plain$settings$robustness_quantile))
stopifnot(isTRUE(all.equal(res_rq$settings$robustness_quantile, rq)))

gain_col <- "usefulness_pmv_gebv"
q_lo_col <- ng_posterior_quantile_col(gain_col, rq)
q_hi_col <- ng_posterior_quantile_col(gain_col, 1 - rq)
expected_direction <- c(yield = "maximize", disease = "minimize")

for (trait in c("yield", "disease")) {
  plain <- res_plain$posterior_predictions[[trait]]
  withq <- res_rq$posterior_predictions[[trait]]
  meta_plain <- attr(plain, "posterior")
  meta_q <- attr(withq, "posterior")
  # Orientation of the ranked value reaches the posterior layer.
  stopifnot(identical(meta_plain$direction, unname(expected_direction[[trait]])))
  stopifnot(identical(meta_q$direction, unname(expected_direction[[trait]])))
  # Both robust tails are cached; the reported interval is untouched by the request.
  stopifnot(all(c(q_lo_col, q_hi_col) %in% names(withq)))
  stopifnot(!any(c(q_lo_col, q_hi_col) %in% names(plain)))
  stopifnot(identical(meta_q$ci_level, 0.95))
  for (col in c(paste0(gain_col, c("_post_mean", "_post_lower", "_post_upper")),
                "pmv_post_lower", "pmv_post_upper")) {
    if (!isTRUE(all.equal(plain[[col]], withq[[col]], tolerance = 0))) {
      stop("robustness_quantile request moved the reported column ", col, " for trait ", trait)
    }
  }
}

# ci_level is an explicit, honest control over the REPORTED interval: widening it moves the
# interval and nothing else.
res_ci80 <- run(ci_level = 0.80)
stopifnot(isTRUE(all.equal(res_ci80$settings$ci_level, 0.80)))
p95 <- res_plain$posterior_predictions$yield
p80 <- res_ci80$posterior_predictions$yield
stopifnot(identical(attr(p80, "posterior")$ci_level, 0.80))
stopifnot(all(p80[[paste0(gain_col, "_post_lower")]] >=
                p95[[paste0(gain_col, "_post_lower")]] - 1e-12))
stopifnot(all(p80[[paste0(gain_col, "_post_upper")]] <=
                p95[[paste0(gain_col, "_post_upper")]] + 1e-12))
stopifnot(isTRUE(all.equal(p80[[paste0(gain_col, "_post_mean")]],
                           p95[[paste0(gain_col, "_post_mean")]], tolerance = 0)))

# End-to-end: the frontend's orchestration (run, then robust-allocate at the breeder's
# quantile) now succeeds for BOTH directions without allow_normal_approximation.
geno_mat <- as.matrix(geno[, -1, drop = FALSE])
rownames(geno_mat) <- geno$NAME
kin <- ng_parent_kinship(geno_mat)
for (trait in c("yield", "disease")) {
  scores <- res_rq$posterior_predictions[[trait]]
  dir_i <- attr(scores, "posterior")$direction
  plan <- ng_optimize_robust_mating_plan(
    posterior_scores = scores, n_crosses = 5L, parent_kinship = kin,
    gain_col = gain_col, robustness_quantile = rq, direction = dir_i,
    max_crosses_per_parent = 3L, method = "greedy_local")
  s <- attr(plan, "summary")
  stopifnot(nrow(plan) == 5L)
  stopifnot(isFALSE(s$robustness_quantile_is_normal_approximation))
  stopifnot(identical(s$robust_direction, dir_i))
  stopifnot(identical(s$robust_quantile_source,
                      if (identical(dir_i, "maximize")) q_lo_col else q_hi_col))
}

bad_rq <- tryCatch(run(robustness_quantile = 1.5), error = function(e) conditionMessage(e))
stopifnot(grepl("robustness_quantile must be one finite probability", bad_rq, fixed = TRUE))
bad_ci <- tryCatch(run(ci_level = 0), error = function(e) conditionMessage(e))
stopifnot(grepl("ci_level must be one finite probability", bad_ci, fixed = TRUE))

cat("runner_robust_quantile: all checks passed\n")
