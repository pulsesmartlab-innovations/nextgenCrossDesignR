root_candidates <- c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  message("AlphaSimR unavailable; skipping multi-trait crop validation test")
  quit(status = 0)
}

scenario <- ng_multitrait_crop_validation_scenario(
  scenario = "compact_selfing",
  n_parents = 8L,
  n_founders = 12L,
  n_chr = 2L,
  seg_sites = 40L,
  snp_per_chr = 12L,
  qtl_per_chr = 3L,
  genome_length_m = 0.75,
  realized_progeny = 4L,
  seed = 31L,
  prediction_noise = 0.10
)

stopifnot(is.data.frame(scenario$scores))
stopifnot(nrow(scenario$scores) == 28L)
stopifnot(is.data.frame(scenario$traits))
stopifnot(nrow(scenario$traits) == 3L)
stopifnot(any(scenario$traits$direction == "maximize"))
stopifnot(any(scenario$traits$direction == "minimize"))
stopifnot(all(c("yield", "disease") %in% scenario$traits$trait))
stopifnot(all(c(
  "parent1", "parent2", "pair_kinship",
  "pred_yield", "pred_disease", "realized_yield", "realized_disease"
) %in% names(scenario$scores)))
stopifnot(all(is.finite(scenario$scores$pair_kinship)))
stopifnot(all(is.finite(scenario$scores$realized_yield)))
stopifnot(all(is.finite(scenario$scores$realized_disease)))
stopifnot(is.matrix(scenario$trait_profile$corA))
stopifnot(!isTRUE(all.equal(scenario$trait_profile$corA, diag(3L))))
stopifnot(identical(scenario$crop_scenario$scenario[[1]], "compact_selfing"))
stopifnot(is.matrix(scenario$parent_genotype))
stopifnot(!any(scenario$parent_genotype == 1, na.rm = TRUE))
stopifnot(is.matrix(scenario$parent_K))
stopifnot(identical(rownames(scenario$parent_K), scenario$parent_values$parent))

ocs_direct <- ng_run_multitrait_validation(
  scores = scenario$scores,
  traits = scenario$traits,
  realized_cols = scenario$realized_cols,
  n_crosses = 3L,
  seed = 31L,
  methods = "auto",
  allocator = "ocs",
  parent_K = scenario$parent_K,
  ocs_lambda_group = 0.20
)
stopifnot("group_coancestry" %in% names(ocs_direct$summary))
stopifnot("lambda_group" %in% names(ocs_direct$summary))
stopifnot(ocs_direct$summary$lambda_group[[1]] == 0.20)
stopifnot(is.finite(ocs_direct$summary$group_coancestry[[1]]))
manual_counts <- ng_parent_counts(ocs_direct$selections, rownames(scenario$parent_K))
manual_group <- ng_group_coancestry(manual_counts, scenario$parent_K)
stopifnot(abs(ocs_direct$summary$group_coancestry[[1]] - manual_group) < 1e-8)

methods <- c("auto", "economic_index", "desired_gain")
grid <- ng_run_multitrait_crop_validation_grid(
  scenarios = c("compact_selfing", "cassava_diploid"),
  parent_sizes = c(6L, 8L),
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 41L,
  methods = methods,
  allocator = "topn",
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L,
  prediction_noise = 0.15
)

summary <- grid$summary
selections <- grid$selections
winners <- grid$winner_summary
stopifnot(is.data.frame(summary))
stopifnot(is.data.frame(selections))
stopifnot(is.data.frame(winners))
stopifnot(nrow(summary) == 2L * 2L * length(methods))
stopifnot(nrow(selections) == 2L * 2L * length(methods) * 2L)
stopifnot(all(c("scenario", "crop", "harness_model", "n_parents", "rep", "seed") %in% names(summary)))
stopifnot(all(summary$method %in% methods))
stopifnot(all(summary$selected_crosses == 2L))
stopifnot(all(is.finite(summary$mean_realized_index)))
stopifnot(all(c("scenario", "n_parents", "metric", "method", "value", "reps") %in% names(winners)))
stopifnot(all(winners$method %in% methods))
stopifnot(all(is.finite(winners$value)))
stopifnot(identical(sort(unique(winners$scenario)), c("cassava_diploid", "compact_selfing")))
stopifnot(identical(sort(unique(winners$n_parents)), c(6L, 8L)))

cassava_alone <- ng_run_multitrait_crop_validation_grid(
  scenarios = "cassava_diploid",
  parent_sizes = 6L,
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 73L,
  methods = "auto",
  allocator = "topn",
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L,
  prediction_noise = 0.15
)
cassava_second <- ng_run_multitrait_crop_validation_grid(
  scenarios = c("compact_selfing", "cassava_diploid"),
  parent_sizes = 6L,
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 73L,
  methods = "auto",
  allocator = "topn",
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L,
  prediction_noise = 0.15
)
cassava_second_summary <- cassava_second$summary[cassava_second$summary$scenario == "cassava_diploid", , drop = FALSE]
stopifnot(all.equal(cassava_alone$summary$mean_realized_index, cassava_second_summary$mean_realized_index))
stopifnot(all.equal(cassava_alone$summary$mean_realized_yield, cassava_second_summary$mean_realized_yield))

ocs_grid <- ng_run_multitrait_crop_validation_grid(
  scenarios = "compact_selfing",
  parent_sizes = 6L,
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 79L,
  methods = c("auto", "economic_index"),
  allocator = "ocs",
  ocs_lambda_group = 0.20,
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L
)
stopifnot(all(is.finite(ocs_grid$summary$group_coancestry)))
stopifnot(all(ocs_grid$summary$allocator == "ocs"))

tmp <- tempfile("ng_multitrait_crop_validation_")
dir.create(tmp, recursive = TRUE)
written <- ng_run_multitrait_crop_validation_grid(
  scenarios = "compact_selfing",
  parent_sizes = 6L,
  reps = 1L,
  n_crosses = 2L,
  realized_progeny = 3L,
  seed = 47L,
  methods = methods,
  allocator = "topn",
  n_founders = 10L,
  n_chr = 2L,
  seg_sites = 30L,
  snp_per_chr = 8L,
  qtl_per_chr = 2L,
  output_dir = tmp,
  prefix = "crop_contract"
)
stopifnot(file.exists(file.path(tmp, "crop_contract_summary.csv")))
stopifnot(file.exists(file.path(tmp, "crop_contract_selections.csv")))
stopifnot(file.exists(file.path(tmp, "crop_contract_winner_summary.csv")))
stopifnot(file.exists(file.path(tmp, "crop_contract_config.csv")))
stopifnot(file.exists(file.path(tmp, "crop_contract_scores.csv")))
stopifnot(nrow(written$winner_summary) > 0L)
stopifnot(is.data.frame(written$scores))
stopifnot(nrow(written$scores) > nrow(written$selections))
stopifnot(all(c("methods", "n_crosses", "alphasimr_threads", "parent_generation_model") %in% names(written$config)))

cat("multi-trait crop validation tests passed\n")
