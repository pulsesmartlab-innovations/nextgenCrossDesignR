helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

tmp <- tempfile("ng_user_cross_prediction_posterior_params_")
dir.create(tmp, recursive = TRUE)

ids <- paste0("P", sprintf("%02d", 1:8))
geno <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 2, 2, 0, 2, 0, 2),
  M02 = c(0, 0, 2, 2, 0, 2, 0, 2),
  M03 = c(2, 0, 2, 0, 2, 0, 2, 0),
  M04 = c(2, 0, 2, 0, 2, 0, 2, 0),
  M05 = c(0, 2, 0, 2, 0, 2, 2, 0),
  M06 = c(2, 2, 0, 0, 2, 0, 0, 2),
  check.names = FALSE
)

phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 60, 53, 68, 55, 64, 59, 63),
  disease = c(4.0, 3.4, 5.1, 2.5, 4.7, 2.9, 3.6, 3.0),
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:6)),
  Chromosome = c(1, 1, 1, 1, 2, 2),
  Position_BP = c(0, 1, 3, 6, 0, 2) * 1e6,
  stringsAsFactors = FALSE
)

direction <- data.frame(
  Trait = c("yield", "disease"),
  PhenotypeColumn = c("yield", "disease"),
  Selection_direction = c("increase", "decrease"),
  stringsAsFactors = FALSE
)

files <- list(
  phenotype = file.path(tmp, "phenotype.csv"),
  genotype = file.path(tmp, "genotype.csv"),
  map = file.path(tmp, "map.csv"),
  direction = file.path(tmp, "trait_direction.csv")
)
write.csv(phenotype, files$phenotype, row.names = FALSE, quote = FALSE)
write.csv(geno, files$genotype, row.names = FALSE, quote = FALSE)
write.csv(marker_map, files$map, row.names = FALSE, quote = FALSE)
write.csv(direction, files$direction, row.names = FALSE, quote = FALSE)

result <- ng_run_cross_prediction(
  phenotype_file = files$phenotype,
  genotype_file = files$genotype,
  map_file = files$map,
  direction_file = files$direction,
  phenotype_id_col = "NAME",
  genotype_id_col = "NAME",
  direction_trait_col = "Trait",
  direction_column_col = "PhenotypeColumn",
  direction_direction_col = "Selection_direction",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,
  prediction_mode = "trait_by_trait",
  trait_value_metric = "usefulness",
  uc_variance_source = "pmv",
  method_varPMV = "full_posterior",
  ril_mode = "infinite",
  run_posterior_prediction = TRUE,
  posterior_method = "closed_form",
  n_iter = 8,
  burn_in = 3,
  use_parallel = FALSE,
  progeny = "DH",
  recomb_model = "haldane",
  selection_prop = 0.20,
  duplicate_action = "none",
  n_crosses = 4,
  max_crosses_per_parent = 3,
  optimizer = "greedy_local",
  allocation_method = "ocs",
  use_ocs = TRUE,
  write_outputs = FALSE,
  write_figures = FALSE,
  assume_inbred = TRUE,
  use_cpp = FALSE,
  seed = 20260623
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$settings$method_varPMV, "full_posterior"))
stopifnot(identical(result$settings$ril_mode, "infinite"))
stopifnot(isTRUE(result$settings$run_posterior_prediction))
stopifnot(identical(result$settings$posterior_method, "closed_form"))
stopifnot(identical(result$settings$n_iter, 8L))
stopifnot(identical(result$settings$burn_in, 3L))
stopifnot(identical(result$settings$posterior_n_draws, 5L))
stopifnot(identical(result$settings$use_parallel, FALSE))

required_cross_cols <- c(
  "yield_pmv_fast",
  "yield_pmv_full_posterior",
  "yield_pmv_used",
  "yield_value",
  "disease_pmv_fast",
  "disease_pmv_full_posterior",
  "disease_pmv_used",
  "disease_value"
)
missing_cols <- setdiff(required_cross_cols, names(result$candidate_crosses))
if (length(missing_cols)) {
  stop("candidate_crosses missing full-posterior PMV columns: ",
       paste(missing_cols, collapse = ", "))
}
stopifnot(all(is.finite(result$candidate_crosses$yield_pmv_full_posterior)))
stopifnot(all(is.finite(result$candidate_crosses$disease_pmv_full_posterior)))
stopifnot(max(abs(result$candidate_crosses$yield_pmv_used -
                    result$candidate_crosses$yield_pmv_full_posterior), na.rm = TRUE) < 1e-10)
stopifnot(max(abs(result$candidate_crosses$disease_pmv_used -
                    result$candidate_crosses$disease_pmv_full_posterior), na.rm = TRUE) < 1e-10)

i <- ng_selection_intensity(0.20)
yield_expected <- result$candidate_crosses$yield_mean +
  i * sqrt(pmax(result$candidate_crosses$yield_pmv_full_posterior, 0))
disease_expected <- result$candidate_crosses$disease_mean -
  i * sqrt(pmax(result$candidate_crosses$disease_pmv_full_posterior, 0))
stopifnot(max(abs(result$candidate_crosses$yield_value - yield_expected), na.rm = TRUE) < 1e-8)
stopifnot(max(abs(result$candidate_crosses$disease_value - disease_expected), na.rm = TRUE) < 1e-8)

pmv_delta <- max(abs(result$candidate_crosses$yield_pmv_full_posterior -
                       result$candidate_crosses$yield_pmv_fast), na.rm = TRUE)
if (!is.finite(pmv_delta) || pmv_delta <= 1e-10) {
  stop("method_varPMV = 'full_posterior' did not change PMV on the correlated-marker fixture")
}

stopifnot(is.list(result$posterior_effects))
stopifnot(is.list(result$posterior_predictions))
stopifnot(identical(sort(names(result$posterior_predictions)), c("disease", "yield")))
for (trait in names(result$posterior_predictions)) {
  posterior_scores <- result$posterior_predictions[[trait]]
  posterior_attr <- attr(posterior_scores, "posterior")
  stopifnot(identical(posterior_attr$n_draws, 5L))
  stopifnot(identical(posterior_attr$method, "closed_form"))
  stopifnot(all(c(
    "pmv_post_mean",
    "pmv_post_lower",
    "pmv_post_upper",
    "usefulness_pmv_gebv_post_mean"
  ) %in% names(posterior_scores)))
}

bad_ril_mode <- tryCatch(
  ng_run_cross_prediction(
    phenotype_file = files$phenotype,
    genotype_file = files$genotype,
    map_file = files$map,
    direction_file = files$direction,
    phenotype_id_col = "NAME",
    genotype_id_col = "NAME",
    direction_trait_col = "Trait",
    direction_column_col = "PhenotypeColumn",
    direction_direction_col = "Selection_direction",
    map_marker_col = "SNP_code",
    map_chr_col = "Chromosome",
    map_pos_bp_col = "Position_BP",
    map_position_unit = "bp",
    trait_value_metric = "pmv",
    method_varPMV = "fast",
    progeny = "RIL",
    ril_mode = "finite",
    duplicate_action = "none",
    n_crosses = 2,
    max_crosses_per_parent = 2,
    write_outputs = FALSE,
    seed = 1
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("ril_mode", bad_ril_mode, fixed = TRUE))

cat("user cross prediction posterior parameter tests passed\n")
