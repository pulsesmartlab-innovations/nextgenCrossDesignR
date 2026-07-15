installed_lib <- Sys.getenv("NG_TEST_INSTALLED_LIB", unset = "")
if (nzchar(installed_lib)) {
  .libPaths(c(
    normalizePath(installed_lib, winslash = "/", mustWork = TRUE),
    if (.Platform$OS.type == "windows") normalizePath(file.path(getwd(), "..", ".Rlib"), winslash = "/", mustWork = FALSE) else NULL,
    .libPaths()
  ))
  library(nextgenCrossDesign)
} else {
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
}

tmp <- tempfile("ng_var_complex_alphamate_")
dir.create(tmp, recursive = TRUE)

ids <- paste0("P", sprintf("%02d", 1:7))
geno <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 2, 2, 0, 2, 0),
  M02 = c(0, 2, 0, 2, 0, 2, 2),
  M03 = c(2, 0, 2, 0, 2, 0, 2),
  M04 = c(2, 2, 0, 0, 2, 0, 0),
  M05 = c(0, 0, 2, 0, 2, 2, 2),
  M06 = c(2, 0, 0, 2, 2, 0, 0),
  M07 = c(0, 2, 2, 2, 0, 0, 2),
  M08 = c(2, 2, 2, 0, 2, 0, 0),
  check.names = FALSE
)
phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 61, 53, 67, 55, 64, 60),
  disease = c(4.0, 3.1, 5.2, 2.6, 4.6, 2.9, 3.4),
  stringsAsFactors = FALSE
)
marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:8)),
  Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
  Position_BP = c(0, 1, 2, 3, 0, 1, 2, 3) * 1e6,
  stringsAsFactors = FALSE
)
direction <- data.frame(
  Trait = c("yield", "disease"),
  Selection_direction = c("increase", "decrease"),
  stringsAsFactors = FALSE
)

files <- list(
  genotype = file.path(tmp, "geno.csv"),
  phenotype = file.path(tmp, "pheno.csv"),
  map = file.path(tmp, "map.csv"),
  direction = file.path(tmp, "direction.csv")
)
write.csv(geno, files$genotype, row.names = FALSE, quote = FALSE)
write.csv(phenotype, files$phenotype, row.names = FALSE, quote = FALSE)
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
  direction_direction_col = "Selection_direction",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,
  trait_value_metric = "var_complex",
  selection_prop = 0.20,
  progeny = "DH",
  duplicate_action = "none",
  n_crosses = 4,
  max_uses_per_parent = 3,
  allocation_method = "alphamate_style",
  alphamate_mode = "ModeOptTarget1",
  alphamate_target_degree = 60,
  alphamate_max_contributions = 3,
  optimizer = "greedy_local",
  use_ocs = TRUE,
  seed = 17
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(result$settings$trait_value_metric, "var_complex"))
stopifnot(identical(result$settings$allocation_method, "alphamate_style"))
stopifnot(nrow(result$selected_crosses) == 4L)
stopifnot(identical(result$plan_summary$alphamate_style, "native_proxy"))
stopifnot(identical(result$plan_summary$alphamate_mode, "ModeOptTarget1"))
stopifnot(identical(result$plan_summary$alphamate_target_degree, 60))
stopifnot(all(c("yield_var_complex", "disease_var_complex") %in% names(result$candidate_crosses)))

i <- {
  z <- stats::qnorm(1 - 0.20)
  stats::dnorm(z) / 0.20
}
yield_expected <- result$candidate_crosses$yield_mean + i * sqrt(pmax(result$candidate_crosses$yield_var_complex, 0))
disease_expected <- result$candidate_crosses$disease_mean - i * sqrt(pmax(result$candidate_crosses$disease_var_complex, 0))
stopifnot(max(abs(result$candidate_crosses$yield_value - yield_expected), na.rm = TRUE) < 1e-8)
stopifnot(max(abs(result$candidate_crosses$disease_value - disease_expected), na.rm = TRUE) < 1e-8)

missing_exe <- tryCatch(
  ng_run_cross_prediction(
    phenotype_file = files$phenotype,
    genotype_file = files$genotype,
    map_file = files$map,
    direction_file = files$direction,
    phenotype_id_col = "NAME",
    genotype_id_col = "NAME",
    direction_trait_col = "Trait",
    direction_direction_col = "Selection_direction",
    map_marker_col = "SNP_code",
    map_chr_col = "Chromosome",
    map_pos_bp_col = "Position_BP",
    map_position_unit = "bp",
    trait_value_metric = "var_complex",
    duplicate_action = "none",
    n_crosses = 2,
    allocation_method = "alphamate_executable",
    alphamate_executable = file.path(tmp, "missing_AlphaMate.exe"),
    seed = 17
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("AlphaMate executable not found", missing_exe, fixed = TRUE))

cat("var_complex and AlphaMate allocation wrapper tests passed\n")
