installed_lib <- Sys.getenv("NG_TEST_INSTALLED_LIB", unset = "")
if (nzchar(installed_lib)) {
  dep_lib <- if (.Platform$OS.type == "windows") normalizePath(file.path("..", ".Rlib"), winslash = "/", mustWork = FALSE) else NULL
  .libPaths(c(normalizePath(installed_lib, winslash = "/", mustWork = TRUE), dep_lib, .libPaths()))
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

tmp <- tempfile("ng_matching_contract_")
dir.create(tmp, recursive = TRUE)

geno_ids <- c("G01", "G02", "G03", "G04", "G05", "G06")
geno <- data.frame(
  TaxaID = geno_ids,
  M01 = c(0, 0, 2, 2, 0, 2),
  M02 = c(0, 2, 0, 2, 0, 2),
  M03 = c(2, 0, 2, 0, 2, 0),
  M04 = c(2, 2, 0, 0, 2, 0),
  M05 = c(0, 0, 2, 0, 2, 2),
  M06 = c(2, 0, 0, 2, 2, 0),
  check.names = FALSE
)

# Phenotype rows are deliberately not in genotype order. The package should
# match by ID and reorder internally before fitting marker effects.
pheno <- data.frame(
  LineID = c("G04", "G02", "G06", "G01", "G05", "G03"),
  Yield_t_ha = c(67, 61, 64, 58, 55, 53),
  Disease_score = c(2.6, 3.1, 2.9, 4.0, 4.6, 5.2),
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  MarkerName = paste0("M", sprintf("%02d", 1:6)),
  Chr = c(1, 1, 1, 2, 2, 2),
  Position_BP = c(0, 5e5, 1e6, 0, 5e5, 1e6),
  stringsAsFactors = FALSE
)

direction <- data.frame(
  TraitName = c("yield", "disease"),
  PhenotypeColumn = c("Yield_t_ha", "Disease_score"),
  SelectionDirection = c("increase", "decrease"),
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(tmp, "phenotype.csv")
genotype_file <- file.path(tmp, "genotype.csv")
map_file <- file.path(tmp, "map.csv")
direction_file <- file.path(tmp, "direction.csv")
write.csv(pheno, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(geno, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(direction, direction_file, row.names = FALSE, quote = FALSE)

result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,
  phenotype_id_col = "LineID",
  genotype_id_col = "TaxaID",
  direction_trait_col = "TraitName",
  direction_column_col = "PhenotypeColumn",
  direction_direction_col = "SelectionDirection",
  map_marker_col = "MarkerName",
  map_chr_col = "Chr",
  map_pos_bp_col = "Position_BP",
  map_position_unit = "bp",
  bp_per_cm = 1e6,
  prediction_mode = "trait_by_trait",
  trait_value_metric = "usefulness",
  uc_variance_source = "pmv",
  multi_trait_method = "auto",
  progeny = "DH",
  duplicate_action = "none",
  n_crosses = 4,
  max_crosses_per_parent = 3,
  optimizer = "auto",
  use_ocs = TRUE,
  assume_inbred = TRUE,
  seed = 13
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
stopifnot(identical(rownames(result$cleaned_data$genotype), geno_ids))
stopifnot(identical(result$input_match_audit$phenotype_id_col, "LineID"))
stopifnot(identical(result$input_match_audit$genotype_id_col, "TaxaID"))
stopifnot(identical(result$input_match_audit$matched_parent_count, 6L))
stopifnot(identical(result$input_match_audit$marker_count, 6L))
stopifnot(identical(result$input_match_audit$map_position_unit, "bp"))
stopifnot(identical(result$input_match_audit$direction_columns$trait, "TraitName"))
stopifnot(identical(result$input_match_audit$direction_columns$column, "PhenotypeColumn"))
stopifnot(identical(result$input_match_audit$direction_columns$direction, "SelectionDirection"))
stopifnot(all(c("yield_value", "disease_value") %in% names(result$candidate_crosses)))
stopifnot(all(c("pos_bp", "pos_cm") %in% names(result$cleaned_data$marker_map)))
stopifnot(identical(result$cleaned_data$marker_map$marker, colnames(result$cleaned_data$genotype)))

bad_map <- marker_map[-6, , drop = FALSE]
write.csv(bad_map, map_file, row.names = FALSE, quote = FALSE)
missing_marker_error <- tryCatch(
  ng_run_cross_prediction(
    phenotype_file = phenotype_file,
    genotype_file = genotype_file,
    map_file = map_file,
    direction_file = direction_file,
    phenotype_id_col = "LineID",
    genotype_id_col = "TaxaID",
    direction_trait_col = "TraitName",
    direction_column_col = "PhenotypeColumn",
    direction_direction_col = "SelectionDirection",
    map_marker_col = "MarkerName",
    map_chr_col = "Chr",
    map_pos_bp_col = "Position_BP",
    map_position_unit = "bp",
    duplicate_action = "none",
    n_crosses = 2
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("marker map", missing_marker_error, ignore.case = TRUE))
stopifnot(grepl("missing", missing_marker_error, ignore.case = TRUE))

write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
bad_direction <- direction
bad_direction$PhenotypeColumn[[1]] <- "not_a_trait_column"
write.csv(bad_direction, direction_file, row.names = FALSE, quote = FALSE)
missing_trait_error <- tryCatch(
  ng_run_cross_prediction(
    phenotype_file = phenotype_file,
    genotype_file = genotype_file,
    map_file = map_file,
    direction_file = direction_file,
    phenotype_id_col = "LineID",
    genotype_id_col = "TaxaID",
    direction_trait_col = "TraitName",
    direction_column_col = "PhenotypeColumn",
    direction_direction_col = "SelectionDirection",
    map_marker_col = "MarkerName",
    map_chr_col = "Chr",
    map_pos_bp_col = "Position_BP",
    map_position_unit = "bp",
    duplicate_action = "none",
    n_crosses = 2
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("phenotype", missing_trait_error, ignore.case = TRUE))
stopifnot(grepl("not_a_trait_column", missing_trait_error, fixed = TRUE))

cat("user cross prediction matching contract tests passed\n")
