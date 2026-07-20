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

make_small_cross_prediction_files <- function() {
  tmp <- tempfile("ng_user_cross_prediction_")
  dir.create(tmp, recursive = TRUE)

  ids <- c("P01", "P02", "P02_copy", "P03", "P04", "P05", "P06")
  geno <- data.frame(
    NAME = ids,
    M01 = c(0, 0, 0, 2, 2, 0, 2),
    M02 = c(0, 2, 2, 0, 2, 0, 2),
    M03 = c(2, 0, 0, 2, 0, 2, 0),
    M04 = c(2, 2, 2, 0, 0, 2, 0),
    M05 = c(0, 0, 0, 2, 0, 2, 2),
    M06 = c(2, 0, 0, 0, 2, 2, 0),
    M07 = c(0, 2, 2, 2, 0, 0, 2),
    M08 = c(2, 2, 2, 0, 2, 0, 0),
    check.names = FALSE
  )
  phenotype <- data.frame(
    NAME = ids,
    yield = c(58, 61, 61, 53, 67, 55, 64),
    disease = c(4.0, 3.1, 3.1, 5.2, 2.6, 4.6, 2.9),
    stringsAsFactors = FALSE
  )
  phenotype$selection_index <- as.numeric(scale(phenotype$yield)) -
    as.numeric(scale(phenotype$disease))
  marker_map <- data.frame(
    SNP_code = paste0("M", sprintf("%02d", 1:8)),
    Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
    Position_BP = c(0, 1, 2, 3, 0, 1, 2, 3) * 1e6,
    stringsAsFactors = FALSE
  )
  direction <- data.frame(
    trait = c("yield", "disease"),
    column = c("yield", "disease"),
    direction = c("increase", "decrease"),
    stringsAsFactors = FALSE
  )

  files <- list(
    dir = tmp,
    genotype = file.path(tmp, "geno.csv"),
    phenotype = file.path(tmp, "pheno.csv"),
    map = file.path(tmp, "map.csv"),
    direction = file.path(tmp, "direction.csv")
  )
  write.csv(geno, files$genotype, row.names = FALSE, quote = FALSE)
  write.csv(phenotype, files$phenotype, row.names = FALSE, quote = FALSE)
  write.csv(marker_map, files$map, row.names = FALSE, quote = FALSE)
  write.csv(direction, files$direction, row.names = FALSE, quote = FALSE)
  files
}

files <- make_small_cross_prediction_files()

trait_result <- ng_run_cross_prediction(
  phenotype_file = files$phenotype,
  genotype_file = files$genotype,
  map_file = files$map,
  direction_file = files$direction,
  id_col = "NAME",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_col = "Position_BP",
  map_pos_cm_divisor = 1e6,
  prediction_mode = "trait_by_trait",
  trait_value_metric = "usefulness",
  uc_variance_source = "pmv",
  multi_trait_method = "auto",
  trait_weights = NULL,
  progeny = "DH",
  duplicate_action = "remove",
  duplicate_threshold = 0.995,
  duplicate_maf_min = 0,
  duplicate_max_missing_prop = 1,
  duplicate_min_compared_markers = 8,
  n_crosses = 5,
  max_crosses_per_parent = 3,
  optimizer = "auto",
  use_ocs = TRUE,
  output_dir = file.path(files$dir, "trait_outputs"),
  write_outputs = FALSE,
  write_figures = TRUE,
  assume_inbred = TRUE,
  seed = 11
)

stopifnot(inherits(trait_result, "ng_cross_prediction_result"))
stopifnot(identical(trait_result$prediction_mode, "trait_by_trait"))
stopifnot(nrow(trait_result$selected_crosses) == 5L)
stopifnot(all(c("multi_trait_score", "priority_rank", "priority_tier") %in% names(trait_result$selected_crosses)))
stopifnot(all(c("yield_value", "disease_value") %in% names(trait_result$candidate_crosses)))
stopifnot(identical(trait_result$objective$diagnostics$method_reason, "equal_weight_rank_default"))
stopifnot("putative_duplicate_genotypes_removed" %in% trait_result$qc$issues$id)
stopifnot(identical(
  trait_result$qc$cleaning$putative_duplicates$removed_parents$removed_parent,
  "P02_copy"
))
stopifnot(!any(trait_result$selected_crosses$parent1 == "P02_copy" | trait_result$selected_crosses$parent2 == "P02_copy"))
stopifnot(file.exists(trait_result$output_files$priority_score_vs_kinship_png))
stopifnot(file.info(trait_result$output_files$priority_score_vs_kinship_png)$size > 1000)

index_result <- ng_run_cross_prediction(
  phenotype_file = files$phenotype,
  genotype_file = files$genotype,
  map_file = files$map,
  id_col = "NAME",
  map_marker_col = "SNP_code",
  map_chr_col = "Chromosome",
  map_pos_col = "Position_BP",
  map_pos_cm_divisor = 1e6,
  prediction_mode = "index_as_trait",
  index_col = "selection_index",
  index_direction = "increase",
  trait_value_metric = "usefulness",
  uc_variance_source = "pmv",
  progeny = "DH",
  duplicate_action = "remove",
  duplicate_threshold = 0.995,
  duplicate_maf_min = 0,
  duplicate_max_missing_prop = 1,
  duplicate_min_compared_markers = 8,
  n_crosses = 4,
  max_crosses_per_parent = 3,
  optimizer = "auto",
  use_ocs = TRUE,
  write_outputs = FALSE,
  assume_inbred = TRUE,
  seed = 11
)

stopifnot(inherits(index_result, "ng_cross_prediction_result"))
stopifnot(identical(index_result$prediction_mode, "index_as_trait"))
stopifnot(nrow(index_result$selected_crosses) == 4L)
stopifnot("selection_index_value" %in% names(index_result$candidate_crosses))
stopifnot(identical(names(index_result$marker_effects), "selection_index"))

cat("user cross prediction workflow tests passed\n")
