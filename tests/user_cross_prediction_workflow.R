# Repository root FIRST. Checking `nextgen_cross_design/` ahead of getwd() made
# these tests load a STALE 0.19.0 copy of the package that sits in the working
# tree under exactly that name -- so they validated a package eleven versions old
# while appearing to cover the current one. The ones that failed were the lucky
# case; the ones that passed gave false assurance. ng_load() already resolves in
# this order; only these hand-rolled preambles inverted it.
root_candidates <- c(
  getwd(),
  file.path(".."),
  file.path(getwd(), "nextgen_cross_design"),
  file.path("..", "nextgen_cross_design")
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

  # A seven-line panel cannot support a genomic prediction: with fewer than ten
  # genotyped-and-phenotyped records the cross-validation cannot run at all, so
  # cv_predictive_r2 is unevaluable and the reliability gate refuses the run. That
  # refusal is correct -- this fixture was demonstrating a workflow on a panel the
  # engine considers too small to draw genomic conclusions from.
  #
  # It is generated rather than hand-written so the panel can be a realistic size while
  # the file stays readable. The structure the test actually depends on is preserved
  # exactly: P02_copy is a byte-identical duplicate of P02 in BOTH genotype and
  # phenotype (so QC must detect and remove it), and yield/disease keep their opposite
  # selection directions and their original scales.
  #
  # 20 markers, not the original 8. With 8 biallelic markers there are only 256 distinct
  # genotypes, so a 60-line panel collides by chance and QC then removes several lines
  # rather than the one planted duplicate -- the assertion below is that P02_copy and
  # ONLY P02_copy is removed. duplicate_min_compared_markers = 8 is a floor and is still
  # satisfied.
  set.seed(4242L)
  n_lines <- 60L
  ids <- c(sprintf("P%02d", seq_len(n_lines)), "P02_copy")
  n_mk <- 20L
  gm <- matrix(2L * rbinom(n_lines * n_mk, 1L, 0.5), nrow = n_lines,
               dimnames = list(sprintf("P%02d", seq_len(n_lines)),
                               paste0("M", sprintf("%02d", seq_len(n_mk)))))
  stopifnot(!anyDuplicated(apply(gm, 1L, paste, collapse = "")))   # only P02_copy duplicates
  gv_y <- as.numeric(gm %*% c(1.2, -0.9, 0.7, 0, 1.0, 0, -0.6, rep(0, n_mk - 7L)))
  gv_d <- as.numeric(gm %*% c(-0.5, 0.8, 0, 0.9, 0, -1.1, 0, 0.6, rep(0, n_mk - 8L)))
  yield_v   <- 60 + 5 * as.numeric(scale(gv_y + rnorm(n_lines, 0, 0.2 * stats::sd(gv_y))))
  disease_v <- 4  + 1 * as.numeric(scale(gv_d + rnorm(n_lines, 0, 0.2 * stats::sd(gv_d))))

  dup <- 2L                                   # P02_copy duplicates P02 exactly
  geno <- data.frame(NAME = ids, rbind(gm, gm[dup, , drop = FALSE]),
                     check.names = FALSE, stringsAsFactors = FALSE)
  phenotype <- data.frame(
    NAME = ids,
    yield = c(yield_v, yield_v[[dup]]),
    disease = c(disease_v, disease_v[[dup]]),
    stringsAsFactors = FALSE
  )
  phenotype$selection_index <- as.numeric(scale(phenotype$yield)) -
    as.numeric(scale(phenotype$disease))
  marker_map <- data.frame(
    SNP_code = paste0("M", sprintf("%02d", seq_len(n_mk))),
    Chromosome = rep(1:2, each = n_mk / 2L),
    Position_BP = rep(seq_len(n_mk / 2L) - 1L, 2) * 1e6,
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
