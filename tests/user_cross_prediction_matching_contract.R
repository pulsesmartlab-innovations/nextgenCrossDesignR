installed_lib <- Sys.getenv("NG_TEST_INSTALLED_LIB", unset = "")
if (nzchar(installed_lib)) {
  dep_lib <- if (.Platform$OS.type == "windows") normalizePath(file.path("..", ".Rlib"), winslash = "/", mustWork = FALSE) else NULL
  .libPaths(c(normalizePath(installed_lib, winslash = "/", mustWork = TRUE), dep_lib, .libPaths()))
  library(nextgenCrossDesign)
} else {
  # This preamble used to try getwd()/nextgen_cross_design BEFORE getwd(), so an
  # untracked stale copy of the package beside the real sources won. In this working
  # tree one exists at 0.19.0 -- eleven releases back, with no reliability gate -- so
  # this file passed locally while testing an obsolete package, and failed on CI where
  # no such copy exists. A green local run was not evidence of anything here.
  #
  # Commit ac35088 fixed this for eight harness tests; this file and
  # user_cross_prediction_var_complex_alphamate.R were missed. Both now use the one
  # shared resolver, which prefers a directory that is ITSELF the package over any copy
  # nested inside it and requires DESCRIPTION to name this package.
  local({
    cands <- file.path(c(".", "..", "../..", "nextgen_cross_design",
                         "../nextgen_cross_design"), "tools", "ng_find_package_root.R")
    hit <- cands[file.exists(cands)]
    if (!length(hit)) stop("cannot locate tools/ng_find_package_root.R", call. = FALSE)
    source(hit[[1L]], local = FALSE)
  })
  root <- normalizePath(ng_find_package_root(getwd()), mustWork = TRUE)
  source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
  source(file.path(root, "R", "load.R"))
  ng_load(root, use_cpp = FALSE, verbose = FALSE)
}

tmp <- tempfile("ng_matching_contract_")
dir.create(tmp, recursive = TRUE)

# Six lines cannot be cross-validated: below ten genotyped-and-phenotyped records
# cv_predictive_r2 is unevaluable, so the reliability gate refuses the run and the
# ID-matching contract this file exists to check never executes. Generated rather than
# hand-written so the panel is a realistic size while the file stays short.
#
# The structural point is preserved and is the whole reason for the fixture: the
# phenotype rows are in a DIFFERENT order from the genotype rows (and use different ID
# and trait column names), so the package must match by ID and reorder internally
# before fitting. The shuffle below is explicit, not incidental.
set.seed(1234L)
n_lines <- 60L; n_mk <- 20L
geno_ids <- sprintf("G%02d", seq_len(n_lines))
gm <- matrix(2L * rbinom(n_lines * n_mk, 1L, 0.5), nrow = n_lines,
             dimnames = list(geno_ids, paste0("M", sprintf("%02d", seq_len(n_mk)))))
geno <- data.frame(TaxaID = geno_ids, gm, check.names = FALSE, stringsAsFactors = FALSE)

gv_y <- as.numeric(gm %*% c(1.2, -0.9, 0.7, 0, 1.0, -0.6, rep(0, n_mk - 6L)))
gv_d <- as.numeric(gm %*% c(-0.5, 0.8, 0, 0.9, -1.1, 0.6, rep(0, n_mk - 6L)))
yield_v   <- 60 + 5 * as.numeric(scale(gv_y + rnorm(n_lines, 0, 0.2 * stats::sd(gv_y))))
disease_v <- 4  + 1 * as.numeric(scale(gv_d + rnorm(n_lines, 0, 0.2 * stats::sd(gv_d))))

# Phenotype rows are deliberately not in genotype order. The package should
# match by ID and reorder internally before fitting marker effects.
ord <- sample.int(n_lines)
stopifnot(!identical(ord, seq_len(n_lines)))        # the shuffle must actually shuffle
pheno <- data.frame(
  LineID = geno_ids[ord],
  Yield_t_ha = yield_v[ord],
  Disease_score = disease_v[ord],
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  MarkerName = colnames(gm),
  Chr = rep(1:2, each = n_mk / 2L),
  Position_BP = rep(seq_len(n_mk / 2L) - 1L, 2) * 5e5,
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
  parent_type = "inbred",
  seed = 13
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
stopifnot(identical(rownames(result$cleaned_data$genotype), geno_ids))
stopifnot(identical(result$input_match_audit$phenotype_id_col, "LineID"))
stopifnot(identical(result$input_match_audit$genotype_id_col, "TaxaID"))
stopifnot(identical(result$input_match_audit$matched_parent_count, n_lines))
stopifnot(identical(result$input_match_audit$marker_count, n_mk))
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
    bp_per_cm = 1e6,
    parent_type = "inbred",
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
    bp_per_cm = 1e6,
    parent_type = "inbred",
    duplicate_action = "none",
    n_crosses = 2
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("phenotype", missing_trait_error, ignore.case = TRUE))
stopifnot(grepl("not_a_trait_column", missing_trait_error, fixed = TRUE))

cat("user cross prediction matching contract tests passed\n")
