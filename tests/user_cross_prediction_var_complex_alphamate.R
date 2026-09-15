installed_lib <- Sys.getenv("NG_TEST_INSTALLED_LIB", unset = "")
if (nzchar(installed_lib)) {
  .libPaths(c(
    normalizePath(installed_lib, winslash = "/", mustWork = TRUE),
    if (.Platform$OS.type == "windows") normalizePath(file.path(getwd(), "..", ".Rlib"), winslash = "/", mustWork = FALSE) else NULL,
    .libPaths()
  ))
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

tmp <- tempfile("ng_var_complex_alphamate_")
dir.create(tmp, recursive = TRUE)

# Seven lines cannot be cross-validated, so the reliability gate refuses the run before
# the var_complex + AlphaMate allocation path this file checks can execute. Generated
# rather than hand-written so the panel is realistic while the file stays short; the
# assertions here are about the allocation wrapper and metric plumbing, not about the
# particular numbers.
set.seed(8675L)
n_lines <- 60L; n_mk <- 20L
ids <- sprintf("P%02d", seq_len(n_lines))
gm <- matrix(2L * rbinom(n_lines * n_mk, 1L, 0.5), nrow = n_lines,
             dimnames = list(ids, paste0("M", sprintf("%02d", seq_len(n_mk)))))
geno <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)

gv_y <- as.numeric(gm %*% c(1.2, -0.9, 0.7, 0, 1.0, -0.6, rep(0, n_mk - 6L)))
gv_d <- as.numeric(gm %*% c(-0.5, 0.8, 0, 0.9, -1.1, 0.6, rep(0, n_mk - 6L)))
phenotype <- data.frame(
  NAME = ids,
  yield   = 60 + 5 * as.numeric(scale(gv_y + rnorm(n_lines, 0, 0.2 * stats::sd(gv_y)))),
  disease = 4  + 1 * as.numeric(scale(gv_d + rnorm(n_lines, 0, 0.2 * stats::sd(gv_d)))),
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  SNP_code = colnames(gm),
  Chromosome = rep(1:2, each = n_mk / 2L),
  Position_BP = rep(seq_len(n_mk / 2L) - 1L, 2) * 1e6,
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
  max_crosses_per_parent = 3,
  allocation_method = "alphamate_style",
  alphamate_mode = "ModeOptTarget1",
  alphamate_target_degree = 60,
  alphamate_max_contributions = 3,
  optimizer = "greedy_local",
  use_ocs = TRUE,
  seed = 17
)

stopifnot(inherits(result, "ng_cross_prediction_result"))
# `var_complex` is a REQUEST, not a resolved method: ng_cp__build_ctx() maps it to
# trait_value_metric = "usefulness" + uc_variance_source = "pmv" (unchanged since at
# least 0.30.0). settings$trait_value_metric therefore carries the resolved token and
# settings$trait_value_metric_input carries what the caller actually typed -- which is
# the pair this assertion should check, not one of them alone.
#
# The old assertion expected "var_complex" back from settings and passed only because
# this file was loading a stale 0.19.0 copy of the package, which predates the
# normalisation. Against the real package it was always wrong.
stopifnot(identical(result$settings$trait_value_metric_input, "var_complex"))
stopifnot(identical(result$settings$trait_value_metric, "usefulness"))
stopifnot(identical(result$effect_summary$trait_value_metric_resolved[[1L]], "usefulness"))
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
    bp_per_cm = 1e6,
    parent_type = "inbred",
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
