helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

scenarios <- ng_crop_genome_scenarios()
stopifnot(all(c(
  "compact_selfing",
  "maize_like",
  "large_diploid_approx",
  "bread_wheat_hexaploid_approx",
  "barley_like",
  "field_pea_like",
  "potato_tetraploid_stress",
  "cassava_diploid",
  "cassava_tetraploid_stress",
  "sugarcane_polyploid_stress"
) %in% scenarios$scenario))
stopifnot(all(c(
  "crop",
  "ploidy_label",
  "harness_model",
  "validation_scope"
) %in% names(scenarios)))
stopifnot(all(is.finite(scenarios$n_chr)))
stopifnot(all(is.finite(scenarios$genome_length_m)))
stopifnot(all(scenarios$genome_length_m > 0))
stopifnot(all(scenarios$snp_per_chr > 0))
stopifnot(all(scenarios$qtl_per_chr > 0))
stopifnot(all(nzchar(scenarios$validation_scope)))

maize <- ng_crop_genome_select("maize_like")
stopifnot(nrow(maize) == 1L)
stopifnot(maize$scenario[[1]] == "maize_like")
stopifnot(maize$n_chr[[1]] == 10L)

wheat <- ng_crop_genome_select("bread_wheat_hexaploid_approx")
stopifnot(nrow(wheat) == 1L)
stopifnot(wheat$crop[[1]] == "wheat")
stopifnot(wheat$n_chr[[1]] == 21L)
stopifnot(grepl("hexaploid", wheat$ploidy_label[[1]], fixed = TRUE))

potato <- ng_crop_genome_select("potato_tetraploid_stress")
stopifnot(nrow(potato) == 1L)
stopifnot(potato$crop[[1]] == "potato")
stopifnot(potato$n_chr[[1]] == 12L)
stopifnot(grepl("approximation", potato$harness_model[[1]], fixed = TRUE))

sugarcane <- ng_crop_genome_select("sugarcane_polyploid_stress")
stopifnot(nrow(sugarcane) == 1L)
stopifnot(sugarcane$crop[[1]] == "sugarcane")
stopifnot(grepl("aneuploid", sugarcane$ploidy_label[[1]], fixed = TRUE))

env <- ng_crop_genome_env(maize)
stopifnot(identical(env[["NG_GRID_N_CHR"]], as.character(maize$n_chr[[1]])))
stopifnot(identical(env[["NG_GRID_GENOME_LENGTH_M"]], as.character(maize$genome_length_m[[1]])))
stopifnot(identical(env[["NG_GRID_SNP_PER_CHR"]], as.character(maize$snp_per_chr[[1]])))
stopifnot(identical(env[["NG_GRID_QTL_PER_CHR"]], as.character(maize$qtl_per_chr[[1]])))
stopifnot(identical(env[["NG_GRID_PHENO_H2"]], as.character(maize$phenotype_h2[[1]])))

bad <- tryCatch(ng_crop_genome_select("unknown_crop"), error = function(e) e)
stopifnot(inherits(bad, "error"))

old <- Sys.getenv("NG_TEST_SCENARIO_ENV", unset = NA_character_)
on.exit({
  if (is.na(old)) {
    Sys.unsetenv("NG_TEST_SCENARIO_ENV")
  } else {
    Sys.setenv(NG_TEST_SCENARIO_ENV = old)
  }
}, add = TRUE)
seen <- ng_with_env_vars(c(NG_TEST_SCENARIO_ENV = "kept_name"), {
  Sys.getenv("NG_TEST_SCENARIO_ENV", unset = NA_character_)
})
stopifnot(identical(seen, "kept_name"))

cat("crop genome scenario tests passed\n")
