ng_crop_genome_scenarios <- function() {
  data.frame(
    scenario = c(
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
    ),
    crop = c(
      "generic_selfing",
      "maize",
      "large_diploid",
      "wheat",
      "barley",
      "field_pea",
      "potato",
      "cassava",
      "cassava",
      "sugarcane"
    ),
    description = c(
      "Compact diploid selfing-crop approximation with moderate heritability.",
      "Medium diploid maize-like approximation with longer genetic map.",
      "Large diploid approximation for marker-dense cereal stress testing.",
      "Bread-wheat-like 21-linkage-group stress test; hexaploidy represented as diploidized chromosomes.",
      "Barley-like diploid cereal approximation with seven long linkage groups.",
      "Field-pea-like diploid pulse approximation with seven large linkage groups.",
      "Cultivated-potato-like tetraploid stress test represented by 12 diploidized homology groups.",
      "Cassava-like diploid clonal-crop approximation with 18 linkage groups.",
      "Cassava autotetraploid stress test represented as diploidized 18-linkage-group dosage pressure.",
      "Modern-sugarcane-like polyploid/aneuploid stress test represented as diploidized homology groups."
    ),
    ploidy_label = c(
      "diploid",
      "diploid",
      "diploid",
      "allohexaploid 2n=6x=42",
      "diploid",
      "diploid 2n=14",
      "autotetraploid 2n=4x=48",
      "diploid 2n=36",
      "autotetraploid 2n=4x=72",
      "complex polyploid/aneuploid 2n~100-120"
    ),
    harness_model = c(
      "diploid DH/RIL",
      "diploid DH/RIL",
      "diploid DH/RIL",
      "diploidized approximation",
      "diploid DH/RIL",
      "diploid DH/RIL",
      "diploidized approximation",
      "diploid DH/RIL",
      "diploidized approximation",
      "diploidized approximation"
    ),
    validation_scope = c(
      "Directly supported by the current diploid DH/RIL harness.",
      "Directly supported by the current diploid DH/RIL harness.",
      "Directly supported by the current diploid DH/RIL harness.",
      "Stress test for chromosome count and map length; not true hexaploid inheritance.",
      "Directly supported by the current diploid DH/RIL harness.",
      "Directly supported by the current diploid DH/RIL harness.",
      "Stress test for potato-like scale and diversity; not true autotetraploid dosage inheritance.",
      "Directly supported by the current diploid DH/RIL harness, but clonal breeding is approximated.",
      "Stress test for cassava tetraploid scale; not true autotetraploid dosage inheritance.",
      "Stress test for sugarcane-like complexity; not true aneuploid polysomic inheritance."
    ),
    n_chr = c(7L, 10L, 21L, 21L, 7L, 7L, 12L, 18L, 18L, 10L),
    genome_length_m = c(0.75, 1.60, 1.20, 2.35, 1.85, 1.80, 0.95, 1.20, 1.20, 1.65),
    seg_sites = c(900L, 1500L, 1800L, 1200L, 1500L, 1200L, 1100L, 1100L, 1200L, 1800L),
    snp_per_chr = c(650L, 700L, 350L, 300L, 700L, 500L, 550L, 350L, 450L, 800L),
    qtl_per_chr = c(24L, 36L, 22L, 18L, 35L, 30L, 35L, 20L, 25L, 40L),
    phenotype_h2 = c(0.55, 0.45, 0.35, 0.40, 0.45, 0.40, 0.35, 0.40, 0.35, 0.30),
    effect_training_n = c(300L, 400L, 500L, 500L, 400L, 350L, 400L, 400L, 450L, 500L),
    n_founders = c(120L, 140L, 160L, 180L, 140L, 120L, 160L, 140L, 160L, 180L),
    stringsAsFactors = FALSE
  )
}

ng_crop_genome_select <- function(scenarios = NULL) {
  available <- ng_crop_genome_scenarios()
  if (is.null(scenarios) || !length(scenarios) ||
      (length(scenarios) == 1L && !nzchar(trimws(as.character(scenarios))))) {
    return(available)
  }
  requested <- unlist(strsplit(paste(as.character(scenarios), collapse = ","), ",", fixed = TRUE),
                      use.names = FALSE)
  requested <- unique(trimws(requested[nzchar(trimws(requested))]))
  missing <- setdiff(requested, available$scenario)
  if (length(missing)) {
    ng_stop("Unknown crop genome scenario(s): ", paste(missing, collapse = ", "))
  }
  available[match(requested, available$scenario), , drop = FALSE]
}

ng_crop_genome_env <- function(scenario) {
  scenario <- as.data.frame(scenario, stringsAsFactors = FALSE)
  if (nrow(scenario) != 1L) ng_stop("scenario must contain exactly one row")
  c(
    NG_CROP_SCENARIO = as.character(scenario$scenario[[1]]),
    NG_CROP = as.character(scenario$crop[[1]]),
    NG_CROP_HARNESS_MODEL = as.character(scenario$harness_model[[1]]),
    NG_CROP_VALIDATION_SCOPE = as.character(scenario$validation_scope[[1]]),
    NG_GRID_N_FOUNDERS = as.character(as.integer(scenario$n_founders[[1]])),
    NG_GRID_N_CHR = as.character(as.integer(scenario$n_chr[[1]])),
    NG_GRID_GENOME_LENGTH_M = as.character(as.numeric(scenario$genome_length_m[[1]])),
    NG_GRID_SEG_SITES = as.character(as.integer(scenario$seg_sites[[1]])),
    NG_GRID_SNP_PER_CHR = as.character(as.integer(scenario$snp_per_chr[[1]])),
    NG_GRID_QTL_PER_CHR = as.character(as.integer(scenario$qtl_per_chr[[1]])),
    NG_GRID_PHENO_H2 = as.character(as.numeric(scenario$phenotype_h2[[1]])),
    NG_GRID_EFFECT_TRAINING_N = as.character(as.integer(scenario$effect_training_n[[1]]))
  )
}

ng_with_env_vars <- function(values, expr) {
  env_names <- names(values)
  if (is.null(env_names) || any(!nzchar(env_names))) {
    ng_stop("environment values must be a fully named vector")
  }
  values <- as.character(values)
  names(values) <- env_names
  old <- Sys.getenv(env_names, unset = NA_character_)
  on.exit({
    for (name in names(old)) {
      if (is.na(old[[name]])) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(as.list(old[[name]]), name))
      }
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(values))
  force(expr)
}
