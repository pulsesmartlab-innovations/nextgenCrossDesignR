# Trait-value metric and posterior parameter guide.
#
# This script is self-contained. It creates small phenotype, genotype, map, and
# trait-direction files, then runs ng_run_cross_prediction() with the main
# trait-value and posterior-related settings that users commonly ask about.
#
# WHICH trait_value_metric SHOULD I USE? (from the package's own simulation studies;
# config-scoped, directional -- see docs/BACKEND_USER_GUIDE.md "Choosing A Trait-Value Metric")
#   * Default: "var_complex" (PMV-based usefulness). Never much worse than the best; clearly
#     best for oligogenic traits with good training + heritability. "pmv" and "vpm" rank
#     crosses almost identically -- treat as interchangeable.
#   * "mean": competitive when the trait is highly polygenic OR training is small / low-h2
#     (the within-family variance term adds little there and can add noise).
#   * "le": a linkage-equilibrium diversity/relatedness proxy, NOT a merit metric -- do not select on it.
#   * Long-term recurrent selection: manage diversity EXPLICITLY (lambda_group, or the
#     strategy dial / target_coancestry) on top of a good merit metric, because pure
#     usefulness metrics exhaust genetic variance and gain plateaus. See example 21.
# Prediction accuracy is similar across these metrics; they differ in variance modelling.

library(nextgenCrossDesign)

required_args <- c(
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "n_iter",
  "burn_in",
  "use_parallel"
)
missing_args <- setdiff(required_args, names(formals(nextgenCrossDesign::ng_run_cross_prediction)))
if (length(missing_args)) {
  stop(
    "This example requires nextgenCrossDesign 0.3.11 or newer. Reinstall the current tarball, restart R, ",
    "and confirm packageVersion('nextgenCrossDesign') >= '0.3.11'. Missing arguments in your active package: ",
    paste(missing_args, collapse = ", "),
    call. = FALSE
  )
}

set.seed(20260623)

out_dir <- file.path(tempdir(), "nextgenCrossDesign_trait_value_parameter_guide")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- paste0("P", sprintf("%02d", 1:8))

genotype <- data.frame(
  NAME = ids,
  M01 = c(0, 0, 2, 2, 0, 2, 0, 2),
  M02 = c(0, 0, 2, 2, 0, 2, 0, 2),
  M03 = c(2, 0, 2, 0, 2, 0, 2, 0),
  M04 = c(2, 0, 2, 0, 2, 0, 2, 0),
  M05 = c(0, 2, 0, 2, 0, 2, 2, 0),
  M06 = c(2, 2, 0, 0, 2, 0, 0, 2),
  M07 = c(0, 2, 2, 2, 0, 0, 2, 0),
  M08 = c(2, 2, 2, 0, 2, 0, 0, 2),
  check.names = FALSE
)

phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 60, 53, 68, 55, 64, 59, 63),
  protein = c(11.0, 10.8, 12.1, 10.5, 11.7, 10.9, 11.2, 10.7),
  disease = c(4.0, 3.4, 5.1, 2.5, 4.7, 2.9, 3.6, 3.0),
  stringsAsFactors = FALSE
)

marker_map <- data.frame(
  SNP_code = paste0("M", sprintf("%02d", 1:8)),
  Chromosome = c(1, 1, 1, 1, 2, 2, 2, 2),
  Position_BP = c(0, 1, 3, 6, 0, 2, 4, 8) * 1e6,
  stringsAsFactors = FALSE
)

trait_direction <- data.frame(
  Trait = c("yield", "protein", "disease"),
  PhenotypeColumn = c("yield", "protein", "disease"),
  Selection_direction = c("increase", "increase", "decrease"),
  stringsAsFactors = FALSE
)

phenotype_file <- file.path(out_dir, "phenotype.csv")
genotype_file <- file.path(out_dir, "genotype.csv")
map_file <- file.path(out_dir, "map.csv")
direction_file <- file.path(out_dir, "trait_direction.csv")

write.csv(phenotype, phenotype_file, row.names = FALSE, quote = FALSE)
write.csv(genotype, genotype_file, row.names = FALSE, quote = FALSE)
write.csv(marker_map, map_file, row.names = FALSE, quote = FALSE)
write.csv(trait_direction, direction_file, row.names = FALSE, quote = FALSE)

# -------------------------------------------------------------------------
# User-facing parameter defaults
# -------------------------------------------------------------------------

trait_value_metric <- "var_complex"
uc_variance_source <- "pmv"
selection_prop <- 0.20
progeny <- "DH"
recomb_model <- "haldane"

method_varPMV <- "fast"
method_varPMV_choices <- c("fast", "full_posterior")

ril_mode <- "infinite"
ril_mode_choices <- c("infinite")

run_posterior_prediction <- FALSE
posterior_method <- "mcmc"
posterior_method_choices <- c("closed_form", "mcmc")
n_iter <- 5000
burn_in <- 500
n_draws <- max(1L, n_iter - burn_in)

use_parallel <- FALSE

parameter_notes <- data.frame(
  parameter = c(
    "trait_value_metric",
    "uc_variance_source",
    "selection_prop",
    "progeny",
    "recomb_model",
    "method_varPMV",
    "ril_mode",
    "run_posterior_prediction",
    "posterior_method",
    "n_iter",
    "burn_in",
    "use_parallel"
  ),
  wrapper_argument = c(
    "trait_value_metric",
    "uc_variance_source",
    "selection_prop",
    "progeny",
    "recomb_model",
    "method_varPMV",
    "ril_mode",
    "run_posterior_prediction",
    "posterior_method",
    "n_iter",
    "burn_in",
    "use_parallel"
  ),
  choices = c(
    "var_complex, usefulness, pmv, vpm, le, mean",
    "pmv, vpm, le",
    "numeric proportion such as 0.10 or 0.20",
    "DH, DHs, RIL, RILs",
    "haldane, kosambi",
    paste(method_varPMV_choices, collapse = ", "),
    paste(ril_mode_choices, collapse = ", "),
    "TRUE, FALSE",
    paste(posterior_method_choices, collapse = ", "),
    "positive integer",
    "non-negative integer",
    "TRUE, FALSE"
  ),
  when_to_use = c(
    "Choose the cross-value formula. var_complex is the practical default; mean ignores within-family variance.",
    "Used only when trait_value_metric = 'usefulness'. Use pmv for marker-effect uncertainty, vpm for recombination variance, or le for a diversity proxy.",
    "Controls selection intensity in mean +/- i * SD. Smaller values emphasize the upper tail more strongly.",
    "Choose DH for doubled haploids or RIL for recombinant inbred lines.",
    "Use haldane as the default. Use kosambi when your genetic map and breeding convention require it.",
    "Use fast for normal screens. Use full_posterior for shortlists when off-diagonal marker-effect covariance should enter PMV.",
    "Current RIL implementation is the infinite-selfing RIL variance approximation.",
    "Set TRUE when you want posterior score intervals returned in result$posterior_predictions.",
    "closed_form is faster. mcmc is slower and integrates over variance-component uncertainty.",
    "For posterior runs, total user iterations; n_draws is max(1, n_iter - burn_in).",
    "For posterior MCMC, discarded burn-in iterations. For closed_form, it only affects the n_draws calculation.",
    "Runs trait fits in parallel where the current R session can start workers; otherwise the package records a serial fallback."
  ),
  stringsAsFactors = FALSE
)

metric_grid <- data.frame(
  run_name = c(
    "var_complex_fast_default",
    "uc_pmv_fast",
    "uc_vpm_fast",
    "pmv_full_posterior",
    "vpm_ril_infinite",
    "le_diversity_proxy",
    "mean_only",
    "posterior_closed_form_small"
  ),
  trait_value_metric = c(
    "var_complex",
    "usefulness",
    "usefulness",
    "pmv",
    "vpm",
    "le",
    "mean",
    "pmv"
  ),
  uc_variance_source = c(
    "pmv",
    "pmv",
    "vpm",
    "pmv",
    "pmv",
    "pmv",
    "pmv",
    "pmv"
  ),
  method_varPMV = c(
    "fast",
    "fast",
    "fast",
    "full_posterior",
    "fast",
    "fast",
    "fast",
    "fast"
  ),
  selection_prop = c(
    0.20,
    0.20,
    0.20,
    0.20,
    0.20,
    0.20,
    0.20,
    0.20
  ),
  progeny = c(
    "DH",
    "DH",
    "DH",
    "DH",
    "RIL",
    "DH",
    "DH",
    "DH"
  ),
  ril_mode = c(
    "infinite",
    "infinite",
    "infinite",
    "infinite",
    "infinite",
    "infinite",
    "infinite",
    "infinite"
  ),
  recomb_model = c(
    "haldane",
    "haldane",
    "haldane",
    "haldane",
    "haldane",
    "haldane",
    "haldane",
    "haldane"
  ),
  run_posterior_prediction = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE
  ),
  posterior_method = c(
    posterior_method,
    posterior_method,
    posterior_method,
    posterior_method,
    posterior_method,
    posterior_method,
    posterior_method,
    "closed_form"
  ),
  n_iter = c(
    n_iter,
    n_iter,
    n_iter,
    n_iter,
    n_iter,
    n_iter,
    n_iter,
    8
  ),
  burn_in = c(
    burn_in,
    burn_in,
    burn_in,
    burn_in,
    burn_in,
    burn_in,
    burn_in,
    3
  ),
  use_parallel = c(
    use_parallel,
    use_parallel,
    use_parallel,
    use_parallel,
    use_parallel,
    use_parallel,
    use_parallel,
    use_parallel
  ),
  note = c(
    "Practical default: PMV-first usefulness with fallback columns.",
    "Usefulness criterion using PMV.",
    "Usefulness criterion using recombination variance only.",
    "Uses dense full beta covariance in PMV; best for shortlists.",
    "Uses RIL infinite-selfing variance approximation.",
    "Relationship-distance proxy, not a true progeny variance.",
    "Ranks crosses by mean only.",
    "Returns posterior interval columns in result$posterior_predictions."
  ),
  stringsAsFactors = FALSE
)

run_metric <- function(row) {
  method_out <- file.path(out_dir, "outputs", row$run_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)

  result <- ng_run_cross_prediction(
    phenotype_file = phenotype_file,
    genotype_file = genotype_file,
    map_file = map_file,
    direction_file = direction_file,

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
    multi_trait_method = "auto",
    trait_weights = NULL,

    trait_value_metric = row$trait_value_metric,
    uc_variance_source = row$uc_variance_source,
    selection_prop = row$selection_prop,
    method_varPMV = row$method_varPMV,
    ril_mode = row$ril_mode,
    run_posterior_prediction = row$run_posterior_prediction,
    posterior_method = row$posterior_method,
    n_iter = row$n_iter,
    burn_in = row$burn_in,
    use_parallel = row$use_parallel,

    progeny = row$progeny,
    recomb_model = row$recomb_model,
    assume_inbred = TRUE,

    duplicate_action = "none",
    n_crosses = 4,
    max_crosses_per_parent = 3,
    optimizer = "greedy_local",
    allocation_method = "ocs",
    use_ocs = TRUE,
    lambda_group = 0.05,
    lambda_mating = 0.02,

    output_dir = method_out,
    write_outputs = FALSE,
    write_figures = TRUE,

    seed = 20260623
  )

  stopifnot(inherits(result, "ng_cross_prediction_result"))
  stopifnot(identical(result$settings$trait_value_metric, row$trait_value_metric))
  stopifnot(identical(result$settings$method_varPMV, row$method_varPMV))
  stopifnot(identical(result$settings$ril_mode, row$ril_mode))
  stopifnot(identical(result$settings$run_posterior_prediction, row$run_posterior_prediction))
  stopifnot(nrow(result$selected_crosses) == 4L)
  stopifnot(identical(rownames(result$cleaned_data$genotype), rownames(result$cleaned_data$phenotype)))
  stopifnot(identical(colnames(result$cleaned_data$genotype), result$cleaned_data$marker_map$marker))

  if (identical(row$method_varPMV, "full_posterior")) {
    stopifnot(all(is.finite(result$candidate_crosses$yield_pmv_full_posterior)))
    stopifnot(max(abs(result$candidate_crosses$yield_pmv_used -
                        result$candidate_crosses$yield_pmv_full_posterior), na.rm = TRUE) < 1e-10)
  }

  if (isTRUE(row$run_posterior_prediction)) {
    stopifnot(length(result$posterior_predictions) > 0L)
    stopifnot(all(c("pmv_post_mean", "usefulness_pmv_gebv_post_mean") %in%
                    names(result$posterior_predictions[[1L]])))
  }

  top <- result$selected_crosses[1, , drop = FALSE]
  data.frame(
    run_name = row$run_name,
    status = "completed",
    trait_value_metric = result$settings$trait_value_metric,
    uc_variance_source = result$settings$uc_variance_source,
    method_varPMV = result$settings$method_varPMV,
    progeny = result$settings$progeny,
    ril_mode = result$settings$ril_mode,
    recomb_model = result$settings$recomb_model,
    run_posterior_prediction = result$settings$run_posterior_prediction,
    posterior_method = result$settings$posterior_method,
    posterior_n_draws = result$settings$posterior_n_draws,
    use_parallel = result$settings$use_parallel,
    parallel_backend = result$settings$parallel_backend,
    selected_crosses = nrow(result$selected_crosses),
    top_parent1 = top$parent1,
    top_parent2 = top$parent2,
    top_multi_trait_score = top$multi_trait_score,
    output_dir = normalizePath(method_out, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

trait_value_results <- lapply(seq_len(nrow(metric_grid)), function(i) {
  run_metric(metric_grid[i, , drop = FALSE])
})

trait_value_summary <- do.call(rbind, trait_value_results)
stopifnot(all(trait_value_summary$status == "completed"))

summary_file <- file.path(out_dir, "trait_value_metric_parameter_summary.csv")
write.csv(trait_value_summary, summary_file, row.names = FALSE)

# Lower-level custom-pipeline equivalents. Most users should use
# ng_run_cross_prediction(); these blocks show what the wrapper calls when a
# developer needs exact control over one trait.
run_lower_level_demonstration <- FALSE
if (run_lower_level_demonstration) {
  geno_matrix <- as.matrix(genotype[, setdiff(names(genotype), "NAME"), drop = FALSE])
  storage.mode(geno_matrix) <- "double"
  rownames(geno_matrix) <- genotype$NAME
  y <- setNames(phenotype$yield, phenotype$NAME)

  effects_full <- ng_fit_ridge_effects(
    geno = geno_matrix,
    y = y,
    ids = rownames(geno_matrix),
    seed = 20260623,
    return_beta_cov_full = TRUE
  )

  marker_map_internal <- data.frame(
    marker = marker_map$SNP_code,
    chr = marker_map$Chromosome,
    pos_cm = marker_map$Position_BP / 1e6,
    stringsAsFactors = FALSE
  )

  full_pmv_scores <- ng_score_crosses(
    geno = geno_matrix,
    effects = effects_full,
    marker_map = marker_map_internal,
    ids = rownames(geno_matrix),
    adjusted_pheno = y,
    target = "DH",
    recomb_model = "haldane",
    posterior_cov_full = effects_full$beta_cov_full
  )

  posterior_effects <- ng_fit_ridge_effects_posterior(
    geno = geno_matrix,
    y = y,
    ids = rownames(geno_matrix),
    n_draws = n_draws,
    method = posterior_method,
    mcmc_burnin = burn_in,
    seed = 20260623
  )

  posterior_scores <- ng_posterior_cross_predict(
    geno = geno_matrix,
    posterior_effects = posterior_effects,
    marker_map = marker_map_internal,
    ids = rownames(geno_matrix),
    adjusted_pheno = y,
    target = "DH",
    recomb_model = "haldane",
    selection_prop = selection_prop
  )
}

parameter_notes
metric_grid
trait_value_summary

message("Trait-value metric parameter guide completed.")
message("Summary file: ", normalizePath(summary_file, winslash = "/", mustWork = TRUE))
