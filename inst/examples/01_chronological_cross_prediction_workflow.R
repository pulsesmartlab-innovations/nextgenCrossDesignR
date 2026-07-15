# Chronological nextgenCrossDesign workflow.
#
# This is written in the order a breeder or analyst normally works:
#   1. Tell the package where the files are.
#   2. Tell the package which columns identify parents, traits, markers, and map positions.
#   3. Choose QC, cross-prediction, optimizer, allocation, and priority settings.
#   4. Run the package workflow.
#   5. Inspect the QC/matching audit and the crossing plan.
#
# The user does not manually intersect phenotype/genotype IDs or marker names.
# The package does that inside QC and stops if the files do not match.

library(nextgenCrossDesign)

required_args <- c(
  "allocation_method",
  "alphamate_mode",
  "alphamate_target_degree",
  "method_varPMV",
  "ril_mode",
  "run_posterior_prediction",
  "posterior_method",
  "nIter",
  "burnIn",
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

# -------------------------------------------------------------------------
# 1. Files
# -------------------------------------------------------------------------
# Self-contained: generate a small, consistent multi-trait dataset and write it to four CSVs,
# then treat those exactly as a user's four input files. To use YOUR OWN data instead, point
# the four *_file paths below at your CSVs (same column names) and delete this block.
.np <- 40L; .nchr <- 5L; .nsnp <- 40L
.ids <- sprintf("L%03d", seq_len(.np))
.markers <- sprintf("M%04d", seq_len(.nchr * .nsnp))
.G <- matrix(2L * rbinom(.np * length(.markers), 1, 0.5), .np,      # inbred (RIL) parents: 0/2
             dimnames = list(.ids, .markers))
.traits <- c("yield", "protein", "disease")                        # disease is a "decrease" trait
.pheno <- data.frame(NAME = .ids, stringsAsFactors = FALSE)
for (.tr in .traits) {
  .gv <- as.numeric(.G %*% rnorm(length(.markers), 0, 0.1))
  .pheno[[.tr]] <- round(60 + .gv + rnorm(.np, 0, stats::sd(.gv)), 3)
}
.genotype  <- data.frame(NAME = .ids, .G, check.names = FALSE, stringsAsFactors = FALSE)
.marker_map <- data.frame(SNP_code = .markers, Chromosome = rep(seq_len(.nchr), each = .nsnp),
                          Position_BP = rep(seq(0, by = 1e6, length.out = .nsnp), .nchr),
                          stringsAsFactors = FALSE)
.direction <- data.frame(Trait = .traits, Selection_direction = c("increase", "increase", "decrease"),
                         stringsAsFactors = FALSE)

data_dir <- file.path(tempdir(), "chronological_data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
output_dir <- file.path(tempdir(), "chronological_user_example")
phenotype_file <- file.path(data_dir, "phenotype.csv")
genotype_file  <- file.path(data_dir, "genotype.csv")
map_file       <- file.path(data_dir, "marker_map.csv")
direction_file <- file.path(data_dir, "trait_direction.csv")
utils::write.csv(.pheno, phenotype_file, row.names = FALSE)
utils::write.csv(.genotype, genotype_file, row.names = FALSE)
utils::write.csv(.marker_map, map_file, row.names = FALSE)
utils::write.csv(.direction, direction_file, row.names = FALSE)

write_workbook <- requireNamespace("openxlsx", quietly = TRUE)
write_figures <- TRUE

# -------------------------------------------------------------------------
# 2. Column names and matching contract
# -------------------------------------------------------------------------

# Parent IDs. These can be different column names in phenotype and genotype.
phenotype_id_col <- "NAME"
genotype_id_col <- "NAME"

# Direction file columns.
# direction_trait_col: report-friendly trait name.
# direction_column_col: phenotype column to model. Use NULL if it is the same
# as the trait name.
# direction_direction_col: increase/decrease, high/low, or max/min.
direction_trait_col <- "Trait"
direction_column_col <- NULL
direction_direction_col <- "Selection_direction"

# Marker map columns. The recommended input unit is base pairs.
map_marker_col <- "SNP_code"
map_chr_col <- "Chromosome"
map_pos_bp_col <- "Position_BP"
map_position_unit <- "bp"
bp_per_cm <- 1e6

# Use NULL to use all traits in the direction file.
traits_to_use <- NULL

# Use NULL when only directions are supplied. If the breeder has weights, use a
# named vector, for example c(YIELD = 0.40, Protein = 0.25).
trait_weights <- NULL

# -------------------------------------------------------------------------
# 3. QC settings
# -------------------------------------------------------------------------

duplicate_action <- "remove"  # "remove", "report", or "none"
duplicate_threshold <- 0.995
duplicate_maf_min <- 0.01
duplicate_max_missing_prop <- 0.40
duplicate_min_compared_markers <- 100

# -------------------------------------------------------------------------
# 4. Cross-prediction method settings
# -------------------------------------------------------------------------

prediction_mode <- "trait_by_trait"
# Use "trait_by_trait" for normal multi-trait parent selection.
# Use "index_as_trait" only when the phenotype file already has a trusted
# selection-index column.

variance_method <- "var_complex"
# Choices:
#   "uc"         mean +/- selection intensity * within-family SD
#   "pmv"        usefulness using PMV as the variance source
#   "vpm"        usefulness using recombination variance as the variance source
#   "var_complex" native PopVar-inspired complex usefulness metric
#   "var_simple" relationship-distance proxy for screening
#   "mean"       cross mean only
trait_value_metric <- variance_method

uc_variance_source <- "pmv"
# Used when variance_method = "uc".
# Choices: "pmv", "vpm", "var_simple"

selection_prop <- 0.10
progeny <- "RILs"                  # "DH", "DHs", "RIL", or "RILs"
recombination_model <- "haldane"   # "haldane" or "kosambi"
assume_inbred <- FALSE
min_effect_reliability <- 0.35

# Optional PMV/posterior knobs are passed directly to ng_run_cross_prediction().
# Keep the fast/default values for normal all-pair screens.
method_varPMV <- "fast"
ril_mode <- "infinite"
run_posterior_prediction <- FALSE
posterior_method <- "mcmc"
nIter <- 5000
burnIn <- 500
use_parallel <- FALSE

# -------------------------------------------------------------------------
# 5. Multi-trait objective settings
# -------------------------------------------------------------------------

multi_trait_method <- "auto"
# Choices:
#   "auto"            equal directional weights when no weights are supplied
#   "weighted"        user-supplied trait_weights
#   "economic_index"  economic weights
#   "desired_gain"    desired-gain targets

threshold_policy <- "soft"   # "soft" or "strict"
threshold_penalty_weight <- 1.0
threshold_penalty_autoscale <- TRUE

# -------------------------------------------------------------------------
# 6. Optimizer, number of crosses, and allocation settings
# -------------------------------------------------------------------------

max_crosses <- 100
n_crosses <- max_crosses
max_uses_per_parent <- 20
min_unique_parents <- NULL
max_pair_kinship <- Inf
min_variance <- NULL

optimizer <- "auto"
# Choices:
#   "auto"
#   "greedy_local"
#   "repair_local"
#   "mip_linear"
#   "mip_contribution"
# Friendly aliases also accepted: "lp", "mip", "ocs", "ocs_qp", "egsi"

solver <- "lpSolve"
use_ocs <- TRUE

allocation_method <- "alphamate_style"
# Choices:
#   "ocs"                  package OCS-style optimizer
#   "alphamate_style"      native package-developed AlphaMate-style allocation
#   "alphamate_executable" call the external AlphaMate executable

lambda_group <- 0.05
lambda_mating <- 0.02
lambda_parent_use <- 0
lambda_parent_use_mode <- "absolute"  # "absolute" or "adaptive"

local_iter <- 2000
ocs_iter <- 5

alphamate_mode <- "ModeOptTarget1"
# Choices:
#   "ModeOptTarget1"       balance criterion and target coancestry degree
#   "ModeMaxCriterion"     emphasize gain/criterion
#   "ModeMinCoancestry"    emphasize diversity

alphamate_target_degree <- 45
alphamate_max_contributions <- max_uses_per_parent
alphamate_number_of_parents <- NULL
alphamate_lambda_group <- NULL
alphamate_lambda_grid <- NULL

# Only used when allocation_method = "alphamate_executable".
alphamate_executable <- NULL
alphamate_runtime_path <- Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")
alphamate_workdir <- NULL
alphamate_keep_files <- FALSE
alphamate_evol_solutions <- 100
alphamate_evol_iterations <- 1000
alphamate_evol_stop <- 200
alphamate_n_threads <- 1

# -------------------------------------------------------------------------
# 7. Priority ranking and optional cross-number sweep
# -------------------------------------------------------------------------

priority_breaks <- c(0.10, 0.35, 0.70, 1.00)
priority_labels <- c("highly_priority", "priority", "medium_priority", "low_priority")
priority_score_weight <- 1.0
priority_kinship_weight <- 0.15
priority_threshold_weight <- 1.0

run_cross_number_sweep <- TRUE
k_range <- seq(50, 200, by = 10)
cross_number_criterion <- "elbow_relative"
relative_threshold <- 0.05
ne_min <- 30
coancestry_max <- 0.05

# -------------------------------------------------------------------------
# 8. Run the package workflow
# -------------------------------------------------------------------------

result <- ng_run_cross_prediction(
  phenotype_file = phenotype_file,
  genotype_file = genotype_file,
  map_file = map_file,
  direction_file = direction_file,

  phenotype_id_col = phenotype_id_col,
  genotype_id_col = genotype_id_col,
  direction_trait_col = direction_trait_col,
  direction_column_col = direction_column_col,
  direction_direction_col = direction_direction_col,
  map_marker_col = map_marker_col,
  map_chr_col = map_chr_col,
  map_pos_bp_col = map_pos_bp_col,
  map_position_unit = map_position_unit,
  bp_per_cm = bp_per_cm,

  prediction_mode = prediction_mode,
  traits_to_use = traits_to_use,
  multi_trait_method = multi_trait_method,
  trait_weights = trait_weights,
  threshold_policy = threshold_policy,
  threshold_penalty_weight = threshold_penalty_weight,
  threshold_penalty_autoscale = threshold_penalty_autoscale,

  trait_value_metric = trait_value_metric,
  uc_variance_source = uc_variance_source,
  selection_prop = selection_prop,
  min_effect_reliability = min_effect_reliability,
  method_varPMV = method_varPMV,
  ril_mode = ril_mode,
  run_posterior_prediction = run_posterior_prediction,
  posterior_method = posterior_method,
  nIter = nIter,
  burnIn = burnIn,
  use_parallel = use_parallel,

  progeny = progeny,
  recombination_model = recombination_model,
  assume_inbred = assume_inbred,

  duplicate_action = duplicate_action,
  duplicate_threshold = duplicate_threshold,
  duplicate_maf_min = duplicate_maf_min,
  duplicate_max_missing_prop = duplicate_max_missing_prop,
  duplicate_min_compared_markers = duplicate_min_compared_markers,

  n_crosses = n_crosses,
  max_uses_per_parent = max_uses_per_parent,
  min_unique_parents = min_unique_parents,
  max_pair_kinship = max_pair_kinship,
  optimizer = optimizer,
  allocation_method = allocation_method,
  use_ocs = use_ocs,
  lambda_group = lambda_group,
  lambda_mating = lambda_mating,
  lambda_parent_use = lambda_parent_use,
  lambda_parent_use_mode = lambda_parent_use_mode,
  local_iter = local_iter,
  ocs_iter = ocs_iter,
  alphamate_mode = alphamate_mode,
  alphamate_target_degree = alphamate_target_degree,
  alphamate_max_contributions = alphamate_max_contributions,
  alphamate_number_of_parents = alphamate_number_of_parents,
  alphamate_lambda_group = alphamate_lambda_group,
  alphamate_lambda_grid = alphamate_lambda_grid,
  alphamate_executable = alphamate_executable,
  alphamate_runtime_path = alphamate_runtime_path,
  alphamate_workdir = alphamate_workdir,
  alphamate_keep_files = alphamate_keep_files,
  alphamate_evol_solutions = alphamate_evol_solutions,
  alphamate_evol_iterations = alphamate_evol_iterations,
  alphamate_evol_stop = alphamate_evol_stop,
  alphamate_n_threads = alphamate_n_threads,

  priority_breaks = priority_breaks,
  priority_labels = priority_labels,
  priority_score_weight = priority_score_weight,
  priority_kinship_weight = priority_kinship_weight,
  priority_threshold_weight = priority_threshold_weight,

  output_dir = output_dir,
  output_file = "breeder_crossing_plan.xlsx",
  write_outputs = write_workbook,
  write_figures = write_figures,

  seed = 20260623
)

# -------------------------------------------------------------------------
# 9. Inspect QC, matching, marker effects, and crossing decisions
# -------------------------------------------------------------------------

result$qc$status
result$qc$issues
result$input_match_audit
result$effect_summary
head(result$candidate_crosses)
head(result$selected_crosses)
result$plan_summary
result$output_files

geno_clean <- result$cleaned_data$genotype
pheno_clean <- result$cleaned_data$phenotype
map_clean <- result$cleaned_data$marker_map

stopifnot(identical(rownames(geno_clean), rownames(pheno_clean)))
stopifnot(identical(colnames(geno_clean), map_clean$marker))
stopifnot(all(c("pos_bp", "pos_cm") %in% names(map_clean)))

# -------------------------------------------------------------------------
# 10. Optional: choose the number of crosses by diminishing returns
# -------------------------------------------------------------------------

if (isTRUE(run_cross_number_sweep)) {
  parent_K <- ng_parent_kinship(geno_clean)
  sweep_optimizer <- switch(
    optimizer,
    lp = "mip_linear",
    mip = "mip_contribution",
    ocs = "mip_contribution",
    ocs_qp = "mip_contribution",
    egsi = "greedy_local",
    optimizer
  )

  curve <- ng_optimize_mating_plan_curve(
    scores = result$candidate_crosses,
    K_range = k_range,
    gain_col = "multi_trait_score",
    parent_K = parent_K,
    max_crosses_per_parent = max_uses_per_parent,
    lambda_group = lambda_group,
    lambda_mating = lambda_mating,
    method = sweep_optimizer,
    criterion = cross_number_criterion,
    relative_threshold = relative_threshold,
    ne_min = ne_min,
    coancestry_max = coancestry_max
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(curve, file.path(output_dir, "cross_number_sweep.csv"), row.names = FALSE)

  if (isTRUE(write_figures)) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ng_plot_diminishing_returns(curve)
      ggplot2::ggsave(
        filename = file.path(output_dir, "diminishing_returns.png"),
        plot = p,
        width = 7,
        height = 4.5,
        dpi = 160
      )
    }
  }

  attr(curve, "elbow_K")
}

# -------------------------------------------------------------------------
# 11. Advanced lower-level API map
# -------------------------------------------------------------------------
#
# The one-call workflow above is the recommended user path. These are the lower
# level functions behind it when a developer needs to build a custom pipeline:
#
# ng_preflight_input_tables()
# ng_score_crosses()
# ng_breeder_selection_objective()
# ng_optimize_breeder_selection_plan()
# ng_alphamate_style_select()
# ng_select_alphamate()
# ng_rank_cross_priority()
# ng_write_cross_priority_workbook()

run_settings <- data.frame(
  setting = c(
    "prediction_mode", "trait_value_metric", "uc_variance_source",
    "multi_trait_method", "optimizer", "solver", "allocation_method", "use_ocs",
    "n_crosses", "max_uses_per_parent", "progeny", "recombination_model",
    "map_position_unit", "bp_per_cm", "method_varPMV", "ril_mode",
    "alphamate_mode", "alphamate_target_degree", "alphamate_max_contributions",
    "run_posterior_prediction", "posterior_method", "nIter", "burnIn",
    "use_parallel"
  ),
  value = as.character(c(
    prediction_mode, trait_value_metric, uc_variance_source,
    multi_trait_method, optimizer, solver, allocation_method, use_ocs,
    n_crosses, max_uses_per_parent, progeny, recombination_model,
    map_position_unit, bp_per_cm, method_varPMV, ril_mode,
    alphamate_mode, alphamate_target_degree, alphamate_max_contributions,
    run_posterior_prediction, posterior_method, nIter, burnIn,
    use_parallel
  )),
  stringsAsFactors = FALSE
)

write.csv(run_settings, file.path(output_dir, "run_settings.csv"), row.names = FALSE)

message("Done.")
message("Output directory: ", normalizePath(output_dir, winslash = "/", mustWork = TRUE))
message("Selected crosses: ", nrow(result$selected_crosses))
message("Candidate crosses scored: ", nrow(result$candidate_crosses))
