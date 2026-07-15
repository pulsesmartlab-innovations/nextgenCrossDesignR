# User-friendly cross prediction workflow.
#
# Use this script when you have four files:
#   1. phenotype file
#   2. genotype file
#   3. marker map file
#   4. trait direction file
#
# The package handles QC, duplicate-genotype removal, phenotype/genotype
# matching, marker-map matching, marker-effect estimation, cross prediction,
# cross allocation, priority ranking, figures, and the workbook.

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

# -------------------------------------------------------------------------
# 1. Files
# -------------------------------------------------------------------------
# Self-contained: generate a small, consistent multi-trait dataset and write it to four CSVs,
# then treat those exactly as a user's four input files. To use YOUR OWN data instead, point
# the four *_file paths below at your CSVs (same column names) and delete this block.
set.seed(20260623)
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

data_dir <- file.path(tempdir(), "user_friendly_data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- file.path(tempdir(), "user_friendly_run")
phenotype_file <- file.path(data_dir, "phenotype.csv")
genotype_file  <- file.path(data_dir, "genotype.csv")
map_file       <- file.path(data_dir, "marker_map.csv")
direction_file <- file.path(data_dir, "trait_direction.csv")
utils::write.csv(.pheno, phenotype_file, row.names = FALSE)
utils::write.csv(.genotype, genotype_file, row.names = FALSE)
utils::write.csv(.marker_map, map_file, row.names = FALSE)
utils::write.csv(.direction, direction_file, row.names = FALSE)

write_workbook <- requireNamespace("openxlsx", quietly = TRUE)
if (!write_workbook) {
  message("openxlsx is not installed, so crossing_plan.xlsx will be skipped.")
}

# -------------------------------------------------------------------------
# 2. Tell the package how the files match
# -------------------------------------------------------------------------

# Parent ID columns. They can have different names in phenotype and genotype.
phenotype_id_col <- "NAME"
genotype_id_col <- "NAME"

# Direction file columns.
# direction_trait_col is the trait name to show in reports.
# direction_column_col is the phenotype column to model.
# Use NULL when trait names are the same as phenotype column names.
# direction_direction_col is increase/decrease, high/low, or max/min.
direction_trait_col <- "Trait"
direction_column_col <- NULL
direction_direction_col <- "Selection_direction"

# Marker-map columns. Use base-pair positions for the user-facing workflow.
map_marker_col <- "SNP_code"
map_chr_col <- "Chromosome"
map_pos_bp_col <- "Position_BP"
map_position_unit <- "bp"
bp_per_cm <- 1e6

# Use all traits listed in the direction file.
traits_to_use <- NULL

# If only directions are supplied, keep weights NULL. The package then uses an
# equal-weight directional objective and reports that decision in diagnostics.
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
# 4. Cross prediction settings
# -------------------------------------------------------------------------

prediction_mode <- "trait_by_trait"

# Native trait-value choices:
#   "var_complex"  native PopVar-inspired complex usefulness metric.
#                  It uses mean +/- i*SD, prefers PMV variance when available,
#                  and falls back to recombination variance or var_simple.
#   "uc"           mean +/- i*SD using uc_variance_source.
#   "pmv"          usefulness using PMV variance.
#   "vpm"          usefulness using recombination variance.
#   "var_simple"   relationship-distance diversity proxy.
#   "mean"         cross mean only.
trait_value_metric <- "var_complex"
uc_variance_source <- "pmv"
selection_prop <- 0.10

progeny <- "RILs"                # "DH", "DHs", "RIL", or "RILs"
recombination_model <- "haldane" # "haldane" or "kosambi"
assume_inbred <- FALSE
min_effect_reliability <- 0.35

# PMV and posterior settings.
# Keep method_varPMV = "fast" for normal all-pair screens. Use
# "full_posterior" only for smaller candidate sets when you want the full
# off-diagonal marker-effect covariance included in PMV.
method_varPMV <- "fast"          # "fast" or "full_posterior"
ril_mode <- "infinite"           # current RIL variance mode
run_posterior_prediction <- FALSE
posterior_method <- "mcmc"        # "closed_form" or "mcmc"
nIter <- 5000
burnIn <- 500
use_parallel <- FALSE

# -------------------------------------------------------------------------
# 5. Multi-trait objective
# -------------------------------------------------------------------------

multi_trait_method <- "auto"
# Choices:
#   "auto"            equal directional weights when no weights are supplied
#   "weighted"        user-supplied trait_weights
#   "economic_index"  economic weights
#   "desired_gain"    desired-gain targets

threshold_policy <- "soft"  # "soft" or "strict"
threshold_penalty_weight <- 1.0
threshold_penalty_autoscale <- TRUE

# -------------------------------------------------------------------------
# 6. Number of crosses, optimizer, and allocation method
# -------------------------------------------------------------------------

n_crosses <- 100
max_uses_per_parent <- 20
min_unique_parents <- NULL
max_pair_kinship <- Inf

optimizer <- "auto"
# Choices:
#   "auto", "greedy_local", "repair_local", "mip_linear", "mip_contribution"
# Friendly aliases:
#   "lp", "mip", "ocs", "ocs_qp", "egsi"

# Allocation choices:
#   "ocs"                  package OCS-style optimizer
#   "alphamate_style"      native package-developed AlphaMate-style allocation
#   "alphamate_executable" call the external AlphaMate executable
allocation_method <- "alphamate_style"
use_ocs <- TRUE

lambda_group <- 0.05
lambda_mating <- 0.02
lambda_parent_use <- 0
lambda_parent_use_mode <- "absolute"
local_iter <- 2000
ocs_iter <- 5

# Used when allocation_method is "alphamate_style" or "alphamate_executable".
alphamate_mode <- "ModeOptTarget1"
alphamate_target_degree <- 45
alphamate_max_contributions <- max_uses_per_parent
alphamate_number_of_parents <- NULL
alphamate_lambda_group <- NULL
alphamate_lambda_grid <- NULL

# Used only when allocation_method = "alphamate_executable".
alphamate_executable <- NULL
alphamate_runtime_path <- Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")
alphamate_workdir <- NULL
alphamate_keep_files <- FALSE
alphamate_evol_solutions <- 100
alphamate_evol_iterations <- 1000
alphamate_evol_stop <- 200
alphamate_n_threads <- 1

# -------------------------------------------------------------------------
# 7. Priority tiers
# -------------------------------------------------------------------------

priority_breaks <- c(0.10, 0.35, 0.70, 1.00)
priority_labels <- c("highly_priority", "priority", "medium_priority", "low_priority")
priority_score_weight <- 1.0
priority_kinship_weight <- 0.15
priority_threshold_weight <- 1.0

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

  output_dir = out_dir,
  output_file = "crossing_plan.xlsx",
  write_outputs = write_workbook,
  write_figures = TRUE,

  seed = 20260623
)

# -------------------------------------------------------------------------
# 9. Inspect results
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
