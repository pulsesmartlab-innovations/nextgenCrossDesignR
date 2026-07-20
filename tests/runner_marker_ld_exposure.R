# ng_run_cross_prediction() must expose the two mate-design capabilities that the
# documented pipeline ng_design_crosses() already has but the runner lacked:
#   (1) LD pruning of the marker matrix, and
#   (2) marker steering (bias allocation toward a target-allele frequency).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(515)
n <- 18L
base <- 12L; dup <- 6L
G <- matrix(2L * rbinom(n * base, 1, 0.5), n, base)
gm <- cbind(G, G[, seq_len(dup), drop = FALSE])         # 6 exact-copy markers -> r2 = 1
colnames(gm) <- sprintf("M%02d", seq_len(ncol(gm)))
ids <- sprintf("P%02d", seq_len(n)); rownames(gm) <- ids
beta <- rnorm(base, 0, 0.2)
y <- as.numeric(G %*% beta) + rnorm(n)
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE)
mm <- data.frame(SNP_code = colnames(gm), Chromosome = 1L,
                 Position_BP = seq_len(ncol(gm)) * 1e6, stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = mm,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP_code", map_chr_col = "Chromosome",
  map_pos_col = "Position_BP", map_pos_cm_divisor = 1e6,
  prediction_mode = "trait_by_trait", trait_value_metric = "usefulness",
  duplicate_action = "none", n_crosses = 8L, max_crosses_per_parent = 4L,
  use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L, ...)

# --- 1. LD pruning reaches the runner and removes the redundant copy markers ---
r_ld <- run(ld_pruning = TRUE, ld_r2_threshold = 0.95)
stopifnot(!is.null(r_ld$ld_pruning_report))
stopifnot(r_ld$ld_pruning_report$markers_after < r_ld$ld_pruning_report$markers_before)
stopifnot(nrow(r_ld$selected_crosses) == 8L)

# --- 2. marker steering reaches the runner and biases allocation ---
mspec <- ng_marker_target_spec("M01", direction = "increase", weight = 1)
r0 <- run()
rm <- run(marker_target_spec = mspec, lambda_marker = 5.0)
stopifnot("marker_target_score" %in% names(rm$candidate_crosses))

# score BOTH plans on the SAME marker-target table (from the steered run) and confirm the
# steered plan selects crosses with higher mean target-allele score.
pkey <- function(df) ifelse(df$parent1 < df$parent2,
                            paste(df$parent1, df$parent2), paste(df$parent2, df$parent1))
ms <- stats::setNames(rm$candidate_crosses$marker_target_score, pkey(rm$candidate_crosses))
plan_score <- function(sel) mean(ms[pkey(sel)], na.rm = TRUE)
stopifnot(plan_score(rm$selected_crosses) > plan_score(r0$selected_crosses))

cat("runner marker-steering + LD-pruning exposure test passed\n")
