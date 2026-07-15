# Enlarging the marker-effect training set with extra individuals that are NOT candidate parents.
#
# Marker effects (used by every metric except var_simple) are sharper when trained on more
# genotyped + phenotyped individuals. `ng_run_cross_prediction()` lets you add extra "training-only"
# individuals via training_genotype / training_phenotype: they enlarge the ridge fit but are never
# crossed. The real parents are exactly the main genotype/phenotype tables; the training-only
# individuals are whoever is in the training_* tables (any parent-ID overlap is dropped).
#
# This script is self-contained (small simulated data) and runs against the installed package.

library(nextgenCrossDesign)

set.seed(20260706)

n_par <- 15L                 # candidate parents (crossed) -- inbred DH lines, dosage 0/2
n_tr  <- 60L                 # training-only individuals -- may be outbred (0/1/2)
m     <- 300L
markers <- sprintf("SNP%03d", seq_len(m))
par_ids <- sprintf("PAR%02d", seq_len(n_par))
tr_ids  <- sprintf("TRN%02d", seq_len(n_tr))

freq  <- runif(m, 0.15, 0.85)
G_par <- vapply(freq, function(p) 2L * rbinom(n_par, 1, p), integer(n_par))   # inbred 0/2
G_tr  <- vapply(freq, function(p) rbinom(n_tr, 2, p), integer(n_tr))          # outbred 0/1/2
colnames(G_par) <- colnames(G_tr) <- markers

# One heritable trait shared by parents and training individuals.
b  <- rnorm(m) * (runif(m) < 0.15)
tv <- function(G) as.numeric(scale(G %*% b))
genotype   <- data.frame(NAME = par_ids, G_par, check.names = FALSE)
phenotype  <- data.frame(NAME = par_ids, yield = 60 + 5 * (tv(G_par) + rnorm(n_par)))
train_geno <- data.frame(NAME = tr_ids, G_tr, check.names = FALSE)
train_phen <- data.frame(NAME = tr_ids, yield = 60 + 5 * (tv(G_tr) + rnorm(n_tr)))
marker_map <- data.frame(SNP_code = markers, Chromosome = rep(1:6, length.out = m),
                         Position_BP = ave(seq_len(m), rep(1:6, length.out = m), FUN = seq_along) * 1e6)
direction  <- data.frame(Trait = "yield", Selection_direction = "increase")

run <- function(...) ng_run_cross_prediction(
  genotype = genotype, phenotype = phenotype, marker_map = marker_map, trait_direction = direction,
  id_col = "NAME", map_position_unit = "bp", bp_per_cm = 1e6,
  trait_value_metric = "var_complex", optimizer = "evolution", allocation_method = "ocs",
  n_crosses = 20L, max_uses_per_parent = 4L, assume_inbred = TRUE,
  write_outputs = FALSE, write_figures = FALSE, seed = 7L, ...)

# --- Parents only vs parents + training set ---
res_parents <- run()
res_trained <- run(training_genotype = train_geno, training_phenotype = train_phen)

cat("Parents only:\n")
cat("  individuals training the effects:", res_parents$input_match_audit$effect_training_n, "\n")
cat("  effect reliability:", round(res_parents$effect_summary$marker_effect_reliability, 3), "\n")

cat("With extra training individuals:\n")
cat("  candidate parents (crossed):     ", res_trained$input_match_audit$matched_parent_count, "\n")
cat("  training-only individuals:       ", res_trained$input_match_audit$training_only_count, "\n")
cat("  individuals training the effects:", res_trained$input_match_audit$effect_training_n, "\n")
cat("  effect reliability:", round(res_trained$effect_summary$marker_effect_reliability, 3), "\n")

# The crossing plan is built ONLY from the real parents -- no training-only individual is crossed.
sel_ids <- unique(c(res_trained$selected_crosses$parent1, res_trained$selected_crosses$parent2))
stopifnot(all(sel_ids %in% par_ids))
stopifnot(!any(tr_ids %in% sel_ids))
stopifnot(res_trained$input_match_audit$training_only_count == n_tr)
stopifnot(res_trained$input_match_audit$effect_training_n == n_par + n_tr)

cat("\nTop crosses (real parents only):\n")
print(utils::head(res_trained$selected_crosses[, c("parent1", "parent2")], 6), row.names = FALSE)
message("Marker-effect training-set example completed.")
