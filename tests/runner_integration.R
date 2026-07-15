# End-to-end integration through the user-friendly runner ng_run_cross_prediction():
# the mate-selection module knobs and lethal guarding must reach the frontend-facing path.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(808)
ids <- sprintf("P%02d", 1:14)
mk <- sprintf("M%02d", 1:10)
# inbred lines (homozygous 0/2) so the DH recombination-variance kernel is valid
gm <- matrix(2L * rbinom(length(ids) * length(mk), 1, 0.5),
             nrow = length(ids), dimnames = list(NULL, mk))
# lethal locus M10: exactly P01 and P02 carry the risk allele (homozygous, dosage 2),
# rest are homozygous-safe (0). Carrier x carrier = P01 x P02.
gm[, "M10"] <- 0L; gm[c(1, 2), "M10"] <- 2L
genotype <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids,
                        yield = rnorm(length(ids), 60, 5),
                        stringsAsFactors = FALSE)
marker_map <- data.frame(SNP_code = mk, Chromosome = rep(1:2, each = 5),
                         Position_BP = rep(c(0, 1, 2, 3, 4) * 1e6, 2),
                         stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)

run <- function(...) {
  ng_run_cross_prediction(
    phenotype = phenotype, genotype = genotype, marker_map = marker_map,
    trait_direction = direction, id_col = "NAME",
    map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_col = "Position_BP", map_pos_cm_divisor = 1e6,
    prediction_mode = "trait_by_trait", trait_value_metric = "uc",
    duplicate_action = "none", n_crosses = 8L, max_uses_per_parent = 4L,
    use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L, ...)
}

# --- 1. strategy dial reaches the runner ---
rs <- run(strategy = "diversity")
stopifnot(inherits(rs, "ng_cross_prediction_result"))
stopifnot(identical(rs$plan_summary$strategy, "diversity"))
stopifnot(is.finite(rs$plan_summary$achieved_emphasis))

# --- 2. min-use-if-used reaches the runner ---
rm <- run(min_crosses_per_parent = 3L)
cnt <- table(c(rm$selected_crosses$parent1, rm$selected_crosses$parent2))
stopifnot(all(as.integer(cnt) >= 3L))

# --- 3. lethal guarding: P01 x P02 (carrier x carrier) excluded from the plan ---
rl <- run(lethal_spec = ng_lethal_recessive_spec("M10", risk_allele = "alt"))
sel <- rl$selected_crosses
stopifnot(!any((sel$parent1 == "P01" & sel$parent2 == "P02") |
                 (sel$parent1 == "P02" & sel$parent2 == "P01")))

# --- 4. progeny-inbreeding emphasis reaches the runner and is recorded ---
rp <- run(lambda_progeny_inbreeding = 20)
stopifnot(is.finite(rp$plan_summary$mean_progeny_inbreeding))
stopifnot(isTRUE(all.equal(rp$plan_summary$lambda_progeny_inbreeding, 20)))

# --- 5. the native evolution optimizer is reachable via optimizer= ---
re <- run(optimizer = "evolution", evol_iterations = 30L, evol_solutions = 20L)
stopifnot(inherits(re, "ng_cross_prediction_result"))
stopifnot(nrow(re$selected_crosses) == 8L)

# --- 6. constrained-OCS (target_coancestry) reaches the runner ---
rtc <- run(target_coancestry = 1.0)   # loose cap -> feasible for this small scenario
stc <- rtc$plan_summary
stopifnot(!is.null(stc$target_coancestry_status))
stopifnot(nrow(rtc$selected_crosses) == 8L)

cat("runner integration test passed\n")
