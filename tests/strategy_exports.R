# Frontend export artifacts: frontier JSON, progeny-inbreeding histogram JSON,
# run-comparison JSON.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(606)
np <- 20L
pp <- sprintf("P%02d", seq_len(np))
cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]],
                     stringsAsFactors = FALSE)
scores$uc_dh_gebv <- rnorm(nrow(scores), 10, 2)
L <- matrix(rnorm(np * np, 0, 0.3), np, np)
G <- crossprod(L) / np; diag(G) <- diag(G) + 1; dimnames(G) <- list(pp, pp)
scores$pair_kinship <- G[cbind(match(scores$parent1, pp), match(scores$parent2, pp))]
scores$expected_progeny_inbreeding <- pmax(0, scores$pair_kinship / 2)
n_crosses <- 12L

# --- 1. frontier payload: schema, points, strategy presets ---
fp <- ng_frontier_export_payload(scores, n_crosses, parent_K = G)
stopifnot(identical(fp$schema, "ng_mating_frontier.v1"))
stopifnot(nrow(fp$points) >= 5L)
stopifnot(all(c("diversity_emphasis", "mean_gain", "group_coancestry",
                "mean_progeny_inbreeding") %in% names(fp$points)))
stopifnot(all(c("high_gain", "balanced", "diversity") %in% names(fp$strategies)))
# high_gain preset should have >= gain of the diversity preset
stopifnot(fp$strategies$high_gain$mean_gain >= fp$strategies$diversity$mean_gain - 1e-8)

# --- 2. progeny-inbreeding histogram payload ---
hp <- ng_progeny_inbreeding_histogram_payload(scores, breaks = 10L)
stopifnot(identical(hp$schema, "ng_progeny_inbreeding_histogram.v1"))
stopifnot(sum(hp$bins$count) == nrow(scores))

# --- 3. run-comparison payload ---
r1 <- ng_optimize_mating_plan(scores, n_crosses, parent_K = G, strategy = "high_gain")
r2 <- ng_optimize_mating_plan(scores, n_crosses, parent_K = G, strategy = "diversity")
cp <- ng_run_comparison_payload(list(high_gain = r1, diversity = r2))
stopifnot(identical(cp$schema, "ng_run_comparison.v1"))
stopifnot(length(cp$runs) == 2L)
stopifnot(identical(cp$runs[[1]]$run, "high_gain"))

# --- 4. JSON writers round-trip through jsonlite ---
if (requireNamespace("jsonlite", quietly = TRUE)) {
  td <- tempfile("ngexport_"); dir.create(td)
  f1 <- ng_write_frontier_json(scores, n_crosses, parent_K = G,
                               output_path = file.path(td, "frontier.json"))
  f2 <- ng_write_progeny_inbreeding_histogram_json(scores,
                               output_path = file.path(td, "hist.json"), breaks = 10L)
  f3 <- ng_write_run_comparison_json(list(high_gain = r1, diversity = r2),
                               output_path = file.path(td, "cmp.json"))
  j1 <- jsonlite::read_json(f1); j2 <- jsonlite::read_json(f2); j3 <- jsonlite::read_json(f3)
  stopifnot(identical(j1$schema, "ng_mating_frontier.v1"))
  stopifnot(identical(j2$schema, "ng_progeny_inbreeding_histogram.v1"))
  stopifnot(identical(j3$schema, "ng_run_comparison.v1"))
}

# --- 5. capability registry advertises the new modules ---
reg <- ng_backend_capability_registry()
need_ids <- c("gain_diversity_balance", "progeny_inbreeding_management",
              "mating_constraints", "marker_steering", "lethal_allele_guarding",
              "cost_logistics")
stopifnot(all(need_ids %in% reg$method_families$id))

# --- 6. frontend integration spec exists and covers the capabilities ---
doc <- c(file.path("docs", "frontend", "MATING_STRATEGY_INTEGRATION.md"),
         file.path("nextgen_cross_design", "docs", "frontend", "MATING_STRATEGY_INTEGRATION.md"),
         file.path("..", "docs", "frontend", "MATING_STRATEGY_INTEGRATION.md"))
doc <- doc[file.exists(doc)]
stopifnot(length(doc) >= 1L)
txt <- paste(readLines(doc[[1L]], warn = FALSE), collapse = "\n")
for (tok in c("diversity_emphasis", "ng_mating_frontier.v1",
              "ng_progeny_inbreeding_histogram.v1", "ng_run_comparison.v1",
              "committed_crosses", "group_permission", "min_crosses_per_parent",
              "marker_target_spec", "lethal_recessive", "budget")) {
  if (!grepl(tok, txt, fixed = TRUE)) stop("integration spec missing token: ", tok)
}

cat("strategy exports test passed\n")
