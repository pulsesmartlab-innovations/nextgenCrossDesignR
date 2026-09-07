helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# DEFECT 5 + DEFECT 7 (multi-trait audit).
#
# 5. ng_posterior_multitrait_cross_predict() (R/32) was exported, documented, tested and
#    release-gated but had ZERO non-test callers, so multi_trait_score -- the merit a
#    multi-trait plan is actually ranked on -- had no posterior at all. Anything wanting an
#    uncertainty on the plan's merit had to reach for posterior_predictions[[1]], i.e. whichever
#    trait sits in row 1 of the breeder's direction file.
# 7. The multi-trait branch of the cross-priority annotation passed neither the posterior SD nor
#    the top-N probability, so prob_top_tier was ABSENT from a multi-trait run's candidate table
#    (present on a single-trait run) and confidence_method always fell back to
#    "midparent_pev_index" even with posterior prediction on.

tmp <- tempfile("ng_runner_mt_index_posterior_")
dir.create(tmp, recursive = TRUE)

ids <- paste0("P", sprintf("%02d", 1:10))
set.seed(4242L)
geno <- data.frame(NAME = ids, matrix(2L * rbinom(10L * 12L, 1L, 0.5), nrow = 10L,
                                      dimnames = list(NULL, sprintf("M%02d", 1:12))),
                   check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(
  NAME = ids,
  yield = c(58, 60, 53, 68, 55, 64, 59, 63, 61, 57),
  disease = c(4.0, 3.4, 5.1, 2.5, 4.7, 2.9, 3.6, 3.0, 3.9, 4.4),
  stringsAsFactors = FALSE
)
marker_map <- data.frame(
  SNP_code = sprintf("M%02d", 1:12),
  Chromosome = rep(1:2, each = 6L),
  Position_BP = rep(c(0, 2, 5, 9, 14, 20), 2) * 1e6,
  stringsAsFactors = FALSE
)
mt_direction <- data.frame(
  Trait = c("yield", "disease"),
  PhenotypeColumn = c("yield", "disease"),
  Selection_direction = c("increase", "decrease"),
  stringsAsFactors = FALSE
)
st_direction <- mt_direction[1L, , drop = FALSE]

files <- list(phenotype = file.path(tmp, "phenotype.csv"),
              genotype = file.path(tmp, "genotype.csv"),
              map = file.path(tmp, "map.csv"),
              direction_mt = file.path(tmp, "trait_direction_mt.csv"),
              direction_st = file.path(tmp, "trait_direction_st.csv"))
write.csv(phenotype, files$phenotype, row.names = FALSE, quote = FALSE)
write.csv(geno, files$genotype, row.names = FALSE, quote = FALSE)
write.csv(marker_map, files$map, row.names = FALSE, quote = FALSE)
write.csv(mt_direction, files$direction_mt, row.names = FALSE, quote = FALSE)
write.csv(st_direction, files$direction_st, row.names = FALSE, quote = FALSE)

run <- function(direction_file, ...) {
  ng_run_cross_prediction(
    phenotype_file = files$phenotype, genotype_file = files$genotype,
    map_file = files$map, direction_file = direction_file,
    phenotype_id_col = "NAME", genotype_id_col = "NAME",
    direction_trait_col = "Trait", direction_column_col = "PhenotypeColumn",
    direction_direction_col = "Selection_direction",
    map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_bp_col = "Position_BP", map_position_unit = "bp", bp_per_cm = 1e6,
    prediction_mode = "trait_by_trait", trait_value_metric = "usefulness",
    uc_variance_source = "pmv", progeny = "DH", recomb_model = "haldane",
    selection_prop = 0.20, duplicate_action = "none",
    posterior_method = "closed_form",
    n_iter = 160L, burn_in = 10L, n_crosses = 5, max_crosses_per_parent = 3,
    optimizer = "greedy_local", write_outputs = FALSE, write_figures = FALSE,
    parent_type = "inbred", use_cpp = FALSE, seed = 909L, ...
  )
}

rq <- 0.25
res_st  <- run(files$direction_st, run_posterior_prediction = TRUE, robustness_quantile = rq)
res_mt  <- run(files$direction_mt, run_posterior_prediction = TRUE, robustness_quantile = rq)
res_off <- run(files$direction_mt, run_posterior_prediction = FALSE)

# ---- DEFECT 7: a multi-trait run carries the columns a single-trait run does -----------------
st_cc <- res_st$candidate_crosses
mt_cc <- res_mt$candidate_crosses
stopifnot("prob_top_tier" %in% names(st_cc))
stopifnot("prob_top_tier" %in% names(mt_cc))
stopifnot(any(is.finite(mt_cc$prob_top_tier)))
stopifnot("prob_top_tier" %in% names(res_mt$selected_crosses))
cat(sprintf("multi-trait prob_top_tier: n_finite = %d, range [%.3f, %.3f]\n",
            sum(is.finite(mt_cc$prob_top_tier)),
            min(mt_cc$prob_top_tier, na.rm = TRUE), max(mt_cc$prob_top_tier, na.rm = TRUE)))

# Posterior-derived confidence, not the mid-parent-PEV fallback, on BOTH run shapes.
cat(sprintf("confidence_method: single-trait = %s, multi-trait = %s, multi-trait posterior off = %s\n",
            st_cc$confidence_method[[1L]], mt_cc$confidence_method[[1L]],
            res_off$candidate_crosses$confidence_method[[1L]]))
stopifnot(identical(mt_cc$confidence_method[[1L]], "posterior_ci"))
stopifnot(isTRUE(res_mt$priority_risk_diagnostics$posterior_used))
# Posterior OFF still falls back to the index PEV path and still emits the column (as NA).
stopifnot("prob_top_tier" %in% names(res_off$candidate_crosses))
stopifnot(all(is.na(res_off$candidate_crosses$prob_top_tier)))
stopifnot(grepl("^midparent_pev_index", res_off$candidate_crosses$confidence_method[[1L]]))

# ---- DEFECT 5: the index posterior is wired and reaches the result ---------------------------
stopifnot(is.null(res_st$posterior_multitrait))    # single trait: no index to build
stopifnot(is.null(res_off$posterior_multitrait))   # posterior off: not paid for
pm <- res_mt$posterior_multitrait
stopifnot(!is.null(pm))
stopifnot(nrow(pm) == nrow(mt_cc))
need <- c("multi_trait_score_post_mean", "multi_trait_score_post_lower",
          "multi_trait_score_post_upper", "multi_trait_score_post_sd")
stopifnot(all(need %in% names(pm)))
stopifnot(all(need %in% names(mt_cc)))            # rides the candidate table
stopifnot(all(is.finite(pm$multi_trait_score_post_mean)))
stopifnot(all(pm$multi_trait_score_post_lower <= pm$multi_trait_score_post_upper))

# Additive: nothing existing was removed or renamed.
stopifnot(all(c("multi_trait_score", "multi_trait_threshold_violation",
                "yield_post_sd", "disease_post_sd") %in% names(mt_cc)))
stopifnot(identical(names(res_mt$posterior_predictions), c("yield", "disease")))

# ---- robustness_quantile: exact empirical quantiles cached from the SAME draws ---------------
q_lo <- ng_posterior_quantile_col("multi_trait_score", rq)
q_hi <- ng_posterior_quantile_col("multi_trait_score", 1 - rq)
stopifnot(all(c(q_lo, q_hi) %in% names(pm)))
stopifnot(all(pm[[q_lo]] <= pm[[q_hi]]))
# The reported credible interval is untouched by the robust request.
meta <- attr(pm, "posterior_multitrait")
stopifnot(identical(meta$ci_level, 0.95))
stopifnot(isTRUE(all.equal(meta$robustness_quantile, rq)))
stopifnot(is.data.frame(meta$posterior_quantiles))
stopifnot(setequal(meta$posterior_quantiles$column, c(q_lo, q_hi)))
res_mt_plain <- run(files$direction_mt, run_posterior_prediction = TRUE)
stopifnot(!any(c(q_lo, q_hi) %in% names(res_mt_plain$posterior_multitrait)))

# ---- direction: the index is normalised higher = better, so it is fixed at "maximize" --------
stopifnot(identical(meta$direction, "maximize"))
stopifnot(identical(attr(pm, "posterior")$direction, "maximize"))
stopifnot(identical(attr(pm, "posterior")$gain_col, "multi_trait_score"))
# EVIDENCE, not assertion: the emitted index really does point up for both trait directions.
# yield is an increase trait, disease a decrease trait.
cat(sprintf("cor(multi_trait_score, yield_mean) = %+.3f ; cor(multi_trait_score, disease_mean) = %+.3f\n",
            cor(mt_cc$multi_trait_score, mt_cc$yield_mean),
            cor(mt_cc$multi_trait_score, mt_cc$disease_mean)))
stopifnot(cor(mt_cc$multi_trait_score, mt_cc$yield_mean) > 0)
stopifnot(cor(mt_cc$multi_trait_score, mt_cc$disease_mean) < 0)

# ---- the per-draw re-standardisation caveat is RECORDED, not silently absorbed ---------------
stopifnot(identical(meta$index_rescaling, "per_draw_restandardized"))
stopifnot(is.character(meta$index_rescaling_note) && nzchar(meta$index_rescaling_note))

# ---- the index posterior is usable by the robust allocator -----------------------------------
# This is the capability the whole fix exists for: robust-allocate on the INDEX the plan ranks
# on, instead of on trait 1's per-trait posterior.
robust <- ng_optimize_robust_mating_plan(
  posterior_scores = pm, n_crosses = 5L,
  gain_col = "multi_trait_score", robustness_quantile = rq, direction = "maximize",
  max_crosses_per_parent = 3, method = "greedy_local")
rs <- attr(robust, "summary")
stopifnot(nrow(robust) == 5L)
stopifnot(identical(rs$robustness_quantile_is_normal_approximation, FALSE))
cat(sprintf("robust plan on the INDEX: quantile source = %s (exact, not a normal approximation)\n",
            rs$robust_quantile_source))
stopifnot(identical(rs$robust_quantile_source, q_lo))

cat("runner_multitrait_index_posterior: OK\n")
