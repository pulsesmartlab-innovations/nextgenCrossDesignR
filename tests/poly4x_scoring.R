helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  message("AlphaSimR unavailable; skipping poly4x scoring test")
  quit(save = "no", status = 0)
}

scenario <- ng_poly4x_select_scenario("potato_autotetraploid_4x")
scenario$n_founders <- 14L
scenario$n_chr <- 2L
scenario$seg_sites <- 24L
scenario$snp_per_chr <- 6L
scenario$qtl_per_chr <- 4L
scenario$genome_length_m <- 0.5
setup <- ng_poly4x_setup_simparam(scenario, seed = 31L, include_digenic = TRUE)
parents <- ng_poly4x_make_parent_pop(setup, n_parents = 8L)
parent_ids <- as.character(parents@id)

duplicate_pairs <- data.frame(
  parent1 = c(parent_ids[[1]], parent_ids[[1]]),
  parent2 = c(parent_ids[[2]], parent_ids[[2]]),
  stringsAsFactors = FALSE
)
duplicate_err <- tryCatch(
  ng_poly4x_score_crosses(
    parent_pop = parents,
    sim_param = setup$sim_param,
    n_score_progeny = 5L,
    candidate_pairs = duplicate_pairs,
    seed = 99L
  ),
  error = function(e) e
)
stopifnot(inherits(duplicate_err, "error"))
stopifnot(grepl("duplicate|repeated cross", conditionMessage(duplicate_err), ignore.case = TRUE))

reciprocal_pairs <- data.frame(
  parent1 = c(parent_ids[[1]], parent_ids[[2]]),
  parent2 = c(parent_ids[[2]], parent_ids[[1]]),
  stringsAsFactors = FALSE
)
reciprocal_err <- tryCatch(
  ng_poly4x_score_crosses(
    parent_pop = parents,
    sim_param = setup$sim_param,
    n_score_progeny = 5L,
    candidate_pairs = reciprocal_pairs,
    seed = 99L
  ),
  error = function(e) e
)
stopifnot(inherits(reciprocal_err, "error"))
stopifnot(grepl("duplicate|repeated cross", conditionMessage(reciprocal_err), ignore.case = TRUE))

scores <- ng_poly4x_score_crosses(
  parent_pop = parents,
  sim_param = setup$sim_param,
  n_score_progeny = 5L,
  selection_prop = 0.20,
  seed = 99L
)

required <- c(
  "parent1", "parent2", "poly4x_mean", "poly4x_var", "poly4x_top10",
  "poly4x_max", "poly4x_usefulness", "poly4x_pair_coancestry",
  "pair_kinship", "poly4x_score_n"
)
stopifnot(all(required %in% names(scores)))
stopifnot(nrow(scores) == 28L)
stopifnot(all(scores$parent1 != scores$parent2))
stopifnot(all(is.finite(scores$poly4x_usefulness)))
stopifnot(all(scores$poly4x_score_n == 5L))

custom_pairs <- data.frame(
  parent1 = c(parent_ids[[1]], parent_ids[[1]], parent_ids[[2]]),
  parent2 = c(parent_ids[[2]], parent_ids[[3]], parent_ids[[4]]),
  stringsAsFactors = FALSE
)
scores_forward <- ng_poly4x_score_crosses(
  parent_pop = parents,
  sim_param = setup$sim_param,
  n_score_progeny = 5L,
  candidate_pairs = custom_pairs,
  seed = 123L
)
scores_reordered <- ng_poly4x_score_crosses(
  parent_pop = parents,
  sim_param = setup$sim_param,
  n_score_progeny = 5L,
  candidate_pairs = custom_pairs[c(3L, 1L, 2L), ],
  seed = 123L
)
pair_key <- function(x) {
  apply(cbind(pmin(x$parent1, x$parent2), pmax(x$parent1, x$parent2)), 1L, paste, collapse = "||")
}
forward_match <- scores_forward[match(pair_key(scores_reordered), pair_key(scores_forward)), ]
score_cols <- c("poly4x_mean", "poly4x_var", "poly4x_top10", "poly4x_max", "poly4x_usefulness")
stopifnot(isTRUE(all.equal(
  forward_match[, score_cols, drop = FALSE],
  scores_reordered[, score_cols, drop = FALSE],
  tolerance = 1e-12,
  check.attributes = FALSE
)))

top_var <- ng_poly4x_var_topn(scores, n_crosses = 3L)
stopifnot(nrow(top_var) == 3L)
stopifnot(identical(unique(top_var$poly4x_method), "ng_poly4x_var_topn"))
stopifnot(!is.null(attr(top_var, "parent_K")))

top_usefulness <- ng_poly4x_usefulness_topn(scores, n_crosses = 3L)
stopifnot(nrow(top_usefulness) == 3L)
stopifnot(identical(unique(top_usefulness$poly4x_method), "ng_poly4x_usefulness_topn"))

parent_K <- attr(scores, "parent_K")
ocs <- ng_poly4x_ocs(
  scores,
  n_crosses = 4L,
  parent_K = parent_K,
  max_crosses_per_parent = 2L,
  lambda_group = 0.5,
  lambda_parent_use = 1.0
)
stopifnot(nrow(ocs) == 4L)
stopifnot(identical(unique(ocs$poly4x_method), "ng_poly4x_ocs"))
counts <- ng_parent_counts(ocs, rownames(parent_K))
stopifnot(max(counts) <= 2L)
stopifnot(all(ocs$parent1 != ocs$parent2))

ocs_from_attr <- ng_poly4x_ocs(
  scores,
  n_crosses = 3L,
  parent_K = NULL,
  max_crosses_per_parent = 2L,
  lambda_group = 0.5,
  lambda_parent_use = 1.0
)
stopifnot(nrow(ocs_from_attr) == 3L)
stopifnot(all(ocs_from_attr$parent1 %in% scores$parent1 | ocs_from_attr$parent1 %in% scores$parent2))
stopifnot(all(ocs_from_attr$parent2 %in% scores$parent1 | ocs_from_attr$parent2 %in% scores$parent2))
stopifnot(identical(unique(ocs_from_attr$poly4x_method), "ng_poly4x_ocs"))

policy_gain <- ng_poly4x_policy(scores, n_crosses = 3L, mode = "gain")
expected_gain <- ng_poly4x_usefulness_topn(scores, n_crosses = 3L)
stopifnot(identical(policy_gain$parent1, expected_gain$parent1))
stopifnot(identical(policy_gain$parent2, expected_gain$parent2))
stopifnot(identical(unique(policy_gain$poly4x_policy_mode), "gain"))
stopifnot(identical(unique(policy_gain$poly4x_policy_scope), "autotetraploid_4x"))

policy_ocs <- ng_poly4x_policy(scores, n_crosses = 4L, mode = "ocs", parent_K = parent_K)
stopifnot(nrow(policy_ocs) == 4L)
stopifnot(max(ng_parent_counts(policy_ocs, rownames(parent_K))) <= 4L)
stopifnot(identical(unique(policy_ocs$poly4x_policy_mode), "ocs"))
stopifnot(identical(attr(policy_ocs, "summary")$poly4x_policy_mode[[1]], "ocs"))

policy_diversity <- ng_poly4x_policy(scores, n_crosses = 4L, mode = "diversity", parent_K = parent_K)
stopifnot(nrow(policy_diversity) == 4L)
stopifnot(max(ng_parent_counts(policy_diversity, rownames(parent_K))) <= 3L)
stopifnot(identical(unique(policy_diversity$poly4x_policy_mode), "diversity"))

invalid_policy <- tryCatch(
  ng_poly4x_policy(scores, n_crosses = 3L, mode = "potato_p20_router"),
  error = function(e) e
)
stopifnot(inherits(invalid_policy, "error"))
stopifnot(grepl("mode", conditionMessage(invalid_policy), ignore.case = TRUE))

poly_gain <- ng_poly4x_policy_select(
  n_parents = 20L,
  crop_scenario = "potato_autotetraploid_4x",
  crop = "potato",
  mode = "gain"
)
stopifnot(identical(poly_gain$method, "poly4x_gain_policy"))
stopifnot(identical(poly_gain$family, "poly4x_policy"))
stopifnot(identical(poly_gain$mode, "gain"))

poly_diversity <- ng_poly4x_policy_select(
  n_parents = 40L,
  crop_scenario = "cassava_autotetraploid_4x",
  crop = "cassava",
  mode = "diversity"
)
stopifnot(identical(poly_diversity$method, "poly4x_diversity_policy"))
stopifnot(identical(poly_diversity$reason, "poly4x_diversity_mode"))

poly_potato_p20 <- ng_poly4x_policy_select(
  n_parents = 20L,
  crop_scenario = "potato_autotetraploid_4x",
  crop = "potato",
  mode = "gain"
)
stopifnot(!identical(poly_potato_p20$method, "poly4x_ocs_policy"))

poly_potato_auto <- ng_poly4x_policy_select(
  n_parents = 40L,
  crop_scenario = "potato_autotetraploid_4x",
  crop = "potato",
  mode = "auto"
)
stopifnot(identical(poly_potato_auto$method, "poly4x_gain_policy"))
stopifnot(identical(poly_potato_auto$mode, "gain"))
stopifnot(grepl("potato", poly_potato_auto$reason, fixed = TRUE))
stopifnot(identical(poly_potato_auto$source, "poly4x_grid_2rep2cycle_20260506"))

poly_potato_p20_auto <- ng_poly4x_policy_select(
  n_parents = 20L,
  crop_scenario = "potato_autotetraploid_4x",
  crop = "potato",
  mode = "auto"
)
stopifnot(identical(poly_potato_p20_auto$method, "poly4x_ocs_policy"))
stopifnot(identical(poly_potato_p20_auto$mode, "ocs"))

poly_cassava_auto <- ng_poly4x_policy_select(
  n_parents = 20L,
  crop_scenario = "cassava_autotetraploid_4x",
  crop = "cassava",
  mode = "auto"
)
stopifnot(identical(poly_cassava_auto$method, "ng_poly4x_usefulness_topn"))
stopifnot(identical(poly_cassava_auto$family, "poly4x_legacy"))
stopifnot(identical(poly_cassava_auto$mode, "legacy_usefulness"))
stopifnot(grepl("cassava", poly_cassava_auto$reason, fixed = TRUE))

poly_cassava_p40_auto <- ng_poly4x_policy_select(
  n_parents = 40L,
  crop_scenario = "cassava_autotetraploid_4x",
  crop = "cassava",
  mode = "auto"
)
stopifnot(identical(poly_cassava_p40_auto$method, "ng_poly4x_ocs"))
stopifnot(identical(poly_cassava_p40_auto$family, "poly4x_legacy"))
stopifnot(identical(poly_cassava_p40_auto$mode, "legacy_ocs"))

poly_auto_fallback <- ng_poly4x_policy_select(
  n_parents = 40L,
  crop_scenario = "potato_autotetraploid_4x",
  crop = "potato",
  mode = "auto",
  available_methods = c("poly4x_diversity_policy")
)
stopifnot(identical(poly_auto_fallback$method, "poly4x_diversity_policy"))
stopifnot(isTRUE(poly_auto_fallback$is_fallback))

cat("poly4x scoring tests passed\n")
