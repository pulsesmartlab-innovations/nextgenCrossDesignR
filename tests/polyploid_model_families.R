helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = function(e) e)
  stopifnot(inherits(err, "error"))
  stopifnot(grepl(pattern, conditionMessage(err), ignore.case = TRUE))
  invisible(err)
}

wheat_decision <- ng_polyploid_model_select(
  crop = "wheat",
  crop_scenario = "bread_wheat_hexaploid_approx",
  inheritance_model = "disomic_subgenome",
  subgenome_names = c("A", "B", "D")
)
stopifnot(identical(wheat_decision$model_family, "allopolyploid_subgenome"))
stopifnot(identical(wheat_decision$inheritance_model, "disomic_subgenome"))
stopifnot(isTRUE(wheat_decision$supported))
stopifnot(identical(wheat_decision$subgenome_names, c("A", "B", "D")))
stopifnot(identical(wheat_decision$ploidy_profile, c(A = 2L, B = 2L, D = 2L)))
stopifnot(identical(wheat_decision$supported_scope, "wheat_like_allopolyploid_disomic"))

potato_decision <- ng_polyploid_model_select(
  crop = "potato",
  crop_scenario = "potato_autotetraploid_4x",
  inheritance_model = "polysomic",
  ploidy = 4L
)
stopifnot(identical(potato_decision$model_family, "autotetraploid_4x"))
stopifnot(identical(potato_decision$supported_scope, "autotetraploid_4x"))
stopifnot(isTRUE(potato_decision$supported))

sugarcane_decision <- ng_polyploid_model_select(
  crop = "sugarcane",
  crop_scenario = "sugarcane_polyploid_stress",
  inheritance_model = "mixed_or_aneuploid"
)
stopifnot(identical(sugarcane_decision$model_family, "complex_polyploid_guard"))
stopifnot(identical(sugarcane_decision$inheritance_model, "mixed_or_aneuploid"))
stopifnot(identical(sugarcane_decision$supported_scope, "requires_empirical_or_external_simulator"))
stopifnot(!isTRUE(sugarcane_decision$supported))

expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "disomic_subgenome"),
  "subgenome"
)
expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "disomic_subgenome", subgenome_names = c("A", "A")),
  "unique"
)
expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "disomic_subgenome", subgenome_names = c("A", "B"), ploidy_profile = c(A = 2L, B = 4L)),
  "diploid"
)
expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "disomic_subgenome", subgenome_names = c("A", "B"), ploidy_profile = c(A = 2.9, B = 2.1)),
  "integer"
)
expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "disomic_subgenome", subgenome_names = c("A", "B"), quad_prob = 0.1),
  "quadrivalent"
)
wheat_null_quad <- ng_polyploid_model_select(
  crop = "wheat",
  inheritance_model = "disomic_subgenome",
  subgenome_names = c("A", "B"),
  quad_prob = NULL
)
stopifnot(identical(wheat_null_quad$model_family, "allopolyploid_subgenome"))
expect_error(
  ng_polyploid_model_select(crop = "wheat", inheritance_model = "unknown", subgenome_names = c("A", "B", "D")),
  "inheritance_model"
)
expect_error(
  ng_polyploid_model_select(crop = "potato", inheritance_model = "polysomic", ploidy = 4.9),
  "integer"
)
expect_error(
  ng_polyploid_model_select(crop = "mystery", inheritance_model = "unknown"),
  "inheritance_model"
)

geno_by_subgenome <- list(
  A = matrix(c(0, 0, 2, 2, 1, 1), nrow = 3, byrow = TRUE,
             dimnames = list(c("P1", "P2", "P3"), c("A_m1", "A_m2"))),
  B = matrix(c(0, 2, 2, 0, 1, 1), nrow = 3, byrow = TRUE,
             dimnames = list(c("P1", "P2", "P3"), c("B_m1", "B_m2")))
)

checked_geno <- ng_polyploid_subgenome_as_dosage_list(geno_by_subgenome)
stopifnot(identical(names(checked_geno), c("A", "B")))
stopifnot(identical(rownames(checked_geno$A), c("P1", "P2", "P3")))
stopifnot(storage.mode(checked_geno$A) == "double")

unnamed <- geno_by_subgenome
names(unnamed) <- NULL
expect_error(ng_polyploid_subgenome_as_dosage_list(unnamed), "named")

missing_subgenome_name <- geno_by_subgenome
names(missing_subgenome_name) <- c("A", NA_character_)
expect_error(ng_polyploid_subgenome_as_dosage_list(missing_subgenome_name), "named")

blank_subgenome_name <- geno_by_subgenome
names(blank_subgenome_name) <- c("A", "  ")
expect_error(ng_polyploid_subgenome_as_dosage_list(blank_subgenome_name), "named")

bad_dosage <- geno_by_subgenome
bad_dosage$A[1, 1] <- 3
expect_error(ng_polyploid_subgenome_as_dosage_list(bad_dosage), "0..2")

bad_rows <- geno_by_subgenome
rownames(bad_rows$B) <- c("P1", "P2", "PX")
expect_error(ng_polyploid_subgenome_as_dosage_list(bad_rows), "parent row names")

na_rows <- geno_by_subgenome
rownames(na_rows$A) <- c("P1", NA_character_, "P3")
expect_error(ng_polyploid_subgenome_as_dosage_list(na_rows), "parent row names")

blank_cols <- geno_by_subgenome
colnames(blank_cols$A) <- c("A_m1", " ")
expect_error(ng_polyploid_subgenome_as_dosage_list(blank_cols), "marker column names")

missing_cols <- geno_by_subgenome
colnames(missing_cols$A) <- NULL
expect_error(ng_polyploid_subgenome_as_dosage_list(missing_cols), "marker column names")

K <- ng_polyploid_subgenome_grm(geno_by_subgenome)
expected_K <- matrix(
  c(1, -1, 0,
    -1, 1, 0,
    0, 0, 0),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("P1", "P2", "P3"), c("P1", "P2", "P3"))
)
stopifnot(isTRUE(all.equal(unclass(K), expected_K, tolerance = 1e-12, check.attributes = FALSE)))
stopifnot(identical(names(attr(K, "subgenome_K")), c("A", "B")))
stopifnot(isTRUE(all.equal(attr(K, "subgenome_weights"), c(A = 1, B = 1), tolerance = 1e-12)))

weighted_K <- ng_polyploid_subgenome_grm(geno_by_subgenome, weights = c(A = 2, B = 1))
stopifnot(isTRUE(all.equal(unclass(weighted_K), expected_K, tolerance = 1e-12, check.attributes = FALSE)))

weighted_fixture <- list(
  A = geno_by_subgenome$A,
  B = matrix(c(0, 0, 1, 1, 2, 2), nrow = 3, byrow = TRUE,
             dimnames = list(c("P1", "P2", "P3"), c("B_m1", "B_m2")))
)
weighted_distinct <- ng_polyploid_subgenome_grm(weighted_fixture, weights = c(A = 2, B = 1))
expected_weighted_distinct <- matrix(
  c(1, -2 / 3, -1 / 3,
    -2 / 3, 2 / 3, 0,
    -1 / 3, 0, 1 / 3),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("P1", "P2", "P3"), c("P1", "P2", "P3"))
)
stopifnot(isTRUE(all.equal(
  unclass(weighted_distinct),
  expected_weighted_distinct,
  tolerance = 1e-12,
  check.attributes = FALSE
)))

bad_weights <- c(A = 1, C = 1)
expect_error(ng_polyploid_subgenome_grm(geno_by_subgenome, weights = bad_weights), "weights")

effects_by_subgenome <- list(
  A = c(A_m1 = 0.50, A_m2 = -0.25),
  B = c(B_m1 = 0.25, B_m2 = 0.75)
)
candidate_pairs <- data.frame(
  parent1 = c("P1", "P1", "P2"),
  parent2 = c("P2", "P3", "P3"),
  stringsAsFactors = FALSE
)

scores <- ng_polyploid_subgenome_score_crosses(
  geno_by_subgenome = geno_by_subgenome,
  effects_by_subgenome = effects_by_subgenome,
  candidate_pairs = candidate_pairs,
  model_decision = wheat_decision,
  selection_prop = 0.10
)
required_score_cols <- c(
  "parent1", "parent2", "poly_gain", "poly_var", "poly_usefulness",
  "pair_kinship", "polyploid_model_family", "inheritance_model",
  "subgenome_count", "subgenome_names", "supported_scope",
  "polyploid_validation_source"
)
stopifnot(all(required_score_cols %in% names(scores)))
stopifnot(nrow(scores) == 3L)
stopifnot(isTRUE(all.equal(scores$poly_gain[[1]], 1.25, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(scores$poly_var[[1]], 0.9375, tolerance = 1e-12)))
stopifnot(all(is.finite(scores$poly_usefulness)))
stopifnot(identical(unique(scores$polyploid_model_family), "allopolyploid_subgenome"))
stopifnot(identical(unique(scores$inheritance_model), "disomic_subgenome"))
stopifnot(identical(unique(scores$subgenome_names), "A|B"))
stopifnot(!is.null(attr(scores, "parent_kinship")))

bad_effects <- effects_by_subgenome
bad_effects$A <- c(A_m1 = 0.5)
expect_error(
  ng_polyploid_subgenome_score_crosses(geno_by_subgenome, bad_effects, candidate_pairs = candidate_pairs),
  "marker"
)

gain_plan <- ng_polyploid_policy(scores, n_crosses = 2L, mode = "gain")
stopifnot(nrow(gain_plan) == 2L)
stopifnot(identical(unique(gain_plan$poly_policy_mode), "gain"))
stopifnot(identical(unique(gain_plan$poly_policy_scope), "wheat_like_allopolyploid_disomic"))
stopifnot(all(diff(gain_plan$poly_usefulness) <= 0))

diversity_plan <- ng_polyploid_policy(scores, n_crosses = 2L, mode = "diversity")
stopifnot(nrow(diversity_plan) == 2L)
stopifnot(identical(unique(diversity_plan$poly_policy_mode), "diversity"))
stopifnot(!is.null(attr(diversity_plan, "summary")))

ocs_plan <- ng_polyploid_policy(scores, n_crosses = 2L, mode = "ocs")
stopifnot(nrow(ocs_plan) == 2L)
stopifnot(identical(unique(ocs_plan$poly_policy_mode), "ocs"))
stopifnot(!is.null(attr(ocs_plan, "summary")))

expect_error(ng_polyploid_policy(scores, n_crosses = 2L, mode = "invalid"), "mode")

cat("polyploid router tests passed\n")
