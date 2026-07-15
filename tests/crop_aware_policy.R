helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

maize <- ng_crop_genome_select("maize_like")
maize_choice <- ng_crop_aware_policy_select(
  n_parents = 20L,
  crop_scenario = maize$scenario[[1]],
  crop = maize$crop[[1]],
  harness_model = maize$harness_model[[1]],
  validation_scope = maize$validation_scope[[1]]
)
stopifnot(identical(maize_choice$method, "ng_frontier_policy_ocs10_lps2"))
stopifnot(identical(maize_choice$family, "frontier_policy"))
stopifnot(identical(maize_choice$primary_method, "ng_frontier_policy_ocs10_lps2"))
stopifnot(identical(maize_choice$primary_family, "frontier_policy"))
stopifnot(identical(maize_choice$reason, "diploid_frontier_validated"))
stopifnot(!isTRUE(maize_choice$is_fallback))

cassava4x <- ng_crop_genome_select("cassava_tetraploid_stress")
cassava4x_choice <- ng_crop_aware_policy_select(
  n_parents = 20L,
  crop_scenario = cassava4x$scenario[[1]],
  crop = cassava4x$crop[[1]],
  harness_model = cassava4x$harness_model[[1]],
  validation_scope = cassava4x$validation_scope[[1]]
)
stopifnot(identical(cassava4x_choice$method, "alphamate_opt60"))
stopifnot(identical(cassava4x_choice$family, "alphamate"))
stopifnot(identical(cassava4x_choice$primary_method, "alphamate_opt60"))
stopifnot(identical(cassava4x_choice$primary_family, "alphamate"))
stopifnot(identical(cassava4x_choice$reason, "crop_smoke_alpha_winner"))
stopifnot("ng_frontier_policy_ocs10_lps2" %in% cassava4x_choice$candidate_methods)

sugarcane <- ng_crop_genome_select("sugarcane_polyploid_stress")
sugarcane_choice <- ng_crop_aware_policy_select(
  n_parents = 20L,
  crop_scenario = sugarcane$scenario[[1]],
  crop = sugarcane$crop[[1]],
  harness_model = sugarcane$harness_model[[1]],
  validation_scope = sugarcane$validation_scope[[1]]
)
stopifnot(identical(sugarcane_choice$method, "alphamate_opt45"))
stopifnot(identical(sugarcane_choice$primary_method, "alphamate_opt45"))
stopifnot(identical(sugarcane_choice$reason, "crop_smoke_alpha_winner"))

available_fallback <- ng_crop_aware_policy_select(
  n_parents = 20L,
  crop_scenario = sugarcane$scenario[[1]],
  crop = sugarcane$crop[[1]],
  harness_model = sugarcane$harness_model[[1]],
  validation_scope = sugarcane$validation_scope[[1]],
  available_methods = "ng_frontier_policy_ocs10_lps2"
)
stopifnot(identical(available_fallback$method, "ng_frontier_policy_ocs10_lps2"))
stopifnot(isTRUE(available_fallback$is_fallback))

env <- ng_crop_genome_env(sugarcane)
stopifnot(identical(env[["NG_CROP_SCENARIO"]], sugarcane$scenario[[1]]))
stopifnot(identical(env[["NG_CROP"]], sugarcane$crop[[1]]))
stopifnot(identical(env[["NG_CROP_HARNESS_MODEL"]], sugarcane$harness_model[[1]]))
stopifnot(identical(env[["NG_CROP_VALIDATION_SCOPE"]], sugarcane$validation_scope[[1]]))

cat("crop aware policy tests passed\n")
