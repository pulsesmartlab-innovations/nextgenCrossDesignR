helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

expected_default <- c(
  `20` = "ng_recomb_gebv_ocs10_lps2",
  `30` = "ng_meta_selector_ocs10_lps2",
  `40` = "ng_pmv_blend_balanced_ocs10_lps2",
  `50` = "ng_meta_router_ocs10_lps2",
  `60` = "popvar_uc_ocs10_lps1",
  `70` = "ng_meta_selector_ocs10_lps2",
  `80` = "ng_meta_portfolio_ocs10_lps2"
)

for (n in names(expected_default)) {
  selected <- ng_frontier_policy_select(as.integer(n))
  stopifnot(identical(selected$method, unname(expected_default[[n]])))
  stopifnot(identical(selected$primary_method, unname(expected_default[[n]])))
  stopifnot(identical(selected$source, "validated_5k_parent_grid_2026_05_03"))
  stopifnot(!isTRUE(selected$is_fallback))
}

external_fallback <- ng_frontier_policy_select(
  60L,
  available_methods = c("simple_usefa_ocs10_lps1", "ng_meta_router_ocs10_lps2")
)
stopifnot(identical(external_fallback$method, "simple_usefa_ocs10_lps1"))
stopifnot(isTRUE(external_fallback$is_fallback))
stopifnot(identical(external_fallback$primary_method, "popvar_uc_ocs10_lps1"))

internal_fallback <- ng_frontier_policy_select(
  60L,
  available_methods = "ng_meta_router_ocs10_lps2"
)
stopifnot(identical(internal_fallback$method, "ng_meta_router_ocs10_lps2"))
stopifnot(isTRUE(internal_fallback$is_fallback))

custom_policy <- ng_frontier_policy_from_spec(
  "35-45=ng_meta_router_ocs10_lps2|var_simple_ocs10_lps2;60=ng_meta_portfolio_ocs10_lps2",
  source = "site_validation_v1"
)
custom <- ng_frontier_policy_select(40L, policy = custom_policy)
stopifnot(identical(custom$method, "ng_meta_router_ocs10_lps2"))
stopifnot(identical(custom$source, "site_validation_v1"))
stopifnot(identical(custom$band_label, "35-45"))

external_need <- ng_external_methods_needed("ng_frontier_policy_ocs10_lps2")
stopifnot(isTRUE(external_need$popvar))
stopifnot(isTRUE(external_need$simplemating_usefa))

cat("frontier policy tests passed\n")
