helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
  message("AlphaSimR unavailable; skipping poly4x simulation test")
  quit(save = "no", status = 0)
}

scenarios <- ng_poly4x_scenarios()
stopifnot(all(c("potato_autotetraploid_4x", "cassava_autotetraploid_4x") %in% scenarios$scenario))

scenario <- ng_poly4x_select_scenario("potato_autotetraploid_4x")
scenario$n_founders <- 12L
scenario$n_chr <- 2L
scenario$seg_sites <- 24L
scenario$snp_per_chr <- 6L
scenario$qtl_per_chr <- 4L
scenario$genome_length_m <- 0.5

seed_err <- tryCatch(
  ng_poly4x_setup_simparam(scenario, seed = c(1L, 2L)),
  error = function(e) e
)
stopifnot(inherits(seed_err, "error"))
stopifnot(grepl("seed", conditionMessage(seed_err), fixed = TRUE))

bad_scenario <- scenario
bad_scenario$n_chr <- 0L
scenario_err <- tryCatch(
  ng_poly4x_setup_simparam(bad_scenario, seed = 17L),
  error = function(e) e
)
stopifnot(inherits(scenario_err, "error"))
stopifnot(grepl("n_chr", conditionMessage(scenario_err), fixed = TRUE))

setup <- ng_poly4x_setup_simparam(scenario, seed = 17L, include_digenic = TRUE)
parents <- ng_poly4x_make_parent_pop(setup, n_parents = 8L)
stopifnot(AlphaSimR::nInd(parents) == 8L)
stopifnot(identical(as.integer(unique(parents@ploidy)), 4L))

dosage <- ng_poly4x_pull_dosage(parents, setup$sim_param)
stopifnot(nrow(dosage) == 8L)
stopifnot(ncol(dosage) == scenario$n_chr[[1]] * scenario$snp_per_chr[[1]])
stopifnot(all(dosage >= 0 & dosage <= 4))
stopifnot(all(abs(dosage - round(dosage)) < 1e-8))

family <- ng_poly4x_make_family(
  parent_pop = parents,
  parent1 = parents@id[[1]],
  parent2 = parents@id[[2]],
  n_progeny = 5L,
  sim_param = setup$sim_param,
  seed = 101L
)
stopifnot(AlphaSimR::nInd(family) == 5L)
stopifnot(identical(as.integer(unique(family@ploidy)), 4L))
family_dosage <- ng_poly4x_pull_dosage(family, setup$sim_param)
stopifnot(all(family_dosage >= 0 & family_dosage <= 4))

self_err <- tryCatch(
  ng_poly4x_make_family(parents, parents@id[[1]], parents@id[[1]], 2L, setup$sim_param),
  error = function(e) e
)
stopifnot(inherits(self_err, "error"))
stopifnot(grepl("self-cross", conditionMessage(self_err), fixed = TRUE))

duplicate_parents <- parents
duplicate_parents@id[[2]] <- duplicate_parents@id[[1]]
duplicate_err <- tryCatch(
  ng_poly4x_make_family(
    duplicate_parents,
    duplicate_parents@id[[1]],
    duplicate_parents@id[[3]],
    2L,
    setup$sim_param
  ),
  error = function(e) e
)
stopifnot(inherits(duplicate_err, "error"))
stopifnot(grepl("unique|duplicate", conditionMessage(duplicate_err), ignore.case = TRUE))

stats <- ng_poly4x_family_stats(family)
stopifnot(all(c("mean_gv", "var_gv", "top10_gv", "max_gv") %in% names(stats)))
stopifnot(all(is.finite(unlist(stats))))

top_prop_err <- tryCatch(
  ng_poly4x_family_stats(family, top_prop = 2),
  error = function(e) e
)
stopifnot(inherits(top_prop_err, "error"))
stopifnot(grepl("top_prop", conditionMessage(top_prop_err), fixed = TRUE))

cat("poly4x simulation tests passed\n")
