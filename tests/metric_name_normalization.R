ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

stopifnot(ng_normalize_metric_token("mid_parent_mean") == "mean")
stopifnot(ng_normalize_metric_token("family_variance") == "vpm")
stopifnot(ng_normalize_metric_token("reliable_family_variance") == "pmv")
stopifnot(ng_normalize_metric_token("USEFULNESS") == "usefulness")   # case/space tolerant
stopifnot(ng_normalize_metric_token(" pmv ") == "pmv")               # canonical passthrough
stopifnot(ng_normalize_metric_token("var_complex") == "var_complex") # passthrough (build_ctx maps it)

# build_ctx: var_complex -> usefulness + pmv (behavior preserving)
cfg <- list(trait_value_metric = "var_complex", uc_variance_source = "vpm",
            prediction_mode = "trait_by_trait", map_position_unit = "bp",
            threshold_policy = "soft", recomb_model = "haldane", grm_method = "vanraden",
            ld_backend = "auto", method_varPMV = "fast", ril_mode = "infinite",
            posterior_method = "closed_form", duplicate_action = "report",
            lambda_parent_use_mode = "absolute", alphamate_mode = "ModeOptTarget1",
            progeny = "DH", optimizer = "ocs", allocation_method = "ocs",
            run_posterior_prediction = FALSE, use_parallel = FALSE,
            n_iter = 100L, burn_in = 10L, n_crosses = 5L,
            # ng_cp__build_ctx() validates a COMPLETE config -- defaults live in
            # ng_run_cross_prediction()'s formals, not in the builder -- so every key a
            # later release added a validation for must appear here. This one was added
            # after the file was written and left it failing on `ci_level must be one
            # finite probability in (0, 1)` before the metric assertions were reached.
            ci_level = 0.95)
ctx <- ng_cp__build_ctx(cfg)
stopifnot(ctx$trait_value_metric == "usefulness", ctx$uc_variance_source == "pmv")

# friendly metric -> canonical through build_ctx
cfg2 <- modifyList(cfg, list(trait_value_metric = "family_variance", uc_variance_source = "reliable_family_variance"))
ctx2 <- ng_cp__build_ctx(cfg2)
stopifnot(ctx2$trait_value_metric == "vpm", ctx2$uc_variance_source == "pmv")

cat("metric name normalization OK\n")
