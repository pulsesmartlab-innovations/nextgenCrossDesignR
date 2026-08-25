ng_design_crosses <- function(geno,
                              y,
                              marker_map = NULL,
                              ids = rownames(geno),
                              adjusted_pheno = NULL,
                              blue = NULL,
                              blup = NULL,
                              n_crosses = 100,
                              max_crosses_per_parent = 6,
                              target = c("DH", "RIL"),
                              parent_type = c("inbred", "dh", "ril"),
                              phased_haplotypes = NULL,
                              recomb_model = c("haldane", "kosambi"),
                              selection_prop = 0.10,
                              min_effect_reliability = 0.35,
                              window_cm = Inf,
                              gain_col = "usefulness_pmv_gebv",
                              method = "auto",
                              lambda_group = 0,
                              lambda_mating = 0,
                              marker_target_spec = NULL,
                              lethal_spec = NULL,
                              marker_ploidy = 2,
                              lambda_marker = 0,
                              drop_lethal_carrier_crosses = TRUE,
                              ld_pruning = FALSE,
                              ld_window = 100L,
                              ld_r2_threshold = 0.9,
                              ld_maf_threshold = 0.01,
                              ld_ploidy = 2,
                              ld_backend = c("auto", "cpp", "r"),
                              use_cpp = TRUE,
                              assume_inbred = NULL,
                              grm_method = c("vanraden", "yang"),
                              seed = 1L,
                              ...) {
  n_crosses_was_missing <- missing(n_crosses)
  ld_backend <- match.arg(ld_backend)
  grm_method <- match.arg(grm_method)
  target <- match.arg(target)
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  recomb_model <- match.arg(recomb_model)
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  ld_pruning_report <- NULL
  if (isTRUE(ld_pruning)) {
    ld_marker_map <- ng_prepare_marker_map(marker_map, colnames(geno), model = recomb_model)
    geno <- ng_ld_prune_geno(
      geno = geno,
      window = ld_window,
      r2_threshold = ld_r2_threshold,
      maf_threshold = ld_maf_threshold,
      ploidy = ld_ploidy,
      backend = ld_backend,
      marker_map = ld_marker_map
    )
    ld_pruning_report <- attr(geno, "ld_pruning_report", exact = TRUE)
  }
  effects <- ng_fit_ridge_effects(geno, y, ids = ids, seed = seed)
  scores <- ng_score_crosses(
    geno = geno,
    effects = effects,
    marker_map = marker_map,
    ids = ids,
    adjusted_pheno = adjusted_pheno,
    blue = blue,
    blup = blup,
    target = target,
    parent_type = parent_type,
    phased_haplotypes = phased_haplotypes,
    recomb_model = recomb_model,
    selection_prop = selection_prop,
    min_effect_reliability = min_effect_reliability,
    window_cm = window_cm,
    use_cpp = use_cpp,
    grm_method = grm_method
  )
  # Marker steering / lethal-allele guarding at the scoring layer (Module 4): attach the
  # per-cross columns, optionally drop carrier x carrier crosses, and blend the steering
  # score into marker_adjusted_gain. When a marker-target emphasis is requested we
  # optimize on that blended column so marker steering flows into the allocation.
  if (!is.null(marker_target_spec) || !is.null(lethal_spec)) {
    scores <- ng_apply_marker_management(
      scores = scores, geno = geno,
      marker_target_spec = marker_target_spec, lethal_spec = lethal_spec,
      ploidy = marker_ploidy, gain_col = gain_col, lambda_marker = lambda_marker,
      drop_lethal_carrier_crosses = drop_lethal_carrier_crosses
    )
    if (is.finite(lambda_marker) && lambda_marker != 0 &&
        "marker_adjusted_gain" %in% names(scores)) {
      gain_col <- "marker_adjusted_gain"
    }
  }
  n_crosses_num <- suppressWarnings(as.numeric(n_crosses))
  if (length(n_crosses_num) != 1L || !is.finite(n_crosses_num) || n_crosses_num < 1 ||
      abs(n_crosses_num - round(n_crosses_num)) > 1e-8) {
    ng_stop("n_crosses must be a positive integer")
  }
  n_crosses <- as.integer(round(n_crosses_num))
  if (isTRUE(n_crosses_was_missing)) {
    n_crosses <- min(n_crosses, max(1L, nrow(scores)))
  }
  parent_kinship <- ng_parent_kinship(geno, method = grm_method)
  # gain_col, method, and ... are forwarded so the documented pipeline can reach every
  # allocator (including method = "evolution") and every optimizer option added in the
  # mate-selection modules: strategy / diversity_emphasis, lambda_progeny_inbreeding,
  # min_crosses_per_parent, committed_crosses, parent_group / group_permission /
  # group_quota, cost_col / budget / lambda_cost, logistic_col / lambda_logistic, etc.
  plan <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = n_crosses,
    gain_col = gain_col,
    parent_kinship = parent_kinship,
    max_crosses_per_parent = max_crosses_per_parent,
    lambda_group = lambda_group,
    lambda_mating = lambda_mating,
    method = method,
    ...
  )
  list(
    effects = effects,
    scores = scores,
    plan = plan,
    plan_summary = attr(plan, "summary"),
    ld_pruning_report = ld_pruning_report
  )
}
