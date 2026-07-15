ng_poly4x_score_crosses <- function(parent_pop,
                                    sim_param,
                                    n_score_progeny = 25L,
                                    candidate_pairs = NULL,
                                    selection_prop = 0.10,
                                    seed = 1L,
                                    top_prop = 0.10) {
  ng_poly4x_require_alphasimr()
  ng_poly4x_assert_pop4x(parent_pop, "parent population")
  ng_poly4x_validate_parent_ids(parent_pop)

  n_score_progeny <- suppressWarnings(as.numeric(n_score_progeny))
  if (length(n_score_progeny) != 1L || !is.finite(n_score_progeny) ||
      n_score_progeny < 2L || abs(n_score_progeny - round(n_score_progeny)) > 1e-8) {
    ng_stop("n_score_progeny must be an integer >= 2")
  }
  n_score_progeny <- as.integer(n_score_progeny)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    ng_stop("seed must be a single finite integer-like numeric value")
  }
  if (abs(seed - round(seed)) > 1e-8) {
    ng_stop("seed must be a single finite integer-like numeric value")
  }
  seed <- round(seed)
  intensity <- ng_selection_intensity(selection_prop)

  ids <- as.character(parent_pop@id)
  if (is.null(candidate_pairs)) {
    candidate_pairs <- ng_make_pairs(ids, include_self = FALSE)
  }
  candidate_pairs <- ng_poly4x_validate_candidate_pairs(candidate_pairs, ids)

  dosage <- ng_poly4x_pull_dosage(parent_pop, sim_param)
  parent_K <- ng_poly4x_parent_relationship(dosage)
  pair_coancestry <- ng_poly4x_pair_coancestry(parent_K, candidate_pairs)
  digenic <- rowMeans(ng_poly4x_digenic_scaled(dosage))
  scoring_parent_pop <- unserialize(serialize(parent_pop, NULL))
  scoring_sim_param <- if (is.function(sim_param$clone)) sim_param$clone(deep = TRUE) else sim_param
  if ("nThreads" %in% names(scoring_sim_param)) scoring_sim_param$nThreads <- 1L

  scoring_order <- order(ng_poly4x_pair_key(candidate_pairs$parent1, candidate_pairs$parent2))
  stats_by_order <- lapply(scoring_order, function(i) {
    family <- ng_poly4x_make_family(
      parent_pop = scoring_parent_pop,
      parent1 = candidate_pairs$parent1[[i]],
      parent2 = candidate_pairs$parent2[[i]],
      n_progeny = n_score_progeny,
      sim_param = scoring_sim_param,
      seed = ng_poly4x_pair_seed(seed, candidate_pairs$parent1[[i]], candidate_pairs$parent2[[i]])
    )
    ng_poly4x_family_stats(family, top_prop = top_prop)
  })
  stats <- do.call(rbind, stats_by_order)
  stats <- stats[match(seq_len(nrow(candidate_pairs)), scoring_order), , drop = FALSE]

  out <- data.frame(
    parent1 = candidate_pairs$parent1,
    parent2 = candidate_pairs$parent2,
    poly4x_mean = stats$mean_gv,
    poly4x_var = stats$var_gv,
    poly4x_top10 = stats$top10_gv,
    poly4x_max = stats$max_gv,
    poly4x_usefulness = stats$mean_gv + intensity * sqrt(pmax(stats$var_gv, 0)),
    poly4x_pair_coancestry = pair_coancestry,
    pair_kinship = pair_coancestry,
    poly4x_digenic_midparent = (digenic[candidate_pairs$parent1] + digenic[candidate_pairs$parent2]) / 2,
    poly4x_score_n = n_score_progeny,
    stringsAsFactors = FALSE
  )
  attr(out, "parent_K") <- parent_K
  out
}

ng_poly4x_validate_candidate_pairs <- function(candidate_pairs, ids) {
  pairs <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("candidate_pairs must contain parent1 and parent2")
  }
  pairs$parent1 <- as.character(pairs$parent1)
  pairs$parent2 <- as.character(pairs$parent2)
  if (anyNA(pairs$parent1) || anyNA(pairs$parent2) ||
      any(!nzchar(pairs$parent1)) || any(!nzchar(pairs$parent2))) {
    ng_stop("candidate_pairs parent1 and parent2 must be non-missing parent IDs")
  }
  self <- which(pairs$parent1 == pairs$parent2)
  if (length(self)) ng_stop("candidate_pairs must not contain self-crosses")
  missing <- setdiff(unique(c(pairs$parent1, pairs$parent2)), ids)
  if (length(missing)) {
    ng_stop("candidate_pairs contains unknown parent IDs: ", paste(missing, collapse = ", "))
  }
  if (!nrow(pairs)) ng_stop("candidate_pairs must contain at least one cross")
  pair_key <- ng_poly4x_pair_key(pairs$parent1, pairs$parent2)
  repeated <- unique(pair_key[duplicated(pair_key)])
  if (length(repeated)) {
    ng_stop("candidate_pairs contains duplicate or repeated cross entries, including reciprocal duplicates")
  }
  pairs[, c("parent1", "parent2"), drop = FALSE]
}

ng_poly4x_pair_key <- function(parent1, parent2) {
  lo <- pmin(parent1, parent2)
  hi <- pmax(parent1, parent2)
  paste0(nchar(lo, type = "bytes"), ":", lo, "|", nchar(hi, type = "bytes"), ":", hi)
}

ng_poly4x_pair_seed <- function(base_seed, parent1, parent2) {
  if (!is.numeric(base_seed) || length(base_seed) != 1L || !is.finite(base_seed) ||
      abs(base_seed - round(base_seed)) > 1e-8) {
    ng_stop("seed must be a single finite integer-like numeric value")
  }
  modulus <- 2147483646
  key <- paste0(round(base_seed), "||", ng_poly4x_pair_key(parent1, parent2))
  bytes <- as.integer(charToRaw(enc2utf8(key)))
  hash <- 0
  for (byte in bytes) {
    hash <- (hash * 4099 + byte + 1L) %% modulus
  }
  as.integer(hash + 1L)
}

ng_poly4x_select_topn <- function(scores, n_crosses, score_col, method) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!(score_col %in% names(scores))) ng_stop("scores missing score column: ", score_col)
  n_crosses <- suppressWarnings(as.numeric(n_crosses))
  if (length(n_crosses) != 1L || !is.finite(n_crosses) ||
      n_crosses < 1L || abs(n_crosses - round(n_crosses)) > 1e-8) {
    ng_stop("n_crosses must be a positive integer")
  }
  n_crosses <- as.integer(n_crosses)
  score <- as.numeric(scores[[score_col]])
  ok <- is.finite(score)
  if (sum(ok) < n_crosses) ng_stop("Not enough finite scores in ", score_col)
  out <- scores[which(ok)[order(score[ok], decreasing = TRUE)][seq_len(n_crosses)], , drop = FALSE]
  out$poly4x_method <- method
  rownames(out) <- NULL
  out
}

ng_poly4x_var_topn <- function(scores, n_crosses) {
  ng_poly4x_select_topn(scores, n_crosses, "poly4x_var", "ng_poly4x_var_topn")
}

ng_poly4x_usefulness_topn <- function(scores, n_crosses) {
  ng_poly4x_select_topn(scores, n_crosses, "poly4x_usefulness", "ng_poly4x_usefulness_topn")
}

ng_poly4x_ocs <- function(scores,
                          n_crosses,
                          parent_K = NULL,
                          max_crosses_per_parent = 4L,
                          lambda_group = 0.5,
                          lambda_mating = 0,
                          lambda_parent_use = 1.0,
                          lambda_parent_use_mode = "adaptive",
                          method = "auto",
                          local_iter = 1000L,
                          ocs_iter = 5L,
                          ...) {
  if (is.null(parent_K)) parent_K <- attr(scores, "parent_K")
  if (is.null(parent_K)) ng_stop("parent_K is required for ng_poly4x_ocs")
  # Forward the shared mate-selection controls (strategy / diversity_emphasis /
  # target_coancestry / committed_crosses / group_permission / group_quota /
  # cost_col / budget / logistic / lambda_progeny_inbreeding / min_crosses_per_parent
  # and method = "evolution") straight through to the generic allocator via `...`.
  # These are the same, already-tested diploid controls; the polyploid path only
  # differs in the usefulness score fed as gain_col, so no new optimizer logic is needed.
  plan <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = n_crosses,
    gain_col = "poly4x_usefulness",
    parent_K = parent_K,
    max_crosses_per_parent = max_crosses_per_parent,
    lambda_group = lambda_group,
    lambda_mating = lambda_mating,
    lambda_parent_use = lambda_parent_use,
    lambda_parent_use_mode = lambda_parent_use_mode,
    method = method,
    local_iter = local_iter,
    ocs_iter = ocs_iter,
    ...
  )
  plan$poly4x_method <- "ng_poly4x_ocs"
  summary <- attr(plan, "summary")
  summary$poly4x_method <- "ng_poly4x_ocs"
  summary$poly4x_gain_col <- "poly4x_usefulness"
  attr(plan, "summary") <- summary
  plan
}

ng_poly4x_policy_modes <- function() {
  c("gain", "diversity", "ocs")
}

ng_poly4x_normalize_policy_mode <- function(mode) {
  mode <- trimws(tolower(as.character(mode)[[1]]))
  if (!nzchar(mode) || is.na(mode) || !(mode %in% ng_poly4x_policy_modes())) {
    ng_stop("mode must be one of: ", paste(ng_poly4x_policy_modes(), collapse = ", "))
  }
  mode
}

ng_poly4x_policy_reason <- function(mode) {
  switch(
    mode,
    gain = "poly4x_gain_usefulness_validated",
    diversity = "poly4x_strongdiv_diversity_validated",
    ocs = "poly4x_default_ocs_available"
  )
}

ng_poly4x_policy_default_source <- function() {
  "poly4x_diagonal10_20260505"
}

ng_poly4x_policy_auto_source <- function() {
  "poly4x_grid_2rep2cycle_20260506"
}

ng_poly4x_add_policy_diagnostics <- function(plan, mode) {
  plan$poly4x_policy_mode <- mode
  plan$poly4x_policy_reason <- ng_poly4x_policy_reason(mode)
  plan$poly4x_policy_scope <- "autotetraploid_4x"
  plan$poly4x_policy_source <- "poly4x_diagonal10_20260505"

  summary <- attr(plan, "summary")
  if (!is.null(summary)) {
    summary$poly4x_policy_mode <- mode
    summary$poly4x_policy_reason <- ng_poly4x_policy_reason(mode)
    summary$poly4x_policy_scope <- "autotetraploid_4x"
    summary$poly4x_policy_source <- "poly4x_diagonal10_20260505"
    attr(plan, "summary") <- summary
  }
  plan
}

ng_poly4x_policy <- function(scores,
                             n_crosses,
                             mode = "gain",
                             parent_K = NULL,
                             method = "auto",
                             local_iter = 1000L,
                             ocs_iter = 5L,
                             ...) {
  mode <- ng_poly4x_normalize_policy_mode(mode)
  if (identical(mode, "gain")) {
    return(ng_poly4x_add_policy_diagnostics(
      ng_poly4x_usefulness_topn(scores, n_crosses = n_crosses),
      mode
    ))
  }

  if (identical(mode, "diversity")) {
    return(ng_poly4x_add_policy_diagnostics(
      ng_poly4x_ocs(
        scores = scores,
        n_crosses = n_crosses,
        parent_K = parent_K,
        max_crosses_per_parent = 3L,
        lambda_group = 1.0,
        lambda_parent_use = 2.0,
        method = method,
        local_iter = local_iter,
        ocs_iter = ocs_iter,
        ...
      ),
      mode
    ))
  }

  ng_poly4x_add_policy_diagnostics(
    ng_poly4x_ocs(
      scores = scores,
      n_crosses = n_crosses,
      parent_K = parent_K,
      max_crosses_per_parent = 4L,
      lambda_group = 0.5,
      lambda_parent_use = 1.0,
      method = method,
      local_iter = local_iter,
      ocs_iter = ocs_iter,
      ...
    ),
    mode
  )
}

ng_poly4x_policy_select <- function(n_parents,
                                    crop_scenario = "",
                                    crop = "",
                                    mode = "gain",
                                    available_methods = NULL,
                                    source = ng_poly4x_policy_default_source()) {
  n_parents <- suppressWarnings(as.integer(n_parents[[1]]))
  requested_mode <- trimws(tolower(as.character(mode)[[1]]))
  crop_scenario_clean <- trimws(tolower(as.character(crop_scenario)[[1]]))
  crop_clean <- trimws(tolower(as.character(crop)[[1]]))

  candidates <- c(
    "poly4x_gain_policy", "poly4x_diversity_policy", "poly4x_ocs_policy",
    "ng_poly4x_usefulness_topn", "ng_poly4x_ocs", "ng_poly4x_var_topn"
  )

  if (identical(requested_mode, "auto")) {
    source <- ng_poly4x_policy_auto_source()
    if (grepl("potato", crop_scenario_clean, fixed = TRUE) || identical(crop_clean, "potato")) {
      if (is.finite(n_parents) && n_parents <= 20L) {
        mode <- "ocs"
        method <- "poly4x_ocs_policy"
        family <- "poly4x_policy"
        reason <- "potato_autotetraploid_4x_p20_poly4x_policy_ocs_top10"
      } else {
        mode <- "gain"
        method <- "poly4x_gain_policy"
        family <- "poly4x_policy"
        reason <- "potato_autotetraploid_4x_p40_poly4x_policy_gain_top10"
      }
    } else if (grepl("cassava", crop_scenario_clean, fixed = TRUE) || identical(crop_clean, "cassava")) {
      if (is.finite(n_parents) && n_parents <= 20L) {
        mode <- "legacy_usefulness"
        method <- "ng_poly4x_usefulness_topn"
        family <- "poly4x_legacy"
        reason <- "cassava_autotetraploid_4x_p20_ng_poly4x_usefulness_topn_top10"
      } else {
        mode <- "legacy_ocs"
        method <- "ng_poly4x_ocs"
        family <- "poly4x_legacy"
        reason <- "cassava_autotetraploid_4x_p40_ng_poly4x_ocs_top10"
      }
    } else {
      mode <- "gain"
      method <- "poly4x_gain_policy"
      family <- "poly4x_policy"
      reason <- "poly4x_auto_unknown_scenario_gain_fallback"
    }
  } else {
    mode <- ng_poly4x_normalize_policy_mode(mode)
    method <- switch(
      mode,
      gain = "poly4x_gain_policy",
      diversity = "poly4x_diversity_policy",
      ocs = "poly4x_ocs_policy"
    )
    family <- "poly4x_policy"
    reason <- switch(
      mode,
      gain = "poly4x_gain_default",
      diversity = "poly4x_diversity_mode",
      ocs = "poly4x_default_ocs_available"
    )
  }
  policy_candidates <- c("poly4x_gain_policy", "poly4x_diversity_policy", "poly4x_ocs_policy")
  selected <- method
  available_methods <- if (is.null(available_methods)) {
    NULL
  } else {
    unique(trimws(as.character(available_methods)))
  }
  if (!is.null(available_methods) && !(selected %in% available_methods)) {
    hits <- candidates[candidates %in% available_methods]
    if (length(hits)) selected <- hits[[1]]
  }

  list(
    method = selected,
    family = if (selected %in% policy_candidates) "poly4x_policy" else if (grepl("^ng_poly4x_", selected)) "poly4x_legacy" else family,
    mode = mode,
    primary_method = method,
    primary_family = family,
    candidate_methods = candidates,
    source = source,
    reason = reason,
    scope = "autotetraploid_4x",
    band_label = if (is.finite(n_parents)) paste0("parents=", n_parents) else NA_character_,
    is_fallback = !identical(selected, method),
    crop_scenario = trimws(as.character(crop_scenario)[[1]]),
    crop = trimws(as.character(crop)[[1]])
  )
}
