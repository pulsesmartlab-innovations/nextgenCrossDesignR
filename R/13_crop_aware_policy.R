ng_crop_aware_policy_default_source <- function() {
  "crop_aware_policy_smoke_2026_05_04"
}

ng_crop_aware_policy_candidates <- function(crop_scenario = "",
                                            crop = "",
                                            harness_model = "",
                                            validation_scope = "") {
  crop_scenario <- trimws(as.character(crop_scenario)[[1]])
  crop <- trimws(as.character(crop)[[1]])
  harness_model <- trimws(as.character(harness_model)[[1]])
  validation_scope <- trimws(as.character(validation_scope)[[1]])

  if (identical(crop_scenario, "cassava_tetraploid_stress")) {
    return(list(
      primary_method = "alphamate_opt60",
      candidate_methods = c(
        "alphamate_opt60",
        "ng_frontier_policy_ocs10_lps2",
        "alphamate_opt45",
        "alphamate_opt30"
      ),
      reason = "crop_smoke_alpha_winner"
    ))
  }
  if (identical(crop_scenario, "sugarcane_polyploid_stress")) {
    return(list(
      primary_method = "alphamate_opt45",
      candidate_methods = c(
        "alphamate_opt45",
        "ng_frontier_policy_ocs10_lps2",
        "alphamate_opt30",
        "alphamate_opt60"
      ),
      reason = "crop_smoke_alpha_winner"
    ))
  }

  stress <- grepl("diploidized", harness_model, fixed = TRUE) ||
    grepl("polyploid", validation_scope, ignore.case = TRUE)
  reason <- if (isTRUE(stress)) {
    "diploidized_stress_frontier_smoke"
  } else {
    "diploid_frontier_validated"
  }
  list(
    primary_method = "ng_frontier_policy_ocs10_lps2",
    candidate_methods = c(
      "ng_frontier_policy_ocs10_lps2",
      "alphamate_opt30",
      "alphamate_opt45",
      "alphamate_opt60"
    ),
    reason = reason
  )
}

ng_crop_aware_policy_select <- function(n_parents,
                                        crop_scenario = Sys.getenv("NG_CROP_SCENARIO", unset = ""),
                                        crop = Sys.getenv("NG_CROP", unset = ""),
                                        harness_model = Sys.getenv("NG_CROP_HARNESS_MODEL", unset = ""),
                                        validation_scope = Sys.getenv("NG_CROP_VALIDATION_SCOPE", unset = ""),
                                        available_methods = NULL,
                                        source = ng_crop_aware_policy_default_source()) {
  n_parents <- suppressWarnings(as.integer(n_parents[[1]]))
  policy <- ng_crop_aware_policy_candidates(
    crop_scenario = crop_scenario,
    crop = crop,
    harness_model = harness_model,
    validation_scope = validation_scope
  )
  candidates <- unique(trimws(as.character(policy$candidate_methods)))
  candidates <- candidates[nzchar(candidates) & !is.na(candidates)]
  if (!length(candidates)) ng_stop("crop-aware policy has no candidate methods")

  available_methods <- if (is.null(available_methods)) {
    NULL
  } else {
    unique(trimws(as.character(available_methods)))
  }
  selected <- candidates[[1]]
  if (!is.null(available_methods)) {
    hits <- candidates[candidates %in% available_methods]
    if (length(hits)) selected <- hits[[1]]
  }

  list(
    method = selected,
    family = ng_frontier_method_family(selected),
    primary_method = policy$primary_method,
    primary_family = ng_frontier_method_family(policy$primary_method),
    candidate_methods = candidates,
    source = source,
    reason = policy$reason,
    band_label = if (is.finite(n_parents)) paste0("parents=", n_parents) else NA_character_,
    is_fallback = !identical(selected, policy$primary_method),
    crop_scenario = trimws(as.character(crop_scenario)[[1]]),
    crop = trimws(as.character(crop)[[1]]),
    harness_model = trimws(as.character(harness_model)[[1]]),
    validation_scope = trimws(as.character(validation_scope)[[1]])
  )
}
