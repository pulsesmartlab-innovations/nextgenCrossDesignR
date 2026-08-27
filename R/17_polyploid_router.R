ng_poly_trim1 <- function(x, default = "") {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  trimws(tolower(as.character(x[[1]])))
}

ng_poly_as_logical1 <- function(x, default = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) return(default)
  if (is.logical(x)) return(isTRUE(x[[1]]))
  ng_poly_trim1(x) %in% c("1", "true", "yes", "y")
}

ng_poly_number1 <- function(x, name, default = NULL, integer = FALSE) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) {
    if (!is.null(default)) return(default)
    ng_stop(name, " must be provided")
  }
  value <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(value)) ng_stop(name, " must be a finite numeric value")
  if (isTRUE(integer) && abs(value - round(value)) > 1e-8) {
    ng_stop(name, " must be an integer value")
  }
  if (isTRUE(integer)) as.integer(round(value)) else value
}

ng_poly_clean_subgenomes <- function(subgenome_names) {
  if (is.null(subgenome_names) || !length(subgenome_names)) {
    ng_stop("subgenome_names are required for disomic allopolyploid models")
  }
  out <- trimws(as.character(subgenome_names))
  if (anyNA(out) || any(!nzchar(out))) {
    ng_stop("subgenome_names must be non-missing non-empty values")
  }
  if (anyDuplicated(out)) ng_stop("subgenome_names must be unique")
  out
}

ng_poly_clean_ploidy_profile <- function(subgenome_names, ploidy_profile = NULL) {
  subgenome_names <- ng_poly_clean_subgenomes(subgenome_names)
  if (is.null(ploidy_profile)) {
    out <- rep.int(2L, length(subgenome_names))
    names(out) <- subgenome_names
    return(out)
  }
  raw <- suppressWarnings(as.numeric(ploidy_profile))
  if (length(raw) != length(subgenome_names) || anyNA(raw) || any(!is.finite(raw))) {
    ng_stop("ploidy_profile must provide one finite numeric value per subgenome")
  }
  if (any(abs(raw - round(raw)) > 1e-8)) {
    ng_stop("ploidy_profile must provide integer values per subgenome")
  }
  out <- as.integer(round(raw))
  profile_names <- names(ploidy_profile)
  if (is.null(profile_names) || any(!nzchar(profile_names))) profile_names <- subgenome_names
  names(out) <- profile_names
  if (!setequal(names(out), subgenome_names)) {
    ng_stop("ploidy_profile names must match subgenome_names")
  }
  out <- out[subgenome_names]
  if (any(out != 2L)) {
    ng_stop("Phase 2A allopolyploid subgenomes must be diploid with ploidy_profile values of 2")
  }
  out
}

ng_poly_model_decision <- function(model_family, inheritance_model, supported,
                                   supported_scope, reason, crop = "",
                                   crop_scenario = "", ploidy = NA_integer_,
                                   ploidy_profile = NULL,
                                   subgenome_names = character()) {
  list(
    model_family = model_family,
    inheritance_model = inheritance_model,
    supported = isTRUE(supported),
    supported_scope = supported_scope,
    reason = reason,
    crop = crop,
    crop_scenario = crop_scenario,
    ploidy = ploidy,
    ploidy_profile = ploidy_profile,
    subgenome_names = subgenome_names
  )
}

ng_polyploid_model_select <- function(crop = "",
                                 crop_scenario = "",
                                 ploidy = NULL,
                                 inheritance_model = "",
                                 subgenome_names = NULL,
                                 ploidy_profile = NULL,
                                 quad_prob = 0,
                                 empirical_generator = FALSE) {
  crop <- ng_poly_trim1(crop)
  crop_scenario <- ng_poly_trim1(crop_scenario)
  inheritance_model <- ng_poly_trim1(inheritance_model)
  empirical_generator <- ng_poly_as_logical1(empirical_generator)
  quad_prob <- ng_poly_number1(quad_prob, "quad_prob", default = 0)

  crop_wheat_like <- crop %in% c("wheat", "bread_wheat") ||
    grepl("wheat", crop_scenario, fixed = TRUE)
  crop_sugarcane_like <- crop %in% c("sugarcane") ||
    grepl("sugarcane", crop_scenario, fixed = TRUE)
  allopolyploid_inheritance <- c("", "disomic_subgenome", "allopolyploid_subgenome", "allopolyploid")

  if (crop_wheat_like && !(inheritance_model %in% allopolyploid_inheritance)) {
    ng_stop("inheritance_model is unsupported for wheat-like disomic subgenome models")
  }

  if (crop_sugarcane_like || inheritance_model %in% c("mixed_or_aneuploid", "complex_polyploid")) {
    reason <- if (empirical_generator) {
      "complex_polyploid_empirical_generator_contract_not_implemented_in_phase2a"
    } else {
      "complex_polyploid_not_supported_without_empirical_generator"
    }
    return(ng_poly_model_decision(
      model_family = "complex_polyploid_guard",
      inheritance_model = "mixed_or_aneuploid",
      supported = FALSE,
      supported_scope = "requires_empirical_or_external_simulator",
      reason = reason,
      crop = crop,
      crop_scenario = crop_scenario
    ))
  }

  if (crop_wheat_like || inheritance_model %in% allopolyploid_inheritance[-1L]) {
    if (quad_prob > 0) {
      ng_stop("quadrivalent pairing is not supported for Phase 2A disomic subgenome models")
    }
    subgenome_names <- ng_poly_clean_subgenomes(subgenome_names)
    ploidy_profile <- ng_poly_clean_ploidy_profile(subgenome_names, ploidy_profile)
    return(ng_poly_model_decision(
      model_family = "allopolyploid_subgenome",
      inheritance_model = "disomic_subgenome",
      supported = TRUE,
      supported_scope = "wheat_like_allopolyploid_disomic",
      reason = "wheat_like_disomic_subgenome_phase2a",
      crop = crop,
      crop_scenario = crop_scenario,
      ploidy = sum(ploidy_profile),
      ploidy_profile = ploidy_profile,
      subgenome_names = subgenome_names
    ))
  }

  ploidy_value <- if (is.null(ploidy)) NA_integer_ else ng_poly_number1(ploidy, "ploidy", integer = TRUE)
  if (identical(inheritance_model, "polysomic") && is.finite(ploidy_value) && ploidy_value == 4L) {
    return(ng_poly_model_decision(
      model_family = "autotetraploid_4x",
      inheritance_model = "polysomic",
      supported = TRUE,
      supported_scope = "autotetraploid_4x",
      reason = "existing_ng_poly4x_model",
      crop = crop,
      crop_scenario = crop_scenario,
      ploidy = 4L
    ))
  }

  if (identical(inheritance_model, "polysomic") && is.finite(ploidy_value) &&
      ploidy_value > 4L && ploidy_value %% 2L == 0L) {
    return(ng_poly_model_decision(
      model_family = "autopolyploid_even_ploidy",
      inheritance_model = "polysomic",
      supported = FALSE,
      supported_scope = "requires_even_ploidy_autopolyploid_validation",
      reason = "autopolyploid_even_ploidy_not_validated",
      crop = crop,
      crop_scenario = crop_scenario,
      ploidy = ploidy_value
    ))
  }

  ng_stop("inheritance_model is unsupported or insufficiently declared for polyploid model selection")
}
