ng_multitrait_default_source <- function() {
  "rank_normalized_ocs_framework_2026_05_06"
}

ng_multitrait_method_family <- function(method) {
  method <- trimws(tolower(as.character(method)[[1]]))
  if (identical(method, "auto")) return("rank_threshold")
  if (identical(method, "threshold")) return("threshold")
  if (identical(method, "weighted")) return("weighted_index")
  if (identical(method, "economic_index")) return("economic_index")
  if (identical(method, "desired_gain")) return("desired_gain")
  "other"
}

ng_multitrait_policy_select <- function(method = "auto",
                                        available_methods = NULL,
                                        source = ng_multitrait_default_source()) {
  method <- trimws(tolower(as.character(method)[[1]]))
  candidates <- c("auto", "threshold", "weighted", "economic_index", "desired_gain")
  if (!(method %in% candidates)) ng_stop("method must be one of: ", paste(candidates, collapse = ", "))
  candidates <- unique(c(method, candidates))
  available_methods <- if (is.null(available_methods)) {
    NULL
  } else {
    unique(trimws(tolower(as.character(available_methods))))
  }
  selected <- candidates[[1]]
  if (!is.null(available_methods)) {
    hits <- candidates[candidates %in% available_methods]
    if (length(hits)) selected <- hits[[1]]
  }
  reason <- if (identical(method, "auto")) {
    "rank_normalized_unknown_weights_default"
  } else {
    paste0("requested_", method)
  }
  list(
    method = selected,
    family = ng_multitrait_method_family(selected),
    primary_method = method,
    primary_family = ng_multitrait_method_family(method),
    candidate_methods = candidates,
    source = source,
    reason = reason,
    is_fallback = !identical(selected, method)
  )
}

ng_multitrait_recycle <- function(x, n, name, default = NA) {
  if (is.null(x)) {
    if (length(default) == 1L) return(rep(default, n))
    if (length(default) == n) return(default)
    ng_stop(name, " default must have length 1 or length ", n)
  }
  if (length(x) == 1L) return(rep(x, n))
  if (length(x) != n) ng_stop(name, " must have length 1 or length ", n)
  x
}

ng_multitrait_direction <- function(direction) {
  direction <- trimws(tolower(as.character(direction)))
  out <- rep(NA_character_, length(direction))
  out[direction %in% c("max", "maximize", "maximise", "increase", "higher", "high", "positive", "+")] <- "maximize"
  out[direction %in% c("min", "minimize", "minimise", "decrease", "lower", "low", "negative", "-")] <- "minimize"
  if (anyNA(out)) {
    bad <- unique(direction[is.na(out)])
    ng_stop("invalid multi-trait direction: ", paste(bad, collapse = ", "))
  }
  out
}

ng_multitrait_spec <- function(trait,
                               column = NULL,
                               direction = NULL,
                               weight = NULL,
                               min_value = NULL,
                               max_value = NULL,
                               threshold_weight = NULL,
                               desired_change = NULL,
                               economic_weight = NULL) {
  if (is.data.frame(trait)) {
    spec <- as.data.frame(trait, stringsAsFactors = FALSE)
    if (!("column" %in% names(spec)) && "trait" %in% names(spec)) spec$column <- spec$trait
    if (!("trait" %in% names(spec)) && "column" %in% names(spec)) spec$trait <- spec$column
    if (!("direction" %in% names(spec))) spec$direction <- "maximize"
    if (!("weight" %in% names(spec))) spec$weight <- NA_real_
    if (!("min_value" %in% names(spec))) spec$min_value <- NA_real_
    if (!("max_value" %in% names(spec))) spec$max_value <- NA_real_
    if (!("threshold_weight" %in% names(spec))) spec$threshold_weight <- 1
    if (!("desired_change" %in% names(spec))) spec$desired_change <- NA_real_
    if (!("economic_weight" %in% names(spec))) spec$economic_weight <- NA_real_
  } else {
    trait <- as.character(trait)
    n <- length(trait)
    if (!n) ng_stop("trait spec must contain at least one trait")
    column <- ng_multitrait_recycle(column, n, "column", trait)
    direction <- ng_multitrait_recycle(direction, n, "direction", "maximize")
    weight <- ng_multitrait_recycle(weight, n, "weight", NA_real_)
    min_value <- ng_multitrait_recycle(min_value, n, "min_value", NA_real_)
    max_value <- ng_multitrait_recycle(max_value, n, "max_value", NA_real_)
    threshold_weight <- ng_multitrait_recycle(threshold_weight, n, "threshold_weight", 1)
    desired_change <- ng_multitrait_recycle(desired_change, n, "desired_change", NA_real_)
    economic_weight <- ng_multitrait_recycle(economic_weight, n, "economic_weight", NA_real_)
    spec <- data.frame(
      trait = trait,
      column = as.character(column),
      direction = as.character(direction),
      weight = weight,
      min_value = min_value,
      max_value = max_value,
      threshold_weight = threshold_weight,
      desired_change = desired_change,
      economic_weight = economic_weight,
      stringsAsFactors = FALSE
    )
  }

  required <- c("trait", "column", "direction", "weight", "min_value", "max_value",
                "threshold_weight", "desired_change", "economic_weight")
  missing <- setdiff(required, names(spec))
  if (length(missing)) ng_stop("multi-trait spec missing columns: ", paste(missing, collapse = ", "))
  spec <- spec[required]
  spec$trait <- trimws(as.character(spec$trait))
  spec$column <- trimws(as.character(spec$column))
  if (any(!nzchar(spec$trait) | is.na(spec$trait))) ng_stop("trait names must be non-empty")
  if (any(!nzchar(spec$column) | is.na(spec$column))) ng_stop("trait columns must be non-empty")
  if (any(duplicated(spec$trait))) ng_stop("trait names must be unique")
  spec$direction <- ng_multitrait_direction(spec$direction)
  spec$weight <- suppressWarnings(as.numeric(spec$weight))
  spec$min_value <- suppressWarnings(as.numeric(spec$min_value))
  spec$max_value <- suppressWarnings(as.numeric(spec$max_value))
  spec$threshold_weight <- suppressWarnings(as.numeric(spec$threshold_weight))
  spec$desired_change <- suppressWarnings(as.numeric(spec$desired_change))
  spec$economic_weight <- suppressWarnings(as.numeric(spec$economic_weight))
  spec$threshold_weight[!is.finite(spec$threshold_weight) | spec$threshold_weight < 0] <- 1
  spec$economic_weight[is.finite(spec$economic_weight) & spec$economic_weight < 0] <- NA_real_
  both <- is.finite(spec$min_value) & is.finite(spec$max_value)
  if (any(both & spec$min_value > spec$max_value)) ng_stop("min_value cannot exceed max_value")
  rownames(spec) <- NULL
  spec
}

ng_breeder_selection_objective <- function(trait,
                                           column = NULL,
                                           direction = NULL,
                                           weight = NULL,
                                           min_value = NULL,
                                           max_value = NULL,
                                           threshold_weight = NULL,
                                           desired_change = NULL,
                                           economic_weight = NULL,
                                           method = "auto",
                                           threshold_policy = c("soft", "strict"),
                                           strict_direction_required = TRUE,
                                           source = "breeder_selection_objective_2026_06_20") {
  threshold_policy <- match.arg(threshold_policy)
  method <- trimws(tolower(as.character(method)[[1]]))
  allowed <- c("auto", "weighted", "economic_index", "desired_gain")
  if (!(method %in% allowed)) ng_stop("method must be one of: ", paste(allowed, collapse = ", "))

  direction_missing <- FALSE
  if (is.data.frame(trait)) {
    direction_missing <- !("direction" %in% names(trait)) ||
      any(!nzchar(trimws(as.character(trait$direction))) | is.na(trait$direction))
  } else {
    direction_missing <- is.null(direction) ||
      any(!nzchar(trimws(as.character(direction))) | is.na(direction))
  }
  if (isTRUE(strict_direction_required) && isTRUE(direction_missing)) {
    ng_stop("direction is required for breeder selection objectives; ",
            "set direction to maximize or minimize for every trait")
  }

  traits <- ng_multitrait_spec(
    trait = trait,
    column = column,
    direction = direction,
    weight = weight,
    min_value = min_value,
    max_value = max_value,
    threshold_weight = threshold_weight,
    desired_change = desired_change,
    economic_weight = economic_weight
  )

  method_requested <- method
  finite_positive <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    any(is.finite(x) & x > 0)
  }
  has_desired <- finite_positive(traits$desired_change)
  has_economic <- finite_positive(traits$economic_weight)
  has_weight <- finite_positive(traits$weight)
  method_reason <- paste0("requested_", method_requested)
  if (identical(method_requested, "auto")) {
    if (has_desired && (has_economic || has_weight)) {
      method <- "desired_gain"
      method_reason <- "auto_promoted_desired_gain"
    } else if (has_economic) {
      method <- "economic_index"
      method_reason <- "auto_promoted_economic_index"
    } else if (has_weight) {
      method <- "weighted"
      method_reason <- "auto_promoted_weighted"
    } else {
      method <- "auto"
      method_reason <- "equal_weight_rank_default"
    }
  }

  weights <- ng_multitrait_resolve_weights(traits, method)
  has_thresholds <- any(is.finite(traits$min_value) | is.finite(traits$max_value))
  out <- list(
    traits = traits,
    method = method,
    method_requested = method_requested,
    strict_direction_required = isTRUE(strict_direction_required),
    threshold_policy = threshold_policy,
    strict_thresholds = identical(threshold_policy, "strict"),
    source = source,
    diagnostics = list(
      method_reason = method_reason,
      resolved_weights = weights,
      has_weights = has_weight,
      has_economic_weights = has_economic,
      has_desired_change = has_desired,
      has_thresholds = has_thresholds,
      threshold_policy = threshold_policy
    )
  )
  class(out) <- c("ng_breeder_selection_objective", "list")
  out
}

ng_multitrait_resolve_weights <- function(spec, method) {
  method <- trimws(tolower(as.character(method)[[1]]))
  raw <- suppressWarnings(as.numeric(spec$weight))
  if (identical(method, "weighted")) {
    if (!any(is.finite(raw) & raw > 0)) raw <- suppressWarnings(as.numeric(spec$economic_weight))
    if (any(!is.finite(raw) | raw <= 0)) ng_stop("weighted multi-trait scoring requires positive finite weights")
  } else if (identical(method, "economic_index")) {
    raw <- suppressWarnings(as.numeric(spec$economic_weight))
    if (!any(is.finite(raw) & raw > 0)) raw <- suppressWarnings(as.numeric(spec$weight))
  } else if (identical(method, "desired_gain")) {
    raw <- suppressWarnings(as.numeric(spec$economic_weight))
    if (!any(is.finite(raw) & raw > 0)) raw <- suppressWarnings(as.numeric(spec$weight))
    if (!any(is.finite(raw) & raw > 0)) raw <- abs(suppressWarnings(as.numeric(spec$desired_change)))
  }

  positive <- is.finite(raw) & raw > 0
  if (!any(positive)) {
    raw <- rep(1, nrow(spec))
  } else if (any(!positive)) {
    raw[!positive] <- stats::median(raw[positive], na.rm = TRUE)
  }
  raw[!is.finite(raw) | raw <= 0] <- 1
  weights <- raw / sum(raw)
  names(weights) <- spec$trait
  weights
}

ng_multitrait_economic_weights <- function(spec) {
  raw <- suppressWarnings(as.numeric(spec$economic_weight))
  if (!any(is.finite(raw) & raw > 0)) raw <- suppressWarnings(as.numeric(spec$weight))
  positive <- is.finite(raw) & raw > 0
  if (!any(positive)) {
    raw <- rep(1, nrow(spec))
  } else if (any(!positive)) {
    raw[!positive] <- stats::median(raw[positive], na.rm = TRUE)
  }
  raw[!is.finite(raw) | raw <= 0] <- 1
  out <- raw / sum(raw)
  names(out) <- spec$trait
  out
}

ng_multitrait_value_scale <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(1)
  scale <- stats::IQR(x, na.rm = TRUE) / 1.349
  if (!is.finite(scale) || scale <= 0) scale <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- diff(range(x, na.rm = TRUE))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  scale
}

ng_multitrait_diag_col <- function(trait, suffix) {
  paste0("multi_trait_", make.names(as.character(trait)), "_", suffix)
}

# Build the covariance matrix used to solve the index. Sources, in order of
# preference:
#   1. `genetic_covariance` (G) supplied by the caller (Smith-Hazel / Pesek-Baker
#      ideal). Use directly if `index_basis = "genetic"` (default when supplied).
#   2. `phenotypic_covariance` (P) supplied by the caller, used as Smith-Hazel
#      P alongside `genetic_covariance` (b = P^{-1} G a).
#   3. Sample Pearson covariance of the oriented candidate-cross trait predictions
#      (legacy "phenotypic proxy"). This is NOT a population-genetic parameter
#      and should be treated as an empirical approximation only.
ng_multitrait_index_covariance <- function(value_z, traits,
                                           phenotypic_covariance = NULL,
                                           genetic_covariance = NULL) {
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  ensure_pxp <- function(M) {
    if (is.null(M)) return(NULL)
    M <- as.matrix(M)
    if (nrow(M) != p || ncol(M) != p) return(NULL)
    M[!is.finite(M)] <- 0
    M
  }
  G <- ensure_pxp(genetic_covariance)
  P <- ensure_pxp(phenotypic_covariance)
  if (!is.null(G) && !is.null(P)) {
    return(list(target_matrix = P, projection = G,
                source = "smith_hazel", solve_form = "b = P^{-1} G a"))
  }
  if (!is.null(G)) {
    return(list(target_matrix = G, projection = NULL,
                source = "genetic_covariance", solve_form = "b = G^{-1} a"))
  }
  if (!is.null(P)) {
    return(list(target_matrix = P, projection = NULL,
                source = "phenotypic_covariance", solve_form = "b = P^{-1} a"))
  }
  complete <- stats::complete.cases(value_z)
  cov_mat <- if (sum(complete) > 1L) {
    stats::cov(value_z[complete, , drop = FALSE])
  } else {
    diag(p)
  }
  cov_mat <- as.matrix(cov_mat)
  cov_mat[!is.finite(cov_mat)] <- 0
  if (nrow(cov_mat) != p || ncol(cov_mat) != p) cov_mat <- diag(p)
  list(target_matrix = cov_mat, projection = NULL,
       source = "phenotypic_proxy", solve_form = "b = Sigma_hat^{-1} a")
}

# Single-ridge solve: solve_mat = target_matrix + ridge * diag_scale * I.
# Previously the diagonal was pre-inflated AND a ridge term was added, applying
# the penalty twice. Returns the coefficients and the predicted response under
# the index theory:
#   - phenotypic_proxy: predicted = Sigma_hat %*% b (sample covariance proxy)
#   - genetic_covariance: predicted = G %*% b (Pesek-Baker realized response)
#   - smith_hazel: predicted = G %*% b (Smith-Hazel realized response)
ng_multitrait_solve_index <- function(target, cov_info, ridge = 1e-6) {
  p <- length(target)
  M <- cov_info$target_matrix
  diag_scale <- mean(diag(M), na.rm = TRUE)
  if (!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
  solve_mat <- M + diag(ridge * diag_scale, p)
  # Smith-Hazel index: b = P^{-1} G a. When a projection (G) is supplied alongside
  # the target matrix (P), the economic weights `target` (a) must be projected
  # through G before the P-solve. Without this, the solve degenerates to
  # b = P^{-1} a and G never enters the scoring coefficients (it would only appear
  # in the reported `predicted` diagnostic). For the single-matrix branches
  # (projection = NULL) rhs stays `target`, preserving b = G^{-1} a / P^{-1} a /
  # Sigma^{-1} a exactly as before.
  rhs <- if (!is.null(cov_info$projection)) as.numeric(cov_info$projection %*% target) else target
  coefficients <- tryCatch(
    as.numeric(qr.solve(solve_mat, rhs)),
    error = function(e) target
  )
  coefficients[!is.finite(coefficients)] <- target[!is.finite(coefficients)]
  if (!any(is.finite(coefficients) & abs(coefficients) > 0)) coefficients <- target
  coefficients <- coefficients / sum(abs(coefficients))
  proj <- if (!is.null(cov_info$projection)) cov_info$projection else cov_info$target_matrix
  predicted <- as.numeric(proj %*% coefficients)
  predicted[!is.finite(predicted)] <- 0
  list(coefficients = coefficients, predicted = predicted)
}

ng_multitrait_desired_gain_fit <- function(value_z, traits, value_scales, ridge = 1e-6,
                                           phenotypic_covariance = NULL,
                                           genetic_covariance = NULL) {
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  if (!p) ng_stop("desired-gain optimizer needs at least one trait")
  economic <- ng_multitrait_economic_weights(traits)
  desired <- abs(suppressWarnings(as.numeric(traits$desired_change)))
  desired[!is.finite(desired) | desired <= 0] <- value_scales[!is.finite(desired) | desired <= 0]
  desired[!is.finite(desired) | desired <= 0] <- 1
  desired_sd <- desired / pmax(value_scales, 1e-8)
  target <- desired_sd * economic
  if (!any(is.finite(target) & target > 0)) target <- rep(1, p)
  target[!is.finite(target) | target < 0] <- 0
  target <- target / sum(target)
  names(target) <- traits$trait

  # External covariances arrive in raw units; map G/P into the value_z space. The desired-gain
  # target is already in value_z scale (desired_sd = |desired| / scale), so it needs no rescaling.
  if (!is.null(phenotypic_covariance) || !is.null(genetic_covariance)) {
    phenotypic_covariance <- ng_multitrait_cov_to_value_z(phenotypic_covariance, traits, value_scales)
    genetic_covariance <- ng_multitrait_cov_to_value_z(genetic_covariance, traits, value_scales)
  }

  cov_info <- ng_multitrait_index_covariance(value_z, traits,
                                             phenotypic_covariance = phenotypic_covariance,
                                             genetic_covariance = genetic_covariance)
  solved <- ng_multitrait_solve_index(target, cov_info, ridge = ridge)
  coefficients <- solved$coefficients
  predicted <- solved$predicted
  names(coefficients) <- traits$trait
  names(predicted) <- traits$trait
  list(
    coefficients = coefficients,
    target = target,
    predicted_response = predicted,
    economic_weights = economic,
    covariance = cov_info$target_matrix,
    cov_source = cov_info$source,
    cov_solve_form = cov_info$solve_form
  )
}

# Map user-supplied EXTERNAL covariances (raw trait units) into the value_z space the index is
# applied to: value_z = sign*(raw - center)/scale, so a raw covariance M maps as M_z = L M L with
# L = diag(sign/scale). Without this, raw-unit index weights multiply standardized/oriented values,
# distorting per-trait weights and flipping minimize-trait cross-covariance signs.
ng_multitrait_cov_to_value_z <- function(M, traits, value_scales) {
  if (is.null(M)) return(NULL)
  signs <- ifelse(traits$direction == "maximize", 1, -1)
  L <- signs / pmax(as.numeric(value_scales), 1e-8)
  L[!is.finite(L)] <- 1
  as.matrix(M) * outer(L, L)
}

ng_multitrait_economic_index_fit <- function(value_z, traits, ridge = 1e-6,
                                             phenotypic_covariance = NULL,
                                             genetic_covariance = NULL,
                                             value_scales = NULL) {
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  if (!p) ng_stop("economic-index optimizer needs at least one trait")
  economic <- ng_multitrait_economic_weights(traits)
  target <- economic
  target[!is.finite(target) | target < 0] <- 0
  if (!any(target > 0)) target <- rep(1 / p, p)
  target <- target / sum(target)
  names(target) <- traits$trait

  # External covariances arrive in raw units: map G/P and the raw-unit economic weights into the
  # value_z space so b'value_z reproduces the raw Smith-Hazel index b_raw'x_raw (a_z = scale * |a|).
  if ((!is.null(phenotypic_covariance) || !is.null(genetic_covariance)) && !is.null(value_scales)) {
    phenotypic_covariance <- ng_multitrait_cov_to_value_z(phenotypic_covariance, traits, value_scales)
    genetic_covariance <- ng_multitrait_cov_to_value_z(genetic_covariance, traits, value_scales)
    target <- target * pmax(as.numeric(value_scales), 1e-8)
    if (sum(abs(target)) > 0) target <- target / sum(abs(target))
    names(target) <- traits$trait
  }

  cov_info <- ng_multitrait_index_covariance(value_z, traits,
                                             phenotypic_covariance = phenotypic_covariance,
                                             genetic_covariance = genetic_covariance)
  solved <- ng_multitrait_solve_index(target, cov_info, ridge = ridge)
  coefficients <- solved$coefficients
  predicted <- solved$predicted
  names(coefficients) <- traits$trait
  names(predicted) <- traits$trait
  list(
    coefficients = coefficients,
    target = target,
    predicted_response = predicted,
    economic_weights = economic,
    covariance = cov_info$target_matrix,
    cov_source = cov_info$source,
    cov_solve_form = cov_info$solve_form
  )
}

ng_add_multitrait_score <- function(scores,
                                    traits,
                                    method = "auto",
                                    out_col = "multi_trait_score",
                                    strict_thresholds = FALSE,
                                    threshold_penalty_weight = 1.0,
                                    threshold_penalty_autoscale = TRUE,
                                    phenotypic_covariance = NULL,
                                    genetic_covariance = NULL,
                                    source = ng_multitrait_default_source()) {
  method <- trimws(tolower(as.character(method)[[1]]))
  methods <- c("auto", "weighted", "threshold", "economic_index", "desired_gain")
  if (!(method %in% methods)) ng_stop("method must be one of: ", paste(methods, collapse = ", "))
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  traits <- ng_multitrait_spec(traits)
  missing_cols <- setdiff(traits$column, names(scores))
  if (length(missing_cols)) ng_stop("scores missing multi-trait columns: ", paste(missing_cols, collapse = ", "))
  weights <- ng_multitrait_resolve_weights(traits, method)
  threshold_penalty_weight <- suppressWarnings(as.numeric(threshold_penalty_weight[[1]]))
  if (!is.finite(threshold_penalty_weight) || threshold_penalty_weight < 0) threshold_penalty_weight <- 1.0
  # Validate optional caller-supplied covariance matrices and align them with the
  # ordering of traits used internally (value_z columns are in traits$trait order).
  ng_multitrait_check_cov <- function(M, name) {
    if (is.null(M)) return(NULL)
    M <- as.matrix(M)
    if (!is.null(dimnames(M)) && all(traits$trait %in% rownames(M)) &&
        all(traits$trait %in% colnames(M))) {
      M <- M[traits$trait, traits$trait, drop = FALSE]
    } else if (nrow(M) != nrow(traits) || ncol(M) != nrow(traits)) {
      ng_stop(name, " must be square with rows/cols equal to the number of traits (", nrow(traits), ")")
    }
    M[!is.finite(M)] <- 0
    M
  }
  phenotypic_covariance <- ng_multitrait_check_cov(phenotypic_covariance, "phenotypic_covariance")
  genetic_covariance <- ng_multitrait_check_cov(genetic_covariance, "genetic_covariance")

  z <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  value_z <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  value_scales <- rep(1, nrow(traits))
  violation <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  colnames(z) <- colnames(value_z) <- colnames(violation) <- traits$trait
  for (i in seq_len(nrow(traits))) {
    x <- suppressWarnings(as.numeric(scores[[traits$column[[i]]]]))
    z[, i] <- ng_rank_normalize(x, bigger_is_better = identical(traits$direction[[i]], "maximize"))
    bad <- !is.finite(x)
    if (any(bad)) {
      finite_z <- z[is.finite(z[, i]), i]
      z[bad, i] <- if (length(finite_z)) min(finite_z, na.rm = TRUE) - 1 else -1
      violation[bad, i] <- violation[bad, i] + 1
    }
    scale <- ng_multitrait_value_scale(x)
    value_scales[[i]] <- scale
    oriented <- if (identical(traits$direction[[i]], "maximize")) x else -x
    finite_oriented <- oriented[is.finite(oriented)]
    center <- if (length(finite_oriented)) stats::median(finite_oriented, na.rm = TRUE) else 0
    value_z[, i] <- (oriented - center) / scale
    bad_value <- !is.finite(value_z[, i])
    if (any(bad_value)) {
      finite_value <- value_z[is.finite(value_z[, i]), i]
      value_z[bad_value, i] <- if (length(finite_value)) min(finite_value, na.rm = TRUE) - 1 else -1
    }
    if (is.finite(traits$min_value[[i]])) {
      violation[, i] <- violation[, i] + pmax(traits$min_value[[i]] - x, 0) / scale
    }
    if (is.finite(traits$max_value[[i]])) {
      violation[, i] <- violation[, i] + pmax(x - traits$max_value[[i]], 0) / scale
    }
    violation[!is.finite(violation[, i]), i] <- 1
    violation[, i] <- violation[, i] * traits$threshold_weight[[i]]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "z")]] <- z[, i]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "desired_z")]] <- value_z[, i]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "violation")]] <- violation[, i]
  }

  desired_gain <- NULL
  economic_index <- NULL
  penalty_weights <- weights
  if (identical(method, "desired_gain")) {
    desired_gain <- ng_multitrait_desired_gain_fit(value_z, traits, value_scales,
                                                   phenotypic_covariance = phenotypic_covariance,
                                                   genetic_covariance = genetic_covariance)
    weighted_score <- as.numeric(value_z %*% desired_gain$coefficients)
    penalty_weights <- desired_gain$target
  } else if (identical(method, "economic_index")) {
    economic_index <- ng_multitrait_economic_index_fit(value_z, traits,
                                                       phenotypic_covariance = phenotypic_covariance,
                                                       genetic_covariance = genetic_covariance,
                                                       value_scales = value_scales)
    weighted_score <- as.numeric(value_z %*% economic_index$coefficients)
    penalty_weights <- economic_index$target
  } else {
    weighted_score <- as.numeric(z %*% weights)
  }
  weighted_violation <- as.numeric(violation %*% penalty_weights)
  # Auto-scale the threshold penalty so a one-SD threshold violation costs
  # ~`threshold_penalty_weight * score_iqr` units of index score. Without
  # this, the multiplicative interaction of per-trait threshold_weight and
  # the global threshold_penalty_weight made it hard for users to tell
  # whether thresholds bound in practice.
  effective_penalty <- threshold_penalty_weight
  if (isTRUE(threshold_penalty_autoscale)) {
    score_iqr <- suppressWarnings(stats::IQR(weighted_score[is.finite(weighted_score)],
                                             na.rm = TRUE))
    if (is.finite(score_iqr) && score_iqr > 0) {
      effective_penalty <- threshold_penalty_weight * score_iqr
    }
  }
  scores$multi_trait_threshold_violation <- rowSums(violation)
  # Design invariant (deliberate; do not "normalize" this away): the emitted index
  # is oriented higher = better for EVERY method. This is guaranteed upstream of the
  # combination, not by rescaling out_col -- each trait is orientation-corrected
  # before standardization (minimize traits are sign-flipped in both `z` and
  # `value_z`) and aggregate index weights are constrained positive, so no method can
  # produce a "lower = better" index. The RAW SCALE of out_col is method-dependent by
  # design (rank-normal units for the weighted/rank-sum path, IQR units for the
  # economic_index/desired_gain solves); this is harmless because any consumer either
  # feeds out_col straight into allocation (ng_gain_scale rescales) or re-ingests it
  # as a trait (the trait loop above re-standardizes any column). Do not add a global
  # out_col normalizer without re-checking the threshold-penalty autoscale below,
  # which keys off IQR(weighted_score).
  scores[[out_col]] <- weighted_score - effective_penalty * weighted_violation
  if (isTRUE(strict_thresholds)) {
    scores[[out_col]][scores$multi_trait_threshold_violation > 1e-12] <- -Inf
  }
  attr(scores, "multi_trait") <- list(
    method = method,
    family = ng_multitrait_method_family(method),
    weights = weights,
    traits = traits,
    out_col = out_col,
    strict_thresholds = isTRUE(strict_thresholds),
    threshold_penalty_weight = threshold_penalty_weight,
    threshold_penalty_autoscale = isTRUE(threshold_penalty_autoscale),
    effective_threshold_penalty = effective_penalty,
    cov_source_phenotypic = if (is.null(phenotypic_covariance)) NA_character_ else "user_supplied",
    cov_source_genetic = if (is.null(genetic_covariance)) NA_character_ else "user_supplied",
    source = source
  )
  if (!is.null(desired_gain)) {
    attr(scores, "multi_trait")$desired_gain_coefficients <- desired_gain$coefficients
    attr(scores, "multi_trait")$desired_gain_target <- desired_gain$target
    attr(scores, "multi_trait")$desired_gain_predicted_response <- desired_gain$predicted_response
    attr(scores, "multi_trait")$economic_weights <- desired_gain$economic_weights
    attr(scores, "multi_trait")$desired_gain_covariance <- desired_gain$covariance
    attr(scores, "multi_trait")$desired_gain_cov_source <- desired_gain$cov_source
    attr(scores, "multi_trait")$desired_gain_cov_solve_form <- desired_gain$cov_solve_form
  }
  if (!is.null(economic_index)) {
    attr(scores, "multi_trait")$economic_index_coefficients <- economic_index$coefficients
    attr(scores, "multi_trait")$economic_index_target <- economic_index$target
    attr(scores, "multi_trait")$economic_index_predicted_response <- economic_index$predicted_response
    attr(scores, "multi_trait")$economic_weights <- economic_index$economic_weights
    attr(scores, "multi_trait")$economic_index_covariance <- economic_index$covariance
    attr(scores, "multi_trait")$economic_index_cov_source <- economic_index$cov_source
    attr(scores, "multi_trait")$economic_index_cov_solve_form <- economic_index$cov_solve_form
  }
  scores
}

ng_multitrait_add_desired_summary <- function(summary, meta) {
  if (!is.null(meta$desired_gain_coefficients)) {
    summary$multitrait_desired_gain_coefficients <- meta$desired_gain_coefficients
    summary$multitrait_desired_gain_target <- meta$desired_gain_target
    summary$multitrait_desired_gain_predicted_response <- meta$desired_gain_predicted_response
    summary$multitrait_economic_weights <- meta$economic_weights
  }
  if (!is.null(meta$economic_index_coefficients)) {
    summary$multitrait_economic_index_coefficients <- meta$economic_index_coefficients
    summary$multitrait_economic_index_target <- meta$economic_index_target
    summary$multitrait_economic_index_predicted_response <- meta$economic_index_predicted_response
    summary$multitrait_economic_weights <- meta$economic_weights
  }
  summary
}

ng_multitrait_summary <- function(plan, out_col, meta) {
  out <- list(
    n_crosses = nrow(plan),
    gain_col = out_col,
    total_gain = sum(plan[[out_col]], na.rm = TRUE),
    mean_gain = mean(plan[[out_col]], na.rm = TRUE),
    multitrait_method = meta$method,
    multitrait_family = meta$family,
    multitrait_trait_count = nrow(meta$traits),
    multitrait_source = meta$source,
    strict_thresholds = meta$strict_thresholds,
    threshold_penalty_weight = meta$threshold_penalty_weight,
    multitrait_weights = meta$weights
  )
  ng_multitrait_add_desired_summary(out, meta)
}

ng_multitrait_select_topn <- function(scores,
                                      traits,
                                      n_crosses,
                                      method = "auto",
                                      out_col = "multi_trait_score",
                                      strict_thresholds = FALSE,
                                      threshold_penalty_weight = 1.0,
                                      source = ng_multitrait_default_source()) {
  n_crosses <- suppressWarnings(as.integer(n_crosses[[1]]))
  if (!is.finite(n_crosses) || n_crosses < 1L) ng_stop("n_crosses must be a positive integer")
  scored <- ng_add_multitrait_score(
    scores = scores,
    traits = traits,
    method = method,
    out_col = out_col,
    strict_thresholds = strict_thresholds,
    threshold_penalty_weight = threshold_penalty_weight,
    source = source
  )
  feasible <- is.finite(scored[[out_col]])
  if (sum(feasible) < n_crosses) ng_stop("Not enough feasible multi-trait candidate crosses")
  tie <- if ("pair_kinship" %in% names(scored)) scored$pair_kinship else rep(0, nrow(scored))
  ord <- order(-scored[[out_col]], tie)
  ord <- ord[feasible[ord]]
  plan <- scored[head(ord, n_crosses), , drop = FALSE]
  attr(plan, "summary") <- ng_multitrait_summary(plan, out_col, attr(scored, "multi_trait"))
  plan
}

ng_score_breeder_objective <- function(scores,
                                       objective,
                                       out_col = "multi_trait_score",
                                       threshold_penalty_weight = 1.0,
                                       threshold_penalty_autoscale = TRUE,
                                       phenotypic_covariance = NULL,
                                       genetic_covariance = NULL) {
  if (!inherits(objective, "ng_breeder_selection_objective")) {
    ng_stop("objective must be created by ng_breeder_selection_objective()")
  }
  scored <- ng_add_multitrait_score(
    scores = scores,
    traits = objective$traits,
    method = objective$method,
    out_col = out_col,
    strict_thresholds = objective$strict_thresholds,
    threshold_penalty_weight = threshold_penalty_weight,
    threshold_penalty_autoscale = threshold_penalty_autoscale,
    phenotypic_covariance = phenotypic_covariance,
    genetic_covariance = genetic_covariance,
    source = objective$source
  )
  scored$multi_trait_threshold_policy <- objective$threshold_policy
  meta <- attr(scored, "multi_trait")
  meta$breeder_objective_method_requested <- objective$method_requested
  meta$breeder_objective_method_reason <- objective$diagnostics$method_reason
  attr(scored, "multi_trait") <- meta
  attr(scored, "breeder_objective") <- objective
  scored
}

ng_optimize_breeder_selection_plan <- function(scores,
                                               objective,
                                               n_crosses,
                                               parent_kinship = NULL,
                                               optimizer_method = "auto",
                                               out_col = "multi_trait_score",
                                               threshold_penalty_weight = 1.0,
                                               threshold_penalty_autoscale = TRUE,
                                               marker_target_spec = NULL,
                                               marker_geno = NULL,
                                               lambda_marker = 0,
                                               marker_ploidy = 2,
                                               ...) {
  scored <- ng_score_breeder_objective(
    scores = scores,
    objective = objective,
    out_col = out_col,
    threshold_penalty_weight = threshold_penalty_weight,
    threshold_penalty_autoscale = threshold_penalty_autoscale
  )
  n_crosses_int <- suppressWarnings(as.integer(n_crosses[[1]]))
  if (!is.finite(n_crosses_int) || n_crosses_int < 1L) {
    ng_stop("n_crosses must be a positive integer")
  }
  feasible <- is.finite(scored[[out_col]])
  if (sum(feasible) < n_crosses_int) {
    ng_stop("breeder objective leaves only ", sum(feasible),
            " feasible crosses for n_crosses = ", n_crosses_int,
            "; relax strict thresholds or request fewer crosses")
  }
  if (isTRUE(objective$strict_thresholds)) {
    scored <- scored[feasible, , drop = FALSE]
  }
  # Optional marker steering (Module 4): blend a target-allele score into the merit column so
  # the allocator biases toward the desired marker frequency, exactly as ng_design_crosses does.
  gain_col_used <- out_col
  if (!is.null(marker_target_spec)) {
    scored <- ng_apply_marker_management(
      scores = scored, geno = marker_geno,
      marker_target_spec = marker_target_spec, ploidy = marker_ploidy,
      gain_col = out_col, lambda_marker = lambda_marker
    )
    if ("marker_adjusted_gain" %in% names(scored)) gain_col_used <- "marker_adjusted_gain"
  }
  plan <- ng_optimize_mating_plan(
    scores = scored,
    n_crosses = n_crosses_int,
    gain_col = gain_col_used,
    parent_kinship = parent_kinship,
    method = optimizer_method,
    ...
  )
  summary <- attr(plan, "summary")
  meta <- attr(scored, "multi_trait")
  summary$breeder_objective_method <- objective$method
  summary$breeder_objective_method_requested <- objective$method_requested
  summary$breeder_objective_method_reason <- objective$diagnostics$method_reason
  summary$breeder_objective_threshold_policy <- objective$threshold_policy
  summary$breeder_objective_zero_violation_selected <- sum(
    plan$multi_trait_threshold_violation <= 1e-12,
    na.rm = TRUE
  )
  summary$breeder_objective_trait_count <- nrow(objective$traits)
  summary$multitrait_method <- meta$method
  summary$multitrait_family <- meta$family
  summary$multitrait_source <- meta$source
  summary$multitrait_weights <- meta$weights
  summary <- ng_multitrait_add_desired_summary(summary, meta)
  attr(plan, "summary") <- summary
  attr(plan, "breeder_objective") <- objective
  plan
}

ng_optimize_multitrait_mating_plan <- function(scores,
                                               traits,
                                               n_crosses,
                                               parent_kinship = NULL,
                                               multitrait_method = "auto",
                                               optimizer_method = "auto",
                                               out_col = "multi_trait_score",
                                               strict_thresholds = FALSE,
                                               threshold_penalty_weight = 1.0,
                                               source = ng_multitrait_default_source(),
                                               ...) {
  scored <- ng_add_multitrait_score(
    scores = scores,
    traits = traits,
    method = multitrait_method,
    out_col = out_col,
    strict_thresholds = strict_thresholds,
    threshold_penalty_weight = threshold_penalty_weight,
    source = source
  )
  meta <- attr(scored, "multi_trait")
  plan <- ng_optimize_mating_plan(
    scores = scored,
    n_crosses = n_crosses,
    gain_col = out_col,
    parent_kinship = parent_kinship,
    method = optimizer_method,
    ...
  )
  summary <- attr(plan, "summary")
  summary$multitrait_method <- meta$method
  summary$multitrait_family <- meta$family
  summary$multitrait_trait_count <- nrow(meta$traits)
  summary$multitrait_source <- meta$source
  summary$strict_thresholds <- meta$strict_thresholds
  summary$threshold_penalty_weight <- meta$threshold_penalty_weight
  summary$multitrait_weights <- meta$weights
  summary <- ng_multitrait_add_desired_summary(summary, meta)
  attr(plan, "summary") <- summary
  plan
}
