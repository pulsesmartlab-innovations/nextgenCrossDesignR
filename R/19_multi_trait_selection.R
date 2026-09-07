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

# DEFECT 2 fix (threshold-fix report): ng_preflight_input_tables() (R/26_data_preflight.R)
# accepts min/minimum as aliases for min_value (max/maximum for max_value) and even
# cross-validates min <= max -- but the scorer's spec intake recognised only the literal
# "min_value"/"max_value" names, so a spec using an alias passed preflight's "thresholds are
# consistent" check and then had its thresholds silently replaced by NA here. Reuse R/26's own
# alias lists (ng_preflight_min_value_aliases / ng_preflight_max_value_aliases) so the two
# cannot drift apart again; do not copy them.
#
# If the spec somehow carries BOTH the canonical column and an alias (e.g. both min_value and
# min), the canonical name wins and a warning is raised -- never a silent pick.
ng_multitrait_apply_threshold_alias <- function(spec, canonical, aliases) {
  hit <- ng_preflight_id_column(spec, aliases)
  if (is.null(hit)) return(spec)
  if (canonical %in% names(spec)) {
    collided <- setdiff(intersect(aliases, names(spec)), canonical)
    if (length(collided)) {
      warning("multi-trait spec has both '", canonical, "' and '", collided[[1L]],
              "'; using '", canonical, "' and ignoring '", collided[[1L]], "'", call. = FALSE)
    }
    return(spec)
  }
  spec[[canonical]] <- spec[[hit]]
  spec
}

ng_multitrait_spec <- function(trait,
                               column = NULL,
                               direction = NULL,
                               weight = NULL,
                               min_value = NULL,
                               max_value = NULL,
                               threshold_weight = NULL,
                               desired_change = NULL,
                               economic_weight = NULL,
                               threshold_column = NULL) {
  if (is.data.frame(trait)) {
    spec <- as.data.frame(trait, stringsAsFactors = FALSE)
    if (!("column" %in% names(spec)) && "trait" %in% names(spec)) spec$column <- spec$trait
    if (!("trait" %in% names(spec)) && "column" %in% names(spec)) spec$trait <- spec$column
    if (!("direction" %in% names(spec))) spec$direction <- "maximize"
    if (!("weight" %in% names(spec))) spec$weight <- NA_real_
    spec <- ng_multitrait_apply_threshold_alias(spec, "min_value", ng_preflight_min_value_aliases)
    spec <- ng_multitrait_apply_threshold_alias(spec, "max_value", ng_preflight_max_value_aliases)
    if (!("min_value" %in% names(spec))) spec$min_value <- NA_real_
    if (!("max_value" %in% names(spec))) spec$max_value <- NA_real_
    if (!("threshold_weight" %in% names(spec))) spec$threshold_weight <- 1
    if (!("desired_change" %in% names(spec))) spec$desired_change <- NA_real_
    if (!("economic_weight" %in% names(spec))) spec$economic_weight <- NA_real_
    # DEFECT 1 fix: a breeder's min_value/max_value are trait-unit thresholds; `column` is
    # whatever metric the run ranks on (e.g. usefulness = mean + i*sqrt(pmv) under the default
    # trait_value_metric), which is NOT trait units. threshold_column gives the threshold
    # comparison its own, independently-controlled basis. Defaulting it to `column` when absent
    # keeps every direct caller of ng_add_multitrait_score()/ng_breeder_selection_objective()
    # working exactly as before: they control both the spec and the `scores` table, so their
    # `column` IS their intended comparison basis. R/39_cross_prediction_runner.R sets
    # threshold_column explicitly to `<trait>_mean`.
    if (!("threshold_column" %in% names(spec))) spec$threshold_column <- spec$column
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
    threshold_column <- ng_multitrait_recycle(threshold_column, n, "threshold_column", column)
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
      threshold_column = as.character(threshold_column),
      stringsAsFactors = FALSE
    )
  }

  required <- c("trait", "column", "direction", "weight", "min_value", "max_value",
                "threshold_weight", "desired_change", "economic_weight", "threshold_column")
  missing <- setdiff(required, names(spec))
  if (length(missing)) ng_stop("multi-trait spec missing columns: ", paste(missing, collapse = ", "))
  spec <- spec[required]
  spec$trait <- trimws(as.character(spec$trait))
  spec$column <- trimws(as.character(spec$column))
  spec$threshold_column <- trimws(as.character(spec$threshold_column))
  if (any(!nzchar(spec$trait) | is.na(spec$trait))) ng_stop("trait names must be non-empty")
  if (any(!nzchar(spec$column) | is.na(spec$column))) ng_stop("trait columns must be non-empty")
  if (any(!nzchar(spec$threshold_column) | is.na(spec$threshold_column))) {
    ng_stop("threshold_column must be non-empty")
  }
  if (any(duplicated(spec$trait))) ng_stop("trait names must be unique")
  spec$direction <- ng_multitrait_direction(spec$direction)
  spec$weight <- suppressWarnings(as.numeric(spec$weight))
  spec$min_value <- suppressWarnings(as.numeric(spec$min_value))
  spec$max_value <- suppressWarnings(as.numeric(spec$max_value))
  spec$threshold_weight <- suppressWarnings(as.numeric(spec$threshold_weight))
  spec$desired_change <- suppressWarnings(as.numeric(spec$desired_change))
  spec$economic_weight <- suppressWarnings(as.numeric(spec$economic_weight))
  if (any(!is.finite(spec$threshold_weight) | spec$threshold_weight < 0)) {
    ng_stop("threshold_weight must be finite and non-negative for every trait")
  }
  if (any(is.finite(spec$weight) & spec$weight < 0)) {
    ng_stop("weight must be non-negative; use direction = 'minimize' for an unfavorable trait")
  }
  if (any(is.finite(spec$economic_weight) & spec$economic_weight < 0)) {
    ng_stop("economic_weight must be non-negative; use direction = 'minimize' for an unfavorable trait")
  }
  if (any(is.finite(spec$desired_change) & spec$desired_change < 0)) {
    ng_stop("desired_change must be non-negative; use direction = 'minimize' for a requested decrease")
  }
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
                                           threshold_column = NULL,
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
    economic_weight = economic_weight,
    threshold_column = threshold_column
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
    if (has_desired) {
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
    if (any(!is.finite(raw) | raw < 0) || !any(raw > 0)) {
      ng_stop("weighted multi-trait scoring requires one finite non-negative weight per trait and at least one positive weight")
    }
  } else if (identical(method, "economic_index")) {
    return(ng_multitrait_economic_weights(spec, required = TRUE))
  } else if (identical(method, "desired_gain")) {
    raw <- suppressWarnings(as.numeric(spec$desired_change))
    if (any(!is.finite(raw) | raw < 0) || !any(raw > 0)) {
      ng_stop("desired_gain requires one finite non-negative desired_change per trait and at least one positive desired_change")
    }
  } else if (any(is.finite(raw))) {
    if (any(!is.finite(raw) | raw < 0) || !any(raw > 0)) {
      ng_stop("partially specified trait weights are not allowed; provide one finite non-negative weight per trait with at least one positive weight")
    }
  }

  if (!any(is.finite(raw))) {
    raw <- rep(1, nrow(spec))
  }
  weights <- raw / sum(raw)
  names(weights) <- spec$trait
  weights
}

ng_multitrait_economic_weights <- function(spec, required = TRUE) {
  raw <- suppressWarnings(as.numeric(spec$economic_weight))
  if (!any(is.finite(raw)) && !isTRUE(required)) {
    out <- rep(NA_real_, nrow(spec))
    names(out) <- spec$trait
    return(out)
  }
  if (any(!is.finite(raw) | raw < 0) || !any(raw > 0)) {
    ng_stop("economic_index requires one finite non-negative economic_weight per trait and at least one positive economic_weight")
  }
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

# --- Label-aware alignment of caller-supplied trait x trait matrices ---------
#
# SAFETY (0.27.0). P and G can now arrive from a GUI over a JSON bridge, and a
# jsonlite matrix round trip DROPS dimnames: jsonlite::fromJSON(toJSON(M))
# returns an unnamed matrix. Every internal use of P/G is positional -- both the
# index solve and the value_z rescaling L = diag(sign / scale) applied by
# ng_multitrait_cov_to_value_z() -- so a matrix whose rows/columns are ordered
# differently from `traits` would silently get the wrong direction sign and the
# wrong scale applied to every element and still return a plausible index.
#
# Contract:
#   * LABELLED input (dimnames on BOTH dimensions) is REORDERED BY NAME into the
#     internal trait order and must name exactly the trait set. A missing, extra
#     or misspelled label is a hard error that names the offending labels.
#     Labels on only one dimension are a hard error too, because the unlabelled
#     dimension cannot be checked against anything.
#   * UNLABELLED input (no dimnames at all) is INTERPRETED POSITIONALLY in
#     `traits$trait` order, and only when the dimension matches. This is the
#     documented behaviour for the JSON bridge; a caller that cannot guarantee
#     column order MUST send dimnames.
ng_multitrait_align_cov <- function(M, trait_names, name) {
  if (is.null(M)) return(NULL)
  M <- as.matrix(M)
  trait_names <- as.character(trait_names)
  p <- length(trait_names)
  rn <- rownames(M)
  cn <- colnames(M)
  if (is.null(rn) && is.null(cn)) {
    if (nrow(M) != p || ncol(M) != p) {
      ng_stop(name, " must be a ", p, " x ", p, " matrix (one row/column per trait); an ",
              "unlabelled matrix is interpreted POSITIONALLY in trait order: ",
              paste(trait_names, collapse = ", "))
    }
    dimnames(M) <- list(trait_names, trait_names)
    return(M)
  }
  if (is.null(rn) || is.null(cn)) {
    ng_stop(name, " has dimnames on only one dimension; supply trait names on BOTH ",
            "rownames and colnames, or on neither (positional interpretation)")
  }
  check_side <- function(labels, side) {
    dups <- unique(labels[duplicated(labels)])
    if (length(dups)) {
      ng_stop(name, " has duplicated ", side, " label(s): ", paste(dups, collapse = ", "))
    }
    missing <- setdiff(trait_names, labels)
    extra <- setdiff(labels, trait_names)
    if (length(missing) || length(extra)) {
      ng_stop(name, " ", side, " labels do not match the trait set. Missing: ",
              if (length(missing)) paste(missing, collapse = ", ") else "<none>",
              "; unexpected: ", if (length(extra)) paste(extra, collapse = ", ") else "<none>",
              ". Expected exactly: ", paste(trait_names, collapse = ", "))
    }
  }
  check_side(rn, "row")
  check_side(cn, "column")
  M <- M[trait_names, trait_names, drop = FALSE]
  dimnames(M) <- list(trait_names, trait_names)
  M
}

# Label-aware alignment PLUS the covariance-matrix validity checks: finite,
# symmetric within tolerance, positive diagonal variances, and positive
# semidefinite by eigenvalue.
ng_multitrait_validate_cov <- function(M, trait_names, name) {
  M <- ng_multitrait_align_cov(M, trait_names, name)
  if (is.null(M)) return(NULL)
  if (any(!is.finite(M))) ng_stop(name, " contains non-finite entries")
  asym <- max(abs(M - t(M)))
  if (!is.finite(asym) || asym > 1e-8) ng_stop(name, " must be symmetric")
  keep <- dimnames(M)
  M <- (M + t(M)) / 2
  dimnames(M) <- keep
  if (any(diag(M) <= 0)) ng_stop(name, " must have positive diagonal variances")
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  tol <- 1e-8 * max(1, max(abs(ev)))
  if (min(ev) < -tol) ng_stop(name, " must be positive semidefinite")
  M
}

# Build the covariance matrices used to solve a formal selection index.
#
# economic_index (Smith-Hazel, b = P^{-1} G a) genuinely requires BOTH P and G.
#
# desired_gain (Pesek-Baker, b = G^{-1} d) requires ONLY G. P never enters the
# coefficient solve: for this purpose target_matrix = G and projection = NULL,
# so ng_multitrait_solve_index() forms b from G and the desired-gain target
# alone. P is used solely to standardise the REPORTED predicted response by the
# index standard deviation sqrt(b' P b). When P is absent the coefficients and
# the emitted index are unchanged and the reported predicted response is NA;
# the names of the unavailable quantities are returned in `unavailable`.
#
# Candidate-cross score covariance is never a substitute for either matrix and
# is not inferred from the scores under any purpose.
ng_multitrait_index_covariance <- function(value_z, traits,
                                           phenotypic_covariance = NULL,
                                           genetic_covariance = NULL,
                                           purpose = c("economic_index", "desired_gain")) {
  purpose <- match.arg(purpose)
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  trait_names <- as.character(traits$trait)
  if (length(trait_names) != p) {
    ng_stop("internal error: value_z has ", p, " columns but the trait spec has ",
            length(trait_names), " traits")
  }
  G <- ng_multitrait_validate_cov(genetic_covariance, trait_names, "genetic_covariance")
  P <- ng_multitrait_validate_cov(phenotypic_covariance, trait_names, "phenotypic_covariance")
  if (identical(purpose, "economic_index")) {
    if (is.null(G) || is.null(P)) {
      ng_stop(purpose, " requires both phenotypic_covariance (P) and genetic_covariance (G) ",
              "for the Smith-Hazel solve b = P^{-1} G a; ",
              "candidate-score covariance is not a substitute for quantitative-genetic covariance")
    }
    return(list(target_matrix = P, projection = G,
                response_G = G, response_P = P,
                source = "smith_hazel", solve_form = "b = P^{-1} G a",
                unavailable = character(0)))
  }
  if (is.null(G)) {
    ng_stop(purpose, " requires genetic_covariance (G) for the Pesek-Baker solve b = G^{-1} d; ",
            "candidate-score covariance is not a substitute for quantitative-genetic covariance")
  }
  list(target_matrix = G, projection = NULL,
       response_G = G, response_P = P,
       source = "pesek_baker", solve_form = "b = G^{-1} d",
       unavailable = if (is.null(P)) "predicted_response" else character(0))
}

# Single-ridge solve: solve_mat = target_matrix + ridge * diag_scale * I.
# Previously the diagonal was pre-inflated AND a ridge term was added, applying
# the penalty twice. Returns the coefficients and the predicted response under
# the index theory:
#   - smith_hazel: predicted = G %*% b (Smith-Hazel realized response)
#   - pesek_baker: predicted = G %*% b (desired-gain response direction)
ng_multitrait_solve_index <- function(target, cov_info, ridge = 1e-6) {
  p <- length(target)
  M <- cov_info$target_matrix
  ridge <- suppressWarnings(as.numeric(ridge))
  if (length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    ng_stop("index ridge must be one finite non-negative number")
  }
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
  # Symmetric Moore-Penrose solve. With ridge > 0 this is the ordinary
  # inverse; with ridge = 0 it gives the uniquely defined minimum-norm
  # solution for a singular PSD P or G. Falling back to `target` would not
  # satisfy either the Smith-Hazel or Pesek-Baker equation.
  es <- eigen(solve_mat, symmetric = TRUE)
  solve_tol <- max(dim(solve_mat)) * .Machine$double.eps * max(1, max(abs(es$values)))
  keep <- es$values > solve_tol
  if (!any(keep)) ng_stop("selection-index covariance has no estimable positive-eigenvalue subspace")
  coefficients <- as.numeric(es$vectors[, keep, drop = FALSE] %*%
    (crossprod(es$vectors[, keep, drop = FALSE], rhs) / es$values[keep]))
  if (any(!is.finite(coefficients)) || !any(abs(coefficients) > solve_tol)) {
    ng_stop("selection-index target has no estimable component in the covariance column space")
  }
  coefficients <- coefficients / sum(abs(coefficients))
  # The predicted response is reported in index-SD units, sigma_I = sqrt(b' P b).
  # P is REQUIRED for that scaling but never for the coefficients: the
  # Pesek-Baker branch (projection = NULL, target_matrix = G) solves b = G^{-1} d
  # without ever touching P. When P is absent the coefficients are unchanged and
  # the predicted response degrades to NA rather than failing the whole solve.
  sigma_i <- if (is.null(cov_info$response_P)) {
    NA_real_
  } else {
    sqrt(max(as.numeric(crossprod(
      coefficients, cov_info$response_P %*% coefficients)), 0))
  }
  predicted <- if (is.na(sigma_i)) {
    rep(NA_real_, p)
  } else if (is.finite(sigma_i) && sigma_i > 0) {
    as.numeric(cov_info$response_G %*% coefficients) / sigma_i
  } else {
    rep(0, p)
  }
  if (!anyNA(predicted)) predicted[!is.finite(predicted)] <- 0
  list(coefficients = coefficients, predicted = predicted, index_sd = sigma_i)
}

ng_multitrait_desired_gain_fit <- function(value_z, traits, value_scales, ridge = 1e-6,
                                           phenotypic_covariance = NULL,
                                           genetic_covariance = NULL) {
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  if (!p) ng_stop("desired-gain optimizer needs at least one trait")
  economic <- ng_multitrait_economic_weights(traits, required = FALSE)
  desired <- suppressWarnings(as.numeric(traits$desired_change))
  if (any(!is.finite(desired)) || any(desired < 0)) {
    ng_stop("desired_gain requires an explicit finite non-negative desired_change for every trait; ",
            "use 0 for a trait with no requested response")
  }
  if (!any(desired > 0)) ng_stop("desired_gain requires at least one positive desired_change")
  desired_sd <- desired / pmax(value_scales, 1e-8)
  # Pesek-Baker desired gains replace economic weights; multiplying desired
  # changes by economic weights would define a different target.
  target <- desired_sd
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
                                             genetic_covariance = genetic_covariance,
                                             purpose = "desired_gain")
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
    cov_solve_form = cov_info$solve_form,
    index_sd = solved$index_sd,
    # Quantities that could not be computed from the covariance matrices that
    # were actually supplied (empty when everything was available).
    unavailable = cov_info$unavailable
  )
}

# Map user-supplied EXTERNAL covariances (raw trait units) into the value_z space the index is
# applied to: value_z = sign*(raw - center)/scale, so a raw covariance M maps as M_z = L M L with
# L = diag(sign/scale). Without this, raw-unit index weights multiply standardized/oriented values,
# distorting per-trait weights and flipping minimize-trait cross-covariance signs.
# `M` is aligned to `traits$trait` BY NAME when it carries dimnames, and is
# interpreted POSITIONALLY (in traits$trait order) only when it carries none --
# see ng_multitrait_align_cov() for the full contract and for why a positional-
# only transform is unsafe for GUI/JSON input.
ng_multitrait_cov_to_value_z <- function(M, traits, value_scales) {
  if (is.null(M)) return(NULL)
  trait_names <- as.character(traits$trait)
  M <- ng_multitrait_align_cov(M, trait_names, "covariance matrix")
  scales <- as.numeric(value_scales)
  # value_scales is built in traits$trait order upstream; when it is named, take
  # the same by-name route so scale and sign can never disagree with M.
  if (!is.null(names(value_scales)) && all(trait_names %in% names(value_scales))) {
    scales <- as.numeric(value_scales[trait_names])
  }
  signs <- ifelse(traits$direction == "maximize", 1, -1)
  L <- signs / pmax(scales, 1e-8)
  L[!is.finite(L)] <- 1
  out <- as.matrix(M) * outer(L, L)
  dimnames(out) <- list(trait_names, trait_names)
  out
}

ng_multitrait_economic_index_fit <- function(value_z, traits, ridge = 1e-6,
                                             phenotypic_covariance = NULL,
                                             genetic_covariance = NULL,
                                             value_scales = NULL) {
  value_z <- as.matrix(value_z)
  p <- ncol(value_z)
  if (!p) ng_stop("economic-index optimizer needs at least one trait")
  economic <- ng_multitrait_economic_weights(traits, required = TRUE)
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
                                             genetic_covariance = genetic_covariance,
                                             purpose = "economic_index")
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
    cov_solve_form = cov_info$solve_form,
    index_sd = solved$index_sd,
    # Quantities that could not be computed from the covariance matrices that
    # were actually supplied (empty when everything was available).
    unavailable = cov_info$unavailable
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
  if (identical(method, "auto")) {
    has_positive <- function(x) {
      x <- suppressWarnings(as.numeric(x))
      any(is.finite(x) & x > 0)
    }
    method <- if (has_positive(traits$desired_change)) {
      "desired_gain"
    } else if (has_positive(traits$economic_weight)) {
      "economic_index"
    } else if (has_positive(traits$weight)) {
      "weighted"
    } else {
      "auto"
    }
  }
  missing_cols <- setdiff(traits$column, names(scores))
  if (length(missing_cols)) ng_stop("scores missing multi-trait columns: ", paste(missing_cols, collapse = ", "))
  # threshold_column naming a column absent from `scores` is a programming error, not a silent
  # fallback (see DEFECT 1 in the threshold-fix report): name the trait and the missing column
  # so the caller can see exactly which threshold_column is wrong.
  missing_threshold <- is.finite(traits$min_value) | is.finite(traits$max_value)
  missing_threshold <- missing_threshold & !(traits$threshold_column %in% names(scores))
  if (any(missing_threshold)) {
    bad <- traits[missing_threshold, c("trait", "threshold_column"), drop = FALSE]
    ng_stop("scores missing threshold_column for trait(s): ",
            paste(sprintf("%s -> %s", bad$trait, bad$threshold_column), collapse = ", "))
  }
  weights <- ng_multitrait_resolve_weights(traits, method)
  threshold_penalty_weight <- suppressWarnings(as.numeric(threshold_penalty_weight[[1]]))
  if (!is.finite(threshold_penalty_weight) || threshold_penalty_weight < 0) threshold_penalty_weight <- 1.0
  # Validate optional caller-supplied covariance matrices and align them with the
  # ordering of traits used internally (value_z columns are in traits$trait order).
  # 0.27.0: this used to fall through to a POSITIONAL read whenever the supplied
  # dimnames did not happen to contain every trait -- so a matrix carrying a
  # misspelled trait, or labels on one dimension only, was silently accepted in
  # whatever order it arrived. ng_multitrait_validate_cov() reorders by name and
  # hard-errors on any label mismatch instead.
  phenotypic_covariance <- ng_multitrait_validate_cov(phenotypic_covariance, traits$trait,
                                                      "phenotypic_covariance")
  genetic_covariance <- ng_multitrait_validate_cov(genetic_covariance, traits$trait,
                                                   "genetic_covariance")

  z <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  value_z <- matrix(0, nrow = nrow(scores), ncol = nrow(traits))
  value_scales <- stats::setNames(rep(1, nrow(traits)), traits$trait)
  # Retained alongside value_scales so a check-reference value (which never appears as a row in
  # `scores`, so it has no rank of its own) can be placed on this SAME affine axis after the
  # fact: oriented <- sign * tau; z <- (oriented - center) / scale. See ng_check_line_value()
  # (R/38_cross_priority_plots.R), which is the sole consumer of these three vectors.
  value_centers <- stats::setNames(rep(0, nrow(traits)), traits$trait)
  value_signs <- stats::setNames(rep(1, nrow(traits)), traits$trait)
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
    sign_i <- if (identical(traits$direction[[i]], "maximize")) 1 else -1
    value_signs[[i]] <- sign_i
    oriented <- if (identical(traits$direction[[i]], "maximize")) x else -x
    finite_oriented <- oriented[is.finite(oriented)]
    center <- if (length(finite_oriented)) stats::median(finite_oriented, na.rm = TRUE) else 0
    value_centers[[i]] <- center
    value_z[, i] <- (oriented - center) / scale
    bad_value <- !is.finite(value_z[, i])
    if (any(bad_value)) {
      finite_value <- value_z[is.finite(value_z[, i]), i]
      value_z[bad_value, i] <- if (length(finite_value)) min(finite_value, na.rm = TRUE) - 1 else -1
    }
    # DEFECT 1 fix: min_value/max_value are trait-unit thresholds set by a breeder ("the family
    # mean must reach 70"), but `x` above is whatever metric this run ranks on -- under the
    # default trait_value_metric = "usefulness" that is mean + i*sqrt(pmv), not the mean, and
    # under "pmv"/"vpm" it is a variance. Comparing min_value/max_value against `x` mixes units
    # and silently biases toward "no violation" in both directions. threshold_column (defaults
    # to `column` -- see ng_multitrait_spec()) gives the threshold its own, correctly-scaled
    # axis; the ranking z-score above is deliberately left on `column`/x.
    x_thresh <- suppressWarnings(as.numeric(scores[[traits$threshold_column[[i]]]]))
    # DEFECT 3 fix (0.26.0): the deficit above is measured on threshold_column, so it must be
    # made dimensionless by the dispersion of THAT column, not by `scale` (the robust spread of
    # the RANKING column `x`). Dividing a trait-unit deficit by the IQR of whatever
    # trait_value_metric produced mixes two different quantities in two different units: under
    # trait_value_metric = "pmv"/"vpm"/"var_complex" the ranking column is a VARIANCE, and the
    # measured effect was a 52x swing in how hard the same thresholds bite (mean penalty 0.35
    # under "usefulness" vs 18.0 under "pmv" on identical data and identical thresholds) driven
    # entirely by a knob that has nothing to do with thresholds. The penalty must be invariant
    # to the ranking metric when the threshold column and the thresholds are unchanged.
    # ng_multitrait_value_scale() never returns 0 or a non-finite value (it falls back
    # IQR -> MAD -> SD -> range -> 1), so this can never divide by zero; when threshold_column
    # IS column (the default for every direct caller) thresh_scale == scale exactly and the
    # penalty is bit-identical to 0.25.0.
    thresh_scale <- if (identical(traits$threshold_column[[i]], traits$column[[i]])) scale else
      ng_multitrait_value_scale(x_thresh)
    if (is.finite(traits$min_value[[i]])) {
      violation[, i] <- violation[, i] + pmax(traits$min_value[[i]] - x_thresh, 0) / thresh_scale
    }
    if (is.finite(traits$max_value[[i]])) {
      violation[, i] <- violation[, i] + pmax(x_thresh - traits$max_value[[i]], 0) / thresh_scale
    }
    violation[!is.finite(violation[, i]), i] <- 1
    violation[, i] <- violation[, i] * traits$threshold_weight[[i]]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "z")]] <- z[, i]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "desired_z")]] <- value_z[, i]
    scores[[ng_multitrait_diag_col(traits$trait[[i]], "violation")]] <- violation[, i]
  }

  desired_gain <- NULL
  economic_index <- NULL
  if (identical(method, "desired_gain")) {
    desired_gain <- ng_multitrait_desired_gain_fit(value_z, traits, value_scales,
                                                   phenotypic_covariance = phenotypic_covariance,
                                                   genetic_covariance = genetic_covariance)
    weighted_score <- as.numeric(value_z %*% desired_gain$coefficients)
  } else if (identical(method, "economic_index")) {
    economic_index <- ng_multitrait_economic_index_fit(value_z, traits,
                                                       phenotypic_covariance = phenotypic_covariance,
                                                       genetic_covariance = genetic_covariance,
                                                       value_scales = value_scales)
    weighted_score <- as.numeric(value_z %*% economic_index$coefficients)
  } else {
    weighted_score <- as.numeric(z %*% weights)
  }
  # Threshold importance is controlled explicitly by threshold_weight in the
  # trait specification. Economic weights or desired-gain targets must not
  # silently weaken a biological/quality threshold (especially when a desired
  # response is exactly zero), so soft violations add independently here.
  weighted_violation <- rowSums(violation)
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
  # combination, not by rescaling out_col: every breeding objective is oriented so
  # selection proceeds upward. Smith-Hazel coefficients themselves need not all be
  # positive because a correlated trait can improve prediction of aggregate merit.
  # The RAW SCALE of out_col is method-dependent by
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
    source = source,
    value_centers = value_centers,
    value_scales = value_scales,
    value_signs = value_signs
  )
  if (!is.null(desired_gain)) {
    attr(scores, "multi_trait")$desired_gain_coefficients <- desired_gain$coefficients
    attr(scores, "multi_trait")$desired_gain_target <- desired_gain$target
    attr(scores, "multi_trait")$desired_gain_predicted_response <- desired_gain$predicted_response
    attr(scores, "multi_trait")$economic_weights <- desired_gain$economic_weights
    attr(scores, "multi_trait")$desired_gain_covariance <- desired_gain$covariance
    attr(scores, "multi_trait")$desired_gain_cov_source <- desired_gain$cov_source
    attr(scores, "multi_trait")$desired_gain_cov_solve_form <- desired_gain$cov_solve_form
    attr(scores, "multi_trait")$desired_gain_index_sd <- desired_gain$index_sd
    # Pesek-Baker never needs P for the coefficients; when P was not supplied the
    # scoring is complete but any P-scaled REPORTED quantity is NA and named here.
    attr(scores, "multi_trait")$desired_gain_unavailable <- desired_gain$unavailable
    attr(scores, "multi_trait")$desired_gain_unavailable_reason <-
      if (length(desired_gain$unavailable)) {
        paste0(paste(desired_gain$unavailable, collapse = ", "),
               " unavailable: reporting it requires phenotypic_covariance (P) for the ",
               "index standard deviation sqrt(b' P b). The desired_gain coefficients ",
               "(b = G^{-1} d) and the emitted index are unaffected.")
      } else {
        NA_character_
      }
  }
  if (!is.null(economic_index)) {
    attr(scores, "multi_trait")$economic_index_coefficients <- economic_index$coefficients
    attr(scores, "multi_trait")$economic_index_target <- economic_index$target
    attr(scores, "multi_trait")$economic_index_predicted_response <- economic_index$predicted_response
    attr(scores, "multi_trait")$economic_weights <- economic_index$economic_weights
    attr(scores, "multi_trait")$economic_index_covariance <- economic_index$covariance
    attr(scores, "multi_trait")$economic_index_cov_source <- economic_index$cov_source
    attr(scores, "multi_trait")$economic_index_cov_solve_form <- economic_index$cov_solve_form
    attr(scores, "multi_trait")$economic_index_index_sd <- economic_index$index_sd
    attr(scores, "multi_trait")$economic_index_unavailable <- economic_index$unavailable
  }
  scores
}

ng_multitrait_add_desired_summary <- function(summary, meta) {
  if (!is.null(meta$desired_gain_coefficients)) {
    summary$multitrait_desired_gain_coefficients <- meta$desired_gain_coefficients
    summary$multitrait_desired_gain_target <- meta$desired_gain_target
    summary$multitrait_desired_gain_predicted_response <- meta$desired_gain_predicted_response
    summary$multitrait_economic_weights <- meta$economic_weights
    summary$multitrait_desired_gain_unavailable <- meta$desired_gain_unavailable
    summary$multitrait_desired_gain_unavailable_reason <- meta$desired_gain_unavailable_reason
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
                                      phenotypic_covariance = NULL,
                                      genetic_covariance = NULL,
                                      source = ng_multitrait_default_source()) {
  n_crosses_num <- suppressWarnings(as.numeric(n_crosses))
  if (length(n_crosses_num) != 1L || !is.finite(n_crosses_num) || n_crosses_num < 1 ||
      abs(n_crosses_num - round(n_crosses_num)) > 1e-8) {
    ng_stop("n_crosses must be a positive integer")
  }
  n_crosses <- as.integer(round(n_crosses_num))
  scored <- ng_add_multitrait_score(
    scores = scores,
    traits = traits,
    method = method,
    out_col = out_col,
    strict_thresholds = strict_thresholds,
    threshold_penalty_weight = threshold_penalty_weight,
    phenotypic_covariance = phenotypic_covariance,
    genetic_covariance = genetic_covariance,
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
  # Preserve the data-frame shape for an empty candidate set. Assigning a
  # scalar to a zero-row data frame raises before the caller can issue the
  # intended "not enough feasible candidates" diagnostic.
  scored$multi_trait_threshold_policy <- rep(objective$threshold_policy, nrow(scored))
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
                                               phenotypic_covariance = NULL,
                                               genetic_covariance = NULL,
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
    threshold_penalty_autoscale = threshold_penalty_autoscale,
    phenotypic_covariance = phenotypic_covariance,
    genetic_covariance = genetic_covariance
  )
  n_crosses_num <- suppressWarnings(as.numeric(n_crosses))
  if (length(n_crosses_num) != 1L || !is.finite(n_crosses_num) || n_crosses_num < 1 ||
      abs(n_crosses_num - round(n_crosses_num)) > 1e-8) {
    ng_stop("n_crosses must be a positive integer")
  }
  n_crosses_int <- as.integer(round(n_crosses_num))
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
                                               phenotypic_covariance = NULL,
                                               genetic_covariance = NULL,
                                               source = ng_multitrait_default_source(),
                                               ...) {
  scored <- ng_add_multitrait_score(
    scores = scores,
    traits = traits,
    method = multitrait_method,
    out_col = out_col,
    strict_thresholds = strict_thresholds,
    threshold_penalty_weight = threshold_penalty_weight,
    phenotypic_covariance = phenotypic_covariance,
    genetic_covariance = genetic_covariance,
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
