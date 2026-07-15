ng_poly_subgenome_names <- function(x) {
  names_x <- names(x)
  if (is.null(names_x) || anyNA(names_x) || any(!nzchar(trimws(names_x)))) {
    ng_stop("geno_by_subgenome must be a named list")
  }
  names_x <- trimws(as.character(names_x))
  if (anyDuplicated(names_x)) ng_stop("subgenome names must be unique")
  names_x
}

ng_poly_subgenome_as_dosage_list <- function(geno_by_subgenome, name = "geno_by_subgenome") {
  if (!is.list(geno_by_subgenome) || !length(geno_by_subgenome)) {
    ng_stop(name, " must be a non-empty named list")
  }
  subgenomes <- ng_poly_subgenome_names(geno_by_subgenome)
  names(geno_by_subgenome) <- subgenomes

  parent_ids <- NULL
  out <- vector("list", length(geno_by_subgenome))
  names(out) <- subgenomes

  for (sg in subgenomes) {
    raw <- geno_by_subgenome[[sg]]
    raw_row_ids <- rownames(raw)
    raw_col_ids <- colnames(raw)
    if (is.null(raw_row_ids) || anyNA(raw_row_ids) || any(!nzchar(trimws(raw_row_ids)))) {
      ng_stop(name, "$", sg, " must have parent row names")
    }
    if (is.null(raw_col_ids) || anyNA(raw_col_ids) || any(!nzchar(trimws(raw_col_ids)))) {
      ng_stop(name, "$", sg, " must have marker column names")
    }
    mat <- ng_as_numeric_matrix(geno_by_subgenome[[sg]], paste0(name, "$", sg))
    row_ids <- rownames(mat)
    col_ids <- colnames(mat)
    if (anyDuplicated(row_ids)) {
      ng_stop(name, "$", sg, " parent row names must be unique")
    }
    if (anyDuplicated(col_ids)) {
      ng_stop(name, "$", sg, " marker column names must be unique")
    }

    bad <- !is.finite(mat) | mat < 0 | mat > 2 | abs(mat - round(mat)) > 1e-8
    if (any(bad)) {
      idx <- which(bad, arr.ind = TRUE)[1, ]
      ng_stop(name, "$", sg, " has diploid dosage outside 0..2 at row ",
              row_ids[idx[[1]]], ", column ", col_ids[idx[[2]]])
    }

    if (is.null(parent_ids)) {
      parent_ids <- row_ids
    } else if (!identical(parent_ids, row_ids)) {
      ng_stop("all subgenomes must share identical parent row names in the same order")
    }

    storage.mode(mat) <- "double"
    out[[sg]] <- mat
  }

  out
}

ng_poly_subgenome_weights <- function(subgenomes, weights = NULL) {
  if (is.null(weights)) {
    out <- rep.int(1, length(subgenomes))
    names(out) <- subgenomes
    return(out)
  }

  out <- suppressWarnings(as.numeric(weights))
  names(out) <- names(weights)
  if (length(out) != length(subgenomes) || any(!is.finite(out)) || any(out <= 0)) {
    ng_stop("weights must provide one positive finite value per subgenome")
  }
  if (is.null(names(out)) || any(!nzchar(names(out))) || !setequal(names(out), subgenomes)) {
    ng_stop("weights names must match subgenome names")
  }

  out[subgenomes]
}

ng_poly_diploid_parent_relationship <- function(geno) {
  X <- geno - 1
  denom <- sum(diag(stats::var(X)))
  if (!is.finite(denom) || denom <= 0) denom <- ncol(X)

  K <- tcrossprod(X) / denom
  rownames(K) <- rownames(geno)
  colnames(K) <- rownames(geno)
  K
}

ng_poly_subgenome_parent_relationship <- function(geno_by_subgenome, weights = NULL) {
  geno_by_subgenome <- ng_poly_subgenome_as_dosage_list(geno_by_subgenome)
  subgenomes <- names(geno_by_subgenome)
  weights <- ng_poly_subgenome_weights(subgenomes, weights)

  component_K <- lapply(geno_by_subgenome, ng_poly_diploid_parent_relationship)
  K <- Reduce(`+`, Map(function(k, w) k * w, component_K, weights)) / sum(weights)
  rownames(K) <- rownames(geno_by_subgenome[[1]])
  colnames(K) <- rownames(geno_by_subgenome[[1]])
  attr(K, "subgenome_K") <- component_K
  attr(K, "subgenome_weights") <- weights
  K
}

ng_poly_subgenome_as_effects_list <- function(effects_by_subgenome, geno_by_subgenome) {
  if (!is.list(effects_by_subgenome) || !length(effects_by_subgenome)) {
    ng_stop("effects_by_subgenome must be a non-empty named list")
  }
  subgenomes <- names(geno_by_subgenome)
  if (is.null(names(effects_by_subgenome)) || !setequal(names(effects_by_subgenome), subgenomes)) {
    ng_stop("effects_by_subgenome names must match geno_by_subgenome names")
  }

  out <- vector("list", length(subgenomes))
  names(out) <- subgenomes
  for (sg in subgenomes) {
    effects <- suppressWarnings(as.numeric(effects_by_subgenome[[sg]]))
    names(effects) <- names(effects_by_subgenome[[sg]])
    if (length(effects) != ncol(geno_by_subgenome[[sg]]) || any(!is.finite(effects))) {
      ng_stop("effects_by_subgenome$", sg, " must provide one finite value per marker")
    }
    if (is.null(names(effects)) || !identical(names(effects), colnames(geno_by_subgenome[[sg]]))) {
      ng_stop("effects_by_subgenome$", sg, " marker names must match genotype marker columns in order")
    }
    out[[sg]] <- effects
  }

  out
}

ng_poly_validate_candidate_pairs <- function(candidate_pairs, parent_ids) {
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
  if (any(pairs$parent1 == pairs$parent2)) ng_stop("candidate_pairs must not contain self-crosses")
  missing <- setdiff(unique(c(pairs$parent1, pairs$parent2)), parent_ids)
  if (length(missing)) ng_stop("candidate_pairs contains unknown parent IDs: ", paste(missing, collapse = ", "))
  if (!nrow(pairs)) ng_stop("candidate_pairs must contain at least one cross")
  key <- ng_poly4x_pair_key(pairs$parent1, pairs$parent2)
  if (anyDuplicated(key)) ng_stop("candidate_pairs contains duplicate or reciprocal duplicate crosses")

  pairs[, c("parent1", "parent2"), drop = FALSE]
}

ng_poly_add_diagnostics <- function(out,
                                    model_decision = NULL,
                                    validation_source = "phase2a_deterministic_core") {
  if (is.null(model_decision)) {
    model_decision <- ng_poly_model_decision(
      model_family = "allopolyploid_subgenome",
      inheritance_model = "disomic_subgenome",
      supported = TRUE,
      supported_scope = "wheat_like_allopolyploid_disomic",
      reason = "phase2a_direct_subgenome_scoring"
    )
  }
  subgenomes <- attr(out, "subgenome_names")
  if (!length(subgenomes)) subgenomes <- model_decision$subgenome_names

  out$polyploid_model_family <- model_decision$model_family
  out$inheritance_model <- model_decision$inheritance_model
  out$subgenome_count <- length(subgenomes)
  out$subgenome_names <- paste(subgenomes, collapse = "|")
  out$supported_scope <- model_decision$supported_scope
  out$polyploid_validation_source <- validation_source
  out
}

ng_poly_subgenome_score_crosses <- function(geno_by_subgenome,
                                            effects_by_subgenome,
                                            candidate_pairs = NULL,
                                            model_decision = NULL,
                                            selection_prop = 0.10,
                                            weights = NULL,
                                            validation_source = "phase2a_deterministic_core") {
  geno_by_subgenome <- ng_poly_subgenome_as_dosage_list(geno_by_subgenome)
  effects_by_subgenome <- ng_poly_subgenome_as_effects_list(effects_by_subgenome, geno_by_subgenome)
  parent_ids <- rownames(geno_by_subgenome[[1]])
  if (is.null(candidate_pairs)) candidate_pairs <- ng_make_pairs(parent_ids, include_self = FALSE)
  candidate_pairs <- ng_poly_validate_candidate_pairs(candidate_pairs, parent_ids)
  parent_K <- ng_poly_subgenome_parent_relationship(geno_by_subgenome, weights = weights)
  intensity <- ng_selection_intensity(selection_prop)

  subgenomes <- names(geno_by_subgenome)
  parent_values <- lapply(subgenomes, function(sg) {
    value <- as.numeric(geno_by_subgenome[[sg]] %*% effects_by_subgenome[[sg]])
    names(value) <- parent_ids
    value
  })
  names(parent_values) <- subgenomes

  gain <- numeric(nrow(candidate_pairs))
  variance <- numeric(nrow(candidate_pairs))
  for (sg in subgenomes) {
    geno <- geno_by_subgenome[[sg]]
    effects <- effects_by_subgenome[[sg]]
    p1 <- candidate_pairs$parent1
    p2 <- candidate_pairs$parent2
    gain <- gain + (parent_values[[sg]][p1] + parent_values[[sg]][p2]) / 2
    contrast <- sweep(geno[p1, , drop = FALSE] - geno[p2, , drop = FALSE], 2, effects, `*`)
    variance <- variance + rowSums(contrast^2) / 4
  }

  out <- data.frame(
    parent1 = candidate_pairs$parent1,
    parent2 = candidate_pairs$parent2,
    poly_gain = as.numeric(gain),
    poly_var = as.numeric(variance),
    poly_usefulness = as.numeric(gain + intensity * sqrt(pmax(variance, 0))),
    pair_kinship = ng_poly4x_pair_coancestry(parent_K, candidate_pairs),
    stringsAsFactors = FALSE
  )
  attr(out, "parent_K") <- parent_K
  attr(out, "subgenome_names") <- subgenomes
  out <- ng_poly_add_diagnostics(out, model_decision = model_decision, validation_source = validation_source)
  attr(out, "parent_K") <- parent_K
  attr(out, "subgenome_names") <- subgenomes
  out
}

ng_poly_policy_modes <- function() {
  c("gain", "diversity", "ocs")
}

ng_poly_normalize_policy_mode <- function(mode) {
  mode <- ng_poly_trim1(mode)
  if (!nzchar(mode) || !(mode %in% ng_poly_policy_modes())) {
    ng_stop("mode must be one of: ", paste(ng_poly_policy_modes(), collapse = ", "))
  }
  mode
}

ng_poly_add_policy_diagnostics <- function(plan, mode, source_scores) {
  plan$poly_policy_mode <- mode
  plan$poly_policy_family <- unique(source_scores$polyploid_model_family)[[1]]
  plan$poly_policy_scope <- unique(source_scores$supported_scope)[[1]]
  plan$poly_policy_source <- unique(source_scores$polyploid_validation_source)[[1]]

  summary <- attr(plan, "summary")
  if (!is.null(summary)) {
    summary$poly_policy_mode <- mode
    summary$poly_policy_family <- plan$poly_policy_family[[1]]
    summary$poly_policy_scope <- plan$poly_policy_scope[[1]]
    summary$poly_policy_source <- plan$poly_policy_source[[1]]
    attr(plan, "summary") <- summary
  }
  plan
}

ng_poly_policy <- function(scores,
                           n_crosses,
                           mode = "gain",
                           parent_K = NULL,
                           max_crosses_per_parent = NULL,
                           method = "auto",
                           local_iter = 1000L,
                           ocs_iter = 5L,
                           ...) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  mode <- ng_poly_normalize_policy_mode(mode)
  if (is.null(parent_K)) parent_K <- attr(scores, "parent_K")

  if (identical(mode, "gain")) {
    plan <- ng_poly4x_select_topn(scores, n_crosses, "poly_usefulness", "ng_poly_gain_topn")
    return(ng_poly_add_policy_diagnostics(plan, mode, scores))
  }

  if (is.null(parent_K)) ng_stop("parent_K is required for OCS polyploid policy modes")
  if (is.null(max_crosses_per_parent)) {
    max_crosses_per_parent <- if (identical(mode, "diversity")) 3L else 4L
  }
  plan <- ng_optimize_mating_plan(
    scores = scores,
    n_crosses = n_crosses,
    gain_col = "poly_usefulness",
    parent_K = parent_K,
    max_crosses_per_parent = max_crosses_per_parent,
    lambda_group = if (identical(mode, "diversity")) 1.0 else 0.5,
    lambda_parent_use = if (identical(mode, "diversity")) 2.0 else 1.0,
    method = method,
    local_iter = local_iter,
    ocs_iter = ocs_iter,
    ...
  )
  plan$poly_method <- paste0("ng_poly_", mode)

  summary <- attr(plan, "summary")
  if (!is.null(summary)) {
    summary$poly_method <- paste0("ng_poly_", mode)
    summary$poly_gain_col <- "poly_usefulness"
    attr(plan, "summary") <- summary
  }

  ng_poly_add_policy_diagnostics(plan, mode, scores)
}
