ng_frontier_default_source <- function() {
  "validated_5k_parent_grid_2026_05_03"
}

ng_frontier_band_label <- function(min_parents, max_parents) {
  min_parents <- as.numeric(min_parents)
  max_parents <- as.numeric(max_parents)
  if (!is.finite(min_parents) && !is.finite(max_parents)) return("*")
  if (!is.finite(min_parents)) return(paste0("<=", as.integer(max_parents)))
  if (!is.finite(max_parents)) return(paste0(as.integer(min_parents), "+"))
  if (identical(as.integer(min_parents), as.integer(max_parents))) {
    return(as.character(as.integer(min_parents)))
  }
  paste(as.integer(min_parents), as.integer(max_parents), sep = "-")
}

ng_frontier_method_family <- function(method) {
  method <- as.character(method)[[1]]
  if (!nzchar(method) || is.na(method)) return(NA_character_)
  if (grepl("^ng_frontier_policy", method)) return("frontier_policy")
  if (grepl("^ng_recomb_gebv", method)) return("recomb_gebv")
  if (grepl("^ng_pmv_.*balanced_ocs", method)) return("pmv_balanced")
  if (grepl("^ng_meta_portfolio", method)) return("portfolio")
  if (grepl("^ng_meta_selector", method)) return("meta_selector")
  if (grepl("^ng_meta_router", method)) return("meta_router")
  if (grepl("^alphamate", method)) return("alphamate")
  if (grepl("^popvar_uc", method)) return("popvar_uc")
  if (grepl("^simple_usefa", method)) return("simple_usefa")
  if (grepl("^var_simple", method)) return("var_simple")
  "other"
}

ng_frontier_split_candidates <- function(x) {
  x <- unlist(strsplit(as.character(x), "\\|", perl = TRUE), use.names = FALSE)
  x <- trimws(x)
  unique(x[nzchar(x) & !is.na(x)])
}

ng_frontier_default_policy <- function(source = ng_frontier_default_source()) {
  # NOTE (v0.1.0): the validated 2026-05-03 frontier was built on the unit-mixed
  # `etk_dh_pmv_scaled_var_blend_cal` score, which was removed in v0.1.0. The
  # `ng_pmv_blend_balanced_ocs10_lps2` slot below is the in-units replacement that
  # uses `etk_dh_pmv_var_blend_cal`. The frontier dispatch is RETAINED so existing
  # users get a sensible default per parent-count band, but the headline claim
  # "ng_frontier_policy_ocs10_lps2 beat AlphaMate in 7/7 bands" MUST be
  # re-validated against the corrected metrics before being repromoted; see
  # VALIDATED_STATE.md.
  primary <- c(
    "ng_recomb_gebv_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2",
    "ng_pmv_blend_balanced_ocs10_lps2",
    "ng_meta_router_ocs10_lps2",
    "popvar_uc_ocs10_lps1",
    "ng_meta_selector_ocs10_lps2",
    "ng_meta_portfolio_ocs10_lps2"
  )
  candidates <- c(
    "ng_recomb_gebv_ocs10_lps2|ng_meta_router_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2|ng_meta_portfolio_ocs10_lps2|ng_meta_router_ocs10_lps2",
    "ng_pmv_blend_balanced_ocs10_lps2|ng_meta_router_ocs10_lps2",
    "ng_meta_router_ocs10_lps2|ng_meta_portfolio_ocs10_lps2",
    "popvar_uc_ocs10_lps1|simple_usefa_ocs10_lps1|ng_meta_router_ocs10_lps2",
    "ng_meta_selector_ocs10_lps2|ng_meta_portfolio_ocs10_lps2|ng_meta_router_ocs10_lps2",
    "ng_meta_portfolio_ocs10_lps2|ng_meta_router_ocs10_lps2"
  )
  out <- data.frame(
    min_parents = c(1, 26, 36, 46, 56, 66, 76),
    max_parents = c(25, 35, 45, 55, 65, 75, Inf),
    validated_parent_count = c(20L, 30L, 40L, 50L, 60L, 70L, 80L),
    primary_method = primary,
    candidate_methods = candidates,
    source = source,
    reason = "validated_parent_size_frontier",
    stringsAsFactors = FALSE
  )
  out$band_label <- vapply(
    seq_len(nrow(out)),
    function(i) ng_frontier_band_label(out$min_parents[[i]], out$max_parents[[i]]),
    character(1)
  )
  out$primary_family <- vapply(out$primary_method, ng_frontier_method_family, character(1))
  out
}

ng_frontier_parse_band <- function(x) {
  x <- trimws(as.character(x)[[1]])
  if (!nzchar(x)) ng_stop("frontier policy band is empty")
  lower <- tolower(x)
  if (lower %in% c("*", "all", "any")) {
    return(list(min_parents = -Inf, max_parents = Inf, label = "*"))
  }
  if (grepl("^>=[[:space:]]*[0-9]+$", x, perl = TRUE)) {
    n <- as.integer(gsub("[^0-9]", "", x))
    return(list(min_parents = n, max_parents = Inf, label = paste0(n, "+")))
  }
  if (grepl("^[0-9]+[[:space:]]*[+]$", x, perl = TRUE)) {
    n <- as.integer(gsub("[^0-9]", "", x))
    return(list(min_parents = n, max_parents = Inf, label = paste0(n, "+")))
  }
  if (grepl("^<=[[:space:]]*[0-9]+$", x, perl = TRUE)) {
    n <- as.integer(gsub("[^0-9]", "", x))
    return(list(min_parents = -Inf, max_parents = n, label = paste0("<=", n)))
  }
  if (grepl("^[0-9]+[[:space:]]*-[[:space:]]*[0-9]+$", x, perl = TRUE)) {
    parts <- as.integer(unlist(strsplit(gsub("[[:space:]]", "", x), "-", fixed = TRUE)))
    if (length(parts) != 2L || any(!is.finite(parts)) || parts[[1]] > parts[[2]]) {
      ng_stop("invalid frontier policy band: ", x)
    }
    return(list(
      min_parents = parts[[1]],
      max_parents = parts[[2]],
      label = paste(parts[[1]], parts[[2]], sep = "-")
    ))
  }
  if (grepl("^[0-9]+$", x, perl = TRUE)) {
    n <- as.integer(x)
    return(list(min_parents = n, max_parents = n, label = as.character(n)))
  }
  ng_stop("invalid frontier policy band: ", x)
}

ng_frontier_policy_from_spec <- function(spec, source = "custom_spec") {
  spec <- paste(as.character(spec), collapse = ",")
  entries <- unlist(strsplit(spec, "[,;\n]+", perl = TRUE), use.names = FALSE)
  entries <- trimws(entries)
  entries <- entries[nzchar(entries) & !is.na(entries)]
  if (!length(entries)) ng_stop("frontier policy spec is empty")

  rows <- lapply(entries, function(entry) {
    sep <- regexpr("=|:", entry, perl = TRUE)
    if (sep[[1]] < 1L) ng_stop("frontier policy entry needs '=' or ':': ", entry)
    lhs <- trimws(substr(entry, 1L, sep[[1]] - 1L))
    rhs <- trimws(substr(entry, sep[[1]] + 1L, nchar(entry)))
    candidates <- ng_frontier_split_candidates(rhs)
    if (!length(candidates)) ng_stop("frontier policy entry has no methods: ", entry)
    band <- ng_frontier_parse_band(lhs)
    data.frame(
      min_parents = band$min_parents,
      max_parents = band$max_parents,
      validated_parent_count = NA_integer_,
      primary_method = candidates[[1]],
      candidate_methods = paste(candidates, collapse = "|"),
      source = source,
      reason = "custom_frontier_policy",
      band_label = band$label,
      primary_family = ng_frontier_method_family(candidates[[1]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ng_frontier_normalize_policy <- function(policy, source = ng_frontier_default_source()) {
  if (is.null(policy)) return(ng_frontier_default_policy(source = source))
  if (is.character(policy) && length(policy) == 1L) {
    policy <- trimws(policy)
    if (!nzchar(policy)) return(ng_frontier_default_policy(source = source))
    return(ng_frontier_policy_from_spec(policy))
  }
  required <- c("min_parents", "max_parents", "primary_method", "candidate_methods")
  missing <- setdiff(required, names(policy))
  if (length(missing)) ng_stop("frontier policy missing columns: ", paste(missing, collapse = ", "))
  policy <- as.data.frame(policy, stringsAsFactors = FALSE)
  if (!("source" %in% names(policy))) policy$source <- source
  if (!("reason" %in% names(policy))) policy$reason <- "frontier_policy"
  if (!("band_label" %in% names(policy))) {
    policy$band_label <- vapply(
      seq_len(nrow(policy)),
      function(i) ng_frontier_band_label(policy$min_parents[[i]], policy$max_parents[[i]]),
      character(1)
    )
  }
  if (!("primary_family" %in% names(policy))) {
    policy$primary_family <- vapply(policy$primary_method, ng_frontier_method_family, character(1))
  }
  policy
}

ng_frontier_policy_select <- function(n_parents,
                                      policy = NULL,
                                      available_methods = NULL,
                                      fallback_method = "ng_meta_router_ocs10_lps2",
                                      source = ng_frontier_default_source()) {
  n_parents <- suppressWarnings(as.integer(n_parents[[1]]))
  policy <- ng_frontier_normalize_policy(policy, source = source)
  available_methods <- if (is.null(available_methods)) {
    NULL
  } else {
    unique(trimws(as.character(available_methods)))
  }
  if (!is.finite(n_parents) || n_parents < 1L) {
    candidates <- unique(c(fallback_method))
    selected <- candidates[[1]]
    return(list(
      method = selected,
      family = ng_frontier_method_family(selected),
      primary_method = NA_character_,
      primary_family = NA_character_,
      candidate_methods = candidates,
      source = source,
      reason = "invalid_parent_count",
      band_label = NA_character_,
      is_fallback = TRUE
    ))
  }

  idx <- which(n_parents >= policy$min_parents & n_parents <= policy$max_parents)
  if (!length(idx)) {
    candidates <- unique(c(fallback_method))
    selected <- candidates[[1]]
    return(list(
      method = selected,
      family = ng_frontier_method_family(selected),
      primary_method = NA_character_,
      primary_family = NA_character_,
      candidate_methods = candidates,
      source = source,
      reason = "no_matching_band",
      band_label = NA_character_,
      is_fallback = TRUE
    ))
  }

  width <- policy$max_parents[idx] - policy$min_parents[idx]
  width[!is.finite(width)] <- Inf
  row <- policy[idx[order(width)][[1]], , drop = FALSE]
  candidates <- ng_frontier_split_candidates(row$candidate_methods[[1]])
  fallback_method <- trimws(as.character(fallback_method)[[1]])
  if (nzchar(fallback_method) && !(fallback_method %in% candidates)) {
    candidates <- c(candidates, fallback_method)
  }
  candidates <- candidates[!grepl("^ng_frontier_policy", candidates)]
  if (!length(candidates)) ng_stop("frontier policy has no non-recursive candidate methods")

  selected <- candidates[[1]]
  if (!is.null(available_methods)) {
    hits <- candidates[candidates %in% available_methods]
    if (length(hits)) {
      selected <- hits[[1]]
    } else if (nzchar(fallback_method) && fallback_method %in% available_methods) {
      selected <- fallback_method
    }
  }

  list(
    method = selected,
    family = ng_frontier_method_family(selected),
    primary_method = row$primary_method[[1]],
    primary_family = row$primary_family[[1]],
    candidate_methods = candidates,
    source = row$source[[1]],
    reason = row$reason[[1]],
    band_label = row$band_label[[1]],
    is_fallback = !identical(selected, candidates[[1]])
  )
}
