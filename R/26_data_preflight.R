ng_preflight_timestamp <- function(generated_at = Sys.time()) {
  format(as.POSIXct(generated_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

ng_preflight_empty_issues <- function() {
  data.frame(
    id = character(),
    severity = character(),
    table = character(),
    field = character(),
    message = character(),
    count = integer(),
    stringsAsFactors = FALSE
  )
}

ng_preflight_issue <- function(id, severity, table, field, message, count = 1L) {
  data.frame(
    id = as.character(id),
    severity = as.character(severity),
    table = as.character(table),
    field = as.character(field),
    message = as.character(message),
    count = as.integer(count),
    stringsAsFactors = FALSE
  )
}

ng_preflight_add_issue <- function(issues, id, severity, table, field, message, count = 1L) {
  rbind(issues, ng_preflight_issue(id, severity, table, field, message, count))
}

ng_preflight_id_column <- function(x, candidates) {
  names_x <- names(x)
  hit <- match(tolower(candidates), tolower(names_x), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) names_x[[hit[[1L]]]] else NULL
}

ng_preflight_table_ids <- function(x, candidates = c("parent", "parent_id", "id", "name", "line", "entry")) {
  if (is.null(x)) return(character())
  if (is.matrix(x)) return(as.character(rownames(x)))
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- ng_preflight_id_column(x, candidates)
  if (!is.null(id_col)) return(trimws(as.character(x[[id_col]])))
  rn <- rownames(x)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(x))))) return(as.character(rn))
  character()
}

ng_preflight_marker_columns <- function(geno) {
  if (is.null(geno)) return(character())
  if (is.matrix(geno)) return(colnames(geno))
  geno <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- ng_preflight_id_column(geno, c("parent", "parent_id", "id", "name", "line", "entry"))
  cols <- names(geno)
  if (!is.null(id_col)) cols <- cols[seq_along(cols) != match(id_col, names(geno))]
  cols
}

ng_preflight_table_summary <- function(name, x) {
  present <- !is.null(x)
  data.frame(
    table = name,
    present = present,
    rows = if (present) nrow(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)) else 0L,
    columns = if (present) ncol(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)) else 0L,
    stringsAsFactors = FALSE
  )
}

ng_preflight_duplicate_values <- function(x) {
  x <- trimws(as.character(x))
  x <- x[nzchar(x) & !is.na(x)]
  unique(x[duplicated(x)])
}

ng_preflight_pair_key <- function(parent1, parent2) {
  parent1 <- trimws(as.character(parent1))
  parent2 <- trimws(as.character(parent2))
  lo <- ifelse(parent1 <= parent2, parent1, parent2)
  hi <- ifelse(parent1 <= parent2, parent2, parent1)
  paste(lo, hi, sep = "\r")
}

ng_preflight_match_arg <- function(value, choices, name) {
  if (is.null(value)) value <- choices[[1L]]
  value <- trimws(tolower(as.character(value[[1L]])))
  hit <- match(value, choices, nomatch = 0L)
  if (!hit) ng_stop(name, " must be one of: ", paste(choices, collapse = ", "))
  choices[[hit]]
}

ng_preflight_subset_table_by_ids <- function(x, keep_ids,
                                             candidates = c("parent", "parent_id", "id", "name", "line", "entry")) {
  if (is.null(x)) return(NULL)
  keep_ids <- trimws(as.character(keep_ids))
  if (is.matrix(x)) {
    ids <- rownames(x)
    if (is.null(ids)) return(x)
    return(x[trimws(as.character(ids)) %in% keep_ids, , drop = FALSE])
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- ng_preflight_id_column(out, candidates)
  if (!is.null(id_col)) {
    return(out[trimws(as.character(out[[id_col]])) %in% keep_ids, , drop = FALSE])
  }
  rn <- rownames(out)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(out))))) {
    return(out[trimws(as.character(rn)) %in% keep_ids, , drop = FALSE])
  }
  out
}

ng_preflight_filter_pairs_by_ids <- function(candidate_pairs, keep_ids) {
  if (is.null(candidate_pairs)) return(NULL)
  out <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("parent1", "parent2") %in% names(out))) return(out)
  keep_ids <- trimws(as.character(keep_ids))
  p1 <- trimws(as.character(out$parent1))
  p2 <- trimws(as.character(out$parent2))
  out[p1 %in% keep_ids & p2 %in% keep_ids, , drop = FALSE]
}

ng_preflight_subset_parent_K <- function(parent_kinship, keep_ids) {
  if (is.null(parent_kinship)) return(NULL)
  rn <- rownames(parent_kinship)
  cn <- colnames(parent_kinship)
  if (is.null(rn) || is.null(cn)) return(parent_kinship)
  keep_ids <- trimws(as.character(keep_ids))
  keep <- intersect(keep_ids, intersect(trimws(as.character(rn)), trimws(as.character(cn))))
  parent_kinship[keep, keep, drop = FALSE]
}

ng_preflight_nrow <- function(x) {
  if (is.null(x)) return(0L)
  nrow(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
}

ng_preflight_rows_removed <- function(original, cleaned, table_name) {
  data.frame(
    table = table_name,
    rows_removed = as.integer(ng_preflight_nrow(original) - ng_preflight_nrow(cleaned)),
    stringsAsFactors = FALSE
  )
}

ng_preflight_duplicate_removal_plan <- function(clusters, geno_ids, dosage) {
  empty_removed <- data.frame(
    cluster_id = character(),
    kept_parent = character(),
    removed_parent = character(),
    kept_missing_prop = numeric(),
    removed_missing_prop = numeric(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  empty_kept <- data.frame(
    cluster_id = character(),
    kept_parent = character(),
    cluster_size = integer(),
    stringsAsFactors = FALSE
  )
  if (!length(clusters)) {
    return(list(kept_parents = empty_kept, removed_parents = empty_removed, removed_ids = character()))
  }
  geno_ids <- trimws(as.character(geno_ids))
  miss <- rowMeans(!is.finite(dosage))
  miss[!is.finite(miss)] <- 1
  names(miss) <- geno_ids
  kept_rows <- list()
  removed_rows <- list()
  kept_i <- 0L
  removed_i <- 0L
  for (i in seq_along(clusters)) {
    members <- trimws(as.character(clusters[[i]]))
    members <- members[members %in% geno_ids]
    if (!length(members)) next
    original_order <- match(members, geno_ids)
    ord <- order(miss[members], original_order, members, na.last = TRUE)
    kept <- members[[ord[[1L]]]]
    cluster_id <- sprintf("dup_cluster_%03d", i)
    kept_i <- kept_i + 1L
    kept_rows[[kept_i]] <- data.frame(
      cluster_id = cluster_id,
      kept_parent = kept,
      cluster_size = length(members),
      stringsAsFactors = FALSE
    )
    removed <- members[ord[-1L]]
    if (length(removed)) {
      for (rm in removed) {
        removed_i <- removed_i + 1L
        removed_rows[[removed_i]] <- data.frame(
          cluster_id = cluster_id,
          kept_parent = kept,
          removed_parent = rm,
          kept_missing_prop = as.numeric(miss[[kept]]),
          removed_missing_prop = as.numeric(miss[[rm]]),
          reason = "kept lowest missingness; ties keep first input row",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  kept_parents <- if (length(kept_rows)) do.call(rbind, kept_rows) else empty_kept
  removed_parents <- if (length(removed_rows)) do.call(rbind, removed_rows) else empty_removed
  list(
    kept_parents = kept_parents,
    removed_parents = removed_parents,
    removed_ids = removed_parents$removed_parent
  )
}

ng_preflight_remove_putative_duplicate_tables <- function(geno,
                                                          phenotype,
                                                          candidate_pairs,
                                                          trait_spec,
                                                          marker_map,
                                                          parent_kinship,
                                                          geno_ids,
                                                          dosage,
                                                          putative_duplicates) {
  plan <- ng_preflight_duplicate_removal_plan(putative_duplicates$clusters, geno_ids, dosage)
  keep_ids <- setdiff(trimws(as.character(geno_ids)), plan$removed_ids)
  cleaned_tables <- list(
    geno = ng_preflight_subset_table_by_ids(geno, keep_ids),
    phenotype = ng_preflight_subset_table_by_ids(phenotype, keep_ids),
    candidate_pairs = ng_preflight_filter_pairs_by_ids(candidate_pairs, keep_ids),
    trait_spec = trait_spec,
    marker_map = marker_map,
    parent_kinship = ng_preflight_subset_parent_K(parent_kinship, keep_ids)
  )
  rows_removed <- do.call(rbind, list(
    ng_preflight_rows_removed(geno, cleaned_tables$geno, "geno"),
    ng_preflight_rows_removed(phenotype, cleaned_tables$phenotype, "phenotype"),
    ng_preflight_rows_removed(candidate_pairs, cleaned_tables$candidate_pairs, "candidate_pairs"),
    ng_preflight_rows_removed(parent_kinship, cleaned_tables$parent_kinship, "parent_kinship")
  ))
  list(
    cleaned_tables = cleaned_tables,
    cleaning = list(
      action = "remove",
      kept_parents = plan$kept_parents,
      removed_parents = plan$removed_parents,
      rows_removed = rows_removed
    )
  )
}

ng_preflight_input_tables <- function(geno = NULL,
                                      phenotype = NULL,
                                      candidate_pairs = NULL,
                                      trait_spec = NULL,
                                      marker_map = NULL,
                                      parent_kinship = NULL,
                                      ploidy = NULL,
                                      putative_duplicate_check = FALSE,
                                      duplicate_threshold = 0.995,
                                      duplicate_maf_min = 0,
                                      duplicate_max_missing_prop = 1,
                                      duplicate_min_compared_markers = 100L,
                                      putative_duplicate_action = c("report", "remove"),
                                      putative_duplicate_return_similarity = FALSE,
                                      generated_at = Sys.time()) {
  putative_duplicate_action <- ng_preflight_match_arg(
    putative_duplicate_action, c("report", "remove"), "putative_duplicate_action"
  )
  issues <- ng_preflight_empty_issues()
  putative_duplicates <- NULL
  cleaned_tables <- NULL
  cleaning <- NULL
  tables <- do.call(rbind, list(
    ng_preflight_table_summary("geno", geno),
    ng_preflight_table_summary("phenotype", phenotype),
    ng_preflight_table_summary("candidate_pairs", candidate_pairs),
    ng_preflight_table_summary("trait_spec", trait_spec),
    ng_preflight_table_summary("marker_map", marker_map),
    ng_preflight_table_summary("parent_kinship", parent_kinship)
  ))

  geno_ids <- ng_preflight_table_ids(geno)
  marker_cols <- ng_preflight_marker_columns(geno)
  if (!is.null(geno)) {
    if (!length(geno_ids)) {
      issues <- ng_preflight_add_issue(issues, "missing_parent_id_column", "blocker", "geno", "parent", "Genotype table needs parent IDs in row names or a parent/id column")
    }
    dup <- ng_preflight_duplicate_values(geno_ids)
    if (length(dup)) {
      issues <- ng_preflight_add_issue(issues, "duplicate_parent_ids", "blocker", "geno", "parent", paste("Duplicate parent IDs:", paste(dup, collapse = ", ")), length(dup))
    }
    dup_markers <- ng_preflight_duplicate_values(marker_cols)
    if (length(dup_markers)) {
      issues <- ng_preflight_add_issue(issues, "duplicate_marker_ids", "blocker", "geno", "marker", paste("Duplicate marker IDs:", paste(dup_markers, collapse = ", ")), length(dup_markers))
    }
    if (!is.null(ploidy) && length(marker_cols)) {
      ploidy <- suppressWarnings(as.numeric(ploidy[[1L]]))
      g <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
      marker_idx <- match(marker_cols, names(g))
      dosage <- suppressWarnings(as.matrix(g[, marker_idx, drop = FALSE]))
      storage.mode(dosage) <- "double"
      bad <- !is.finite(dosage) | dosage < 0 | dosage > ploidy | abs(dosage - round(dosage)) > 1e-8
      if (is.finite(ploidy) && any(bad, na.rm = TRUE)) {
        issues <- ng_preflight_add_issue(issues, "invalid_dosage_range", "blocker", "geno", "markers", paste0("Dosage values must be integer values within 0..", ploidy), sum(bad, na.rm = TRUE))
      }
    }
    if (isTRUE(putative_duplicate_check) && length(geno_ids) && !anyDuplicated(geno_ids) && length(marker_cols)) {
      g <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
      marker_idx <- match(marker_cols, names(g))
      dosage <- suppressWarnings(as.matrix(g[, marker_idx, drop = FALSE]))
      storage.mode(dosage) <- "double"
      rownames(dosage) <- geno_ids
      putative_duplicates <- tryCatch(
        ng_detect_putative_duplicates(
          dosage,
          ids = geno_ids,
          duplicate_threshold = duplicate_threshold,
          maf_min = duplicate_maf_min,
          max_missing_prop = duplicate_max_missing_prop,
          min_compared_markers = duplicate_min_compared_markers,
          return_similarity = isTRUE(putative_duplicate_return_similarity)
        ),
        error = function(e) {
          issues <<- ng_preflight_add_issue(
            issues,
            "putative_duplicate_check_failed",
            "warning",
            "geno",
            "markers",
            paste("Putative duplicate check could not be completed:", conditionMessage(e)),
            1L
          )
          NULL
        }
      )
      if (!is.null(putative_duplicates) && nrow(putative_duplicates$pairs)) {
        examples <- paste(
          utils::head(paste(putative_duplicates$pairs$parent1,
                            putative_duplicates$pairs$parent2,
                            sep = "/"), 4L),
          collapse = ", "
        )
        if (identical(putative_duplicate_action, "remove")) {
          cleanup <- ng_preflight_remove_putative_duplicate_tables(
            geno = geno,
            phenotype = phenotype,
            candidate_pairs = candidate_pairs,
            trait_spec = trait_spec,
            marker_map = marker_map,
            parent_kinship = parent_kinship,
            geno_ids = geno_ids,
            dosage = dosage,
            putative_duplicates = putative_duplicates
          )
          cleaned_tables <- cleanup$cleaned_tables
          cleaning <- list(putative_duplicates = cleanup$cleaning)
          geno <- cleaned_tables$geno
          phenotype <- cleaned_tables$phenotype
          candidate_pairs <- cleaned_tables$candidate_pairs
          parent_kinship <- cleaned_tables$parent_kinship
          geno_ids <- ng_preflight_table_ids(geno)
          marker_cols <- ng_preflight_marker_columns(geno)
          n_removed <- nrow(cleanup$cleaning$removed_parents)
          issues <- ng_preflight_add_issue(
            issues,
            "putative_duplicate_genotypes_removed",
            "warning",
            "geno",
            "markers",
            paste0(
              "Putative duplicate genotype profiles exceeded IBS threshold ",
              duplicate_threshold,
              "; pairs=", nrow(putative_duplicates$pairs),
              ", clusters=", length(putative_duplicates$clusters),
              ", removed_parents=", n_removed,
              ", examples: ", examples
            ),
            n_removed
          )
        } else {
          issues <- ng_preflight_add_issue(
            issues,
            "putative_duplicate_genotypes",
            "warning",
            "geno",
            "markers",
            paste0(
              "Putative duplicate genotype profiles exceeded IBS threshold ",
              duplicate_threshold,
              "; pairs=", nrow(putative_duplicates$pairs),
              ", clusters=", length(putative_duplicates$clusters),
              ", examples: ", examples
            ),
            nrow(putative_duplicates$pairs)
          )
        }
      } else if (identical(putative_duplicate_action, "remove")) {
        cleaned_tables <- list(
          geno = geno,
          phenotype = phenotype,
          candidate_pairs = candidate_pairs,
          trait_spec = trait_spec,
          marker_map = marker_map,
          parent_kinship = parent_kinship
        )
        cleaning <- list(
          putative_duplicates = list(
            action = "remove",
            kept_parents = data.frame(cluster_id = character(), kept_parent = character(), cluster_size = integer(), stringsAsFactors = FALSE),
            removed_parents = data.frame(cluster_id = character(), kept_parent = character(), removed_parent = character(), kept_missing_prop = numeric(), removed_missing_prop = numeric(), reason = character(), stringsAsFactors = FALSE),
            rows_removed = data.frame(table = c("geno", "phenotype", "candidate_pairs", "parent_kinship"), rows_removed = 0L, stringsAsFactors = FALSE)
          )
        )
      }
    }
  }

  pheno_ids <- ng_preflight_table_ids(phenotype)
  if (!is.null(phenotype)) {
    p <- as.data.frame(phenotype, stringsAsFactors = FALSE, check.names = FALSE)
    parent_col <- ng_preflight_id_column(p, c("parent", "parent_id", "id", "name", "line", "entry"))
    trait_col <- ng_preflight_id_column(p, c("trait", "trait_name", "variable"))
    if (!is.null(parent_col)) {
      key <- trimws(as.character(p[[parent_col]]))
      if (!is.null(trait_col)) key <- paste(key, trimws(as.character(p[[trait_col]])), sep = "\r")
      dup <- unique(key[duplicated(key) & nzchar(key)])
      if (length(dup)) {
        issues <- ng_preflight_add_issue(issues, "duplicate_phenotype_rows", "warning", "phenotype", parent_col, "Phenotype table has repeated parent/trait records", length(dup))
      }
    }
    if (length(geno_ids) && length(pheno_ids)) {
      missing_in_geno <- setdiff(pheno_ids, geno_ids)
      missing_in_pheno <- setdiff(geno_ids, pheno_ids)
      if (length(missing_in_geno) || length(missing_in_pheno)) {
        issues <- ng_preflight_add_issue(
          issues,
          "genotype_phenotype_id_mismatch",
          "blocker",
          "phenotype",
          "parent",
          paste0("Genotype/phenotype IDs differ; phenotype-only=", length(missing_in_geno), ", genotype-only=", length(missing_in_pheno)),
          length(unique(c(missing_in_geno, missing_in_pheno)))
        )
      }
    }
  }

  if (!is.null(candidate_pairs)) {
    pairs <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("parent1", "parent2") %in% names(pairs))) {
      issues <- ng_preflight_add_issue(issues, "missing_candidate_pair_columns", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs must contain parent1 and parent2")
    } else {
      p1 <- trimws(as.character(pairs$parent1))
      p2 <- trimws(as.character(pairs$parent2))
      missing_pair_ids <- !nzchar(p1) | !nzchar(p2) | is.na(p1) | is.na(p2)
      if (any(missing_pair_ids)) {
        issues <- ng_preflight_add_issue(issues, "missing_candidate_parent", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs contain missing parent IDs", sum(missing_pair_ids))
      }
      self <- p1 == p2 & nzchar(p1)
      if (any(self)) {
        issues <- ng_preflight_add_issue(issues, "self_cross", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs contain self-crosses", sum(self))
      }
      if (length(geno_ids)) {
        unknown <- setdiff(unique(c(p1, p2)), geno_ids)
        unknown <- unknown[nzchar(unknown) & !is.na(unknown)]
        if (length(unknown)) {
          issues <- ng_preflight_add_issue(issues, "unknown_candidate_parent", "blocker", "candidate_pairs", "parent1,parent2", paste("Unknown candidate parents:", paste(unknown, collapse = ", ")), length(unknown))
        }
      }
      direct_key <- paste(p1, p2, sep = "\r")
      dup_direct <- unique(direct_key[duplicated(direct_key) & nzchar(direct_key)])
      if (length(dup_direct)) {
        issues <- ng_preflight_add_issue(issues, "duplicate_candidate_crosses", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs contain repeated parent1/parent2 rows", length(dup_direct))
      }
      canonical <- ng_preflight_pair_key(p1, p2)
      dup_canonical <- unique(canonical[duplicated(canonical) & nzchar(canonical)])
      if (length(dup_canonical) && !length(dup_direct)) {
        issues <- ng_preflight_add_issue(issues, "duplicate_candidate_crosses", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs contain repeated unordered crosses", length(dup_canonical))
      }
      reciprocal_count <- length(setdiff(dup_canonical, dup_direct))
      if (reciprocal_count > 0L || length(dup_canonical) > length(dup_direct)) {
        issues <- ng_preflight_add_issue(issues, "reciprocal_candidate_crosses", "blocker", "candidate_pairs", "parent1,parent2", "Candidate pairs contain reciprocal duplicate crosses", max(1L, reciprocal_count))
      }
    }
  }

  if (!is.null(trait_spec)) {
    spec <- as.data.frame(trait_spec, stringsAsFactors = FALSE, check.names = FALSE)
    trait_col <- ng_preflight_id_column(spec, c("trait", "trait_name", "name"))
    if (is.null(trait_col)) {
      issues <- ng_preflight_add_issue(issues, "missing_trait_column", "blocker", "trait_spec", "trait", "Trait specification must contain a trait column")
    } else {
      dup_traits <- ng_preflight_duplicate_values(spec[[trait_col]])
      if (length(dup_traits)) {
        issues <- ng_preflight_add_issue(issues, "duplicate_trait_names", "blocker", "trait_spec", trait_col, paste("Duplicate trait names:", paste(dup_traits, collapse = ", ")), length(dup_traits))
      }
    }
    direction_col <- ng_preflight_id_column(spec, c("direction", "trait_direction"))
    if (!is.null(direction_col)) {
      valid <- c("max", "maximize", "maximise", "increase", "higher", "high", "positive", "+",
                 "min", "minimize", "minimise", "decrease", "lower", "low", "negative", "-")
      direction <- trimws(tolower(as.character(spec[[direction_col]])))
      bad <- nzchar(direction) & !(direction %in% valid)
      if (any(bad, na.rm = TRUE)) {
        issues <- ng_preflight_add_issue(issues, "invalid_trait_direction", "blocker", "trait_spec", direction_col, paste("Invalid trait directions:", paste(unique(direction[bad]), collapse = ", ")), sum(bad, na.rm = TRUE))
      }
    }
    min_col <- ng_preflight_id_column(spec, c("min_value", "minimum", "min"))
    max_col <- ng_preflight_id_column(spec, c("max_value", "maximum", "max"))
    if (!is.null(min_col) && !is.null(max_col)) {
      min_value <- suppressWarnings(as.numeric(spec[[min_col]]))
      max_value <- suppressWarnings(as.numeric(spec[[max_col]]))
      bad <- is.finite(min_value) & is.finite(max_value) & min_value > max_value
      if (any(bad)) {
        issues <- ng_preflight_add_issue(issues, "invalid_trait_threshold", "blocker", "trait_spec", "min_value,max_value", "Trait minimum thresholds cannot exceed maximum thresholds", sum(bad))
      }
    }
  }

  if (!is.null(marker_map) && length(marker_cols)) {
    map <- as.data.frame(marker_map, stringsAsFactors = FALSE, check.names = FALSE)
    marker_col <- ng_preflight_id_column(map, c("marker", "marker_id", "snp_code", "id"))
    if (is.null(marker_col)) {
      issues <- ng_preflight_add_issue(issues, "missing_marker_map_column", "blocker", "marker_map", "marker", "Marker map must contain a marker column")
    } else {
      map_markers <- trimws(as.character(map[[marker_col]]))
      dup_map <- ng_preflight_duplicate_values(map_markers)
      if (length(dup_map)) {
        issues <- ng_preflight_add_issue(issues, "duplicate_marker_map_entries", "blocker", "marker_map", marker_col, paste("Duplicate marker-map entries:", paste(dup_map, collapse = ", ")), length(dup_map))
      }
      missing_markers <- setdiff(unique(marker_cols), map_markers)
      if (length(missing_markers)) {
        issues <- ng_preflight_add_issue(issues, "missing_marker_map_entries", "blocker", "marker_map", marker_col, paste("Marker map is missing genotype markers:", paste(missing_markers, collapse = ", ")), length(missing_markers))
      }
    }
  }

  if (!is.null(parent_kinship)) {
    rn <- rownames(parent_kinship)
    cn <- colnames(parent_kinship)
    dup <- unique(c(ng_preflight_duplicate_values(rn), ng_preflight_duplicate_values(cn)))
    if (length(dup)) {
      issues <- ng_preflight_add_issue(issues, "parent_relationship_duplicate_ids", "blocker", "parent_kinship", "rownames,colnames", paste("Parent relationship matrix has duplicate IDs:", paste(dup, collapse = ", ")), length(dup))
    }
    if (length(geno_ids) && (!setequal(geno_ids, rn) || !setequal(geno_ids, cn))) {
      issues <- ng_preflight_add_issue(issues, "parent_relationship_id_mismatch", "blocker", "parent_kinship", "rownames,colnames", "Parent relationship matrix IDs must match genotype parent IDs", 1L)
    }
  }

  blockers <- sum(issues$severity == "blocker")
  warnings <- sum(issues$severity == "warning")
  status <- if (blockers > 0L) "blocker" else if (warnings > 0L) "warning" else "pass"
  rownames(issues) <- NULL
  rownames(tables) <- NULL
  out <- list(
    schema_version = "ng_data_preflight.v1",
    generated_at = ng_preflight_timestamp(generated_at),
    status = status,
    counts = list(
      blockers = as.integer(blockers),
      warnings = as.integer(warnings),
      issues = as.integer(nrow(issues))
    ),
    tables = tables,
    issues = issues,
    putative_duplicates = putative_duplicates
  )
  if (!is.null(cleaned_tables)) {
    out$cleaned_tables <- cleaned_tables
    out$cleaned_tables_summary <- do.call(rbind, list(
      ng_preflight_table_summary("geno", cleaned_tables$geno),
      ng_preflight_table_summary("phenotype", cleaned_tables$phenotype),
      ng_preflight_table_summary("candidate_pairs", cleaned_tables$candidate_pairs),
      ng_preflight_table_summary("trait_spec", cleaned_tables$trait_spec),
      ng_preflight_table_summary("marker_map", cleaned_tables$marker_map),
      ng_preflight_table_summary("parent_kinship", cleaned_tables$parent_kinship)
    ))
    rownames(out$cleaned_tables_summary) <- NULL
  }
  if (!is.null(cleaning)) out$cleaning <- cleaning
  out
}

ng_write_data_preflight_json <- function(output_path,
                                         geno = NULL,
                                         phenotype = NULL,
                                         candidate_pairs = NULL,
                                         trait_spec = NULL,
                                         marker_map = NULL,
                                         parent_kinship = NULL,
                                         ploidy = NULL,
                                         putative_duplicate_check = FALSE,
                                         duplicate_threshold = 0.995,
                                         duplicate_maf_min = 0,
                                         duplicate_max_missing_prop = 1,
                                         duplicate_min_compared_markers = 100L,
                                         putative_duplicate_action = c("report", "remove"),
                                         putative_duplicate_return_similarity = FALSE,
                                         generated_at = Sys.time()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to write data preflight JSON")
  }
  payload <- ng_preflight_input_tables(
    geno = geno,
    phenotype = phenotype,
    candidate_pairs = candidate_pairs,
    trait_spec = trait_spec,
    marker_map = marker_map,
    parent_kinship = parent_kinship,
    ploidy = ploidy,
    putative_duplicate_check = putative_duplicate_check,
    duplicate_threshold = duplicate_threshold,
    duplicate_maf_min = duplicate_maf_min,
    duplicate_max_missing_prop = duplicate_max_missing_prop,
    duplicate_min_compared_markers = duplicate_min_compared_markers,
    putative_duplicate_action = putative_duplicate_action,
    putative_duplicate_return_similarity = putative_duplicate_return_similarity,
    generated_at = generated_at
  )
  payload$cleaned_tables <- NULL
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
