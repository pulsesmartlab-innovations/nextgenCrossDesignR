ng_cpw_clean_sheet_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", as.character(x))
  x <- gsub("_+", "_", x)
  substr(x, 1L, 31L)
}

ng_cpw_first_col <- function(data, candidates) {
  hit <- match(tolower(candidates), tolower(names(data)), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) names(data)[hit[[1L]]] else NULL
}

ng_cpw_tier_sheet_name <- function(tier) {
  key <- tolower(trimws(as.character(tier)))
  if (identical(key, "highly_priority")) return("Highly_Priority")
  if (identical(key, "priority")) return("Priority")
  if (identical(key, "medium_priority")) return("Medium_Priority")
  if (identical(key, "low_priority")) return("Low_Priority")
  ng_cpw_clean_sheet_name(tier)
}

ng_cpw_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

ng_cpw_format_value <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- ifelse(is.finite(x), format(round(x, 4L), trim = TRUE, scientific = FALSE), "")
  out
}

ng_cpw_trait_table <- function(trait_directions, data) {
  if (is.null(trait_directions)) return(data.frame())
  trait_directions <- as.data.frame(trait_directions, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(trait_directions)) return(data.frame())
  trait_col <- ng_cpw_first_col(trait_directions, c("trait", "name"))
  column_col <- ng_cpw_first_col(trait_directions, c("column", "score_col", "trait_column"))
  direction_col <- ng_cpw_first_col(trait_directions, c("direction", "selection_directionn", "selection_direction"))
  weight_col <- ng_cpw_first_col(trait_directions, c("weight", "economic_weight"))
  if (is.null(trait_col)) trait_col <- names(trait_directions)[[1L]]
  trait <- as.character(trait_directions[[trait_col]])
  column <- if (!is.null(column_col)) as.character(trait_directions[[column_col]]) else trait
  column <- ifelse(column %in% names(data), column, paste0("pred_", trait))
  direction <- if (!is.null(direction_col)) as.character(trait_directions[[direction_col]]) else "maximize"
  direction <- ifelse(tolower(direction) %in% c("min", "minimize", "minimise", "low", "lower", "decrease", "-"),
                      "minimize", "maximize")
  weight <- if (!is.null(weight_col)) ng_cpw_numeric(trait_directions[[weight_col]]) else rep(1, length(trait))
  weight[!is.finite(weight)] <- 1
  out <- data.frame(
    trait = trait,
    column = column,
    direction = direction,
    weight = weight,
    available = column %in% names(data),
    stringsAsFactors = FALSE
  )
  out[nzchar(out$trait), , drop = FALSE]
}

ng_cpw_trait_evidence <- function(crosses, trait_info, n_traits = 3L) {
  empty <- rep("", nrow(crosses))
  if (!nrow(crosses) || is.null(trait_info) || !nrow(trait_info)) {
    return(list(favorable = empty, risk = empty))
  }
  trait_info <- trait_info[trait_info$available, , drop = FALSE]
  if (!nrow(trait_info)) return(list(favorable = empty, risk = empty))
  components <- matrix(0, nrow = nrow(crosses), ncol = nrow(trait_info))
  raw_values <- matrix(NA_real_, nrow = nrow(crosses), ncol = nrow(trait_info))
  for (i in seq_len(nrow(trait_info))) {
    x <- ng_cpw_numeric(crosses[[trait_info$column[[i]]]])
    raw_values[, i] <- x
    components[, i] <- ng_rank_normalize(x, bigger_is_better = identical(trait_info$direction[[i]], "maximize"))
  }
  favorable <- character(nrow(crosses))
  risk <- character(nrow(crosses))
  for (r in seq_len(nrow(crosses))) {
    ok <- is.finite(raw_values[r, ]) & is.finite(components[r, ])
    if (!any(ok)) {
      favorable[[r]] <- ""
      risk[[r]] <- ""
    } else {
      ord_fav <- order(components[r, ], decreasing = TRUE, na.last = NA)
      ord_risk <- order(components[r, ], decreasing = FALSE, na.last = NA)
      ord_fav <- ord_fav[seq_len(min(length(ord_fav), as.integer(n_traits)))]
      ord_risk <- ord_risk[seq_len(min(length(ord_risk), as.integer(n_traits)))]
      favorable[[r]] <- paste(
        paste0(trait_info$trait[ord_fav], "=", ng_cpw_format_value(raw_values[r, ord_fav])),
        collapse = "; "
      )
      risk[[r]] <- paste(
        paste0(trait_info$trait[ord_risk], "=", ng_cpw_format_value(raw_values[r, ord_risk])),
        collapse = "; "
      )
    }
  }
  list(favorable = favorable, risk = risk)
}

ng_cpw_duplicate_parent_set <- function(duplicate_pairs) {
  if (is.null(duplicate_pairs)) return(character())
  duplicate_pairs <- as.data.frame(duplicate_pairs, stringsAsFactors = FALSE)
  p1 <- ng_cpw_first_col(duplicate_pairs, c("parent1", "sample1", "id1", "line1"))
  p2 <- ng_cpw_first_col(duplicate_pairs, c("parent2", "sample2", "id2", "line2"))
  if (is.null(p1) || is.null(p2) || !nrow(duplicate_pairs)) return(character())
  unique(c(as.character(duplicate_pairs[[p1]]), as.character(duplicate_pairs[[p2]])))
}

ng_cpw_parent_use_flag <- function(crosses, parent_use) {
  out <- rep("OK", nrow(crosses))
  if (is.null(parent_use) || !nrow(crosses)) return(out)
  parent_use <- as.data.frame(parent_use, stringsAsFactors = FALSE)
  parent_col <- ng_cpw_first_col(parent_use, c("parent", "line", "id", "NAME"))
  count_col <- ng_cpw_first_col(parent_use, c("crosses_selected", "Actual_Participation", "count", "uses"))
  if (is.null(parent_col) || is.null(count_col)) return(out)
  counts <- ng_cpw_numeric(parent_use[[count_col]])
  names(counts) <- as.character(parent_use[[parent_col]])
  threshold <- stats::quantile(counts[is.finite(counts)], probs = 0.90, na.rm = TRUE, names = FALSE, type = 1)
  if (!is.finite(threshold)) return(out)
  p1 <- as.character(crosses$parent1)
  p2 <- as.character(crosses$parent2)
  max_use <- pmax(counts[p1], counts[p2], na.rm = TRUE)
  out[is.finite(max_use) & max_use >= threshold & max_use > 1] <- "High parent use"
  out
}

# Per-trait mid-parent GEBV columns ride on the cross table as `<trait>_mean_gebv` (R/39). They
# are excluded from the workbook by default to keep it lean; `include_trait_gebv = TRUE` (R/39
# runner formal, threaded through ng_write_cross_priority_workbook) opts them back in, relabelled
# so breeders see them as mid-parent GEBV rather than the internal column name.
ng_cpw_gebv_cols <- function(data) {
  grep("_mean_gebv$", names(data), value = TRUE)
}

ng_cpw_gebv_label <- function(col) {
  sub("_mean_gebv$", "_mid_parent_gebv", col)
}

# Portfolio + risk columns for the workbook, in breeder-reading order: where the cross sits, how
# much to trust it, and (multi-trait only) which trait is driving that doubt. Returns an empty
# list when the run carries no portfolio annotation at all.
#
# risk_bin / cross_confidence are WITHIN-RUN quantities (tertiles and a min-max normalization of
# the full post-filter candidate pool), so the header carries that caveat -- a spreadsheet
# outlives the app session that explains it, and these labels are not comparable across runs.
ng_cpw_portfolio_cols <- function(crosses, n) {
  has <- function(col) col %in% names(crosses) && !all(is.na(crosses[[col]]))
  out <- list()
  if (has("portfolio_profile")) out[["portfolio_profile"]] <- as.character(crosses$portfolio_profile)
  if (has("cross_level"))       out[["cross_level"]]       <- ng_cpw_numeric(crosses$cross_level)
  if (has("cross_upside"))      out[["cross_upside"]]      <- ng_cpw_numeric(crosses$cross_upside)
  if (has("risk_bin"))          out[["risk_bin_within_candidate_pool"]] <- as.character(crosses$risk_bin)
  if (has("cross_confidence"))  out[["cross_confidence_within_candidate_pool"]] <- ng_cpw_numeric(crosses$cross_confidence)
  # Multi-trait only: names the component trait contributing most of the index prediction error.
  if (has("risk_driver_trait")) out[["risk_driver_trait"]] <- as.character(crosses$risk_driver_trait)
  if (has("risk_driver_share")) out[["risk_driver_share"]] <- ng_cpw_numeric(crosses$risk_driver_share)
  # portfolio_basis is constant per run, but a spreadsheet gets separated from its run notes --
  # carrying it per row is what stops a rank-index quadrant being read as a decomposition.
  if (has("portfolio_basis"))   out[["portfolio_basis"]]   <- as.character(crosses$portfolio_basis)
  lapply(out, function(v) if (length(v) == n) v else rep(NA, n))
}

# The `Checks` sheet: the reference itself, one row per trait. Breeders read the check first
# and the crosses second, so it gets its own sheet rather than being inferred from columns.
ng_cpw_checks_sheet <- function(trait_check_reference, crosses) {
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  if (!nrow(spec)) return(data.frame())
  # n_wrong_side (R/51_check_reference.R::ng_attach_check_reference) is counted once, over the
  # full candidate pool, before allocation/selection narrows it down to a crossing plan. Its
  # denominator must therefore be that SAME pool (diagnostics$n_candidates) -- not nrow(crosses),
  # which at the workbook layer is typically the much smaller SELECTED plan and would pair a
  # candidate-pool numerator with a selected-plan denominator. Falls back to nrow(crosses) only
  # if diagnostics does not carry the count.
  n_total <- trait_check_reference$diagnostics$n_candidates
  if (is.null(n_total)) n_total <- nrow(crosses)
  wrong <- trait_check_reference$diagnostics$n_wrong_side
  not_eval <- trait_check_reference$diagnostics$n_not_evaluable
  value <- vapply(seq_len(nrow(spec)), function(k) {
    suppressWarnings(as.numeric(trait_check_reference$values[[spec$trait[[k]]]][[spec$check[[k]]]]))
  }, numeric(1))
  data.frame(
    trait = spec$trait,
    check_id = spec$check,
    direction = ifelse(spec$reject_if == "below", "want above", "want below"),
    value = value,
    source = as.character(trait_check_reference$source[spec$trait]),
    # A check whose own `value` never resolved (NA) was never comparable to ANY cross for this
    # trait: rendering "0 / N" here would be an affirmative false claim ("0 crosses on the wrong
    # side") when in truth zero crosses were ever evaluated. "not evaluable" is the only honest
    # string for that case -- "0 / 66" must be unreachable for a check that was never compared.
    # When the check itself DID resolve but some individual crosses' own mean/variance did not
    # (diagnostics$n_not_evaluable > 0 for this trait), that partial gap is surfaced alongside the
    # wrong/total ratio rather than silently folded into either count.
    n_crosses_on_wrong_side = vapply(seq_len(nrow(spec)), function(k) {
      tr <- spec$trait[[k]]
      if (!is.finite(value[[k]])) return("not evaluable")
      w <- wrong[[tr]]
      ne <- not_eval[[tr]]
      base <- sprintf("%d / %d", if (is.null(w)) NA_integer_ else as.integer(w), n_total)
      if (!is.null(ne) && length(ne) == 1L && is.finite(ne) && ne > 0) {
        paste0(base, sprintf(" (%d not evaluable)", as.integer(ne)))
      } else {
        base
      }
    }, character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
}

# Insert named columns immediately before an existing column, preserving order.
ng_cpw_insert_before <- function(df, before, cols) {
  if (!length(cols)) return(df)
  for (nm in names(cols)) df[[nm]] <- cols[[nm]]
  at <- match(before, names(df))
  if (is.na(at)) return(df)
  added <- names(cols)
  keep <- setdiff(names(df), added)
  at <- match(before, keep)
  df[, append(keep, added, after = at - 1L), drop = FALSE]
}

ng_cpw_make_selected <- function(crosses, trait_info, parent_use, duplicate_pairs, block_size,
                                 include_trait_gebv = FALSE) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("priority_rank" %in% names(crosses)) || !("priority_tier" %in% names(crosses))) {
    crosses <- ng_rank_cross_priority(crosses)
  }
  crosses <- crosses[order(crosses$priority_rank), , drop = FALSE]
  block_size <- as.integer(block_size[[1L]])
  if (!is.finite(block_size) || block_size < 1L) block_size <- 10L
  evidence <- ng_cpw_trait_evidence(crosses, trait_info)
  duplicate_parents <- ng_cpw_duplicate_parent_set(duplicate_pairs)
  dup_flag <- ifelse(crosses$parent1 %in% duplicate_parents | crosses$parent2 %in% duplicate_parents,
                     "Parent in putative duplicate pair", "OK")
  threshold_violation <- if ("multi_trait_threshold_violation" %in% names(crosses)) {
    ng_cpw_numeric(crosses$multi_trait_threshold_violation)
  } else {
    rep(NA_real_, nrow(crosses))
  }
  out <- data.frame(
    Cross_Rank = seq_len(nrow(crosses)),
    Cross_ID = sprintf("C%03d", seq_len(nrow(crosses))),
    Cross_Block = sprintf("Block_%02d", ceiling(seq_len(nrow(crosses)) / block_size)),
    Block_Position = ((seq_len(nrow(crosses)) - 1L) %% block_size) + 1L,
    Parent_1 = as.character(crosses$parent1),
    Parent_2 = as.character(crosses$parent2),
    priority_rank = as.integer(crosses$priority_rank),
    priority_tier = as.character(crosses$priority_tier),
    priority_index = if ("priority_index" %in% names(crosses)) ng_cpw_numeric(crosses$priority_index) else NA_real_,
    multi_trait_score = if ("multi_trait_score" %in% names(crosses)) ng_cpw_numeric(crosses$multi_trait_score) else NA_real_,
    pair_kinship = if ("pair_kinship" %in% names(crosses)) ng_cpw_numeric(crosses$pair_kinship) else NA_real_,
    threshold_violation = threshold_violation,
    top_favorable_traits = evidence$favorable,
    top_risk_traits = evidence$risk,
    parent_use_flag = ng_cpw_parent_use_flag(crosses, parent_use),
    duplicate_qc_flag = dup_flag,
    package_evidence_note = ifelse(
      is.finite(threshold_violation) & threshold_violation > 0,
      "Ranked with soft threshold penalty; review risk traits and breeder notes",
      "Ranked from supplied trait directions, weights, kinship, and parent-use constraints"
    ),
    Breeder_Rationale = "",
    Breeder_Notes = "",
    Final_Decision = "",
    Crossing_Status = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  # Portfolio + risk (R/37_cross_portfolio_risk.R). The workbook is what leaves the building, so
  # a plan exported here should carry the same risk/portfolio read the app shows -- otherwise a
  # breeder working from the spreadsheet sees priority tiers with no indication of which crosses
  # are speculative. Inserted before the free-text breeder columns so the decision fields stay
  # rightmost. Every column is optional: runs predating the layer (or where the annotation could
  # not be resolved) simply omit it rather than emitting an all-NA column.
  out <- ng_cpw_insert_before(out, "Breeder_Rationale",
                              ng_cpw_portfolio_cols(crosses, nrow(out)))
  trait_cols <- trait_info$column[trait_info$available]
  trait_cols <- trait_cols[trait_cols %in% names(crosses)]
  for (col in trait_cols) out[[col]] <- crosses[[col]]
  if (isTRUE(include_trait_gebv)) {
    gebv_cols <- ng_cpw_gebv_cols(crosses)
    for (col in gebv_cols) out[[ng_cpw_gebv_label(col)]] <- ng_cpw_numeric(crosses[[col]])
  }
  # Per-trait check-reference columns (R/51_check_reference.R / R/39 runner): so a breeder
  # scanning a cross sees its margin against the check without cross-referencing the Checks
  # sheet. No `keep` vector exists in this function (unlike ng_cpw_candidate_table), so the
  # check columns are appended the same way trait/GEBV columns above are -- copied straight
  # from `crosses` when present, absent entirely on a no-checks run.
  check_cols <- c(grep("_check_id$|_check_value$|_vs_check$|_check_ok$|_p_beat_check$",
                       names(crosses), value = TRUE),
                  intersect(c("checks_all_ok", "p_beat_all_checks"), names(crosses)))
  for (col in check_cols) out[[col]] <- crosses[[col]]
  out
}

ng_cpw_candidate_table <- function(scored, selected, include_trait_gebv = FALSE) {
  if (is.null(scored)) return(data.frame())
  scored <- as.data.frame(scored, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(scored)) return(data.frame())
  score_col <- if ("multi_trait_score" %in% names(scored)) "multi_trait_score" else ng_cpw_first_col(scored, c("priority_index", "usefulness_pmv_gebv"))
  if (!is.null(score_col)) {
    ord <- order(-ng_cpw_numeric(scored[[score_col]]), scored$parent1, scored$parent2, na.last = TRUE)
    scored <- scored[ord, , drop = FALSE]
  }
  key <- paste(selected$Parent_1, selected$Parent_2, sep = "\r")
  cand_key <- paste(scored$parent1, scored$parent2, sep = "\r")
  out <- data.frame(
    Candidate_Rank = seq_len(nrow(scored)),
    Selected_in_plan = ifelse(cand_key %in% key, "Yes", "No"),
    Selection_Status = ifelse(cand_key %in% key, selected$priority_tier[match(cand_key, key)], "Not selected"),
    Parent_1 = as.character(scored$parent1),
    Parent_2 = as.character(scored$parent2),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  keep <- c("multi_trait_score", "pair_kinship", "multi_trait_threshold_violation", "priority_rank", "priority_tier")
  pred_cols <- grep("^pred_", names(scored), value = TRUE)
  keep <- unique(c(keep, pred_cols))
  # Per-trait check-reference columns ride along on the candidate table too, same as Selected_All.
  keep <- c(keep, grep("_check_id$|_check_value$|_vs_check$|_check_ok$|_p_beat_check$",
                       names(scored), value = TRUE),
            "checks_all_ok", "p_beat_all_checks")
  keep <- intersect(unique(keep), names(scored))
  out <- cbind(out, scored[, keep, drop = FALSE])
  if (isTRUE(include_trait_gebv)) {
    gebv_cols <- ng_cpw_gebv_cols(scored)
    for (col in gebv_cols) out[[ng_cpw_gebv_label(col)]] <- ng_cpw_numeric(scored[[col]])
  }
  out
}

ng_cpw_dashboard <- function(selected, scored, trait_info, parent_use, duplicate_pairs, n_crosses_requested) {
  tier_counts <- table(selected$priority_tier)
  tier_text <- paste(paste(names(tier_counts), as.integer(tier_counts), sep = "="), collapse = "; ")
  data.frame(
    Item = c(
      "Requested crosses",
      "Selected crosses",
      "Candidate crosses",
      "Priority tiers",
      "Traits included",
      "Parents used",
      "Putative duplicate pairs",
      "Breeder rationale"
    ),
    Value = c(
      if (is.null(n_crosses_requested)) nrow(selected) else as.integer(n_crosses_requested[[1L]]),
      nrow(selected),
      if (is.null(scored)) NA_integer_ else nrow(as.data.frame(scored)),
      tier_text,
      nrow(trait_info[trait_info$available, , drop = FALSE]),
      length(unique(c(selected$Parent_1, selected$Parent_2))),
      if (is.null(duplicate_pairs)) 0L else nrow(as.data.frame(duplicate_pairs)),
      "Editable by user"
    ),
    Notes = c(
      "User-requested crossing capacity.",
      "Rows included in the final selected crossing plan.",
      "Candidate parent pairs scored before final selection.",
      "Priority tiers are assigned after optimization for practical execution.",
      "Only available trait columns are used in evidence summaries.",
      "Unique parents represented in selected crosses.",
      "Pairs identified during putative duplicate genotype QC.",
      "Breeder_Rationale, Breeder_Notes, Final_Decision, and Crossing_Status are intentionally blank."
    ),
    stringsAsFactors = FALSE
  )
}

ng_cpw_scoring_method <- function() {
  data.frame(
    Section = c(
      "Selection engine",
      "Priority tiering",
      "Evidence columns",
      "Editable breeder columns",
      "Breeder rationale policy",
      "Duplicate and parent-use QC"
    ),
    Details = c(
      "Crosses are ranked from supplied trait directions, trait weights, kinship, threshold penalties, and parent-use constraints.",
      "Priority tiers divide the selected plan into practical execution groups; users can change tier breaks and crossing capacity.",
      "top_favorable_traits and top_risk_traits summarize objective trait evidence from the selected cross table.",
      "Breeder_Rationale, Breeder_Notes, Final_Decision, and Crossing_Status are editable fields for program-specific decisions.",
      "The package does not write one-size-fits-all breeding rationale text because rationale differs by breeder, market, nursery, and cycle.",
      "parent_use_flag and duplicate_qc_flag surface review points without replacing breeder judgment."
    ),
    stringsAsFactors = FALSE
  )
}

ng_cross_priority_workbook_tables <- function(crosses,
                                              scored = NULL,
                                              trait_directions = NULL,
                                              parent_use = NULL,
                                              duplicate_pairs = NULL,
                                              figures = NULL,
                                              n_crosses_requested = NULL,
                                              block_size = 10L,
                                              include_trait_gebv = FALSE,
                                              trait_check_reference = NULL) {
  crosses <- as.data.frame(crosses, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(crosses)) ng_stop("crosses must contain at least one selected cross")
  if (!all(c("parent1", "parent2") %in% names(crosses))) {
    ng_stop("crosses must contain parent1 and parent2 columns")
  }
  trait_info <- ng_cpw_trait_table(trait_directions, crosses)
  selected <- ng_cpw_make_selected(crosses, trait_info, parent_use, duplicate_pairs, block_size,
                                   include_trait_gebv = include_trait_gebv)
  out <- list(
    Dashboard = ng_cpw_dashboard(selected, scored, trait_info, parent_use, duplicate_pairs, n_crosses_requested),
    Scoring_Method = ng_cpw_scoring_method(),
    Trait_Directions = trait_info,
    Selected_All = selected
  )
  if (!is.null(trait_check_reference)) {
    out$Checks <- ng_cpw_checks_sheet(trait_check_reference, selected)
  }
  tier_order <- c("highly_priority", "priority", "medium_priority", "low_priority")
  for (tier in tier_order) {
    sheet <- ng_cpw_tier_sheet_name(tier)
    idx <- tolower(selected$priority_tier) == tier
    out[[sheet]] <- selected[idx, , drop = FALSE]
  }
  out$Candidate_Crosses <- ng_cpw_candidate_table(scored, selected, include_trait_gebv = include_trait_gebv)
  out$Parent_Use_QC <- if (is.null(parent_use)) data.frame() else as.data.frame(parent_use, stringsAsFactors = FALSE)
  out$Duplicate_QC <- if (is.null(duplicate_pairs)) data.frame() else as.data.frame(duplicate_pairs, stringsAsFactors = FALSE)
  if (!is.null(figures)) out$Figure_Index <- as.data.frame(figures, stringsAsFactors = FALSE)
  out
}

ng_cpw_write_sheet <- function(wb, sheet, data) {
  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
  title_style <- openxlsx::createStyle(fontSize = 13, textDecoration = "bold", fontColour = "#1F4E78")
  header_style <- openxlsx::createStyle(textDecoration = "bold", fontColour = "#FFFFFF", fgFill = "#1F4E78")
  openxlsx::writeData(wb, sheet, gsub("_", " ", sheet), startRow = 1L, startCol = 1L)
  openxlsx::addStyle(wb, sheet, title_style, rows = 1L, cols = 1L, stack = TRUE)
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(data)) {
    openxlsx::writeData(wb, sheet, "No rows available", startRow = 3L, startCol = 1L)
    return(invisible(FALSE))
  }
  openxlsx::writeDataTable(
    wb, sheet, data, startRow = 3L, startCol = 1L,
    tableStyle = "TableStyleMedium2", withFilter = TRUE
  )
  openxlsx::addStyle(wb, sheet, header_style, rows = 3L, cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 4L)
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
  editable <- which(names(data) %in% c("Breeder_Rationale", "Breeder_Notes", "Final_Decision", "Crossing_Status"))
  if (length(editable)) {
    editable_style <- openxlsx::createStyle(fgFill = "#FFF2CC")
    openxlsx::addStyle(
      wb, sheet, editable_style,
      rows = seq.int(4L, 3L + nrow(data)), cols = editable,
      gridExpand = TRUE, stack = TRUE
    )
  }
  invisible(TRUE)
}

ng_write_cross_priority_workbook <- function(output_path,
                                             crosses,
                                             scored = NULL,
                                             trait_directions = NULL,
                                             parent_use = NULL,
                                             duplicate_pairs = NULL,
                                             figures = NULL,
                                             n_crosses_requested = NULL,
                                             block_size = 10L,
                                             include_trait_gebv = FALSE,
                                             trait_check_reference = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    ng_stop("openxlsx is required to write cross priority workbooks")
  }
  if (missing(output_path) || is.null(output_path) || !nzchar(as.character(output_path[[1L]]))) {
    ng_stop("output_path is required")
  }
  tables <- ng_cross_priority_workbook_tables(
    crosses = crosses,
    scored = scored,
    trait_directions = trait_directions,
    parent_use = parent_use,
    duplicate_pairs = duplicate_pairs,
    figures = figures,
    n_crosses_requested = n_crosses_requested,
    block_size = block_size,
    include_trait_gebv = include_trait_gebv,
    trait_check_reference = trait_check_reference
  )
  wb <- openxlsx::createWorkbook(creator = "nextgenCrossDesign")
  for (sheet in names(tables)) {
    ng_cpw_write_sheet(wb, ng_cpw_clean_sheet_name(sheet), tables[[sheet]])
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
