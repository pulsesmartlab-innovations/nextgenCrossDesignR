ng_crossing_plan_find_col <- function(data, candidates, required = TRUE, label = "column") {
  names_data <- names(data)
  hit <- match(tolower(candidates), tolower(names_data), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) return(names_data[[hit[[1L]]]])
  if (isTRUE(required)) ng_stop("Could not find ", label, ": ", paste(candidates, collapse = ", "))
  NULL
}

ng_crossing_plan_family_root <- function(cross) {
  cross <- trimws(as.character(cross))
  out <- sub("-[^-]+$", "", cross)
  out[!nzchar(out) | is.na(out)] <- cross[!nzchar(out) | is.na(out)]
  out
}

ng_crossing_plan_direction <- function(direction) {
  x <- trimws(tolower(as.character(direction)))
  increase <- c("max", "maximize", "maximise", "increase", "higher", "high", "positive", "+")
  decrease <- c("min", "minimize", "minimise", "decrease", "lower", "low", "negative", "-")
  out <- ifelse(x %in% decrease, "decrease", ifelse(x %in% increase, "increase", NA_character_))
  if (anyNA(out)) ng_stop("Invalid crossing-plan trait direction: ", paste(unique(x[is.na(out)]), collapse = ", "))
  out
}

ng_crossing_plan_criteria <- function(candidates,
                                      criteria = NULL,
                                      exclude_cols = character()) {
  if (is.null(criteria)) {
    numeric_cols <- names(candidates)[vapply(candidates, is.numeric, logical(1L))]
    numeric_cols <- setdiff(numeric_cols, exclude_cols)
    if (!length(numeric_cols)) ng_stop("criteria is required when no numeric trait columns can be inferred")
    criteria <- data.frame(
      trait = numeric_cols,
      label = numeric_cols,
      weight = rep(1 / length(numeric_cols), length(numeric_cols)),
      direction = rep("increase", length(numeric_cols)),
      meaning = paste("Higher", numeric_cols, "is preferred"),
      stringsAsFactors = FALSE
    )
  } else {
    criteria <- as.data.frame(criteria, stringsAsFactors = FALSE, check.names = FALSE)
    trait_col <- ng_crossing_plan_find_col(criteria, c("trait", "column", "name"), label = "criteria trait column")
    criteria$trait <- as.character(criteria[[trait_col]])
    if (!("label" %in% names(criteria))) criteria$label <- criteria$trait
    if (!("weight" %in% names(criteria))) criteria$weight <- 1
    if (!("direction" %in% names(criteria))) criteria$direction <- "increase"
    if (!("meaning" %in% names(criteria))) criteria$meaning <- ""
    criteria <- criteria[, c("trait", "label", "weight", "direction", "meaning"), drop = FALSE]
  }
  criteria$trait <- trimws(as.character(criteria$trait))
  criteria$label <- trimws(as.character(criteria$label))
  criteria$weight <- suppressWarnings(as.numeric(criteria$weight))
  criteria$direction <- ng_crossing_plan_direction(criteria$direction)
  criteria$meaning <- as.character(criteria$meaning)
  criteria <- criteria[nzchar(criteria$trait) & criteria$trait %in% names(candidates), , drop = FALSE]
  if (!nrow(criteria)) ng_stop("No criteria traits matched candidate columns")
  criteria$weight[!is.finite(criteria$weight) | criteria$weight < 0] <- 0
  if (sum(criteria$weight) <= 0) criteria$weight <- rep(1 / nrow(criteria), nrow(criteria))
  criteria$weight <- criteria$weight / sum(criteria$weight)
  rownames(criteria) <- NULL
  criteria
}

ng_crossing_plan_score_trait <- function(x, direction) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  r <- range(x[ok], na.rm = TRUE)
  if (!all(is.finite(r)) || abs(diff(r)) < .Machine$double.eps) {
    out[ok] <- 1
  } else {
    out[ok] <- (x[ok] - r[[1L]]) / diff(r)
    if (identical(direction, "decrease")) out[ok] <- 1 - out[ok]
  }
  out
}

ng_crossing_plan_prepare <- function(candidates,
                                     criteria = NULL,
                                     line_col = NULL,
                                     cross_col = NULL,
                                     color_col = NULL,
                                     family_col = NULL,
                                     top_n = 5L) {
  candidates <- as.data.frame(candidates, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(line_col)) {
    line_col <- ng_crossing_plan_find_col(candidates, c("line", "entry", "parent", "parent_id", "id"), label = "line ID column")
  }
  if (is.null(cross_col)) {
    cross_col <- ng_crossing_plan_find_col(candidates, c("cross", "family", "pedigree"), label = "cross/family column")
  }
  if (is.null(color_col)) {
    color_col <- ng_crossing_plan_find_col(candidates, c("seed_color", "seed color", "color", "colour"), required = FALSE)
  }
  if (is.null(family_col)) {
    family_col <- ng_crossing_plan_find_col(candidates, c("family_root", "family root", "family_id"), required = FALSE)
  }

  exclude <- c(line_col, cross_col, color_col, family_col)
  criteria <- ng_crossing_plan_criteria(candidates, criteria = criteria, exclude_cols = exclude)

  out <- data.frame(
    Cross = as.character(candidates[[cross_col]]),
    Line = as.character(candidates[[line_col]]),
    Seed.Color = if (!is.null(color_col)) as.character(candidates[[color_col]]) else NA_character_,
    Family.Root = if (!is.null(family_col)) as.character(candidates[[family_col]]) else ng_crossing_plan_family_root(candidates[[cross_col]]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  index <- rep(0, nrow(out))
  for (i in seq_len(nrow(criteria))) {
    trait <- criteria$trait[[i]]
    label <- criteria$label[[i]]
    score_label <- paste(label, "Score")
    raw <- suppressWarnings(as.numeric(candidates[[trait]]))
    score <- ng_crossing_plan_score_trait(raw, criteria$direction[[i]])
    out[[label]] <- raw
    out[[score_label]] <- score
    index <- index + criteria$weight[[i]] * ifelse(is.finite(score), score, 0)
  }
  out$Selection.Index <- index
  out$Family.Count <- stats::ave(out$Line, out$Family.Root, FUN = length)
  out$Within.Family.Rank <- stats::ave(
    -out$Selection.Index,
    out$Family.Root,
    FUN = function(x) rank(x, ties.method = "first")
  )
  out$Within.Family.Rank <- as.integer(out$Within.Family.Rank)
  out$Full.Sib.Excluded <- ifelse(out$Within.Family.Rank > 1L, "Yes", "No")

  eligible <- which(out$Within.Family.Rank == 1L & is.finite(out$Selection.Index))
  selected <- eligible[order(-out$Selection.Index[eligible], out$Line[eligible], method = "radix")]
  selected <- utils::head(selected, as.integer(top_n))
  out$Selected.Top.N <- "No"
  out$Selected.Top.N[selected] <- "Yes"

  list(scored = out, criteria = criteria, selected_rows = selected)
}

ng_crossing_plan_breeding_use <- function(row, rank, criteria) {
  pieces <- character()
  if (rank == 1L) pieces <- c(pieces, "Anchor parent: strongest selection index among unique families")
  harvest <- grep("harvest|lodg", criteria$trait, ignore.case = TRUE)
  if (length(harvest)) {
    label <- criteria$label[[harvest[[1L]]]]
    score_col <- paste(label, "Score")
    if (score_col %in% names(row) && is.finite(row[[score_col]]) && row[[score_col]] >= 0.65) {
      pieces <- c(pieces, "upright/harvestability donor")
    } else {
      pieces <- c(pieces, "monitor lodging or harvestability")
    }
  }
  if (length(grep("seed|tsw|weight|quality", criteria$trait, ignore.case = TRUE))) {
    pieces <- c(pieces, "seed-size/quality contributor")
  }
  if (!length(pieces)) pieces <- "high-index unique-family parent"
  paste(unique(pieces), collapse = "; ")
}

ng_crossing_plan_top_parents <- function(prepared) {
  scored <- prepared$scored
  selected <- scored[prepared$selected_rows, , drop = FALSE]
  selected <- selected[order(selected$Selection.Index, decreasing = TRUE), , drop = FALSE]
  selected$Rank <- seq_len(nrow(selected))
  selected$Full.sib.exclusion.status <- "Selected - unique family"
  selected$Breeding.Use <- vapply(seq_len(nrow(selected)), function(i) {
    ng_crossing_plan_breeding_use(selected[i, , drop = FALSE], selected$Rank[[i]], prepared$criteria)
  }, character(1L))

  trait_labels <- prepared$criteria$label
  score_labels <- paste(trait_labels, "Score")
  cols <- c(
    "Rank", "Line", "Seed.Color", "Cross", "Family.Root",
    trait_labels, score_labels, "Selection.Index",
    "Full.sib.exclusion.status", "Breeding.Use"
  )
  selected[, cols[cols %in% names(selected)], drop = FALSE]
}

ng_crossing_plan_rationale <- function(p1, p2, criteria) {
  pieces <- character()
  if (p1$Rank[[1L]] == 1L || p2$Rank[[1L]] == 1L) pieces <- c(pieces, "uses top index anchor")
  if (!is.na(p1$Seed.Color[[1L]]) && !is.na(p2$Seed.Color[[1L]]) && !identical(p1$Seed.Color[[1L]], p2$Seed.Color[[1L]])) {
    pieces <- c(pieces, paste(p1$Seed.Color[[1L]], "x", p2$Seed.Color[[1L]], "diversity"))
  }
  harvest <- grep("harvest|lodg", criteria$trait, ignore.case = TRUE)
  if (length(harvest)) {
    score_col <- paste(criteria$label[[harvest[[1L]]]], "Score")
    if (score_col %in% names(p1) && is.finite(p1[[score_col]]) && is.finite(p2[[score_col]]) &&
        p1[[score_col]] >= 0.65 && p2[[score_col]] >= 0.65) {
      pieces <- c(pieces, "both parents upright/harvestable")
    } else {
      pieces <- c(pieces, "monitor lodging in progeny")
    }
  }
  if (!length(pieces)) pieces <- "high-index complementary parents"
  paste(unique(pieces), collapse = "; ")
}

ng_crossing_plan_recommended_crosses <- function(top_parents,
                                                 criteria,
                                                 cross_n = 10L) {
  if (nrow(top_parents) < 2L) return(data.frame())
  rows <- list()
  k <- 0L
  score_cols <- paste(criteria$label, "Score")
  score_cols <- score_cols[score_cols %in% names(top_parents)]
  for (i in seq_len(nrow(top_parents) - 1L)) {
    for (j in (i + 1L):nrow(top_parents)) {
      p1 <- top_parents[i, , drop = FALSE]
      p2 <- top_parents[j, , drop = FALSE]
      if (identical(p1$Family.Root[[1L]], p2$Family.Root[[1L]])) next
      mean_index <- mean(c(p1$Selection.Index[[1L]], p2$Selection.Index[[1L]]), na.rm = TRUE)
      complement <- if (length(score_cols)) {
        mean(abs(as.numeric(p1[1, score_cols]) - as.numeric(p2[1, score_cols])), na.rm = TRUE)
      } else {
        0
      }
      color_bonus <- if (!is.na(p1$Seed.Color[[1L]]) && !is.na(p2$Seed.Color[[1L]]) &&
                           !identical(p1$Seed.Color[[1L]], p2$Seed.Color[[1L]])) 1 else 0
      cross_score <- 0.85 * mean_index + 0.10 * min(p1$Selection.Index[[1L]], p2$Selection.Index[[1L]]) +
        0.03 * complement + 0.02 * color_bonus
      k <- k + 1L
      rows[[k]] <- data.frame(
        Parent.1 = p1$Line[[1L]],
        Parent.1.Cross = p1$Cross[[1L]],
        Parent.1.Color = p1$Seed.Color[[1L]],
        Parent.2 = p2$Line[[1L]],
        Parent.2.Cross = p2$Cross[[1L]],
        Parent.2.Color = p2$Seed.Color[[1L]],
        Mean.Parent.Index = mean_index,
        Cross.Priority.Score = cross_score,
        Breeding.Rationale = ng_crossing_plan_rationale(p1, p2, criteria),
        Full.Sib.Check = "Pass - different family roots",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out <- out[order(out$Cross.Priority.Score, decreasing = TRUE), , drop = FALSE]
  out <- utils::head(out, as.integer(cross_n))
  out$Cross.Rank <- seq_len(nrow(out))
  breaks <- stats::quantile(out$Cross.Priority.Score, probs = c(1 / 3, 2 / 3), na.rm = TRUE, names = FALSE, type = 1)
  out$Tier <- ifelse(out$Cross.Priority.Score >= breaks[[2L]], "Tier 1",
                     ifelse(out$Cross.Priority.Score >= breaks[[1L]], "Tier 2", "Tier 3"))
  out[, c(
    "Cross.Rank", "Parent.1", "Parent.1.Cross", "Parent.1.Color",
    "Parent.2", "Parent.2.Cross", "Parent.2.Color",
    "Mean.Parent.Index", "Cross.Priority.Score", "Tier",
    "Breeding.Rationale", "Full.Sib.Check"
  ), drop = FALSE]
}

ng_crossing_plan_criteria_sheet <- function(criteria) {
  trait_rows <- data.frame(
    Rule...Parameter = paste(criteria$label, "weight"),
    Value = format(criteria$weight, digits = 4, trim = TRUE),
    Meaning = criteria$meaning,
    Editable. = "Yes",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rules <- data.frame(
    Rule...Parameter = c("Full-sib avoidance", "Harvestability interpretation"),
    Value = c("One line per family root", "1 = upright; 5 = pods on ground"),
    Meaning = c(
      "Family root is the cross name before the final dash; lower-ranked same-root lines are excluded from top parent selection.",
      "Harvestability/lodging traits are treated as decrease traits when included in criteria."
    ),
    Editable. = c("No", "No"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rbind(trait_rows, rules)
}

ng_crossing_plan_write_sheet <- function(wb, sheet, title, subtitle, data, start_row = 3L) {
  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
  title_style <- openxlsx::createStyle(fontSize = 15, textDecoration = "bold", fontColour = "#183A37")
  subtitle_style <- openxlsx::createStyle(fontSize = 10, fontColour = "#52605E", wrapText = TRUE)
  openxlsx::writeData(wb, sheet, title, startRow = 1, startCol = 1)
  openxlsx::addStyle(wb, sheet, title_style, rows = 1, cols = 1, stack = TRUE)
  if (!is.null(subtitle) && nzchar(subtitle)) {
    openxlsx::writeData(wb, sheet, subtitle, startRow = 2, startCol = 1)
    openxlsx::addStyle(wb, sheet, subtitle_style, rows = 2, cols = 1, stack = TRUE)
  }
  if (nrow(data)) {
    openxlsx::writeDataTable(wb, sheet, data, startRow = start_row, startCol = 1, tableStyle = "TableStyleMedium2")
    openxlsx::freezePane(wb, sheet, firstActiveRow = start_row + 1L)
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = "auto")
  }
}

ng_crossing_plan_insert_plot <- function(wb, sheet, top_parents) {
  if (!nrow(top_parents)) return(invisible(FALSE))
  img <- tempfile("ng_crossing_plan_", fileext = ".png")
  grDevices::png(img, width = 960, height = 420, res = 130)
  old <- graphics::par(no.readonly = TRUE)
  tryCatch({
    graphics::par(mar = c(7, 4, 2, 1))
    graphics::barplot(
      top_parents$Selection.Index,
      names.arg = top_parents$Line,
      las = 2,
      col = "#2F7D6D",
      border = NA,
      ylab = "Selection index",
      main = "Top parent selection index"
    )
  }, finally = {
    graphics::par(old)
    grDevices::dev.off()
  })
  openxlsx::insertImage(wb, sheet, img, startRow = max(10L, nrow(top_parents) + 6L), startCol = 1, width = 7.6, height = 3.3)
  invisible(TRUE)
}

ng_write_crossing_plan_workbook <- function(output_path,
                                            candidates,
                                            criteria = NULL,
                                            top_n = 5L,
                                            cross_n = 10L,
                                            title = NULL,
                                            line_col = NULL,
                                            cross_col = NULL,
                                            color_col = NULL,
                                            family_col = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    ng_stop("openxlsx is required to write crossing plan workbooks")
  }
  if (missing(output_path) || is.null(output_path) || !nzchar(as.character(output_path[[1L]]))) {
    ng_stop("output_path is required")
  }
  if (is.null(title)) title <- paste0("Top ", as.integer(top_n), " Crossing Candidates with Crossing Plan")
  prepared <- ng_crossing_plan_prepare(
    candidates = candidates,
    criteria = criteria,
    line_col = line_col,
    cross_col = cross_col,
    color_col = color_col,
    family_col = family_col,
    top_n = top_n
  )
  top_parents <- ng_crossing_plan_top_parents(prepared)
  crosses <- ng_crossing_plan_recommended_crosses(top_parents, prepared$criteria, cross_n = cross_n)
  scored <- prepared$scored
  names(scored)[names(scored) == "Selected.Top.N"] <- paste0("Selected.Top.", as.integer(top_n))
  family <- scored[order(scored$Family.Root, scored$Within.Family.Rank, scored$Line), , drop = FALSE]
  family <- family[, c("Family.Root", "Cross", "Line", "Seed.Color", prepared$criteria$label,
                       "Selection.Index", "Within.Family.Rank", "Full.Sib.Excluded"), drop = FALSE]
  criteria_sheet <- ng_crossing_plan_criteria_sheet(prepared$criteria)

  wb <- openxlsx::createWorkbook(creator = "nextgenCrossDesign")
  ng_crossing_plan_write_sheet(
    wb, "Top Parents", title,
    "Selection enforces one representative per family root to avoid selecting likely full-sibs.",
    top_parents,
    start_row = 3L
  )
  ng_crossing_plan_write_sheet(
    wb, "Recommended Crosses", "Recommended Cross Combinations",
    "Crosses are ranked by parent index, complementarity, and full-sib family-root checks.",
    crosses,
    start_row = 3L
  )
  ng_crossing_plan_write_sheet(
    wb, "Selection Index - All Lines", "Selection Index for All Candidate Lines",
    "Rows marked Full-sib excluded are lower-ranked representatives within the same family root.",
    scored,
    start_row = 3L
  )
  ng_crossing_plan_write_sheet(
    wb, "Criteria", "Selection Criteria and Interpretation", NULL, criteria_sheet, start_row = 2L
  )
  ng_crossing_plan_write_sheet(
    wb, "Family Screening", "Full-Sib / Family Screening", NULL, family, start_row = 2L
  )
  ng_crossing_plan_insert_plot(wb, "Top Parents", top_parents)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
