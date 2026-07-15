ng_hhvr_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&#39;", x, fixed = TRUE)
  x
}

ng_hhvr_humanize <- function(x) {
  x <- gsub("_", " ", as.character(x), fixed = TRUE)
  x[is.na(x)] <- ""
  x
}

ng_hhvr_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "")
}

ng_hhvr_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), formatC(x, format = "f", digits = 0), "")
}

ng_hhvr_first <- function(x, default = "") {
  x <- as.character(x)
  x <- x[nzchar(x) & !is.na(x)]
  if (length(x)) x[[1]] else default
}

ng_hhvr_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

ng_hhvr_metric_label <- function(metric) {
  labels <- c(
    mean_realized_index = "Selection index",
    mean_realized_yield = "Yield",
    mean_realized_disease = "Disease",
    mean_realized_quality = "Quality",
    unique_parents = "Unique parents",
    max_parent_use = "Max parent use",
    mean_pair_kinship = "Pair kinship",
    group_coancestry = "Group coancestry"
  )
  out <- labels[as.character(metric)]
  out[is.na(out)] <- ng_hhvr_humanize(metric[is.na(out)])
  unname(out)
}

ng_hhvr_metric_direction <- function(metric) {
  directions <- c(
    mean_realized_index = "maximize",
    mean_realized_yield = "maximize",
    mean_realized_disease = "minimize",
    mean_realized_quality = "maximize",
    unique_parents = "maximize",
    max_parent_use = "minimize",
    mean_pair_kinship = "minimize",
    group_coancestry = "minimize"
  )
  out <- directions[[as.character(metric)]]
  if (is.null(out)) "maximize" else out
}

ng_hhvr_scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  finite <- is.finite(x)
  out <- rep(0.5, length(x))
  if (!any(finite)) return(out)
  r <- range(x[finite], na.rm = TRUE)
  if (!all(is.finite(r)) || abs(diff(r)) < .Machine$double.eps) return(out)
  out[finite] <- (x[finite] - r[[1]]) / diff(r)
  out
}

ng_hhvr_read_csv <- function(results_dir, prefix, suffix, required = FALSE) {
  if (is.null(results_dir) || is.null(prefix)) {
    if (isTRUE(required)) ng_stop("results_dir and prefix are required to read ", suffix)
    return(data.frame(stringsAsFactors = FALSE))
  }
  path <- file.path(results_dir, paste0(prefix, "_", suffix, ".csv"))
  if (!file.exists(path)) {
    if (isTRUE(required)) ng_stop("visual report input not found: ", path)
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

ng_hhvr_method_summary <- function(summary) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  required <- c("method")
  missing <- setdiff(required, names(summary))
  if (length(missing)) ng_stop("visual report summary missing columns: ", paste(missing, collapse = ", "))

  metric_cols <- intersect(c(
    "mean_realized_index",
    "mean_realized_yield",
    "mean_realized_disease",
    "mean_realized_quality",
    "unique_parents",
    "max_parent_use",
    "mean_pair_kinship",
    "group_coancestry"
  ), names(summary))
  methods <- unique(as.character(summary$method))
  rows <- vector("list", length(methods))
  for (i in seq_along(methods)) {
    part <- summary[summary$method == methods[[i]], , drop = FALSE]
    row <- data.frame(
      method = methods[[i]],
      benchmark_role = if ("benchmark_role" %in% names(part)) ng_hhvr_first(part$benchmark_role) else "",
      external_tool = if ("external_tool" %in% names(part)) ng_hhvr_first(part$external_tool) else "",
      implementation = if ("implementation" %in% names(part)) ng_hhvr_first(part$implementation) else "",
      exact_external_status = if ("exact_external_status" %in% names(part)) ng_hhvr_first(part$exact_external_status) else "",
      records = nrow(part),
      stringsAsFactors = FALSE
    )
    for (metric in metric_cols) row[[metric]] <- ng_hhvr_mean(part[[metric]])
    rows[[i]] <- row
  }
  out <- ng_bind_rows_fill(rows)
  if ("mean_realized_index" %in% names(out)) {
    out <- out[order(-suppressWarnings(as.numeric(out$mean_realized_index)), out$method), , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

ng_hhvr_recommended_method <- function(method_summary) {
  if (!nrow(method_summary)) return("")
  pool <- method_summary
  if ("benchmark_role" %in% names(pool) && any(pool$benchmark_role == "candidate")) {
    pool <- pool[pool$benchmark_role == "candidate", , drop = FALSE]
  }
  if ("mean_realized_index" %in% names(pool)) {
    values <- suppressWarnings(as.numeric(pool$mean_realized_index))
    finite <- is.finite(values)
    if (any(finite)) return(as.character(pool$method[which.max(ifelse(finite, values, -Inf))]))
  }
  as.character(pool$method[[1]])
}

ng_hhvr_kpi_card <- function(label, value, note = "") {
  paste0(
    "<div class=\"kpi\"><div class=\"kpi-label\">", ng_hhvr_escape(label), "</div>",
    "<div class=\"kpi-value\">", ng_hhvr_escape(value), "</div>",
    "<div class=\"kpi-note\">", ng_hhvr_escape(note), "</div></div>"
  )
}

ng_hhvr_role_class <- function(role) {
  role <- tolower(trimws(as.character(role)))
  ifelse(role == "candidate", "candidate", ifelse(role == "baseline", "baseline", "neutral"))
}

ng_hhvr_role_badge <- function(role, implementation = "") {
  role <- ng_hhvr_first(role, "method")
  impl <- ng_hhvr_humanize(ng_hhvr_first(implementation, ""))
  label <- if (nzchar(impl)) paste(role, impl, sep = " | ") else role
  paste0("<span class=\"badge ", ng_hhvr_role_class(role), "\">", ng_hhvr_escape(label), "</span>")
}

ng_hhvr_section <- function(title, body, id = NULL, subtitle = NULL) {
  attr <- if (is.null(id)) "" else paste0(" data-panel=\"", ng_hhvr_escape(id), "\"")
  paste0(
    "<section class=\"panel\"", attr, ">",
    "<div class=\"panel-head\"><h2>", ng_hhvr_escape(title), "</h2>",
    if (!is.null(subtitle)) paste0("<p>", ng_hhvr_escape(subtitle), "</p>") else "",
    "</div>",
    body,
    "</section>"
  )
}

ng_hhvr_index_svg <- function(method_summary) {
  if (!("mean_realized_index" %in% names(method_summary)) || !nrow(method_summary)) {
    return("<p class=\"empty\">Selection index response is unavailable.</p>")
  }
  dat <- method_summary[is.finite(suppressWarnings(as.numeric(method_summary$mean_realized_index))), , drop = FALSE]
  if (!nrow(dat)) return("<p class=\"empty\">Selection index response is unavailable.</p>")
  dat <- dat[order(suppressWarnings(as.numeric(dat$mean_realized_index))), , drop = FALSE]
  n <- nrow(dat)
  height <- 86 + 34 * n
  x0 <- 285
  width <- 900
  plot_w <- 520
  values <- suppressWarnings(as.numeric(dat$mean_realized_index))
  norm <- ng_hhvr_scale01(values)
  lines <- character(n)
  for (i in seq_len(n)) {
    y <- 48 + i * 34
    bar_w <- max(6, round(norm[[i]] * plot_w))
    role <- ng_hhvr_role_class(dat$benchmark_role[[i]])
    lines[[i]] <- paste0(
      "<text x=\"18\" y=\"", y + 5, "\" class=\"svg-label\">",
      ng_hhvr_escape(dat$method[[i]]), "</text>",
      "<rect x=\"", x0, "\" y=\"", y - 14, "\" width=\"", bar_w,
      "\" height=\"18\" rx=\"3\" class=\"bar ", role, "\"><title>",
      ng_hhvr_escape(paste(dat$method[[i]], ng_hhvr_num(values[[i]], 3))), "</title></rect>",
      "<text x=\"", x0 + bar_w + 10, "\" y=\"", y + 1, "\" class=\"svg-value\">",
      ng_hhvr_num(values[[i]], 3), "</text>"
    )
  }
  paste0(
    "<svg class=\"chart index-chart\" viewBox=\"0 0 ", width, " ", height,
    "\" role=\"img\" aria-label=\"Selection index response by method\">",
    "<text x=\"18\" y=\"28\" class=\"svg-title\">Mean realized selection index</text>",
    "<line x1=\"", x0, "\" y1=\"46\" x2=\"", x0 + plot_w, "\" y2=\"46\" class=\"axis\"/>",
    paste(lines, collapse = ""),
    "</svg>"
  )
}

ng_hhvr_trait_panel <- function(method_summary) {
  trait_metrics <- data.frame(
    metric = c("mean_realized_yield", "mean_realized_disease", "mean_realized_quality"),
    action = c("Increase yield", "Decrease disease", "Increase quality"),
    stringsAsFactors = FALSE
  )
  trait_metrics <- trait_metrics[trait_metrics$metric %in% names(method_summary), , drop = FALSE]
  if (!nrow(trait_metrics)) return("<p class=\"empty\">Trait response columns are unavailable.</p>")

  cards <- vector("list", nrow(trait_metrics))
  for (i in seq_len(nrow(trait_metrics))) {
    metric <- trait_metrics$metric[[i]]
    values <- suppressWarnings(as.numeric(method_summary[[metric]]))
    direction <- ng_hhvr_metric_direction(metric)
    oriented <- if (identical(direction, "minimize")) -values else values
    width <- pmax(4, round(100 * ng_hhvr_scale01(oriented)))
    ord <- order(-oriented, method_summary$method)
    rows <- character(length(ord))
    for (j in seq_along(ord)) {
      idx <- ord[[j]]
      rows[[j]] <- paste0(
        "<div class=\"trait-row\"><div class=\"trait-method\">",
        ng_hhvr_escape(method_summary$method[[idx]]), "</div>",
        "<div class=\"trait-bar\"><span style=\"width:", width[[idx]], "%\"></span></div>",
        "<div class=\"trait-value\">", ng_hhvr_num(values[[idx]], 2), "</div></div>"
      )
    }
    cards[[i]] <- paste0(
      "<div class=\"trait-card\"><div class=\"trait-action\">",
      ng_hhvr_escape(trait_metrics$action[[i]]), "</div>",
      paste(rows, collapse = ""),
      "</div>"
    )
  }
  paste0("<div class=\"trait-grid\">", paste(cards, collapse = ""), "</div>")
}

ng_hhvr_frontier_svg <- function(method_summary) {
  x_col <- if ("group_coancestry" %in% names(method_summary)) "group_coancestry" else "mean_pair_kinship"
  y_col <- "mean_realized_index"
  if (!(x_col %in% names(method_summary)) || !(y_col %in% names(method_summary))) {
    return("<p class=\"empty\">Gain-diversity frontier needs selection index and coancestry columns.</p>")
  }
  x <- suppressWarnings(as.numeric(method_summary[[x_col]]))
  y <- suppressWarnings(as.numeric(method_summary[[y_col]]))
  ok <- is.finite(x) & is.finite(y)
  dat <- method_summary[ok, , drop = FALSE]
  x <- x[ok]
  y <- y[ok]
  if (!nrow(dat)) return("<p class=\"empty\">Gain-diversity frontier is unavailable.</p>")

  width <- 900
  height <- 420
  left <- 88
  right <- 28
  top <- 42
  bottom <- 70
  plot_w <- width - left - right
  plot_h <- height - top - bottom
  sx <- left + ng_hhvr_scale01(x) * plot_w
  sy <- top + (1 - ng_hhvr_scale01(y)) * plot_h
  points <- character(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    role <- ng_hhvr_role_class(dat$benchmark_role[[i]])
    points[[i]] <- paste0(
      "<g class=\"point-wrap\"><circle cx=\"", round(sx[[i]], 1), "\" cy=\"", round(sy[[i]], 1),
      "\" r=\"8\" class=\"point ", role, "\"><title>",
      ng_hhvr_escape(paste(dat$method[[i]], "index", ng_hhvr_num(y[[i]], 3), x_col, ng_hhvr_num(x[[i]], 3))),
      "</title></circle>",
      "<text x=\"", round(sx[[i]] + 12, 1), "\" y=\"", round(sy[[i]] + 4, 1),
      "\" class=\"point-label\">", ng_hhvr_escape(dat$method[[i]]), "</text></g>"
    )
  }
  paste0(
    "<svg class=\"chart frontier-chart\" viewBox=\"0 0 ", width, " ", height,
    "\" role=\"img\" aria-label=\"Gain and diversity frontier\">",
    "<rect x=\"", left, "\" y=\"", top, "\" width=\"", plot_w, "\" height=\"", plot_h, "\" class=\"plot-bg\"/>",
    "<line x1=\"", left, "\" y1=\"", top + plot_h, "\" x2=\"", left + plot_w, "\" y2=\"", top + plot_h, "\" class=\"axis\"/>",
    "<line x1=\"", left, "\" y1=\"", top, "\" x2=\"", left, "\" y2=\"", top + plot_h, "\" class=\"axis\"/>",
    "<text x=\"", left, "\" y=\"24\" class=\"svg-title\">Higher gain is better; lower coancestry is left</text>",
    "<text x=\"", left + plot_w / 2, "\" y=\"392\" class=\"svg-axis\">",
    ng_hhvr_escape(ng_hhvr_metric_label(x_col)), " risk</text>",
    "<text x=\"18\" y=\"210\" transform=\"rotate(-90 18 210)\" class=\"svg-axis\">Selection index</text>",
    paste(points, collapse = ""),
    "<g class=\"legend\"><circle cx=\"670\" cy=\"24\" r=\"6\" class=\"point candidate\"/>",
    "<text x=\"682\" y=\"28\" class=\"legend-text\">candidate</text>",
    "<circle cx=\"760\" cy=\"24\" r=\"6\" class=\"point baseline\"/>",
    "<text x=\"772\" y=\"28\" class=\"legend-text\">baseline</text></g>",
    "</svg>"
  )
}

ng_hhvr_parent_counts <- function(selections, method = NULL) {
  selections <- as.data.frame(selections, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(selections)) || !nrow(selections)) {
    return(data.frame(parent = character(), uses = integer(), stringsAsFactors = FALSE))
  }
  if (!is.null(method) && "method" %in% names(selections)) {
    selections <- selections[selections$method == method, , drop = FALSE]
  }
  parents <- sort(unique(c(as.character(selections$parent1), as.character(selections$parent2))))
  parents <- parents[nzchar(parents) & !is.na(parents)]
  counts <- stats::setNames(integer(length(parents)), parents)
  for (i in seq_len(nrow(selections))) {
    p <- c(as.character(selections$parent1[[i]]), as.character(selections$parent2[[i]]))
    p <- p[p %in% names(counts)]
    counts[p] <- counts[p] + 1L
  }
  out <- data.frame(parent = names(counts), uses = as.integer(counts), stringsAsFactors = FALSE)
  out <- out[order(-out$uses, out$parent), , drop = FALSE]
  rownames(out) <- NULL
  out
}

ng_hhvr_parent_contribution_panel <- function(selections, method) {
  counts <- ng_hhvr_parent_counts(selections, method = method)
  if (!nrow(counts)) return("<p class=\"empty\">Parent contribution is unavailable.</p>")
  counts <- counts[seq_len(min(12L, nrow(counts))), , drop = FALSE]
  max_use <- max(counts$uses, na.rm = TRUE)
  rows <- character(nrow(counts))
  for (i in seq_len(nrow(counts))) {
    width <- if (is.finite(max_use) && max_use > 0) round(100 * counts$uses[[i]] / max_use) else 0
    rows[[i]] <- paste0(
      "<div class=\"parent-row\"><div class=\"parent-name\">", ng_hhvr_escape(counts$parent[[i]]),
      "</div><div class=\"parent-use\"><span style=\"width:", width, "%\"></span></div>",
      "<div class=\"parent-count\">", counts$uses[[i]], "</div></div>"
    )
  }
  paste0(
    "<div class=\"focus-note\">Recommended method shown: <strong>", ng_hhvr_escape(method), "</strong></div>",
    "<div class=\"parent-list\">", paste(rows, collapse = ""), "</div>"
  )
}

ng_hhvr_pair_key <- function(p1, p2) {
  p <- sort(c(as.character(p1), as.character(p2)))
  paste(p, collapse = "||")
}

ng_hhvr_matrix_panel <- function(selections, method, max_parents = 12L) {
  selections <- as.data.frame(selections, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(selections)) || !nrow(selections)) {
    return("<p class=\"empty\">Mate allocation matrix is unavailable.</p>")
  }
  if ("method" %in% names(selections)) selections <- selections[selections$method == method, , drop = FALSE]
  if (!nrow(selections)) return("<p class=\"empty\">No selected crosses are available for the recommended method.</p>")
  counts <- ng_hhvr_parent_counts(selections)
  parents <- counts$parent[seq_len(min(max_parents, nrow(counts)))]
  pair_counts <- list()
  for (i in seq_len(nrow(selections))) {
    key <- ng_hhvr_pair_key(selections$parent1[[i]], selections$parent2[[i]])
    pair_counts[[key]] <- if (is.null(pair_counts[[key]])) 1L else pair_counts[[key]] + 1L
  }
  max_pair <- max(unlist(pair_counts, use.names = FALSE), na.rm = TRUE)
  header <- paste0("<th></th>", paste0("<th>", ng_hhvr_escape(parents), "</th>", collapse = ""))
  rows <- character(length(parents))
  for (i in seq_along(parents)) {
    cells <- character(length(parents))
    for (j in seq_along(parents)) {
      if (i == j) {
        cells[[j]] <- "<td class=\"diag\"></td>"
      } else {
        key <- ng_hhvr_pair_key(parents[[i]], parents[[j]])
        val <- if (is.null(pair_counts[[key]])) 0L else pair_counts[[key]]
        alpha <- if (max_pair > 0) 0.10 + 0.70 * val / max_pair else 0
        cells[[j]] <- paste0(
          "<td class=\"matrix-cell\" style=\"background:rgba(43, 125, 107, ", sprintf("%.2f", alpha),
          ")\">", if (val > 0L) val else "", "</td>"
        )
      }
    }
    rows[[i]] <- paste0("<tr><th>", ng_hhvr_escape(parents[[i]]), "</th>", paste(cells, collapse = ""), "</tr>")
  }
  paste0(
    "<div class=\"matrix-wrap\"><table class=\"matrix\"><thead><tr>", header, "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

ng_hhvr_ranked_decisions_panel <- function(selections, method, n = 12L) {
  selections <- as.data.frame(selections, stringsAsFactors = FALSE)
  required <- c("parent1", "parent2")
  if (!all(required %in% names(selections)) || !nrow(selections)) {
    return("<p class=\"empty\">Ranked crossing decisions are unavailable.</p>")
  }
  if ("method" %in% names(selections)) selections <- selections[selections$method == method, , drop = FALSE]
  if (!nrow(selections)) return("<p class=\"empty\">No selected crosses are available for the recommended method.</p>")
  order_cols <- intersect(c("scenario", "n_parents", "rep", "selection_rank"), names(selections))
  if (length(order_cols)) {
    selections <- selections[do.call(order, selections[order_cols]), , drop = FALSE]
  }
  selections <- selections[seq_len(min(n, nrow(selections))), , drop = FALSE]
  headers <- c("Rank", "Cross", "Score", "Index", "Yield", "Disease", "Quality", "Kinship")
  rows <- character(nrow(selections))
  for (i in seq_len(nrow(selections))) {
    rank <- if ("selection_rank" %in% names(selections)) selections$selection_rank[[i]] else i
    cross <- paste(selections$parent1[[i]], selections$parent2[[i]], sep = " x ")
    value <- function(col, digits = 2) if (col %in% names(selections)) ng_hhvr_num(selections[[col]][[i]], digits) else ""
    cells <- c(
      ng_hhvr_escape(rank),
      ng_hhvr_escape(cross),
      value("multi_trait_score", 2),
      value("realized_index", 2),
      value("realized_yield", 1),
      value("realized_disease", 1),
      value("realized_quality", 1),
      value("pair_kinship", 3)
    )
    rows[[i]] <- paste0("<tr>", paste0("<td>", cells, "</td>", collapse = ""), "</tr>")
  }
  paste0(
    "<div class=\"table-wrap\"><table class=\"decision-table\"><thead><tr>",
    paste0("<th>", ng_hhvr_escape(headers), "</th>", collapse = ""),
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

ng_hhvr_winner_panel <- function(winner_summary) {
  winner_summary <- as.data.frame(winner_summary, stringsAsFactors = FALSE)
  if (!all(c("metric", "method") %in% names(winner_summary)) || !nrow(winner_summary)) {
    return("<p class=\"empty\">Winner summary is unavailable.</p>")
  }
  rows <- character(nrow(winner_summary))
  for (i in seq_len(nrow(winner_summary))) {
    parent_size <- if ("n_parents" %in% names(winner_summary)) winner_summary$n_parents[[i]] else ""
    value <- if ("value" %in% names(winner_summary)) ng_hhvr_num(winner_summary$value[[i]], 3) else ""
    rows[[i]] <- paste0(
      "<tr><td>", ng_hhvr_escape(parent_size), "</td><td>",
      ng_hhvr_escape(ng_hhvr_metric_label(winner_summary$metric[[i]])), "</td><td>",
      ng_hhvr_escape(ng_hhvr_metric_direction(winner_summary$metric[[i]])), "</td><td>",
      ng_hhvr_escape(winner_summary$method[[i]]), "</td><td>", value, "</td></tr>"
    )
  }
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>",
    "<th>Parents</th><th>Decision metric</th><th>Direction</th><th>Best method</th><th>Value</th>",
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

ng_hhvr_comparison_panel <- function(comparisons) {
  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)
  if (!all(c("method", "baseline_method", "better_than_baseline") %in% names(comparisons)) || !nrow(comparisons)) {
    return("<p class=\"empty\">Candidate-versus-baseline comparisons are unavailable.</p>")
  }
  keys <- unique(comparisons[, c("method", "baseline_method"), drop = FALSE])
  rows <- character(nrow(keys))
  for (i in seq_len(nrow(keys))) {
    part <- comparisons[
      comparisons$method == keys$method[[i]] & comparisons$baseline_method == keys$baseline_method[[i]],
      , drop = FALSE
    ]
    better <- mean(as.logical(part$better_than_baseline), na.rm = TRUE)
    median_delta <- if ("delta" %in% names(part)) stats::median(suppressWarnings(as.numeric(part$delta)), na.rm = TRUE) else NA_real_
    rows[[i]] <- paste0(
      "<tr><td>", ng_hhvr_escape(keys$method[[i]]), "</td><td>",
      ng_hhvr_escape(keys$baseline_method[[i]]), "</td><td>",
      ng_hhvr_num(100 * better, 0), "%</td><td>", ng_hhvr_num(median_delta, 3), "</td></tr>"
    )
  }
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>",
    "<th>Candidate</th><th>Baseline</th><th>Metrics better</th><th>Median signed delta</th>",
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

ng_hhvr_caveat_panel <- function(method_registry, comparisons) {
  method_registry <- as.data.frame(method_registry, stringsAsFactors = FALSE)
  rows <- character()
  if (nrow(method_registry) && all(c("label", "implementation", "exact_external_status") %in% names(method_registry))) {
    style_proxy <- method_registry[method_registry$implementation == "style_proxy", , drop = FALSE]
    if (nrow(style_proxy)) {
      rows <- character(nrow(style_proxy))
      for (i in seq_len(nrow(style_proxy))) {
        reason <- if ("fallback_reason" %in% names(style_proxy)) style_proxy$fallback_reason[[i]] else ""
        tool <- if ("external_tool" %in% names(style_proxy)) style_proxy$external_tool[[i]] else ""
        rows[[i]] <- paste0(
          "<tr><td>", ng_hhvr_escape(style_proxy$label[[i]]), "</td><td>",
          ng_hhvr_escape(tool), "</td><td>style proxy</td><td>",
          ng_hhvr_escape(ng_hhvr_humanize(style_proxy$exact_external_status[[i]])), "</td><td>",
          ng_hhvr_escape(reason), "</td></tr>"
        )
      }
    }
  }
  if (!length(rows) && nrow(comparisons)) {
    baseline_cols <- intersect(c(
      "baseline_method", "baseline_external_tool", "baseline_implementation",
      "baseline_exact_external_status", "baseline_fallback_reason"
    ), names(comparisons))
    if (length(baseline_cols) >= 3L) {
      baselines <- unique(comparisons[, baseline_cols, drop = FALSE])
      rows <- character(nrow(baselines))
      for (i in seq_len(nrow(baselines))) {
        rows[[i]] <- paste0(
          "<tr><td>", ng_hhvr_escape(baselines$baseline_method[[i]]), "</td><td>",
          ng_hhvr_escape(if ("baseline_external_tool" %in% names(baselines)) baselines$baseline_external_tool[[i]] else ""),
          "</td><td>", ng_hhvr_escape(ng_hhvr_humanize(if ("baseline_implementation" %in% names(baselines)) baselines$baseline_implementation[[i]] else "")),
          "</td><td>", ng_hhvr_escape(ng_hhvr_humanize(if ("baseline_exact_external_status" %in% names(baselines)) baselines$baseline_exact_external_status[[i]] else "")),
          "</td><td>", ng_hhvr_escape(if ("baseline_fallback_reason" %in% names(baselines)) baselines$baseline_fallback_reason[[i]] else ""),
          "</td></tr>"
        )
      }
    }
  }
  if (!length(rows)) return("<p class=\"empty\">No benchmark caveats were reported.</p>")
  paste0(
    "<div class=\"table-wrap\"><table><thead><tr>",
    "<th>Method</th><th>Tool</th><th>Evidence type</th><th>Status</th><th>Reason</th>",
    "</tr></thead><tbody>", paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

ng_hhvr_context_value <- function(summary, col, empty = "not reported") {
  if (!(col %in% names(summary))) return(empty)
  values <- unique(as.character(summary[[col]]))
  values <- values[nzchar(values) & !is.na(values)]
  if (!length(values)) return(empty)
  paste(values, collapse = ", ")
}

ng_hhvr_css <- function() {
  paste0(
    "<style>",
    ":root{--ink:#1f2933;--muted:#64707d;--line:#d7dee6;--panel:#ffffff;--bg:#f5f7fa;--candidate:#2b7d6b;--baseline:#b85c38;--accent:#2f5f98;--warn:#b7791f;}",
    "*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font-family:Arial,Helvetica,sans-serif;line-height:1.42;}",
    ".shell{max-width:1280px;margin:0 auto;padding:24px 22px 42px}.hero{display:grid;grid-template-columns:minmax(0,1.7fr) minmax(280px,.8fr);gap:18px;align-items:end;margin-bottom:18px;}",
    "h1{font-size:30px;line-height:1.1;margin:0 0 10px}h2{font-size:18px;margin:0}p{margin:0}.sub{color:var(--muted);max-width:780px}.stamp{color:var(--muted);font-size:13px;text-align:right}",
    ".kpis{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));gap:10px;margin:14px 0 18px}.kpi{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:12px}.kpi-label{font-size:12px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}.kpi-value{font-size:22px;font-weight:700;margin-top:4px}.kpi-note{font-size:12px;color:var(--muted);margin-top:2px}",
    ".context-strip{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:10px;background:#fff;border:1px solid var(--line);border-radius:8px;padding:10px;margin:-4px 0 18px}.context-item span{display:block;font-size:11px;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}.context-item strong{display:block;margin-top:2px;font-size:13px;overflow-wrap:anywhere}",
    ".grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:14px}.panel{background:var(--panel);border:1px solid var(--line);border-radius:8px;padding:16px;min-width:0}.panel.wide{grid-column:1/-1}.panel-head{display:flex;align-items:flex-start;justify-content:space-between;gap:12px;margin-bottom:12px}.panel-head p{font-size:13px;color:var(--muted);max-width:640px}",
    ".badge{display:inline-block;border-radius:999px;padding:2px 8px;font-size:12px;border:1px solid var(--line);white-space:nowrap}.badge.candidate{background:#e7f3ef;color:#155a4a}.badge.baseline{background:#f7ece6;color:#873f23}.badge.neutral{background:#edf1f5;color:#425160}",
    ".chart{display:block;width:100%;height:auto}.axis{stroke:#8794a3;stroke-width:1}.plot-bg{fill:#fbfcfd;stroke:#e3e8ee}.svg-title{font-size:15px;font-weight:700;fill:#263442}.svg-label{font-size:12px;fill:#263442}.svg-value,.svg-axis,.legend-text{font-size:12px;fill:#5d6875}.bar.candidate{fill:var(--candidate)}.bar.baseline{fill:var(--baseline)}.bar.neutral{fill:var(--accent)}.point{stroke:white;stroke-width:2}.point.candidate{fill:var(--candidate)}.point.baseline{fill:var(--baseline)}.point.neutral{fill:var(--accent)}.point-label{font-size:11px;fill:#33404d}",
    ".trait-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px}.trait-card{border:1px solid var(--line);border-radius:7px;padding:10px;background:#fbfcfd}.trait-action{font-weight:700;margin-bottom:8px}.trait-row,.parent-row{display:grid;grid-template-columns:minmax(120px,1fr) minmax(80px,1.3fr) 48px;gap:8px;align-items:center;margin:6px 0}.trait-method,.parent-name{font-size:12px;overflow-wrap:anywhere}.trait-bar,.parent-use{height:10px;background:#e8edf2;border-radius:999px;overflow:hidden}.trait-bar span{display:block;height:100%;background:var(--accent)}.parent-use span{display:block;height:100%;background:var(--candidate)}.trait-value,.parent-count{font-size:12px;text-align:right;color:var(--muted)}",
    ".focus-note{font-size:13px;color:var(--muted);margin-bottom:8px}.matrix-wrap,.table-wrap{overflow:auto}.matrix{border-collapse:collapse;min-width:520px}.matrix th,.matrix td{border:1px solid #e0e6ed;min-width:34px;height:30px;text-align:center;font-size:12px}.matrix th{background:#f3f6f9;font-weight:600}.matrix .diag{background:#f0f3f6}.matrix-cell{font-weight:700;color:#143e35}",
    "table{width:100%;border-collapse:collapse}th,td{border-bottom:1px solid #e4e9ef;padding:8px;text-align:left;font-size:13px;vertical-align:top}th{color:#4b5866;background:#f5f7fa;font-weight:700}.decision-table td:first-child{font-weight:700}.empty{color:var(--muted);font-size:13px;padding:8px 0}",
    "@media(max-width:900px){.hero,.grid,.kpis,.trait-grid,.context-strip{grid-template-columns:1fr}.stamp{text-align:left}.panel-head{display:block}.trait-row,.parent-row{grid-template-columns:minmax(100px,1fr) minmax(90px,1.2fr) 44px}}",
    "</style>"
  )
}

ng_hhvr_html <- function(summary,
                         selections,
                         comparisons,
                         winner_summary,
                         method_registry,
                         generated_at = Sys.time()) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  selections <- as.data.frame(selections, stringsAsFactors = FALSE)
  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)
  winner_summary <- as.data.frame(winner_summary, stringsAsFactors = FALSE)
  method_registry <- as.data.frame(method_registry, stringsAsFactors = FALSE)
  method_summary <- ng_hhvr_method_summary(summary)
  recommended <- ng_hhvr_recommended_method(method_summary)

  n_scenarios <- if ("scenario" %in% names(summary)) length(unique(summary$scenario)) else NA_integer_
  n_parent_sizes <- if ("n_parents" %in% names(summary)) length(unique(summary$n_parents)) else NA_integer_
  n_methods <- length(unique(summary$method))
  n_comparisons <- nrow(comparisons)
  better_rate <- if ("better_than_baseline" %in% names(comparisons) && nrow(comparisons)) {
    paste0(ng_hhvr_num(100 * mean(as.logical(comparisons$better_than_baseline), na.rm = TRUE), 0), "%")
  } else {
    ""
  }

  recommended_row <- method_summary[method_summary$method == recommended, , drop = FALSE]
  recommended_badge <- if (nrow(recommended_row)) {
    ng_hhvr_role_badge(recommended_row$benchmark_role[[1]], recommended_row$implementation[[1]])
  } else {
    ""
  }

  kpis <- paste0(
    "<div class=\"kpis\">",
    ng_hhvr_kpi_card("Scenarios", ng_hhvr_int(n_scenarios), "crop/genome contexts"),
    ng_hhvr_kpi_card("Parent sizes", ng_hhvr_int(n_parent_sizes), "selection scale"),
    ng_hhvr_kpi_card("Methods", ng_hhvr_int(n_methods), "candidate and baseline"),
    ng_hhvr_kpi_card("Comparisons", ng_hhvr_int(n_comparisons), "metric-level contrasts"),
    ng_hhvr_kpi_card("Better rate", better_rate, "candidate vs baseline"),
    "</div>"
  )
  context <- paste0(
    "<div class=\"context-strip\">",
    "<div class=\"context-item\"><span>Scenarios</span><strong>",
    ng_hhvr_escape(ng_hhvr_context_value(summary, "scenario")), "</strong></div>",
    "<div class=\"context-item\"><span>Crop</span><strong>",
    ng_hhvr_escape(ng_hhvr_context_value(summary, "crop")), "</strong></div>",
    "<div class=\"context-item\"><span>Harness model</span><strong>",
    ng_hhvr_escape(ng_hhvr_context_value(summary, "harness_model")), "</strong></div>",
    "<div class=\"context-item\"><span>Parent sizes</span><strong>",
    ng_hhvr_escape(ng_hhvr_context_value(summary, "n_parents")), "</strong></div>",
    "</div>"
  )

  paste0(
    "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
    "<title>Head-to-Head Multi-Trait Decision Report</title>",
    ng_hhvr_css(),
    "</head><body><main class=\"shell\">",
    "<div class=\"hero\"><div><h1>Head-to-Head Multi-Trait Decision Report</h1>",
    "<p class=\"sub\">Breeder-facing view of multi-trait gain, disease reduction, quality response, diversity risk, parent use, and selected mate allocations. The backend remains the quantitative genetics package; this report is the frontend decision layer over the benchmark outputs.</p></div>",
    "<div class=\"stamp\"><div>Generated: ", ng_hhvr_escape(format(generated_at, "%Y-%m-%d %H:%M:%S %Z")), "</div>",
    "<div>Recommended focus: <strong>", ng_hhvr_escape(recommended), "</strong> ", recommended_badge, "</div></div></div>",
    kpis,
    context,
    "<div class=\"grid\">",
    ng_hhvr_section(
      "Selection Index Response",
      ng_hhvr_index_svg(method_summary),
      id = "selection-index-response",
      subtitle = "Mean realized selection index by method. Higher is better."
    ),
    ng_hhvr_section(
      "Trait Direction Response",
      ng_hhvr_trait_panel(method_summary),
      id = "trait-direction-response",
      subtitle = "Each trait is oriented in the breeder's desired direction."
    ),
    ng_hhvr_section(
      "Gain-Diversity Frontier",
      ng_hhvr_frontier_svg(method_summary),
      id = "gain-diversity-frontier",
      subtitle = "Methods in the upper-left combine high realized gain with lower coancestry risk."
    ),
    ng_hhvr_section(
      "Parent Contribution",
      ng_hhvr_parent_contribution_panel(selections, recommended),
      id = "parent-contribution",
      subtitle = "Parent-use pressure in the recommended mating plan."
    ),
    ng_hhvr_section(
      "Mate Allocation Matrix",
      ng_hhvr_matrix_panel(selections, recommended),
      id = "mate-allocation-matrix",
      subtitle = "Selected crosses among the most-used parents for the recommended method."
    ),
    ng_hhvr_section(
      "Ranked Crossing Decisions",
      ng_hhvr_ranked_decisions_panel(selections, recommended),
      id = "ranked-crossing-decisions",
      subtitle = "Top selected crosses with realized trait response diagnostics."
    ),
    "<section class=\"panel wide\">",
    "<div class=\"panel-head\"><h2>Winner Summary</h2><p>Best method by metric and parent-size setting.</p></div>",
    ng_hhvr_winner_panel(winner_summary),
    "</section>",
    "<section class=\"panel wide\">",
    "<div class=\"panel-head\"><h2>Candidate vs Baseline Evidence</h2><p>Signed deltas follow the metric direction: positive means the candidate is better than the baseline.</p></div>",
    ng_hhvr_comparison_panel(comparisons),
    "</section>",
    "<section class=\"panel wide\">",
    "<div class=\"panel-head\"><h2>Benchmark Evidence Caveats</h2><p>External-tool baselines can be exact runs or style proxies depending on what is available in the local benchmark environment.</p></div>",
    ng_hhvr_caveat_panel(method_registry, comparisons),
    "</section>",
    "</div></main></body></html>"
  )
}

ng_write_head_to_head_visual_report <- function(output_path,
                                                summary = NULL,
                                                selections = NULL,
                                                comparisons = NULL,
                                                winner_summary = NULL,
                                                method_registry = NULL,
                                                results_dir = NULL,
                                                prefix = NULL,
                                                generated_at = Sys.time()) {
  if (missing(output_path) || is.null(output_path) || !nzchar(as.character(output_path[[1]]))) {
    ng_stop("output_path is required")
  }
  output_path <- as.character(output_path[[1]])
  if (is.null(summary)) summary <- ng_hhvr_read_csv(results_dir, prefix, "summary", required = TRUE)
  if (is.null(selections)) selections <- ng_hhvr_read_csv(results_dir, prefix, "selections", required = FALSE)
  if (is.null(comparisons)) comparisons <- ng_hhvr_read_csv(results_dir, prefix, "comparisons", required = FALSE)
  if (is.null(winner_summary)) winner_summary <- ng_hhvr_read_csv(results_dir, prefix, "winner_summary", required = FALSE)
  if (is.null(method_registry)) method_registry <- ng_hhvr_read_csv(results_dir, prefix, "method_registry", required = FALSE)

  html <- ng_hhvr_html(
    summary = summary,
    selections = selections,
    comparisons = comparisons,
    winner_summary = winner_summary,
    method_registry = method_registry,
    generated_at = generated_at
  )
  parent <- dirname(output_path)
  if (nzchar(parent) && !dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_path, useBytes = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
