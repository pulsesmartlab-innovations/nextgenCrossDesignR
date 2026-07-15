local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

read_optional_csv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.csv(path, stringsAsFactors = FALSE)
}

fmt_num <- function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "")
}

code <- function(x) {
  x <- as.character(x)
  ifelse(nzchar(x), paste0("`", x, "`"), "")
}

md_table <- function(headers, rows) {
  if (!length(rows)) return(character())
  body <- vapply(rows, function(row) {
    paste0("| ", paste(as.character(row), collapse = " | "), " |")
  }, character(1))
  c(
    paste0("| ", paste(headers, collapse = " | "), " |"),
    paste0("| ", paste(rep("---", length(headers)), collapse = " | "), " |"),
    body
  )
}

path_for <- function(results_dir, prefix, suffix) {
  file.path(results_dir, paste0(prefix, suffix))
}

status_sentence <- function(status) {
  required <- c("popvar_status", "simple_status", "gms_status")
  if (!nrow(status)) return("Family status: not available.")
  if (all(required %in% names(status)) &&
      all(status$popvar_status == "ok") &&
      all(status$simple_status == "ok") &&
      all(status$gms_status == "ok")) {
    return(sprintf(
      "Family status: all %d calibration runs completed with PopVar, SimpleMating, and GMS status ok.",
      nrow(status)
    ))
  }

  summarize_one <- function(col) {
    if (!(col %in% names(status))) return(paste0(col, "=missing"))
    counts <- table(status[[col]], useNA = "ifany")
    paste(paste0(names(counts), ":", as.integer(counts)), collapse = ",")
  }
  paste(
    "Family status:",
    paste(vapply(required, summarize_one, character(1)), collapse = "; ")
  )
}

winner_rows <- function(winners, metric) {
  if (!nrow(winners) || !all(c("n_parents", "metric", "method", "value") %in% names(winners))) {
    return(list())
  }
  d <- winners[winners$metric == metric, , drop = FALSE]
  if (!nrow(d)) return(list())
  d$n_parents <- suppressWarnings(as.integer(d$n_parents))
  d <- d[order(d$n_parents), , drop = FALSE]
  lapply(seq_len(nrow(d)), function(i) {
    c(d$n_parents[[i]], code(d$method[[i]]), fmt_num(d$value[[i]]), d$effect_training_n[[i]])
  })
}

all_metric_winner_rows <- function(winners) {
  if (!nrow(winners) || !all(c("n_parents", "metric", "method", "value") %in% names(winners))) {
    return(list())
  }
  keep <- winners$metric %in% c("mean_gv", "top10_gv", "max_gv")
  d <- winners[keep, , drop = FALSE]
  if (!nrow(d)) return(list())
  d$n_parents <- suppressWarnings(as.integer(d$n_parents))
  d <- d[order(d$n_parents, match(d$metric, c("mean_gv", "top10_gv", "max_gv"))), , drop = FALSE]
  lapply(seq_len(nrow(d)), function(i) {
    c(d$n_parents[[i]], d$metric[[i]], code(d$method[[i]]), fmt_num(d$value[[i]]))
  })
}

router_gap_rows <- function(overall, winners) {
  needed <- c("n_parents", "method", "top10_gv")
  if (!nrow(overall) || !all(needed %in% names(overall))) return(list())
  frontier <- winners[winners$metric == "top10_gv", , drop = FALSE]
  if (!nrow(frontier)) return(list())
  overall$n_parents <- suppressWarnings(as.integer(overall$n_parents))
  overall$top10_gv <- suppressWarnings(as.numeric(overall$top10_gv))
  frontier$n_parents <- suppressWarnings(as.integer(frontier$n_parents))
  frontier$value <- suppressWarnings(as.numeric(frontier$value))

  rows <- list()
  for (n in sort(unique(frontier$n_parents))) {
    f <- frontier[frontier$n_parents == n, , drop = FALSE]
    d <- overall[overall$n_parents == n & grepl("^ng_meta_router", overall$method), , drop = FALSE]
    if (!nrow(f) || !nrow(d)) next
    d <- d[order(-d$top10_gv), , drop = FALSE]
    gap <- d$top10_gv[[1]] - f$value[[1]]
    rows[[length(rows) + 1L]] <- c(
      n,
      code(f$method[[1]]),
      fmt_num(f$value[[1]]),
      code(d$method[[1]]),
      fmt_num(d$top10_gv[[1]]),
      fmt_num(gap)
    )
  }
  rows
}

calibration_rows <- function(metrics) {
  if (!nrow(metrics) || !("target" %in% names(metrics)) || !("spearman" %in% names(metrics))) {
    return(list())
  }
  score_col <- if ("score_col" %in% names(metrics)) "score_col" else if ("score" %in% names(metrics)) "score" else NA_character_
  if (is.na(score_col) || !("n_parents" %in% names(metrics))) return(list())
  d <- metrics[metrics$target %in% c("realized_top10", "realized_var"), , drop = FALSE]
  if (!nrow(d)) return(list())
  d$n_parents <- suppressWarnings(as.integer(d$n_parents))
  d$spearman <- suppressWarnings(as.numeric(d$spearman))
  if ("top_overlap" %in% names(d)) {
    d$top_overlap <- suppressWarnings(as.numeric(d$top_overlap))
  } else {
    d$top_overlap <- NA_real_
  }

  rows <- list()
  for (n in sort(unique(d$n_parents))) {
    for (target in c("realized_top10", "realized_var")) {
      dt <- d[d$n_parents == n & d$target == target & is.finite(d$spearman), , drop = FALSE]
      if (!nrow(dt)) next
      idx <- which.max(dt$spearman)
      rows[[length(rows) + 1L]] <- c(
        n,
        target,
        code(dt[[score_col]][[idx]]),
        fmt_num(dt$spearman[[idx]]),
        fmt_num(dt$top_overlap[[idx]])
      )
    }
  }
  rows
}

results_dir <- env_chr("NG_VALIDATION_REPORT_RESULTS_DIR", "results")
family_prefix <- env_chr("NG_VALIDATION_REPORT_FAMILY_PREFIX", "diagnostic_family_5k")
allocator_prefix <- env_chr("NG_VALIDATION_REPORT_ALLOCATOR_PREFIX", "diagnostic_allocator_5k")
grid_prefix <- env_chr("NG_VALIDATION_REPORT_GRID_PREFIX", "diagnostic_parent_grid_5k")
out_path <- env_chr("NG_VALIDATION_REPORT_OUT", file.path(results_dir, paste0(grid_prefix, "_validation_report.md")))

family_status <- read_optional_csv(path_for(results_dir, family_prefix, "_status.csv"))
family_metrics <- read_optional_csv(path_for(results_dir, family_prefix, "_metric_summary_avg.csv"))
allocator_winners <- read_optional_csv(path_for(results_dir, allocator_prefix, "_winner_summary.csv"))
grid_winners <- read_optional_csv(path_for(results_dir, grid_prefix, "_winner_summary.csv"))
grid_overall <- read_optional_csv(path_for(results_dir, grid_prefix, "_overall_avg.csv"))

if (!nrow(grid_winners)) {
  stop("Missing parent-size grid winner summary for prefix: ", grid_prefix, call. = FALSE)
}

lines <- c(
  "# Framework Validation Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Results directory: `", results_dir, "`"),
  paste0("Family prefix: `", family_prefix, "`"),
  paste0("Allocator prefix: `", allocator_prefix, "`"),
  paste0("Parent grid prefix: `", grid_prefix, "`"),
  "",
  status_sentence(family_status),
  "",
  "## Parent-size top10 frontier",
  "",
  md_table(
    c("Parents", "Top10 winner", "Top10", "Effect training N"),
    winner_rows(grid_winners, "top10_gv")
  ),
  "",
  "## Mean/top/max winners",
  "",
  md_table(
    c("Parents", "Metric", "Winner", "Value"),
    all_metric_winner_rows(grid_winners)
  )
)

gap <- router_gap_rows(grid_overall, grid_winners)
if (length(gap)) {
  lines <- c(lines, "", "## Router context", "", md_table(
    c("Parents", "Top10 winner", "Winner top10", "Best router", "Router top10", "Router gap"),
    gap
  ))
}

alloc <- winner_rows(allocator_winners, "top10_gv")
if (length(alloc)) {
  lines <- c(lines, "", "## Allocator crosscheck top10 winners", "", md_table(
    c("Parents", "Top10 winner", "Top10", "Effect training N"),
    alloc
  ))
}

cal <- calibration_rows(family_metrics)
if (length(cal)) {
  lines <- c(lines, "", "## Family calibration signal", "", md_table(
    c("Parents", "Target", "Best score", "Spearman", "Top overlap"),
    cal
  ))
}

lines <- c(
  lines,
  "",
  "## Interpretation",
  "",
  "- Do not promote a universal winner from mixed parent-size evidence.",
  "- Treat the parent-size top10 frontier as the current empirical decision rule until a new method beats it in replicated validation.",
  "- Prioritize method changes where the router gap is negative or where allocator crosscheck winners disagree with the multi-cycle grid.",
  "- Keep PopVar, SimpleMating, var-simple, PMV-balanced, recombination-GEBV, meta-selector, and meta-portfolio as controls in the next validation run."
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
writeLines(lines, out_path, useBytes = TRUE)
message("Wrote validation report: ", out_path)
