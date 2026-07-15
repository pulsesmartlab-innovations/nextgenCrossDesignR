ng_hhdj_required_columns <- function(data, cols, name) {
  missing <- setdiff(cols, names(data))
  if (length(missing)) ng_stop(name, " missing dashboard JSON columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

ng_hhdj_unique_values <- function(data, col) {
  if (!(col %in% names(data))) return(character())
  values <- unique(as.character(data[[col]]))
  values[nzchar(values) & !is.na(values)]
}

ng_hhdj_context <- function(summary) {
  list(
    scenarios = ng_hhdj_unique_values(summary, "scenario"),
    crops = ng_hhdj_unique_values(summary, "crop"),
    harness_models = ng_hhdj_unique_values(summary, "harness_model"),
    parent_sizes = sort(unique(suppressWarnings(as.integer(summary$n_parents[is.finite(summary$n_parents)])))),
    reps = sort(unique(suppressWarnings(as.integer(summary$rep[is.finite(summary$rep)]))))
  )
}

ng_hhdj_source_files <- function(prefix) {
  suffixes <- c("summary", "selections", "comparisons", "winner_summary", "method_registry", "config")
  stats::setNames(paste0(prefix, "_", suffixes, ".csv"), suffixes)
}

ng_hhdj_sanitize <- function(x) {
  if (is.data.frame(x)) {
    out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
    for (name in names(out)) out[[name]] <- ng_hhdj_sanitize(out[[name]])
    rownames(out) <- NULL
    return(out)
  }
  if (is.list(x)) return(lapply(x, ng_hhdj_sanitize))
  if (is.numeric(x)) {
    x[!is.finite(x)] <- NA_real_
    return(x)
  }
  x
}

ng_head_to_head_dashboard_payload <- function(summary,
                                              selections = NULL,
                                              comparisons = NULL,
                                              winner_summary = NULL,
                                              method_registry = NULL,
                                              config = NULL,
                                              prefix = "head_to_head_benchmark",
                                              generated_at = Sys.time()) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  selections <- if (is.null(selections)) data.frame(stringsAsFactors = FALSE) else as.data.frame(selections, stringsAsFactors = FALSE)
  comparisons <- if (is.null(comparisons)) data.frame(stringsAsFactors = FALSE) else as.data.frame(comparisons, stringsAsFactors = FALSE)
  winner_summary <- if (is.null(winner_summary)) data.frame(stringsAsFactors = FALSE) else as.data.frame(winner_summary, stringsAsFactors = FALSE)
  method_registry <- if (is.null(method_registry)) data.frame(stringsAsFactors = FALSE) else as.data.frame(method_registry, stringsAsFactors = FALSE)
  config <- if (is.null(config)) data.frame(stringsAsFactors = FALSE) else as.data.frame(config, stringsAsFactors = FALSE)

  ng_hhdj_required_columns(summary, c("method", "benchmark_role"), "summary")
  if (nrow(selections)) ng_hhdj_required_columns(selections, c("method", "parent1", "parent2"), "selections")
  if (nrow(comparisons)) {
    ng_hhdj_required_columns(
      comparisons,
      c("metric", "direction", "method", "baseline_method", "value", "baseline_value", "delta"),
      "comparisons"
    )
  }
  if (nrow(winner_summary)) ng_hhdj_required_columns(winner_summary, c("metric", "method"), "winner_summary")
  if (nrow(method_registry)) ng_hhdj_required_columns(method_registry, c("label", "role", "implementation"), "method_registry")

  method_summary <- ng_hhvr_method_summary(summary)
  payload <- list(
    schema_version = "ng_head_to_head_dashboard.v1",
    generated_at = format(as.POSIXct(generated_at, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    source = list(
      prefix = as.character(prefix[[1]]),
      files = as.list(ng_hhdj_source_files(prefix))
    ),
    metric_directions = as.list(ng_head_to_head_metric_directions()),
    context = ng_hhdj_context(summary),
    recommended_method = ng_hhvr_recommended_method(method_summary),
    method_summary = method_summary,
    summary = summary,
    selections = selections,
    comparisons = comparisons,
    winner_summary = winner_summary,
    method_registry = method_registry,
    config = config
  )
  ng_hhdj_sanitize(payload)
}

ng_write_head_to_head_dashboard_json <- function(output_path,
                                                 summary = NULL,
                                                 selections = NULL,
                                                 comparisons = NULL,
                                                 winner_summary = NULL,
                                                 method_registry = NULL,
                                                 config = NULL,
                                                 results_dir = NULL,
                                                 prefix = "head_to_head_benchmark",
                                                 generated_at = Sys.time()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to export dashboard JSON; install jsonlite or skip JSON export")
  }
  if (missing(output_path) || is.null(output_path) || !nzchar(as.character(output_path[[1]]))) {
    ng_stop("output_path is required")
  }
  output_path <- as.character(output_path[[1]])
  if (is.null(summary)) summary <- ng_hhvr_read_csv(results_dir, prefix, "summary", required = TRUE)
  if (is.null(selections)) selections <- ng_hhvr_read_csv(results_dir, prefix, "selections", required = FALSE)
  if (is.null(comparisons)) comparisons <- ng_hhvr_read_csv(results_dir, prefix, "comparisons", required = FALSE)
  if (is.null(winner_summary)) winner_summary <- ng_hhvr_read_csv(results_dir, prefix, "winner_summary", required = FALSE)
  if (is.null(method_registry)) method_registry <- ng_hhvr_read_csv(results_dir, prefix, "method_registry", required = FALSE)
  if (is.null(config)) config <- ng_hhvr_read_csv(results_dir, prefix, "config", required = FALSE)

  payload <- ng_head_to_head_dashboard_payload(
    summary = summary,
    selections = selections,
    comparisons = comparisons,
    winner_summary = winner_summary,
    method_registry = method_registry,
    config = config,
    prefix = prefix,
    generated_at = generated_at
  )
  parent <- dirname(output_path)
  if (nzchar(parent) && !dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(payload, output_path, dataframe = "rows", auto_unbox = TRUE, na = "null", pretty = TRUE)
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
