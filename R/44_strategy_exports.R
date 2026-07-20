# Export artifacts for the decision-workbench frontend (a separate project). These emit
# schema-versioned JSON the UI can render directly: the gain-vs-diversity frontier with
# strategy presets, the progeny-inbreeding histogram, and a side-by-side run comparison.
# The backend does not render any UI; docs/frontend/MATING_STRATEGY_INTEGRATION.md maps
# each artifact/field to the control it drives.

# Build the gain-vs-diversity frontier payload: one point per lambda on the native sweep
# (annotated with its achieved diversity_emphasis and mean progeny inbreeding), plus the
# three strategy presets resolved to their operating points.
ng_frontier_export_payload <- function(scores,
                                       n_crosses,
                                       parent_kinship,
                                       gain_col = "usefulness_pmv_gebv",
                                       lambdas = 10 ^ seq(-2, 3, length.out = 10),
                                       strategies = c("high_gain", "balanced", "diversity"),
                                       generated_at = Sys.time(),
                                       ...) {
  sweep <- ng_pareto_mate_allocation(scores = scores, n_crosses = n_crosses,
                                     gain_col = gain_col, parent_kinship = parent_kinship,
                                     lambdas = lambdas, ...)
  fr <- sweep$frontier
  gain_range <- range(suppressWarnings(as.numeric(fr$mean_gain)), na.rm = TRUE, finite = TRUE)
  progeny_f <- vapply(sweep$plans, function(p) {
    v <- attr(p, "summary")$mean_progeny_inbreeding
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  points <- data.frame(
    lambda_group = fr$lambda_group,
    diversity_emphasis = vapply(seq_len(nrow(fr)),
      function(i) ng_frontier_achieved_emphasis(fr$mean_gain[[i]], gain_range), numeric(1)),
    mean_gain = fr$mean_gain,
    group_coancestry = fr$group_coancestry,
    mean_progeny_inbreeding = progeny_f,
    unique_parents = fr$unique_parents,
    max_parent_use = fr$max_parent_use,
    stringsAsFactors = FALSE
  )
  strat <- lapply(strategies, function(st) {
    plan <- ng_select_by_strategy(scores = scores, n_crosses = n_crosses,
                                  parent_kinship = parent_kinship, gain_col = gain_col,
                                  strategy = st, lambdas = lambdas, ...)
    s <- attr(plan, "summary")
    list(strategy = st,
         diversity_emphasis = s$diversity_emphasis,
         achieved_emphasis = s$achieved_emphasis,
         selected_lambda_group = s$selected_lambda_group,
         mean_gain = s$mean_gain,
         group_coancestry = s$group_coancestry,
         mean_progeny_inbreeding = s$mean_progeny_inbreeding,
         unique_parents = s$unique_parents)
  })
  names(strat) <- strategies
  list(
    schema = "ng_mating_frontier.v1",
    generated_at = format(generated_at, tz = "UTC"),
    n_crosses = n_crosses,
    gain_col = gain_col,
    emphasis_scale = "0 = maximum genetic gain; 100 = maximum diversity",
    points = points,
    strategies = strat
  )
}

ng_write_frontier_json <- function(scores, n_crosses, parent_kinship, output_path = NULL,
                                   gain_col = "usefulness_pmv_gebv", generated_at = Sys.time(), ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to write frontier JSON")
  }
  if (is.null(output_path) || !nzchar(output_path)) {
    output_path <- file.path("results", "mating_frontier.json")
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  payload <- ng_frontier_export_payload(scores, n_crosses, parent_kinship, gain_col = gain_col,
                                        generated_at = generated_at, ...)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null",
                       dataframe = "rows")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}

ng_progeny_inbreeding_histogram_payload <- function(x, breaks = 20L, generated_at = Sys.time()) {
  h <- ng_progeny_inbreeding_histogram(x, breaks = breaks)
  list(
    schema = "ng_progeny_inbreeding_histogram.v1",
    generated_at = format(generated_at, tz = "UTC"),
    n = attr(h, "n"),
    mean_progeny_inbreeding = attr(h, "mean_progeny_inbreeding"),
    max_progeny_inbreeding = attr(h, "max_progeny_inbreeding"),
    bins = h
  )
}

ng_write_progeny_inbreeding_histogram_json <- function(x, output_path = NULL, breaks = 20L,
                                                       generated_at = Sys.time()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to write progeny-inbreeding histogram JSON")
  }
  if (is.null(output_path) || !nzchar(output_path)) {
    output_path <- file.path("results", "progeny_inbreeding_histogram.json")
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  payload <- ng_progeny_inbreeding_histogram_payload(x, breaks = breaks, generated_at = generated_at)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null",
                       dataframe = "rows")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}

# Side-by-side comparison of several optimizer runs (the "compare runs" view). `runs` is a
# named list of plans produced by ng_optimize_mating_plan (each carrying a summary attr).
ng_run_comparison_payload <- function(runs, generated_at = Sys.time()) {
  if (!length(runs)) ng_stop("runs is empty")
  labels <- names(runs)
  if (is.null(labels)) labels <- paste0("run", seq_along(runs))
  fields <- c("n_crosses", "mean_gain", "total_gain", "group_coancestry",
              "mean_progeny_inbreeding", "max_progeny_inbreeding", "mean_pair_kinship",
              "unique_parents", "max_parent_use", "diversity_emphasis", "achieved_emphasis",
              "total_cost", "n_committed", "constraints_reduced_plan")
  rows <- lapply(seq_along(runs), function(i) {
    s <- attr(runs[[i]], "summary")
    if (is.null(s)) s <- list()
    row <- list(run = labels[[i]])
    for (f in fields) row[[f]] <- if (is.null(s[[f]])) NA else s[[f]]
    row
  })
  list(
    schema = "ng_run_comparison.v1",
    generated_at = format(generated_at, tz = "UTC"),
    fields = fields,
    runs = rows
  )
}

ng_write_run_comparison_json <- function(runs, output_path = NULL, generated_at = Sys.time()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    ng_stop("jsonlite is required to write run-comparison JSON")
  }
  if (is.null(output_path) || !nzchar(output_path)) {
    output_path <- file.path("results", "run_comparison.json")
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  payload <- ng_run_comparison_payload(runs, generated_at = generated_at)
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE, pretty = TRUE, na = "null",
                       dataframe = "rows")
  normalizePath(output_path, winslash = "/", mustWork = TRUE)
}
