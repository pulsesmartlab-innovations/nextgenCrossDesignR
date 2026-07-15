# Render publication-grade figures for the nextgenCrossDesign package.
#
# Usage: Rscript tools/render_publication_figures.R [fig_id...]
#   No args  -> renders all six figures (F1..F6)
#   "f1 f4"  -> renders only F1 and F4 (case-insensitive)
#
# Outputs (under docs/figures/, both formats per figure):
#   F<n>_<slug>.png   300 DPI, ~7"x5" raster for README / release notes
#   F<n>_<slug>.pdf   vector copy for publication / poster
#
# Data sources (read from <repo>/../results/ for head-to-head + multi-trait,
# and <repo>/results/ for poly4x + cpp_speedup_bench). If an artifact is
# missing the loader invokes the matching runner once.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ---- Repo + path resolution -----------------------------------------------

find_repo_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(start, dirname(start), dirname(dirname(start)),
                         file.path(start, "nextgen_cross_design"),
                         file.path(dirname(start), "nextgen_cross_design")))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "R", "load.R")) &&
        file.exists(file.path(cand, "DESCRIPTION"))) {
      return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Could not locate nextgen_cross_design package root.")
}

repo_root <- find_repo_root()
external_results <- normalizePath(file.path(repo_root, "..", "results"),
                                  winslash = "/", mustWork = FALSE)
internal_results <- normalizePath(file.path(repo_root, "results"),
                                  winslash = "/", mustWork = FALSE)
figures_dir <- normalizePath(file.path(repo_root, "docs", "figures"),
                             winslash = "/", mustWork = FALSE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

# ---- Shared theme + palette (mirrors R/23_head_to_head_visual_report) -----

ng_pub_palette <- c(
  candidate = "#2b7d6b",   # frontier-policy / nextgen methods
  baseline  = "#b85c38",   # AlphaMate / SimpleMating / external
  neutral   = "#2f5f98",   # greedy / random
  warn      = "#b7791f"    # style-proxy or low-evidence
)

ng_publication_theme <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2,
                                margin = margin(b = 6)),
      plot.subtitle = element_text(colour = "#4b5866",
                                   margin = margin(b = 10), size = base_size - 1),
      plot.caption = element_text(colour = "#64707d", hjust = 0,
                                  size = base_size - 2, margin = margin(t = 8)),
      axis.title = element_text(colour = "#1f2933"),
      axis.text = element_text(colour = "#4b5866"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e4e9ef"),
      strip.text = element_text(face = "bold", colour = "#1f2933"),
      strip.background = element_rect(fill = "#f5f7fa", colour = NA),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 1)
    )
}

# Map a method label to the candidate/baseline/neutral role used by the
# HTML report so the publication figures color-match it.
ng_method_role <- function(method, role = NULL) {
  if (!is.null(role)) {
    r <- tolower(trimws(as.character(role)))
    out <- ifelse(r == "candidate", "candidate",
                  ifelse(r == "baseline", "baseline", "neutral"))
    return(out)
  }
  m <- tolower(as.character(method))
  out <- rep("neutral", length(m))
  # Candidate first, then baseline — baseline regex wins on overlap because
  # the next assignment overrides. `style` / `style_proxy` suffixes are
  # baseline regardless of which external tool they emulate.
  out[grepl("nextgen|frontier|^ng_", m)] <- "candidate"
  out[grepl("alphamate|simplemating|simplemate|popvar|external|style", m)] <- "baseline"
  out
}

# ---- Loaders --------------------------------------------------------------

read_csv_required <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path))
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

invoke_runner <- function(script, env = list(), label) {
  message(sprintf("[loader] %s missing, invoking %s ...", label, script))
  # Cross-platform env passing: set in parent process before Rscript fork.
  if (length(env)) {
    prior <- Sys.getenv(names(env), names = TRUE, unset = NA)
    do.call(Sys.setenv, env)
    on.exit({
      for (nm in names(prior)) {
        if (is.na(prior[[nm]])) Sys.unsetenv(nm) else Sys.setenv(setNames(list(prior[[nm]]), nm))
      }
    }, add = TRUE)
  }
  status <- system2("Rscript", args = shQuote(script), wait = TRUE)
  if (status != 0L) stop(sprintf("Runner failed (exit %d) for %s", status, script))
}

load_head_to_head <- function(prefix = "head_to_head_pub",
                              results_dir = external_results) {
  summary_path <- file.path(results_dir, paste0(prefix, "_summary.csv"))
  if (!file.exists(summary_path)) {
    invoke_runner(
      file.path(repo_root, "tools", "run_head_to_head_benchmark.R"),
      env = list(
        NG_HEAD_TO_HEAD_PREFIX = prefix,
        NG_HEAD_TO_HEAD_PARENT_SIZES = "20,40,60",
        NG_HEAD_TO_HEAD_REPS = "3",
        NG_HEAD_TO_HEAD_CROSSES = "6",
        NG_HEAD_TO_HEAD_REALIZED_PROGENY = "20",
        NG_HEAD_TO_HEAD_SCENARIOS = "compact_selfing,cassava_diploid"
      ),
      label = "head-to-head benchmark outputs"
    )
  }
  list(
    summary   = read_csv_required(summary_path, "head-to-head summary"),
    selections = read_csv_required(
      file.path(results_dir, paste0(prefix, "_selections.csv")), "selections"),
    method_registry = tryCatch(
      read_csv_required(file.path(results_dir, paste0(prefix, "_method_registry.csv")),
                        "method_registry"),
      error = function(e) data.frame())
  )
}

load_cpp_bench <- function() {
  path <- file.path(external_results, "cpp_speedup_bench.csv")
  if (!file.exists(path)) {
    invoke_runner(file.path(repo_root, "tools", "cpp_speedup_bench.R"),
                  env = list(NG_USE_CPP = "1"),
                  label = "C++ speedup bench CSV")
  }
  read_csv_required(path, "cpp speedup bench")
}

load_multitrait <- function(prefix = "multitrait_grid_pub",
                            results_dir = external_results) {
  summary_path <- file.path(results_dir, paste0(prefix, "_summary.csv"))
  if (!file.exists(summary_path)) {
    invoke_runner(
      file.path(repo_root, "tools", "run_multitrait_validation_grid.R"),
      env = list(
        NG_MULTITRAIT_GRID_PREFIX = prefix,
        NG_MULTITRAIT_GRID_PARENT_SIZES = "20,40,60",
        NG_MULTITRAIT_GRID_REPS = "3",
        NG_MULTITRAIT_GRID_CROSSES = "6"
      ),
      label = "multitrait grid outputs"
    )
  }
  read_csv_required(summary_path, "multitrait grid summary")
}

load_poly4x <- function(results_dir = internal_results,
                        pattern = "^poly4x_controlled_runner_smoke_.+_metrics\\.csv$") {
  files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No poly4x metrics CSVs found in ", results_dir)
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE,
                        check.names = FALSE))
}

# ---- Helpers used by builders --------------------------------------------

# 95% bootstrap CI for the mean of x (BCa-free, percentile method).
bootstrap_mean_ci <- function(x, R = 2000L, seed = 1L) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(c(mean = if (length(x)) x else NA_real_,
                               lower = NA_real_, upper = NA_real_))
  set.seed(seed)
  draws <- replicate(R, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x),
    lower = unname(quantile(draws, 0.025, na.rm = TRUE)),
    upper = unname(quantile(draws, 0.975, na.rm = TRUE)))
}

summarize_by <- function(data, group_cols, value_col) {
  parts <- split(data, lapply(group_cols, function(g) data[[g]]), drop = TRUE)
  out <- do.call(rbind, lapply(parts, function(part) {
    ci <- bootstrap_mean_ci(part[[value_col]])
    row <- part[1L, group_cols, drop = FALSE]
    row$mean  <- ci[["mean"]]
    row$lower <- ci[["lower"]]
    row$upper <- ci[["upper"]]
    row$n     <- nrow(part)
    row
  }))
  rownames(out) <- NULL
  out
}

# ---- Builders -------------------------------------------------------------

# F1. Gain-diversity Pareto frontier
build_f1 <- function(payload) {
  summary <- payload$summary
  agg <- summarize_by(summary, c("method", "n_parents"), "mean_realized_index")
  agg$coancestry <- vapply(seq_len(nrow(agg)), function(i) {
    keep <- summary$method == agg$method[i] & summary$n_parents == agg$n_parents[i]
    mean(summary$group_coancestry[keep], na.rm = TRUE)
  }, numeric(1))
  agg$role <- ng_method_role(agg$method,
                             role = vapply(seq_len(nrow(agg)), function(i) {
                               keep <- summary$method == agg$method[i]
                               if (!"benchmark_role" %in% names(summary) || !any(keep)) return(NA_character_)
                               as.character(summary$benchmark_role[keep][1L])
                             }, character(1)))

  # Compute upper-left Pareto frontier (max y, min x).
  ord <- order(agg$coancestry, -agg$mean)
  front_idx <- integer()
  best_y <- -Inf
  for (i in ord) {
    if (agg$mean[i] > best_y) { front_idx <- c(front_idx, i); best_y <- agg$mean[i] }
  }
  frontier <- agg[front_idx, , drop = FALSE]
  frontier <- frontier[order(frontier$coancestry), , drop = FALSE]

  ggplot(agg, aes(x = coancestry, y = mean, color = role,
                  shape = factor(n_parents))) +
    geom_line(data = frontier, aes(x = coancestry, y = mean),
              inherit.aes = FALSE, color = "#1f2933",
              linetype = "dashed", linewidth = 0.4) +
    geom_point(size = 3, stroke = 0.7) +
    ggrepel_or_text(agg) +
    scale_color_manual(values = ng_pub_palette, name = "Role",
                       labels = c(candidate = "Candidate (NextGen)",
                                  baseline  = "Baseline (external)",
                                  neutral   = "Greedy / random")) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8),
                       name = "Parent-pool size") +
    labs(
      title = "Gain-diversity Pareto frontier",
      subtitle = "Higher selection index and lower group coancestry are preferred (upper-left).",
      x = "Group coancestry (lower is better)",
      y = "Mean realized selection index (higher is better)",
      caption = sprintf("Source: head-to-head benchmark (%d methods x %d parent sizes).",
                        length(unique(agg$method)), length(unique(agg$n_parents)))
    ) +
    ng_publication_theme()
}

# Light-weight label placement that doesn't require ggrepel.
ggrepel_or_text <- function(agg) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    return(ggrepel::geom_text_repel(aes(label = method), size = 2.7,
                                    max.overlaps = 50, segment.alpha = 0.4))
  }
  geom_text(aes(label = method), size = 2.7, hjust = -0.15, vjust = -0.4)
}

# F2. Per-band ranking with bootstrap CI
build_f2 <- function(payload) {
  summary <- payload$summary
  agg <- summarize_by(summary, c("method", "n_parents"), "mean_realized_index")
  agg$role <- ng_method_role(agg$method)
  # Order methods within each band by mean (highest at top).
  agg <- agg[order(agg$n_parents, agg$mean), , drop = FALSE]
  agg$method <- factor(agg$method, levels = unique(agg$method))

  ggplot(agg, aes(x = mean, y = method, fill = role)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25,
                  color = "#1f2933", linewidth = 0.35,
                  orientation = "y") +
    facet_wrap(~ paste("n =", n_parents), scales = "free_y", ncol = 1) +
    scale_fill_manual(values = ng_pub_palette, name = "Role") +
    labs(
      title = "Per-band method ranking",
      subtitle = "Mean realized selection index with 95% bootstrap CI across replicates.",
      x = "Mean realized selection index",
      y = NULL,
      caption = "Wider CIs == low evidence. Compare against VALIDATED_STATE.md before claiming a band-level win."
    ) +
    ng_publication_theme() +
    theme(panel.grid.major.y = element_blank())
}

# F3. Per-method distribution across reps
build_f3 <- function(payload) {
  summary <- payload$summary
  summary$role <- ng_method_role(summary$method)
  # Drop methods seen in only one rep within all bands.
  ggplot(summary, aes(x = mean_realized_index, y = method, fill = role)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.85, color = "#1f2933",
                 linewidth = 0.3) +
    geom_jitter(height = 0.18, size = 0.5, color = "#1f2933", alpha = 0.45) +
    facet_wrap(~ paste("n =", n_parents), scales = "free_y", ncol = 1) +
    scale_fill_manual(values = ng_pub_palette, name = "Role") +
    labs(
      title = "Realized index distribution across replicates",
      subtitle = "Boxes are quartiles; jitter shows per-rep realizations. Wider boxes reveal noise hidden by F2 means.",
      x = "Realized selection index per rep",
      y = NULL,
      caption = sprintf("n = %d replicates x %d scenarios x %d methods.",
                        length(unique(summary$rep)),
                        length(unique(summary$scenario)),
                        length(unique(summary$method)))
    ) +
    ng_publication_theme() +
    theme(panel.grid.major.y = element_blank())
}

# F4. C++ speedup bar chart
build_f4 <- function(bench) {
  long <- data.frame(
    kernel  = rep(bench$kernel, 2L),
    backend = rep(c("R", "C++"), each = nrow(bench)),
    seconds = c(bench$r_seconds, bench$cpp_seconds),
    stringsAsFactors = FALSE
  )
  long$backend <- factor(long$backend, levels = c("R", "C++"))
  bench$label <- sprintf("%.0fx", bench$speedup_x)
  bench$ypos <- pmax(bench$r_seconds, bench$cpp_seconds)

  p <- ggplot(long, aes(x = kernel, y = seconds, fill = backend)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(data = bench, aes(x = kernel, y = ypos, label = label),
              inherit.aes = FALSE, vjust = -0.6, fontface = "bold",
              size = 3.5, color = "#1f2933") +
    scale_y_log10(labels = scales::label_number(suffix = "s"),
                  expand = expansion(mult = c(0.02, 0.2))) +
    scale_fill_manual(values = c(R = ng_pub_palette[["baseline"]],
                                  `C++` = ng_pub_palette[["candidate"]]),
                       name = "Backend") +
    labs(
      title = "v0.3.0 C++ kernel speedups",
      subtitle = "Wall time per call on n=60 parents, m=1500 markers. Speedup multiplier annotated.",
      x = NULL, y = "Wall time (log scale)",
      caption = "Skipped after merit check: BLAS-backed full-posterior (current C++ beats BLAS via column-skip + fused loop); MCMC Gibbs (~20s typical, modest 2-3x port gain, defer)."
    ) +
    ng_publication_theme() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  p
}

# F5. Multi-trait index method comparison
build_f5 <- function(grid) {
  metric_cols <- c(index   = "mean_realized_index",
                   yield   = "mean_realized_yield",
                   disease = "mean_realized_disease",
                   quality = "mean_realized_quality")
  metric_cols <- metric_cols[metric_cols %in% names(grid)]
  if (!length(metric_cols)) stop("No expected multi-trait metric columns found")
  long <- do.call(rbind, lapply(names(metric_cols), function(short) {
    data.frame(method = grid$method,
               family = grid$method_family,
               n_parents = grid$n_parents,
               rep = grid$rep,
               metric = short,
               value = grid[[metric_cols[[short]]]],
               stringsAsFactors = FALSE)
  }))
  long$metric <- factor(long$metric, levels = c("index", "yield", "disease", "quality"),
                        labels = c("Selection index (max)",
                                   "Yield (max)",
                                   "Disease (min)",
                                   "Quality (max)"))
  agg <- do.call(rbind, by(long, list(long$method, long$metric, long$n_parents),
                           function(part) {
    ci <- bootstrap_mean_ci(part$value)
    data.frame(method = part$method[1L], metric = part$metric[1L],
               n_parents = part$n_parents[1L],
               family = part$family[1L],
               mean = ci[["mean"]], lower = ci[["lower"]], upper = ci[["upper"]],
               stringsAsFactors = FALSE)
  }, simplify = FALSE))

  ggplot(agg, aes(x = factor(n_parents), y = mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.8), width = 0.2,
                  color = "#1f2933", linewidth = 0.3) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Set2", name = "Method") +
    labs(
      title = "Multi-trait scoring method comparison",
      subtitle = "Realized response per trait and per scoring method. Error bars: 95% bootstrap CI across reps.",
      x = "Parent-pool size",
      y = "Mean realized value",
      caption = "Methods: auto = rank-threshold, weighted, economic_index = Smith-Hazel, desired_gain = Pesek-Baker, threshold = soft-penalty."
    ) +
    ng_publication_theme()
}

# F6. Polyploid (autotetraploid) scenario
build_f6 <- function(metrics) {
  metrics$method <- gsub("^ng_poly4x_", "", metrics$method)
  agg <- summarize_by(metrics, c("method"), "mean_gv")
  agg <- agg[order(agg$mean), , drop = FALSE]
  metrics$method <- factor(metrics$method, levels = agg$method)

  p_box <- ggplot(metrics, aes(x = mean_gv, y = method)) +
    geom_boxplot(fill = ng_pub_palette[["candidate"]], alpha = 0.85,
                 outlier.size = 0.6, color = "#1f2933", linewidth = 0.3) +
    geom_jitter(height = 0.18, size = 0.5, color = "#1f2933", alpha = 0.4) +
    labs(
      title = "Autotetraploid realized mean GV per method",
      subtitle = "Boxes: quartiles across replicates / cycles / configs. Higher is better.",
      x = "Mean realized genetic value", y = NULL
    ) +
    ng_publication_theme() +
    theme(panel.grid.major.y = element_blank())

  p_scatter <- ggplot(metrics, aes(x = var_gv, y = top10_gv, color = method)) +
    geom_point(size = 2.2, alpha = 0.8) +
    scale_color_brewer(palette = "Dark2", name = "Method") +
    labs(
      title = "Diversity x top-10 trade-off",
      subtitle = "Within-family variance vs top-10 GV - the polyploid frontier shape.",
      x = "Within-family variance of GV",
      y = "Top-10 GV"
    ) +
    ng_publication_theme()

  p_box / p_scatter +
    plot_annotation(
      caption = sprintf("Source: %d data points across %d configs of the poly4x controlled runner.",
                        nrow(metrics), length(unique(metrics$config)))
    )
}

# ---- Registry + dispatch --------------------------------------------------

figure_registry <- list(
  f1 = list(slug = "gain_diversity_pareto",
            builder = build_f1, loader = function() load_head_to_head(),
            width = 8.5, height = 6),
  f2 = list(slug = "per_band_ranking",
            builder = build_f2, loader = function() load_head_to_head(),
            width = 7.5, height = 8),
  f3 = list(slug = "per_method_distribution",
            builder = build_f3, loader = function() load_head_to_head(),
            width = 7.5, height = 8),
  f4 = list(slug = "cpp_speedup",
            builder = build_f4, loader = load_cpp_bench,
            width = 8, height = 5),
  f5 = list(slug = "multitrait_method_comparison",
            builder = build_f5, loader = load_multitrait,
            width = 8.5, height = 6.5),
  f6 = list(slug = "polyploid_autotetraploid",
            builder = build_f6, loader = load_poly4x,
            width = 8, height = 8)
)

render_one <- function(fig_id) {
  entry <- figure_registry[[fig_id]]
  if (is.null(entry)) stop(sprintf("Unknown figure id: %s", fig_id))
  message(sprintf("[%s] loading data ...", fig_id))
  payload <- entry$loader()
  message(sprintf("[%s] building plot ...", fig_id))
  plot_obj <- entry$builder(payload)
  base <- file.path(figures_dir,
                    sprintf("F%s_%s", toupper(sub("^f", "", fig_id)), entry$slug))
  ggsave(paste0(base, ".png"), plot_obj, width = entry$width, height = entry$height,
         dpi = 300, units = "in")
  ggsave(paste0(base, ".pdf"), plot_obj, width = entry$width, height = entry$height,
         units = "in", device = grDevices::cairo_pdf)
  message(sprintf("[%s] wrote %s.{png,pdf}", fig_id, base))
  invisible(base)
}

# ---- CLI -----------------------------------------------------------------

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  selected <- if (!length(args)) names(figure_registry) else tolower(args)
  unknown <- setdiff(selected, names(figure_registry))
  if (length(unknown)) stop(sprintf("Unknown figure id(s): %s",
                                    paste(unknown, collapse = ", ")))
  cat(sprintf("Rendering %d figure(s) -> %s\n",
              length(selected), figures_dir))
  for (id in selected) render_one(id)
  cat("done.\n")
}

if (sys.nframe() == 0L) main()
