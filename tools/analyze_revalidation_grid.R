# Analyze a v0.2.0 re-validation grid. Configurable via env vars:
#   NG_ANALYSIS_PREFIX   default "revalidation_v0_2_0_full"
#   NG_ANALYSIS_SIZES    default "20,30,40,50,60,70,80"
#   NG_ANALYSIS_CYCLE    default "max" (final cycle only); "all" pools all
#                                       cycles into a single multi-cycle mean
#   NG_ANALYSIS_RESULTS  default "../results"
# Emits a markdown report to stdout. The replicated v0.2.0 smoke grid can
# also be analyzed by setting NG_ANALYSIS_PREFIX=revalidation_v0_2_0_replicated
# and NG_ANALYSIS_SIZES=20,40,60,80.

env_chr <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || !nzchar(v)) default else v
}

prefix <- env_chr("NG_ANALYSIS_PREFIX", "revalidation_v0_2_0_full")
sizes  <- as.integer(strsplit(env_chr("NG_ANALYSIS_SIZES", "20,30,40,50,60,70,80"), ",", fixed = TRUE)[[1L]])
cycle_mode <- env_chr("NG_ANALYSIS_CYCLE", "max")
results_dir <- normalizePath(env_chr("NG_ANALYSIS_RESULTS", file.path("..", "results")))

read_size <- function(p) {
  f <- file.path(results_dir, sprintf("%s_%dp_metrics_by_rep.csv", prefix, p))
  if (!file.exists(f)) {
    warning("missing ", f); return(NULL)
  }
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$n_parents <- p
  d
}
raw <- do.call(rbind, lapply(sizes, read_size))
raw <- raw[!is.na(raw$top10_gv), , drop = FALSE]
cycles_avail <- sort(unique(raw$cycle))
nonzero_cycles <- cycles_avail[cycles_avail > 0]
target_cycle <- if (identical(cycle_mode, "all")) nonzero_cycles else max(nonzero_cycles)
raw <- raw[raw$cycle %in% target_cycle, , drop = FALSE]

agg <- do.call(rbind, lapply(sizes, function(p) {
  sub <- raw[raw$n_parents == p, , drop = FALSE]
  by_method <- split(sub$top10_gv, sub$method)
  out <- data.frame(
    n_parents = p,
    method    = names(by_method),
    n_obs     = vapply(by_method, length, integer(1L)),
    mean_top10 = vapply(by_method, mean, numeric(1L)),
    sd_top10  = vapply(by_method, function(x) if (length(x) > 1L) sd(x) else NA_real_, numeric(1L)),
    stringsAsFactors = FALSE
  )
  out[order(-out$mean_top10), , drop = FALSE]
}))

cat(sprintf("# Re-validation grid analysis — `%s`\n\n", prefix))
cat(sprintf("Cycle(s): %s. Sizes: %s. Results dir: `%s`.\n\n",
            paste(target_cycle, collapse = ","),
            paste(sizes, collapse = ","), results_dir))

cat("## Multi-rep realized top-10 GV (mean ± SD per parent size)\n")
for (p in sizes) {
  sub <- agg[agg$n_parents == p, , drop = FALSE]
  if (!nrow(sub)) next
  cat(sprintf("\n### p=%d (n=%d obs/method)\n\n", p, sub$n_obs[1L]))
  cat("| Method | mean top10 | SD | Rank |\n|---|---:|---:|---:|\n")
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("| `%s` | %.3f | %.3f | %d |\n",
                sub$method[i], sub$mean_top10[i],
                ifelse(is.finite(sub$sd_top10[i]), sub$sd_top10[i], 0),
                i))
  }
}

cat("\n## Frontier policy vs best non-frontier per parent size\n\n")
cat("| Parents | Frontier mean ± SD | Best non-frontier | Δ % | Welch p |\n")
cat("|---:|---|---|---:|---:|\n")
for (p in sizes) {
  rep_p <- raw[raw$n_parents == p, , drop = FALSE]
  frontier <- rep_p[rep_p$method == "ng_frontier_policy_ocs10_lps2", "top10_gv"]
  others   <- rep_p[rep_p$method != "ng_frontier_policy_ocs10_lps2", , drop = FALSE]
  if (!length(frontier) || !nrow(others)) next
  by_m <- split(others$top10_gv, others$method)
  means <- vapply(by_m, mean, numeric(1L))
  best_name <- names(means)[which.max(means)]
  best <- by_m[[best_name]]
  delta_pct <- 100 * (mean(best) - mean(frontier)) / mean(frontier)
  welch <- tryCatch(stats::t.test(best, frontier, var.equal = FALSE)$p.value,
                    error = function(e) NA_real_)
  cat(sprintf("| %d | %.3f ± %.3f | `%s` %.3f ± %.3f | %+0.2f%% | %.3g |\n",
              p, mean(frontier), if (length(frontier) > 1L) sd(frontier) else 0,
              best_name, mean(best), if (length(best) > 1L) sd(best) else 0,
              delta_pct, welch))
}

alphamate_methods <- grep("^alphamate_", unique(raw$method), value = TRUE)
if (length(alphamate_methods)) {
  cat("\n## Frontier policy vs real AlphaMate per parent size\n\n")
  cat("Positive Δ = frontier beats AlphaMate on average. The v0.0.x headline was 'frontier beats best AlphaMate in 7/7 bands'.\n\n")
  cat("| Parents | Frontier mean ± SD | AlphaMate method | AlphaMate mean ± SD | Δ % | Welch p |\n")
  cat("|---:|---|---|---|---:|---:|\n")
  for (p in sizes) {
    rep_p <- raw[raw$n_parents == p, , drop = FALSE]
    frontier <- rep_p[rep_p$method == "ng_frontier_policy_ocs10_lps2", "top10_gv"]
    if (!length(frontier)) next
    for (am in alphamate_methods) {
      am_vals <- rep_p[rep_p$method == am, "top10_gv"]
      if (!length(am_vals)) next
      delta_pct <- 100 * (mean(frontier) - mean(am_vals)) / mean(am_vals)
      welch <- tryCatch(stats::t.test(frontier, am_vals, var.equal = FALSE)$p.value,
                        error = function(e) NA_real_)
      cat(sprintf("| %d | %.3f ± %.3f | `%s` | %.3f ± %.3f | %+0.2f%% | %.3g |\n",
                  p,
                  mean(frontier), if (length(frontier) > 1L) sd(frontier) else 0,
                  am,
                  mean(am_vals), if (length(am_vals) > 1L) sd(am_vals) else 0,
                  delta_pct, welch))
    }
  }
}

cat("\n## Consistency: top-3 finishes across parent sizes\n\n")
methods <- unique(raw$method)
top3 <- vapply(methods, function(m) {
  sum(vapply(sizes, function(p) {
    sub <- agg[agg$n_parents == p, , drop = FALSE]
    if (!nrow(sub)) return(FALSE)
    m %in% sub$method[seq_len(min(3L, nrow(sub)))]
  }, logical(1L)))
}, integer(1L))
ord <- order(-top3)
cat("| Method | Top-3 finishes |\n|---|---:|\n")
for (i in ord) {
  cat(sprintf("| `%s` | %d / %d |\n", methods[i], top3[i], length(sizes)))
}
