local({ .h <- file.path("tools", "ng_project_libpath.R"); if (file.exists(.h)) { source(.h); ng_prepend_project_lib(".Rlib") } else .libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths())) })

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) return(default)
  value
}

bind_rows_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    for (m in missing) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  do.call(rbind, x)
}

prefixes <- strsplit(env_chr(
  "NG_COMPARE_GRID_PREFIXES",
  "balanced_external_exact_20_40_3rep,balanced_external_shortlist_50_80_3rep"
), ",", fixed = TRUE)[[1]]
prefixes <- trimws(prefixes[nzchar(prefixes)])
tier_names <- strsplit(env_chr("NG_COMPARE_GRID_TIERS", ""), ",", fixed = TRUE)[[1]]
tier_names <- trimws(tier_names[nzchar(tier_names)])
if (length(tier_names) != length(prefixes)) {
  tier_names <- prefixes
}
out_prefix <- env_chr("NG_COMPARE_OUTPUT_PREFIX", "balanced_external_parent_grid")

overall <- list()
for (i in seq_along(prefixes)) {
  path <- file.path("results", paste0(prefixes[[i]], "_overall_avg.csv"))
  if (!file.exists(path)) next
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$tier <- tier_names[[i]]
  overall[[prefixes[[i]]]] <- d
}
overall <- bind_rows_fill(overall)
if (!nrow(overall)) stop("No grid overall summary files found.", call. = FALSE)

overall$group <- ifelse(
  grepl("^(popvar_|simple_)", overall$method),
  "external",
  ifelse(grepl("^var_simple", overall$method), "var_simple", "local")
)

metrics <- c("mean_gv", "top10_gv", "max_gv")
rows <- list()
for (n in sort(unique(overall$n_parents))) {
  dn <- overall[overall$n_parents == n, , drop = FALSE]
  for (metric in metrics) {
    for (group in c("local", "external", "var_simple")) {
      dg <- dn[dn$group == group, , drop = FALSE]
      if (!nrow(dg)) next
      idx <- which.max(dg[[metric]])
      rows[[length(rows) + 1L]] <- data.frame(
        n_parents = n,
        tier = dg$tier[[idx]],
        metric = metric,
        group = group,
        method = dg$method[[idx]],
        value = dg[[metric]][[idx]],
        top_crosses = dg$top_crosses[[idx]],
        effect_training_n = dg$effect_training_n[[idx]],
        stringsAsFactors = FALSE
      )
    }
  }
}
group_best <- bind_rows_fill(rows)

wide <- reshape(group_best, idvar = c("n_parents", "tier", "metric"),
                timevar = "group", direction = "wide")
wide$local_minus_external <- wide$value.local - wide$value.external
wide$local_minus_var_simple <- wide$value.local - wide$value.var_simple

write.csv(overall, file.path("results", paste0(out_prefix, "_overall_avg_combined.csv")), row.names = FALSE)
write.csv(group_best, file.path("results", paste0(out_prefix, "_group_best_long.csv")), row.names = FALSE)
write.csv(wide, file.path("results", paste0(out_prefix, "_group_comparison.csv")), row.names = FALSE)

print(wide[order(wide$n_parents, wide$metric), ], row.names = FALSE)
