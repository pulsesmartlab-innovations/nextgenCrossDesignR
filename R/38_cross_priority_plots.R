ng_priority_plot_pair_key <- function(data, parent1_col, parent2_col) {
  paste(as.character(data[[parent1_col]]), as.character(data[[parent2_col]]), sep = "\r")
}

ng_priority_plot_tier_colors <- function(labels = c("highly_priority", "priority", "medium_priority", "low_priority")) {
  cols <- c(
    highly_priority = "#1b6b54",
    priority = "#4f8f7a",
    medium_priority = "#d6a23c",
    low_priority = "#b85c38"
  )
  labels <- as.character(labels)
  missing <- setdiff(labels, names(cols))
  if (length(missing)) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    cols <- c(cols, extra)
  }
  cols[labels]
}

ng_plot_priority_score_vs_kinship <- function(scored,
                                              selected = NULL,
                                              output_path = NULL,
                                              score_col = "multi_trait_score",
                                              kinship_col = "pair_kinship",
                                              tier_col = "priority_tier",
                                              parent1_col = "parent1",
                                              parent2_col = "parent2",
                                              title = "Selected priority tiers versus all candidate crosses",
                                              width = 8,
                                              height = 5,
                                              res = 150) {
  scored <- as.data.frame(scored, stringsAsFactors = FALSE)
  required <- c(parent1_col, parent2_col, score_col, kinship_col)
  miss <- setdiff(required, names(scored))
  if (length(miss)) ng_stop("scored missing columns: ", paste(miss, collapse = ", "))
  if (is.null(selected)) selected <- scored
  selected <- as.data.frame(selected, stringsAsFactors = FALSE)
  required_selected <- c(parent1_col, parent2_col, score_col, kinship_col, tier_col)
  miss_selected <- setdiff(required_selected, names(selected))
  if (length(miss_selected)) ng_stop("selected missing columns: ", paste(miss_selected, collapse = ", "))
  if (!nrow(scored)) ng_stop("scored must contain at least one candidate cross")
  if (!nrow(selected)) ng_stop("selected must contain at least one selected cross")

  owns_device <- !is.null(output_path)
  dev_no <- NULL
  if (isTRUE(owns_device)) {
    output_path <- as.character(output_path[[1L]])
    dev_no <- ng_plot_open_device(output_path, width = width, height = height, res = res)
    on.exit(ng_plot_close_device(dev_no), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  # Restore par only while a device is still open; after we close an owned device the
  # restore would otherwise spawn a stray default device (Rplots.pdf).
  on.exit(if (length(grDevices::dev.list())) graphics::par(old_par), add = TRUE)

  score <- suppressWarnings(as.numeric(scored[[score_col]]))
  kinship <- suppressWarnings(as.numeric(scored[[kinship_col]]))
  selected_score <- suppressWarnings(as.numeric(selected[[score_col]]))
  selected_kinship <- suppressWarnings(as.numeric(selected[[kinship_col]]))
  ok <- is.finite(score) & is.finite(kinship)
  selected_ok <- is.finite(selected_score) & is.finite(selected_kinship)
  if (!any(ok)) ng_stop("scored contains no finite score/kinship pairs")
  if (!any(selected_ok)) ng_stop("selected contains no finite score/kinship pairs")

  graphics::plot(
    kinship[ok], score[ok],
    pch = 16,
    cex = 0.42,
    col = grDevices::adjustcolor("#5f6b73", alpha.f = 0.22),
    xlab = "Pair kinship",
    ylab = "Multi-trait score",
    main = title
  )

  selected <- selected[selected_ok, , drop = FALSE]
  selected_score <- selected_score[selected_ok]
  selected_kinship <- selected_kinship[selected_ok]
  if ("priority_rank" %in% names(selected)) {
    ord <- order(suppressWarnings(as.integer(selected$priority_rank)), na.last = TRUE)
    selected <- selected[ord, , drop = FALSE]
    selected_score <- selected_score[ord]
    selected_kinship <- selected_kinship[ord]
  }
  tiers <- as.character(selected[[tier_col]])
  tier_levels <- if (is.factor(selected[[tier_col]])) levels(selected[[tier_col]]) else unique(tiers)
  tier_levels <- tier_levels[tier_levels %in% tiers]
  tier_cols <- ng_priority_plot_tier_colors(tier_levels)
  graphics::points(
    selected_kinship,
    selected_score,
    pch = 16,
    cex = 0.9,
    col = tier_cols[tiers]
  )
  graphics::legend(
    "topright",
    legend = c("All candidate crosses", gsub("_", " ", tier_levels)),
    col = c(grDevices::adjustcolor("#5f6b73", alpha.f = 0.35), tier_cols),
    pch = 16,
    pt.cex = c(0.65, rep(0.9, length(tier_levels))),
    bty = "n",
    cex = 0.82
  )
  if (isTRUE(owns_device)) {
    ng_plot_close_device(dev_no)
    invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
  } else {
    invisible(NULL)
  }
}
