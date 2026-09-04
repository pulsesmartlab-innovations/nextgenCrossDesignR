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


# Look up the ok/P(beat check) pair that a scored table carries for its ACTIVE check(s), without
# requiring the caller to pass trait_check_reference through: the joint columns
# (checks_all_ok / p_beat_all_checks) win when present (multi-trait), and a single unambiguous
# per-trait pair (<trait>_check_ok / <trait>_p_beat_check) is used otherwise. Returns NULL when
# neither is resolvable so callers can fall back to unstyled markers -- this is what keeps a
# no-checks run's plot byte-for-byte the same as before this feature existed.
ng_priority_plot_check_style <- function(df) {
  ok <- df[["checks_all_ok"]]
  p <- df[["p_beat_all_checks"]]
  if (is.null(ok) || is.null(p)) {
    ok_cols <- grep("_check_ok$", names(df), value = TRUE)
    p_cols <- grep("_p_beat_check$", names(df), value = TRUE)
    if (length(ok_cols) != 1L || length(p_cols) != 1L) return(NULL)
    ok <- df[[ok_cols]]
    p <- df[[p_cols]]
  }
  list(ok = as.logical(ok), p = suppressWarnings(as.numeric(p)))
}

# Per-point alpha blend. grDevices::adjustcolor()'s alpha.f is a SCALAR -- it builds one
# 4x4 transform matrix via diag(), so passing a vector errors ("d == c(4L, 4L) are not all
# TRUE"). Opacity-by-P(beat check) needs one alpha per point, so blend manually via
# col2rgb()/rgb() instead, which are properly vectorized over both col and alpha.
ng_priority_plot_alpha_col <- function(col, alpha) {
  alpha <- pmax(0, pmin(1, alpha))
  rgb_mat <- grDevices::col2rgb(col)
  grDevices::rgb(rgb_mat[1L, ] / 255, rgb_mat[2L, ] / 255, rgb_mat[3L, ] / 255, alpha = alpha)
}

# Apply the check-reference marker styling: a point on the wrong side of its check is
# de-emphasised (flagged grey) rather than recolored to look like a normal tier point, and
# P(beat check) is encoded as opacity so a high-variance cross sitting below the line -- one
# that can still throw a superior tail -- stays visibly distinct from one that genuinely cannot.
# `style` is the (possibly NULL) result of ng_priority_plot_check_style(); base_col is the
# color each point would use with no check active.
ng_priority_plot_apply_check_style <- function(base_col, style) {
  if (is.null(style)) return(base_col)
  col <- ifelse(style$ok %in% FALSE, "#9AA0A6", base_col)
  alpha <- ifelse(is.finite(style$p), pmax(0.15, pmin(1, style$p)), 1)
  ng_priority_plot_alpha_col(col, alpha)
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
                                              res = 150,
                                              check_line = NULL,
                                              check_label = NULL) {
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
  if (!is.null(check_line) && is.finite(check_line)) {
    graphics::abline(h = check_line, lty = 2, lwd = 2, col = "#B00020")
    graphics::mtext(sprintf("check: %s", if (is.null(check_label)) "reference" else check_label),
                    side = 4, at = check_line, las = 1, cex = 0.7, col = "#B00020")
  }

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
  point_col <- ng_priority_plot_apply_check_style(tier_cols[tiers],
                                                  ng_priority_plot_check_style(selected))
  graphics::points(
    selected_kinship,
    selected_score,
    pch = 16,
    cex = 0.9,
    col = point_col
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

# The y value at which to draw the check reference line. Single-trait runs plot the trait mean
# itself, so the check's value is the line. A LINEAR index (weighted / economic) is a fixed
# combination of trait values, so the check's index value is exact. A RANK-based index is a
# function of the candidate distribution, and a check has no rank because it is not a cross --
# there is no honest line, so we return NA and the caller says so rather than drawing a number
# that looks authoritative and is not.
#
# THE AXIS IS NEVER RAW TRAIT UNITS. ng_score_breeder_objective() -> ng_add_multitrait_score()
# (R/19_multi_trait_selection.R) z-normalises multi_trait_score for EVERY method, including a
# single-trait run (one trait still goes through the same standardization loop). A check has no
# row in that computation and therefore no z of its own, so the only honest way to place it on
# this axis is to push its raw value through the SAME affine transform the candidates went
# through: oriented <- sign * tau; z <- (oriented - center) / scale. ng_add_multitrait_score()
# retains those per-trait parameters in multi_trait_meta$value_centers / $value_scales /
# $value_signs (named by trait) for exactly this purpose. Without them (an older result, or a
# multi_trait_meta the caller built by hand without these fields) there is no honest line, so we
# return NA -- the same discipline already applied to the rank-based-index case, now applied to
# every case, because drawing a raw-units number on a z-normalised axis is not "conservative
# single-trait math", it is simply wrong (it was invisible off the plotted range in real runs).
ng_check_line_value <- function(trait_check_reference, multi_trait_meta = NULL, trait = NULL) {
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  if (!nrow(spec)) return(NA_real_)
  val <- function(tr) {
    ck <- spec$check[[match(tr, spec$trait)]]
    suppressWarnings(as.numeric(trait_check_reference$values[[tr]][[ck]]))
  }
  centers <- multi_trait_meta$value_centers
  scales <- multi_trait_meta$value_scales
  signs <- multi_trait_meta$value_signs
  has_transform <- !is.null(centers) && !is.null(scales) && !is.null(signs)
  # The check's raw value, mapped onto the SAME standardized axis the candidates were scored on.
  z_of <- function(tr) {
    tau <- val(tr)
    if (!length(tau) || !is.finite(tau)) return(NA_real_)
    if (!has_transform || !(tr %in% names(centers)) || !(tr %in% names(scales)) ||
        !(tr %in% names(signs))) {
      return(NA_real_)
    }
    sign_i <- suppressWarnings(as.numeric(signs[[tr]]))
    center_i <- suppressWarnings(as.numeric(centers[[tr]]))
    scale_i <- suppressWarnings(as.numeric(scales[[tr]]))
    if (!is.finite(sign_i) || !is.finite(center_i) || !is.finite(scale_i) || scale_i <= 0) {
      return(NA_real_)
    }
    oriented <- if (sign_i > 0) tau else -tau
    (oriented - center_i) / scale_i
  }
  if (!is.null(trait) || nrow(spec) == 1L) {
    tr <- if (is.null(trait)) spec$trait[[1L]] else trait
    if (!(tr %in% spec$trait)) return(NA_real_)
    return(z_of(tr))
  }
  method <- as.character(if (is.null(multi_trait_meta$method)) "" else multi_trait_meta$method)
  if (!(method %in% c("weighted", "economic_index", "desired_gain"))) return(NA_real_)
  w <- multi_trait_meta$weights
  if (is.null(w) || !length(w)) return(NA_real_)
  tr <- intersect(spec$trait, names(w))
  if (!length(tr)) return(NA_real_)
  v <- vapply(tr, z_of, numeric(1))
  if (any(!is.finite(v))) return(NA_real_)
  sum(v * as.numeric(w[tr]))
}

# Label for the check reference line. `active$check[[1L]]` alone is only correct when every
# active trait shares the same check genotype -- the common case, but not the only one a linear
# index can legitimately combine. When the traits contributing to the line use DIFFERENT checks,
# naming just the first would misrepresent the line as a single check's value rather than the
# blend it actually is.
ng_check_line_label <- function(active) {
  active <- as.data.frame(active, stringsAsFactors = FALSE)
  cks <- unique(as.character(active$check))
  if (length(cks) <= 1L) return(if (length(cks)) cks[[1L]] else NA_character_)
  paste0("blended (", paste(cks, collapse = " + "), ")")
}

# Per-trait check panels: one facet per trait that has a check, y = that trait's mid-parent
# mean, x = pair kinship, with the trait's own check line. This is where multi-trait checks
# live, because a single index axis cannot carry several check lines on different scales.
# Each panel also encodes P(beat check) as marker opacity when that trait's <key>_p_beat_check
# column is present, for the same reason the main scatter does: a point below the line with a
# high P(beat check) can still throw a superior progeny and must stay visually distinct from
# one that genuinely cannot.
ng_plot_check_panels <- function(scored, trait_check_reference, output_path = NULL,
                                 kinship_col = "pair_kinship", width = 10, height = 4,
                                 res = 150) {
  if (is.null(trait_check_reference)) return(invisible(NULL))
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  keep <- paste0(key, "_mean") %in% names(scored)
  traits <- spec$trait[keep]; key <- key[keep]
  if (!length(traits)) return(invisible(NULL))
  owns_device <- !is.null(output_path)
  dev_no <- NULL
  if (isTRUE(owns_device)) {
    output_path <- as.character(output_path[[1L]])
    dev_no <- ng_plot_open_device(output_path, width = width, height = height, res = res)
    on.exit(ng_plot_close_device(dev_no), add = TRUE)
  }
  op <- graphics::par(mfrow = c(1L, length(traits)), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  for (i in seq_along(traits)) {
    tr <- traits[[i]]; kk <- key[[i]]
    y <- suppressWarnings(as.numeric(scored[[paste0(kk, "_mean")]]))
    x <- suppressWarnings(as.numeric(scored[[kinship_col]]))
    ok <- scored[[paste0(kk, "_check_ok")]]
    base_col <- ifelse(ok %in% FALSE, "#BBBBBB", "#1F4E78")
    p_col <- paste0(kk, "_p_beat_check")
    p <- if (p_col %in% names(scored)) suppressWarnings(as.numeric(scored[[p_col]])) else NA_real_
    alpha <- ifelse(is.finite(p), pmax(0.15, pmin(1, p)), 1)
    graphics::plot(x, y, pch = 19, col = ng_priority_plot_alpha_col(base_col, alpha),
                   xlab = "Pair kinship", ylab = paste(tr, "mid-parent"), main = tr)
    tau <- suppressWarnings(as.numeric(scored[[paste0(kk, "_check_value")]][[1L]]))
    if (length(tau) && is.finite(tau)) {
      graphics::abline(h = tau, lty = 2, lwd = 2, col = "#B00020")
    }
  }
  if (isTRUE(owns_device)) {
    ng_plot_close_device(dev_no)
    invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
  } else {
    invisible(output_path)
  }
}
