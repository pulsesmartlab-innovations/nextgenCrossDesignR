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

# Draw a check reference line on whichever axis carries the check's mean value. A check has a
# value on the mean axis and NONE on any other, so the line is perpendicular to whichever axis
# carries the mean: horizontal ("y") when the mean is plotted on y, vertical ("x") when the mean
# is plotted on x (e.g. a diversity-on-y-vs-mean-GEBV-on-x scatter). A plot with NEITHER axis
# carrying a mean has no honest way to place the check at all, so `mean_axis = NULL` (the
# default) draws nothing -- there is deliberately no default that guesses "y", because that
# would work silently on every plot built so far and fail silently on the first one that puts
# the mean on x. Label placement follows orientation so it never collides with the axis it sits
# on: side 4 (right margin) for a horizontal line, side 3 (top margin) for a vertical one.
ng_plot_check_reference_line <- function(value, label = NULL, mean_axis = NULL) {
  if (is.null(mean_axis) || is.null(value) || !is.finite(value)) return(invisible(NULL))
  mean_axis <- match.arg(mean_axis, c("y", "x"))
  if (identical(mean_axis, "y")) {
    graphics::abline(h = value, lty = 2, lwd = 2, col = "#B00020")
    if (!is.null(label)) {
      graphics::mtext(sprintf("check: %s", label), side = 4, at = value, las = 1, cex = 0.7,
                      col = "#B00020")
    }
  } else {
    graphics::abline(v = value, lty = 2, lwd = 2, col = "#B00020")
    if (!is.null(label)) {
      graphics::mtext(sprintf("check: %s", label), side = 3, at = value, las = 1, cex = 0.7,
                      col = "#B00020")
    }
  }
  invisible(NULL)
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
  ng_plot_check_reference_line(check_line, label = if (is.null(check_line)) NULL else
                               (if (is.null(check_label)) "reference" else check_label),
                               mean_axis = "y")

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

# Mirrors ng_rank_normalize() (R/04_optimizers.R) exactly, so a check can be placed on the SAME
# standardized axis WITHOUT perturbing the candidates it is being compared against. The
# candidates' z is: orient by direction -> rank (ties = average) -> (r-0.5)/n -> qnorm -> zero-
# fill non-finite rows -> standardize (mean/sd) over ALL rows including the zero-filled ones.
# The check never joins that ranked vector (it is not a cross); instead its fractional position
# is read off the candidates' empirical CDF with a Blom-style continuity correction --
# (#candidates <= check + 0.5) / n -- then pushed through qnorm and the SAME standardization
# (mean/sd) the candidates' vector used, so it lands on their exact axis.
ng_check_rank_axis_z <- function(x, bigger_is_better, tau_raw) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  n_ok <- sum(ok)
  if (n_ok <= 1L || !is.finite(tau_raw)) return(NA_real_)
  value <- if (isTRUE(bigger_is_better)) x[ok] else -x[ok]
  tau_oriented <- if (isTRUE(bigger_is_better)) tau_raw else -tau_raw
  r <- rank(value, ties.method = "average")
  p <- (r - 0.5) / length(r)
  out <- rep(0, length(x))
  out[ok] <- stats::qnorm(p)
  bad <- !is.finite(out)
  if (any(bad)) out[bad] <- min(out[!bad], na.rm = TRUE) - 1
  mu <- mean(out)
  sd_out <- stats::sd(out)
  if (!is.finite(sd_out) || sd_out <= 0) return(NA_real_)
  p_check <- (sum(value <= tau_oriented) + 0.5) / length(value)
  (stats::qnorm(p_check) - mu) / sd_out
}

# The y value at which to draw the check reference line. multi_trait_score is z-normalised by
# ng_add_multitrait_score() (R/19_multi_trait_selection.R), but WHICH standardization applies
# depends on the run's actual method family, and a check must be placed on the axis that
# actually exists rather than a plausible-looking substitute:
#
#   - "rank_threshold" (auto), "threshold", "weighted_index" (weighted): the score is
#     z %*% weights where z = ng_rank_normalize() -- a RANK transform of the candidates. A check
#     has no rank of its own (it never enters that vector), but its VALUE can be positioned
#     within the candidate distribution by quantile lookup -- see ng_check_rank_axis_z(). This
#     needs the candidates' own raw trait values, so `candidate_scores` (the scored/candidate
#     table) is required for this family; without it, no line.
#   - "economic_index", "desired_gain": the score is value_z %*% <solved coefficients>. value_z
#     is a plain affine standardization (center/scale/sign, retained in multi_trait_meta), so
#     the check maps through the identical transform. The combining vector must be the SOLVED
#     coefficients actually used for scoring (multi_trait_meta$economic_index_coefficients /
#     $desired_gain_coefficients) -- NOT the raw input weights, which are a different vector in
#     general.
#   - anything else, or any case whose required parameters are missing: NA_real_. No plausible-
#     looking number is drawn for a case that cannot be placed exactly.
ng_check_line_value <- function(trait_check_reference, multi_trait_meta = NULL, trait = NULL,
                                candidate_scores = NULL, trait_value_metric = NULL) {
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  if (!nrow(spec)) return(NA_real_)
  val <- function(tr) {
    ck <- spec$check[[match(tr, spec$trait)]]
    suppressWarnings(as.numeric(trait_check_reference$values[[tr]][[ck]]))
  }
  tr_list <- if (is.null(trait)) spec$trait else trait
  if (!all(tr_list %in% spec$trait)) return(NA_real_)

  family <- as.character(if (is.null(multi_trait_meta$family)) "" else multi_trait_meta$family)

  # value_z placement, for the LINEAR families (economic_index / desired_gain) only.
  centers <- multi_trait_meta$value_centers
  scales <- multi_trait_meta$value_scales
  signs <- multi_trait_meta$value_signs
  has_value_z_transform <- !is.null(centers) && !is.null(scales) && !is.null(signs)
  z_of_value_z <- function(tr) {
    tau <- val(tr)
    if (!length(tau) || !is.finite(tau)) return(NA_real_)
    if (!has_value_z_transform || !(tr %in% names(centers)) || !(tr %in% names(scales)) ||
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

  # Quantile placement onto the candidates' rank-normalized axis, for the RANK families only.
  # The check's tau is on the MEAN scale (ng_check_reference_value() resolves it onto whatever
  # source produced the cross MEANS -- R/51:69-85), never on `traits_df$column` (typically
  # `<trait>_value`, i.e. mean +/- i*sigma under the default trait_value_metric = "usefulness",
  # a pure variance for pmv/vpm, or a relationship distance for parent_distance). Ranking the
  # check against `_value` compares it to a different quantity than the one it actually is,
  # which silently shifts which side of the check a cross appears to fall on. The fix ranks
  # against the trait's `_mean` column instead: the rank-normal z's mu/sd depend only on n and
  # tie structure, not on which variable was ranked, so this is an exact placement, not an
  # approximation. `spec$column_key` is the sanitised lookup key ng_attach_check_reference() (and
  # the cross table itself) already use; a hand-built trait_check_reference without a
  # `column_key` column falls back to the raw trait name, exactly like R/51's own `key` fallback.
  traits_df <- multi_trait_meta$traits
  z_of_rank <- function(tr) {
    tau <- val(tr)
    if (!length(tau) || !is.finite(tau)) return(NA_real_)
    if (is.null(candidate_scores) || is.null(traits_df) || is.null(signs) ||
        !(tr %in% names(signs))) {
      return(NA_real_)
    }
    idx <- match(tr, traits_df$trait)
    if (is.na(idx)) return(NA_real_)
    key <- if ("column_key" %in% names(spec)) spec$column_key[[match(tr, spec$trait)]] else tr
    col <- paste0(key, "_mean")
    if (is.null(col) || !(col %in% names(candidate_scores))) return(NA_real_)
    sign_i <- suppressWarnings(as.numeric(signs[[tr]]))
    if (!is.finite(sign_i)) return(NA_real_)
    x <- suppressWarnings(as.numeric(candidate_scores[[col]]))
    ng_check_rank_axis_z(x, bigger_is_better = sign_i > 0, tau_raw = tau)
  }

  if (family %in% c("economic_index", "desired_gain")) {
    # value_z's center/scale (R/19:577-585) are computed from `traits$column`, i.e. `_value`
    # under the run's trait_value_metric. `_value = mean_value` ONLY when trait_value_metric ==
    # "mean" (R/39_cross_prediction_runner.R ng_run_cp_trait_value(): every other metric adds a
    # sign*i*sqrt(variance) term, or is a variance/distance in its own right, none of which is an
    # affine function of the mean). The check's tau is always mean-scale, so remapping it through
    # value_z's center/scale is only exact when trait_value_metric == "mean"; there is no
    # approximate substitute, so every other metric refuses rather than drawing a plausible-but-
    # wrong line.
    metric <- tolower(as.character(if (is.null(trait_value_metric)) "" else trait_value_metric[[1L]]))
    if (!identical(metric, "mean")) return(NA_real_)
    coef_vec <- if (identical(family, "economic_index")) multi_trait_meta$economic_index_coefficients
                else multi_trait_meta$desired_gain_coefficients
    if (is.null(coef_vec) || !length(coef_vec)) return(NA_real_)
    tr <- intersect(tr_list, names(coef_vec))
    if (!length(tr)) return(NA_real_)
    # A check on only SOME of the index's traits cannot be plotted: the omitted traits' terms
    # would be silently dropped from the sum, which is mathematically identical to imputing each
    # missing trait's value_z at exactly the population median (value_z's own center) -- a
    # plausible-looking substitute presented as if it were the check's real position. Refuse
    # unless every index trait is represented among the checked (+ requested) traits.
    if (length(tr) < length(coef_vec)) return(NA_real_)
    v <- vapply(tr, z_of_value_z, numeric(1))
    if (any(!is.finite(v))) return(NA_real_)
    return(sum(v * as.numeric(coef_vec[tr])))
  }
  if (family %in% c("rank_threshold", "threshold", "weighted_index")) {
    # D4: candidates' axis position for this family is z = ng_rank_normalize(<trait>_value)
    # (R/19_multi_trait_selection.R:568). `<trait>_value` is USEFULNESS (mean +/- i*sigma) under
    # the default trait_value_metric, a DIFFERENT distribution from the check, which has no
    # usefulness at all (a check is never crossed, so it has no variance -- usefulness needs
    # one). ng_check_rank_axis_z() below only ever ranks the check's (mean-scale) tau against the
    # candidates' `_mean` column, which is only the SAME axis the candidates were ranked on when
    # trait_value_metric == "mean" (then `_value` IS `_mean`, so rank(_value) == rank(_mean)
    # exactly). Every other metric draws no line, matching the gate the economic branch above
    # already applies (line 315) for the identical reason.
    metric <- tolower(as.character(if (is.null(trait_value_metric)) "" else trait_value_metric[[1L]]))
    if (!identical(metric, "mean")) return(NA_real_)
    w <- multi_trait_meta$weights
    if (is.null(w) || !length(w)) return(NA_real_)
    tr <- intersect(tr_list, names(w))
    if (!length(tr)) return(NA_real_)
    # D5: mirrors the economic branch's refusal (line 326) -- a check on only SOME of the
    # weighted traits would silently drop the omitted traits' terms from the sum, which is
    # mathematically identical to imputing each missing trait's z at exactly 0 (the rank-normal
    # axis's own mean), a plausible-looking substitute for a position the check does not have.
    if (length(tr) < length(w)) return(NA_real_)
    v <- vapply(tr, z_of_rank, numeric(1))
    if (any(!is.finite(v))) return(NA_real_)
    return(sum(v * as.numeric(w[tr])))
  }
  NA_real_
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
  # traits, key AND checks must all be filtered by `keep` together -- indexing the
  # unfiltered spec inside the loop mislabels every panel after a dropped trait.
  traits <- spec$trait[keep]; key <- key[keep]; checks <- as.character(spec$check)[keep]
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
    ng_plot_check_reference_line(if (length(tau)) tau[[1L]] else NA_real_,
                                 label = checks[[i]], mean_axis = "y")
  }
  if (isTRUE(owns_device)) {
    ng_plot_close_device(dev_no)
    invisible(normalizePath(output_path, winslash = "/", mustWork = TRUE))
  } else {
    invisible(output_path)
  }
}
