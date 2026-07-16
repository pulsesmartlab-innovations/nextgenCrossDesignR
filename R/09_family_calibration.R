ng_bind_rows_fill <- function(x) {
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

ng_cheap_cross_screen <- function(geno,
                                  effects,
                                  ids = rownames(geno),
                                  adjusted_pheno = NULL,
                                  blue = NULL,
                                  blup = NULL,
                                  selection_prop = 0.10,
                                  min_effect_reliability = 0.35,
                                  include_self = FALSE) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- as.character(ids)
  geno <- ng_check_same_ids(geno, ids, "geno")
  pairs <- ng_make_pairs(ids, include_self = include_self)
  mean_source <- ng_choose_mean_source(
    geno = geno,
    effects = effects,
    adjusted_pheno = adjusted_pheno,
    blue = blue,
    blup = blup,
    ids = ids,
    min_reliability = min_effect_reliability
  )
  K <- ng_parent_kinship(geno)
  rel_var <- ng_pair_relationship_variance(pairs, K)
  p1 <- match(pairs$parent1, ids)
  p2 <- match(pairs$parent2, ids)
  gebv <- setNames(ng_predict_gebv(geno, effects), ids)
  i <- ng_selection_intensity(selection_prop)

  out <- pairs
  out$mean_source <- mean_source$source
  out$effect_reliability <- mean_source$reliability
  out$cross_mean <- 0.5 * (mean_source$value[p1] + mean_source$value[p2])
  out$mpv <- 0.5 * (gebv[p1] + gebv[p2])
  out$var_simple <- rel_var$var_simple
  out$pair_kinship <- rel_var$pair_kinship
  out$uc_var_simple <- out$cross_mean + i * sqrt(pmax(out$var_simple, 0))
  if (!is.null(adjusted_pheno)) {
    adj <- ng_match_vector(adjusted_pheno, ids, "adjusted_pheno")
    out$cross_mean_adjusted_pheno <- 0.5 * (adj[p1] + adj[p2])
  }
  out
}

ng_select_calibration_pair_rows <- function(cheap_scores,
                                            n_families = 100L,
                                            random_n = NULL,
                                            top_per_metric = 5L,
                                            bottom_per_metric = 2L,
                                            metric_cols = c("cross_mean", "mpv", "var_simple", "uc_var_simple"),
                                            seed = NULL) {
  cheap_scores <- as.data.frame(cheap_scores, stringsAsFactors = FALSE)
  n <- nrow(cheap_scores)
  if (!n) ng_stop("cheap_scores has no rows")
  n_families <- max(1L, min(n, as.integer(n_families)))
  if (is.null(random_n)) random_n <- ceiling(n_families / 2)
  random_n <- max(0L, min(n_families, as.integer(random_n)))
  top_per_metric <- max(0L, as.integer(top_per_metric))
  bottom_per_metric <- max(0L, as.integer(bottom_per_metric))
  metric_cols <- intersect(metric_cols, names(cheap_scores))
  if (!length(metric_cols)) metric_cols <- "uc_var_simple"

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  if (!is.null(seed)) set.seed(seed)
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  chosen <- integer(0)
  if (random_n > 0L) chosen <- sample.int(n, random_n)
  for (col in metric_cols) {
    x <- cheap_scores[[col]]
    ok <- which(is.finite(x))
    if (!length(ok)) next
    if (top_per_metric > 0L) {
      chosen <- union(chosen, head(ok[order(x[ok], decreasing = TRUE)], top_per_metric))
    }
    if (bottom_per_metric > 0L) {
      chosen <- union(chosen, head(ok[order(x[ok], decreasing = FALSE)], bottom_per_metric))
    }
  }
  if (length(chosen) < n_families) {
    rest <- setdiff(seq_len(n), chosen)
    add_n <- min(length(rest), n_families - length(chosen))
    if (add_n > 0L) chosen <- c(chosen, sample(rest, add_n))
  }
  if (length(chosen) > n_families) {
    protected <- chosen
    chosen <- protected[seq_len(n_families)]
  }
  as.integer(chosen)
}

ng_add_gms_vpm_scores <- function(scores,
                                  geno,
                                  effects,
                                  marker_map,
                                  row_index = NULL,
                                  max_markers = Inf,
                                  ncores = 1L,
                                  nBLASthreads = NULL) {
  scores$gms_vpm <- NA_real_
  scores$gms_status <- "not_run"
  if (!ng_has_optional_pkg("genomicMateSelectR")) {
    scores$gms_status <- "package_unavailable"
    return(scores)
  }

  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- rownames(geno)
  if (ncol(geno) > max_markers) {
    scores$gms_status <- "skipped_marker_limit"
    return(scores)
  }
  if (is.null(row_index)) row_index <- seq_len(nrow(scores))
  row_index <- as.integer(row_index)
  row_index <- row_index[is.finite(row_index) & row_index >= 1L & row_index <= nrow(scores)]
  if (!length(row_index)) return(scores)

  marker_map <- ng_prepare_marker_map(marker_map, colnames(geno))
  marker_map <- marker_map[match(colnames(geno), marker_map$marker), , drop = FALSE]
  recomb_decay <- ng_recomb_decay_matrix(marker_map)
  beta <- effects$beta[colnames(geno)]
  beta[!is.finite(beta)] <- 0

  hap <- as.matrix(geno / 2)
  hap[!is.finite(hap)] <- 0
  hap <- pmin(1, pmax(0, hap))
  dim(hap) <- c(length(ids), ncol(geno))
  dimnames(hap) <- dimnames(geno)
  hap_n <- length(ids)
  haplo_mat <- matrix(NA_real_, nrow = 2 * hap_n, ncol = ncol(hap))
  colnames(haplo_mat) <- colnames(hap)
  rownames(haplo_mat) <- as.vector(rbind(paste0(ids, "_HapA"), paste0(ids, "_HapB")))
  haplo_mat[seq(1L, nrow(haplo_mat), by = 2L), ] <- hap
  haplo_mat[seq(2L, nrow(haplo_mat), by = 2L), ] <- hap

  if (!is.null(nBLASthreads) && requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(nBLASthreads)
  }
  for (idx in row_index) {
    sire <- as.character(scores$parent1[idx])
    dam <- as.character(scores$parent2[idx])
    f1_id <- paste0("F1_", idx)
    pair_hap <- rbind(
      haplo_mat[paste0(sire, "_HapA"), , drop = FALSE],
      haplo_mat[paste0(dam, "_HapA"), , drop = FALSE]
    )
    rownames(pair_hap) <- paste0(f1_id, c("_HapA", "_HapB"))
    x <- colSums(pair_hap)
    seg <- names(x[x > 0 & x < 2])
    if (!length(seg)) {
      scores$gms_vpm[idx] <- 0
      scores$gms_status[idx] <- "ok"
      next
    }
    val <- tryCatch({
      progeny_ld <- ng_optional_pkg_fun("genomicMateSelectR", "calcCrossLD")(
        f1_id,
        f1_id,
        recombFreqMat = recomb_decay[seg, seg, drop = FALSE],
        haploMat = pair_hap[, seg, drop = FALSE]
      )
      2 * as.numeric(ng_optional_pkg_fun("genomicMateSelectR", "quadform")(
        D = progeny_ld,
        x = beta[seg],
        y = beta[seg]
      ))
    }, error = function(e) {
      attr(scores, "gms_error") <- conditionMessage(e)
      NA_real_
    })
    if (is.finite(val)) {
      scores$gms_vpm[idx] <- pmax(val, 0)
      scores$gms_status[idx] <- "ok"
    } else {
      scores$gms_status[idx] <- "failed"
    }
  }
  scores
}

ng_family_metric_registry <- function(scores) {
  mean_cols <- c(
    "cross_mean", "cross_mean_gebv", "cross_mean_adjusted_pheno", "cross_mean_blend", "mpv",
    "popvar_mu", "simple_mpv", "simple_usefa_mean"
  )
  variance_cols <- c(
    "var_simple", "dh_recomb_var", "dh_pmv_var",
    "var_simple_cal", "dh_recomb_var_cal", "dh_pmv_var_cal",
    "ng_portfolio_var", "popvar_varG", "simple_usefa_var", "gms_vpm"
  )
  usefulness_cols <- c(
    "uc_recomb", "uc_dh",
    "uc_recomb_gebv", "uc_dh_gebv",
    "uc_recomb_adj", "uc_dh_adj",
    "uc_recomb_blend", "uc_dh_blend",
    "etk_var_simple_cal", "etk_dh_recomb_var_cal", "etk_dh_pmv_var_cal",
    "etk_var_simple_gebv_cal", "etk_dh_recomb_var_gebv_cal", "etk_dh_pmv_var_gebv_cal",
    "etk_var_simple_adj_cal", "etk_dh_recomb_var_adj_cal", "etk_dh_pmv_var_adj_cal",
    "etk_var_simple_blend_cal", "etk_dh_recomb_var_blend_cal", "etk_dh_pmv_var_blend_cal",
    "ng_portfolio_score",
    "popvar_uc", "popvar_musp_high", "simple_usefa"
  )
  registry <- data.frame(
    score_col = c(mean_cols, variance_cols, usefulness_cols),
    metric_group = c(
      rep("mean", length(mean_cols)),
      rep("variance", length(variance_cols)),
      rep("usefulness", length(usefulness_cols))
    ),
    stringsAsFactors = FALSE
  )
  registry <- registry[registry$score_col %in% names(scores), , drop = FALSE]
  registry
}

ng_family_metric_summary <- function(families,
                                     registry = NULL,
                                     target_cols = c("realized_mean", "realized_var", "realized_top10", "realized_max"),
                                     top_prop = 0.10) {
  families <- as.data.frame(families, stringsAsFactors = FALSE)
  if (is.null(registry)) registry <- ng_family_metric_registry(families)
  target_cols <- intersect(target_cols, names(families))
  if (!nrow(registry) || !length(target_cols)) return(data.frame())

  out <- list()
  k <- 1L
  for (i in seq_len(nrow(registry))) {
    score_col <- registry$score_col[i]
    x <- as.numeric(families[[score_col]])
    for (target_col in target_cols) {
      y <- as.numeric(families[[target_col]])
      ok <- is.finite(x) & is.finite(y)
      n <- sum(ok)
      row <- data.frame(
        score_col = score_col,
        metric_group = registry$metric_group[i],
        target = target_col,
        n = n,
        pearson = NA_real_,
        spearman = NA_real_,
        slope = NA_real_,
        intercept = NA_real_,
        linear_rmse = NA_real_,
        raw_rmse = NA_real_,
        top_n = NA_integer_,
        top_pred_mean = NA_real_,
        all_mean = if (n) mean(y[ok]) else NA_real_,
        top_pred_delta = NA_real_,
        top_overlap = NA_real_,
        primary_target = FALSE,
        stringsAsFactors = FALSE
      )
      row$primary_target <- (row$metric_group == "mean" && target_col == "realized_mean") ||
        (row$metric_group == "variance" && target_col == "realized_var") ||
        (row$metric_group == "usefulness" && target_col %in% c("realized_top10", "realized_max"))
      if (n >= 3L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
        fit <- stats::lm(y[ok] ~ x[ok])
        pred <- stats::fitted(fit)
        row$pearson <- suppressWarnings(stats::cor(x[ok], y[ok], method = "pearson"))
        row$spearman <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
        row$slope <- unname(stats::coef(fit)[2])
        row$intercept <- unname(stats::coef(fit)[1])
        row$linear_rmse <- sqrt(mean((y[ok] - pred)^2))
        row$raw_rmse <- sqrt(mean((y[ok] - x[ok])^2))
        top_n <- max(1L, ceiling(n * top_prop))
        ord_x <- order(x[ok], decreasing = TRUE)
        ord_y <- order(y[ok], decreasing = TRUE)
        top_x <- ord_x[seq_len(top_n)]
        top_y <- ord_y[seq_len(top_n)]
        row$top_n <- top_n
        row$top_pred_mean <- mean(y[ok][top_x])
        row$top_pred_delta <- row$top_pred_mean - row$all_mean
        row$top_overlap <- length(intersect(top_x, top_y)) / top_n
      }
      out[[k]] <- row
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}
