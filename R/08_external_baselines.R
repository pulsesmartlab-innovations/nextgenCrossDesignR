ng_external_methods_needed <- function(methods) {
  methods <- as.character(methods)
  meta <- any(grepl("^ng_meta_(portfolio|selector|router)", methods, perl = TRUE))
  frontier <- any(grepl("^ng_frontier_policy", methods, perl = TRUE))
  list(
    popvar = any(grepl("^popvar_", methods)) || meta || frontier,
    simplemating_mpv = any(grepl("^simple_mpv", methods)),
    simplemating_usefa = any(grepl("^simple_usefa", methods)) || meta || frontier
  )
}

ng_add_external_baseline_scores <- function(scores,
                                            geno,
                                            effects,
                                            marker_map,
                                            adjusted_pheno = NULL,
                                            methods = character(),
                                            selection_prop = 0.10,
                                            popvar_tail_p = selection_prop,
                                            simplemating_generation = 1L,
                                            simplemating_threads = 1L,
                                            n_crosses = NULL,
                                            shortlist_n = NULL,
                                            shortlist_multiplier = NULL,
                                            shortlist_score_col = "etk_pmv_cal",
                                            popvar_engine = "auto",
                                            simplemating_engine = "auto") {
  need <- ng_external_methods_needed(methods)
  row_index <- ng_external_shortlist_indices(
    scores = scores,
    n_crosses = n_crosses,
    shortlist_n = shortlist_n,
    shortlist_multiplier = shortlist_multiplier,
    score_col = shortlist_score_col
  )
  if (isTRUE(need$popvar)) {
    scores <- ng_add_popvar_scores(
      scores = scores,
      geno = geno,
      effects = effects,
      marker_map = marker_map,
      adjusted_pheno = adjusted_pheno,
      tail_p = popvar_tail_p,
      row_index = row_index,
      engine = popvar_engine
    )
  }
  if (isTRUE(need$simplemating_mpv) || isTRUE(need$simplemating_usefa)) {
    scores <- ng_add_simplemating_scores(
      scores = scores,
      geno = geno,
      effects = effects,
      marker_map = marker_map,
      adjusted_pheno = adjusted_pheno,
      prop_sel = selection_prop,
      generation = simplemating_generation,
      n_threads = simplemating_threads,
      include_mpv = need$simplemating_mpv,
      include_usefa = need$simplemating_usefa,
      row_index = row_index,
      engine = simplemating_engine
    )
  }
  scores
}

ng_external_shortlist_indices <- function(scores,
                                          n_crosses = NULL,
                                          shortlist_n = NULL,
                                          shortlist_multiplier = NULL,
                                          score_col = "etk_pmv_cal") {
  n <- nrow(scores)
  if (!is.null(shortlist_n) && is.finite(shortlist_n) && shortlist_n > 0) {
    k <- min(n, as.integer(shortlist_n))
  } else if (!is.null(shortlist_multiplier) && is.finite(shortlist_multiplier) &&
             shortlist_multiplier > 0 && !is.null(n_crosses) &&
             is.finite(n_crosses) && n_crosses > 0) {
    k <- min(n, max(as.integer(n_crosses), ceiling(shortlist_multiplier * n_crosses)))
  } else {
    return(seq_len(n))
  }
  score_cols <- ng_external_shortlist_score_cols(score_col)
  valid_cols <- score_cols[score_cols %in% names(scores)]
  if (!length(valid_cols)) {
    fallback <- c("usefulness_pmv_gebv", "usefulness_pmv", "usefulness_vpm_gebv", "usefulness_vpm")
    valid_cols <- fallback[fallback %in% names(scores)]
  }
  if (!length(valid_cols)) return(seq_len(k))

  selected <- integer()
  for (col in valid_cols) {
    x <- scores[[col]]
    ok <- is.finite(x)
    if (!any(ok)) next
    nk <- min(k, sum(ok))
    selected <- unique(c(
      selected,
      which(ok)[order(x[ok], decreasing = TRUE)[seq_len(nk)]]
    ))
  }
  if (!length(selected)) return(seq_len(k))
  selected
}

ng_external_shortlist_score_cols <- function(score_col) {
  x <- unlist(strsplit(as.character(score_col), ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  unique(x[nzchar(x) & !is.na(x)])
}

ng_add_popvar_scores <- function(scores,
                                 geno,
                                 effects,
                                 marker_map,
                                 adjusted_pheno = NULL,
                                 tail_p = 0.10,
                                 self_gen = 0,
                                 dh = TRUE,
                                 row_index = NULL,
                                 fallback = c("native_proxy", "none"),
                                 engine = c("auto", "native", "external")) {
  fallback <- match.arg(fallback)
  engine <- match.arg(engine)
  scores$popvar_mu <- NA_real_
  scores$popvar_varG <- NA_real_
  scores$popvar_musp_low <- NA_real_
  scores$popvar_musp_high <- NA_real_
  scores$popvar_uc <- NA_real_
  scores$popvar_status <- "not_run"

  if (identical(engine, "native")) {
    return(ng_apply_popvar_native_proxy(scores, tail_p = tail_p, status = "native_proxy"))
  }

  if (!ng_has_optional_pkg("PopVar")) {
    if (identical(engine, "auto") && identical(fallback, "native_proxy")) {
      return(ng_apply_popvar_native_proxy(scores, tail_p = tail_p, status = "native_proxy_package_unavailable"))
    }
    scores$popvar_status <- "package_unavailable"
    return(scores)
  }

  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- rownames(geno)
  marker_map <- ng_prepare_marker_map(marker_map, colnames(geno))
  marker_map <- marker_map[match(colnames(geno), marker_map$marker), , drop = FALSE]
  if (is.null(row_index)) row_index <- seq_len(nrow(scores))
  row_index <- as.integer(row_index)
  row_index <- row_index[is.finite(row_index) & row_index >= 1L & row_index <= nrow(scores)]
  if (!length(row_index)) return(scores)

  if (any(!stats::complete.cases(geno))) {
    if (identical(engine, "auto") && identical(fallback, "native_proxy")) {
      return(ng_apply_popvar_native_proxy(scores, tail_p = tail_p, status = "native_proxy_missing_genotypes"))
    }
    scores$popvar_status <- "missing_genotypes"
    return(scores)
  }
  if (any(!(as.vector(geno) %in% c(0, 2)))) {
    if (identical(engine, "auto") && identical(fallback, "native_proxy")) {
      return(ng_apply_popvar_native_proxy(scores, tail_p = tail_p, status = "native_proxy_requires_inbred_0_2_genotypes"))
    }
    scores$popvar_status <- "requires_inbred_0_2_genotypes"
    return(scores)
  }

  beta <- effects$beta[colnames(geno)]
  beta[!is.finite(beta)] <- 0
  marker_effects <- data.frame(
    marker = colnames(geno),
    Trait = as.numeric(beta),
    stringsAsFactors = FALSE
  )
  y <- if (is.null(adjusted_pheno)) {
    rep(NA_real_, length(ids))
  } else {
    ng_match_vector(adjusted_pheno, ids, "adjusted_pheno")
  }
  y_in <- data.frame(Entry = ids, Trait = as.numeric(y), stringsAsFactors = FALSE)
  map_in <- data.frame(
    marker = marker_map$marker,
    chr = marker_map$chr,
    pos = marker_map$pos_cm,
    stringsAsFactors = FALSE
  )
  crossing_table <- data.frame(
    parent1 = as.character(scores$parent1[row_index]),
    parent2 = as.character(scores$parent2[row_index]),
    stringsAsFactors = FALSE
  )

  out <- tryCatch({
    ng_optional_pkg_fun("PopVar", "pop_predict2")(
      M = geno - 1,
      y.in = y_in,
      marker.effects = marker_effects,
      map.in = map_in,
      crossing.table = crossing_table,
      tail.p = tail_p,
      self.gen = self_gen,
      DH = dh,
      models = "rrBLUP"
    )
  }, error = function(e) {
    attr(scores, "popvar_error") <- conditionMessage(e)
    NULL
  })
  if (is.null(out) || !nrow(out)) {
    if (identical(engine, "auto") && identical(fallback, "native_proxy")) {
      return(ng_apply_popvar_native_proxy(scores, tail_p = tail_p, status = "native_proxy_exact_failed"))
    }
    scores$popvar_status <- "failed"
    return(scores)
  }

  idx <- match(ng_pair_key(scores$parent1, scores$parent2), ng_pair_key(out$parent1, out$parent2))
  hit <- !is.na(idx)
  scores$popvar_mu[hit] <- out$pred_mu[idx[hit]]
  scores$popvar_varG[hit] <- pmax(out$pred_varG[idx[hit]], 0)
  scores$popvar_musp_low[hit] <- out$pred_musp_low[idx[hit]]
  scores$popvar_musp_high[hit] <- out$pred_musp_high[idx[hit]]
  i <- ng_selection_intensity(tail_p)
  scores$popvar_uc <- scores$popvar_mu + i * sqrt(pmax(scores$popvar_varG, 0))
  scores$popvar_status[hit] <- "ok"
  scores
}

ng_add_simplemating_scores <- function(scores,
                                       geno,
                                       effects,
                                       marker_map,
                                       adjusted_pheno = NULL,
                                       prop_sel = 0.10,
                                       type = "DH",
                                       generation = 1L,
                                       n_threads = 1L,
                                       include_mpv = TRUE,
                                       include_usefa = TRUE,
                                       row_index = NULL,
                                       fallback = c("native_proxy", "none"),
                                       engine = c("auto", "native", "external")) {
  fallback <- match.arg(fallback)
  engine <- match.arg(engine)
  scores$simple_mpv <- NA_real_
  scores$simple_usefa_mean <- NA_real_
  scores$simple_usefa_var <- NA_real_
  scores$simple_usefa_sd <- NA_real_
  scores$simple_usefa <- NA_real_
  scores$simple_status <- "not_run"

  if (identical(engine, "native")) {
    return(ng_apply_simplemating_native_proxy(scores, prop_sel = prop_sel, status = "native_proxy"))
  }

  if (!ng_has_optional_pkg("SimpleMating")) {
    if (identical(engine, "auto") && identical(fallback, "native_proxy")) {
      return(ng_apply_simplemating_native_proxy(scores, prop_sel = prop_sel, status = "native_proxy_package_unavailable"))
    }
    scores$simple_status <- "package_unavailable"
    return(scores)
  }

  geno <- ng_as_numeric_matrix(geno, "geno")
  ids <- rownames(geno)
  marker_map <- ng_prepare_marker_map(marker_map, colnames(geno))
  marker_map <- marker_map[match(colnames(geno), marker_map$marker), , drop = FALSE]
  if (is.null(row_index)) row_index <- seq_len(nrow(scores))
  row_index <- as.integer(row_index)
  row_index <- row_index[is.finite(row_index) & row_index >= 1L & row_index <= nrow(scores)]
  if (!length(row_index)) return(scores)
  beta <- effects$beta[colnames(geno)]
  beta[!is.finite(beta)] <- 0

  K <- ng_parent_kinship(geno)
  mate_plan <- data.frame(
    Parent1 = as.character(scores$parent1[row_index]),
    Parent2 = as.character(scores$parent2[row_index]),
    stringsAsFactors = FALSE
  )
  criterion <- if (!is.null(adjusted_pheno)) {
    setNames(ng_match_vector(adjusted_pheno, ids, "adjusted_pheno"), ids)
  } else {
    setNames(ng_predict_gebv(geno, effects), ids)
  }
  crit_df <- data.frame(Id = ids, Criterion = as.numeric(criterion[ids]), stringsAsFactors = FALSE)

  mpv_error <- NULL
  usefa_error <- NULL
  if (isTRUE(include_mpv)) {
    mid_parent_value <- tryCatch({
      utils::capture.output({
        tmp <- ng_optional_pkg_fun("SimpleMating", "getMPV")(MatePlan = mate_plan, Criterion = crit_df, K = K)
      })
      tmp
    }, error = function(e) {
      mpv_error <<- conditionMessage(e)
      NULL
    })
    if (!is.null(mid_parent_value) && nrow(mid_parent_value)) {
      idx <- match(ng_pair_key(scores$parent1, scores$parent2), ng_pair_key(mid_parent_value$Parent1, mid_parent_value$Parent2))
      hit <- !is.na(idx)
      scores$simple_mpv[hit] <- mid_parent_value$Y[idx[hit]]
    }
  }

  if (isTRUE(include_usefa)) {
    map_in <- data.frame(
      chr = marker_map$chr,
      pos = marker_map$pos_cm,
      marker = marker_map$marker,
      stringsAsFactors = FALSE
    )
    # SimpleMating getUsefA requires the parent Markers to be homozygous-coded (0, 2, or
    # NA) because a non-phased genotype cannot resolve haplotype phase. RIL parents are
    # mostly homozygous but can carry a few unfixed heterozygous loci; we treat those
    # individual calls as missing data (NA) -- which getUsefA explicitly supports -- rather
    # than forcing whole-genome homozygosity (which would turn RILs into DH lines). If a
    # LARGE fraction is heterozygous the parents are genuinely outbred and need a
    # phased-haplotype variance path, so we warn.
    markers_usefa <- geno
    het_mask <- is.finite(markers_usefa) & !(markers_usefa %in% c(0, 2))
    het_frac <- if (length(markers_usefa)) mean(het_mask) else 0
    if (any(het_mask)) {
      markers_usefa[het_mask] <- NA_real_
      if (het_frac > 0.25) {
        warning(sprintf(
          "SimpleMating getUsefA: %.1f%% of parent genotype calls are heterozygous and were treated as missing; for substantially outbred parents use a phased-haplotype variance path.",
          100 * het_frac), call. = FALSE)
      }
    }
    attr(scores, "simple_usefa_het_frac") <- het_frac
    usefa <- tryCatch({
      utils::capture.output({
        tmp <- ng_optional_pkg_fun("SimpleMating", "getUsefA")(
          MatePlan = mate_plan,
          Markers = markers_usefa,
          addEff = as.numeric(beta),
          K = K,
          Map.In = map_in,
          propSel = prop_sel,
          Type = type,
          Generation = generation,
          n_threads = n_threads,
          display_progress = FALSE
        )
      })
      tmp
    }, error = function(e) {
      usefa_error <<- conditionMessage(e)
      NULL
    })
    if (!is.null(usefa) && length(usefa) >= 1L && nrow(usefa[[1]])) {
      u <- usefa[[1]]
      idx <- match(ng_pair_key(scores$parent1, scores$parent2), ng_pair_key(u$Parent1, u$Parent2))
      hit <- !is.na(idx)
      scores$simple_usefa_mean[hit] <- u$Mean[idx[hit]]
      scores$simple_usefa_var[hit] <- pmax(u$Variance[idx[hit]], 0)
      scores$simple_usefa_sd[hit] <- u$sd[idx[hit]]
      scores$simple_usefa[hit] <- u$Usefulness[idx[hit]]
    }
  }
  if (!is.null(mpv_error)) attr(scores, "simple_mpv_error") <- mpv_error
  if (!is.null(usefa_error)) attr(scores, "simple_usefa_error") <- usefa_error
  ok <- is.finite(scores$simple_mpv) | is.finite(scores$simple_usefa)
  if (!any(ok) && identical(engine, "auto") && identical(fallback, "native_proxy")) {
    proxy <- ng_apply_simplemating_native_proxy(scores, prop_sel = prop_sel, status = "native_proxy_exact_failed")
    if (!is.null(usefa_error)) attr(proxy, "simple_usefa_error") <- usefa_error
    if (!is.null(mpv_error)) attr(proxy, "simple_mpv_error") <- mpv_error
    return(proxy)
  }
  scores$simple_status[ok] <- "ok"
  scores
}

ng_parent_criterion_from_cross_mean <- function(scores, criterion_col = "cross_mean") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!(criterion_col %in% names(scores))) ng_stop("scores missing AlphaMate criterion column: ", criterion_col)
  required <- c("parent1", "parent2")
  miss <- setdiff(required, names(scores))
  if (length(miss)) ng_stop("scores missing columns: ", paste(miss, collapse = ", "))
  parents <- sort(unique(c(as.character(scores$parent1), as.character(scores$parent2))))
  y <- suppressWarnings(as.numeric(scores[[criterion_col]]))
  ok <- is.finite(y)
  if (!any(ok)) ng_stop("AlphaMate criterion column has no finite values: ", criterion_col)
  x <- matrix(0, nrow = sum(ok), ncol = length(parents), dimnames = list(NULL, parents))
  p1 <- match(as.character(scores$parent1[ok]), parents)
  p2 <- match(as.character(scores$parent2[ok]), parents)
  for (i in seq_len(nrow(x))) {
    if (identical(p1[[i]], p2[[i]])) {
      x[i, p1[[i]]] <- 1
    } else {
      x[i, p1[[i]]] <- x[i, p1[[i]]] + 0.5
      x[i, p2[[i]]] <- x[i, p2[[i]]] + 0.5
    }
  }
  fit <- stats::lm.fit(x = x, y = y[ok])
  coef <- as.numeric(fit$coefficients)
  names(coef) <- parents
  if (anyNA(coef) || fit$rank < length(parents)) {
    fallback <- rowsum(rep(y[ok], 2), c(as.character(scores$parent1[ok]), as.character(scores$parent2[ok])))
    counts <- rowsum(rep(1, 2 * sum(ok)), c(as.character(scores$parent1[ok]), as.character(scores$parent2[ok])))
    approx <- as.numeric(fallback[, 1] / counts[, 1])
    names(approx) <- rownames(fallback)
    coef[is.na(coef)] <- approx[names(coef)[is.na(coef)]]
  }
  if (any(!is.finite(coef))) ng_stop("Could not reconstruct finite AlphaMate parent criterion")
  coef
}

ng_alphamate_default_executable <- function() {
  env <- Sys.getenv("NG_ALPHAMATE_EXE", unset = "")
  if (nzchar(env)) return(env)
  # AlphaMate ships as `AlphaMate.exe` on Windows and `AlphaMate` (no extension) on
  # macOS/Linux. Probe both names so discovery is OS-agnostic, preferring the
  # native name for the running platform.
  exe_names <- if (.Platform$OS.type == "windows") {
    c("AlphaMate.exe", "AlphaMate")
  } else {
    c("AlphaMate", "AlphaMate.exe")
  }
  bases <- c(
    file.path("external", "AlphaMate", "binaries"),
    file.path("..", "external", "AlphaMate", "binaries")
  )
  candidates <- as.vector(t(outer(bases, exe_names, file.path)))
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[[1]], winslash = "/", mustWork = FALSE) else exe_names[[1L]]
}

ng_alphamate_runtime_paths <- function(runtime_path = Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")) {
  # The bundled torch runtime hint is Windows-specific; only offer it on Windows so
  # the default carries no OS-specific assumptions elsewhere. Non-Windows users
  # supply any runtime dirs via NG_ALPHAMATE_RUNTIME_PATH.
  extra <- if (.Platform$OS.type == "windows") "C:/Python/Lib/site-packages/torch/lib" else character(0)
  candidates <- c(
    unlist(strsplit(runtime_path, .Platform$path.sep, fixed = TRUE), use.names = FALSE),
    extra
  )
  candidates <- unique(trimws(candidates))
  candidates[nzchar(candidates) & dir.exists(candidates)]
}

ng_write_alphamate_table <- function(path, ids, values) {
  ids <- as.character(ids)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(ids)) {
    writeLines(paste(ids[[i]], format(as.numeric(values[[i]]), digits = 17, scientific = TRUE)), con)
  }
}

ng_write_alphamate_matrix <- function(path, K) {
  ids <- rownames(K)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(ids)) {
    vals <- format(as.numeric(K[i, ]), digits = 17, scientific = TRUE)
    writeLines(paste(c(ids[[i]], vals), collapse = " "), con)
  }
}

ng_write_alphamate_spec <- function(path,
                                    n_crosses,
                                    target_degree = 45,
                                    number_of_parents = NULL,
                                    max_contributions = NULL,
                                    evol_solutions = 100L,
                                    evol_iterations = 1000L,
                                    evol_stop = 200L,
                                    n_threads = 1L,
                                    selfing_weight = -1e6) {
  lines <- c(
    "Seed, 15",
    "NrmMatrixFile, Nrm.txt",
    "SelCriterionFile, Criterion.txt",
    paste("NumberOfMatings,", as.integer(n_crosses)),
    if (!is.null(number_of_parents) && is.finite(number_of_parents)) {
      paste("NumberOfParents,", as.integer(number_of_parents))
    },
    if (!is.null(max_contributions) && is.finite(max_contributions)) {
      paste("LimitContributionsMax,", as.integer(max_contributions))
    },
    paste("TargetDegree,", as.numeric(target_degree)),
    "AllowRepeatedMatings, no",
    "AllowSelfing, no",
    paste("SelfingWeight,", as.numeric(selfing_weight)),
    paste("EvolAlgNumberOfSolutions,", as.integer(evol_solutions)),
    paste("EvolAlgNumberOfIterations,", as.integer(evol_iterations)),
    paste("EvolAlgNumberOfIterationsStop,", as.integer(evol_stop)),
    paste("NumberOfThreads,", as.integer(n_threads)),
    "Stop"
  )
  writeLines(lines[!vapply(lines, is.null, logical(1))], path)
}

ng_run_alphamate <- function(executable,
                             spec_file = "AlphaMateSpec.txt",
                             workdir,
                             runtime_path = "",
                             stdout_file = "AlphaMate.stdout.txt",
                             stderr_file = "AlphaMate.stderr.txt") {
  executable <- normalizePath(executable, winslash = "/", mustWork = TRUE)
  workdir <- normalizePath(workdir, winslash = "/", mustWork = TRUE)
  runtime_paths <- ng_alphamate_runtime_paths(runtime_path)
  old_path <- Sys.getenv("PATH", unset = "")
  if (length(runtime_paths)) {
    Sys.setenv(PATH = paste(c(runtime_paths, old_path), collapse = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  stdout_path <- file.path(workdir, stdout_file)
  stderr_path <- file.path(workdir, stderr_file)
  old_wd <- getwd()
  setwd(workdir)
  on.exit(setwd(old_wd), add = TRUE)
  status <- system2(
    command = executable,
    args = spec_file,
    stdout = stdout_path,
    stderr = stderr_path,
    wait = TRUE
  )
  if (is.null(status)) status <- 0L
  list(
    status = as.integer(status),
    stdout = stdout_path,
    stderr = stderr_path,
    runtime_paths = runtime_paths,
    executable = executable
  )
}

ng_parse_alphamate_mating_plan <- function(path) {
  if (!file.exists(path)) ng_stop("AlphaMate mating-plan output not found: ", path)
  out <- utils::read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Parent1", "Parent2")
  miss <- setdiff(required, names(out))
  if (length(miss)) ng_stop("AlphaMate mating-plan output missing columns: ", paste(miss, collapse = ", "))
  if (!("MatingsCount" %in% names(out))) out$MatingsCount <- 1L
  out$MatingsCount <- as.integer(out$MatingsCount)
  out$MatingsCount[!is.finite(out$MatingsCount) | out$MatingsCount < 1L] <- 1L
  out[rep(seq_len(nrow(out)), out$MatingsCount), , drop = FALSE]
}

ng_alphamate_pair_labels <- function(parent1, parent2, max_n = 5L) {
  labels <- paste(as.character(parent1), as.character(parent2), sep = " x ")
  labels <- unique(labels)
  if (length(labels) > max_n) {
    labels <- c(utils::head(labels, max_n), paste0("...", length(labels) - max_n, " more"))
  }
  paste(labels, collapse = ", ")
}

ng_match_alphamate_plan <- function(scores, am_plan, n_crosses) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  am_plan <- as.data.frame(am_plan, stringsAsFactors = FALSE)
  score_required <- c("parent1", "parent2")
  plan_required <- c("Parent1", "Parent2")
  score_miss <- setdiff(score_required, names(scores))
  plan_miss <- setdiff(plan_required, names(am_plan))
  if (length(score_miss)) ng_stop("scores missing columns: ", paste(score_miss, collapse = ", "))
  if (length(plan_miss)) ng_stop("AlphaMate mating plan missing columns: ", paste(plan_miss, collapse = ", "))

  p1 <- as.character(am_plan$Parent1)
  p2 <- as.character(am_plan$Parent2)
  score_key <- paste(as.character(scores$parent1), as.character(scores$parent2), sep = "\r")
  score_key_rev <- paste(as.character(scores$parent2), as.character(scores$parent1), sep = "\r")
  selected_key <- paste(p1, p2, sep = "\r")
  idx <- match(selected_key, score_key)
  missing <- is.na(idx)
  if (any(missing)) idx[missing] <- match(selected_key[missing], score_key_rev)

  self <- identical(length(p1), length(p2)) & p1 == p2
  canonical <- ifelse(p1 <= p2, paste(p1, p2, sep = "\r"), paste(p2, p1, sep = "\r"))
  repeated <- duplicated(canonical)
  unmatched <- is.na(idx)
  valid <- !(self | repeated | unmatched)
  details <- character()
  if (any(self)) {
    details <- c(details, paste0("self-cross output: ", ng_alphamate_pair_labels(p1[self], p2[self])))
  }
  if (any(repeated)) {
    details <- c(details, paste0("repeated output pairs: ", ng_alphamate_pair_labels(p1[repeated], p2[repeated])))
  }
  if (any(unmatched)) {
    details <- c(details, paste0("unmatched output pairs: ", ng_alphamate_pair_labels(p1[unmatched], p2[unmatched])))
  }
  if (sum(valid) < n_crosses) {
    if (!length(details)) {
      details <- paste0("only ", sum(valid), " valid matched crosses for requested ", n_crosses)
    }
    ng_stop("AlphaMate selected crosses could not be matched back to scores: ", paste(details, collapse = "; "))
  }
  idx[valid][seq_len(n_crosses)]
}

ng_select_alphamate <- function(scores,
                                criterion_col = "cross_mean",
                                n_crosses,
                                parent_kinship,
                                executable = ng_alphamate_default_executable(),
                                runtime_path = Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = ""),
                                target_degree = 45,
                                max_contributions = NULL,
                                number_of_parents = NULL,
                                evol_solutions = 100L,
                                evol_iterations = 1000L,
                                evol_stop = 200L,
                                n_threads = 1L,
                                selfing_weight = -1e6,
                                workdir = tempfile("ng_alphamate_"),
                                keep_files = FALSE,
                                mode = "ModeOptTarget1") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!file.exists(executable)) ng_stop("AlphaMate executable not found: ", executable)
  if (is.null(parent_kinship)) ng_stop("AlphaMate selection requires parent_kinship")
  parents <- sort(unique(c(as.character(scores$parent1), as.character(scores$parent2))))
  parent_kinship <- parent_kinship[parents, parents, drop = FALSE]
  criterion <- ng_parent_criterion_from_cross_mean(scores, criterion_col)
  criterion <- criterion[parents]
  aliases <- sprintf("AM%06d", seq_along(parents))
  names(aliases) <- parents
  parent_by_alias <- setNames(parents, aliases)
  scores_for_am <- scores
  scores_for_am$parent1 <- unname(aliases[as.character(scores_for_am$parent1)])
  scores_for_am$parent2 <- unname(aliases[as.character(scores_for_am$parent2)])
  if (anyNA(scores_for_am$parent1) || anyNA(scores_for_am$parent2)) {
    ng_stop("Could not alias AlphaMate score parent IDs")
  }
  parent_K_for_am <- parent_kinship
  rownames(parent_K_for_am) <- colnames(parent_K_for_am) <- aliases[parents]
  criterion_for_am <- criterion
  names(criterion_for_am) <- aliases[names(criterion_for_am)]
  if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(workdir)) ng_stop("Could not create AlphaMate workdir: ", workdir)
  if (!isTRUE(keep_files)) on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  ng_write_alphamate_table(file.path(workdir, "Criterion.txt"), aliases[parents], criterion_for_am[aliases[parents]])
  ng_write_alphamate_matrix(file.path(workdir, "Nrm.txt"), parent_K_for_am)
  ng_write_alphamate_spec(
    path = file.path(workdir, "AlphaMateSpec.txt"),
    n_crosses = n_crosses,
    target_degree = target_degree,
    number_of_parents = number_of_parents,
    max_contributions = max_contributions,
    evol_solutions = evol_solutions,
    evol_iterations = evol_iterations,
    evol_stop = evol_stop,
    n_threads = n_threads,
    selfing_weight = selfing_weight
  )
  run <- ng_run_alphamate(
    executable = executable,
    spec_file = "AlphaMateSpec.txt",
    workdir = workdir,
    runtime_path = runtime_path
  )
  if (!identical(run$status, 0L)) {
    stderr <- if (file.exists(run$stderr)) paste(readLines(run$stderr, warn = FALSE), collapse = "\n") else ""
    ng_stop("AlphaMate failed with exit code ", run$status, if (nzchar(stderr)) paste0(": ", stderr) else "")
  }

  plan_file <- switch(
    mode,
    ModeMaxCriterion = "MatingPlanModeMaxCriterion.txt",
    ModeMinCoancestry = "MatingPlanModeMinCoancestry.txt",
    ModeOptTarget1 = "MatingPlanModeOptTarget1.txt",
    ng_stop("Unsupported AlphaMate output mode: ", mode)
  )
  am_plan <- ng_parse_alphamate_mating_plan(file.path(workdir, plan_file))
  idx <- ng_match_alphamate_plan(scores_for_am, am_plan, n_crosses)
  selected <- ng_plan_summary(
    selected = idx,
    scores = scores_for_am,
    gain_col = criterion_col,
    parent_kinship = parent_K_for_am,
    lambda_group = 0,
    lambda_mating = 0,
    lambda_parent_use = 0,
    lambda_parent_use_input = 0,
    lambda_parent_use_mode = "alphamate",
    score_scale = NA_real_
  )
  selected$parent1 <- unname(parent_by_alias[as.character(selected$parent1)])
  selected$parent2 <- unname(parent_by_alias[as.character(selected$parent2)])
  if (anyNA(selected$parent1) || anyNA(selected$parent2)) {
    ng_stop("Could not restore AlphaMate selected parent IDs from aliases")
  }
  selected$alphamate_mode <- mode
  selected$alphamate_target_degree <- target_degree
  selected$alphamate_rank <- seq_len(nrow(selected))
  s <- attr(selected, "summary")
  s$alphamate_executable <- run$executable
  s$alphamate_runtime_path <- paste(run$runtime_paths, collapse = .Platform$path.sep)
  s$alphamate_workdir <- normalizePath(workdir, winslash = "/", mustWork = FALSE)
  s$alphamate_mode <- mode
  s$alphamate_target_degree <- target_degree
  s$alphamate_exit_code <- run$status
  s$alphamate_criterion_col <- criterion_col
  attr(selected, "summary") <- s
  selected
}

ng_select_simplemating <- function(scores,
                                   score_col,
                                   n_crosses,
                                   parent_kinship,
                                   max_crosses_per_parent = 4L,
                                   min_crosses_per_parent = 1L,
                                   max_crosses_to_search = 1e5,
                                   culling_pairwise_k = NULL) {
  if (!ng_has_optional_pkg("SimpleMating")) {
    ng_stop("SimpleMating package is not available")
  }
  if (!(score_col %in% names(scores))) ng_stop("scores missing ", score_col)
  ok <- is.finite(scores[[score_col]]) & is.finite(scores$pair_kinship)
  if (sum(ok) < n_crosses) ng_stop("Not enough finite SimpleMating score rows for ", score_col)
  data <- data.frame(
    Parent1 = as.character(scores$parent1[ok]),
    Parent2 = as.character(scores$parent2[ok]),
    Y = as.numeric(scores[[score_col]][ok]),
    K = as.numeric(scores$pair_kinship[ok]),
    stringsAsFactors = FALSE
  )
  out <- tryCatch({
    utils::capture.output({
      tmp <- ng_optional_pkg_fun("SimpleMating", "selectCrosses")(
        data = data,
        n.cross = n_crosses,
        max.cross = max_crosses_per_parent,
        min.cross = min_crosses_per_parent,
        max.cross.to.search = max_crosses_to_search,
        culling.pairwise.k = culling_pairwise_k
      )
    })
    tmp
  }, error = function(e) {
    ng_stop("SimpleMating selectCrosses failed: ", conditionMessage(e))
  })
  idx <- match(ng_pair_key(out$plan$Parent1, out$plan$Parent2), ng_pair_key(scores$parent1, scores$parent2))
  idx <- idx[!is.na(idx)]
  if (length(idx) < n_crosses) ng_stop("SimpleMating selected crosses could not be matched back to scores")
  ng_plan_summary(
    scores = scores,
    selected = idx[seq_len(n_crosses)],
    gain_col = score_col,
    parent_kinship = parent_kinship,
    lambda_group = 0,
    lambda_mating = 0,
    lambda_parent_use = 0,
    lambda_parent_use_input = 0,
    lambda_parent_use_mode = "simplemating_select",
    score_scale = NA_real_
  )
}
