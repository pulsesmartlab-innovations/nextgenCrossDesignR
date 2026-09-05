ng_run_cp_read_csv <- function(path, table_name) {
  if (is.null(path)) return(NULL)
  path <- as.character(path[[1L]])
  if (!nzchar(path) || !file.exists(path)) ng_stop(table_name, " file not found: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

ng_run_cp_first_col <- function(x, candidates, table_name) {
  hit <- match(tolower(candidates), tolower(names(x)), nomatch = 0L)
  hit <- hit[hit > 0L]
  if (!length(hit)) {
    ng_stop(table_name, " is missing one of these columns: ", paste(candidates, collapse = ", "))
  }
  names(x)[[hit[[1L]]]]
}

ng_run_cp_id_col <- function(x, id_col, table_name) {
  if (!is.null(id_col)) {
    if (!(id_col %in% names(x))) ng_stop(table_name, " is missing id_col: ", id_col)
    return(id_col)
  }
  ng_run_cp_first_col(x, c("parent", "parent_id", "id", "name", "line", "entry", "NAME"), table_name)
}

ng_run_cp_canonical_id_table <- function(x, id_col, table_name, canonical_col = "parent") {
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  original_id_col <- ng_run_cp_id_col(x, id_col, table_name)
  if (!identical(original_id_col, canonical_col)) {
    if (canonical_col %in% names(x)) {
      ng_stop(table_name, " already contains reserved internal ID column: ", canonical_col)
    }
    names(x)[match(original_id_col, names(x))] <- canonical_col
  }
  attr(x, "ng_original_id_col") <- original_id_col
  attr(x, "ng_internal_id_col") <- canonical_col
  x
}

ng_run_cp_geno_matrix <- function(geno, id_col = NULL) {
  if (is.matrix(geno)) {
    out <- ng_as_numeric_matrix(geno, "genotype")
    return(out)
  }
  geno <- as.data.frame(geno, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- ng_run_cp_id_col(geno, id_col, "genotype")
  ids <- trimws(as.character(geno[[id_col]]))
  if (any(!nzchar(ids) | is.na(ids))) ng_stop("genotype contains missing parent IDs")
  marker_cols <- setdiff(names(geno), id_col)
  if (!length(marker_cols)) ng_stop("genotype has no marker columns after id_col")
  out <- suppressWarnings(as.matrix(geno[, marker_cols, drop = FALSE]))
  storage.mode(out) <- "double"
  rownames(out) <- ids
  out
}

ng_run_cp_pheno_frame <- function(phenotype, id_col = NULL) {
  phenotype <- as.data.frame(phenotype, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- ng_run_cp_id_col(phenotype, id_col, "phenotype")
  ids <- trimws(as.character(phenotype[[id_col]]))
  if (any(!nzchar(ids) | is.na(ids))) ng_stop("phenotype contains missing parent IDs")
  rownames(phenotype) <- ids
  attr(phenotype, "ng_id_col") <- id_col
  phenotype
}

# Build the optional marker-effect TRAINING augmentation. Returns a genotype matrix aligned to
# the parents' markers whose rows are training-ONLY individuals (any parent overlap is dropped,
# so parents are never double-counted), plus their phenotypes. These individuals enlarge the
# ridge fit only and are never candidate parents. Returns NULL when no training input is given.
ng_run_cp_training_set <- function(training_genotype, training_phenotype,
                                   training_genotype_file, training_phenotype_file,
                                   training_genotype_id_col, training_phenotype_id_col,
                                   parent_ids, parent_markers, trait_columns) {
  have_geno <- !is.null(training_genotype) || !is.null(training_genotype_file)
  have_pheno <- !is.null(training_phenotype) || !is.null(training_phenotype_file)
  if (!have_geno && !have_pheno) return(NULL)
  if (have_geno != have_pheno) {
    ng_stop("training_genotype and training_phenotype must BOTH be supplied: extra training ",
            "individuals need genotypes AND phenotypes to inform marker effects")
  }
  geno_raw <- if (is.null(training_genotype)) {
    ng_run_cp_read_csv(training_genotype_file, "training_genotype")
  } else training_genotype
  pheno_raw <- if (is.null(training_phenotype)) {
    ng_run_cp_read_csv(training_phenotype_file, "training_phenotype")
  } else training_phenotype

  g <- ng_run_cp_geno_matrix(geno_raw, id_col = training_genotype_id_col)
  missing_markers <- setdiff(parent_markers, colnames(g))
  if (length(missing_markers)) {
    ng_stop("training_genotype is missing marker(s) present in the parents: ",
            paste(utils::head(missing_markers, 5L), collapse = ", "),
            if (length(missing_markers) > 5L) ", ..." else "")
  }
  g <- g[, parent_markers, drop = FALSE]         # align to parent markers, same order
  storage.mode(g) <- "double"

  p <- ng_run_cp_pheno_frame(pheno_raw, id_col = training_phenotype_id_col)
  missing_traits <- setdiff(trait_columns, names(p))
  if (length(missing_traits)) {
    ng_stop("training_phenotype is missing trait column(s) needed for effect estimation: ",
            paste(missing_traits, collapse = ", "))
  }

  common <- setdiff(intersect(rownames(g), rownames(p)), parent_ids)   # training-only, geno+pheno
  if (!length(common)) return(NULL)
  g <- g[common, , drop = FALSE]
  p <- p[common, trait_columns, drop = FALSE]
  if (anyNA(g)) {                                # mean-impute so the ridge fit is well-defined
    for (j in seq_len(ncol(g))) {
      col <- g[, j]
      if (anyNA(col)) { col[is.na(col)] <- mean(col, na.rm = TRUE); g[, j] <- col }
    }
  }
  list(geno = g, pheno = as.data.frame(p, stringsAsFactors = FALSE), ids = common)
}

ng_run_cp_marker_map <- function(marker_map,
                                 map_marker_col = NULL,
                                 map_chr_col = NULL,
                                 map_pos_col = NULL,
                                 map_pos_cm_col = NULL,
                                 map_pos_bp_col = NULL,
                                 map_position_unit = c("bp", "cM"),
                                 bp_per_cm = NULL,
                                 map_pos_cm_divisor = 1) {
  if (is.null(marker_map)) return(NULL)
  map_position_unit <- match.arg(map_position_unit)
  marker_map <- as.data.frame(marker_map, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(map_marker_col)) {
    map_marker_col <- ng_run_cp_first_col(marker_map, c("marker", "marker_id", "snp_code", "SNP_code", "id"), "marker_map")
  }
  if (!(map_marker_col %in% names(marker_map))) ng_stop("marker_map is missing map_marker_col: ", map_marker_col)
  if (is.null(map_chr_col)) {
    chr_hit <- match(tolower(c("chr", "chromosome", "Chromosome")), tolower(names(marker_map)), nomatch = 0L)
    chr_hit <- chr_hit[chr_hit > 0L]
    map_chr_col <- if (length(chr_hit)) names(marker_map)[[chr_hit[[1L]]]] else NULL
  }

  if (identical(map_position_unit, "bp")) {
    if (is.null(bp_per_cm)) {
      legacy_divisor <- suppressWarnings(as.numeric(map_pos_cm_divisor[[1L]]))
      if (is.finite(legacy_divisor) && legacy_divisor > 0 && !identical(legacy_divisor, 1)) {
        warning("map_pos_cm_divisor is deprecated; use bp_per_cm", call. = FALSE)
        bp_per_cm <- legacy_divisor
      } else {
        ng_stop("bp_per_cm must be supplied explicitly when map_position_unit = 'bp'; ",
                "physical base-pair positions cannot be treated as centimorgans")
      }
    }
    bp_per_cm <- suppressWarnings(as.numeric(bp_per_cm[[1L]]))
    if (!is.finite(bp_per_cm) || bp_per_cm <= 0) ng_stop("bp_per_cm must be a positive number")
    if (is.null(map_pos_bp_col)) map_pos_bp_col <- map_pos_col
    if (is.null(map_pos_bp_col)) {
      pos_hit <- match(
        tolower(c("pos_bp", "position_bp", "Position_BP", "bp", "physical_pos", "position", "pos")),
        tolower(names(marker_map)),
        nomatch = 0L
      )
      pos_hit <- pos_hit[pos_hit > 0L]
      if (length(pos_hit)) map_pos_bp_col <- names(marker_map)[[pos_hit[[1L]]]]
    }
    if (is.null(map_pos_bp_col) || !(map_pos_bp_col %in% names(marker_map))) {
      ng_stop("marker_map needs map_pos_bp_col, or map_pos_col with map_position_unit = 'bp'")
    }
    pos_bp <- suppressWarnings(as.numeric(marker_map[[map_pos_bp_col]]))
    if (any(!is.finite(pos_bp) | pos_bp < 0, na.rm = TRUE)) {
      ng_stop("marker_map base-pair positions must be finite non-negative numbers")
    }
    pos_cm <- pos_bp / bp_per_cm
    position_col <- map_pos_bp_col
  } else {
    bp_per_cm <- NA_real_
    if (is.null(map_pos_cm_col)) map_pos_cm_col <- map_pos_col
    if (is.null(map_pos_cm_col)) {
      pos_hit <- match(
        tolower(c("pos_cm", "cm", "cM", "genetic_pos", "position_cm", "pos", "position")),
        tolower(names(marker_map)),
        nomatch = 0L
      )
      pos_hit <- pos_hit[pos_hit > 0L]
      if (length(pos_hit)) map_pos_cm_col <- names(marker_map)[[pos_hit[[1L]]]]
    }
    if (is.null(map_pos_cm_col) || !(map_pos_cm_col %in% names(marker_map))) {
      ng_stop("marker_map needs map_pos_cm_col, or map_pos_col with map_position_unit = 'cM'")
    }
    pos_cm <- suppressWarnings(as.numeric(marker_map[[map_pos_cm_col]]))
    if (any(!is.finite(pos_cm) | pos_cm < 0, na.rm = TRUE)) {
      ng_stop("marker_map cM positions must be finite non-negative numbers")
    }
    pos_bp <- if (!is.null(map_pos_bp_col)) {
      if (!(map_pos_bp_col %in% names(marker_map))) ng_stop("marker_map is missing map_pos_bp_col: ", map_pos_bp_col)
      suppressWarnings(as.numeric(marker_map[[map_pos_bp_col]]))
    } else {
      rep(NA_real_, length(pos_cm))
    }
    position_col <- map_pos_cm_col
  }

  out <- data.frame(
    marker = trimws(as.character(marker_map[[map_marker_col]])),
    chr = if (is.null(map_chr_col)) 1L else marker_map[[map_chr_col]],
    pos_bp = pos_bp,
    pos_cm = pos_cm,
    stringsAsFactors = FALSE
  )
  if (any(!nzchar(out$marker) | is.na(out$marker))) ng_stop("marker_map contains missing marker IDs")
  attr(out, "map_marker_col") <- map_marker_col
  attr(out, "map_chr_col") <- map_chr_col
  attr(out, "map_position_col") <- position_col
  attr(out, "map_position_unit") <- map_position_unit
  attr(out, "bp_per_cm") <- bp_per_cm
  out
}

ng_run_cp_align_marker_map <- function(geno, marker_map) {
  marker_map_attrs <- attributes(marker_map)[
    c("map_marker_col", "map_chr_col", "map_position_col", "map_position_unit", "bp_per_cm")
  ]
  marker_map <- as.data.frame(marker_map, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("marker" %in% names(marker_map))) ng_stop("marker map must contain an internal marker column")
  if (!("pos_cm" %in% names(marker_map))) ng_stop("marker map must contain internal pos_cm values")
  if (!("pos_bp" %in% names(marker_map))) ng_stop("marker map must contain internal pos_bp values")
  marker_map$marker <- trimws(as.character(marker_map$marker))
  dup_map <- unique(marker_map$marker[duplicated(marker_map$marker) & nzchar(marker_map$marker)])
  if (length(dup_map)) ng_stop("marker map contains duplicate marker(s): ", paste(dup_map, collapse = ", "))
  geno_markers <- colnames(geno)
  missing_markers <- setdiff(geno_markers, marker_map$marker)
  if (length(missing_markers)) {
    ng_stop(
      "marker map is missing ", length(missing_markers),
      " genotype marker(s): ", paste(utils::head(missing_markers, 20L), collapse = ", ")
    )
  }
  extra_markers <- setdiff(marker_map$marker, geno_markers)
  if (length(extra_markers)) {
    ng_stop(
      "marker map contains ", length(extra_markers),
      " marker(s) not present in genotype: ", paste(utils::head(extra_markers, 20L), collapse = ", ")
    )
  }
  out <- marker_map[match(geno_markers, marker_map$marker), , drop = FALSE]
  rownames(out) <- NULL
  for (nm in names(marker_map_attrs)) attr(out, nm) <- marker_map_attrs[[nm]]
  out
}

ng_run_cp_align_parent_tables <- function(geno, phenotype) {
  geno_ids <- rownames(geno)
  pheno_ids <- rownames(phenotype)
  missing_in_genotype <- setdiff(pheno_ids, geno_ids)
  missing_in_phenotype <- setdiff(geno_ids, pheno_ids)
  if (length(missing_in_genotype) || length(missing_in_phenotype)) {
    ng_stop(
      "genotype and phenotype parent IDs must match exactly after QC; phenotype-only=",
      length(missing_in_genotype), ", genotype-only=", length(missing_in_phenotype)
    )
  }
  phenotype <- phenotype[match(geno_ids, pheno_ids), , drop = FALSE]
  rownames(phenotype) <- geno_ids
  list(geno = geno, phenotype = phenotype)
}

ng_run_cp_direction_table <- function(direction,
                                      direction_trait_col = NULL,
                                      direction_column_col = NULL,
                                      direction_direction_col = NULL) {
  direction <- as.data.frame(direction, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(direction_trait_col)) {
    direction_trait_col <- ng_run_cp_first_col(direction, c("trait", "trait_name", "Trait", "TraitName", "name"), "direction_file")
  }
  if (!(direction_trait_col %in% names(direction))) {
    ng_stop("direction_file is missing direction_trait_col: ", direction_trait_col)
  }
  if (is.null(direction_column_col)) {
    column_hit <- match(
      tolower(c("column", "phenotype_column", "PhenotypeColumn", "pheno_col", "trait_column")),
      tolower(names(direction)),
      nomatch = 0L
    )
    column_hit <- column_hit[column_hit > 0L]
    direction_column_col <- if (length(column_hit)) names(direction)[[column_hit[[1L]]]] else direction_trait_col
  }
  if (!(direction_column_col %in% names(direction))) {
    ng_stop("direction_file is missing direction_column_col: ", direction_column_col)
  }
  if (is.null(direction_direction_col)) {
    direction_direction_col <- ng_run_cp_first_col(
      direction,
      c("direction", "trait_direction", "Selection_direction", "SelectionDirection", "selection_direction"),
      "direction_file"
    )
  }
  if (!(direction_direction_col %in% names(direction))) {
    ng_stop("direction_file is missing direction_direction_col: ", direction_direction_col)
  }

  other_cols <- setdiff(names(direction), c(direction_trait_col, direction_column_col, direction_direction_col))
  out <- data.frame(
    trait = direction[[direction_trait_col]],
    column = direction[[direction_column_col]],
    direction = direction[[direction_direction_col]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (length(other_cols)) out <- cbind(out, direction[, other_cols, drop = FALSE])
  attr(out, "direction_columns") <- list(
    trait = direction_trait_col,
    column = direction_column_col,
    direction = direction_direction_col
  )
  out
}

ng_run_cp_direction <- function(direction) {
  ng_multitrait_direction(direction)
}

ng_run_cp_target <- function(progeny) {
  token <- trimws(tolower(as.character(progeny[[1L]])))
  if (token %in% c("dh", "dhs", "doubled_haploid", "doubled_haploids")) return("DH")
  if (token %in% c("ril", "rils")) return("RIL")
  ng_stop("progeny must be DH, DHs, RIL, or RILs")
}

ng_run_cp_optimizer <- function(optimizer) {
  token <- trimws(tolower(as.character(optimizer[[1L]])))
  aliases <- c(
    lp = "mip_linear",
    mip = "mip_contribution",
    ocs = "mip_contribution",
    ocs_qp = "mip_contribution",
    egsi = "greedy_local",
    ga = "evolution",
    de = "evolution",
    memetic = "evolution"
  )
  if (token %in% names(aliases)) token <- aliases[[token]]
  allowed <- c("auto", "greedy_local", "repair_local", "mip_linear", "mip_contribution", "evolution")
  if (!(token %in% allowed)) ng_stop("optimizer must be one of: ", paste(c(allowed, names(aliases)), collapse = ", "))
  token
}

ng_run_cp_allocation_method <- function(allocation_method) {
  token <- trimws(tolower(as.character(allocation_method[[1L]])))
  aliases <- c(
    ocs = "ocs",
    package_ocs = "ocs",
    nextgen_ocs = "ocs",
    alphamate_style = "alphamate_style",
    alphamate_native = "alphamate_style",
    alpha_style = "alphamate_style",
    alphamate_executable = "alphamate_executable",
    alphamate_exact = "alphamate_executable",
    alpha_executable = "alphamate_executable"
  )
  if (!(token %in% names(aliases))) {
    ng_stop("allocation_method must be one of: ", paste(names(aliases), collapse = ", "))
  }
  aliases[[token]]
}

ng_run_cp_method_varPMV <- function(method_varPMV) {
  token <- trimws(tolower(as.character(method_varPMV[[1L]])))
  aliases <- c(
    fast = "fast",
    diagonal = "fast",
    diagonal_posterior = "fast",
    full = "full_posterior",
    full_posterior = "full_posterior",
    full_pmv = "full_posterior"
  )
  if (!(token %in% names(aliases))) {
    ng_stop("method_varPMV must be one of: fast, full_posterior")
  }
  aliases[[token]]
}

ng_run_cp_ril_mode <- function(ril_mode) {
  token <- trimws(tolower(as.character(ril_mode[[1L]])))
  aliases <- c(
    infinite = "infinite",
    infinite_selfing = "infinite",
    fixed = "infinite"
  )
  if (!(token %in% names(aliases))) {
    ng_stop("ril_mode currently supports only 'infinite' for progeny = 'RIL'")
  }
  aliases[[token]]
}

ng_run_cp_integer <- function(x, name, min_value = 0L) {
  raw <- suppressWarnings(as.numeric(x))
  if (length(raw) != 1L || !is.finite(raw) || raw < min_value ||
      abs(raw - round(raw)) > 1e-8) {
    ng_stop(name, " must be an integer >= ", min_value)
  }
  as.integer(round(raw))
}

ng_run_cp_logical <- function(x, name) {
  if (length(x) != 1L || is.na(x)) ng_stop(name, " must be TRUE or FALSE")
  isTRUE(x)
}

ng_run_cp_pmv_col <- function(scored_trait, method_varPMV = "fast") {
  method_varPMV <- ng_run_cp_method_varPMV(method_varPMV)
  if (identical(method_varPMV, "fast")) return("pmv")
  col <- "pmv_full_posterior"
  if (!(col %in% names(scored_trait))) {
    ng_stop("method_varPMV = 'full_posterior' requires score column: ", col)
  }
  vals <- suppressWarnings(as.numeric(scored_trait[[col]]))
  if (!any(is.finite(vals))) {
    ng_stop("method_varPMV = 'full_posterior' produced no finite full-posterior PMV values")
  }
  col
}

ng_run_cp_parallel_cores <- function(n_threads, n_jobs) {
  if (!is.null(n_threads)) {
    cores <- ng_run_cp_integer(n_threads, "n_threads", min_value = 1L)
  } else {
    detected <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    cores <- if (is.finite(detected) && detected > 1L) detected - 1L else 1L
  }
  as.integer(max(1L, min(cores, n_jobs)))
}

ng_run_cp_apply <- function(jobs, fun, use_parallel = FALSE, n_threads = NULL) {
  if (!isTRUE(use_parallel) || length(jobs) < 2L) {
    out <- lapply(jobs, fun)
    attr(out, "parallel_backend") <- "serial"
    attr(out, "n_threads") <- 1L
    return(out)
  }
  cores <- ng_run_cp_parallel_cores(n_threads, length(jobs))
  if (cores < 2L) {
    out <- lapply(jobs, fun)
    attr(out, "parallel_backend") <- "serial"
    attr(out, "n_threads") <- 1L
    return(out)
  }
  if (!identical(.Platform$OS.type, "windows")) {
    out <- parallel::mclapply(jobs, fun, mc.cores = cores)
    attr(out, "parallel_backend") <- "mclapply"
    attr(out, "n_threads") <- cores
    return(out)
  }
  out <- tryCatch({
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, jobs, fun)
  }, error = function(e) {
    warning(
      "use_parallel = TRUE could not start the Windows PSOCK trait workers; ",
      "falling back to serial trait processing. Reason: ", conditionMessage(e),
      call. = FALSE
    )
    NULL
  })
  if (is.null(out)) {
    out <- lapply(jobs, fun)
    attr(out, "parallel_backend") <- "serial_after_parallel_error"
    attr(out, "n_threads") <- 1L
  } else {
    attr(out, "parallel_backend") <- "psock"
    attr(out, "n_threads") <- cores
  }
  out
}

ng_run_cp_variance_col <- function(trait_value_metric,
                                   uc_variance_source,
                                   scored_trait = NULL,
                                   method_varPMV = "fast") {
  metric <- trimws(tolower(as.character(trait_value_metric[[1L]])))
  source <- trimws(tolower(as.character(uc_variance_source[[1L]])))
  if (identical(metric, "usefulness")) metric <- source
  if (metric %in% c("parent_distance", "le") &&
      identical(trimws(tolower(as.character(trait_value_metric[[1L]]))), "usefulness")) {
    ng_stop("uc_variance_source cannot be parent_distance: genomic distance is not a trait variance")
  }
  if (identical(metric, "pmv")) {
    if (is.null(scored_trait)) return("pmv")
    return(ng_run_cp_pmv_col(scored_trait, method_varPMV))
  }
  if (identical(metric, "vpm")) return("vpm")
  # "parent_distance" is the canonical name for the parental genomic-distance
  # variance; "le" is kept as a deprecated alias for callers/results from <= 0.9.0.
  if (metric %in% c("parent_distance", "le")) return("parent_distance")
  if (identical(metric, "var_complex")) {
    if (is.null(scored_trait)) return("pmv")
    return(ng_run_cp_var_complex_col(scored_trait, method_varPMV))
  }
  if (identical(metric, "mean")) return(NA_character_)
  ng_stop("trait_value_metric must be one of: usefulness, pmv, vpm, parent_distance, var_complex, mean")
}

ng_run_cp_var_complex_col <- function(scored_trait, method_varPMV = "fast") {
  pmv_col <- ng_run_cp_pmv_col(scored_trait, method_varPMV)
  candidates <- c(pmv_col, "vpm")
  hit <- candidates[candidates %in% names(scored_trait)]
  if (!length(hit)) ng_stop("scored trait table is missing a usable var_complex variance column")
  hit[[1L]]
}

ng_run_cp_trait_value <- function(scored_trait,
                                  direction,
                                  trait_value_metric = "usefulness",
                                  uc_variance_source = "pmv",
                                  selection_prop = 0.10,
                                  method_varPMV = "fast") {
  direction <- ng_run_cp_direction(direction)[[1L]]
  metric <- trimws(tolower(as.character(trait_value_metric[[1L]])))
  mean_value <- suppressWarnings(as.numeric(scored_trait$cross_mean_blend))
  if (identical(metric, "mean")) return(mean_value)
  if (metric %in% c("parent_distance", "le")) {
    if (!("parent_distance" %in% names(scored_trait))) {
      ng_stop("scored trait table is missing parent_distance")
    }
    return(suppressWarnings(as.numeric(scored_trait$parent_distance)))
  }
  # Pure-variance objectives rank crosses by the predicted family variance
  # itself (PopVar-style), not by mean + i*SD. The usefulness path
  # (trait_value_metric = "usefulness") still combines this variance with the
  # mean. Variance magnitude is direction-agnostic, so no sign is applied.
  if (metric %in% c("vpm", "pmv")) {
    var_col <- ng_run_cp_variance_col(metric, uc_variance_source, scored_trait, method_varPMV)
    if (!(var_col %in% names(scored_trait))) {
      ng_stop("scored trait table is missing variance column: ", var_col)
    }
    return(pmax(suppressWarnings(as.numeric(scored_trait[[var_col]])), 0))
  }
  var_col <- if (identical(metric, "var_complex")) {
    ng_run_cp_var_complex_col(scored_trait, method_varPMV)
  } else {
    ng_run_cp_variance_col(metric, uc_variance_source, scored_trait, method_varPMV)
  }
  if (!(var_col %in% names(scored_trait))) ng_stop("scored trait table is missing variance column: ", var_col)
  sign <- if (identical(direction, "maximize")) 1 else -1
  mean_value + sign * ng_selection_intensity(selection_prop) * sqrt(pmax(suppressWarnings(as.numeric(scored_trait[[var_col]])), 0))
}

ng_run_cp_clean_trait_name <- function(x) {
  out <- make.names(as.character(x))
  out <- gsub("[.]+", "_", out)
  out <- gsub("^_|_$", "", out)
  ifelse(nzchar(out), out, "trait")
}

ng_run_cp_trait_spec <- function(direction,
                                 traits_to_use = NULL,
                                 trait_weights = NULL) {
  direction <- as.data.frame(direction, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("trait" %in% names(direction)) && "column" %in% names(direction)) direction$trait <- direction$column
  if (!("column" %in% names(direction)) && "trait" %in% names(direction)) direction$column <- direction$trait
  if (!("direction" %in% names(direction))) ng_stop("direction_file must contain a direction column")
  required <- c("trait", "column", "direction")
  missing <- setdiff(required, names(direction))
  if (length(missing)) ng_stop("direction_file missing columns: ", paste(missing, collapse = ", "))
  out <- direction
  out$trait <- trimws(as.character(out$trait))
  out$column <- trimws(as.character(out$column))
  out$direction <- ng_run_cp_direction(out$direction)
  if (!is.null(traits_to_use)) {
    keep <- out$trait %in% traits_to_use | out$column %in% traits_to_use
    out <- out[keep, , drop = FALSE]
  }
  if (!nrow(out)) ng_stop("No traits remain after applying traits_to_use")
  if (!("weight" %in% names(out))) out$weight <- NA_real_
  if (!is.null(trait_weights)) {
    tw <- suppressWarnings(as.numeric(trait_weights))
    if (!is.null(names(trait_weights)) && any(nzchar(names(trait_weights)))) {
      names(tw) <- names(trait_weights)
      mapped <- tw[out$trait]
      missing <- !is.finite(mapped)
      mapped[missing] <- tw[out$column[missing]]
      out$weight <- as.numeric(mapped)
    } else {
      if (length(tw) != nrow(out)) ng_stop("unnamed trait_weights must match the number of selected traits")
      out$weight <- tw
    }
  }
  out
}

ng_run_cp_index_spec <- function(index_col, index_direction) {
  if (is.null(index_col) || !nzchar(as.character(index_col[[1L]]))) {
    ng_stop("index_col is required when prediction_mode = 'index_as_trait'")
  }
  data.frame(
    trait = "selection_index",
    column = as.character(index_col[[1L]]),
    direction = ng_run_cp_direction(index_direction)[[1L]],
    weight = NA_real_,
    stringsAsFactors = FALSE
  )
}

ng_run_cp_parent_use <- function(plan) {
  if (!nrow(plan)) {
    return(data.frame(parent = character(), crosses = integer(), stringsAsFactors = FALSE))
  }
  parents <- c(as.character(plan$parent1), as.character(plan$parent2))
  tab <- sort(table(parents), decreasing = TRUE)
  data.frame(parent = names(tab), crosses = as.integer(tab), stringsAsFactors = FALSE)
}

ng_run_cp_output_files <- function(output_dir,
                                   output_file,
                                   selected_crosses,
                                   candidate_crosses,
                                   trait_spec,
                                   qc,
                                   write_outputs,
                                   write_figures,
                                   n_crosses,
                                   include_trait_gebv = FALSE,
                                   trait_check_reference = NULL,
                                   multi_trait_meta = NULL,
                                   trait_value_metric = NULL) {
  files <- list()
  if (!isTRUE(write_outputs) && !isTRUE(write_figures)) return(files)
  if (is.null(output_dir) || !nzchar(as.character(output_dir[[1L]]))) {
    ng_stop("output_dir is required when write_outputs or write_figures is TRUE")
  }
  output_dir <- as.character(output_dir[[1L]])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figures <- NULL
  if (isTRUE(write_figures)) {
    chk_line <- if (is.null(trait_check_reference)) NA_real_ else
      ng_check_line_value(trait_check_reference, multi_trait_meta,
                          candidate_scores = candidate_crosses,
                          trait_value_metric = trait_value_metric)
    plot_path <- file.path(output_dir, "priority_score_vs_kinship.png")
    ng_plot_priority_score_vs_kinship(
      scored = candidate_crosses,
      selected = selected_crosses,
      output_path = plot_path,
      check_line = if (is.finite(chk_line)) chk_line else NULL,
      check_label = if (is.finite(chk_line)) ng_check_line_label(trait_check_reference$active) else NULL
    )
    files$priority_score_vs_kinship_png <- normalizePath(plot_path, winslash = "/", mustWork = TRUE)
    figures <- data.frame(
      figure = "priority_score_vs_kinship",
      path = files$priority_score_vs_kinship_png,
      stringsAsFactors = FALSE
    )
    if (!is.null(trait_check_reference) && nrow(trait_check_reference$active) > 1L) {
      panel_path <- file.path(output_dir, "check_panels.png")
      ng_plot_check_panels(candidate_crosses, trait_check_reference, output_path = panel_path)
      files$check_panels_png <- normalizePath(panel_path, winslash = "/", mustWork = TRUE)
      figures <- rbind(figures, data.frame(figure = "check_panels",
                                           path = files$check_panels_png,
                                           stringsAsFactors = FALSE))
    }
  }
  if (isTRUE(write_outputs)) {
    workbook_path <- file.path(output_dir, if (is.null(output_file)) "crossing_plan.xlsx" else as.character(output_file[[1L]]))
    files$workbook <- ng_write_cross_priority_workbook(
      output_path = workbook_path,
      crosses = selected_crosses,
      scored = candidate_crosses,
      trait_directions = trait_spec,
      parent_use = ng_run_cp_parent_use(selected_crosses),
      duplicate_pairs = if (!is.null(qc$putative_duplicates)) qc$putative_duplicates$pairs else NULL,
      figures = figures,
      n_crosses_requested = n_crosses,
      include_trait_gebv = isTRUE(include_trait_gebv),
      trait_check_reference = trait_check_reference
    )
  }
  files
}


# ---------------------------------------------------------------------------
# Staged pipeline: ng_run_cross_prediction is decomposed into internal stage
# functions over a shared `ctx` list. The public runner is a thin driver that
# builds ctx and calls the stages in order; a later resumable ng_run_stage
# calls the SAME stage functions so staged execution is byte-identical.
# ---------------------------------------------------------------------------

ng_cp_stage_order <- function() c("qc", "predict", "index", "allocate", "rank")

# Map breeder-facing metric names to the canonical tokens the pipeline uses.
# Friendly names are additive aliases; canonical tokens keep working unchanged.
ng_normalize_metric_token <- function(x) {
  token <- trimws(tolower(as.character(x[[1L]])))
  switch(token,
    mid_parent_mean          = "mean",
    family_variance          = "vpm",
    reliable_family_variance = "pmv",
    token)
}

ng_cp__build_ctx <- function(config) {
  ctx <- config
  ctx$prediction_mode <- match.arg(ctx$prediction_mode, c("trait_by_trait", "index_as_trait"))
  ctx$map_position_unit <- match.arg(ctx$map_position_unit, c("bp", "cM"))
  # Accept breeder-facing metric names as aliases of the canonical tokens.
  # var_complex is a deprecated alias for usefulness + reliable (pmv) variance,
  # so it maps to BOTH fields to preserve its exact (UC-PMV) behavior.
  # Preserve the user's ORIGINAL choices so result$settings echoes what they
  # supplied (e.g. "var_complex"/"family_variance"), not the normalized token.
  ctx$trait_value_metric_input <- ctx$trait_value_metric
  ctx$uc_variance_source_input <- ctx$uc_variance_source
  .tv <- ng_normalize_metric_token(ctx$trait_value_metric)
  if (identical(.tv, "var_complex")) {
    ctx$trait_value_metric <- "usefulness"
    ctx$uc_variance_source <- "pmv"
  } else {
    ctx$trait_value_metric <- .tv
    ctx$uc_variance_source <- ng_normalize_metric_token(ctx$uc_variance_source)
  }
  ctx$trait_value_metric <- match.arg(ctx$trait_value_metric, c("usefulness", "pmv", "vpm", "parent_distance", "le", "var_complex", "mean"))
  ctx$uc_variance_source <- match.arg(ctx$uc_variance_source, c("pmv", "vpm", "parent_distance", "le"))
  ctx$threshold_policy <- match.arg(ctx$threshold_policy, c("soft", "strict"))
  ctx$recomb_model <- match.arg(ctx$recomb_model, c("haldane", "kosambi"))
  ctx$grm_method <- match.arg(ctx$grm_method, c("vanraden", "yang"))
  ctx$ld_backend <- match.arg(ctx$ld_backend, c("auto", "cpp", "r"))
  ctx$method_varPMV <- ng_run_cp_method_varPMV(ctx$method_varPMV)
  ctx$ril_mode <- ng_run_cp_ril_mode(ctx$ril_mode)
  ctx$posterior_method <- match.arg(ctx$posterior_method, c("mcmc", "closed_form"))
  ctx$duplicate_action <- match.arg(ctx$duplicate_action, c("remove", "report", "none"))
  ctx$lambda_parent_use_mode <- match.arg(ctx$lambda_parent_use_mode, c("absolute", "adaptive"))
  ctx$alphamate_mode <- match.arg(ctx$alphamate_mode, c("ModeOptTarget1", "ModeMaxCriterion", "ModeMinCoancestry"))
  ctx$target <- ng_run_cp_target(ctx$progeny)
  ctx$optimizer_method <- ng_run_cp_optimizer(ctx$optimizer)
  ctx$allocation_method <- ng_run_cp_allocation_method(ctx$allocation_method)
  ctx$run_posterior_prediction <- ng_run_cp_logical(ctx$run_posterior_prediction, "run_posterior_prediction")
  ctx$use_parallel <- ng_run_cp_logical(ctx$use_parallel, "use_parallel")
  ctx$n_iter <- ng_run_cp_integer(ctx$n_iter, "n_iter", min_value = 1L)
  ctx$burn_in <- ng_run_cp_integer(ctx$burn_in, "burn_in", min_value = 0L)
  ctx$posterior_n_draws <- max(1L, ctx$n_iter - ctx$burn_in)
  ctx$n_crosses <- ng_run_cp_integer(ctx$n_crosses, "n_crosses", min_value = 1L)
  ctx
}

ng_cp__stage_qc <- function(ctx) {
  list2env(ctx, environment())
  direction_canonical <- ctx$direction_canonical
  phenotype_raw <- if (is.null(phenotype)) ng_run_cp_read_csv(phenotype_file, "phenotype") else phenotype
  genotype_raw <- if (is.null(genotype)) ng_run_cp_read_csv(genotype_file, "genotype") else genotype
  map_raw <- if (is.null(marker_map)) ng_run_cp_read_csv(map_file, "marker_map") else marker_map
  direction_raw <- if (is.null(trait_direction)) ng_run_cp_read_csv(direction_file, "direction") else trait_direction
  if (is.null(phenotype_raw)) ng_stop("phenotype_file or phenotype is required")
  if (is.null(genotype_raw)) ng_stop("genotype_file or genotype is required")
  if (is.null(map_raw)) ng_stop("map_file or marker_map is required")

  if (is.null(phenotype_id_col)) phenotype_id_col <- id_col
  if (is.null(genotype_id_col)) genotype_id_col <- id_col
  phenotype_input <- ng_run_cp_canonical_id_table(phenotype_raw, phenotype_id_col, "phenotype")
  genotype_input <- ng_run_cp_canonical_id_table(genotype_raw, genotype_id_col, "genotype")
  phenotype_id_col_used <- attr(phenotype_input, "ng_original_id_col")
  genotype_id_col_used <- attr(genotype_input, "ng_original_id_col")

  trait_spec <- if (identical(prediction_mode, "trait_by_trait")) {
    if (is.null(direction_raw)) ng_stop("direction_file or trait_direction is required when prediction_mode = 'trait_by_trait'")
    direction_canonical <- ng_run_cp_direction_table(
      direction_raw,
      direction_trait_col = direction_trait_col,
      direction_column_col = direction_column_col,
      direction_direction_col = direction_direction_col
    )
    direction_columns <- attr(direction_canonical, "direction_columns")
    ng_run_cp_trait_spec(direction_canonical, traits_to_use = traits_to_use, trait_weights = trait_weights)
  } else {
    direction_columns <- list(
      trait = NA_character_,
      column = as.character(index_col[[1L]]),
      direction = "index_direction"
    )
    ng_run_cp_index_spec(index_col, index_direction)
  }

  marker_map_std <- ng_run_cp_marker_map(
    map_raw,
    map_marker_col = map_marker_col,
    map_chr_col = map_chr_col,
    map_pos_col = map_pos_col,
    map_pos_cm_col = map_pos_cm_col,
    map_pos_bp_col = map_pos_bp_col,
    map_position_unit = map_position_unit,
    bp_per_cm = bp_per_cm,
    map_pos_cm_divisor = map_pos_cm_divisor
  )

  qc <- ng_preflight_input_tables(
    geno = genotype_input,
    phenotype = phenotype_input,
    trait_spec = trait_spec,
    marker_map = marker_map_std,
    ploidy = 2L,
    putative_duplicate_check = !identical(duplicate_action, "none"),
    putative_duplicate_action = if (identical(duplicate_action, "none")) "report" else duplicate_action,
    duplicate_threshold = duplicate_threshold,
    duplicate_maf_min = duplicate_maf_min,
    duplicate_max_missing_prop = duplicate_max_missing_prop,
    duplicate_min_compared_markers = duplicate_min_compared_markers
  )
  # NOTE: on a blocker, the ORIGINAL monolith threw right here, before the
  # clean/align block below ever ran. The staged path must not throw (the
  # frontend needs to persist qc.json with status="blocker" + issues), so we
  # set ctx$qc and return early instead -- clean/align is still skipped on a
  # blocker, exactly as before. The one-shot driver (ng_run_cross_prediction)
  # re-introduces the throw, with this exact message, immediately after this
  # stage returns -- see the stage loop below.
  ctx <- ng_ctx_put(ctx, qc = qc)
  if (any(qc$issues$severity == "blocker")) {
    return(ctx)
  }

  cleaned <- if (!is.null(qc$cleaned_tables)) qc$cleaned_tables else list(
    geno = genotype_input,
    phenotype = phenotype_input,
    marker_map = marker_map_std,
    trait_spec = trait_spec
  )
  geno <- ng_run_cp_geno_matrix(cleaned$geno, id_col = "parent")
  pheno <- ng_run_cp_pheno_frame(cleaned$phenotype, id_col = "parent")
  aligned <- ng_run_cp_align_parent_tables(geno, pheno)
  geno <- aligned$geno
  pheno <- aligned$phenotype
  ids <- rownames(geno)
  if (length(ids) < 2L) ng_stop("At least two matched parents are required after QC")
  missing_trait_cols <- setdiff(trait_spec$column, names(pheno))
  if (length(missing_trait_cols)) {
    ng_stop(
      "phenotype is missing trait/index column(s) from direction file: ",
      paste(missing_trait_cols, collapse = ", ")
    )
  }
  marker_map_std <- ng_run_cp_align_marker_map(geno, cleaned$marker_map)
  ctx <- ng_ctx_put(
    ctx,
    phenotype_id_col_used = phenotype_id_col_used,
    genotype_id_col_used = genotype_id_col_used,
    trait_spec = trait_spec,
    direction_columns = direction_columns,
    direction_canonical = direction_canonical,
    marker_map_std = marker_map_std,
    qc = qc,
    geno = geno,
    pheno = pheno,
    ids = ids
  )
  ctx
}

ng_cp__stage_predict <- function(ctx) {
  list2env(ctx, environment())
  # Optional LD pruning of the marker matrix (same capability as ng_design_crosses): drop
  # redundant/low-MAF markers before scoring. Per-trait effects are fit from `geno` inside
  # the loop below, so pruning here keeps effects, geno, and the map aligned automatically.
  ld_pruning_report <- NULL
  if (isTRUE(ld_pruning)) {
    geno <- ng_ld_prune_geno(
      geno = geno, window = ld_window, r2_threshold = ld_r2_threshold,
      maf_threshold = ld_maf_threshold, ploidy = ld_ploidy, backend = ld_backend,
      marker_map = marker_map_std
    )
    ld_pruning_report <- attr(geno, "ld_pruning_report", exact = TRUE)
    # subset the already-aligned map to the retained markers, preserving its column attrs
    map_attrs <- attributes(marker_map_std)[
      c("map_marker_col", "map_chr_col", "map_position_col", "map_position_unit", "bp_per_cm")]
    marker_map_std <- marker_map_std[match(colnames(geno), marker_map_std$marker), , drop = FALSE]
    rownames(marker_map_std) <- NULL
    for (nm in names(map_attrs)) attr(marker_map_std, nm) <- map_attrs[[nm]]
    ids <- rownames(geno)
  }

  # Optional marker-effect TRAINING augmentation: extra individuals ("others") that enlarge the
  # ridge fit but are NOT candidate parents. Aligned to the (possibly LD-pruned) parent markers,
  # and any id that is also a parent is dropped so parents are never double-counted.
  tg_idcol <- if (is.null(training_genotype_id_col)) id_col else training_genotype_id_col
  tp_idcol <- if (is.null(training_phenotype_id_col)) id_col else training_phenotype_id_col
  training_set <- ng_run_cp_training_set(
    training_genotype = training_genotype, training_phenotype = training_phenotype,
    training_genotype_file = training_genotype_file, training_phenotype_file = training_phenotype_file,
    training_genotype_id_col = tg_idcol, training_phenotype_id_col = tp_idcol,
    parent_ids = ids, parent_markers = colnames(geno), trait_columns = trait_spec$column
  )
  if (!is.null(training_set) && trait_value_metric[[1L]] %in% c("parent_distance", "le")) {
    warning("trait_value_metric = 'parent_distance' does not use marker effects; the supplied ",
            "training set will not change the selected metric or the crossing plan.", call. = FALSE)
  }
  training_only_count <- if (is.null(training_set)) 0L else length(training_set$ids)
  training_ids <- if (is.null(training_set)) character(0L) else training_set$ids

  effects_list <- list()
  trait_scores <- list()
  posterior_effects_list <- list()
  posterior_predictions_list <- list()
  cross_table <- NULL
  set.seed(seed)

  run_trait_job <- function(i) {
    trait <- trait_spec$trait[[i]]
    column <- trait_spec$column[[i]]
    if (!(column %in% names(pheno))) ng_stop("phenotype is missing trait/index column: ", column)
    y <- suppressWarnings(as.numeric(pheno[[column]]))
    names(y) <- rownames(pheno)
    if (sum(is.finite(y)) < 2L) ng_stop("trait/index column has fewer than two finite values: ", column)
    # Enlarge the training set with the extra individuals that have a finite phenotype for this
    # trait. Effects are fit on parents + training; crosses are still scored on parents only.
    fit_geno <- geno; fit_y <- y; fit_ids <- ids
    if (!is.null(training_set)) {
      tr_y <- suppressWarnings(as.numeric(training_set$pheno[[column]]))
      names(tr_y) <- training_set$ids
      keep <- is.finite(tr_y)
      if (any(keep)) {
        fit_geno <- rbind(geno, training_set$geno[keep, , drop = FALSE])
        fit_y <- c(y, tr_y[keep])
        fit_ids <- c(ids, training_set$ids[keep])
      }
    }
    n_effect_training <- length(fit_ids)
    # Also request the full beta covariance (used for the mid-parent PEV / cross-priority
    # risk annotation, Task 3) whenever the marker count is tractable, even outside the
    # full_posterior PMV path. Multi-trait runs need it per trait (the index PEV is a
    # weighted sum of the per-trait mid-parent PEVs), so the marker cap -- not the trait
    # count -- is what gates it; each job reduces Sigma_beta to an n_pairs PEV vector and
    # drops the dense matrix below, so peak memory is one Sigma_beta per parallel worker.
    want_beta_cov <- identical(method_varPMV, "full_posterior") || ncol(geno) <= 6000L
    fit <- ng_fit_ridge_effects(
      geno = fit_geno,
      y = fit_y,
      ids = fit_ids,
      seed = seed + i - 1L,
      return_beta_cov_full = want_beta_cov
    )
    posterior_cov_full <- if (identical(method_varPMV, "full_posterior")) fit$beta_cov_full else NULL
    scored_trait <- ng_score_crosses(
      geno = geno,
      effects = fit,
      marker_map = marker_map_std,
      ids = ids,
      adjusted_pheno = y,
      target = target,
      selection_prop = selection_prop,
      min_effect_reliability = min_effect_reliability,
      recomb_model = recomb_model,
      use_cpp = use_cpp,
      parent_type = parent_type,
      phased_haplotypes = phased_haplotypes,
      grm_method = grm_method,
      posterior_cov_full = posterior_cov_full
    )
    value <- ng_run_cp_trait_value(
      scored_trait = scored_trait,
      direction = trait_spec$direction[[i]],
      trait_value_metric = trait_value_metric,
      uc_variance_source = uc_variance_source,
      selection_prop = selection_prop,
      method_varPMV = method_varPMV
    )
    posterior_effects <- NULL
    posterior_scores <- NULL
    if (isTRUE(run_posterior_prediction)) {
      posterior_effects <- ng_fit_ridge_effects_posterior(
        geno = fit_geno,
        y = fit_y,
        ids = fit_ids,
        n_draws = posterior_n_draws,
        method = posterior_method,
        mcmc_burnin = burn_in,
        seed = seed + 1000L + i - 1L
      )
      posterior_scores <- ng_posterior_cross_predict(
        geno = geno,
        posterior_effects = posterior_effects,
        marker_map = marker_map_std,
        ids = ids,
        adjusted_pheno = y,
        target = target,
        selection_prop = selection_prop,
        min_effect_reliability = min_effect_reliability,
        recomb_model = recomb_model,
        use_cpp = use_cpp,
        parent_type = parent_type,
        # Summarize the SELECTED metric, not the hardcoded usefulness column: a posterior
        # interval computed on usefulness must never be presented as the uncertainty of a
        # `mean` / `var_complex` run (design doc F2).
        # method_varPMV = "fast" INSIDE the draw loop is deliberate, not a fallback.
        # "full_posterior" inflates the variance by the marker-effect uncertainty; inside the
        # posterior loop that uncertainty is already carried by the draws themselves (R/30
        # scores each draw with beta_var = 0), so requesting it here would double-count it --
        # and the per-draw tables carry no pmv_full_posterior column to begin with.
        value_fun = function(sc) ng_run_cp_trait_value(
          scored_trait = sc, direction = trait_spec$direction[[i]],
          trait_value_metric = trait_value_metric, uc_variance_source = uc_variance_source,
          selection_prop = selection_prop, method_varPMV = "fast"),
        # include the plan size so prob_top_tier answers "top-N where N is the plan"
        top_n_targets = unique(as.integer(c(n_crosses, min(n_crosses, 10L), 10L, 20L, 50L)))
      )
    }
    clean_trait <- ng_run_cp_clean_trait_name(trait)
    pmv_used_col <- ng_run_cp_pmv_col(scored_trait, method_varPMV)
    var_complex_col <- ng_run_cp_var_complex_col(scored_trait, method_varPMV)
    # Cross-priority risk: reduce the dense Sigma_beta to the per-cross mid-parent PEV here,
    # while it is still in scope, so only an n_pairs vector travels back from the job. The
    # PEV then rides the cross table as `<trait>_midparent_pev` and survives candidate
    # filtering and allocation by construction (rather than being recomputed downstream from
    # a Sigma_beta that, for multi-trait runs, only the last trait would have supplied).
    midparent_pev <- ng_midparent_pev(geno, scored_trait[, c("parent1", "parent2")],
                                      fit$beta_cov_full, fit$marker_mean)
    # Posterior spread of the ranked value + P(top-N) on that same value. Reduced to vectors
    # here so they ride the cross table and survive candidate filtering / allocation, exactly
    # like the mid-parent PEV above.
    post_sd <- NULL; post_topn <- NULL
    if (!is.null(posterior_scores)) {
      if ("ranked_value_post_sd" %in% names(posterior_scores)) {
        post_sd <- suppressWarnings(as.numeric(posterior_scores$ranked_value_post_sd))
      }
      tn <- paste0("posterior_topn_prob_", as.integer(n_crosses))
      if (tn %in% names(posterior_scores)) {
        post_topn <- suppressWarnings(as.numeric(posterior_scores[[tn]]))
      }
    }
    if (!identical(method_varPMV, "full_posterior")) fit$beta_cov_full <- NULL
    list(
      index = i,
      trait = trait,
      column = column,
      clean_trait = clean_trait,
      y = y,
      fit = fit,
      midparent_pev = midparent_pev,
      post_sd = post_sd,
      post_topn = post_topn,
      scored_trait = scored_trait,
      value = value,
      pmv_used_col = pmv_used_col,
      var_complex_col = var_complex_col,
      posterior_effects = posterior_effects,
      posterior_scores = posterior_scores,
      effect_summary = data.frame(
        trait = trait,
        column = column,
        direction = trait_spec$direction[[i]],
        value_column = paste0(clean_trait, "_value"),
        marker_effect_reliability = NA_real_,
        cv_predictive_r2 = fit$cv_predictive_r2,
        marker_effect_training_n = as.integer(n_effect_training),
        ridge_lambda = fit$lambda,
        method_varPMV = method_varPMV,
        pmv_column_used = pmv_used_col,
        posterior_prediction = isTRUE(run_posterior_prediction),
        posterior_method = if (isTRUE(run_posterior_prediction)) posterior_method else NA_character_,
        posterior_draws = if (isTRUE(run_posterior_prediction)) posterior_n_draws else NA_integer_,
        stringsAsFactors = FALSE
      )
    )
  }

  trait_results <- ng_run_cp_apply(
    jobs = seq_len(nrow(trait_spec)),
    fun = run_trait_job,
    use_parallel = use_parallel,
    n_threads = n_threads
  )
  parallel_backend <- attr(trait_results, "parallel_backend")
  parallel_cores_used <- attr(trait_results, "n_threads")
  trait_results <- trait_results[order(vapply(trait_results, `[[`, integer(1L), "index"))]
  effect_summary <- vector("list", length(trait_results))
  trait_mean_source <- list()

  for (j in seq_along(trait_results)) {
    item <- trait_results[[j]]
    trait <- item$trait
    scored_trait <- item$scored_trait
    clean_trait <- item$clean_trait
    if (is.null(cross_table)) {
      pair_cols <- intersect(c("parent1", "parent2", "pair_relationship", "pair_kinship"),
                             names(scored_trait))
      cross_table <- scored_trait[, pair_cols, drop = FALSE]
    }
    cross_table[[paste0(clean_trait, "_value")]] <- item$value
    cross_table[[paste0(clean_trait, "_mean")]] <- scored_trait$cross_mean_blend
    if (!exists("trait_mean_source", inherits = FALSE)) trait_mean_source <- list()
    trait_mean_source[[trait]] <- scored_trait$mean_source[[1L]]
    cross_table[[paste0(clean_trait, "_pmv")]] <- scored_trait$pmv
    cross_table[[paste0(clean_trait, "_pmv_fast")]] <- scored_trait$pmv
    cross_table[[paste0(clean_trait, "_pmv_full_posterior")]] <- scored_trait$pmv_full_posterior
    cross_table[[paste0(clean_trait, "_pmv_used")]] <- scored_trait[[item$pmv_used_col]]
    cross_table[[paste0(clean_trait, "_vpm")]] <- scored_trait$vpm
    cross_table[[paste0(clean_trait, "_mean_gebv")]] <- scored_trait$cross_mean_gebv
    cross_table[[paste0(clean_trait, "_parent_distance")]] <- scored_trait$parent_distance
    cross_table[[paste0(clean_trait, "_var_complex")]] <- scored_trait[[item$var_complex_col]]
    cross_table[[paste0(clean_trait, "_reliability")]] <- scored_trait$effect_reliability
    cross_table[[paste0(clean_trait, "_cv_predictive_r2")]] <- scored_trait$effect_cv_predictive_r2
    # Per-trait mid-parent PEV (cross-priority risk). NA when Sigma_beta was not requested.
    cross_table[[paste0(clean_trait, "_midparent_pev")]] <- item$midparent_pev
    if (!is.null(item$post_sd))   cross_table[[paste0(clean_trait, "_post_sd")]]   <- item$post_sd
    if (!is.null(item$post_topn)) cross_table[[paste0(clean_trait, "_post_topn")]] <- item$post_topn
    effects_list[[trait]] <- item$fit
    trait_scores[[trait]] <- scored_trait
    if (isTRUE(run_posterior_prediction)) {
      posterior_effects_list[[trait]] <- item$posterior_effects
      posterior_predictions_list[[trait]] <- item$posterior_scores
    }
    effect_summary[[j]] <- item$effect_summary
  }

  # EXACT within-family cross-trait covariance (multi-trait only): Cov(t, s | cross) = a_t' R a_s,
  # the recombination-aware generalization of the per-trait a'Ra used for vpm. This is what makes
  # the multi-trait portfolio upside sqrt(w' S w) correct for genetically correlated traits --
  # summing per-trait variances would overstate the index spread whenever the traits are
  # antagonistic within the family (the classic yield/protein case). The wf_var_/wf_cov_ columns
  # ride the cross table so they survive candidate filtering and allocation by row alignment.
  # Cost scales with T (one kernel pass per trait); the T(T+1)/2 trait pairs are cheap colSums.
  # Hoisted to NULL so it is always defined -- including for the check-reference joint
  # probability computed later in ng_cp__stage_index, which falls back to its G_hat proxy
  # when this is NULL (single-trait runs, or nrow(cross_table) == 0).
  ctc <- NULL
  if (nrow(trait_spec) > 1L && !is.null(cross_table) && nrow(cross_table) > 0L) {
    betas <- vapply(trait_spec$trait, function(tr) {
      b <- suppressWarnings(as.numeric(effects_list[[tr]]$beta))
      b[!is.finite(b)] <- 0
      b
    }, numeric(ncol(geno)))
    betas <- matrix(betas, nrow = ncol(geno), ncol = nrow(trait_spec),
                    dimnames = list(colnames(geno), trait_spec$trait))
    # Chunk the PAIR dimension. ng_cross_trait_within_family_cov materializes an n_pairs x m
    # contrast matrix plus one m x n_pairs matrix per trait, so a full candidate pool at
    # production scale (200 parents = 19,900 pairs x 5,000 markers x 3 traits) would allocate
    # several GB in one go. Crosses are row-independent here, so chunking is numerically
    # identical and caps peak memory at roughly chunk x markers x 8 bytes x (1 + 2T).
    all_pairs <- cross_table[, c("parent1", "parent2"), drop = FALSE]
    chunk <- max(200L, as.integer(5e6 / max(ncol(geno), 1L)))
    starts <- seq(1L, nrow(all_pairs), by = chunk)
    ctc <- do.call(rbind, lapply(starts, function(s) {
      rows <- seq(s, min(s + chunk - 1L, nrow(all_pairs)))
      ng_cross_trait_within_family_cov(
        geno = geno, betas = betas, marker_map = marker_map_std, ids = ids,
        pairs = all_pairs[rows, , drop = FALSE],
        target = target, recomb_model = recomb_model)
    }))
    for (cc in setdiff(names(ctc), c("parent1", "parent2"))) cross_table[[cc]] <- ctc[[cc]]
  }

  # Carried into ng_cp__stage_index for the check-reference block: the per-trait mean_source
  # (so a check resolves onto the same scale as the cross means it is compared against) and the
  # exact within-family cross-trait covariance (for the multi-check joint probability).
  ctx <- ng_ctx_put(
    ctx,
    ld_pruning_report = ld_pruning_report,
    geno = geno,
    marker_map_std = marker_map_std,
    ids = ids,
    training_only_count = training_only_count,
    training_ids = training_ids,
    effects_list = effects_list,
    trait_scores = trait_scores,
    posterior_effects_list = posterior_effects_list,
    posterior_predictions_list = posterior_predictions_list,
    cross_table = cross_table,
    trait_mean_source = trait_mean_source,
    ctc = ctc,
    parallel_backend = parallel_backend,
    parallel_cores_used = parallel_cores_used,
    effect_summary = effect_summary
  )
  ctx
}

ng_cp__stage_index <- function(ctx) {
  list2env(ctx, environment())
  direction_canonical <- ctx$direction_canonical
  # ctc (exact within-family cross-trait covariance, R/39 stage_predict) is deliberately assigned
  # NULL rather than omitted when it cannot be computed (single-trait runs, or an empty candidate
  # pool) -- ng_attach_joint_check_probability() below is written to treat NULL as "fall back to
  # the G_hat proxy". Plain `ctx$ctc <- ctc` would silently DROP the key when ctc is NULL
  # (assigning NULL into a list element deletes it), so list2env() above would never bind a bare
  # `ctc` in that case; the upstream write now goes through ng_ctx_put() (R/00_utils.R), which
  # keeps the key. Read it explicitly via `ctx$` anyway (returns NULL for an absent key, and costs
  # nothing when the key is present) so this path never regresses to "object 'ctc' not found" if a
  # future write here is ever done the plain way again.
  ctc <- ctx$ctc
  # check_pheno is NULL on almost every run and is never reassigned on ctx by any earlier
  # stage, so list2env() above already binds it correctly -- but it is read via `ctx$` anyway,
  # the same defensive pattern as `ctc` just above, since this is exactly the shape of field
  # ("usually NULL") that has silently gone unbound here before.
  check_pheno <- ctx$check_pheno
  # Per-cross COST / LOGISTICS: candidate crosses are generated internally, so a breeder who
  # wants cost/budget/logistic-aware allocation supplies cross_cost -- a data frame with
  # parent1/parent2 + one numeric column per factor (e.g. cost, distance). Joined onto the
  # candidate crosses by UNORDERED parent pair; cost_col / logistic_col then reference those
  # columns. Without this the runner cannot attach per-cross cost, so the knobs would no-op.
  if (!is.null(cross_cost)) {
    cc <- as.data.frame(cross_cost, stringsAsFactors = FALSE)
    if (!all(c("parent1", "parent2") %in% names(cc))) {
      ng_stop("cross_cost must have parent1 and parent2 columns")
    }
    val_cols <- setdiff(names(cc), c("parent1", "parent2"))
    if (!length(val_cols)) {
      ng_stop("cross_cost must have at least one value column (e.g. cost, distance)")
    }
    key_ct <- ng_group_pair_key(cross_table$parent1, cross_table$parent2)
    key_cc <- ng_group_pair_key(as.character(cc$parent1), as.character(cc$parent2))
    idx <- match(key_ct, key_cc)
    for (vc in val_cols) cross_table[[vc]] <- suppressWarnings(as.numeric(cc[[vc]][idx]))
    n_unmatched <- sum(is.na(idx))
    if (n_unmatched) {
      message(sprintf("cross_cost: %d of %d candidate crosses had no matching cost row (left as NA).",
                      n_unmatched, nrow(cross_table)))
    }
  }
  # A cost/logistic column requested for allocation MUST exist on the candidate crosses --
  # otherwise the constraint would be silently ignored. Fail loudly with a fix hint.
  for (req_col in c(cost_col, logistic_col)) {
    if (!is.null(req_col) && !(req_col %in% names(cross_table))) {
      ng_stop("cost_col/logistic_col '", req_col, "' is not a column on the candidate crosses; ",
              "supply it via cross_cost = data.frame(parent1, parent2, ", req_col, " = ...).")
    }
  }

  objective_traits <- trait_spec
  objective_traits$column <- vapply(objective_traits$trait, function(trait) {
    paste0(ng_run_cp_clean_trait_name(trait), "_value")
  }, character(1L))
  objective <- ng_breeder_selection_objective(
    trait = objective_traits,
    method = multi_trait_method,
    threshold_policy = threshold_policy
  )

  # Lethal-allele guarding (Module 4): drop carrier x carrier matings at the nominated
  # deleterious recessive loci before scoring/allocation, so they never enter the plan.
  # Record the candidate count before filtering so the number dropped is recoverable for
  # the run diagnostics (ng_apply_marker_management returns the filtered frame only).
  n_candidates_pre_lethal <- nrow(cross_table)
  if (!is.null(lethal_spec)) {
    cross_table <- ng_apply_marker_management(
      scores = cross_table, geno = geno, lethal_spec = lethal_spec,
      ploidy = marker_ploidy, drop_lethal_carrier_crosses = drop_lethal_carrier_crosses)
  }

  # Check lines as REFERENCES (Module: check reference). Attaches per-trait reference columns
  # and never touches the candidate set: no filtering, no penalty, no reordering. Checks are
  # kept out of geno entirely, so QC, LD, the GRM, and ng_make_pairs() are all untouched and a
  # run with checks is numerically identical to the same run without them.
  trait_check_reference <- NULL
  if (!is.null(trait_checks)) {
    if (is.null(check_geno)) {
      ng_stop("trait_checks needs check_geno: check lines are supplied in their own genotype ",
              "matrix, separate from the candidate parents, and are never crossed")
    }
    if (!identical(prediction_mode, "trait_by_trait")) {
      ng_stop("trait_checks requires prediction_mode = 'trait_by_trait' (checks are keyed by trait)")
    }
    if (!(inherits(trait_checks, "data.frame") && all(c("trait", "check") %in% names(trait_checks)))) {
      ng_stop("trait_checks must be a data.frame with trait + check columns")
    }
    # Progeny per family is the breeder's number, not ours. It scales P(beat check) directly,
    # so there is no default: a made-up family size would silently drive a reported probability.
    kp <- if (is.null(check_progeny_size)) NA_integer_ else
      suppressWarnings(as.integer(check_progeny_size[[1L]]))
    if (!length(kp) || is.na(kp) || kp < 1L) {
      ng_stop("check_progeny_size is required with trait_checks: give the number of progeny ",
              "you will raise per family. It sets P(beat check) -- the chance a cross throws a ",
              "line beating the check -- so it must be your program's figure, not a default.")
    }
    check_geno <- ng_align_check_geno(check_geno, colnames(geno), ploidy = marker_ploidy)
    clash <- intersect(rownames(check_geno), rownames(geno))
    if (length(clash)) {
      ng_stop("check id(s) also a candidate parent: ", paste(clash, collapse = ", "),
              ". A check is a benchmark, not breeding material; give it a distinct id or ",
              "remove it from the parent genotypes.")
    }
    tdir <- stats::setNames(direction_canonical$direction, direction_canonical$trait)
    tc_direction <- if (is.null(trait_checks$direction)) NA else trait_checks$direction
    tc_spec <- ng_trait_check_spec(trait_checks$trait, trait_checks$check,
                                   direction = tc_direction, trait_direction = tdir)
    # cross_table columns are named with the sanitised trait name; carry it as the lookup key
    # so a trait like "Days to flower" resolves to Days_to_flower_mean rather than erroring.
    tc_spec$column_key <- ng_run_cp_clean_trait_name(tc_spec$trait)
    missing_chk <- setdiff(tc_spec$check, rownames(check_geno))
    if (length(missing_chk)) {
      ng_stop("trait_checks names check line(s) absent from check_geno: ",
              paste(missing_chk, collapse = ", "))
    }
    check_values <- list(); check_source <- character(0)
    # Resolved LAZILY, only the first time some trait actually needs check_pheno -- see below.
    check_pheno_id_col <- NULL
    for (tr in unique(tc_spec$trait)) {
      if (!(tr %in% names(effects_list))) {
        ng_stop("trait_checks references a trait not present in trait_direction/effects: ", tr)
      }
      # The per-trait mean_source stamped on the scored table is the authority: the check must
      # land on whatever source produced the cross means for this trait, or the reference line
      # would sit on a different scale from the axis it is drawn on.
      src <- as.character(trait_mean_source[[tr]])
      check_source[[tr]] <- src
      # Explicit check_records wins where both are supplied: check_pheno only fills in the
      # per-trait records the caller did NOT already give via check_records. A GEBV source is
      # never filled from check_pheno -- ng_check_records_from_pheno() itself returns NULL for
      # a GEBV source, since the value must then come from the check's own markers.
      recs <- check_records[[tr]]
      # A GEBV-sourced trait never consults check_pheno at all -- ng_check_records_from_pheno()
      # returns NULL immediately for a GEBV source (the value comes from the check's own markers
      # instead) -- so it must never trigger id-column resolution either. Gating on
      # `is.null(recs)` alone was NOT sufficient: for a GEBV trait `recs` is also NULL (no
      # check_records supplied), so the branch below still ran, called ng_run_cp_id_col() on
      # check_pheno, and hard-errored on an unrecognised id column even though nothing in this
      # trait's evaluation needed check_pheno at all.
      if (is.null(recs) && !is.null(check_pheno) && !startsWith(src, "GEBV")) {
        # check_pheno's id column resolves exactly the way every other table's does, and ONLY
        # once we know some trait actually needs check_pheno (this gate) -- resolving eagerly,
        # above the loop, would fail a run over a check_pheno table no trait ever consults (e.g.
        # every trait supplied via check_records, or GEBV-sourced): attaching an input nothing
        # reads must never turn a working run into a failing one. An explicit id_col passes
        # through unchanged so a genuinely missing explicit column still hits
        # ng_check_records_from_pheno()'s own check and its existing message, unaltered; an
        # omitted id_col auto-detects via ng_run_cp_id_col()'s own candidate list (parent,
        # parent_id, id, name, line, entry, NAME) -- the SAME list/function phenotype and
        # genotype use, not a second hand-rolled copy of it. Resolved once (cached in the loop
        # variable) and reused across every trait that shares this same check_pheno table.
        if (is.null(check_pheno_id_col)) {
          check_pheno_id_col <- if (is.null(id_col)) {
            ng_run_cp_id_col(check_pheno, NULL, "check_pheno")
          } else {
            id_col
          }
        }
        col <- trait_spec$column[match(tr, trait_spec$trait)]
        conv <- ng_check_records_from_pheno(check_pheno, check_pheno_id_col,
                                            stats::setNames(list(col), tr), src)
        recs <- if (is.null(conv)) NULL else conv[[tr]]
      }
      check_values[[tr]] <- ng_check_reference_value(src, check_geno, effects_list[[tr]],
                                                     check_records = recs)
    }
    # Design section 5's "not-evaluable rule" promises "NA with a diagnostic" -- until now there
    # was no diagnostic anywhere in the run: a phenotypic source with no matching check_pheno/
    # check_records record silently produced NA with no warning at all. Warn once per trait whose
    # check value came back non-finite, naming the trait, the resolved source, and pointing the
    # breeder at check_pheno (or check_records) as the fix.
    for (tr in unique(tc_spec$trait)) {
      ck <- tc_spec$check[match(tr, tc_spec$trait)]
      cv <- check_values[[tr]]
      val <- if (!is.null(cv) && ck %in% names(cv)) cv[[ck]] else NA_real_
      if (!length(val) || !is.finite(val)) {
        warning(sprintf(
          "check reference for trait '%s' (check '%s') is not evaluable: the resolved mean ",
          tr, ck),
          sprintf("source '%s' has no record for this check. Supply a matching entry via ",
                  check_source[[tr]]),
          "check_pheno (or check_records) so this check can be compared against the cross ",
          "means; its reference columns and 'checks_all_ok' will report NA for this trait ",
          "until then.", call. = FALSE)
      }
    }
    cross_table <- ng_attach_check_reference(cross_table, tc_spec, trait_values = NULL,
                                             check_values = check_values,
                                             k_progeny = kp)
    # Multi-trait only: P(a progeny beats every check at once). Skipped for a single check,
    # where p_beat_all_checks would just duplicate the per-trait column.
    if (nrow(tc_spec) > 1L) {
      cross_table <- ng_attach_joint_check_probability(
        cross_table, tc_spec, check_values, k_progeny = kp,
        cross_trait_cov = ctc)
    }
    trait_check_reference <- list(
      active = tc_spec, values = check_values, source = check_source,
      progeny_size = kp,
      diagnostics = attr(cross_table, "check_reference_diagnostics"))
    attr(cross_table, "check_reference_diagnostics") <- NULL
  }

  scored_crosses <- ng_score_breeder_objective(
    cross_table,
    objective,
    threshold_penalty_weight = threshold_penalty_weight,
    threshold_penalty_autoscale = threshold_penalty_autoscale,
    phenotypic_covariance = phenotypic_covariance,
    genetic_covariance = genetic_covariance
  )
  # Marker steering (Module 4): attach the target-allele score (and, when lambda_marker != 0,
  # the blended merit) to the reported candidate table and drive the criterion the native /
  # AlphaMate-style allocators optimize. The OCS path performs the same blend internally.
  # The index metadata (method + resolved weights / solved coefficients) is what the multi-trait
  # portfolio axes are built from; capture it here because `attr` does not survive the row
  # subsetting the allocator performs on the way to the plan.
  ctx <- ng_ctx_put(ctx, multi_trait_meta = attr(scored_crosses, "multi_trait"))
  allocation_criterion_col <- "multi_trait_score"
  if (!is.null(marker_target_spec)) {
    scored_crosses <- ng_apply_marker_management(
      scores = scored_crosses, geno = geno,
      marker_target_spec = marker_target_spec, ploidy = marker_ploidy,
      gain_col = "multi_trait_score", lambda_marker = lambda_marker
    )
    if ("marker_adjusted_gain" %in% names(scored_crosses)) {
      allocation_criterion_col <- "marker_adjusted_gain"
    }
  }
  ctx <- ng_ctx_put(
    ctx,
    cross_table = cross_table,
    objective = objective,
    n_candidates_pre_lethal = n_candidates_pre_lethal,
    trait_check_reference = trait_check_reference,
    scored_crosses = scored_crosses,
    allocation_criterion_col = allocation_criterion_col
  )
  ctx
}

ng_cp__stage_allocate <- function(ctx) {
  list2env(ctx, environment())
  parent_kinship <- if (isTRUE(use_ocs) || !identical(allocation_method, "ocs")) ng_parent_kinship(geno, method = grm_method) else NULL
  if (identical(allocation_method, "ocs")) {
    # The mate-selection module knobs (strategy / diversity_emphasis, progeny inbreeding,
    # breeder constraints, cost/logistics) are forwarded through ng_optimize_breeder_
    # selection_plan's ... into ng_optimize_mating_plan, so the user-friendly runner
    # reaches the same capabilities as the low-level API.
    plan <- ng_optimize_breeder_selection_plan(
      scores = cross_table,
      objective = objective,
      n_crosses = n_crosses,
      parent_kinship = parent_kinship,
      optimizer_method = optimizer_method,
      threshold_penalty_weight = threshold_penalty_weight,
      threshold_penalty_autoscale = threshold_penalty_autoscale,
      phenotypic_covariance = phenotypic_covariance,
      genetic_covariance = genetic_covariance,
      marker_target_spec = marker_target_spec,
      marker_geno = geno,
      lambda_marker = lambda_marker,
      marker_ploidy = marker_ploidy,
      max_crosses_per_parent = max_crosses_per_parent,
      min_crosses_per_parent = min_crosses_per_parent,
      min_unique_parents = min_unique_parents,
      max_pair_kinship = max_pair_kinship,
      committed_crosses = committed_crosses,
      parent_group = parent_group,
      group_permission = group_permission,
      group_quota = group_quota,
      cost_col = cost_col,
      budget = budget,
      lambda_cost = lambda_cost,
      logistic_col = logistic_col,
      lambda_logistic = lambda_logistic,
      strategy = strategy,
      diversity_emphasis = diversity_emphasis,
      target_coancestry = target_coancestry,
      lambda_group = if (isTRUE(use_ocs)) lambda_group else 0,
      lambda_mating = if (isTRUE(use_ocs)) lambda_mating else 0,
      lambda_progeny_inbreeding = if (isTRUE(use_ocs)) lambda_progeny_inbreeding else 0,
      mate_relatedness = if (isTRUE(use_ocs)) match.arg(mate_relatedness, c("off", "avoid_inbreeding", "favor_complementarity")) else "off",
      mate_relatedness_weight = if (isTRUE(use_ocs)) mate_relatedness_weight else 0,
      lambda_parent_use = if (isTRUE(use_ocs)) lambda_parent_use else 0,
      lambda_parent_use_mode = lambda_parent_use_mode,
      local_iter = local_iter,
      ocs_iter = ocs_iter,
      evol_solutions = evol_solutions,
      evol_iterations = evol_iterations,
      evol_stop = evol_stop,
      evol_seed = evol_seed
    )
  } else if (identical(allocation_method, "alphamate_style")) {
    alpha_max_contrib <- if (is.null(alphamate_max_contributions)) max_crosses_per_parent else alphamate_max_contributions
    plan <- ng_alphamate_style_select(
      scores = scored_crosses,
      criterion_col = allocation_criterion_col,
      n_crosses = n_crosses,
      parent_kinship = parent_kinship,
      mode = alphamate_mode,
      target_degree = alphamate_target_degree,
      max_contributions = alpha_max_contrib,
      lambda_group = alphamate_lambda_group,
      lambda_grid = alphamate_lambda_grid,
      method = optimizer_method,
      local_iter = local_iter
    )
  } else if (identical(allocation_method, "alphamate_executable")) {
    alpha_max_contrib <- if (is.null(alphamate_max_contributions)) max_crosses_per_parent else alphamate_max_contributions
    alpha_executable <- if (is.null(alphamate_executable)) ng_alphamate_default_executable() else alphamate_executable
    alpha_workdir <- if (is.null(alphamate_workdir)) tempfile("ng_alphamate_run_") else alphamate_workdir
    plan <- ng_select_alphamate(
      scores = scored_crosses,
      criterion_col = allocation_criterion_col,
      n_crosses = n_crosses,
      parent_kinship = parent_kinship,
      executable = alpha_executable,
      runtime_path = alphamate_runtime_path,
      target_degree = alphamate_target_degree,
      max_contributions = alpha_max_contrib,
      number_of_parents = alphamate_number_of_parents,
      evol_solutions = alphamate_evol_solutions,
      evol_iterations = alphamate_evol_iterations,
      evol_stop = alphamate_evol_stop,
      n_threads = alphamate_n_threads,
      workdir = alpha_workdir,
      keep_files = alphamate_keep_files,
      mode = alphamate_mode
    )
  } else {
    ng_stop("Unsupported allocation_method: ", allocation_method)
  }
  plan_summary <- attr(plan, "summary")
  plan_summary$allocation_method <- allocation_method
  attr(plan, "summary") <- plan_summary
  ctx <- ng_ctx_put(ctx, plan = plan)
  ctx
}

ng_cp__stage_rank <- function(ctx) {
  list2env(ctx, environment())
  multi_trait_meta <- ctx$multi_trait_meta
  selected <- ng_rank_cross_priority(
    crosses = plan,
    score_col = "multi_trait_score",
    pair_kinship_col = "pair_kinship",
    threshold_violation_col = "multi_trait_threshold_violation",
    score_weight = priority_score_weight,
    kinship_weight = priority_kinship_weight,
    threshold_weight = priority_threshold_weight,
    breaks = priority_breaks,
    labels = priority_labels,
    sort = TRUE
  )
  attr(selected, "summary") <- attr(plan, "summary")

  # Cross-priority portfolio + risk annotation (Tasks 1-4): attaches cross_level/cross_upside/
  # cross_confidence/risk_bin/confidence_method/portfolio_profile to the selected and candidate
  # tables, plus a summary diagnostics block. Single-trait resolves level/upside on the trait
  # itself; multi-trait resolves them on the SELECTION INDEX (see R/37) -- same six columns and
  # the same quadrant semantics either way, so the reporting surface does not fork.
  priority_risk_diagnostics <- NULL
  traits_clean <- vapply(trait_spec$trait, ng_run_cp_clean_trait_name, character(1L),
                         USE.NAMES = FALSE)
  # The merit's sqrt(X) term is effect-based (i.e. draws on marker-effect estimation error,
  # not just the mid-parent mean) for every metric except mean/parent_distance, and for
  # usefulness only when its variance source isn't the effect-free parent_distance proxy.
  effect_based_x <- !(trait_value_metric %in% c("mean", "parent_distance", "le")) &&
    !(identical(trait_value_metric, "usefulness") &&
      uc_variance_source %in% c("parent_distance", "le"))
  lvl_cols <- paste0(traits_clean, "_mean_gebv")
  vpm_cols <- paste0(traits_clean, "_vpm")
  pev_cols <- paste0(traits_clean, "_midparent_pev")

  if (length(trait_spec$column) == 1L) {
    sd_col <- paste0(traits_clean, "_post_sd")
    tn_col <- paste0(traits_clean, "_post_topn")
    ann_one <- function(tbl) {
      if (!nrow(tbl)) return(tbl)
      if (!all(c(lvl_cols, vpm_cols, "parent1", "parent2") %in% names(tbl))) return(tbl)
      pev <- if (pev_cols %in% names(tbl)) suppressWarnings(as.numeric(tbl[[pev_cols]])) else NULL
      psd <- if (sd_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[sd_col]])) else NULL
      ptn <- if (tn_col %in% names(tbl)) suppressWarnings(as.numeric(tbl[[tn_col]])) else NULL
      ng_annotate_cross_priority(tbl, level = tbl[[lvl_cols]], vpm = tbl[[vpm_cols]],
                                 pev = pev, effect_based_x = effect_based_x,
                                 post_sd = psd, prob_top_tier = ptn)
    }
  } else {
    # Index coefficients: the linear solve for the economic_index / desired_gain families, the
    # resolved (normalized, positive) trait weights for the weighted / rank-sum family.
    mt_method <- if (is.null(multi_trait_meta$method)) "weighted" else multi_trait_meta$method
    mt_coef <- multi_trait_meta$economic_index_coefficients
    if (is.null(mt_coef)) mt_coef <- multi_trait_meta$desired_gain_coefficients
    if (is.null(mt_coef)) mt_coef <- multi_trait_meta$weights
    wf_cols <- paste0("wf_var_", trait_spec$trait)
    mt_ready <- !is.null(mt_coef) &&
      all(c(lvl_cols, vpm_cols, wf_cols, "parent1", "parent2") %in% names(scored_crosses))
    # One basis for the whole run, resolved on the candidate pool: the per-trait scale is an
    # IQR of the rows it sees, so deriving it per table would put the plan and the candidate
    # pool on different index axes.
    mt_basis <- if (mt_ready && nrow(scored_crosses) > 0L) {
      ng_multitrait_index_reference(scored_crosses, trait_order = trait_spec$trait,
                                    mean_gebv_cols = lvl_cols,
                                    directions = trait_spec$direction, coefficients = mt_coef)
    } else NULL
    ann_one <- function(tbl) {
      if (!nrow(tbl) || is.null(mt_basis)) return(tbl)
      if (!all(c(lvl_cols, vpm_cols, wf_cols, "parent1", "parent2") %in% names(tbl))) return(tbl)
      ng_annotate_cross_priority_multitrait(
        tbl, trait_order = trait_spec$trait, mean_gebv_cols = lvl_cols, vpm_cols = vpm_cols,
        directions = trait_spec$direction, coefficients = mt_coef,
        pev_cols = if (all(pev_cols %in% names(tbl))) pev_cols else NULL,
        effect_based_x = effect_based_x, index_method = mt_method, basis = mt_basis)
    }
  }
  # Resolve every relative cut (confidence min/max, risk tertiles, portfolio
  # medians) once on the full post-filter candidate pool, then copy those exact
  # annotations to the selected subset. Otherwise one cross can receive two
  # different labels in the two tables returned by the same run.
  scored_crosses <- ann_one(scored_crosses)
  selected <- ng_copy_cross_priority_annotations(selected, scored_crosses)
  if (nrow(selected) > 0L && "confidence_method" %in% names(selected)) {
    cm1 <- selected$confidence_method[[1L]]
    priority_risk_diagnostics <- list(
      confidence_method = if (is.null(cm1)) NA_character_ else cm1,
      basis = if (length(trait_spec$column) == 1L) "single_trait" else "multi_trait_index",
      n_traits = nrow(trait_spec),
      reference_population = "candidate_crosses_after_filters_before_allocation",
      n_reference_crosses = nrow(scored_crosses),
      n_by_risk    = as.list(table(selected$risk_bin)),
      n_by_profile = as.list(table(selected$portfolio_profile)),
      top_tier_high_risk = sum(as.character(selected$priority_tier) ==
                                 levels(selected$priority_tier)[[1L]] &
                               as.character(selected$risk_bin) == "high", na.rm = TRUE),
      posterior_used = identical(cm1, "posterior_ci")
    )
    if (length(trait_spec$column) > 1L) {
      ib <- attr(selected, "index_basis")
      priority_risk_diagnostics$index_method <- multi_trait_meta$method
      priority_risk_diagnostics$index_weights <- if (is.null(ib)) NULL else as.list(ib$w)
      priority_risk_diagnostics$upside_method <- "exact_within_family_cov"
      # "linearized_rank_index" means multi_trait_score is a RANK index, so the quadrant is
      # indicative rather than a decomposition of the ranked merit -- the frontend must badge
      # it. See ng_multitrait_portfolio_basis.
      pb <- selected$portfolio_basis[[1L]]
      priority_risk_diagnostics$portfolio_basis <- if (is.null(pb)) NA_character_ else pb
      # Per-trait accounting for the index, over the CANDIDATE pool (a plan of n_crosses is too
      # small a sample to characterise a trait). Answers the two questions a breeder asks of a
      # multi-trait run without opening the cross table: which trait is carrying the risk, and
      # is any trait effectively absent from the index?
      if (!is.null(ib) && nrow(scored_crosses) > 0L) {
        grab_rt <- function(cols) {
          if (!all(cols %in% names(scored_crosses))) return(NULL)
          matrix(vapply(cols, function(cc)
            suppressWarnings(as.numeric(scored_crosses[[cc]])), numeric(nrow(scored_crosses))),
            nrow = nrow(scored_crosses), dimnames = list(NULL, trait_spec$trait))
        }
        col_mean <- function(M) {
          if (is.null(M)) return(rep(NA_real_, nrow(trait_spec)))
          v <- colMeans(M, na.rm = TRUE)
          v[!is.finite(v)] <- NA_real_
          v
        }
        var_share <- col_mean(ng_multitrait_variance_shares(
          scored_crosses, trait_spec$trait, ib$w, vpm = grab_rt(vpm_cols)))
        pm <- grab_rt(pev_cols)
        pev_share <- col_mean(if (is.null(pm)) NULL else
          ng_multitrait_pev_shares(pm, ib$w, trait_spec$trait))
        pred_r2 <- col_mean(grab_rt(paste0(traits_clean, "_cv_predictive_r2")))
        priority_risk_diagnostics$index_traits <- lapply(seq_len(nrow(trait_spec)), function(k)
          list(trait = trait_spec$trait[[k]],
               direction = trait_spec$direction[[k]],
               weight = unname(ib$w[[k]]),
               marker_effect_reliability = NA_real_,
               cv_predictive_r2 = unname(pred_r2[[k]]),
               mean_variance_share = unname(var_share[[k]]),
               mean_pev_share = unname(pev_share[[k]])))
        # Traits that buy more uncertainty than opportunity: carrying above an equal share of
        # the index PEV, and more than twice as much of the index's risk as of its spread. This
        # is the actionable signal (drop the trait, or phenotype it better) and it is what
        # actually catches a poorly-estimated trait -- a bare "contributes little" threshold
        # misses one that contributes little gain but a great deal of error.
        equal_share <- 1 / nrow(trait_spec)
        flagged <- is.finite(pev_share) & pev_share > equal_share &
          pev_share > 2 * pmax(var_share, 0, na.rm = FALSE)
        flagged[is.na(flagged)] <- FALSE
        priority_risk_diagnostics$risk_disproportionate_traits <-
          as.list(trait_spec$trait[flagged])

        # Is the index confidence actually an INDEX quantity? sum_k w_k^2 PEV_k is exact given
        # the fitted models, but the per-trait PEVs are only comparable across traits when
        # their residual variances are -- and sigma_e^2 is estimated in sample. A trait whose
        # ridge lambda lands on the grid floor can interpolate its training rows, collapsing
        # sigma_e^2 and reporting a near-zero PEV regardless of how well the trait is really
        # predicted. When one trait then carries essentially all of the index PEV, risk_bin is
        # a single-trait statement wearing an index label, and must not be read otherwise.
        conc <- if (any(is.finite(pev_share))) max(pev_share, na.rm = TRUE) else NA_real_
        priority_risk_diagnostics$pev_concentration <- conc
        priority_risk_diagnostics$pev_concentration_note <-
          if (is.finite(conc) && conc > 0.9) sprintf(paste0(
            "index confidence is dominated by '%s' (%.0f%% of the index prediction-error ",
            "variance), so risk_bin ranks crosses on that trait rather than on the index. ",
            "Per-trait PEVs are only comparable when their residual variances are; check ",
            "ridge_lambda and cv_predictive_r2 per trait (a lambda at the grid floor ",
            "interpolates its training rows and reports a near-zero PEV)."),
            trait_spec$trait[[which.max(replace(pev_share, !is.finite(pev_share), -Inf))]],
            100 * conc) else NA_character_
      }
      priority_risk_diagnostics$portfolio_basis_note <-
        if (identical(pb, "linearized_rank_index")) paste0(
          "index method '", multi_trait_meta$method, "' combines rank-normalized traits, so no ",
          "linear index exists in genetic units; cross_level/cross_upside use the trait weights ",
          "as standardized-unit coefficients and are indicative, not a decomposition of ",
          "multi_trait_score. Supply economic_weight (or desired_change) per trait for a solved ",
          "linear index.") else NA_character_
    }
  }

  # Constraint / marker-management / cost diagnostics: a single structured record of what
  # the breeder knobs actually did to the plan, so the frontend can surface visible run
  # notes (plan shrink, min-unique relaxation, lethal-carrier drops, marker steering, budget
  # binding) instead of leaving these effects buried in stderr messages.
  ps <- attr(plan, "summary")
  `%||%` <- function(a, b) if (is.null(a)) b else a
  is_ocs <- identical(allocation_method, "ocs")
  n_delivered <- nrow(selected)
  lethal_dropped <- if (!is.null(lethal_spec)) max(0L, n_candidates_pre_lethal - nrow(cross_table)) else 0L
  cost_vals <- if (!is.null(cost_col) && cost_col %in% names(selected)) suppressWarnings(as.numeric(selected[[cost_col]])) else NULL
  constraint_diagnostics <- list(
    allocation_method = allocation_method,
    # (2) breeder mating constraints
    n_crosses_requested = as.integer(n_crosses),
    n_crosses_delivered = as.integer(n_delivered),
    constraints_reduced_plan = isTRUE(ps$constraints_reduced_plan),
    n_committed = as.integer(ps$n_committed %||% 0L),
    min_use_if_used = ps$min_use_if_used %||% NA_integer_,
    min_unique_requested = ps$balanced_min_unique_requested %||% NA_integer_,
    min_unique_used = ps$balanced_min_unique_used %||% NA_integer_,
    min_unique_relaxation_attempts = as.integer(ps$balanced_min_unique_relaxation_attempts %||% 0L),
    group_quota_active = !is.null(group_quota) && length(group_quota) > 0L,
    group_permission_active = !is.null(group_permission),
    max_pair_kinship = if (is.null(max_pair_kinship)) NA_real_ else as.numeric(max_pair_kinship),
    # (3) marker steering & lethal guarding
    lethal_active = !is.null(lethal_spec),
    lethal_n_loci = if (!is.null(lethal_spec)) nrow(as.data.frame(lethal_spec)) else 0L,
    lethal_candidates_pre = as.integer(n_candidates_pre_lethal),
    lethal_dropped = as.integer(lethal_dropped),
    lethal_dropped_from_plan = isTRUE(drop_lethal_carrier_crosses),
    marker_steering_active = !is.null(marker_target_spec),
    marker_n_loci = if (!is.null(marker_target_spec)) nrow(as.data.frame(marker_target_spec)) else 0L,
    marker_lambda = as.numeric(lambda_marker),
    marker_blended = identical(allocation_criterion_col, "marker_adjusted_gain"),
    # (4) cost & logistics
    budget_active = is.finite(budget) && is_ocs,
    budget = as.numeric(budget),
    plan_total_cost = if (!is.null(cost_vals)) sum(cost_vals[is.finite(cost_vals)]) else NA_real_,
    cost_emphasis = as.numeric(lambda_cost),
    logistic_emphasis = as.numeric(lambda_logistic),
    logistic_active = !is.null(logistic_col)
  )

  output_files <- ng_run_cp_output_files(
    output_dir = output_dir,
    output_file = output_file,
    selected_crosses = selected,
    candidate_crosses = scored_crosses,
    trait_spec = trait_spec,
    qc = qc,
    write_outputs = write_outputs,
    write_figures = write_figures,
    n_crosses = n_crosses,
    include_trait_gebv = include_trait_gebv,
    trait_check_reference = ctx$trait_check_reference,
    multi_trait_meta = ctx$multi_trait_meta,
    trait_value_metric = ctx$trait_value_metric
  )
  ctx <- ng_ctx_put(
    ctx,
    selected = selected,
    scored_crosses = scored_crosses,
    priority_risk_diagnostics = priority_risk_diagnostics,
    constraint_diagnostics = constraint_diagnostics,
    output_files = output_files
  )
  ctx
}

# The ng_cp__* stage functions unpack the shared ctx list into their local
# environment via list2env(); codetools cannot see those injected bindings, so
# declare the union of ctx fields they read as package globals to keep R CMD
# check clean. Regenerate with codetools::findGlobals over the ng_cp__* fns if
# ctx gains fields.
utils::globalVariables(c(
  "allocation_criterion_col", "allocation_method", "alphamate_evol_iterations", "alphamate_evol_solutions",
  "alphamate_evol_stop", "alphamate_executable", "alphamate_keep_files", "alphamate_lambda_grid",
  "alphamate_lambda_group", "alphamate_max_contributions", "alphamate_mode", "alphamate_n_threads",
  "alphamate_number_of_parents", "alphamate_runtime_path", "alphamate_target_degree", "alphamate_workdir",
  "assume_inbred", "parent_type", "phased_haplotypes", "bp_per_cm", "budget", "burn_in",
  "check_geno", "check_pheno", "check_progeny_size", "check_records",
  "committed_crosses", "constraint_diagnostics", "cost_col",
  "cross_cost", "cross_table", "ctc", "direction_column_col", "direction_columns",
  "direction_direction_col", "direction_file", "direction_trait_col", "diversity_emphasis",
  "drop_lethal_carrier_crosses", "duplicate_action", "duplicate_maf_min", "duplicate_max_missing_prop",
  "duplicate_min_compared_markers", "duplicate_threshold", "effect_summary", "effects_list",
  "evol_iterations", "evol_seed", "evol_solutions", "evol_stop",
  "geno", "genotype", "genotype_file",
  "genotype_id_col_used", "grm_method", "group_permission", "group_quota",
  "id_col", "ids", "include_trait_gebv", "index_col", "phenotypic_covariance", "genetic_covariance",
  "index_direction", "lambda_cost", "lambda_group", "lambda_logistic",
  "lambda_marker", "lambda_mating", "lambda_parent_use", "lambda_parent_use_mode",
  "lambda_progeny_inbreeding", "ld_backend", "ld_maf_threshold", "ld_ploidy",
  "ld_pruning", "ld_r2_threshold", "ld_window", "lethal_spec",
  "local_iter", "logistic_col", "map_chr_col", "map_file",
  "map_marker_col", "map_pos_bp_col", "map_pos_cm_col", "map_pos_cm_divisor",
  "map_pos_col", "map_position_unit", "marker_map", "marker_map_std",
  "marker_ploidy", "marker_target_spec", "mate_relatedness", "mate_relatedness_weight",
  "max_crosses_per_parent", "max_pair_kinship", "method_varPMV", "min_crosses_per_parent",
  "min_effect_reliability", "min_unique_parents", "multi_trait_method", "n_candidates_pre_lethal",
  "n_crosses", "n_iter", "n_threads", "objective",
  "ocs_iter", "optimizer", "optimizer_method", "output_dir",
  "output_file", "output_files", "parallel_backend", "parallel_cores_used",
  "parent_group", "pheno", "phenotype", "phenotype_file",
  "phenotype_id_col_used", "plan", "posterior_effects_list", "posterior_method",
  "posterior_n_draws", "posterior_predictions_list", "prediction_mode", "priority_breaks",
  "priority_kinship_weight", "priority_labels", "priority_score_weight", "priority_threshold_weight",
  "qc", "recomb_model", "ril_mode", "run_posterior_prediction",
  "scored_crosses", "seed", "selected", "selection_prop",
  "strategy", "target", "target_coancestry", "threshold_penalty_autoscale",
  "threshold_penalty_weight", "threshold_policy", "training_genotype", "training_genotype_file",
  "training_genotype_id_col", "training_ids", "training_only_count", "training_phenotype",
  "training_phenotype_file", "training_phenotype_id_col",
  "trait_check_reference", "trait_checks", "trait_direction", "trait_mean_source",
  "trait_scores", "trait_spec", "trait_value_metric", "trait_value_metric_input",
  "uc_variance_source_input", "trait_weights",
  "traits_to_use", "uc_variance_source", "use_cpp", "use_ocs",
  "use_parallel", "write_figures", "write_outputs"
))

ng_cp__assemble_result <- function(ctx) {
  list2env(ctx, environment())
  ld_pruning_report <- ctx$ld_pruning_report
  trait_check_reference <- ctx$trait_check_reference
  priority_risk_diagnostics <- ctx$priority_risk_diagnostics
  result <- list(
    prediction_mode = prediction_mode,
    qc = qc,
    cleaned_data = list(
      phenotype = pheno,
      genotype = geno,
      marker_map = marker_map_std
    ),
    ld_pruning_report = ld_pruning_report,
    trait_direction = trait_spec,
    input_match_audit = list(
      phenotype_id_col = phenotype_id_col_used,
      genotype_id_col = genotype_id_col_used,
      internal_id_col = "parent",
      matched_parent_count = as.integer(length(ids)),
      training_only_count = as.integer(training_only_count),
      training_ids = training_ids,
      effect_training_n = as.integer(length(ids) + training_only_count),
      parent_order = ids,
      marker_count = as.integer(ncol(geno)),
      marker_order = colnames(geno),
      map_marker_col = attr(marker_map_std, "map_marker_col"),
      map_chr_col = attr(marker_map_std, "map_chr_col"),
      map_position_col = attr(marker_map_std, "map_position_col"),
      map_position_unit = attr(marker_map_std, "map_position_unit"),
      bp_per_cm = attr(marker_map_std, "bp_per_cm"),
      direction_columns = direction_columns,
      trait_columns = trait_spec$column
    ),
    effect_summary = do.call(rbind, effect_summary),
    marker_effects = effects_list,
    trait_scores = trait_scores,
    posterior_effects = posterior_effects_list,
    posterior_predictions = posterior_predictions_list,
    candidate_crosses = scored_crosses,
    selected_crosses = selected,
    objective = objective,
    plan_summary = attr(plan, "summary"),
    constraint_diagnostics = constraint_diagnostics,
    trait_check_reference = trait_check_reference,
    priority_risk_diagnostics = priority_risk_diagnostics,
    output_files = output_files,
    settings = list(
      trait_value_metric = trait_value_metric_input,
      uc_variance_source = uc_variance_source_input,
      method_varPMV = method_varPMV,
      multi_trait_method = multi_trait_method,
      progeny = target,
      ril_mode = ril_mode,
      recomb_model = recomb_model,
      selection_prop = selection_prop,
      run_posterior_prediction = run_posterior_prediction,
      posterior_method = posterior_method,
      n_iter = n_iter,
      burn_in = burn_in,
      posterior_n_draws = if (isTRUE(run_posterior_prediction)) posterior_n_draws else 0L,
      use_parallel = use_parallel,
      parallel_backend = parallel_backend,
      n_threads = parallel_cores_used,
      optimizer = optimizer,
      optimizer_method = optimizer_method,
      allocation_method = allocation_method,
      alphamate_mode = alphamate_mode,
      alphamate_target_degree = alphamate_target_degree,
      alphamate_max_contributions = if (is.null(alphamate_max_contributions)) max_crosses_per_parent else alphamate_max_contributions,
      use_ocs = isTRUE(use_ocs),
      n_crosses = n_crosses,
      max_crosses_per_parent = max_crosses_per_parent
    )
  )
  class(result) <- c("ng_cross_prediction_result", "list")
  result
}

ng_cp_pipeline <- list(
  qc = ng_cp__stage_qc,
  predict = ng_cp__stage_predict,
  index = ng_cp__stage_index,
  allocate = ng_cp__stage_allocate,
  rank = ng_cp__stage_rank
)

ng_run_cross_prediction <- function(phenotype_file = NULL,
                                    genotype_file = NULL,
                                    map_file = NULL,
                                    direction_file = NULL,
                                    phenotype = NULL,
                                    genotype = NULL,
                                    marker_map = NULL,
                                    trait_direction = NULL,
                                    id_col = NULL,
                                    training_genotype = NULL,
                                    training_phenotype = NULL,
                                    training_genotype_file = NULL,
                                    training_phenotype_file = NULL,
                                    training_genotype_id_col = NULL,
                                    training_phenotype_id_col = NULL,
                                    phenotype_id_col = NULL,
                                    genotype_id_col = NULL,
                                    direction_trait_col = NULL,
                                    direction_column_col = NULL,
                                    direction_direction_col = NULL,
                                    map_marker_col = NULL,
                                    map_chr_col = NULL,
                                    map_pos_col = NULL,
                                    map_pos_cm_col = NULL,
                                    map_pos_bp_col = NULL,
                                    map_position_unit = c("bp", "cM"),
                                    bp_per_cm = NULL,
                                    map_pos_cm_divisor = 1,
                                    prediction_mode = c("trait_by_trait", "index_as_trait"),
                                    traits_to_use = NULL,
                                    index_col = NULL,
                                    index_direction = "increase",
                                    trait_value_metric = c("usefulness", "pmv", "vpm", "parent_distance", "le", "var_complex", "mean"),
                                    uc_variance_source = c("pmv", "vpm", "parent_distance", "le"),
                                    multi_trait_method = "auto",
                                    trait_weights = NULL,
                                    phenotypic_covariance = NULL,
                                    genetic_covariance = NULL,
                                    threshold_policy = c("soft", "strict"),
                                    threshold_penalty_weight = 1.0,
                                    threshold_penalty_autoscale = TRUE,
                                    progeny = "DH",
                                    recomb_model = c("haldane", "kosambi"),
                                    selection_prop = 0.10,
                                    min_effect_reliability = 0.35,
                                    grm_method = c("vanraden", "yang"),
                                    method_varPMV = c("fast", "full_posterior"),
                                    ril_mode = "infinite",
                                    run_posterior_prediction = FALSE,
                                    posterior_method = c("mcmc", "closed_form"),
                                    n_iter = 5000L,
                                    burn_in = 500L,
                                    use_parallel = FALSE,
                                    n_threads = NULL,
                                    duplicate_action = c("remove", "report", "none"),
                                    duplicate_threshold = 0.995,
                                    duplicate_maf_min = 0.01,
                                    duplicate_max_missing_prop = 0.40,
                                    duplicate_min_compared_markers = 100L,
                                    ld_pruning = FALSE,
                                    ld_window = 100L,
                                    ld_r2_threshold = 0.9,
                                    ld_maf_threshold = 0.01,
                                    ld_ploidy = 2,
                                    ld_backend = c("auto", "cpp", "r"),
                                    n_crosses = 100,
                                    max_crosses_per_parent = 6,
                                    min_unique_parents = NULL,
                                    max_pair_kinship = Inf,
                                    optimizer = "auto",
                                    allocation_method = "ocs",
                                    use_ocs = TRUE,
                                    lambda_group = 0.05,
                                    lambda_mating = 0,
                                    lambda_parent_use = 0,
                                    lambda_parent_use_mode = c("absolute", "adaptive"),
                                    lambda_progeny_inbreeding = 0,
                                    mate_relatedness = c("off", "avoid_inbreeding", "favor_complementarity"),
                                    mate_relatedness_weight = 0,
                                    strategy = NULL,
                                    diversity_emphasis = NULL,
                                    target_coancestry = NULL,
                                    min_crosses_per_parent = 0,
                                    committed_crosses = NULL,
                                    parent_group = NULL,
                                    group_permission = NULL,
                                    group_quota = NULL,
                                    cross_cost = NULL,
                                    cost_col = NULL,
                                    budget = Inf,
                                    lambda_cost = 0,
                                    logistic_col = NULL,
                                    lambda_logistic = 0,
                                    lethal_spec = NULL,
                                    trait_checks = NULL,
                                    check_geno = NULL,
                                    check_records = NULL,
                                    check_pheno = NULL,
                                    check_progeny_size = NULL,
                                    include_trait_gebv = FALSE,
                                    marker_target_spec = NULL,
                                    lambda_marker = 0,
                                    marker_ploidy = 2,
                                    drop_lethal_carrier_crosses = TRUE,
                                    local_iter = 2000,
                                    ocs_iter = 5L,
                                    evol_solutions = 100L,
                                    evol_iterations = 200L,
                                    evol_stop = 40L,
                                    evol_seed = NULL,
                                    alphamate_mode = c("ModeOptTarget1", "ModeMaxCriterion", "ModeMinCoancestry"),
                                    alphamate_target_degree = 45,
                                    alphamate_max_contributions = NULL,
                                    alphamate_number_of_parents = NULL,
                                    alphamate_lambda_group = NULL,
                                    alphamate_lambda_grid = NULL,
                                    alphamate_executable = NULL,
                                    alphamate_runtime_path = Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = ""),
                                    alphamate_workdir = NULL,
                                    alphamate_keep_files = FALSE,
                                    alphamate_evol_solutions = 100L,
                                    alphamate_evol_iterations = 1000L,
                                    alphamate_evol_stop = 200L,
                                    alphamate_n_threads = 1L,
                                    priority_breaks = c(0.10, 0.35, 0.70, 1.00),
                                    priority_labels = c("highly_priority", "priority", "medium_priority", "low_priority"),
                                    priority_score_weight = 1.0,
                                    priority_kinship_weight = 0.15,
                                    priority_threshold_weight = 1.0,
                                    output_dir = NULL,
                                    output_file = "crossing_plan.xlsx",
                                    write_outputs = FALSE,
                                    write_figures = FALSE,
                                    parent_type = c("inbred", "dh", "ril"),
                                    assume_inbred = NULL,
                                    phased_haplotypes = NULL,
                                    use_cpp = TRUE,
                                    seed = 1L) {
  config <- mget(names(formals()))
  # parent_type governs the heterozygosity audit: 'inbred'/'dh' (fully fixed
  # lines) block het as a data error; 'ril' accepts residual het. Legacy
  # assume_inbred is deprecated (reconciled with a one-time warning). Canonicalise
  # once here so the staged config and inner calls carry a single clean value.
  parent_type <- ng_reconcile_parent_type(parent_type, assume_inbred)
  assume_inbred <- NULL
  config$parent_type <- parent_type
  config$assume_inbred <- NULL
  ctx <- ng_cp__build_ctx(config)
  for (s in ng_cp_stage_order()) {
    ctx <- ng_cp_pipeline[[s]](ctx)
    # The qc stage (R/39 ng_cp__stage_qc) no longer throws on a blocker --
    # it returns ctx with ctx$qc$status == "blocker" so the staged runner
    # (R/45 ng_run_stage) can persist qc.json for the frontend. The one-shot
    # driver reinstates the original throw-on-blocker behavior here, at the
    # exact original point (right after qc, before predict), with the exact
    # original message.
    if (identical(s, "qc") && identical(ctx$qc$status, "blocker")) {
      blockers <- ctx$qc$issues$message[ctx$qc$issues$severity == "blocker"]
      ng_stop(
        "Input QC found blocker issues before cross prediction: ",
        paste(utils::head(blockers, 5L), collapse = " | ")
      )
    }
  }
  ng_cp__assemble_result(ctx)
}
