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

  bp_per_cm <- if (is.null(bp_per_cm)) map_pos_cm_divisor else bp_per_cm
  bp_per_cm <- suppressWarnings(as.numeric(bp_per_cm[[1L]]))
  if (!is.finite(bp_per_cm) || bp_per_cm <= 0) ng_stop("bp_per_cm must be a positive number")

  if (identical(map_position_unit, "bp")) {
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
      pos_cm * bp_per_cm
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
  out <- suppressWarnings(as.integer(x[[1L]]))
  if (!is.finite(out) || out < min_value) {
    ng_stop(name, " must be an integer >= ", min_value)
  }
  out
}

ng_run_cp_logical <- function(x, name) {
  if (length(x) != 1L || is.na(x)) ng_stop(name, " must be TRUE or FALSE")
  isTRUE(x)
}

ng_run_cp_pmv_col <- function(scored_trait, method_varPMV = "fast") {
  method_varPMV <- ng_run_cp_method_varPMV(method_varPMV)
  if (identical(method_varPMV, "fast")) return("dh_pmv_var")
  col <- "dh_pmv_var_full_posterior"
  if (!(col %in% names(scored_trait))) {
    ng_stop("method_varPMV = 'full_posterior' requires score column: ", col)
  }
  vals <- suppressWarnings(as.numeric(scored_trait[[col]]))
  if (!any(is.finite(vals))) {
    ng_stop("method_varPMV = 'full_posterior' produced no finite full-posterior PMV values")
  }
  col
}

ng_run_cp_parallel_cores <- function(parallel_cores, n_jobs) {
  if (!is.null(parallel_cores)) {
    cores <- ng_run_cp_integer(parallel_cores, "parallel_cores", min_value = 1L)
  } else {
    detected <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    cores <- if (is.finite(detected) && detected > 1L) detected - 1L else 1L
  }
  as.integer(max(1L, min(cores, n_jobs)))
}

ng_run_cp_apply <- function(jobs, fun, use_parallel = FALSE, parallel_cores = NULL) {
  if (!isTRUE(use_parallel) || length(jobs) < 2L) {
    out <- lapply(jobs, fun)
    attr(out, "parallel_backend") <- "serial"
    attr(out, "parallel_cores") <- 1L
    return(out)
  }
  cores <- ng_run_cp_parallel_cores(parallel_cores, length(jobs))
  if (cores < 2L) {
    out <- lapply(jobs, fun)
    attr(out, "parallel_backend") <- "serial"
    attr(out, "parallel_cores") <- 1L
    return(out)
  }
  if (!identical(.Platform$OS.type, "windows")) {
    out <- parallel::mclapply(jobs, fun, mc.cores = cores)
    attr(out, "parallel_backend") <- "mclapply"
    attr(out, "parallel_cores") <- cores
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
    attr(out, "parallel_cores") <- 1L
  } else {
    attr(out, "parallel_backend") <- "psock"
    attr(out, "parallel_cores") <- cores
  }
  out
}

ng_run_cp_variance_col <- function(trait_value_metric,
                                   uc_variance_source,
                                   scored_trait = NULL,
                                   method_varPMV = "fast") {
  metric <- trimws(tolower(as.character(trait_value_metric[[1L]])))
  source <- trimws(tolower(as.character(uc_variance_source[[1L]])))
  if (identical(metric, "uc")) metric <- source
  if (identical(metric, "pmv")) {
    if (is.null(scored_trait)) return("dh_pmv_var")
    return(ng_run_cp_pmv_col(scored_trait, method_varPMV))
  }
  if (identical(metric, "vpm")) return("dh_recomb_var")
  if (identical(metric, "var_simple")) return("var_simple")
  if (identical(metric, "var_complex")) {
    if (is.null(scored_trait)) return("dh_pmv_var")
    return(ng_run_cp_var_complex_col(scored_trait, method_varPMV))
  }
  if (identical(metric, "mean")) return(NA_character_)
  ng_stop("trait_value_metric must be one of: uc, pmv, vpm, var_simple, var_complex, mean")
}

ng_run_cp_var_complex_col <- function(scored_trait, method_varPMV = "fast") {
  pmv_col <- ng_run_cp_pmv_col(scored_trait, method_varPMV)
  candidates <- c(pmv_col, "dh_recomb_var", "var_simple")
  hit <- candidates[candidates %in% names(scored_trait)]
  if (!length(hit)) ng_stop("scored trait table is missing a usable var_complex variance column")
  hit[[1L]]
}

ng_run_cp_trait_value <- function(scored_trait,
                                  direction,
                                  trait_value_metric = "uc",
                                  uc_variance_source = "pmv",
                                  selection_prop = 0.10,
                                  method_varPMV = "fast") {
  direction <- ng_run_cp_direction(direction)[[1L]]
  metric <- trimws(tolower(as.character(trait_value_metric[[1L]])))
  mean_value <- suppressWarnings(as.numeric(scored_trait$cross_mean_blend))
  if (identical(metric, "mean")) return(mean_value)
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
                                   n_crosses) {
  files <- list()
  if (!isTRUE(write_outputs) && !isTRUE(write_figures)) return(files)
  if (is.null(output_dir) || !nzchar(as.character(output_dir[[1L]]))) {
    ng_stop("output_dir is required when write_outputs or write_figures is TRUE")
  }
  output_dir <- as.character(output_dir[[1L]])
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figures <- NULL
  if (isTRUE(write_figures)) {
    plot_path <- file.path(output_dir, "priority_score_vs_kinship.png")
    ng_plot_priority_score_vs_kinship(
      scored = candidate_crosses,
      selected = selected_crosses,
      output_path = plot_path
    )
    files$priority_score_vs_kinship_png <- normalizePath(plot_path, winslash = "/", mustWork = TRUE)
    figures <- data.frame(
      figure = "priority_score_vs_kinship",
      path = files$priority_score_vs_kinship_png,
      stringsAsFactors = FALSE
    )
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
      n_crosses_requested = n_crosses
    )
  }
  files
}

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
                                    trait_value_metric = c("uc", "pmv", "vpm", "var_simple", "var_complex", "mean"),
                                    uc_variance_source = c("pmv", "vpm", "var_simple"),
                                    multi_trait_method = "auto",
                                    trait_weights = NULL,
                                    threshold_policy = c("soft", "strict"),
                                    threshold_penalty_weight = 1.0,
                                    threshold_penalty_autoscale = TRUE,
                                    progeny = "DH",
                                    recombination_model = c("haldane", "kosambi"),
                                    selection_prop = 0.10,
                                    min_effect_reliability = 0.35,
                                    grm_method = c("vanraden", "yang"),
                                    method_varPMV = c("fast", "full_posterior"),
                                    ril_mode = "infinite",
                                    run_posterior_prediction = FALSE,
                                    posterior_method = c("mcmc", "closed_form"),
                                    nIter = 5000L,
                                    burnIn = 500L,
                                    use_parallel = FALSE,
                                    parallel_cores = NULL,
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
                                    max_uses_per_parent = 6,
                                    min_unique_parents = NULL,
                                    max_pair_kinship = Inf,
                                    optimizer = "auto",
                                    allocation_method = "ocs",
                                    use_ocs = TRUE,
                                    lambda_group = 0.05,
                                    lambda_mating = 0.02,
                                    lambda_parent_use = 0,
                                    lambda_parent_use_mode = c("absolute", "adaptive"),
                                    lambda_progeny_inbreeding = 0,
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
                                    assume_inbred = TRUE,
                                    use_cpp = TRUE,
                                    seed = 1L) {
  prediction_mode <- match.arg(prediction_mode)
  map_position_unit <- match.arg(map_position_unit)
  trait_value_metric <- match.arg(trait_value_metric)
  uc_variance_source <- match.arg(uc_variance_source)
  threshold_policy <- match.arg(threshold_policy)
  recombination_model <- match.arg(recombination_model)
  grm_method <- match.arg(grm_method)
  ld_backend <- match.arg(ld_backend)
  method_varPMV <- ng_run_cp_method_varPMV(method_varPMV)
  ril_mode <- ng_run_cp_ril_mode(ril_mode)
  posterior_method <- match.arg(posterior_method)
  duplicate_action <- match.arg(duplicate_action)
  lambda_parent_use_mode <- match.arg(lambda_parent_use_mode)
  alphamate_mode <- match.arg(alphamate_mode)
  target <- ng_run_cp_target(progeny)
  optimizer_method <- ng_run_cp_optimizer(optimizer)
  allocation_method <- ng_run_cp_allocation_method(allocation_method)
  run_posterior_prediction <- ng_run_cp_logical(run_posterior_prediction, "run_posterior_prediction")
  use_parallel <- ng_run_cp_logical(use_parallel, "use_parallel")
  nIter <- ng_run_cp_integer(nIter, "nIter", min_value = 1L)
  burnIn <- ng_run_cp_integer(burnIn, "burnIn", min_value = 0L)
  posterior_n_draws <- max(1L, nIter - burnIn)
  n_crosses <- suppressWarnings(as.integer(n_crosses[[1L]]))
  if (!is.finite(n_crosses) || n_crosses < 1L) ng_stop("n_crosses must be a positive integer")

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
  if (any(qc$issues$severity == "blocker")) {
    blockers <- qc$issues$message[qc$issues$severity == "blocker"]
    ng_stop(
      "Input QC found blocker issues before cross prediction: ",
      paste(utils::head(blockers, 5L), collapse = " | ")
    )
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

  # Optional LD pruning of the marker matrix (same capability as ng_design_crosses): drop
  # redundant/low-MAF markers before scoring. Per-trait effects are fit from `geno` inside
  # the loop below, so pruning here keeps effects, geno, and the map aligned automatically.
  ld_pruning_report <- NULL
  if (isTRUE(ld_pruning)) {
    geno <- ng_ld_prune_geno(
      geno = geno, window = ld_window, r2_threshold = ld_r2_threshold,
      maf_threshold = ld_maf_threshold, ploidy = ld_ploidy, backend = ld_backend
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
  if (!is.null(training_set) && identical(trait_value_metric, "var_simple")) {
    warning("trait_value_metric = 'var_simple' does not use marker effects; the supplied ",
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
    fit <- ng_fit_ridge_effects(
      geno = fit_geno,
      y = fit_y,
      ids = fit_ids,
      seed = seed + i - 1L,
      return_beta_cov_full = identical(method_varPMV, "full_posterior")
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
      recomb_model = recombination_model,
      use_cpp = use_cpp,
      assume_inbred = assume_inbred,
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
        mcmc_burnin = burnIn,
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
        recomb_model = recombination_model,
        use_cpp = use_cpp,
        assume_inbred = assume_inbred,
        top_n_targets = unique(as.integer(c(min(n_crosses, 10L), 10L, 20L, 50L)))
      )
    }
    clean_trait <- ng_run_cp_clean_trait_name(trait)
    pmv_used_col <- ng_run_cp_pmv_col(scored_trait, method_varPMV)
    var_complex_col <- ng_run_cp_var_complex_col(scored_trait, method_varPMV)
    list(
      index = i,
      trait = trait,
      column = column,
      clean_trait = clean_trait,
      y = y,
      fit = fit,
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
        marker_effect_reliability = fit$reliability,
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
    parallel_cores = parallel_cores
  )
  parallel_backend <- attr(trait_results, "parallel_backend")
  parallel_cores_used <- attr(trait_results, "parallel_cores")
  trait_results <- trait_results[order(vapply(trait_results, `[[`, integer(1L), "index"))]
  effect_summary <- vector("list", length(trait_results))

  for (j in seq_along(trait_results)) {
    item <- trait_results[[j]]
    trait <- item$trait
    scored_trait <- item$scored_trait
    clean_trait <- item$clean_trait
    if (is.null(cross_table)) {
      cross_table <- scored_trait[, c("parent1", "parent2", "pair_kinship"), drop = FALSE]
    }
    cross_table[[paste0(clean_trait, "_value")]] <- item$value
    cross_table[[paste0(clean_trait, "_mean")]] <- scored_trait$cross_mean_blend
    cross_table[[paste0(clean_trait, "_pmv")]] <- scored_trait$dh_pmv_var
    cross_table[[paste0(clean_trait, "_pmv_fast")]] <- scored_trait$dh_pmv_var
    cross_table[[paste0(clean_trait, "_pmv_full_posterior")]] <- scored_trait$dh_pmv_var_full_posterior
    cross_table[[paste0(clean_trait, "_pmv_used")]] <- scored_trait[[item$pmv_used_col]]
    cross_table[[paste0(clean_trait, "_vpm")]] <- scored_trait$dh_recomb_var
    cross_table[[paste0(clean_trait, "_var_simple")]] <- scored_trait$var_simple
    cross_table[[paste0(clean_trait, "_var_complex")]] <- scored_trait[[item$var_complex_col]]
    cross_table[[paste0(clean_trait, "_reliability")]] <- scored_trait$effect_reliability
    effects_list[[trait]] <- item$fit
    trait_scores[[trait]] <- scored_trait
    if (isTRUE(run_posterior_prediction)) {
      posterior_effects_list[[trait]] <- item$posterior_effects
      posterior_predictions_list[[trait]] <- item$posterior_scores
    }
    effect_summary[[j]] <- item$effect_summary
  }

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
  if (!is.null(lethal_spec)) {
    cross_table <- ng_apply_marker_management(
      scores = cross_table, geno = geno, lethal_spec = lethal_spec,
      ploidy = marker_ploidy, drop_lethal_carrier_crosses = drop_lethal_carrier_crosses)
  }
  scored_crosses <- ng_score_breeder_objective(
    cross_table,
    objective,
    threshold_penalty_weight = threshold_penalty_weight,
    threshold_penalty_autoscale = threshold_penalty_autoscale
  )
  # Marker steering (Module 4): attach the target-allele score (and, when lambda_marker != 0,
  # the blended merit) to the reported candidate table and drive the criterion the native /
  # AlphaMate-style allocators optimize. The OCS path performs the same blend internally.
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
  parent_K <- if (isTRUE(use_ocs) || !identical(allocation_method, "ocs")) ng_parent_kinship(geno, method = grm_method) else NULL
  if (identical(allocation_method, "ocs")) {
    # The mate-selection module knobs (strategy / diversity_emphasis, progeny inbreeding,
    # breeder constraints, cost/logistics) are forwarded through ng_optimize_breeder_
    # selection_plan's ... into ng_optimize_mating_plan, so the user-friendly runner
    # reaches the same capabilities as the low-level API.
    plan <- ng_optimize_breeder_selection_plan(
      scores = cross_table,
      objective = objective,
      n_crosses = n_crosses,
      parent_K = parent_K,
      optimizer_method = optimizer_method,
      threshold_penalty_weight = threshold_penalty_weight,
      threshold_penalty_autoscale = threshold_penalty_autoscale,
      marker_target_spec = marker_target_spec,
      marker_geno = geno,
      lambda_marker = lambda_marker,
      marker_ploidy = marker_ploidy,
      max_crosses_per_parent = max_uses_per_parent,
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
    alpha_max_contrib <- if (is.null(alphamate_max_contributions)) max_uses_per_parent else alphamate_max_contributions
    plan <- ng_alphamate_style_select(
      scores = scored_crosses,
      criterion_col = allocation_criterion_col,
      n_crosses = n_crosses,
      parent_K = parent_K,
      mode = alphamate_mode,
      target_degree = alphamate_target_degree,
      max_contributions = alpha_max_contrib,
      lambda_group = alphamate_lambda_group,
      lambda_grid = alphamate_lambda_grid,
      method = optimizer_method,
      local_iter = local_iter
    )
  } else if (identical(allocation_method, "alphamate_executable")) {
    alpha_max_contrib <- if (is.null(alphamate_max_contributions)) max_uses_per_parent else alphamate_max_contributions
    alpha_executable <- if (is.null(alphamate_executable)) ng_alphamate_default_executable() else alphamate_executable
    alpha_workdir <- if (is.null(alphamate_workdir)) tempfile("ng_alphamate_run_") else alphamate_workdir
    plan <- ng_select_alphamate(
      scores = scored_crosses,
      criterion_col = allocation_criterion_col,
      n_crosses = n_crosses,
      parent_K = parent_K,
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

  output_files <- ng_run_cp_output_files(
    output_dir = output_dir,
    output_file = output_file,
    selected_crosses = selected,
    candidate_crosses = scored_crosses,
    trait_spec = trait_spec,
    qc = qc,
    write_outputs = write_outputs,
    write_figures = write_figures,
    n_crosses = n_crosses
  )

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
    output_files = output_files,
    settings = list(
      trait_value_metric = trait_value_metric,
      uc_variance_source = uc_variance_source,
      method_varPMV = method_varPMV,
      multi_trait_method = multi_trait_method,
      progeny = target,
      ril_mode = ril_mode,
      recombination_model = recombination_model,
      selection_prop = selection_prop,
      run_posterior_prediction = run_posterior_prediction,
      posterior_method = posterior_method,
      nIter = nIter,
      burnIn = burnIn,
      posterior_n_draws = if (isTRUE(run_posterior_prediction)) posterior_n_draws else 0L,
      use_parallel = use_parallel,
      parallel_backend = parallel_backend,
      parallel_cores = parallel_cores_used,
      optimizer = optimizer,
      optimizer_method = optimizer_method,
      allocation_method = allocation_method,
      alphamate_mode = alphamate_mode,
      alphamate_target_degree = alphamate_target_degree,
      alphamate_max_contributions = if (is.null(alphamate_max_contributions)) max_uses_per_parent else alphamate_max_contributions,
      use_ocs = isTRUE(use_ocs),
      n_crosses = n_crosses,
      max_uses_per_parent = max_uses_per_parent
    )
  )
  class(result) <- c("ng_cross_prediction_result", "list")
  result
}
