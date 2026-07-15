ng_poly4x_scenarios <- function() {
  data.frame(
    scenario = c("potato_autotetraploid_4x", "cassava_autotetraploid_4x"),
    crop = c("potato", "cassava"),
    description = c(
      "Autotetraploid potato 4x simulation with polysomic inheritance.",
      "Autotetraploid cassava 4x simulation with polysomic inheritance."
    ),
    ploidy = c(4L, 4L),
    n_chr = c(12L, 18L),
    genome_length_m = c(0.95, 1.20),
    seg_sites = c(1100L, 1200L),
    snp_per_chr = c(550L, 450L),
    qtl_per_chr = c(35L, 25L),
    n_founders = c(160L, 160L),
    phenotype_h2 = c(0.35, 0.35),
    quad_prob = c(0.10, 0.10),
    digenic_var = c(0.20, 0.20),
    stringsAsFactors = FALSE
  )
}

ng_poly4x_select_scenario <- function(scenarios = NULL) {
  available <- ng_poly4x_scenarios()
  if (is.null(scenarios) || !length(scenarios) ||
      (length(scenarios) == 1L && !nzchar(trimws(as.character(scenarios))))) {
    return(available)
  }
  requested <- unlist(strsplit(paste(as.character(scenarios), collapse = ","), ",", fixed = TRUE),
                      use.names = FALSE)
  requested <- unique(trimws(requested[nzchar(trimws(requested))]))
  missing <- setdiff(requested, available$scenario)
  if (length(missing)) {
    ng_stop("Unknown poly4x scenario(s): ", paste(missing, collapse = ", "))
  }
  available[match(requested, available$scenario), , drop = FALSE]
}

ng_poly4x_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    ng_stop("seed must be a single finite numeric value")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

ng_poly4x_scenario_number <- function(scenario, field, integer = FALSE, lower = -Inf,
                                      upper = Inf, lower_inclusive = TRUE,
                                      upper_inclusive = TRUE) {
  if (!field %in% names(scenario)) {
    ng_stop("scenario ", field, " must be provided")
  }
  value <- scenario[[field]]
  if (length(value) != 1L) {
    ng_stop("scenario ", field, " must be a single finite numeric value")
  }
  number <- suppressWarnings(as.numeric(value[[1]]))
  if (!is.finite(number)) {
    ng_stop("scenario ", field, " must be a single finite numeric value")
  }
  if (isTRUE(integer) && number != as.integer(number)) {
    ng_stop("scenario ", field, " must be an integer value")
  }
  lower_ok <- if (lower_inclusive) number >= lower else number > lower
  upper_ok <- if (upper_inclusive) number <= upper else number < upper
  if (!lower_ok || !upper_ok) {
    ng_stop("scenario ", field, " is outside the allowed range")
  }
  if (isTRUE(integer)) as.integer(number) else number
}

ng_poly4x_validate_parent_ids <- function(parent_pop) {
  ids <- as.character(parent_pop@id)
  duplicate_idx <- anyDuplicated(ids)
  if (duplicate_idx != 0L) {
    ng_stop("parent population ids must be unique; duplicate id: ", ids[[duplicate_idx]])
  }
  invisible(TRUE)
}

ng_poly4x_require_alphasimr <- function() {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    ng_stop("AlphaSimR is required for the poly4x runner")
  }
  invisible(TRUE)
}

ng_poly4x_setup_simparam <- function(scenario, seed = 1L, include_digenic = FALSE) {
  ng_poly4x_require_alphasimr()
  scenario <- as.data.frame(scenario, stringsAsFactors = FALSE)
  if (nrow(scenario) != 1L) ng_stop("scenario must contain exactly one row")
  ploidy <- ng_poly4x_scenario_number(scenario, "ploidy", integer = TRUE, lower = 1L)
  if (ploidy != 4L) ng_stop("scenario ploidy must be 4")
  n_founders <- ng_poly4x_scenario_number(scenario, "n_founders", integer = TRUE, lower = 1L)
  n_chr <- ng_poly4x_scenario_number(scenario, "n_chr", integer = TRUE, lower = 1L)
  seg_sites <- ng_poly4x_scenario_number(scenario, "seg_sites", integer = TRUE, lower = 1L)
  genome_length_m <- ng_poly4x_scenario_number(scenario, "genome_length_m", lower = 0, lower_inclusive = FALSE)
  quad_prob <- ng_poly4x_scenario_number(scenario, "quad_prob", lower = 0, upper = 1)
  qtl_per_chr <- ng_poly4x_scenario_number(scenario, "qtl_per_chr", integer = TRUE, lower = 1L)
  digenic_var <- ng_poly4x_scenario_number(scenario, "digenic_var", lower = 0)
  snp_per_chr <- ng_poly4x_scenario_number(scenario, "snp_per_chr", integer = TRUE, lower = 1L)

  ng_poly4x_with_seed(seed, {
    founder <- AlphaSimR::quickHaplo(
      nInd = n_founders,
      nChr = n_chr,
      segSites = seg_sites,
      genLen = genome_length_m,
      ploidy = 4L,
      inbred = FALSE
    )
    sim_param <- AlphaSimR::SimParam$new(founder)
    sim_param$quadProb <- quad_prob
    if (isTRUE(include_digenic)) {
      sim_param$addTraitAD(
        nQtlPerChr = qtl_per_chr,
        var = 1,
        meanDD = 0,
        varDD = digenic_var
      )
    } else {
      sim_param$addTraitA(nQtlPerChr = qtl_per_chr, var = 1)
    }
    sim_param$addSnpChip(nSnpPerChr = snp_per_chr)
    list(founder = founder, sim_param = sim_param, scenario = scenario)
  })
}

ng_poly4x_make_parent_pop <- function(setup, n_parents) {
  ng_poly4x_require_alphasimr()
  n_parents <- as.integer(n_parents)
  if (length(n_parents) != 1L || is.na(n_parents) || n_parents < 2L) {
    ng_stop("n_parents must be at least 2")
  }
  parent_pop <- AlphaSimR::newPop(setup$founder, simParam = setup$sim_param)
  if (AlphaSimR::nInd(parent_pop) < n_parents) ng_stop("not enough founders for requested parents")
  parent_pop <- parent_pop[seq_len(n_parents)]
  parent_pop@id <- paste0("poly4x_P", seq_len(n_parents))
  ng_poly4x_assert_pop4x(parent_pop, "parent population")
  parent_pop
}

ng_poly4x_assert_pop4x <- function(pop, context = "population") {
  ploidy <- unique(as.integer(pop@ploidy))
  if (!identical(ploidy, 4L)) ng_stop(context, " must have ploidy 4")
  invisible(TRUE)
}

ng_poly4x_pull_dosage <- function(pop, sim_param, marker_names = NULL) {
  ng_poly4x_require_alphasimr()
  ng_poly4x_assert_pop4x(pop)
  dosage <- AlphaSimR::pullSnpGeno(pop, simParam = sim_param)
  if (is.null(marker_names)) {
    marker_names <- paste0("snp", seq_len(ncol(dosage)))
  }
  if (length(marker_names) != ncol(dosage)) ng_stop("marker_names length must match SNP count")
  colnames(dosage) <- as.character(marker_names)
  rownames(dosage) <- pop@id
  ng_poly4x_as_dosage_matrix(dosage, ploidy = 4L, name = "poly4x dosage")
}

ng_poly4x_parent_index <- function(parent_pop, parent_id) {
  ng_poly4x_validate_parent_ids(parent_pop)
  parent_id <- as.character(parent_id)
  if (length(parent_id) != 1L || is.na(parent_id)) ng_stop("parent_id must be a single value")
  idx <- match(parent_id, parent_pop@id)
  if (is.na(idx)) ng_stop("parent not found: ", parent_id)
  idx
}

ng_poly4x_make_family <- function(parent_pop, parent1, parent2, n_progeny, sim_param, seed = NULL) {
  ng_poly4x_require_alphasimr()
  ng_poly4x_assert_pop4x(parent_pop, "parent population")
  parent1 <- as.character(parent1)
  parent2 <- as.character(parent2)
  if (identical(parent1, parent2)) ng_stop("self-cross is not supported for poly4x families")
  n_progeny <- as.integer(n_progeny)
  if (length(n_progeny) != 1L || is.na(n_progeny) || n_progeny < 1L) {
    ng_stop("n_progeny must be at least 1")
  }
  idx1 <- ng_poly4x_parent_index(parent_pop, parent1)
  idx2 <- ng_poly4x_parent_index(parent_pop, parent2)
  family <- ng_poly4x_with_seed(seed, {
    AlphaSimR::makeCross(
      parent_pop,
      crossPlan = matrix(c(idx1, idx2), nrow = 1L),
      nProgeny = n_progeny,
      simParam = sim_param
    )
  })
  family@id <- paste0(parent1, "x", parent2, "_", seq_len(AlphaSimR::nInd(family)))
  ng_poly4x_assert_pop4x(family, "family")
  family
}

ng_poly4x_family_stats <- function(pop, top_prop = 0.10) {
  ng_poly4x_require_alphasimr()
  gv <- as.numeric(AlphaSimR::gv(pop)[, 1])
  if (!is.numeric(top_prop) || length(top_prop) != 1L ||
      !is.finite(top_prop) || top_prop <= 0 || top_prop > 1) {
    ng_stop("top_prop must be a single finite numeric value in (0, 1]")
  }
  top_n <- min(length(gv), max(1L, ceiling(length(gv) * top_prop)))
  data.frame(
    mean_gv = mean(gv),
    var_gv = if (length(gv) > 1L) stats::var(gv) else 0,
    top10_gv = mean(sort(gv, decreasing = TRUE)[seq_len(top_n)]),
    max_gv = max(gv),
    stringsAsFactors = FALSE
  )
}
