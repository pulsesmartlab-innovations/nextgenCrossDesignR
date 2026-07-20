ng_multitrait_crop_validation_int <- function(value, default, name, min_value = 1L) {
  if (is.null(value)) value <- default
  value <- suppressWarnings(as.integer(value[[1]]))
  if (!is.finite(value) || value < min_value) {
    ng_stop(name, " must be an integer >= ", min_value)
  }
  value
}

ng_multitrait_crop_validation_num <- function(value, default, name, min_value = NULL) {
  if (is.null(value)) value <- default
  value <- suppressWarnings(as.numeric(value[[1]]))
  if (!is.finite(value)) ng_stop(name, " must be finite")
  if (!is.null(min_value) && value < min_value) ng_stop(name, " must be >= ", min_value)
  value
}

ng_multitrait_crop_validation_seed <- function(seed, scenario, n_parents, rep, offset = 0L) {
  seed <- suppressWarnings(as.integer(seed[[1]]))
  if (!is.finite(seed)) seed <- 1L
  key <- utf8ToInt(paste(as.character(scenario), collapse = "|"))
  hash <- 0L
  if (length(key)) {
    hash <- sum((seq_along(key) + 17L) * key) %% 900000L
  }
  out <- as.numeric(seed) + as.numeric(hash) + as.numeric(n_parents) * 1000 + as.numeric(rep) + as.numeric(offset)
  out <- ((out - 1) %% 2147483000) + 1
  as.integer(out)
}

ng_multitrait_crop_make_dh_parents <- function(founder_pop, n_parents, sim_param) {
  n_parents <- suppressWarnings(as.integer(n_parents[[1]]))
  if (!is.finite(n_parents) || n_parents < 1L) ng_stop("n_parents must be positive")
  n_dh_per_founder <- ceiling(n_parents / AlphaSimR::nInd(founder_pop))
  dh <- AlphaSimR::makeDH(
    founder_pop,
    nDH = n_dh_per_founder,
    keepParents = FALSE,
    simParam = sim_param
  )
  dh[seq_len(n_parents)]
}

ng_multitrait_crop_trait_profile <- function(scenario = "compact_selfing") {
  crop_scenario <- ng_crop_genome_select(scenario)
  if (nrow(crop_scenario) != 1L) ng_stop("scenario must select exactly one crop scenario")
  crop <- as.character(crop_scenario$crop[[1]])
  scenario_name <- as.character(crop_scenario$scenario[[1]])

  if (crop %in% c("cassava", "potato")) {
    weight <- c(1.5, 4.0, 2.5)
    economic_weight <- c(1.5, 4.0, 2.5)
    desired_change <- c(4.0, 25.0, 3.0)
    corA <- matrix(c(
      1.00, 0.35, -0.45,
      0.35, 1.00, -0.25,
      -0.45, -0.25, 1.00
    ), nrow = 3L, byrow = TRUE)
    quality_label <- "dry_matter_proxy"
  } else if (crop %in% c("wheat", "barley", "field_pea")) {
    weight <- c(1.5, 3.5, 2.0)
    economic_weight <- c(1.5, 3.5, 2.0)
    desired_change <- c(3.5, 22.0, 2.5)
    corA <- matrix(c(
      1.00, 0.40, -0.30,
      0.40, 1.00, -0.25,
      -0.30, -0.25, 1.00
    ), nrow = 3L, byrow = TRUE)
    quality_label <- "grain_quality_proxy"
  } else if (crop %in% c("sugarcane")) {
    weight <- c(1.5, 3.0, 3.0)
    economic_weight <- c(1.5, 3.0, 3.0)
    desired_change <- c(5.0, 20.0, 3.5)
    corA <- matrix(c(
      1.00, 0.30, -0.35,
      0.30, 1.00, -0.20,
      -0.35, -0.20, 1.00
    ), nrow = 3L, byrow = TRUE)
    quality_label <- "sucrose_quality_proxy"
  } else {
    weight <- c(1.0, 4.0, 2.0)
    economic_weight <- c(1.0, 4.0, 2.0)
    desired_change <- c(4.0, 25.0, 2.0)
    corA <- matrix(c(
      1.00, 0.45, -0.25,
      0.45, 1.00, -0.35,
      -0.25, -0.35, 1.00
    ), nrow = 3L, byrow = TRUE)
    quality_label <- "quality_proxy"
  }

  trait <- c("yield", "disease", "quality")
  traits <- ng_multitrait_spec(
    trait = trait,
    column = paste0("pred_", trait),
    direction = c("maximize", "minimize", "maximize"),
    weight = weight,
    economic_weight = economic_weight,
    desired_change = desired_change,
    max_value = c(NA, 40, NA),
    threshold_weight = c(1, 4, 1)
  )
  rownames(corA) <- colnames(corA) <- trait
  list(
    scenario = scenario_name,
    crop = crop,
    quality_label = quality_label,
    traits = traits,
    corA = corA,
    mean = stats::setNames(c(100, 40, 35), trait),
    var = stats::setNames(c(64, 100, 9), trait)
  )
}

ng_multitrait_crop_validation_scenario <- function(scenario = "compact_selfing",
                                                   n_parents = 20L,
                                                   n_founders = NULL,
                                                   n_chr = NULL,
                                                   seg_sites = NULL,
                                                   snp_per_chr = NULL,
                                                   qtl_per_chr = NULL,
                                                   genome_length_m = NULL,
                                                   realized_progeny = 20L,
                                                   seed = 1L,
                                                   prediction_noise = 0.20,
                                                   alphasimr_threads = 1L) {
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    ng_stop("AlphaSimR is required for multi-trait crop validation")
  }
  crop_scenario <- ng_crop_genome_select(scenario)
  if (nrow(crop_scenario) != 1L) ng_stop("scenario must select exactly one crop scenario")
  profile <- ng_multitrait_crop_trait_profile(crop_scenario$scenario[[1]])
  traits <- profile$traits
  realized_cols <- ng_multitrait_validation_realized_cols(traits)

  n_parents <- ng_multitrait_crop_validation_int(n_parents, 20L, "n_parents", min_value = 4L)
  n_founders <- ng_multitrait_crop_validation_int(
    n_founders, max(as.integer(crop_scenario$n_founders[[1]]), n_parents), "n_founders", min_value = n_parents
  )
  n_chr <- ng_multitrait_crop_validation_int(n_chr, crop_scenario$n_chr[[1]], "n_chr", min_value = 1L)
  seg_sites <- ng_multitrait_crop_validation_int(seg_sites, crop_scenario$seg_sites[[1]], "seg_sites", min_value = 4L)
  snp_per_chr <- ng_multitrait_crop_validation_int(snp_per_chr, crop_scenario$snp_per_chr[[1]], "snp_per_chr")
  qtl_per_chr <- ng_multitrait_crop_validation_int(qtl_per_chr, crop_scenario$qtl_per_chr[[1]], "qtl_per_chr")
  genome_length_m <- ng_multitrait_crop_validation_num(
    genome_length_m, crop_scenario$genome_length_m[[1]], "genome_length_m", min_value = .Machine$double.eps
  )
  realized_progeny <- ng_multitrait_crop_validation_int(realized_progeny, 20L, "realized_progeny")
  seed <- ng_multitrait_crop_validation_int(seed, 1L, "seed", min_value = 1L)
  alphasimr_threads <- ng_multitrait_crop_validation_int(alphasimr_threads, 1L, "alphasimr_threads")
  prediction_noise <- ng_multitrait_crop_validation_num(prediction_noise, 0.20, "prediction_noise", min_value = 0)

  set.seed(seed)
  founder <- AlphaSimR::quickHaplo(n_founders, n_chr, seg_sites, genLen = genome_length_m)
  sim_param <- AlphaSimR::SimParam$new(founder)
  sim_param$nThreads <- max(1L, alphasimr_threads)
  sim_param$addTraitA(
    nQtlPerChr = qtl_per_chr,
    mean = as.numeric(profile$mean),
    var = as.numeric(profile$var),
    corA = profile$corA,
    name = traits$trait
  )
  sim_param$addSnpChip(nSnpPerChr = snp_per_chr)

  base_pop <- AlphaSimR::newPop(founder, simParam = sim_param)
  parent_pop <- ng_multitrait_crop_make_dh_parents(base_pop, n_parents, sim_param)
  parent_ids <- sprintf("P%02d", seq_len(n_parents))
  parent_pop@id <- parent_ids

  parent_gv <- as.matrix(AlphaSimR::gv(parent_pop))
  colnames(parent_gv) <- traits$trait
  rownames(parent_gv) <- parent_ids
  geno <- AlphaSimR::pullSnpGeno(parent_pop, simParam = sim_param)
  rownames(geno) <- parent_ids
  if (is.null(colnames(geno))) colnames(geno) <- sprintf("M%05d", seq_len(ncol(geno)))
  parent_kinship <- ng_parent_kinship(geno)

  pairs <- ng_make_pairs(parent_ids, include_self = FALSE)
  p1 <- match(pairs$parent1, parent_ids)
  p2 <- match(pairs$parent2, parent_ids)
  scores <- pairs
  for (trait in traits$trait) {
    mean_value <- 0.5 * (parent_gv[p1, trait] + parent_gv[p2, trait])
    scores[[paste0("pred_", trait)]] <- mean_value +
      stats::rnorm(nrow(scores), sd = sqrt(profile$var[[trait]]) * prediction_noise)
  }

  realized <- matrix(NA_real_, nrow = nrow(scores), ncol = nrow(traits))
  colnames(realized) <- traits$trait
  for (i in seq_len(nrow(scores))) {
    set.seed(ng_multitrait_crop_validation_seed(
      seed, paste0(profile$scenario, "|", scores$parent1[[i]], "|", scores$parent2[[i]]),
      realized_progeny, i, offset = 100000L
    ))
    f1 <- AlphaSimR::makeCross(
      parent_pop,
      matrix(c(p1[[i]], p2[[i]]), ncol = 2L),
      nProgeny = 1L,
      simParam = sim_param
    )
    dh <- AlphaSimR::makeDH(f1, nDH = realized_progeny, keepParents = FALSE, simParam = sim_param)
    family_gv <- as.matrix(AlphaSimR::gv(dh))
    colnames(family_gv) <- traits$trait
    realized[i, ] <- colMeans(family_gv)
  }
  for (trait in traits$trait) {
    scores[[realized_cols[[trait]]]] <- realized[, trait]
  }
  rel <- ng_pair_relationship_variance(scores[, c("parent1", "parent2")], parent_kinship)
  scores$pair_kinship <- rel$pair_kinship

  parent_values <- data.frame(parent = parent_ids, parent_gv, stringsAsFactors = FALSE)
  list(
    scores = scores,
    traits = traits,
    realized_cols = realized_cols,
    parent_values = parent_values,
    parent_genotype = geno,
    parent_kinship = parent_kinship,
    trait_profile = profile,
    crop_scenario = crop_scenario,
    config = data.frame(
      scenario = as.character(crop_scenario$scenario[[1]]),
      crop = as.character(crop_scenario$crop[[1]]),
      harness_model = as.character(crop_scenario$harness_model[[1]]),
      parent_generation_model = "DH",
      n_parents = n_parents,
      n_founders = n_founders,
      n_chr = n_chr,
      seg_sites = seg_sites,
      snp_per_chr = snp_per_chr,
      qtl_per_chr = qtl_per_chr,
      genome_length_m = genome_length_m,
      realized_progeny = realized_progeny,
      prediction_noise = prediction_noise,
      alphasimr_threads = alphasimr_threads,
      seed = seed,
      stringsAsFactors = FALSE
    )
  )
}

ng_multitrait_crop_validation_winner_summary <- function(summary) {
  summary <- as.data.frame(summary, stringsAsFactors = FALSE)
  required <- c("scenario", "crop", "harness_model", "n_parents", "rep")
  missing <- setdiff(required, names(summary))
  if (length(missing)) ng_stop("summary missing columns: ", paste(missing, collapse = ", "))
  keys <- unique(summary[, c("scenario", "crop", "harness_model", "n_parents"), drop = FALSE])
  out <- list()
  for (i in seq_len(nrow(keys))) {
    key <- keys[i, , drop = FALSE]
    rows <- summary$scenario == key$scenario[[1]] &
      summary$n_parents == key$n_parents[[1]]
    part <- summary[rows, , drop = FALSE]
    winners <- ng_multitrait_validation_winner_summary(part)
    if (!nrow(winners)) next
    winners$scenario <- key$scenario[[1]]
    winners$crop <- key$crop[[1]]
    winners$harness_model <- key$harness_model[[1]]
    out[[length(out) + 1L]] <- winners[, c(
      "scenario", "crop", "harness_model", "n_parents",
      "metric", "direction", "method", "tied_methods", "tied_method_count",
      "value", "reps"
    ), drop = FALSE]
  }
  if (!length(out)) {
    return(data.frame(
      scenario = character(),
      crop = character(),
      harness_model = character(),
      n_parents = integer(),
      metric = character(),
      direction = character(),
      method = character(),
      tied_methods = character(),
      tied_method_count = integer(),
      value = numeric(),
      reps = integer(),
      stringsAsFactors = FALSE
    ))
  }
  combined <- do.call(rbind, out)
  rownames(combined) <- NULL
  combined
}

ng_run_multitrait_crop_validation_grid <- function(scenarios = c("compact_selfing", "maize_like", "cassava_diploid"),
                                                   parent_sizes = c(12L, 20L),
                                                   reps = 1L,
                                                   n_crosses = 4L,
                                                   realized_progeny = 8L,
                                                   seed = 1L,
                                                   methods = ng_multitrait_validation_default_methods(),
                                                   allocator = c("topn", "ocs"),
                                                   n_founders = NULL,
                                                   n_chr = NULL,
                                                   seg_sites = NULL,
                                                   snp_per_chr = NULL,
                                                   qtl_per_chr = NULL,
                                                   genome_length_m = NULL,
                                                   prediction_noise = 0.20,
                                                   alphasimr_threads = 1L,
                                                   ocs_lambda_group = 0.05,
                                                   ocs_lambda_mating = 0,
                                                   output_dir = NULL,
                                                   prefix = "multitrait_crop_grid") {
  allocator <- match.arg(allocator)
  crop_scenarios <- ng_crop_genome_select(scenarios)
  parent_sizes <- unique(suppressWarnings(as.integer(parent_sizes)))
  parent_sizes <- parent_sizes[is.finite(parent_sizes)]
  if (!length(parent_sizes)) ng_stop("parent_sizes must contain at least one size")
  if (any(parent_sizes < 4L)) ng_stop("parent_sizes must all be at least 4")
  reps <- ng_multitrait_crop_validation_int(reps, 1L, "reps")
  seed <- ng_multitrait_crop_validation_int(seed, 1L, "seed")

  summaries <- list()
  selections <- list()
  scores_out <- list()
  configs <- list()
  k <- 0L
  for (s in seq_len(nrow(crop_scenarios))) {
    scenario_name <- as.character(crop_scenarios$scenario[[s]])
    for (n_parents in parent_sizes) {
      for (rep in seq_len(reps)) {
        run_seed <- ng_multitrait_crop_validation_seed(seed, scenario_name, n_parents, rep)
        scenario_obj <- ng_multitrait_crop_validation_scenario(
          scenario = scenario_name,
          n_parents = n_parents,
          n_founders = n_founders,
          n_chr = n_chr,
          seg_sites = seg_sites,
          snp_per_chr = snp_per_chr,
          qtl_per_chr = qtl_per_chr,
          genome_length_m = genome_length_m,
          realized_progeny = realized_progeny,
          seed = run_seed,
          prediction_noise = prediction_noise,
          alphasimr_threads = alphasimr_threads
        )
        run <- ng_run_multitrait_validation(
          scores = scenario_obj$scores,
          traits = scenario_obj$traits,
          realized_cols = scenario_obj$realized_cols,
          n_crosses = n_crosses,
          seed = run_seed,
          methods = methods,
          allocator = allocator,
          parent_kinship = scenario_obj$parent_kinship,
          ocs_lambda_group = ocs_lambda_group,
          ocs_lambda_mating = ocs_lambda_mating
        )
        k <- k + 1L
        run$summary$scenario <- scenario_name
        run$summary$crop <- as.character(crop_scenarios$crop[[s]])
        run$summary$harness_model <- as.character(crop_scenarios$harness_model[[s]])
        run$summary$n_parents <- n_parents
        run$summary$rep <- rep
        run$summary$seed <- run_seed
        run$summary$allocator <- allocator
        run$summary$realized_progeny <- realized_progeny
        run$selections$scenario <- scenario_name
        run$selections$crop <- as.character(crop_scenarios$crop[[s]])
        run$selections$harness_model <- as.character(crop_scenarios$harness_model[[s]])
        run$selections$n_parents <- n_parents
        run$selections$rep <- rep
        run$selections$seed <- run_seed
        run$selections$allocator <- allocator
        run$selections$realized_progeny <- realized_progeny
        summaries[[k]] <- run$summary
        selections[[k]] <- run$selections
        score_rows <- scenario_obj$scores
        score_rows$scenario <- scenario_name
        score_rows$crop <- as.character(crop_scenarios$crop[[s]])
        score_rows$harness_model <- as.character(crop_scenarios$harness_model[[s]])
        score_rows$n_parents <- n_parents
        score_rows$rep <- rep
        score_rows$seed <- run_seed
        score_rows$realized_progeny <- realized_progeny
        scores_out[[k]] <- score_rows
        cfg <- scenario_obj$config
        cfg$rep <- rep
        cfg$allocator <- allocator
        cfg$n_crosses <- n_crosses
        cfg$methods <- paste(methods, collapse = ",")
        cfg$ocs_lambda_group <- ocs_lambda_group
        cfg$ocs_lambda_mating <- ocs_lambda_mating
        configs[[k]] <- cfg
      }
    }
  }

  summary <- do.call(rbind, summaries)
  selections <- do.call(rbind, selections)
  scores <- do.call(rbind, scores_out)
  config <- do.call(rbind, configs)
  rownames(summary) <- NULL
  rownames(selections) <- NULL
  rownames(scores) <- NULL
  rownames(config) <- NULL
  winner_summary <- ng_multitrait_crop_validation_winner_summary(summary)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(summary, file.path(output_dir, paste0(prefix, "_summary.csv")), row.names = FALSE)
    write.csv(selections, file.path(output_dir, paste0(prefix, "_selections.csv")), row.names = FALSE)
    write.csv(scores, file.path(output_dir, paste0(prefix, "_scores.csv")), row.names = FALSE)
    write.csv(winner_summary, file.path(output_dir, paste0(prefix, "_winner_summary.csv")), row.names = FALSE)
    write.csv(config, file.path(output_dir, paste0(prefix, "_config.csv")), row.names = FALSE)
  }

  list(
    summary = summary,
    selections = selections,
    scores = scores,
    winner_summary = winner_summary,
    config = config,
    scenarios = crop_scenarios,
    parent_sizes = parent_sizes,
    reps = reps,
    n_crosses = n_crosses,
    realized_progeny = realized_progeny,
    seed = seed,
    allocator = allocator
  )
}
