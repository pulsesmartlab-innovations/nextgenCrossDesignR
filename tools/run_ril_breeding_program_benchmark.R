# Reproducible RIL breeding-program benchmark for nextgenCrossDesign.
#
# Purpose:
#   Compare nextgenCrossDesign against known breeding/genomic mating software in
#   a simulation loop that looks like a breeder-run RIL program: parents ->
#   planned crosses -> F2/F3/F4/F5/F6 advancement -> selected RIL parents for
#   the next cycle. This is intentionally different from a one-shot score-table
#   benchmark.
#
# Default run is review/smoke scale. Increase NG_RIL_* settings for publication
# scale. PopVar and SimpleMating are required exact comparison branches for the
# default benchmark and are installed from GitHub when missing. The script
# records package/binary availability and exact_external_status fields so setup
# failures are not confused with style-proxy baselines.
#
# Multi-trait design:
#   - nextgenCrossDesign methods use trait_by_trait scoring with a trait
#     direction file equivalent: yield/protein increase, disease decrease by
#     default, plus weights/economic weights/desired changes.
#   - PopVar and SimpleMating single-index comparison branches receive an
#     already-oriented selection index with index_direction = increase, because
#     those tools do not model multiple trait directions directly.

find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    start,
    file.path(start, "nextgen_cross_design"),
    dirname(start),
    file.path(dirname(start), "nextgen_cross_design"),
    file.path(dirname(dirname(start)), "nextgen_cross_design")
  ), winslash = "/", mustWork = FALSE))

  for (candidate in candidates) {
    desc <- file.path(candidate, "DESCRIPTION")
    if (file.exists(file.path(candidate, "R", "load.R")) &&
        file.exists(desc) &&
        any(grepl("^Package:\\s*nextgenCrossDesign\\s*$", readLines(desc, n = 20L, warn = FALSE)))) {
      return(candidate)
    }
  }
  stop("Could not locate nextgenCrossDesign root containing DESCRIPTION and R/load.R", call. = FALSE)
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (!is.finite(out)) as.integer(default) else out
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) as.numeric(default) else out
}

env_bool <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(value)) return(isTRUE(default))
  value %in% c("1", "true", "yes", "y")
}

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  x[nzchar(x) & !is.na(x)]
}

env_csv <- function(name, default) {
  split_csv(env_chr(name, default))
}

env_num_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(split_csv(value)))
  out <- out[is.finite(out)]
  if (!length(out)) as.numeric(default) else out
}

setup_benchmark_library <- function(cfg) {
  lib <- normalizePath(cfg$r_lib, winslash = "/", mustWork = FALSE)
  if (!dir.exists(lib)) dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(lib)) stop("Could not create R library for benchmark dependencies: ", lib, call. = FALSE)
  .libPaths(unique(c(lib, .libPaths())))
  lib
}

package_record <- function(package,
                           source,
                           repository = NA_character_,
                           ref = NA_character_,
                           subdir = NA_character_,
                           status,
                           reason = "") {
  version <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
  path <- if (requireNamespace(package, quietly = TRUE)) {
    normalizePath(find.package(package), winslash = "/", mustWork = FALSE)
  } else {
    NA_character_
  }
  data.frame(
    package = package,
    source = source,
    repository = repository,
    ref = ref,
    subdir = subdir,
    status = status,
    version = version,
    path = path,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

ensure_r_cran_package <- function(package, cfg) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(package_record(package, "CRAN", status = "already_installed"))
  }
  if (!isTRUE(cfg$install_missing_dependencies)) {
    stop(package, " is required but is not installed. Set NG_RIL_INSTALL_MISSING_DEPS=TRUE to install it.", call. = FALSE)
  }
  lib <- setup_benchmark_library(cfg)
  message("Installing required benchmark package from CRAN: ", package)
  tryCatch(
    utils::install.packages(package, lib = lib, repos = cfg$cran_repos, dependencies = TRUE),
    error = function(e) {
      stop("Could not install ", package, " from CRAN: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Installed ", package, " but it is still not loadable from .libPaths().", call. = FALSE)
  }
  package_record(package, "CRAN", status = "installed")
}

ensure_r_github_package <- function(package,
                                    repo,
                                    cfg,
                                    ref = "",
                                    subdir = "") {
  repo <- trimws(as.character(repo))
  ref <- trimws(as.character(ref))
  subdir <- trimws(as.character(subdir))
  if (!nzchar(repo)) stop("Missing GitHub repository for required package: ", package, call. = FALSE)
  if (requireNamespace(package, quietly = TRUE) && !isTRUE(cfg$force_github_install)) {
    return(package_record(
      package = package,
      source = "GitHub",
      repository = repo,
      ref = if (nzchar(ref)) ref else "HEAD",
      subdir = if (nzchar(subdir)) subdir else NA_character_,
      status = "already_installed"
    ))
  }
  if (!isTRUE(cfg$install_missing_dependencies) && !isTRUE(cfg$force_github_install)) {
    stop(package, " is required but is not installed. Set NG_RIL_INSTALL_MISSING_DEPS=TRUE to install it.", call. = FALSE)
  }
  lib <- setup_benchmark_library(cfg)
  if (!requireNamespace("remotes", quietly = TRUE)) {
    message("Installing remotes from CRAN so GitHub benchmark dependencies can be installed")
    utils::install.packages("remotes", lib = lib, repos = cfg$cran_repos, dependencies = TRUE)
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("The remotes package is required to install ", package, " from GitHub.", call. = FALSE)
  }
  message("Installing required benchmark package from GitHub: ", package, " <- ", repo)
  args <- list(
    repo = repo,
    ref = if (nzchar(ref)) ref else "HEAD",
    subdir = if (nzchar(subdir)) subdir else NULL,
    lib = lib,
    dependencies = TRUE,
    upgrade = cfg$github_upgrade,
    build_vignettes = FALSE,
    quiet = FALSE
  )
  tryCatch(
    do.call(remotes::install_github, args),
    error = function(e) {
      stop("Could not install ", package, " from GitHub repository ", repo, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Installed ", package, " from GitHub but it is still not loadable from .libPaths().", call. = FALSE)
  }
  package_record(
    package = package,
    source = "GitHub",
    repository = repo,
    ref = if (nzchar(ref)) ref else "HEAD",
    subdir = if (nzchar(subdir)) subdir else NA_character_,
    status = "installed"
  )
}

benchmark_methods_need_tool <- function(methods, tool) {
  registry <- ril_benchmark_method_registry()
  rule <- registry[match(methods, registry$method), , drop = FALSE]
  any(stats::na.omit(rule$external_tool) == tool)
}

install_required_benchmark_dependencies <- function(cfg) {
  setup_benchmark_library(cfg)
  # AlphaSimR is required (it provides the ground-truth simulation). The external
  # comparators PopVar/SimpleMating are OPTIONAL: if their GitHub build fails (e.g. a
  # missing transitive dependency such as sommer), warn and record the failure so
  # method_availability() skips those methods -- the benchmark must not halt because
  # an optional comparator will not compile on this machine.
  records <- list(ensure_r_cran_package("AlphaSimR", cfg))
  install_optional_external <- function(package, repo, ref, subdir) {
    tryCatch(
      ensure_r_github_package(package = package, repo = repo, ref = ref, subdir = subdir, cfg = cfg),
      error = function(e) {
        warning("Optional external comparator ", package, " is unavailable and will be skipped: ",
                conditionMessage(e), call. = FALSE)
        package_record(
          package = package, source = "GitHub", repository = repo,
          ref = if (nzchar(ref)) ref else "HEAD",
          subdir = if (nzchar(subdir)) subdir else NA_character_,
          status = "install_failed"
        )
      }
    )
  }
  if (benchmark_methods_need_tool(cfg$methods, "PopVar")) {
    records[[length(records) + 1L]] <- install_optional_external(
      "PopVar", cfg$popvar_github_repo, cfg$popvar_github_ref, cfg$popvar_github_subdir
    )
  }
  if (benchmark_methods_need_tool(cfg$methods, "SimpleMating")) {
    records[[length(records) + 1L]] <- install_optional_external(
      "SimpleMating", cfg$simplemating_github_repo, cfg$simplemating_github_ref, cfg$simplemating_github_subdir
    )
  }
  bind_rows_fill(records)
}

bind_rows_fill <- function(x) {
  x <- Filter(function(z) !is.null(z) && is.data.frame(z) && nrow(z) > 0L, x)
  if (!length(x)) return(data.frame())
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    for (m in missing) z[[m]] <- NA
    z[, cols, drop = FALSE]
  })
  out <- do.call(rbind, x)
  rownames(out) <- NULL
  out
}

with_rng_seed <- function(seed, expr) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  if (!is.null(seed) && is.finite(seed)) set.seed(as.integer(seed))
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

stable_seed <- function(..., seed = 1L) {
  txt <- paste(..., collapse = "|")
  bytes <- utf8ToInt(txt)
  h <- as.numeric(seed %% 2147483647L)
  for (b in bytes) h <- (h * 131 + b) %% 2147483647
  as.integer(max(1, h))
}

top_prop_mean <- function(x, prop = 0.10) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  n <- max(1L, ceiling(length(x) * prop))
  mean(utils::head(sort(x, decreasing = TRUE), n))
}

resolve_benchmark_trait_spec <- function(cfg) {
  traits <- as.character(cfg$trait_names)
  n <- length(traits)
  if (!n) stop("NG_RIL_TRAITS must define at least one trait", call. = FALSE)
  recycle_to_traits <- function(x, default, name) {
    if (is.null(x) || !length(x)) x <- default
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) stop(name, " must have length 1 or match NG_RIL_TRAITS", call. = FALSE)
    x
  }
  directions <- recycle_to_traits(cfg$trait_directions, "increase", "NG_RIL_TRAIT_DIRECTIONS")
  directions <- tolower(trimws(directions))
  direction_sign <- rep(NA_real_, n)
  direction_sign[directions %in% c("increase", "max", "maximize", "maximise", "higher", "+")] <- 1
  direction_sign[directions %in% c("decrease", "min", "minimize", "minimise", "lower", "-")] <- -1
  if (any(!is.finite(direction_sign))) {
    stop("NG_RIL_TRAIT_DIRECTIONS must use increase/decrease style values", call. = FALSE)
  }
  directions <- ifelse(direction_sign > 0, "increase", "decrease")
  weights <- suppressWarnings(as.numeric(recycle_to_traits(cfg$trait_weights, 1, "NG_RIL_TRAIT_WEIGHTS")))
  weights[!is.finite(weights) | weights < 0] <- 0
  if (!any(weights > 0)) weights <- rep(1, n)
  weights <- weights / sum(weights)
  economic <- suppressWarnings(as.numeric(recycle_to_traits(
    cfg$trait_economic_weights, weights, "NG_RIL_TRAIT_ECONOMIC_WEIGHTS"
  )))
  economic[!is.finite(economic) | economic < 0] <- NA_real_
  desired <- suppressWarnings(as.numeric(recycle_to_traits(
    cfg$trait_desired_change, NA_real_, "NG_RIL_TRAIT_DESIRED_CHANGE"
  )))
  variances <- suppressWarnings(as.numeric(recycle_to_traits(
    cfg$trait_variances, 1, "NG_RIL_TRAIT_VARIANCES"
  )))
  variances[!is.finite(variances) | variances <= 0] <- 1
  data.frame(
    trait = traits,
    column = traits,
    direction = directions,
    Selection_direction = directions,
    direction_sign = direction_sign,
    weight = weights,
    economic_weight = economic,
    desired_change = desired,
    simulated_variance = variances,
    stringsAsFactors = FALSE
  )
}

gv_matrix <- function(pop, cfg) {
  g <- AlphaSimR::gv(pop)
  if (is.null(dim(g))) g <- matrix(as.numeric(g), ncol = 1L)
  g <- as.matrix(g)
  n_traits <- nrow(cfg$trait_spec)
  if (ncol(g) < n_traits) stop("AlphaSimR returned fewer traits than configured", call. = FALSE)
  g <- g[, seq_len(n_traits), drop = FALSE]
  colnames(g) <- cfg$trait_spec$trait
  g
}

# The aggregate breeding-goal index H = sum_t weight_t * direction_sign_t * gv_t on
# TRUE genetic values (normalized weights; minimize traits sign-flipped so higher H is
# better). This fixed LINEAR index is the merit ALL methods are graded on. NOTE: the
# gain comparison is objective-aligned only when a method's own aggregation is also
# this linear index (the default multi_trait_method="weighted"). Under non-linear
# modes (economic_index / desired_gain / threshold) the graded merit (linear H) is not
# exactly the composite a method optimizes, so those gain deltas mix index choice with
# method quality -- interpret accordingly.
selection_index_from_trait_values <- function(values, cfg) {
  values <- as.matrix(values)
  if (ncol(values) < nrow(cfg$trait_spec)) {
    stop("trait value matrix has fewer columns than the benchmark trait spec", call. = FALSE)
  }
  values <- values[, cfg$trait_spec$trait, drop = FALSE]
  oriented <- sweep(values, 2, cfg$trait_spec$direction_sign, `*`)
  as.numeric(oriented %*% cfg$trait_spec$weight)
}

gv_vec <- function(pop, cfg) {
  selection_index_from_trait_values(gv_matrix(pop, cfg), cfg)
}

# Per-line composite index built with the package's OWN multi-trait objective, so
# single-trait external tools (PopVar/SimpleMating) can be trained on EXACTLY the
# index the package optimizes in trait_by_trait mode (apples-to-apples target,
# per the benchmark design decision to feed externals the package's exact index).
package_composite_index <- function(trait_values, cfg) {
  values <- as.matrix(trait_values)
  values <- values[, cfg$trait_spec$trait, drop = FALSE]
  df <- as.data.frame(values, stringsAsFactors = FALSE)
  names(df) <- cfg$trait_spec$trait
  spec <- cfg$trait_spec
  spec$column <- spec$trait
  objective <- ng_breeder_selection_objective(
    trait = spec,
    method = cfg$multi_trait_method,
    threshold_policy = cfg$threshold_policy
  )
  scored <- ng_score_breeder_objective(
    df, objective,
    out_col = ".composite_index",
    threshold_penalty_weight = cfg$threshold_penalty_weight,
    threshold_penalty_autoscale = cfg$threshold_penalty_autoscale
  )
  as.numeric(scored[[".composite_index"]])
}

# The per-line phenotype/target that single-trait engines (PopVar/SimpleMating) and
# the package's index_as_trait path are trained on. In trait_by_trait mode this is
# the package's composite index; in index_as_trait mode it is the oriented linear
# selection index (which is itself the breeding goal in that mode).
external_training_index <- function(trait_values, cfg) {
  if (identical(cfg$prediction_mode, "trait_by_trait")) {
    package_composite_index(trait_values, cfg)
  } else {
    selection_index_from_trait_values(trait_values, cfg)
  }
}

# Genetic-diversity summary of a breeding pool: what breeders watch erode over
# cycles. `add_var_index` is the population additive genetic variance of the index
# (var(gv); additive traits) INCLUDING LD (Bulmer), not the genic variance; it and
# He/polymorphism/MAF decline as alleles fix, while mean parental RELATIONSHIP
# (VanRaden G off-diagonal ~ 2 x coancestry coefficient) rises.
population_diversity <- function(geno, index_gv, parent_K = NULL) {
  geno <- as.matrix(geno)
  p <- colMeans(geno, na.rm = TRUE) / 2
  p <- p[is.finite(p)]
  He <- if (length(p)) mean(2 * p * (1 - p)) else NA_real_
  poly <- if (length(p)) mean(p > 1e-9 & p < 1 - 1e-9) else NA_real_
  maf <- if (length(p)) mean(pmin(p, 1 - p)) else NA_real_
  if (is.null(parent_K)) parent_K <- ng_parent_kinship(geno)
  offdiag <- parent_K[upper.tri(parent_K)]
  list(
    # Population additive genetic variance of the index INCLUDING linkage
    # disequilibrium (Bulmer effect) -- traits are purely additive so this is Va of
    # the pool, but it is NOT the LD-free genic variance sum(2p(1-p)a^2).
    add_var_index = if (length(index_gv) > 1L) stats::var(index_gv) else 0,
    expected_heterozygosity = He,
    prop_polymorphic_markers = poly,
    mean_maf = maf,
    # Off-diagonals of the VanRaden genomic relationship matrix (G ~ numerator
    # relationship A ~ 2 x kinship coefficient). Reported as "relationship", NOT the
    # coancestry coefficient in [0,1]; it rises as the pool loses diversity.
    mean_parent_relationship = if (length(offdiag)) mean(offdiag, na.rm = TRUE) else NA_real_,
    max_parent_relationship = if (length(offdiag)) max(offdiag, na.rm = TRUE) else NA_real_
  )
}

cor_safe <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  if (stats::sd(x[ok]) <= 0 || stats::sd(y[ok]) <= 0) return(NA_real_)
  suppressWarnings(tryCatch(stats::cor(x[ok], y[ok], method = method), error = function(e) NA_real_))
}

# Realize a set of crosses to RIL families (single-seed descent) to obtain the
# ground-truth realized family performance used to score prediction accuracy.
# No selection is applied: this is a pure predicted-vs-realized calibration set.
realize_cross_families <- function(pop, cross_plan, sim_param, cfg, seed, family_size = NULL) {
  if (is.null(family_size)) family_size <- cfg$accuracy_family_size
  n_cross <- nrow(cross_plan)
  f1 <- with_rng_seed(seed, AlphaSimR::makeCross(pop, crossPlan = cross_plan, nProgeny = 1, simParam = sim_param))
  lines <- with_rng_seed(seed + 1L, AlphaSimR::self(f1, nProgeny = family_size, simParam = sim_param))
  self_gens <- max(0L, as.integer(cfg$accuracy_self_generations) - 1L)
  for (g in seq_len(self_gens)) {
    lines <- with_rng_seed(seed + 10L + g, AlphaSimR::self(lines, nProgeny = 1, simParam = sim_param))
  }
  family <- rep(seq_len(n_cross), each = family_size)
  idx <- selection_index_from_trait_values(gv_matrix(lines, cfg), cfg)
  rows <- lapply(seq_len(n_cross), function(i) {
    vals <- idx[family == i]
    data.frame(
      parent1 = cross_plan[i, 1],
      parent2 = cross_plan[i, 2],
      realized_family_mean = mean(vals),
      realized_family_var = if (length(vals) > 1L) stats::var(vals) else 0,
      realized_family_top10 = top_prop_mean(vals, 0.10),
      realized_family_max = max(vals),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

trait_phenotype_matrix <- function(pop, h2, seed, cfg) {
  g <- gv_matrix(pop, cfg)
  h2 <- max(min(as.numeric(h2), 0.999), 0.001)
  out <- g
  with_rng_seed(seed, {
    for (j in seq_len(ncol(g))) {
      vg <- stats::var(g[, j])
      if (!is.finite(vg) || vg <= 0) vg <- cfg$trait_spec$simulated_variance[[j]]
      ve <- vg * (1 - h2) / h2
      out[, j] <- g[, j] + stats::rnorm(nrow(g), sd = sqrt(ve))
    }
  })
  colnames(out) <- cfg$trait_spec$trait
  out
}

phenotype_values <- function(pop, h2, seed, cfg) {
  selection_index_from_trait_values(trait_phenotype_matrix(pop, h2, seed, cfg), cfg)
}

select_top_individuals <- function(pop, n_keep, h2, seed, cfg) {
  n_keep <- min(as.integer(n_keep), AlphaSimR::nInd(pop))
  if (!is.finite(n_keep) || n_keep < 1L) stop("n_keep must be positive", call. = FALSE)
  pheno <- phenotype_values(pop, h2 = h2, seed = seed, cfg = cfg)
  keep <- order(pheno, decreasing = TRUE)[seq_len(n_keep)]
  list(pop = pop[keep], index = keep, phenotype = pheno[keep])
}

select_top_within_groups <- function(pop, groups, n_per_group, h2, seed, cfg) {
  groups <- as.character(groups)
  if (length(groups) != AlphaSimR::nInd(pop)) {
    stop("group vector length does not match population size", call. = FALSE)
  }
  n_per_group <- max(1L, as.integer(n_per_group))
  pheno <- phenotype_values(pop, h2 = h2, seed = seed, cfg = cfg)
  keep <- unlist(lapply(split(seq_along(groups), groups), function(idx) {
    idx[order(pheno[idx], decreasing = TRUE)[seq_len(min(n_per_group, length(idx)))]]
  }), use.names = FALSE)
  keep <- sort(unique(keep))
  list(pop = pop[keep], groups = groups[keep], index = keep, phenotype = pheno[keep])
}

trait_mean_columns <- function(values, prefix, cfg) {
  values <- as.matrix(values)
  cols <- lapply(seq_len(ncol(values)), function(i) mean(values[, i], na.rm = TRUE))
  names(cols) <- paste0(prefix, "_", make.names(colnames(values)))
  as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
}

stage_summary <- function(pop, stage, rep, method, cycle, h2, cfg, selected = NA) {
  trait_gv <- gv_matrix(pop, cfg)
  g <- selection_index_from_trait_values(trait_gv, cfg)
  cbind(data.frame(
    rep = rep,
    method = method,
    cycle = cycle,
    stage = stage,
    h2 = h2,
    selected = selected,
    prediction_mode = cfg$prediction_mode,
    multi_trait_method = cfg$multi_trait_method,
    n_ind = length(g),
    mean_index_gv = mean(g),
    max_index_gv = max(g),
    top10_index_gv = top_prop_mean(g, 0.10),
    var_index_gv = stats::var(g),
    stringsAsFactors = FALSE
  ), trait_mean_columns(trait_gv, "mean_gv", cfg))
}

family_summary <- function(pop, family_id, rep, method, cycle, selected_crosses, cfg) {
  family_id <- as.integer(family_id)
  trait_gv <- gv_matrix(pop, cfg)
  g <- selection_index_from_trait_values(trait_gv, cfg)
  rows <- lapply(sort(unique(family_id)), function(fid) {
    idx <- which(family_id == fid)
    cross <- selected_crosses[fid, , drop = FALSE]
    cbind(data.frame(
      rep = rep,
      method = method,
      cycle = cycle,
      cross_id = fid,
      parent1 = as.character(cross$parent1[[1]]),
      parent2 = as.character(cross$parent2[[1]]),
      n_lines = length(idx),
      realized_index_mean = mean(g[idx]),
      realized_index_var = if (length(idx) > 1L) stats::var(g[idx]) else 0,
      realized_index_top10 = top_prop_mean(g[idx], 0.10),
      realized_index_max = max(g[idx]),
      stringsAsFactors = FALSE
    ), trait_mean_columns(trait_gv[idx, , drop = FALSE], "realized_mean_gv", cfg))
  })
  bind_rows_fill(rows)
}

advance_ril_pipeline <- function(parent_pop,
                                 selected_crosses,
                                 sim_param,
                                 cfg,
                                 rep,
                                 method,
                                 cycle) {
  h2 <- cfg$stage_heritabilities
  cross_plan <- as.matrix(selected_crosses[, c("parent1", "parent2"), drop = FALSE])
  f1 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F1", seed = cfg$seed),
    AlphaSimR::makeCross(parent_pop, crossPlan = cross_plan, nProgeny = 1, simParam = sim_param)
  )

  f2 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F2", seed = cfg$seed),
    AlphaSimR::self(f1, nProgeny = cfg$f2_per_cross, simParam = sim_param)
  )
  f2_family <- rep(seq_len(nrow(selected_crosses)), each = cfg$f2_per_cross)
  f2_keep_n <- max(1L, round(cfg$f2_per_cross * cfg$within_family_selection_prop))
  f2_sel <- select_top_within_groups(
    f2, f2_family, n_per_group = f2_keep_n, h2 = h2[[1]],
    seed = stable_seed(rep, method, cycle, "select_F2", seed = cfg$seed),
    cfg = cfg
  )

  f3 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F3", seed = cfg$seed),
    AlphaSimR::self(f2_sel$pop, nProgeny = cfg$bulk_size, simParam = sim_param)
  )
  f3_family <- rep(f2_sel$groups, each = cfg$bulk_size)
  f3_sel <- select_top_within_groups(
    f3, f3_family, n_per_group = 1L, h2 = h2[[2]],
    seed = stable_seed(rep, method, cycle, "select_F3", seed = cfg$seed),
    cfg = cfg
  )

  f4 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F4", seed = cfg$seed),
    AlphaSimR::self(f3_sel$pop, nProgeny = cfg$bulk_size, simParam = sim_param)
  )
  f4_family <- rep(f3_sel$groups, each = cfg$bulk_size)
  f4_keep_n <- min(
    AlphaSimR::nInd(f4),
    max(cfg$n_parents, round(AlphaSimR::nInd(f4) * cfg$between_family_selection_prop))
  )
  f4_sel <- select_top_individuals(
    f4, n_keep = f4_keep_n, h2 = h2[[3]],
    seed = stable_seed(rep, method, cycle, "select_F4", seed = cfg$seed),
    cfg = cfg
  )
  f4_sel_family <- f4_family[f4_sel$index]

  f5 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F5", seed = cfg$seed),
    AlphaSimR::self(f4_sel$pop, nProgeny = 1, simParam = sim_param)
  )
  f5_keep_n <- min(
    AlphaSimR::nInd(f5),
    max(cfg$n_parents, round(AlphaSimR::nInd(f5) * cfg$ayt_selection_prop))
  )
  f5_sel <- select_top_individuals(
    f5, n_keep = f5_keep_n, h2 = h2[[4]],
    seed = stable_seed(rep, method, cycle, "select_F5", seed = cfg$seed),
    cfg = cfg
  )
  f5_sel_family <- f4_sel_family[f5_sel$index]

  f6 <- with_rng_seed(
    stable_seed(rep, method, cycle, "F6", seed = cfg$seed),
    AlphaSimR::self(f5_sel$pop, nProgeny = 1, simParam = sim_param)
  )
  f6_family <- f5_sel_family
  parent_sel <- select_top_individuals(
    f6, n_keep = cfg$n_parents, h2 = h2[[5]],
    seed = stable_seed(rep, method, cycle, "select_F6_next_parents", seed = cfg$seed),
    cfg = cfg
  )
  next_parents <- parent_sel$pop
  # RILs are NOT fully homozygous -- they carry ~(1/2)^(g-1) residual heterozygosity per
  # locus, especially when parents are recycled early (here at F6). We therefore keep the
  # RIL parents as-is (do NOT force them to DH homozygosity, which would change the
  # breeding system). The optional flag below finishes parents to doubled-haploid ONLY for
  # users explicitly modelling a DH program; it is off by default so the default study
  # stays a true RIL program.
  if (isTRUE(cfg$make_dh_parents)) {
    next_parents <- with_rng_seed(
      stable_seed(rep, method, cycle, "makeDH_parents", seed = cfg$seed),
      AlphaSimR::makeDH(next_parents, nDH = 1, simParam = sim_param)
    )
  }
  next_parents@id <- paste0(method, "_R", rep, "_C", cycle, "_P", seq_len(AlphaSimR::nInd(next_parents)))

  stages <- bind_rows_fill(list(
    stage_summary(f2, "F2", rep, method, cycle, h2[[1]], cfg = cfg, selected = FALSE),
    stage_summary(f2_sel$pop, "F2_selected", rep, method, cycle, h2[[1]], cfg = cfg, selected = TRUE),
    stage_summary(f3, "F3_bulk", rep, method, cycle, h2[[2]], cfg = cfg, selected = FALSE),
    stage_summary(f3_sel$pop, "F3_selected", rep, method, cycle, h2[[2]], cfg = cfg, selected = TRUE),
    stage_summary(f4, "F4_PYT_candidates", rep, method, cycle, h2[[3]], cfg = cfg, selected = FALSE),
    stage_summary(f4_sel$pop, "F4_PYT_selected", rep, method, cycle, h2[[3]], cfg = cfg, selected = TRUE),
    stage_summary(f5, "F5_AYT_candidates", rep, method, cycle, h2[[4]], cfg = cfg, selected = FALSE),
    stage_summary(f5_sel$pop, "F5_AYT_selected", rep, method, cycle, h2[[4]], cfg = cfg, selected = TRUE),
    stage_summary(f6, "F6_EYT_candidates", rep, method, cycle, h2[[5]], cfg = cfg, selected = FALSE),
    stage_summary(next_parents, "F6_next_cycle_parents", rep, method, cycle, h2[[5]], cfg = cfg, selected = TRUE)
  ))

  list(
    next_parents = next_parents,
    stage_summary = stages,
    family_summary = family_summary(f6, f6_family, rep, method, cycle, selected_crosses, cfg)
  )
}

make_auxiliary_ril_training_lines <- function(parent_pop, target_n, sim_param, cfg, seed) {
  if (target_n <= AlphaSimR::nInd(parent_pop)) return(parent_pop)
  extra_n <- target_n - AlphaSimR::nInd(parent_pop)
  n_par <- AlphaSimR::nInd(parent_pop)
  p1 <- with_rng_seed(seed, sample.int(n_par, extra_n, replace = TRUE))
  p2 <- with_rng_seed(seed + 1L, sample.int(n_par, extra_n, replace = TRUE))
  same <- p1 == p2
  # Re-draw self-pairs (p1 == p2) so training lines come from true crosses. The seed
  # MUST advance every iteration: a seed keyed to sum(same) reproduces the identical
  # resample whenever the collision count is unchanged, which is an infinite loop
  # (needs >= 2 parents to be resolvable at all). Use an iteration counter and a
  # safety cap; fall back to a deterministic offset partner if any remain.
  if (n_par >= 2L) {
    iter <- 0L
    while (any(same) && iter < 1000L) {
      iter <- iter + 1L
      p2[same] <- with_rng_seed(seed + 2L + iter, sample.int(n_par, sum(same), replace = TRUE))
      same <- p1 == p2
    }
  }
  if (any(same)) p2[same] <- ((p1[same]) %% n_par) + 1L  # guaranteed distinct partner
  f1 <- with_rng_seed(seed + 10L, AlphaSimR::makeCross(parent_pop, cbind(p1, p2), nProgeny = 1, simParam = sim_param))
  ril <- f1
  for (g in seq_len(cfg$training_self_generations)) {
    ril <- with_rng_seed(seed + 100L + g, AlphaSimR::self(ril, nProgeny = 1, simParam = sim_param))
  }
  out <- c(parent_pop, ril)
  out[seq_len(target_n)]
}

prepare_marker_assets <- function(sim_param) {
  snp_map <- AlphaSimR::getSnpMap(simParam = sim_param)
  marker_names <- paste0("Chr", snp_map$chr, "_", snp_map$id)
  marker_map <- data.frame(
    marker = marker_names,
    chr = snp_map$chr,
    pos_cm = snp_map$pos * 100,
    stringsAsFactors = FALSE
  )
  list(marker_names = marker_names, marker_map = marker_map)
}

population_genotypes <- function(pop, sim_param, marker_names) {
  geno <- AlphaSimR::pullSnpGeno(pop, simParam = sim_param)
  colnames(geno) <- marker_names
  rownames(geno) <- pop@id
  geno
}

select_top_by_score <- function(scores, score_col, n_crosses) {
  if (!(score_col %in% names(scores))) stop("scores missing score column: ", score_col, call. = FALSE)
  values <- suppressWarnings(as.numeric(scores[[score_col]]))
  ok <- is.finite(values)
  if (sum(ok) < n_crosses) {
    stop("Not enough finite scores in ", score_col, " for requested crosses", call. = FALSE)
  }
  idx <- which(ok)[order(values[ok], decreasing = TRUE)[seq_len(n_crosses)]]
  scores[idx, , drop = FALSE]
}

score_parent_crosses <- function(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle) {
  rule <- method_score_rule(method)
  ids <- parent_pop@id
  parent_geno <- population_genotypes(parent_pop, sim_param, marker_assets$marker_names)
  # Kinship of the parents is trait-independent: compute the VanRaden G once and reuse
  # it for every per-trait ng_score_crosses call and as the returned parent_K, instead
  # of recomputing the O(n^2 * m) matrix 3x (per trait) + 1x below.
  parent_K_cache <- ng_parent_kinship(parent_geno)
  training_pop <- make_auxiliary_ril_training_lines(
    parent_pop = parent_pop,
    target_n = cfg$effect_training_n,
    sim_param = sim_param,
    cfg = cfg,
    seed = stable_seed(rep, method, cycle, "training_lines", seed = cfg$seed)
  )
  training_pop@id <- paste0(method, "_R", rep, "_C", cycle, "_T", seq_len(AlphaSimR::nInd(training_pop)))
  training_geno <- population_genotypes(training_pop, sim_param, marker_assets$marker_names)
  training_trait_y <- trait_phenotype_matrix(
    training_pop,
    h2 = cfg$effect_training_h2,
    seed = stable_seed(rep, method, cycle, "training_pheno", seed = cfg$seed),
    cfg = cfg
  )
  parent_trait_pheno <- trait_phenotype_matrix(
    parent_pop,
    h2 = cfg$parent_observed_h2,
    seed = stable_seed(rep, method, cycle, "parent_pheno", seed = cfg$seed),
    cfg = cfg
  )

  native_multitrait <- identical(cfg$prediction_mode, "trait_by_trait") &&
    identical(rule$external_tool, "nextgenCrossDesign")

  if (isTRUE(native_multitrait)) {
    cross_table <- NULL
    effects_list <- list()
    reliability <- numeric(nrow(cfg$trait_spec))
    for (i in seq_len(nrow(cfg$trait_spec))) {
      trait <- cfg$trait_spec$trait[[i]]
      clean_trait <- ng_run_cp_clean_trait_name(trait)
      y <- as.numeric(training_trait_y[, trait])
      names(y) <- training_pop@id
      effects_i <- ng_fit_ridge_effects(
        geno = training_geno,
        y = y,
        ids = training_pop@id,
        h2_prior = cfg$effect_h2_prior,
        kfold = cfg$effect_kfold,
        seed = stable_seed(rep, method, cycle, "ridge", trait, seed = cfg$seed),
        return_beta_cov_full = identical(rule$method_varPMV, "full_posterior")
      )
      parent_pheno_i <- as.numeric(parent_trait_pheno[, trait])
      names(parent_pheno_i) <- ids
      scored_trait <- tryCatch(
        ng_score_crosses(
          geno = parent_geno,
          parent_K = parent_K_cache,
          effects = effects_i,
          marker_map = marker_assets$marker_map,
          ids = ids,
          adjusted_pheno = parent_pheno_i,
          target = "RIL",
          selection_prop = cfg$selection_prop,
          min_effect_reliability = cfg$min_effect_reliability,
          recomb_model = cfg$recombination_model,
          use_cpp = cfg$use_cpp,
          assume_inbred = cfg$assume_inbred_for_scoring,
          posterior_cov_full = effects_i$beta_cov_full
        ),
        error = function(e) {
          if (!isTRUE(cfg$retry_scoring_as_noninbred)) stop(e)
          warning(conditionMessage(e), call. = FALSE)
          ng_score_crosses(
            geno = parent_geno,
            parent_K = parent_K_cache,
            effects = effects_i,
            marker_map = marker_assets$marker_map,
            ids = ids,
            adjusted_pheno = parent_pheno_i,
            target = "RIL",
            selection_prop = cfg$selection_prop,
            min_effect_reliability = cfg$min_effect_reliability,
            recomb_model = cfg$recombination_model,
            use_cpp = cfg$use_cpp,
            assume_inbred = FALSE,
            posterior_cov_full = effects_i$beta_cov_full
          )
        }
      )
      value <- ng_run_cp_trait_value(
        scored_trait = scored_trait,
        direction = cfg$trait_spec$direction[[i]],
        trait_value_metric = rule$trait_value_metric,
        uc_variance_source = rule$uc_variance_source,
        selection_prop = cfg$selection_prop,
        method_varPMV = rule$method_varPMV
      )
      if (is.null(cross_table)) {
        cross_table <- scored_trait[, c("parent1", "parent2", "pair_kinship"), drop = FALSE]
      }
      cross_table[[paste0(clean_trait, "_value")]] <- value
      cross_table[[paste0(clean_trait, "_mean")]] <- scored_trait$cross_mean_blend
      cross_table[[paste0(clean_trait, "_pmv")]] <- scored_trait$dh_pmv_var
      cross_table[[paste0(clean_trait, "_vpm")]] <- scored_trait$dh_recomb_var
      effects_list[[trait]] <- effects_i
      reliability[[i]] <- effects_i$reliability
    }
    objective_traits <- cfg$trait_spec
    objective_traits$column <- vapply(objective_traits$trait, function(trait) {
      paste0(ng_run_cp_clean_trait_name(trait), "_value")
    }, character(1L))
    objective <- ng_breeder_selection_objective(
      trait = objective_traits,
      method = cfg$multi_trait_method,
      threshold_policy = cfg$threshold_policy
    )
    scores <- ng_score_breeder_objective(
      cross_table,
      objective,
      out_col = ".package_method_score",
      threshold_penalty_weight = cfg$threshold_penalty_weight,
      threshold_penalty_autoscale = cfg$threshold_penalty_autoscale
    )
    return(list(
      scores = scores,
      parent_geno = parent_geno,
      parent_K = parent_K_cache,
      effects = effects_list,
      parent_pheno = selection_index_from_trait_values(parent_trait_pheno, cfg),
      parent_trait_pheno = parent_trait_pheno,   # per-trait phenotypes for the user-API method
      effect_reliability = mean(reliability, na.rm = TRUE),
      scoring_basis = "trait_by_trait_multi_trait_directional"
    ))
  }

  # Single-trait engines (PopVar/SimpleMating) and the package's index_as_trait
  # path are trained on external_training_index(): the package's composite index in
  # trait_by_trait mode (so externals target the EXACT index the package optimizes),
  # or the oriented linear selection index in index_as_trait mode.
  training_y <- external_training_index(training_trait_y, cfg)
  names(training_y) <- training_pop@id
  effects <- ng_fit_ridge_effects(
    geno = training_geno,
    y = training_y,
    ids = training_pop@id,
    h2_prior = cfg$effect_h2_prior,
    kfold = cfg$effect_kfold,
    seed = stable_seed(rep, method, cycle, "ridge", "selection_index", seed = cfg$seed),
    return_beta_cov_full = identical(rule$method_varPMV, "full_posterior")
  )
  parent_pheno <- external_training_index(parent_trait_pheno, cfg)
  names(parent_pheno) <- ids

  scores <- tryCatch(
    ng_score_crosses(
      geno = parent_geno,
      effects = effects,
      marker_map = marker_assets$marker_map,
      ids = ids,
      adjusted_pheno = parent_pheno,
      target = "RIL",
      selection_prop = cfg$selection_prop,
      min_effect_reliability = cfg$min_effect_reliability,
      recomb_model = cfg$recombination_model,
      use_cpp = cfg$use_cpp,
      assume_inbred = cfg$assume_inbred_for_scoring,
      posterior_cov_full = effects$beta_cov_full
    ),
    error = function(e) {
      if (!isTRUE(cfg$retry_scoring_as_noninbred)) stop(e)
      warning(conditionMessage(e), call. = FALSE)
      ng_score_crosses(
        geno = parent_geno,
        effects = effects,
        marker_map = marker_assets$marker_map,
        ids = ids,
        adjusted_pheno = parent_pheno,
        target = "RIL",
        selection_prop = cfg$selection_prop,
        min_effect_reliability = cfg$min_effect_reliability,
        recomb_model = cfg$recombination_model,
        use_cpp = cfg$use_cpp,
        assume_inbred = FALSE,
        posterior_cov_full = effects$beta_cov_full
      )
    }
  )

  list(
    scores = scores,
    parent_geno = parent_geno,
    parent_K = ng_parent_kinship(parent_geno),
    effects = effects,
    parent_pheno = parent_pheno,
    effect_reliability = effects$reliability,
    scoring_basis = if (identical(cfg$prediction_mode, "trait_by_trait")) {
      "package_composite_index_single_trait"
    } else {
      "index_as_trait_increase"
    }
  )
}

ril_benchmark_method_registry <- function() {
  base <- data.frame(
    method = c(
      "random",
      "var_complex_ocs",
      "var_complex_alphamate_style",
      "var_complex_alphamate_executable",
      "uc_pmv_ocs",
      "uc_vpm_ocs",
      "pmv_ocs",
      "vpm_ocs",
      "var_simple_ocs",
      "mean_ocs",
      "popvar_uc_topn",
      "simple_usefa_topn",
      "simple_usefa_select"
    ),
    source = c(
      "breeder_control",
      "inst/examples/05_ocs_user_run.R; inst/examples/03_variance_method_comparison.R",
      "inst/examples/18_allocation_method_comparison.R; inst/examples/03_variance_method_comparison.R",
      "inst/examples/04_alphamate_executable_user_run.R; inst/examples/18_allocation_method_comparison.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/07_trait_value_metric_parameter_guide.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/07_trait_value_metric_parameter_guide.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/20_full_posterior_pmv_shortlist.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/07_trait_value_metric_parameter_guide.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/07_trait_value_metric_parameter_guide.R",
      "inst/examples/03_variance_method_comparison.R; inst/examples/07_trait_value_metric_parameter_guide.R",
      "R/08_external_baselines.R; R/25_backend_capability_registry.R",
      "R/08_external_baselines.R; R/25_backend_capability_registry.R",
      "R/08_external_baselines.R; R/25_backend_capability_registry.R"
    ),
    trait_value_metric = c(
      NA,
      "var_complex",
      "var_complex",
      "var_complex",
      "uc",
      "uc",
      "pmv",
      "vpm",
      "var_simple",
      "mean",
      NA,
      NA,
      NA
    ),
    uc_variance_source = c(
      NA,
      "pmv",
      "pmv",
      "pmv",
      "pmv",
      "vpm",
      "pmv",
      "pmv",
      "pmv",
      "pmv",
      NA,
      NA,
      NA
    ),
    method_varPMV = c(
      NA,
      "fast",
      "fast",
      "fast",
      "fast",
      "fast",
      "full_posterior",
      "fast",
      "fast",
      "fast",
      NA,
      NA,
      NA
    ),
    allocation_method = c(
      "random",
      "ocs",
      "alphamate_style",
      "alphamate_executable",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "external_score_topn",
      "external_score_topn",
      "simplemating_select"
    ),
    allocator = c(
      "random",
      "ocs",
      "alphamate_style",
      "alphamate_executable",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "ocs",
      "topn",
      "topn",
      "simplemating_select"
    ),
    score_col = c(
      NA,
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      ".package_method_score",
      "popvar_uc",
      "simple_usefa",
      "simple_usefa"
    ),
    implementation = c(
      "native_random",
      "native_nextgen",
      "native_nextgen",
      "exact_external",
      "native_nextgen",
      "native_nextgen",
      "native_nextgen",
      "native_nextgen",
      "native_nextgen",
      "native_nextgen",
      "exact_external",
      "exact_external",
      "exact_external"
    ),
    external_tool = c(
      "none",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "AlphaMate",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "nextgenCrossDesign",
      "PopVar",
      "SimpleMating",
      "SimpleMating"
    ),
    default = c(
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE
    ),
    stringsAsFactors = FALSE
  )
  # Complete the PopVar/SimpleMating comparison across ALL criteria the package
  # supports (VALIDATED_STATE.md documents popvar_mu/var/uc/musp and
  # simple_mpv/usefa as validated external branches). Restricting the benchmark to
  # a single external criterion would undersell the breadth of the comparison. All
  # remain exact_external calls (engine = "external") with provenance columns.
  extra <- data.frame(
    method = c("popvar_mu_topn", "popvar_var_topn", "popvar_musp_topn",
               "simple_mpv_topn", "simple_mpv_select"),
    source = rep("R/08_external_baselines.R; R/25_backend_capability_registry.R", 5),
    trait_value_metric = rep(NA_character_, 5),
    uc_variance_source = rep(NA_character_, 5),
    method_varPMV = rep(NA_character_, 5),
    allocation_method = c("external_score_topn", "external_score_topn", "external_score_topn",
                          "external_score_topn", "simplemating_select"),
    allocator = c("topn", "topn", "topn", "topn", "simplemating_select"),
    score_col = c("popvar_mu", "popvar_varG", "popvar_musp_high", "simple_mpv", "simple_mpv"),
    implementation = rep("exact_external", 5),
    external_tool = c("PopVar", "PopVar", "PopVar", "SimpleMating", "SimpleMating"),
    default = rep(TRUE, 5),
    stringsAsFactors = FALSE
  )
  # Native evolutionary (memetic GA) allocator on the flagship var_complex score, so
  # the study can compare the new optimizer head-to-head against the OCS/MIP path.
  evo <- data.frame(
    method = "var_complex_evolution",
    source = "R/40_evolutionary_optimizer.R; R/04_optimizers.R",
    trait_value_metric = "var_complex", uc_variance_source = "pmv", method_varPMV = "fast",
    allocation_method = "evolution", allocator = "evolution", score_col = ".package_method_score",
    implementation = "native_nextgen", external_tool = "nextgenCrossDesign", default = TRUE,
    stringsAsFactors = FALSE
  )
  # Gain-diversity balancing dial (R/41): the same var_complex score allocated at three
  # points on the frontier so the study can see the realized gain-vs-diversity trade-off.
  strat <- data.frame(
    method = c("strategy_high_gain", "strategy_balanced", "strategy_diversity"),
    source = rep("R/41_strategy_frontier.R; R/04_optimizers.R", 3),
    trait_value_metric = rep("var_complex", 3), uc_variance_source = rep("pmv", 3),
    method_varPMV = rep("fast", 3),
    allocation_method = rep("strategy", 3),
    allocator = c("strategy_high_gain", "strategy_balanced", "strategy_diversity"),
    score_col = rep(".package_method_score", 3),
    implementation = rep("native_nextgen", 3), external_tool = rep("nextgenCrossDesign", 3),
    default = rep(FALSE, 3),
    stringsAsFactors = FALSE
  )
  # The package driven through its TOP-LEVEL user API (ng_run_cross_prediction), exactly
  # as an end user would call it (data frames in, selected crosses out) -- the same
  # black-box way SimpleMating is invoked -- so the study compares user-API to user-API,
  # not the internal ng_score_crosses/ng_optimize_mating_plan calls the other package
  # methods use.
  usr <- data.frame(
    method = "package_user_api",
    source = "R/39_cross_prediction_runner.R (ng_run_cross_prediction)",
    trait_value_metric = "var_complex", uc_variance_source = "pmv", method_varPMV = "fast",
    allocation_method = "user_api", allocator = "user_api", score_col = ".package_method_score",
    implementation = "native_user_api", external_tool = "nextgenCrossDesign", default = FALSE,
    stringsAsFactors = FALSE
  )
  rbind(base, extra, evo, strat, usr)
}

ril_benchmark_default_methods <- function() {
  registry <- ril_benchmark_method_registry()
  registry$method[registry$default]
}

method_availability <- function(methods, cfg) {
  registry <- ril_benchmark_method_registry()
  rows <- lapply(methods, function(method) {
    rule <- registry[match(method, registry$method), , drop = FALSE]
    if (!nrow(rule) || is.na(rule$method[[1]])) {
      return(data.frame(
        method = method,
        active = FALSE,
        external_tool = "unknown",
        exact_external_status = "unsupported_method",
        reason = "Method is not in ril_benchmark_method_registry()",
        stringsAsFactors = FALSE
      ))
    }
    active <- TRUE
    status <- "not_applicable"
    reason <- ""
    if (identical(rule$external_tool[[1]], "PopVar")) {
      active <- requireNamespace("PopVar", quietly = TRUE)
      status <- if (active) "available" else "package_unavailable"
      reason <- if (active) "" else "PopVar is not installed"
    } else if (identical(rule$external_tool[[1]], "SimpleMating")) {
      active <- requireNamespace("SimpleMating", quietly = TRUE)
      status <- if (active) "available" else "package_unavailable"
      reason <- if (active) "" else "SimpleMating is not installed"
    } else if (identical(rule$external_tool[[1]], "AlphaMate")) {
      active <- nzchar(cfg$alphamate_executable) && file.exists(cfg$alphamate_executable)
      status <- if (active) "available" else "executable_unavailable"
      reason <- if (active) "" else "AlphaMate executable is not configured"
    } else {
      status <- "not_applicable"
    }
    data.frame(
      method = method,
      active = active,
      external_tool = rule$external_tool[[1]],
      implementation = rule$implementation[[1]],
      source = rule$source[[1]],
      trait_value_metric = rule$trait_value_metric[[1]],
      uc_variance_source = rule$uc_variance_source[[1]],
      method_varPMV = rule$method_varPMV[[1]],
      allocation_method = rule$allocation_method[[1]],
      exact_external_status = status,
      reason = reason,
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

method_score_rule <- function(method) {
  registry <- ril_benchmark_method_registry()
  idx <- match(method, registry$method)
  if (!is.finite(idx) || is.na(idx)) stop("Unsupported benchmark method: ", method, call. = FALSE)
  as.list(registry[idx, , drop = FALSE])
}

add_package_method_score <- function(scores, rule, cfg) {
  if (".package_method_score" %in% names(scores)) {
    return(scores)
  }
  if (!is.character(rule$trait_value_metric) || is.na(rule$trait_value_metric) ||
      !nzchar(rule$trait_value_metric)) {
    return(scores)
  }
  scores$.package_method_score <- ng_run_cp_trait_value(
    scored_trait = scores,
    direction = "increase",
    trait_value_metric = rule$trait_value_metric,
    uc_variance_source = rule$uc_variance_source,
    selection_prop = cfg$selection_prop,
    method_varPMV = rule$method_varPMV
  )
  scores
}

# Score all candidate crosses for a method WITHOUT allocation. Shared by
# select_crosses_for_method() (which then allocates) and evaluate_shared_accuracy()
# (which correlates the predicted cross values against realized family performance).
# For PopVar/SimpleMating the exact external engines populate every criterion column
# (popvar_mu/varG/musp/uc, simple_mpv/usefa); the method's score_col selects which.
compute_method_scores <- function(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle) {
  scoring <- score_parent_crosses(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle)
  scores <- scoring$scores
  rule <- method_score_rule(method)
  scores <- add_package_method_score(scores, rule, cfg)
  exact_status <- "not_applicable"
  fallback_reason <- ""

  if (grepl("^popvar_", method)) {
    exact_status <- "ok"
    scores <- tryCatch(
      ng_add_popvar_scores(
        scores = scores,
        geno = scoring$parent_geno,
        effects = scoring$effects,
        marker_map = marker_assets$marker_map,
        adjusted_pheno = scoring$parent_pheno,
        tail_p = cfg$selection_prop,
        engine = "external",
        fallback = "none"
      ),
      error = function(e) {
        exact_status <<- "failed"
        fallback_reason <<- conditionMessage(e)
        scores
      }
    )
    if ("popvar_status" %in% names(scores)) {
      exact_status <- paste(unique(as.character(scores$popvar_status)), collapse = "|")
    }
  }

  if (grepl("^simple_", method)) {
    exact_status <- "ok"
    scores <- tryCatch(
      ng_add_simplemating_scores(
        scores = scores,
        geno = scoring$parent_geno,
        effects = scoring$effects,
        marker_map = marker_assets$marker_map,
        adjusted_pheno = scoring$parent_pheno,
        prop_sel = cfg$selection_prop,
        # Match SimpleMating's progeny-variance formula to the breeding system the study
        # models: RIL by default (the package's own scoring uses target = "RIL"), or DH
        # when parents are finished as doubled haploids. getUsefA's Type must agree with
        # target so the two engines score the SAME kind of cross.
        type = if (isTRUE(cfg$make_dh_parents)) "DH" else "RIL",
        generation = cfg$simplemating_generation,
        n_threads = cfg$simplemating_threads,
        engine = "external",
        fallback = "none"
      ),
      error = function(e) {
        exact_status <<- "failed"
        fallback_reason <<- conditionMessage(e)
        scores
      }
    )
    if ("simple_status" %in% names(scores)) {
      exact_status <- paste(unique(as.character(scores$simple_status)), collapse = "|")
    }
  }

  list(scores = scores, rule = rule, scoring = scoring,
       exact_status = exact_status, fallback_reason = fallback_reason)
}

select_crosses_for_method <- function(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle) {
  computed <- compute_method_scores(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle)
  scores <- computed$scores
  rule <- computed$rule
  scoring <- computed$scoring
  exact_status <- computed$exact_status
  fallback_reason <- computed$fallback_reason

  selected <- if (identical(rule$allocator, "random")) {
    idx <- with_rng_seed(
      stable_seed(rep, method, cycle, "random_crosses", seed = cfg$seed),
      sample.int(nrow(scores), min(cfg$n_crosses, nrow(scores)))
    )
    scores[idx, , drop = FALSE]
  } else if (identical(rule$allocator, "topn")) {
    select_top_by_score(scores, rule$score_col, cfg$n_crosses)
  } else if (identical(rule$allocator, "ocs")) {
    ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$n_crosses,
      gain_col = rule$score_col,
      parent_K = scoring$parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = cfg$lambda_parent_use,
      lambda_parent_use_mode = cfg$lambda_parent_use_mode,
      method = cfg$optimizer,
      local_iter = cfg$local_iter,
      ocs_iter = cfg$ocs_iter
    )
  } else if (identical(rule$allocator, "evolution")) {
    ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$n_crosses,
      gain_col = rule$score_col,
      parent_K = scoring$parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_group = cfg$lambda_group,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = cfg$lambda_parent_use,
      lambda_parent_use_mode = cfg$lambda_parent_use_mode,
      method = "evolution",
      local_iter = cfg$local_iter,
      evol_solutions = cfg$evol_solutions,
      evol_iterations = cfg$evol_iterations,
      evol_stop = cfg$evol_stop,
      evol_seed = stable_seed(rep, method, cycle, "evolution", seed = cfg$seed)
    )
  } else if (grepl("^strategy_", rule$allocator)) {
    # Gain-diversity balancing: allocate at the strategy's frontier point. The strategy
    # dispatch sweeps its own lambda_group, so we do not pass one here.
    ng_optimize_mating_plan(
      scores = scores,
      n_crosses = cfg$n_crosses,
      gain_col = rule$score_col,
      parent_K = scoring$parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      max_pair_kinship = cfg$max_pair_kinship,
      lambda_mating = cfg$lambda_mating,
      lambda_parent_use = cfg$lambda_parent_use,
      lambda_parent_use_mode = cfg$lambda_parent_use_mode,
      strategy = sub("^strategy_", "", rule$allocator),
      method = cfg$optimizer,
      local_iter = cfg$local_iter
    )
  } else if (identical(rule$allocator, "user_api")) {
    # Drive the package through its top-level user API exactly as an end user would:
    # hand ng_run_cross_prediction() the parent genotypes + per-trait phenotypes + map +
    # directions as data frames, and take back the selected crosses. This exercises the
    # real user-facing workflow (QC -> effect estimation -> DH/RIL usefulness -> multi-trait
    # index -> OCS), the same black-box way SimpleMating is called, rather than the internal
    # function calls the other package methods use.
    api_ids <- rownames(scoring$parent_geno)
    geno_df <- data.frame(NAME = api_ids, scoring$parent_geno, check.names = FALSE, stringsAsFactors = FALSE)
    pheno_df <- data.frame(NAME = api_ids, as.data.frame(scoring$parent_trait_pheno, stringsAsFactors = FALSE),
                           check.names = FALSE, stringsAsFactors = FALSE)
    mmap <- marker_assets$marker_map
    map_df <- data.frame(SNP = mmap$marker, Chr = mmap$chr, PosCM = mmap$pos_cm, stringsAsFactors = FALSE)
    dir_df <- data.frame(trait = cfg$trait_spec$trait, column = cfg$trait_spec$trait,
                         direction = cfg$trait_spec$direction, stringsAsFactors = FALSE)
    api <- ng_run_cross_prediction(
      phenotype = pheno_df, genotype = geno_df, marker_map = map_df, trait_direction = dir_df,
      id_col = "NAME", map_marker_col = "SNP", map_chr_col = "Chr", map_pos_cm_col = "PosCM",
      map_position_unit = "cM", prediction_mode = "trait_by_trait",
      trait_value_metric = "var_complex", uc_variance_source = "pmv",
      multi_trait_method = cfg$multi_trait_method, trait_weights = cfg$trait_weights,
      progeny = "RIL", ril_mode = "infinite",
      n_crosses = cfg$n_crosses, max_uses_per_parent = cfg$max_crosses_per_parent,
      use_ocs = TRUE, lambda_group = cfg$lambda_group, lambda_mating = cfg$lambda_mating,
      lambda_parent_use = cfg$lambda_parent_use, optimizer = cfg$optimizer,
      # Recurrent RIL parents carry residual heterozygosity from finite selfing; a user
      # scores them with assume_inbred = FALSE (the top-level API has no internal
      # het-retry, unlike the benchmark's native scoring path).
      duplicate_action = "none", assume_inbred = FALSE,
      use_cpp = cfg$use_cpp, write_outputs = FALSE, write_figures = FALSE,
      seed = stable_seed(rep, method, cycle, "user_api", seed = cfg$seed))
    sel <- api$selected_crosses
    # match the API's selected crosses back to the candidate scores rows (order-insensitive)
    ck <- ng_pair_key(scores$parent1, scores$parent2)
    fwd <- match(ng_pair_key(sel$parent1, sel$parent2), ck)
    rev <- match(ng_pair_key(sel$parent2, sel$parent1), ck)
    idx <- ifelse(is.na(fwd), rev, fwd)
    scores[idx[!is.na(idx)], , drop = FALSE]
  } else if (identical(rule$allocator, "alphamate_style")) {
    ng_alphamate_style_select(
      scores = scores,
      criterion_col = rule$score_col,
      n_crosses = cfg$n_crosses,
      parent_K = scoring$parent_K,
      mode = cfg$alphamate_mode,
      target_degree = cfg$alphamate_target_degree,
      max_contributions = cfg$alphamate_max_contributions,
      lambda_group = cfg$alphamate_lambda_group,
      method = cfg$optimizer,
      local_iter = cfg$local_iter
    )
  } else if (identical(rule$allocator, "simplemating_select")) {
    ng_select_simplemating(
      scores = scores,
      score_col = rule$score_col,
      n_crosses = cfg$n_crosses,
      parent_K = scoring$parent_K,
      max_crosses_per_parent = cfg$max_crosses_per_parent,
      min_crosses_per_parent = cfg$simplemating_min_cross,
      culling_pairwise_k = if (is.finite(cfg$simplemating_culling_k)) cfg$simplemating_culling_k else NULL
    )
  } else if (identical(rule$allocator, "alphamate_executable")) {
    exact_status <- "ok"
    ng_select_alphamate(
      scores = scores,
      criterion_col = rule$score_col,
      n_crosses = cfg$n_crosses,
      parent_K = scoring$parent_K,
      executable = cfg$alphamate_executable,
      runtime_path = cfg$alphamate_runtime_path,
      target_degree = cfg$alphamate_target_degree,
      max_contributions = cfg$alphamate_max_contributions,
      number_of_parents = cfg$alphamate_number_of_parents,
      evol_solutions = cfg$alphamate_evol_solutions,
      evol_iterations = cfg$alphamate_evol_iterations,
      evol_stop = cfg$alphamate_evol_stop,
      n_threads = cfg$alphamate_threads,
      keep_files = cfg$alphamate_keep_files
    )
  } else {
    stop("Unsupported allocator: ", rule$allocator, call. = FALSE)
  }

  plan_summary <- attr(selected, "summary")
  selected <- as.data.frame(selected, stringsAsFactors = FALSE)
  selected <- selected[seq_len(min(nrow(selected), cfg$n_crosses)), , drop = FALSE]
  selected$cross_id <- seq_len(nrow(selected))
  selected$method <- method
  selected$cycle <- cycle
  selected$rep <- rep
  selected$score_col <- rule$score_col
  selected$allocator <- rule$allocator
  selected$implementation <- rule$implementation
  selected$source <- rule$source
  selected$trait_value_metric <- rule$trait_value_metric
  selected$uc_variance_source <- rule$uc_variance_source
  selected$method_varPMV <- rule$method_varPMV
  selected$allocation_method <- rule$allocation_method
  selected$prediction_mode <- cfg$prediction_mode
  selected$multi_trait_method <- cfg$multi_trait_method
  selected$trait_direction <- paste(cfg$trait_spec$trait, cfg$trait_spec$direction, sep = ":", collapse = ";")
  selected$index_direction <- "increase"
  selected$scoring_basis <- scoring$scoring_basis
  selected$exact_external_status <- exact_status
  selected$fallback_reason <- fallback_reason

  if (is.null(plan_summary)) plan_summary <- list()
  mean_gain <- if (!is.null(plan_summary$mean_gain)) {
    plan_summary$mean_gain
  } else if (is.character(rule$score_col) && length(rule$score_col) == 1L &&
             !is.na(rule$score_col) && rule$score_col %in% names(selected)) {
    mean(selected[[rule$score_col]], na.rm = TRUE)
  } else {
    NA_real_
  }
  selection_summary <- data.frame(
    rep = rep,
    method = method,
    cycle = cycle,
    selected_crosses = nrow(selected),
    score_col = rule$score_col,
    allocator = rule$allocator,
    implementation = rule$implementation,
    source = rule$source,
    trait_value_metric = rule$trait_value_metric,
    uc_variance_source = rule$uc_variance_source,
    method_varPMV = rule$method_varPMV,
    allocation_method = rule$allocation_method,
    prediction_mode = cfg$prediction_mode,
    multi_trait_method = cfg$multi_trait_method,
    trait_direction = paste(cfg$trait_spec$trait, cfg$trait_spec$direction, sep = ":", collapse = ";"),
    index_direction = "increase",
    scoring_basis = scoring$scoring_basis,
    exact_external_status = exact_status,
    fallback_reason = fallback_reason,
    effect_reliability = scoring$effect_reliability,
    mean_gain = mean_gain,
    unique_parents = length(unique(c(selected$parent1, selected$parent2))),
    max_parent_use = max(tabulate(match(c(selected$parent1, selected$parent2), unique(c(selected$parent1, selected$parent2))))),
    group_coancestry = if (!is.null(plan_summary$group_coancestry)) plan_summary$group_coancestry else NA_real_,
    mean_pair_kinship = mean(selected$pair_kinship, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  list(selected = selected, selection_summary = selection_summary)
}

# Pick the predicted-variance column a method exposes, for the variance-prediction
# accuracy axis (recombination-aware segregation variance is the package's core
# differentiator, so the comparison must include it, not just the mean/usefulness).
method_variance_col <- function(rule, score_names) {
  if (identical(rule$external_tool, "PopVar") && "popvar_varG" %in% score_names) return("popvar_varG")
  if (identical(rule$external_tool, "SimpleMating") && "simple_usefa_var" %in% score_names) return("simple_usefa_var")
  if ("dh_pmv_var" %in% score_names) return("dh_pmv_var")
  if ("dh_recomb_var" %in% score_names) return("dh_recomb_var")
  NA_character_
}

# Unbiased prediction accuracy: on the shared base parent pool (identical for all
# methods at rep start), every method scores the SAME candidate crosses; a common
# random sample of those crosses is realized to RIL families (shared ground truth);
# each method's predicted cross value/variance is correlated with realized family
# performance. Not range-restricted, so accuracies are comparable across methods.
evaluate_shared_accuracy <- function(base_parents, sim_param, marker_assets, cfg, rep, active_methods) {
  if (!isTRUE(cfg$accuracy_enable)) return(data.frame())
  n <- AlphaSimR::nInd(base_parents)
  if (n < 4L) return(data.frame())
  ids <- base_parents@id
  all_pairs <- t(utils::combn(n, 2))
  n_cand <- min(nrow(all_pairs), cfg$accuracy_n_candidates)
  cand_rows <- with_rng_seed(stable_seed(rep, "acc_cand", seed = cfg$seed),
                             sort(sample.int(nrow(all_pairs), n_cand)))
  cand <- all_pairs[cand_rows, , drop = FALSE]
  n_real <- min(nrow(cand), cfg$accuracy_n_realized)
  real_rows <- with_rng_seed(stable_seed(rep, "acc_real", seed = cfg$seed),
                             sort(sample.int(nrow(cand), n_real)))
  real_pairs <- cand[real_rows, , drop = FALSE]
  cross_plan <- cbind(ids[real_pairs[, 1]], ids[real_pairs[, 2]])
  realized <- realize_cross_families(base_parents, cross_plan, sim_param, cfg,
                                     seed = stable_seed(rep, "acc_realize", seed = cfg$seed))
  key <- function(a, b) paste(a, b, sep = "__")
  realized$.key <- key(realized$parent1, realized$parent2)

  aligned <- function(s, col) {
    fwd <- stats::setNames(as.numeric(s[[col]]), key(as.character(s$parent1), as.character(s$parent2)))
    rev <- stats::setNames(as.numeric(s[[col]]), key(as.character(s$parent2), as.character(s$parent1)))
    out <- fwd[realized$.key]
    miss <- is.na(out)
    if (any(miss)) out[miss] <- rev[realized$.key[miss]]
    as.numeric(out)
  }

  rows <- lapply(active_methods, function(method) {
    rule <- method_score_rule(method)
    sc <- rule$score_col
    if (!is.character(sc) || is.na(sc) || !nzchar(sc)) return(NULL)  # random: no prediction
    computed <- tryCatch(
      compute_method_scores(base_parents, sim_param, marker_assets, cfg, rep, method, cycle = 0L),
      error = function(e) NULL
    )
    if (is.null(computed) || !(sc %in% names(computed$scores))) return(NULL)
    s <- computed$scores
    predicted <- aligned(s, sc)
    var_col <- method_variance_col(rule, names(s))
    pred_var <- if (!is.na(var_col)) aligned(s, var_col) else rep(NA_real_, nrow(realized))
    data.frame(
      rep = rep, method = method,
      external_tool = rule$external_tool, implementation = rule$implementation,
      score_col = sc, variance_col = ifelse(is.na(var_col), "", var_col),
      exact_external_status = computed$exact_status,
      n_realized = sum(is.finite(predicted) & is.finite(realized$realized_family_mean)),
      acc_pearson_mean = cor_safe(predicted, realized$realized_family_mean, "pearson"),
      acc_spearman_mean = cor_safe(predicted, realized$realized_family_mean, "spearman"),
      acc_pearson_top10 = cor_safe(predicted, realized$realized_family_top10, "pearson"),
      acc_spearman_top10 = cor_safe(predicted, realized$realized_family_top10, "spearman"),
      acc_pearson_max = cor_safe(predicted, realized$realized_family_max, "pearson"),
      acc_var_pearson = cor_safe(pred_var, realized$realized_family_var, "pearson"),
      acc_var_spearman = cor_safe(pred_var, realized$realized_family_var, "spearman"),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

# Range-restricted, in-program accuracy for context: within each rep/method/cycle,
# correlate the predicted selection score of the crosses a method actually made
# against their realized family performance. Reuses existing benchmark outputs.
selected_cross_accuracy <- function(selected_df, family_df) {
  if (!nrow(selected_df) || !nrow(family_df)) return(data.frame())
  groups <- unique(selected_df[, c("rep", "method", "cycle"), drop = FALSE])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    gk <- groups[i, ]
    sc <- method_score_rule(gk$method)$score_col
    if (!is.character(sc) || is.na(sc) || !nzchar(sc)) return(NULL)
    s <- selected_df[selected_df$rep == gk$rep & selected_df$method == gk$method &
                       selected_df$cycle == gk$cycle, , drop = FALSE]
    f <- family_df[family_df$rep == gk$rep & family_df$method == gk$method &
                     family_df$cycle == gk$cycle, , drop = FALSE]
    if (!(sc %in% names(s)) || !nrow(f)) return(NULL)
    m <- merge(s[, c("cross_id", sc)], f[, c("cross_id", "realized_index_mean", "realized_index_top10")],
               by = "cross_id")
    if (nrow(m) < 3L) return(NULL)
    data.frame(
      rep = gk$rep, method = gk$method, cycle = gk$cycle, score_col = sc, n_selected = nrow(m),
      sel_pearson_mean = cor_safe(m[[sc]], m$realized_index_mean, "pearson"),
      sel_spearman_mean = cor_safe(m[[sc]], m$realized_index_mean, "spearman"),
      sel_pearson_top10 = cor_safe(m[[sc]], m$realized_index_top10, "pearson"),
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

make_base_parents <- function(founder, sim_param, cfg, rep) {
  base_pop <- AlphaSimR::newPop(founder, simParam = sim_param)
  n_per_founder <- ceiling(cfg$n_parents / AlphaSimR::nInd(base_pop))
  dh <- with_rng_seed(
    stable_seed(rep, "base_dh", seed = cfg$seed),
    AlphaSimR::makeDH(base_pop, nDH = n_per_founder, keepParents = FALSE, simParam = sim_param)
  )
  parents <- dh[seq_len(cfg$n_parents)]
  parents@id <- paste0("base_R", rep, "_P", seq_len(AlphaSimR::nInd(parents)))
  parents
}

software_versions <- function(cfg) {
  pkg_version <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(pkg))
  }
  data.frame(
    software = c("R", "nextgenCrossDesign", "AlphaSimR", "PopVar", "SimpleMating", "AlphaMate"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      tryCatch(as.character(utils::packageVersion("nextgenCrossDesign")), error = function(e) "source_tree"),
      pkg_version("AlphaSimR"),
      pkg_version("PopVar"),
      pkg_version("SimpleMating"),
      NA_character_
    ),
    status = c(
      "available",
      "available",
      if (requireNamespace("AlphaSimR", quietly = TRUE)) "available" else "package_unavailable",
      if (requireNamespace("PopVar", quietly = TRUE)) "available" else "package_unavailable",
      if (requireNamespace("SimpleMating", quietly = TRUE)) "available" else "package_unavailable",
      if (nzchar(cfg$alphamate_executable) && file.exists(cfg$alphamate_executable)) "available" else "executable_unavailable"
    ),
    path = c(
      R.home(),
      cfg$root,
      NA_character_,
      NA_character_,
      NA_character_,
      cfg$alphamate_executable
    ),
    stringsAsFactors = FALSE
  )
}

add_dependency_provenance <- function(software_df, dependency_df) {
  software_df$source <- NA_character_
  software_df$repository <- NA_character_
  software_df$ref <- NA_character_
  software_df$dependency_status <- NA_character_
  if (!nrow(dependency_df)) return(software_df)
  for (i in seq_len(nrow(dependency_df))) {
    idx <- match(dependency_df$package[[i]], software_df$software)
    if (is.na(idx)) next
    software_df$source[[idx]] <- dependency_df$source[[i]]
    software_df$repository[[idx]] <- dependency_df$repository[[i]]
    software_df$ref[[idx]] <- dependency_df$ref[[i]]
    software_df$dependency_status[[idx]] <- dependency_df$status[[i]]
    if ("path" %in% names(dependency_df) && nzchar(as.character(dependency_df$path[[i]]))) {
      software_df$path[[idx]] <- dependency_df$path[[i]]
    }
  }
  software_df
}

ril_config_from_env <- function(root = find_project_root()) {
  h2 <- env_num_vec("NG_RIL_STAGE_H2", c(0.03, 0.15, 0.40, 0.60, 0.70))
  if (length(h2) < 5L) h2 <- rep(h2, length.out = 5L)
  list(
    root = root,
    # Platform-specific default so a Windows-built `.Rlib` (e.g. synced via
    # OneDrive) is never reused or polluted by macOS/Linux installs. Windows keeps
    # the historical `.Rlib`; other platforms get their own tagged library. Override
    # with NG_RIL_R_LIB.
    r_lib = env_chr("NG_RIL_R_LIB", file.path(
      dirname(root),
      if (.Platform$OS.type == "windows") ".Rlib" else paste0(".Rlib-", R.version$os, "-", R.version$arch)
    )),
    install_missing_dependencies = env_bool("NG_RIL_INSTALL_MISSING_DEPS", TRUE),
    force_github_install = env_bool("NG_RIL_FORCE_GITHUB_INSTALL", FALSE),
    cran_repos = env_chr("NG_RIL_CRAN_REPOS", "https://cloud.r-project.org"),
    github_upgrade = env_chr("NG_RIL_GITHUB_UPGRADE", "never"),
    popvar_github_repo = env_chr("NG_POPVAR_GITHUB_REPO", "UMN-BarleyOatSilphium/PopVar"),
    popvar_github_ref = env_chr("NG_POPVAR_GITHUB_REF", "HEAD"),
    popvar_github_subdir = env_chr("NG_POPVAR_GITHUB_SUBDIR", ""),
    simplemating_github_repo = env_chr("NG_SIMPLEMATING_GITHUB_REPO", "Resende-Lab/SimpleMating"),
    simplemating_github_ref = env_chr("NG_SIMPLEMATING_GITHUB_REF", "HEAD"),
    simplemating_github_subdir = env_chr("NG_SIMPLEMATING_GITHUB_SUBDIR", ""),
    prediction_mode = env_chr("NG_RIL_PREDICTION_MODE", "trait_by_trait"),
    index_direction = env_chr("NG_RIL_INDEX_DIRECTION", "increase"),
    multi_trait_method = env_chr("NG_RIL_MULTI_TRAIT_METHOD", "weighted"),
    trait_names = env_csv("NG_RIL_TRAITS", "yield,protein,disease"),
    trait_directions = env_csv("NG_RIL_TRAIT_DIRECTIONS", "increase,increase,decrease"),
    trait_weights = env_num_vec("NG_RIL_TRAIT_WEIGHTS", c(0.50, 0.20, 0.30)),
    trait_economic_weights = env_num_vec("NG_RIL_TRAIT_ECONOMIC_WEIGHTS", c(1.00, 0.45, 1.20)),
    trait_desired_change = env_num_vec("NG_RIL_TRAIT_DESIRED_CHANGE", c(0.50, 0.20, 0.40)),
    trait_variances = env_num_vec("NG_RIL_TRAIT_VARIANCES", c(1.00, 0.60, 0.80)),
    trait_corA = env_num_vec("NG_RIL_TRAIT_CORA", c(
      1.00, 0.25, -0.35,
      0.25, 1.00, -0.20,
      -0.35, -0.20, 1.00
    )),
    threshold_policy = env_chr("NG_RIL_THRESHOLD_POLICY", "soft"),
    threshold_penalty_weight = env_num("NG_RIL_THRESHOLD_PENALTY_WEIGHT", 1.0),
    threshold_penalty_autoscale = env_bool("NG_RIL_THRESHOLD_PENALTY_AUTOSCALE", TRUE),
    seed = env_int("NG_RIL_SEED", 20260630L),
    # Breeders want many recurrent-selection cycles (to watch gain accrue and
    # diversity erode) and enough replicates for stable means. Reps are independent
    # and run in parallel; cycles build on each other and stay sequential.
    reps = env_int("NG_RIL_REPS", 10L),
    cycles = env_int("NG_RIL_CYCLES", 25L),
    # Parallelism (fork-based on macOS/Linux; serial fallback on Windows).
    #   scope "reps"    -> parallelize the independent replicate loop (default)
    #   scope "methods" -> reps serial, parallelize methods within each cycle
    #   scope "none"    -> fully serial
    parallel_scope = env_chr("NG_RIL_PARALLEL_SCOPE", "reps"),
    parallel_workers = env_int("NG_RIL_PARALLEL_WORKERS", 0L),
    n_founders = env_int("NG_RIL_N_FOUNDERS", 80L),
    n_parents = env_int("NG_RIL_N_PARENTS", 30L),
    n_chr = env_int("NG_RIL_N_CHR", 5L),
    seg_sites = env_int("NG_RIL_SEG_SITES", 300L),
    snp_per_chr = env_int("NG_RIL_SNP_PER_CHR", 200L),
    qtl_per_chr = env_int("NG_RIL_QTL_PER_CHR", 20L),
    genome_length_m = env_num("NG_RIL_GENOME_LENGTH_M", 1.0),
    n_crosses = env_int("NG_RIL_N_CROSSES", 20L),
    f2_per_cross = env_int("NG_RIL_F2_PER_CROSS", 80L),
    bulk_size = env_int("NG_RIL_BULK_SIZE", 5L),
    within_family_selection_prop = env_num("NG_RIL_WITHIN_FAMILY_SELECTION_PROP", 0.10),
    between_family_selection_prop = env_num("NG_RIL_BETWEEN_FAMILY_SELECTION_PROP", 0.20),
    ayt_selection_prop = env_num("NG_RIL_AYT_SELECTION_PROP", 0.70),
    # RIL parents keep their natural residual heterozygosity (mostly homozygous with a few
    # unfixed outlier loci). Set TRUE only to model a DH program (finish parents to full
    # homozygosity); OFF by default so the default study is a true RIL program. Residual-het
    # loci are handled at scoring time (treated as missing for tools that require 0/2), not
    # by forcing whole-genome homozygosity here.
    make_dh_parents = env_bool("NG_RIL_MAKE_DH_PARENTS", FALSE),
    stage_heritabilities = h2[seq_len(5L)],
    effect_training_n = env_int("NG_RIL_EFFECT_TRAINING_N", 120L),
    training_self_generations = env_int("NG_RIL_TRAINING_SELF_GENERATIONS", 5L),
    effect_training_h2 = env_num("NG_RIL_EFFECT_TRAINING_H2", 0.50),
    effect_h2_prior = env_num("NG_RIL_EFFECT_H2_PRIOR", 0.50),
    effect_kfold = env_int("NG_RIL_EFFECT_KFOLD", 0L),
    parent_observed_h2 = env_num("NG_RIL_PARENT_OBSERVED_H2", 0.50),
    selection_prop = env_num("NG_RIL_SELECTION_PROP", 0.10),
    min_effect_reliability = env_num("NG_RIL_MIN_EFFECT_RELIABILITY", 0.20),
    recombination_model = env_chr("NG_RIL_RECOMBINATION_MODEL", "haldane"),
    use_cpp = env_bool("NG_RIL_USE_CPP", FALSE),
    alphasimr_threads = env_int("NG_RIL_ALPHASIMR_THREADS", 1L),
    assume_inbred_for_scoring = env_bool("NG_RIL_ASSUME_INBRED_FOR_SCORING", TRUE),
    retry_scoring_as_noninbred = env_bool("NG_RIL_RETRY_SCORING_AS_NONINBRED", TRUE),
    max_crosses_per_parent = env_int("NG_RIL_MAX_CROSSES_PER_PARENT", 6L),
    max_pair_kinship = env_num("NG_RIL_MAX_PAIR_KINSHIP", Inf),
    lambda_group = env_num("NG_RIL_LAMBDA_GROUP", 0.05),
    lambda_mating = env_num("NG_RIL_LAMBDA_MATING", 0.02),
    lambda_parent_use = env_num("NG_RIL_LAMBDA_PARENT_USE", 0.0),
    lambda_parent_use_mode = env_chr("NG_RIL_LAMBDA_PARENT_USE_MODE", "absolute"),
    # greedy_local is the default OCS solver for the study: it optimizes the same
    # objective as the MIP, is fast and robust, and is within tolerance of MIP in the
    # bake-off. This avoids lpSolve branch-and-bound blow-up on the large, degenerate
    # (inbred late-cycle) score landscapes a recurrent program produces. Set
    # NG_RIL_OPTIMIZER=mip_contribution to force the exact MIP (guarded/time-capped).
    optimizer = env_chr("NG_RIL_OPTIMIZER", "greedy_local"),
    local_iter = env_int("NG_RIL_LOCAL_ITER", 1000L),
    ocs_iter = env_int("NG_RIL_OCS_ITER", 5L),
    evol_solutions = env_int("NG_RIL_EVOL_SOLUTIONS", 100L),
    evol_iterations = env_int("NG_RIL_EVOL_ITERATIONS", 150L),
    evol_stop = env_int("NG_RIL_EVOL_STOP", 30L),
    simplemating_generation = env_int("NG_RIL_SIMPLEMATING_GENERATION", 6L),
    simplemating_threads = env_int("NG_RIL_SIMPLEMATING_THREADS", 1L),
    simplemating_min_cross = env_int("NG_RIL_SIMPLEMATING_MIN_CROSS", 1L),
    simplemating_culling_k = env_num("NG_RIL_SIMPLEMATING_CULLING_K", NA_real_),
    alphamate_mode = env_chr("NG_RIL_ALPHAMATE_MODE", "ModeOptTarget1"),
    alphamate_target_degree = env_num("NG_RIL_ALPHAMATE_TARGET_DEGREE", 45),
    alphamate_lambda_group = env_num("NG_RIL_ALPHAMATE_LAMBDA_GROUP", NA_real_),
    alphamate_executable = env_chr("NG_RIL_ALPHAMATE_EXE", Sys.getenv("NG_ALPHAMATE_EXE", unset = "")),
    alphamate_runtime_path = env_chr("NG_RIL_ALPHAMATE_RUNTIME_PATH", Sys.getenv("NG_ALPHAMATE_RUNTIME_PATH", unset = "")),
    alphamate_max_contributions = env_int("NG_RIL_ALPHAMATE_MAX_CONTRIBUTIONS", NA_integer_),
    alphamate_number_of_parents = env_int("NG_RIL_ALPHAMATE_NUMBER_OF_PARENTS", NA_integer_),
    alphamate_evol_solutions = env_int("NG_RIL_ALPHAMATE_EVOL_SOLUTIONS", 100L),
    alphamate_evol_iterations = env_int("NG_RIL_ALPHAMATE_EVOL_ITERATIONS", 1000L),
    alphamate_evol_stop = env_int("NG_RIL_ALPHAMATE_EVOL_STOP", 200L),
    alphamate_threads = env_int("NG_RIL_ALPHAMATE_THREADS", 1L),
    alphamate_keep_files = env_bool("NG_RIL_ALPHAMATE_KEEP_FILES", FALSE),
    methods = unique(env_csv(
      "NG_RIL_METHODS",
      paste(ril_benchmark_default_methods(), collapse = ",")
    )),
    baseline_method = env_chr("NG_RIL_BASELINE_METHOD", "var_simple_ocs"),
    accuracy_enable = env_bool("NG_RIL_ACCURACY_ENABLE", TRUE),
    accuracy_n_candidates = env_int("NG_RIL_ACCURACY_N_CANDIDATES", 60L),
    accuracy_n_realized = env_int("NG_RIL_ACCURACY_N_REALIZED", 40L),
    accuracy_family_size = env_int("NG_RIL_ACCURACY_FAMILY_SIZE", 30L),
    accuracy_self_generations = env_int("NG_RIL_ACCURACY_SELF_GENERATIONS", 5L),
    output_dir = env_chr("NG_RIL_OUTPUT_DIR", file.path(root, "results")),
    output_prefix = env_chr("NG_RIL_OUTPUT_PREFIX", "ril_breeding_program_benchmark")
  )
}

validate_ril_config <- function(cfg) {
  if (cfg$reps < 1L) stop("NG_RIL_REPS must be >= 1", call. = FALSE)
  if (cfg$cycles < 1L) stop("NG_RIL_CYCLES must be >= 1", call. = FALSE)
  if (cfg$n_parents < 4L) stop("NG_RIL_N_PARENTS must be >= 4", call. = FALSE)
  if (cfg$n_crosses < 1L) stop("NG_RIL_N_CROSSES must be >= 1", call. = FALSE)
  if (cfg$f2_per_cross < 2L) stop("NG_RIL_F2_PER_CROSS must be >= 2", call. = FALSE)
  if (!identical(cfg$recombination_model, "haldane") && !identical(cfg$recombination_model, "kosambi")) {
    stop("NG_RIL_RECOMBINATION_MODEL must be haldane or kosambi", call. = FALSE)
  }
  cfg$r_lib <- normalizePath(cfg$r_lib, winslash = "/", mustWork = FALSE)
  cfg$github_upgrade <- match.arg(tolower(trimws(cfg$github_upgrade)), c("never", "ask", "always"))
  cfg$prediction_mode <- match.arg(tolower(trimws(cfg$prediction_mode)), c("trait_by_trait", "index_as_trait"))
  cfg$index_direction <- match.arg(tolower(trimws(cfg$index_direction)), c("increase", "decrease"))
  if (identical(cfg$prediction_mode, "index_as_trait") && !identical(cfg$index_direction, "increase")) {
    stop("NG_RIL_INDEX_DIRECTION must be increase when using index_as_trait; direct indexes must already be oriented so larger is better", call. = FALSE)
  }
  cfg$multi_trait_method <- match.arg(tolower(trimws(cfg$multi_trait_method)), c("auto", "weighted", "economic_index", "desired_gain"))
  cfg$threshold_policy <- match.arg(tolower(trimws(cfg$threshold_policy)), c("soft", "strict"))
  cfg$trait_spec <- resolve_benchmark_trait_spec(cfg)
  n_traits <- nrow(cfg$trait_spec)
  corA <- suppressWarnings(as.numeric(cfg$trait_corA))
  if (length(corA) == 1L) {
    corA <- matrix(corA, n_traits, n_traits)
    diag(corA) <- 1
  } else if (length(corA) == n_traits * n_traits) {
    corA <- matrix(corA, n_traits, n_traits, byrow = TRUE)
  } else {
    stop("NG_RIL_TRAIT_CORA must be length 1 or n_traits^2", call. = FALSE)
  }
  corA[!is.finite(corA)] <- 0
  corA <- (corA + t(corA)) / 2
  diag(corA) <- 1
  eig <- eigen(corA, symmetric = TRUE, only.values = TRUE)$values
  if (min(eig) <= 1e-8) {
    corA <- corA + diag(abs(min(eig)) + 1e-6, n_traits)
    D <- sqrt(diag(corA))
    corA <- corA / tcrossprod(D)
    diag(corA) <- 1
  }
  dimnames(corA) <- list(cfg$trait_spec$trait, cfg$trait_spec$trait)
  cfg$trait_corA_matrix <- corA
  cfg$lambda_parent_use_mode <- match.arg(tolower(trimws(cfg$lambda_parent_use_mode)), c("absolute", "adaptive"))
  cfg$alphamate_mode <- match.arg(cfg$alphamate_mode, c("ModeOptTarget1", "ModeMaxCriterion", "ModeMinCoancestry"))
  if (!is.finite(cfg$alphasimr_threads) || cfg$alphasimr_threads < 1L) cfg$alphasimr_threads <- 1L
  if (!is.finite(cfg$alphamate_lambda_group)) cfg$alphamate_lambda_group <- NULL
  if (!is.finite(cfg$alphamate_max_contributions)) cfg$alphamate_max_contributions <- NULL
  if (!is.finite(cfg$alphamate_number_of_parents)) cfg$alphamate_number_of_parents <- NULL
  cfg$parallel_scope <- match.arg(tolower(trimws(cfg$parallel_scope)), c("reps", "methods", "none"))
  if (!is.finite(cfg$parallel_workers) || cfg$parallel_workers < 0L) cfg$parallel_workers <- 0L
  cfg
}

# Fork-based parallel map on macOS/Linux; deterministic serial fallback on Windows
# (and whenever workers <= 1). Each element's work is independent and RNG is keyed by
# stable_seed(), so results are identical regardless of worker count or OS.
parallel_lapply <- function(x, fun, workers) {
  workers <- as.integer(workers)
  if (!is.finite(workers) || workers <= 1L || length(x) <= 1L) return(lapply(x, fun))
  if (.Platform$OS.type == "windows" || !requireNamespace("parallel", quietly = TRUE)) {
    return(lapply(x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = min(workers, length(x)), mc.preschedule = FALSE)
}

resolve_workers <- function(cfg, n_units) {
  w <- cfg$parallel_workers
  avail <- tryCatch(parallel::detectCores(), error = function(e) 1L)
  if (!is.finite(avail) || avail < 1L) avail <- 1L
  if (!is.finite(w) || w <= 0L) w <- max(1L, avail - 1L)  # leave one core free
  min(as.integer(w), as.integer(n_units), as.integer(avail))
}

compare_with_baseline <- function(metrics, baseline_method) {
  if (!nrow(metrics) || !(baseline_method %in% metrics$method)) return(data.frame())
  overall <- stats::aggregate(
    cbind(mean_index_gv, max_index_gv, top10_index_gv, var_index_gv) ~ method + cycle,
    metrics,
    mean
  )
  baseline <- overall[overall$method == baseline_method, , drop = FALSE]
  rows <- lapply(seq_len(nrow(overall)), function(i) {
    b <- baseline[baseline$cycle == overall$cycle[[i]], , drop = FALSE]
    if (nrow(b) != 1L || identical(overall$method[[i]], baseline_method)) return(NULL)
    data.frame(
      method = overall$method[[i]],
      baseline_method = baseline_method,
      cycle = overall$cycle[[i]],
      delta_mean_index_gv = overall$mean_index_gv[[i]] - b$mean_index_gv[[1]],
      delta_top10_index_gv = overall$top10_index_gv[[i]] - b$top10_index_gv[[1]],
      delta_max_index_gv = overall$max_index_gv[[i]] - b$max_index_gv[[1]],
      delta_var_index_gv = overall$var_index_gv[[i]] - b$var_index_gv[[1]],
      stringsAsFactors = FALSE
    )
  })
  bind_rows_fill(rows)
}

write_outputs <- function(result, cfg) {
  dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- file.path(cfg$output_dir, cfg$output_prefix)
  write.csv(result$config, paste0(prefix, "_config.csv"), row.names = FALSE)
  write.csv(result$trait_spec, paste0(prefix, "_trait_spec.csv"), row.names = FALSE)
  write.csv(result$dependency_setup, paste0(prefix, "_dependency_setup.csv"), row.names = FALSE)
  write.csv(result$software_versions, paste0(prefix, "_software_versions.csv"), row.names = FALSE)
  write.csv(result$method_availability, paste0(prefix, "_method_availability.csv"), row.names = FALSE)
  write.csv(result$metrics, paste0(prefix, "_metrics_by_cycle.csv"), row.names = FALSE)
  write.csv(result$stage_summary, paste0(prefix, "_stage_summary.csv"), row.names = FALSE)
  write.csv(result$family_summary, paste0(prefix, "_family_summary.csv"), row.names = FALSE)
  write.csv(result$selected_crosses, paste0(prefix, "_selected_crosses.csv"), row.names = FALSE)
  write.csv(result$selection_summary, paste0(prefix, "_selection_summary.csv"), row.names = FALSE)
  write.csv(result$comparison, paste0(prefix, "_comparison_vs_baseline.csv"), row.names = FALSE)
  if (!is.null(result$prediction_accuracy) && nrow(result$prediction_accuracy)) {
    write.csv(result$prediction_accuracy, paste0(prefix, "_prediction_accuracy.csv"), row.names = FALSE)
  }
  if (!is.null(result$selected_cross_accuracy) && nrow(result$selected_cross_accuracy)) {
    write.csv(result$selected_cross_accuracy, paste0(prefix, "_selected_cross_accuracy.csv"), row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()), paste0(prefix, "_session_info.txt"))
  normalizePath(cfg$output_dir, winslash = "/", mustWork = TRUE)
}

# Per-cycle breeder metrics: genetic gain (cumulative delta from the cycle-0 base),
# population mean, and the genetic-diversity signals breeders watch erode over cycles.
# var_index_gv is the population additive genetic variance of the index INCLUDING LD
# (Bulmer), not the LD-free genic variance. He/polymorphism/MAF fall as alleles fix;
# mean_parent_relationship (VanRaden G off-diagonal ~ 2 x coancestry coefficient)
# rises as relatedness accumulates.
ril_metrics_row <- function(pop, rep, method, cycle, cfg, base_mean, base_best, marker_assets, sim_param) {
  trait_gv <- gv_matrix(pop, cfg)
  g <- selection_index_from_trait_values(trait_gv, cfg)
  geno <- population_genotypes(pop, sim_param, marker_assets$marker_names)
  div <- population_diversity(geno, g)
  cbind(data.frame(
    rep = rep,
    method = method,
    cycle = cycle,
    prediction_mode = cfg$prediction_mode,
    multi_trait_method = cfg$multi_trait_method,
    index_direction = "increase",
    n_ind = length(g),
    mean_index_gv = mean(g),
    max_index_gv = max(g),
    top10_index_gv = top_prop_mean(g, 0.10),
    var_index_gv = stats::var(g),
    genetic_gain_index = mean(g) - base_mean,
    max_gain_index = max(g) - base_mean,
    transgressive_rate = mean(g > base_best),
    expected_heterozygosity = div$expected_heterozygosity,
    prop_polymorphic_markers = div$prop_polymorphic_markers,
    mean_maf = div$mean_maf,
    mean_parent_relationship = div$mean_parent_relationship,
    max_parent_relationship = div$max_parent_relationship,
    stringsAsFactors = FALSE
  ), trait_mean_columns(trait_gv, "mean_gv", cfg))
}

# One independent replicate: founder haplotypes -> DH base parents -> recurrent
# selection cycles (sequential; each builds on the last). Methods within a cycle are
# independent and may run in parallel. Returns per-rep result frames plus the shared
# prediction-accuracy table for this replicate.
run_single_rep <- function(rep, cfg, active_methods) {
  # Pin BLAS/OpenMP to a single thread inside each forked replicate worker. Reps run
  # in parallel via mclapply and each does dense matrix ops (scoring); without this a
  # multithreaded BLAS spawns (workers)^2 threads that thrash the cores and serialize
  # the run. AlphaSimR is already pinned via sim_param$nThreads.
  Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1L), silent = TRUE)
  }
  founder <- with_rng_seed(
    stable_seed(rep, "founder", seed = cfg$seed),
    AlphaSimR::quickHaplo(cfg$n_founders, cfg$n_chr, cfg$seg_sites, genLen = cfg$genome_length_m)
  )
  sim_param <- AlphaSimR::SimParam$new(founder)
  sim_param$nThreads <- max(1L, as.integer(cfg$alphasimr_threads))
  sim_param$addTraitA(
    nQtlPerChr = cfg$qtl_per_chr,
    mean = rep(0, nrow(cfg$trait_spec)),
    var = cfg$trait_spec$simulated_variance,
    corA = cfg$trait_corA_matrix,
    name = cfg$trait_spec$trait
  )
  sim_param$addSnpChip(nSnpPerChr = cfg$snp_per_chr)
  marker_assets <- prepare_marker_assets(sim_param)

  base_parents <- make_base_parents(founder, sim_param, cfg, rep)
  branches <- stats::setNames(lapply(active_methods, function(method) {
    p <- base_parents
    p@id <- paste0(method, "_R", rep, "_C0_P", seq_len(AlphaSimR::nInd(p)))
    p
  }), active_methods)

  base_idx <- gv_vec(base_parents, cfg)
  base_mean <- mean(base_idx)
  base_best <- max(base_idx)

  metrics_out <- list()
  stages_out <- list()
  families_out <- list()
  selected_out <- list()
  selection_summary_out <- list()

  # cycle-0 baseline (identical starting pool) so gain/diversity trajectories show
  # the starting point every method departs from.
  for (method in active_methods) {
    metrics_out[[paste(rep, method, 0, sep = "_")]] <-
      ril_metrics_row(branches[[method]], rep, method, 0L, cfg, base_mean, base_best, marker_assets, sim_param)
  }

  method_workers <- if (identical(cfg$parallel_scope, "methods")) resolve_workers(cfg, length(active_methods)) else 1L

  for (cycle in seq_len(cfg$cycles)) {
    run_method <- function(method) {
      parent_pop <- branches[[method]]
      selection <- tryCatch(
        select_crosses_for_method(parent_pop, sim_param, marker_assets, cfg, rep, method, cycle),
        error = function(e) {
          warning("Skipping method ", method, " at rep ", rep, " cycle ", cycle, ": ", conditionMessage(e), call. = FALSE)
          NULL
        }
      )
      if (is.null(selection)) return(NULL)
      advanced <- advance_ril_pipeline(
        parent_pop = parent_pop,
        selected_crosses = selection$selected,
        sim_param = sim_param, cfg = cfg, rep = rep, method = method, cycle = cycle
      )
      list(
        method = method,
        next_parents = advanced$next_parents,
        metrics = ril_metrics_row(advanced$next_parents, rep, method, cycle, cfg, base_mean, base_best, marker_assets, sim_param),
        stage = advanced$stage_summary,
        family = advanced$family_summary,
        selected = selection$selected,
        selection_summary = selection$selection_summary
      )
    }
    results <- parallel_lapply(active_methods, run_method, method_workers)
    # Methods are independent within a cycle: collect all, THEN advance branches so
    # every method saw the same start-of-cycle pool.
    for (res in results) {
      if (is.null(res)) next
      k <- paste(rep, res$method, cycle, sep = "_")
      branches[[res$method]] <- res$next_parents
      metrics_out[[k]] <- res$metrics
      stages_out[[k]] <- res$stage
      families_out[[k]] <- res$family
      selected_out[[k]] <- res$selected
      selection_summary_out[[k]] <- res$selection_summary
    }
  }

  accuracy <- tryCatch(
    evaluate_shared_accuracy(base_parents, sim_param, marker_assets, cfg, rep, active_methods),
    error = function(e) {
      warning("prediction-accuracy evaluation failed at rep ", rep, ": ", conditionMessage(e), call. = FALSE)
      data.frame()
    }
  )

  list(
    metrics = bind_rows_fill(metrics_out),
    stage_summary = bind_rows_fill(stages_out),
    family_summary = bind_rows_fill(families_out),
    selected_crosses = bind_rows_fill(selected_out),
    selection_summary = bind_rows_fill(selection_summary_out),
    prediction_accuracy = accuracy
  )
}

run_ril_breeding_program_benchmark <- function(cfg = ril_config_from_env()) {
  # SimpleMating imports rgl, which aborts at load time on headless machines with no
  # OpenGL/X11 (common on servers and CI). rgl's null-device mode loads without a
  # display; enable it by default so the SimpleMating comparator is usable headless.
  # Forked replicate workers inherit this environment.
  if (!nzchar(Sys.getenv("RGL_USE_NULL"))) Sys.setenv(RGL_USE_NULL = "TRUE")
  cfg <- validate_ril_config(cfg)
  dependency_df <- install_required_benchmark_dependencies(cfg)
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) {
    stop("AlphaSimR is required for the RIL breeding-program benchmark", call. = FALSE)
  }
  .libPaths(c(normalizePath(cfg$r_lib, mustWork = FALSE), .libPaths()))
  source(file.path(cfg$root, "R", "load.R"))
  ng_load(cfg$root, use_cpp = cfg$use_cpp, verbose = FALSE)
  suppressPackageStartupMessages(library(AlphaSimR))

  method_info <- method_availability(cfg$methods, cfg)
  active_methods <- method_info$method[method_info$active]
  if (!length(active_methods)) stop("No active benchmark methods remain after dependency checks", call. = FALSE)

  config_value <- function(x) {
    if (is.data.frame(x)) return(paste(utils::capture.output(utils::str(x, give.attr = FALSE)), collapse = " "))
    if (is.matrix(x)) return(paste(as.numeric(x), collapse = ","))
    paste(x, collapse = ",")
  }
  config_df <- data.frame(
    name = names(cfg),
    value = vapply(cfg, config_value, character(1)),
    stringsAsFactors = FALSE
  )
  version_df <- add_dependency_provenance(software_versions(cfg), dependency_df)

  rep_workers <- if (identical(cfg$parallel_scope, "reps")) resolve_workers(cfg, cfg$reps) else 1L
  method_workers <- if (identical(cfg$parallel_scope, "methods")) resolve_workers(cfg, length(active_methods)) else 1L
  message(sprintf(
    "RIL benchmark: %d reps x %d cycles x %d methods | parallel scope=%s (rep workers=%d, method workers=%d)",
    cfg$reps, cfg$cycles, length(active_methods), cfg$parallel_scope, rep_workers, method_workers
  ))

  rep_results <- parallel_lapply(seq_len(cfg$reps), function(rep) {
    message("  replicate ", rep, " of ", cfg$reps)
    run_single_rep(rep, cfg, active_methods)
  }, rep_workers)

  gather <- function(field) bind_rows_fill(lapply(rep_results, function(r) r[[field]]))
  metrics <- gather("metrics")
  selected_crosses <- gather("selected_crosses")
  family_summary <- gather("family_summary")

  result <- list(
    config = config_df,
    trait_spec = cfg$trait_spec,
    dependency_setup = dependency_df,
    software_versions = version_df,
    method_availability = method_info,
    metrics = metrics,
    stage_summary = gather("stage_summary"),
    family_summary = family_summary,
    selected_crosses = selected_crosses,
    selection_summary = gather("selection_summary"),
    comparison = compare_with_baseline(metrics, cfg$baseline_method),
    prediction_accuracy = gather("prediction_accuracy"),
    selected_cross_accuracy = selected_cross_accuracy(selected_crosses, family_summary)
  )
  out_dir <- write_outputs(result, cfg)
  message("Wrote RIL breeding-program benchmark outputs to: ", out_dir)
  message("Output prefix: ", cfg$output_prefix)
  result
}

is_this_script <- function() {
  args <- commandArgs(FALSE)
  any(grepl("run_ril_breeding_program_benchmark\\.R$", args))
}

if (is_this_script()) {
  invisible(run_ril_breeding_program_benchmark())
}
