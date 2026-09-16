# One field answers "is this column an index?" -- and it is value_kind.
#
# Two parameters asserted the same fact. `prediction_mode = "index_as_trait"` said it at
# the RUN level (and did nothing but rename the trait to "selection_index"), while
# value_kind says it PER TRAIT. Two sources of truth that can disagree is the defect
# this project already fixed once: posterior_predictions$mean_source and
# effect_summary$mean_source both answered "which basis?" and gave opposite answers, and
# the fix was not to validate them against each other -- it was to make one field answer
# the question.
#
# So index_as_trait is deprecated and MAPS onto value_kind rather than co-existing with
# it. Mapping beats guarding: downstream there is only one mechanism, so there is nothing
# left to contradict and no contradiction check to maintain. Same pattern the package
# already uses for assume_inbred -> parent_type and "le" -> parent_distance.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11); n <- 40L; m <- 40L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
gv <- as.numeric(g %*% c(rnorm(8, 0, 1), rep(0, m - 8)))
idx <- 500 + 37 * as.numeric(scale(gv + rnorm(n, 0, 0.2 * stats::sd(gv))))
pheno <- data.frame(NAME = ids, my_index = idx, stringsAsFactors = FALSE)
mm <- data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                 bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6)
base <- list(
  phenotype = pheno, genotype = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = mm, id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr",
  map_pos_col = "bp", bp_per_cm = 1e6, n_crosses = 8L, write_outputs = FALSE,
  write_figures = FALSE, run_posterior_prediction = FALSE, seed = 5L)

# ---- 1. the deprecated mode warns, in the PARENT process -----------------
# A warning raised inside an mclapply worker never reaches the caller (see the 0.32.0
# advisory finding), so a deprecation emitted during fitting would be silent. This one
# is raised at config time and must be catchable here.
w <- NULL
r_dep <- withCallingHandlers(
  do.call(ng_run_cross_prediction, c(base, list(
    trait_direction = data.frame(trait = "my_index", column = "my_index",
                                 direction = "increase", stringsAsFactors = FALSE),
    prediction_mode = "index_as_trait", index_col = "my_index"))),
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
stopifnot(!is.null(w))
dep <- grep("index_as_trait", w, value = TRUE)
stopifnot(length(dep) > 0L)
stopifnot(grepl("deprecat", dep[[1L]], ignore.case = TRUE))
stopifnot(grepl("value_kind", dep[[1L]], fixed = TRUE))   # names the replacement

# ---- 2. it MAPS: the deprecated path now transforms and reports as an index
es <- r_dep$effect_summary
stopifnot(identical(as.character(es$value_kind[[1L]]), "index"))
stopifnot(identical(as.character(es$index_transform[[1L]]), "inverse_normal_rank"))
stopifnot(identical(as.character(es$trait[[1L]]), "selection_index"))  # rename preserved

# ---- 3. the mapping is FAITHFUL, not approximate -------------------------
# The whole claim of a deprecation-by-mapping is that the old spelling produces exactly
# what the new one does. Compare the numbers, since the trait NAME differs by design.
r_new <- do.call(ng_run_cross_prediction, c(base, list(
  trait_direction = data.frame(trait = "my_index", column = "my_index",
                               direction = "increase", value_kind = "index",
                               stringsAsFactors = FALSE))))
stopifnot(identical(paste(r_dep$selected_crosses$parent1, r_dep$selected_crosses$parent2),
                    paste(r_new$selected_crosses$parent1, r_new$selected_crosses$parent2)))
stopifnot(isTRUE(all.equal(r_dep$candidate_crosses$selection_index_mean,
                           r_new$candidate_crosses$my_index_mean)))
stopifnot(isTRUE(all.equal(r_dep$candidate_crosses$selection_index_value,
                           r_new$candidate_crosses$my_index_value)))

# ---- 4. the modern spelling does NOT warn -------------------------------
w2 <- NULL
invisible(withCallingHandlers(
  do.call(ng_run_cross_prediction, c(base, list(
    trait_direction = data.frame(trait = "my_index", column = "my_index",
                                 direction = "increase", value_kind = "index",
                                 stringsAsFactors = FALSE)))),
  warning = function(cond) { w2 <<- c(w2, conditionMessage(cond)); invokeRestart("muffleWarning") }))
stopifnot(!length(grep("index_as_trait", w2 %||% character())))

# ---- 5. a direction table with no value_kind column still runs ----------
# Every config written before this feature existed must be unaffected.
r_old <- do.call(ng_run_cross_prediction, c(base, list(
  trait_direction = data.frame(trait = "my_index", column = "my_index",
                               direction = "increase", stringsAsFactors = FALSE))))
stopifnot(identical(as.character(r_old$effect_summary$value_kind[[1L]]), "trait"))
stopifnot(is.na(r_old$effect_summary$index_transform[[1L]]))

cat("index_as_trait_deprecated: PASS\n")
