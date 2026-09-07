helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 0.28.0 -- ORDER INDEPENDENCE OF PER-TRAIT RESULTS.
#
# The bug: ng_cp__stage_predict()'s per-trait loop derived both of its seeds from the trait's
# ROW POSITION `i` in the breeder's direction file:
#   * `seed + i - 1L` reached ng_choose_ridge_lambda() -> set.seed(); sample(), i.e. the k-fold
#     CV partition behind ridge lambda selection. A different partition can select a different
#     lambda, which changes the fitted marker effects and hence the progeny variance.
#   * `seed + 1000L + i - 1L` drove the per-trait posterior draws.
# So simply moving a trait from row 1 to row 2 of the direction file changed its science.
#
# The panel below is chosen (data seed 346) so the second trait's lambda genuinely flips between
# consecutive seeds. Measured on 0.27.0, `disease_vpm` averaged 2.5124748 when disease was listed
# FIRST and 3.197428e-07 when it was listed SECOND -- seven orders of magnitude -- and the
# selected crossing plan differed. This test pins that shut.

set.seed(346L)
n <- 24L; m <- 40L
X <- matrix(2L * rbinom(n * m, 1L, 0.5), nrow = n, ncol = m,
            dimnames = list(sprintf("P%02d", seq_len(n)), sprintf("M%02d", seq_len(m))))
disease <- as.numeric(X %*% rnorm(m, 0, 0.25)) + rnorm(n, 0, 2)
set.seed(11L)
yield <- 60 + as.numeric(X %*% rnorm(m, 0, 0.20)) + rnorm(n, 0, 2)

tmp <- tempfile("ng_trait_order_"); dir.create(tmp, recursive = TRUE)
files <- list(pheno = file.path(tmp, "phenotype.csv"), geno = file.path(tmp, "genotype.csv"),
              map = file.path(tmp, "map.csv"))
write.csv(data.frame(NAME = rownames(X), X, check.names = FALSE, stringsAsFactors = FALSE),
          files$geno, row.names = FALSE, quote = FALSE)
write.csv(data.frame(NAME = rownames(X), yield = yield, disease = disease,
                     stringsAsFactors = FALSE),
          files$pheno, row.names = FALSE, quote = FALSE)
write.csv(data.frame(SNP_code = colnames(X), Chromosome = rep(1:2, each = m / 2L),
                     Position_BP = rep(seq_len(m / 2L) * 3e6, 2), stringsAsFactors = FALSE),
          files$map, row.names = FALSE, quote = FALSE)

direction <- data.frame(Trait = c("yield", "disease"), PhenotypeColumn = c("yield", "disease"),
                        Selection_direction = c("increase", "decrease"),
                        stringsAsFactors = FALSE)
write_direction <- function(rows, name) {
  path <- file.path(tmp, name)
  write.csv(direction[rows, , drop = FALSE], path, row.names = FALSE, quote = FALSE)
  path
}
dir_yd <- write_direction(c(1L, 2L), "dir_yield_disease.csv")
dir_dy <- write_direction(c(2L, 1L), "dir_disease_yield.csv")
dir_y  <- write_direction(1L, "dir_yield.csv")
dir_d  <- write_direction(2L, "dir_disease.csv")

run <- function(direction_file, ...) {
  ng_run_cross_prediction(
    phenotype_file = files$pheno, genotype_file = files$geno, map_file = files$map,
    direction_file = direction_file,
    phenotype_id_col = "NAME", genotype_id_col = "NAME",
    direction_trait_col = "Trait", direction_column_col = "PhenotypeColumn",
    direction_direction_col = "Selection_direction",
    map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_bp_col = "Position_BP",
    map_position_unit = "bp", bp_per_cm = 1e6,
    prediction_mode = "trait_by_trait", trait_value_metric = "usefulness",
    uc_variance_source = "pmv", progeny = "DH", recomb_model = "haldane",
    selection_prop = 0.20, duplicate_action = "none", posterior_method = "closed_form",
    n_iter = 120L, burn_in = 10L, n_crosses = 5, max_crosses_per_parent = 3,
    optimizer = "greedy_local", write_outputs = FALSE, write_figures = FALSE,
    parent_type = "inbred", use_cpp = FALSE, seed = 909L, ...)
}

pair_key <- function(x) paste(x$parent1, x$parent2, sep = "|")
align_to <- function(x, reference) x[match(pair_key(reference), pair_key(x)), , drop = FALSE]

# Every numeric column that belongs to `trait` must survive a permutation of the direction file
# bit-for-bit. tolerance = 0 -- the requirement is identity, not closeness.
assert_same_trait <- function(a, b, trait, label) {
  cols <- grep(paste0("^", trait, "_"), names(a), value = TRUE)
  cols <- cols[vapply(cols, function(cc) is.numeric(a[[cc]]), logical(1L))]
  stopifnot(length(cols) > 0L)
  a <- align_to(a, b)
  stopifnot(!anyNA(pair_key(a)))
  bad <- character(0L)
  for (cc in cols) {
    x <- a[[cc]]; y <- b[[cc]]
    if (all(is.na(x)) && all(is.na(y))) next
    if (!isTRUE(all.equal(x, y, tolerance = 0))) bad <- c(bad, cc)
  }
  if (length(bad)) {
    stop(sprintf("%s: %s columns depend on the trait's row position: %s",
                 label, trait, paste(bad, collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("  %-46s %2d %s columns identical (tolerance = 0)\n", label, length(cols), trait))
}

# ---- 1. posterior OFF: permuting the direction file changes nothing -------------------------
res_yd <- run(dir_yd)
res_dy <- run(dir_dy)
cat("posterior off:\n")
assert_same_trait(res_yd$candidate_crosses, res_dy$candidate_crosses, "disease", "candidate_crosses")
assert_same_trait(res_yd$candidate_crosses, res_dy$candidate_crosses, "yield",   "candidate_crosses")
cat(sprintf("  disease_vpm mean: yield-first = %.10g, disease-first = %.10g\n",
            mean(res_yd$candidate_crosses$disease_vpm),
            mean(res_dy$candidate_crosses$disease_vpm)))

# The plan itself, and the index the plan is ranked on, are order-invariant too.
plan_key <- function(r) sort(pair_key(r$selected_crosses))
stopifnot(identical(plan_key(res_yd), plan_key(res_dy)))
mt_yd <- align_to(res_yd$candidate_crosses, res_dy$candidate_crosses)
stopifnot(isTRUE(all.equal(mt_yd$multi_trait_score, res_dy$candidate_crosses$multi_trait_score,
                           tolerance = 0)))

# ---- 2. posterior ON: the posterior columns are order-invariant as well ----------------------
post_yd <- run(dir_yd, run_posterior_prediction = TRUE)
post_dy <- run(dir_dy, run_posterior_prediction = TRUE)
cat("posterior on:\n")
assert_same_trait(post_yd$candidate_crosses, post_dy$candidate_crosses, "disease", "candidate_crosses")
assert_same_trait(post_yd$candidate_crosses, post_dy$candidate_crosses, "yield",   "candidate_crosses")

# the whole per-trait posterior table, not just the two columns that ride the cross table
for (trait in c("yield", "disease")) {
  a <- post_yd$posterior_predictions[[trait]]
  b <- post_dy$posterior_predictions[[trait]]
  stopifnot(!is.null(a), !is.null(b))
  a <- align_to(a, b)
  num <- names(b)[vapply(b, is.numeric, logical(1L))]
  bad <- num[!vapply(num, function(cc) isTRUE(all.equal(a[[cc]], b[[cc]], tolerance = 0)),
                     logical(1L))]
  if (length(bad)) {
    stop(sprintf("posterior_predictions[['%s']] depends on row position: %s",
                 trait, paste(bad, collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("  %-46s %2d %s columns identical (tolerance = 0)\n",
              "posterior_predictions", length(num), trait))
}

# The INDEX posterior (R/32, ng_posterior_multitrait_cross_predict) had the same defect
# (`seed + j`) and is now keyed on the trait name.
pm_yd <- align_to(post_yd$posterior_multitrait, post_dy$posterior_multitrait)
pm_dy <- post_dy$posterior_multitrait
stopifnot(!is.null(pm_dy))
pm_num <- names(pm_dy)[vapply(pm_dy, is.numeric, logical(1L))]
pm_bad <- pm_num[!vapply(pm_num, function(cc) isTRUE(all.equal(pm_yd[[cc]], pm_dy[[cc]],
                                                               tolerance = 0)), logical(1L))]
if (length(pm_bad)) {
  stop("posterior_multitrait depends on trait order: ", paste(pm_bad, collapse = ", "),
       call. = FALSE)
}
cat(sprintf("  %-46s %2d index-posterior columns identical (tolerance = 0)\n",
            "posterior_multitrait", length(pm_num)))

# ---- 3. a trait's DETERMINISTIC results do not depend on its neighbours ----------------------
# Sharing one CV fold partition across traits (rather than one per row position) also makes a
# multi-trait run agree with the corresponding single-trait run, and keeps single-trait runs
# bit-identical to 0.27.0.
solo_d <- run(dir_d)
solo_y <- run(dir_y)
for (trait in c("disease", "yield")) {
  solo <- if (identical(trait, "disease")) solo_d else solo_y
  cols <- grep(paste0("^", trait, "_"), names(solo$candidate_crosses), value = TRUE)
  cols <- setdiff(cols, grep("_post", cols, value = TRUE))
  cols <- cols[vapply(cols, function(cc) is.numeric(solo$candidate_crosses[[cc]]), logical(1L))]
  a <- align_to(solo$candidate_crosses, res_dy$candidate_crosses)
  bad <- cols[!vapply(cols, function(cc) {
    x <- a[[cc]]; y <- res_dy$candidate_crosses[[cc]]
    (all(is.na(x)) && all(is.na(y))) || isTRUE(all.equal(x, y, tolerance = 0))
  }, logical(1L))]
  if (length(bad)) {
    stop(sprintf("single-trait vs multi-trait %s mismatch: %s", trait,
                 paste(bad, collapse = ", ")), call. = FALSE)
  }
  cat(sprintf("  single-trait run reproduces multi-trait %-8s %2d columns (tolerance = 0)\n",
              trait, length(cols)))
}

# ---- 4. the seed helper itself ---------------------------------------------------------------
# Stable, explicit, and not derived from any hashing internal that could change between R
# versions or platforms. These literals are the values of the documented byte polynomial
# h <- (h * 131 + byte) %% (2^31 - 1) and must not drift.
stopifnot(identical(ng_name_hash32(""), 0))
stopifnot(identical(ng_name_hash32("a"), 97))
stopifnot(identical(ng_name_hash32("ab"), 97 * 131 + 98))
stopifnot(identical(ng_name_hash32("abc"), (97 * 131 + 98) * 131 + 99))
stopifnot(length(ng_name_hash32(c("yield", "disease"))) == 2L)
h_yield <- ng_name_hash32("yield")
h_disease <- ng_name_hash32("disease")
stopifnot(h_yield != h_disease)
stopifnot(h_yield >= 0, h_yield < 2147483647)
# a long name must not overflow into a non-integer or a negative seed
long <- paste(rep("trait_with_a_very_long_name", 50L), collapse = "_")
s_long <- ng_trait_rng_seed(2147483000, long, salt = 1000L)
stopifnot(is.integer(s_long), !is.na(s_long), s_long >= 0L)
# distinct traits get distinct streams; the same trait always gets the same one
stopifnot(ng_trait_rng_seed(909L, "yield", 1000L) != ng_trait_rng_seed(909L, "disease", 1000L))
stopifnot(identical(ng_trait_rng_seed(909L, "yield", 1000L),
                    ng_trait_rng_seed(909L, "yield", 1000L)))
# the salt separates independent uses of the same trait
stopifnot(ng_trait_rng_seed(909L, "yield", 0L) != ng_trait_rng_seed(909L, "yield", 1000L))
# non-ASCII names are keyed off their UTF-8 bytes, so the value does not depend on the locale
stopifnot(ng_name_hash32("proteïne") > 0)
cat(sprintf("  ng_name_hash32: yield = %.0f, disease = %.0f (stable across sessions)\n",
            h_yield, h_disease))

cat("trait_order_invariance: OK\n")
