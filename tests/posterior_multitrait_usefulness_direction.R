helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# DEFECT 6 (multi-trait audit, R/32_posterior_multitrait.R): under
# value_mode = "usefulness" the per-trait, per-draw value was built as
#   mean + i * sigma
# with NO direction sign, so a MINIMIZE trait (disease, lodging) was scored on the
# UNFAVOURABLE tail of its family -- within-family variance charged as a liability instead of
# credited as an opportunity -- and then value_z applied -1 on top of that. The wired
# single-trait convention is ng_run_cp_trait_value() (R/39): sign = +1 for maximize, -1 for
# minimize, applied ONLY to the i*sigma term, NEVER to the mean.
#
# The corrected term lives inside a draw loop, so this test pins it from the outside: with ONE
# trait and ONE draw, multi_trait_score is exactly the rank-normalization of that draw's
# per-trait usefulness value, so the function's output identifies the expression it used.

set.seed(404)
n_ind <- 30L
m <- 60L
geno <- matrix(sample(c(0, 2), n_ind * m, replace = TRUE), nrow = n_ind)
rownames(geno) <- paste0("L", seq_len(n_ind))
colnames(geno) <- paste0("M", seq_len(m))
marker_map <- data.frame(
  marker = colnames(geno),
  chr = rep(1:3, each = m / 3),
  pos_cm = rep(seq(0, 90, length.out = m / 3), times = 3),
  stringsAsFactors = FALSE
)
beta_d <- rnorm(m, sd = 0.6)
disease <- as.numeric(scale(geno %*% beta_d)) * 5 + 40 + rnorm(n_ind, sd = 1)
Y <- cbind(disease = disease)
rownames(Y) <- rownames(geno)

pairs <- ng_make_pairs(rownames(geno))[1:20, , drop = FALSE]

traits <- data.frame(trait = "disease", column = "disease", direction = "minimize",
                     weight = 1, stringsAsFactors = FALSE)

SEED <- 7L
N_DRAWS <- 1L
res <- suppressWarnings(ng_posterior_multitrait_cross_predict(
  geno = geno, Y = Y, traits = traits, marker_map = marker_map,
  pairs = pairs, n_draws = N_DRAWS, posterior_method = "closed_form",
  index_method = "weighted", value_mode = "usefulness",
  selection_prop = 0.10, target = "DH", seed = SEED
))

# ---- Rebuild the single draw the function used, outside the package -------------------------
# Same call shape R/32 makes for trait j = 1: seed = seed + j.
intensity <- ng_selection_intensity(0.10)
post_d <- ng_fit_ridge_effects_posterior(geno = geno, y = as.numeric(Y[, 1L]),
                                         ids = rownames(geno), lambda = NULL, kfold = 5L,
                                         n_draws = N_DRAWS, method = "closed_form",
                                         seed = SEED + 1L)
beta_s <- post_d$beta_draws[, 1L]
map_s <- ng_prepare_marker_map(marker_map, colnames(geno), model = "haldane")
sorted <- ng_sort_by_map(geno, beta_s, rep(0, ncol(geno)), marker_map = map_s)
scored_pair <- ng_dh_recomb_variance_pairs(
  geno = sorted$geno, beta = sorted$effects, beta_var = sorted$beta_var,
  marker_map = sorted$marker_map, ids = rownames(geno), pairs = pairs,
  window_cm = Inf, use_cpp = TRUE, recomb_model = "haldane", target = "DH")
fit_d <- post_d$fit
gebv <- fit_d$intercept + sum(fit_d$marker_mean * fit_d$beta) -
  sum(fit_d$marker_mean * beta_s) + as.numeric(geno %*% beta_s)
p1 <- match(pairs$parent1, rownames(geno)); p2 <- match(pairs$parent2, rownames(geno))
mp <- 0.5 * (gebv[p1] + gebv[p2])
sigma <- sqrt(pmax(scored_pair$pmv, 0))

value_fixed  <- mp - intensity * sigma   # correct for a MINIMIZE trait
value_broken <- mp + intensity * sigma   # what 0.25.0 produced regardless of direction

cat(sprintf("i*sigma range [%.5f, %.5f]; |fixed - broken| range [%.5f, %.5f]\n",
            min(intensity * sigma), max(intensity * sigma),
            min(abs(value_fixed - value_broken)), max(abs(value_fixed - value_broken))))
stopifnot(all(intensity * sigma > 1e-6))   # the sign has something to bite on

# The mean must NOT be signed: flipping the direction moves the value by exactly 2*i*sigma.
stopifnot(isTRUE(all.equal(value_broken - value_fixed, 2 * intensity * sigma)))

# The two expressions must ORDER the crosses differently, otherwise this test proves nothing.
stopifnot(!identical(order(value_fixed), order(value_broken)))

expect_fixed  <- as.numeric(ng_rank_normalize(value_fixed,  bigger_is_better = FALSE))
expect_broken <- as.numeric(ng_rank_normalize(value_broken, bigger_is_better = FALSE))
got <- as.numeric(res$multi_trait_score_post_mean)

stopifnot(isTRUE(all.equal(got, expect_fixed, tolerance = 1e-10)))
stopifnot(!isTRUE(all.equal(got, expect_broken, tolerance = 1e-6)))

# A MAXIMIZE trait must be unchanged by the fix: sign = +1 there.
traits_max <- traits; traits_max$direction <- "maximize"
res_max <- suppressWarnings(ng_posterior_multitrait_cross_predict(
  geno = geno, Y = Y, traits = traits_max, marker_map = marker_map,
  pairs = pairs, n_draws = N_DRAWS, posterior_method = "closed_form",
  index_method = "weighted", value_mode = "usefulness",
  selection_prop = 0.10, target = "DH", seed = SEED
))
stopifnot(isTRUE(all.equal(as.numeric(res_max$multi_trait_score_post_mean),
                           as.numeric(ng_rank_normalize(value_broken, bigger_is_better = TRUE)),
                           tolerance = 1e-10)))

cat("posterior_multitrait_usefulness_direction: OK\n")
