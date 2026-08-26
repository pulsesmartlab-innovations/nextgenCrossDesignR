helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# v0.3.0 preview — end-to-end multi-trait posterior cross prediction.
# Verifies that ng_posterior_multitrait_cross_predict() joins per-trait
# beta posteriors + posterior G + multi-trait index in a single call,
# returns per-cross multi-trait index CIs, and the rank-stability
# probabilities sum to N for each requested top_n_target.

set.seed(20260523L)
n <- 60L
m <- 120L
geno <- matrix(2L * rbinom(n * m, 1L, 0.45), nrow = n, ncol = m)
ids <- paste0("L", seq_len(n))
markers <- paste0("M", seq_len(m))
rownames(geno) <- ids; colnames(geno) <- markers

# True G for 3 traits; simulate phenotypes Y = geno %*% beta + noise.
G_true <- matrix(c(1.0, 0.4, -0.3,
                   0.4, 1.2,  0.1,
                  -0.3, 0.1,  0.6), nrow = 3L, byrow = TRUE)
B_true <- matrix(stats::rnorm(m * 3L), nrow = m) %*% chol(G_true)
E <- matrix(stats::rnorm(n * 3L, sd = 0.5), nrow = n)
Y <- geno %*% B_true + E
colnames(Y) <- c("yield", "disease", "lodging")
rownames(Y) <- ids

traits <- ng_multitrait_spec(
  trait = c("yield", "disease", "lodging"),
  direction = c("maximize", "minimize", "minimize"),
  economic_weight = c(2, 1, 1)
)

mk <- data.frame(
  marker = markers,
  chr = rep(1:3, length.out = m),
  pos_cm = rep(seq(0, 80, length.out = m / 3L), 3L)
)

# Use a small subset of crosses to keep the test fast.
pair_ids <- ids[seq_len(20L)]
pairs <- ng_make_pairs(pair_ids)
geno_use <- geno[pair_ids, , drop = FALSE]
Y_use <- Y[pair_ids, , drop = FALSE]

# Default mode = "mean" (fast); 30 draws is enough for a smoke check.
post_scores <- ng_posterior_multitrait_cross_predict(
  geno = geno_use, Y = Y_use, traits = traits, marker_map = mk,
  ids = rownames(geno_use), pairs = pairs,
  n_draws = 30L, posterior_method = "closed_form",
  genetic_covariance_method = "beta_posterior",
  genetic_covariance = G_true,
  index_method = "economic_index", value_mode = "mean",
  selection_prop = 0.20, target = "DH",
  recomb_model = "haldane", use_cpp = FALSE,
  top_n_targets = c(3L, 5L), seed = 1L
)

# Required columns
expected_cols <- c(
  "multi_trait_score_post_mean",
  "multi_trait_score_post_lower",
  "multi_trait_score_post_upper",
  "multitrait_posterior_topn_prob_3",
  "multitrait_posterior_topn_prob_5"
)
miss <- setdiff(expected_cols, names(post_scores))
if (length(miss)) stop("missing columns: ", paste(miss, collapse = ", "))
stopifnot(nrow(post_scores) == nrow(pairs))

# Posterior CI brackets the posterior mean.
stopifnot(all(post_scores$multi_trait_score_post_lower <=
              post_scores$multi_trait_score_post_mean + 1e-6))
stopifnot(all(post_scores$multi_trait_score_post_mean <=
              post_scores$multi_trait_score_post_upper + 1e-6))
# CI width is non-trivial (the posterior actually moves across draws).
ci_width <- post_scores$multi_trait_score_post_upper -
  post_scores$multi_trait_score_post_lower
if (max(ci_width) < 1e-3) {
  stop("multi_trait_score posterior CI widths are all near zero — sampler not propagating uncertainty")
}

# posterior_topn_prob columns are in [0, 1] and the row sums equal N.
for (N in c(3L, 5L)) {
  col <- paste0("multitrait_posterior_topn_prob_", N)
  stopifnot(all(post_scores[[col]] >= 0 & post_scores[[col]] <= 1))
  if (abs(sum(post_scores[[col]], na.rm = TRUE) - N) > 0.05) {
    stop(sprintf("Sum of multitrait_posterior_topn_prob_%d should equal %d, got %.3f",
                 N, N, sum(post_scores[[col]], na.rm = TRUE)))
  }
}

meta <- attr(post_scores, "posterior_multitrait")
stopifnot(identical(meta$index_method, "economic_index"))
stopifnot(identical(meta$value_mode, "mean"))
stopifnot(!is.null(meta$G_mean))
stopifnot(!is.null(meta$P_hat))

cat("posterior_multitrait_cross_predict: 5/5 checks passed\n")
cat(sprintf("  n_pairs=%d  n_draws=%d  index_method=%s\n",
            nrow(post_scores), meta$n_draws, meta$index_method))
cat(sprintf("  posterior CI width: median=%.3f  max=%.3f\n",
            stats::median(ci_width), max(ci_width)))
