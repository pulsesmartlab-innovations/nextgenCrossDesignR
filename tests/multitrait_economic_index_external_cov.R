# Smith-Hazel economic index with USER-SUPPLIED external covariances. The external G/P arrive in
# RAW trait units, but the index is applied to value_z = sign*(raw-center)/scale. The fix maps the
# external G/P (and the economic weights) into the value_z space so the package's economic-index
# ranking reproduces the true raw Smith-Hazel index b = P^{-1} G a exactly. Regression: without an
# external covariance the (self-consistent) internal path is unchanged.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(21)
K <- 60L
scores <- data.frame(
  parent1 = sprintf("A%02d", seq_len(K)), parent2 = sprintf("B%02d", seq_len(K)),
  pred_yield   = rnorm(K, 60, 6),
  pred_disease = rnorm(K, 3, 1.2),
  stringsAsFactors = FALSE)
traits <- data.frame(trait = c("yield", "disease"), column = c("pred_yield", "pred_disease"),
                     direction = c("increase", "decrease"), economic_weight = c(2, 3),
                     stringsAsFactors = FALSE)
lab <- c("yield", "disease")
G <- matrix(c(4, -1.2, -1.2, 1.0), 2, 2, dimnames = list(lab, lab))   # genetic covariance (raw units)
P <- matrix(c(6, -0.8, -0.8, 1.6), 2, 2, dimnames = list(lab, lab))   # phenotypic covariance (raw units)

scored <- ng_add_multitrait_score(scores, traits, method = "economic_index",
                                  genetic_covariance = G, phenotypic_covariance = P, out_col = "mt")

# True raw Smith-Hazel: economic weights carry the direction sign (a decrease trait has NEGATIVE
# marginal economic value); b = P^{-1} G a; index = predicted trait matrix %*% b.
a_raw     <- c(2, -3)                                        # yield increase (+), disease decrease (-)
b_raw     <- as.numeric(solve(P, G %*% a_raw))
raw_index <- as.numeric(as.matrix(scores[, c("pred_yield", "pred_disease")]) %*% b_raw)

rho <- cor(scored$mt, raw_index, method = "spearman")
cat(sprintf("  spearman(package economic_index, raw Smith-Hazel) = %.4f\n", rho))
stopifnot(rho > 0.999)

# Sanity on orientation: the top-ranked cross should beat the bottom on the economically-weighted,
# direction-oriented merit (high yield and/or low disease).
top <- scores[which.max(scored$mt), ]; bot <- scores[which.min(scored$mt), ]
merit <- function(r) 2 * r$pred_yield - 3 * r$pred_disease      # oriented economic merit
stopifnot(merit(top) > merit(bot))

# Regression: economic_index WITHOUT an external covariance still runs and ranks (internal path).
scored0 <- ng_add_multitrait_score(scores, traits, method = "economic_index", out_col = "mt")
stopifnot(is.finite(scored0$mt), length(unique(rank(scored0$mt))) > 1L)

cat("multitrait_economic_index_external_cov.R: PASS\n")
