# The rank summation index (Mulamba & Mock 1978), weighted.
#
# The package had five multi-trait methods and none of them was this one. The confusion
# is understandable: `auto` resolves to a family called "rank_threshold" and runs through
# ng_rank_normalize(), so the word "rank" is everywhere. But that function does
#
#     r <- rank(value); p <- (r - 0.5)/n; z <- qnorm(p); z <- standardize(z)
#
# which is a van der Waerden / Blom NORMAL-SCORE transform. It is rank-BASED, and it is
# not a rank summation index. The qnorm step deliberately re-imposes Gaussian spacing --
# the gap between the top two entries stretches, the gaps in the middle compress --
# which is precisely the spacing an RSI exists to discard.
#
# THE PROPERTY THAT MAKES IT WORTH HAVING, asserted below: under RSI a trait's extreme
# outlier contributes exactly one rank step, no matter how extreme. Under a normal-score
# index it contributes a stretched tail value. A breeder choosing RSI is choosing that
# insensitivity; if the implementation did not have it, the method would be RSI in name
# only.
#
# Weighted, with equal weights the default -- and with equal weights sum(1 * rank_i) is
# exactly the classic unweighted index, so Mulamba-Mock is the default case rather than
# a separate mode.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---- 1. the method exists, validates, and has its own family --------------
stopifnot("rank_sum" %in% ng_multitrait_methods())
stopifnot(identical(ng_multitrait_validate_method("rank_sum"), "rank_sum"))
stopifnot(identical(ng_multitrait_method_family("rank_sum"), "rank_sum"))
# a caller who types the literature name is told the spelling, not just rejected
e <- tryCatch(ng_multitrait_validate_method("mulamba_mock"), error = function(e) conditionMessage(e))
stopifnot(is.character(e), grepl("rank_sum", e, fixed = TRUE))

# ---- 2. equal weights reproduce the classic unweighted rank sum -----------
n <- 12L
sc <- data.frame(
  parent1 = sprintf("P%02d", seq_len(n)), parent2 = sprintf("Q%02d", seq_len(n)),
  yield   = c(10, 14, 12, 18, 11, 16, 13, 15, 17, 9, 19, 20),
  disease = c(5, 2, 7, 3, 6, 1, 8, 4, 2.5, 9, 3.5, 6.5),
  stringsAsFactors = FALSE)
td <- data.frame(trait = c("yield","disease"), column = c("yield","disease"),
                 direction = c("maximize","minimize"), weight = c(1, 1),
                 stringsAsFactors = FALSE)
out <- ng_add_multitrait_score(sc, td, method = "rank_sum")
got <- out$multi_trait_score

# hand-rolled Mulamba-Mock: rank each trait with 1 = BEST, sum, smallest wins.
r_yield   <- rank(-sc$yield,   ties.method = "average")   # maximize -> negate
r_disease <- rank( sc$disease, ties.method = "average")   # minimize -> as-is
rsi <- r_yield + r_disease
# the package reports higher = better (every other family does), so the ORDER must be
# the exact reverse of the classic smallest-sum-wins ordering.
stopifnot(identical(order(got, decreasing = TRUE), order(rsi, decreasing = FALSE)))
stopifnot(isTRUE(all.equal(cor(got, rsi, method = "spearman"), -1)))

# ---- 3. THE ROBUSTNESS PROPERTY: an outlier moves the score one rank step --
# Make the best yield astronomically better. Under RSI nothing changes: it was already
# rank 1 and still is. Under a normal-score index the top z stretches and the whole
# index shifts.
sc_out <- sc; sc_out$yield[12] <- 1e6
got_out <- ng_add_multitrait_score(sc_out, td, method = "rank_sum")$multi_trait_score
stopifnot(isTRUE(all.equal(got, got_out)))          # RSI: entirely unmoved

z_before <- ng_add_multitrait_score(sc,     td, method = "weighted")$multi_trait_score
z_after  <- ng_add_multitrait_score(sc_out, td, method = "weighted")$multi_trait_score
# Not asserting the normal-score index is WRONG -- only that the two genuinely differ,
# which is why having both is worth something.
stopifnot(!isTRUE(all.equal(got, z_before)))

# ---- 4. weights do what weights should ------------------------------------
td0 <- td; td0$weight <- c(1, 0)          # disease switched off entirely
g0 <- ng_add_multitrait_score(sc, td0, method = "rank_sum")$multi_trait_score
stopifnot(isTRUE(all.equal(cor(g0, r_yield, method = "spearman"), -1)))   # yield alone

td2 <- td; td2$weight <- c(3, 1)          # yield dominant
g2 <- ng_add_multitrait_score(sc, td2, method = "rank_sum")$multi_trait_score
stopifnot(cor(g2, -r_yield, method="spearman") > cor(got, -r_yield, method="spearman"))

# ---- 5. direction is honoured, not assumed --------------------------------
td_flip <- td; td_flip$direction <- c("minimize", "minimize")
gf <- ng_add_multitrait_score(sc, td_flip, method = "rank_sum")$multi_trait_score
stopifnot(cor(gf, sc$yield, method = "spearman") < 0)   # now LOW yield is good

# ---- 6. end to end through the runner -------------------------------------
set.seed(9); N <- 40L; m <- 40L; ids <- sprintf("P%02d", seq_len(N))
g <- matrix(2L * rbinom(N * m, 1, 0.5), N, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
y1 <- as.numeric(g %*% c(rnorm(8,0,1), rep(0, m-8)))
y2 <- as.numeric(g %*% c(rep(0,8), rnorm(8,0,1), rep(0, m-16)))
r <- ng_run_cross_prediction(
  phenotype = data.frame(NAME = ids, yield = y1 + rnorm(N,0,0.2*sd(y1)),
                                     protein = y2 + rnorm(N,0,0.2*sd(y2))),
  genotype = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                          bp = rep(seq(0,100,length.out=m/2),2)[seq_len(m)]*1e6),
  trait_direction = data.frame(trait = c("yield","protein"), column = c("yield","protein"),
                               direction = c("increase","increase")),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 8L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L, multi_trait_method = "rank_sum")
stopifnot(inherits(r, "ng_cross_prediction_result"))
stopifnot(nrow(r$selected_crosses) == 8L)
mt <- attr(r$candidate_crosses, "multi_trait")
stopifnot(identical(mt$family, "rank_sum"))
stopifnot(all(is.finite(r$candidate_crosses$multi_trait_score)))
# the run REPORTS the method it used, not the one requested-but-ignored
stopifnot(identical(r$settings$multi_trait_method, "rank_sum"))

# ---- 7. it is declarable, so a frontend can offer it ----------------------
ctls <- ng_backend_controls()
mtm <- Filter(function(c) identical(c$id, "multi_trait_method"), ctls)[[1L]]
stopifnot("rank_sum" %in% vapply(mtm$choices, function(x) x$value, character(1)))

cat("rank_summation_index: PASS\n")
