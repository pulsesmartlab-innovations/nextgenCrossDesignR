# A supplied column can be a measured trait or an index the breeder computed elsewhere,
# and the engine had no way to tell.
#
# prediction_mode = "index_as_trait" looks like the distinction but is not: it renames
# the trait to "selection_index" and does nothing else. It appears in R/39 and the
# capability registry and NOWHERE else in R/ -- it reaches no scoring decision, no
# validation and no output field, and settings does not even report it. And it is
# run-level, so a multi-trait run cannot say "these three are traits and that one is an
# index I already built".
#
# WHY IT MATTERS, and why a label alone is not enough: an externally computed index has
# an ARBITRARY scale and origin. It is a linear combination with weights the breeder
# chose, so "12.4" carries no units anyone can check, and a threshold, a check value or
# a covariance entry stated against it means something different from the same number
# against a real trait. Declaring a column an index makes the engine standardise it, so
# everything downstream is in SD units of that index and says so.
#
# Per-trait, on the direction table -- not a run mode -- so real traits and a supplied
# index can appear in the same run.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(42); n <- 40L; m <- 40L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
gv <- as.numeric(g %*% c(rnorm(8, 0, 1), rep(0, m - 8)))
# A deliberately awkward scale: a breeder's index on an arbitrary origin and spread.
idx <- 500 + 37 * as.numeric(scale(gv + rnorm(n, 0, 0.2 * stats::sd(gv))))
pheno <- data.frame(NAME = ids, my_index = idx, stringsAsFactors = FALSE)
mm <- data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                 bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6)
run <- function(dirtab, ...) ng_run_cross_prediction(
  phenotype = pheno, genotype = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = mm, trait_direction = dirtab,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 8L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L, ...)

d_trait <- data.frame(trait = "my_index", column = "my_index", direction = "increase",
                      stringsAsFactors = FALSE)
d_index <- cbind(d_trait, value_kind = "index", stringsAsFactors = FALSE)

# ---- 1. default is "trait": nothing changes for existing configs ----------
r_t <- run(d_trait)
stopifnot(identical(as.character(r_t$effect_summary$value_kind[[1L]]), "trait"))
# a measured trait is left on its own scale
stopifnot(abs(mean(r_t$candidate_crosses$my_index_mean) - mean(idx)) < 5)

# ---- 2. declaring it an index applies the inverse normal transform -------
r_i <- run(d_index)
stopifnot(identical(as.character(r_i$effect_summary$value_kind[[1L]]), "index"))
stopifnot(identical(as.character(r_i$effect_summary$index_transform[[1L]]),
                    "inverse_normal_rank"))
# the reported values are now normal scores, not the breeder's arbitrary scale
stopifnot(abs(mean(r_i$candidate_crosses$my_index_mean)) < 0.5)
stopifnot(abs(mean(r_t$candidate_crosses$my_index_mean)) > 100)   # the untouched scale

# ---- 3. the transform itself: rank-based, clamped, monotone --------------
# Asserted directly on the function, because these are the properties that make it the
# right choice for a quantity whose units nobody can check.
set.seed(1); xx <- c(rlnorm(200), NA, NA)          # heavily skewed, with gaps
tt <- ng_inverse_normal_transform(xx)
stopifnot(length(tt) == length(xx))
stopifnot(all(is.na(tt[is.na(xx)])))                        # NA in, NA out
stopifnot(all(is.finite(tt[!is.na(xx)])))                   # eps clamp: never +/-Inf
stopifnot(identical(order(tt[!is.na(xx)]), order(xx[!is.na(xx)])))   # MONOTONE
fin <- tt[!is.na(tt)]
stopifnot(abs(mean(fin)) < 0.05, abs(stats::sd(fin) - 1) < 0.1)     # ~ standard normal
# invariant to ANY monotone re-expression -- the property centre-and-scale does not have
stopifnot(isTRUE(all.equal(tt, ng_inverse_normal_transform(exp(xx / 3)))))
stopifnot(isTRUE(all.equal(tt, ng_inverse_normal_transform(1e6 * xx + 500))))
# ties share a value (average ranks), rather than being broken by input order
tie <- ng_inverse_normal_transform(c(5, 5, 1, 9))
stopifnot(isTRUE(all.equal(tie[[1L]], tie[[2L]])))

# ---- 4. declaring an index CHANGES the selection, and that is the point ---
# The inverse normal transform is monotone but NOT affine, so unlike centre-and-scale it
# does not leave the ranking alone: the mid-parent of transformed values is not the
# transform of the mid-parent, and the within-family variance moves non-proportionally.
#
# This is a real consequence a breeder must be told, not a side effect to hide. Ranking
# on a raw supplied index means ranking on a scale whose spacing came from weights
# nobody can check; the transform replaces that with a defensible one. Asserting the
# difference here keeps it an explicit, tested decision.
sel_t <- paste(r_t$selected_crosses$parent1, r_t$selected_crosses$parent2)
sel_i <- paste(r_i$selected_crosses$parent1, r_i$selected_crosses$parent2)
stopifnot(!identical(sel_t, sel_i))
# ...and the transformed run is still a well formed plan, not merely a different one
stopifnot(nrow(r_i$selected_crosses) == nrow(r_t$selected_crosses))
stopifnot(all(is.finite(r_i$selected_crosses$my_index_value)))
stopifnot(all(sel_i %in% paste(r_i$candidate_crosses$parent1,
                               r_i$candidate_crosses$parent2)))

# PARENT-level order is preserved, which is what "monotone" buys: the transform never
# reorders the lines themselves, only the family-level arithmetic built on them.
stopifnot(identical(order(ng_inverse_normal_transform(idx)), order(idx)))

# ---- 5. an unknown value_kind is refused at entry, before any fitting ----
bad <- cbind(d_trait, value_kind = "composite", stringsAsFactors = FALSE)
e <- tryCatch({ run(bad); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(e), grepl("value_kind", e, fixed = TRUE))
stopifnot(grepl("trait", e, fixed = TRUE), grepl("index", e, fixed = TRUE))

# ---- 6. traits and a supplied index can coexist in one run ---------------
# This is the case a run-level mode cannot express at all.
pheno2 <- data.frame(NAME = ids, yield = gv + rnorm(n, 0, 0.2 * stats::sd(gv)),
                     my_index = idx, stringsAsFactors = FALSE)
d_mix <- data.frame(trait = c("yield", "my_index"), column = c("yield", "my_index"),
                    direction = c("increase", "increase"),
                    value_kind = c("trait", "index"), stringsAsFactors = FALSE)
r_m <- ng_run_cross_prediction(
  phenotype = pheno2, genotype = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = mm, trait_direction = d_mix,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 8L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L)
es <- r_m$effect_summary
stopifnot(identical(as.character(es$value_kind[match("yield", es$trait)]), "trait"))
stopifnot(identical(as.character(es$value_kind[match("my_index", es$trait)]), "index"))
stopifnot(is.na(es$index_transform[match("yield", es$trait)]))          # untouched
stopifnot(identical(as.character(es$index_transform[match("my_index", es$trait)]),
                    "inverse_normal_rank"))

# ---- 7. DIRECTION is honoured for an index, because it is not always "higher" -----
# A supplied index may itself be a rank summation index, where LOWER is better, or a
# Smith-Hazel index, where higher is. The engine must not assume. Standardisation uses
# a positive scale (an SD), so it preserves order and cannot silently flip this.
d_idx_dn <- data.frame(trait = "my_index", column = "my_index", direction = "decrease",
                       value_kind = "index", stringsAsFactors = FALSE)
r_dn <- run(d_idx_dn)
r_up <- run(d_index)
stopifnot(identical(as.character(r_dn$effect_summary$value_kind[[1L]]), "index"))
# same transform either way -- direction is orientation, not a rescaling
stopifnot(identical(as.character(r_dn$effect_summary$index_transform[[1L]]),
                    as.character(r_up$effect_summary$index_transform[[1L]])))
# but the merit points the other way: with "decrease", usefulness sits BELOW the mean
i10 <- ng_selection_intensity(r_dn$settings$selection_prop)
cu <- r_up$candidate_crosses; cd <- r_dn$candidate_crosses
stopifnot(all(cu$my_index_value >= cu$my_index_mean - 1e-9))   # increase: mean + i*sd
stopifnot(all(cd$my_index_value <= cd$my_index_mean + 1e-9))   # decrease: mean - i*sd
# and the two runs genuinely choose different crosses
stopifnot(!identical(paste(r_up$selected_crosses$parent1, r_up$selected_crosses$parent2),
                     paste(r_dn$selected_crosses$parent1, r_dn$selected_crosses$parent2)))

# the transform never flips direction: a monotone map cannot invert "higher is better"
stopifnot(cor(cu$my_index_mean, r_t$candidate_crosses$my_index_mean,
              method = "spearman") > 0.95)

cat("declared_index_is_normalized: PASS\n")
