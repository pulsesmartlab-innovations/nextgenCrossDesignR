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

# ---- 2. declaring it an index standardises it, and records the transform --
r_i <- run(d_index)
stopifnot(identical(as.character(r_i$effect_summary$value_kind[[1L]]), "index"))
ctr <- r_i$effect_summary$index_center[[1L]]
scl <- r_i$effect_summary$index_scale[[1L]]
stopifnot(is.finite(ctr), is.finite(scl), scl > 0)
stopifnot(abs(ctr - mean(idx)) < 1e-8)
stopifnot(abs(scl - stats::sd(idx)) < 1e-8)
# the cross means are now in SD units of the index, centred near zero
stopifnot(abs(mean(r_i$candidate_crosses$my_index_mean)) < 0.5)
stopifnot(abs(mean(r_t$candidate_crosses$my_index_mean)) > 100)   # the untouched scale

# ---- 3. the transform is recoverable: it is a stated affine map ----------
# A breeder must be able to get back to their own units from the reported numbers.
back <- r_i$candidate_crosses$my_index_mean * scl + ctr
stopifnot(max(abs(sort(back) - sort(r_t$candidate_crosses$my_index_mean))) < 1e-6)

# ---- 4. normalising does NOT perturb the ranking -------------------------
# Standardisation is affine and the merit is affine in the mean, so a single-trait
# ordering is invariant. Worth asserting: it tells a breeder the declaration changes
# the REPORTED UNITS, not which crosses are chosen.
stopifnot(identical(paste(r_t$selected_crosses$parent1, r_t$selected_crosses$parent2),
                    paste(r_i$selected_crosses$parent1, r_i$selected_crosses$parent2)))

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
stopifnot(!is.finite(es$index_scale[match("yield", es$trait)]))   # untouched
stopifnot(is.finite(es$index_scale[match("my_index", es$trait)]))

# ---- 7. DIRECTION is honoured for an index, because it is not always "higher" -----
# A supplied index may itself be a rank summation index, where LOWER is better, or a
# Smith-Hazel index, where higher is. The engine must not assume. Standardisation uses
# a positive scale (an SD), so it preserves order and cannot silently flip this.
d_idx_dn <- data.frame(trait = "my_index", column = "my_index", direction = "decrease",
                       value_kind = "index", stringsAsFactors = FALSE)
r_dn <- run(d_idx_dn)
r_up <- run(d_index)
stopifnot(identical(as.character(r_dn$effect_summary$value_kind[[1L]]), "index"))
# same normalisation either way -- direction is orientation, not scale
stopifnot(isTRUE(all.equal(r_dn$effect_summary$index_scale[[1L]],
                           r_up$effect_summary$index_scale[[1L]])))
# but the merit points the other way: with "decrease", usefulness sits BELOW the mean
i10 <- ng_selection_intensity(r_dn$settings$selection_prop)
cu <- r_up$candidate_crosses; cd <- r_dn$candidate_crosses
stopifnot(all(cu$my_index_value >= cu$my_index_mean - 1e-9))   # increase: mean + i*sd
stopifnot(all(cd$my_index_value <= cd$my_index_mean + 1e-9))   # decrease: mean - i*sd
# and the two runs genuinely choose different crosses
stopifnot(!identical(paste(r_up$selected_crosses$parent1, r_up$selected_crosses$parent2),
                     paste(r_dn$selected_crosses$parent1, r_dn$selected_crosses$parent2)))

# normalisation never flips sign: the standardised column must preserve the raw order
stopifnot(isTRUE(all.equal(order(cu$my_index_mean),
                           order(r_t$candidate_crosses$my_index_mean))))

cat("declared_index_is_normalized: PASS\n")
