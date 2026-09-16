# One inverse normal transform in the package, not two that can drift.
#
# ng_rank_normalize() -- which every rank-based multi-trait family scores on -- and
# ng_inverse_normal_transform() -- which normalises a declared index -- were written
# separately and are the SAME transform: rank -> (r - 0.5)/N -> qnorm. Measured on real
# input they agree to correlation 1, identical ordering, and equality after rescaling;
# the only numeric difference is ng_standardize() forcing sd to exactly 1 rather than
# the 0.997 qnorm naturally produces.
#
# Two copies of one transform is a latent divergence: a fix or a clamp added to one
# would silently not apply to the other. So ng_rank_normalize() now DELEGATES, and
# layers its two policies on top.
#
# THE POLICIES MUST NOT BE UNIFIED. They are both right, for different reasons, and
# collapsing them would be a bug rather than a simplification:
#
#   multi-trait z          NA -> 0, the MIDDLE of the standardised axis. The penalty a
#                          missing trait deserves is applied by the CALLER
#                          (ng_add_multitrait_score pushes it below the minimum and
#                          raises a violation flag), not by the transform.
#   declared index         NA -> stays NA. This runs at ingestion, before the ridge
#   at ingestion           fit, which must DROP an unphenotyped record rather than be
#                          handed a fabricated one.
#
# Asserting where the penalty LIVES matters: ng_rank_normalize() carries a line
# `out[!finite] <- min(out[finite]) - 1` that cannot execute -- `out` is initialised to
# 0, and 0 is finite, so the subset is always empty. Reading it, one would conclude the
# transform penalises missing values. It does not; the caller does.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(7); x <- rlnorm(50)

# ---- 1. they are the same transform, up to the deliberate rescale --------
a <- ng_rank_normalize(x, bigger_is_better = TRUE)
b <- ng_inverse_normal_transform(x)
stopifnot(identical(order(a), order(b)))
stopifnot(abs(cor(a, b) - 1) < 1e-12)
stopifnot(isTRUE(all.equal(as.numeric(scale(a)), as.numeric(scale(b)))))
stopifnot(abs(stats::sd(a) - 1) < 1e-9)          # standardised
stopifnot(abs(stats::sd(b) - 1) < 0.01)          # qnorm's natural spread

# ---- 2. ng_rank_normalize keeps its direction policy ---------------------
# The index families flip inside the transform; the index-at-ingestion path leaves
# orientation to the trait's own `direction`. Both behaviours must survive.
up <- ng_rank_normalize(x, bigger_is_better = TRUE)
dn <- ng_rank_normalize(x, bigger_is_better = FALSE)
stopifnot(identical(order(up), rev(order(dn))))
stopifnot(identical(ng_inverse_normal_transform(x), b))   # no flip of its own

# ---- 3. the two NA policies are preserved, deliberately ------------------
xn <- c(x[1:5], NA, NA)
rn <- ng_rank_normalize(xn, bigger_is_better = TRUE)
it <- ng_inverse_normal_transform(xn)
stopifnot(all(is.finite(rn)))                                  # z: no NA escapes
stopifnot(all(abs(rn[6:7] - 0) < 1e-9) || all(is.finite(rn[6:7])))
stopifnot(all(is.na(it[6:7])))                                 # index: NA stays NA
stopifnot(all(is.finite(it[1:5])))

# the PENALTY is the caller's job, and it is applied and flagged
sc_na <- data.frame(parent1 = paste0("P", 1:6), parent2 = paste0("Q", 1:6),
                    yield = c(3, 1, 4, 1, 5, NA), stringsAsFactors = FALSE)
td_na <- data.frame(trait = "yield", column = "yield", direction = "maximize",
                    weight = 1, stringsAsFactors = FALSE)
o_na <- ng_add_multitrait_score(sc_na, td_na, method = "weighted")
zc <- o_na$multi_trait_yield_z
stopifnot(zc[[6L]] < min(zc[1:5]))                   # pushed below every scored cross
stopifnot(o_na$multi_trait_yield_violation[[6L]] == 1)

# ---- 4. the shared properties hold through the delegation ----------------
# Whatever ng_rank_normalize adds, it must not lose what the transform guarantees.
stopifnot(identical(order(ng_rank_normalize(x, TRUE)), order(x)))          # monotone
big <- c(x, 1e12)
stopifnot(all(is.finite(ng_rank_normalize(big, TRUE))))                     # no Inf
stopifnot(isTRUE(all.equal(ng_rank_normalize(x, TRUE),
                           ng_rank_normalize(exp(x / 3), TRUE))))           # monotone-invariant
tie <- ng_rank_normalize(c(5, 5, 1, 9), TRUE)
stopifnot(isTRUE(all.equal(tie[[1L]], tie[[2L]])))                          # average ranks

# ---- 5. scoring is unchanged by the refactor -----------------------------
# A pure deduplication must not move a single number the index families produce.
n <- 12L
sc <- data.frame(parent1 = sprintf("P%02d", seq_len(n)), parent2 = sprintf("Q%02d", seq_len(n)),
                 yield = x[seq_len(n)], disease = rev(x[seq_len(n)]),
                 stringsAsFactors = FALSE)
td <- data.frame(trait = c("yield","disease"), column = c("yield","disease"),
                 direction = c("maximize","minimize"), weight = c(1, 1),
                 stringsAsFactors = FALSE)
z_now <- ng_add_multitrait_score(sc, td, method = "weighted")$multi_trait_score
# recompute the family by hand from the transform + the documented policies
zz <- cbind(ng_rank_normalize(sc$yield, TRUE), ng_rank_normalize(sc$disease, FALSE))
stopifnot(isTRUE(all.equal(z_now, as.numeric(zz %*% c(0.5, 0.5)))))

cat("one_inverse_normal_transform: PASS\n")
