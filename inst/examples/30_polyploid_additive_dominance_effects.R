# Polyploid additive + dominance genomic prediction: breeding value vs total genotypic value.
#
# WHY DOMINANCE? For CLONAL crops (cassava, potato, sugarcane) you deploy the propagated clone, so
# the relevant merit is TOTAL GENOTYPIC value = breeding value (additive) + dominance deviation.
# ng_polyploid_fit_effects(model = "additive_dominance") estimates both marker-effect components;
# ng_polyploid_predict_value(type = "genotypic" | "breeding") returns either.
#   * select CLONES for release/propagation on "genotypic" value;
#   * rank PARENTS for expected crossing gain on "breeding" value.
# Additive-only (model = "additive") stays the robust default when non-additive variance is small or
# poorly estimated.

library(nextgenCrossDesign)

if (!"ng_polyploid_fit_effects" %in% getNamespaceExports("nextgenCrossDesign")) {
  stop("This example needs a build with ng_polyploid_fit_effects. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ploidy <- 4L
n <- 240L; m <- 150L
ids <- sprintf("clone%03d", seq_len(n))
freq <- runif(m, 0.2, 0.8)
dosage <- matrix(rbinom(n * m, ploidy, rep(freq, each = n)), n, m,
                 dimnames = list(ids, sprintf("snp%03d", seq_len(m))))
# a trait with real additive AND dominance genetic variance
p <- colMeans(dosage) / ploidy
W <- sweep(dosage, 2, ploidy * p)
H <- dosage * (ploidy - dosage); D <- sweep(H, 2, colMeans(H))
true_gv <- as.numeric(W %*% rnorm(m)) + as.numeric(D %*% rnorm(m, 0, 1.2))
names(true_gv) <- ids
phenotype <- true_gv + rnorm(n, 0, sd(true_gv) * 0.5)     # h2 ~ 0.8 on genotypic value
names(phenotype) <- ids

train <- ids[1:180]; test <- ids[181:240]

# --- fit additive-only vs additive + dominance on the training clones ---
fit_a <- ng_polyploid_fit_effects(dosage[train, ], phenotype[train], ploidy = ploidy, model = "additive")
fit_ad <- ng_polyploid_fit_effects(dosage[train, ], phenotype[train], ploidy = ploidy,
                                   model = "additive_dominance")

# --- predict on held-out clones ---
gv_a <- ng_polyploid_predict_value(fit_a, dosage[test, ], type = "genotypic")   # additive-only
gv_ad <- ng_polyploid_predict_value(fit_ad, dosage[test, ], type = "genotypic") # additive + dominance
bv_ad <- ng_polyploid_predict_value(fit_ad, dosage[test, ], type = "breeding")  # breeding value only

cat(sprintf("Held-out genotypic-value accuracy: additive-only %.3f  |  additive+dominance %.3f\n",
            cor(gv_a, true_gv[test]), cor(gv_ad, true_gv[test])))
cat(sprintf("Breeding value != genotypic value (dominance present): cor = %.3f\n", cor(bv_ad, gv_ad)))

# top clones to PROPAGATE are ranked on genotypic value; top PARENTS to cross on breeding value
top_clones <- names(sort(gv_ad, decreasing = TRUE))[1:5]
top_parents <- names(sort(bv_ad, decreasing = TRUE))[1:5]
cat("Top clones to propagate (genotypic):", paste(top_clones, collapse = ", "), "\n")
cat("Top parents to cross (breeding)   :", paste(top_parents, collapse = ", "), "\n")

# (modelling dominance usually raises genotypic-value accuracy when non-additive variance is real
# and well-estimated; the exact margin depends on training size and dominance magnitude.)
stopifnot(is.finite(cor(gv_ad, true_gv[test])), length(top_clones) == 5L)
message("Polyploid additive+dominance effects example completed.")
