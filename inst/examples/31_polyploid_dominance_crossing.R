# Dominance-aware mate design for clonal / heterosis crops (cassava, sugarcane, potato), even ploidy.
#
# For these crops a cross's value is the distribution of TOTAL genotypic value among its progeny
# clones -- so both the heterosis-inclusive cross mean AND the within-family variance matter. From
# additive+dominance marker effects, ng_polyploid_score_crosses_dominance() predicts, per cross:
#   cross_mean    = mid-parent breeding value + expected progeny dominance (HETEROSIS)
#   cross_var     = additive segregation + dominance segregation (+ their covariance)
#   cross_usefulness = mean + i*sqrt(var)
# ng_polyploid_design_crosses(dominance = TRUE) does fit -> score -> allocate in one call, and forwards
# every native control (strategy dial, target_coancestry, committed matings, ...). Nonzero
# double_reduction is currently restricted to the autotetraploid single-IBD-pair model.
# Additive-only is the default; dominance requires an explicit experimental opt-in.

library(nextgenCrossDesign)

if (!"dominance" %in% names(formals(ng_polyploid_design_crosses))) {
  stop("This example needs a build with the dominance-aware polyploid design. Reinstall the tarball.",
       call. = FALSE)
}

set.seed(20260703)
ploidy <- 4L
n <- 24L; m <- 160L
ids <- sprintf("cassava%02d", seq_len(n))
freq <- runif(m, 0.2, 0.8)
dosage <- matrix(rbinom(n * m, ploidy, rep(freq, each = n)), n, m,
                 dimnames = list(ids, sprintf("snp%03d", seq_len(m))))
p <- colMeans(dosage) / ploidy
W <- sweep(dosage, 2, ploidy * p); H <- dosage * (ploidy - dosage); Dd <- sweep(H, 2, colMeans(H))
phenotype <- as.numeric(W %*% rnorm(m)) + as.numeric(Dd %*% rnorm(m, 0, 1.2)) + rnorm(n, 0, 3)
names(phenotype) <- ids

# --- lower-level: fit A+D effects, then score crosses (with double reduction) ---
fit <- ng_polyploid_fit_effects(
  dosage, phenotype, ploidy = ploidy, model = "additive_dominance",
  allow_experimental_dominance = TRUE
)
scores <- ng_polyploid_score_crosses_dominance(fit, dosage, selection_prop = 0.10, double_reduction = 0.08)
cat(sprintf("Scored %d crosses. Heterosis range [%.2f, %.2f]; usefulness range [%.2f, %.2f]\n",
            nrow(scores), min(scores$heterosis), max(scores$heterosis),
            min(scores$cross_usefulness), max(scores$cross_usefulness)))
print(utils::head(scores[order(-scores$cross_usefulness),
                         c("parent1", "parent2", "cross_mean", "heterosis", "cross_var", "cross_usefulness")], 4))

# --- one-call: QC -> A+D fit -> dominance cross scoring -> allocation with a strategy dial ---
plan <- ng_polyploid_design_crosses(
  dosage = dosage, n_crosses = 12L, ploidy = ploidy, phenotype = phenotype,
  dominance = TRUE, allow_experimental_dominance = TRUE,
  gain = "usefulness",                            # optimize dominance-aware usefulness
  double_reduction = 0.08,                        # autotetraploid IBD-pair coefficient
  grm_method = "vanraden",                        # or "yang"
  strategy = "balanced",                          # gain-vs-diversity dial (native control)
  max_crosses_per_parent = 4L, method = "greedy_local")

stopifnot(nrow(plan) == 12L, isTRUE(attr(plan, "summary")$dominance),
          all(c("heterosis", "cross_usefulness") %in% names(plan)))
cat(sprintf("\nDesigned 12 crosses (dominance-aware). Mean heterosis of chosen crosses: %.3f\n",
            mean(plan$heterosis)))
message("Polyploid dominance crossing example completed.")
