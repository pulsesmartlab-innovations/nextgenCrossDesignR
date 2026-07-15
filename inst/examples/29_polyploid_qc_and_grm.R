# Polyploid data QC and genomic relationship matrices (any ploidy).
#
# WHICH GRM METHOD? Both are correct allele-frequency-based estimators generalized to ploidy
# (NOT the old diploid-midpoint shortcut). Choice is a weighting preference, not a performance race:
#   * "vanraden" (default): single overall scaling; weights rare-allele markers more.
#   * "yang" (GCTA): each marker standardized to unit variance; weights markers equally.
# Use ng_polyploid_dominance_grm() for a DOMINANCE relationship (needed for clonal/heterosis crops).

library(nextgenCrossDesign)

if (!"ng_polyploid_grm" %in% getNamespaceExports("nextgenCrossDesign")) {
  stop("This example needs a build with the polyploid GRM/QC functions. Reinstall the current tarball.",
       call. = FALSE)
}

set.seed(20260703)
ploidy <- 4L                                  # autotetraploid; set 2 / 6 for diploid / hexaploid
n <- 60L; m <- 200L
ids <- sprintf("clone%03d", seq_len(n))
freq <- runif(m, 0.1, 0.9)
dosage <- matrix(rbinom(n * m, ploidy, rep(freq, each = n)), n, m,
                 dimnames = list(ids, sprintf("snp%03d", seq_len(m))))
# inject a few QC problems: a monomorphic marker, some missing, a duplicate id
dosage[, 1] <- ploidy
dosage[sample(length(dosage), 200)] <- NA
dosage[2, ] <- dosage[1, ]; rownames(dosage)[2] <- rownames(dosage)[1]

# --- ploidy-aware QC: flags out-of-range/missing/monomorphic/duplicates, returns a clean matrix ---
qc <- ng_polyploid_qc(dosage, ploidy = ploidy, max_missing_marker = 0.20,
                      max_missing_sample = 0.50, min_maf = 0.05, drop_monomorphic = TRUE)
cat("QC summary:\n"); str(qc$summary, give.attr = FALSE)
clean <- qc$clean
stopifnot(qc$summary$n_monomorphic >= 1L, length(qc$duplicate_samples) >= 1L)

# de-duplicate ids for the relationship demo (QC flags duplicates; here we simply drop one)
clean <- clean[!duplicated(rownames(clean)), , drop = FALSE]

# --- additive GRM: VanRaden vs Yang (both diagonal ~1, centered) ---
G_vr <- ng_polyploid_grm(clean, ploidy = ploidy, method = "vanraden")
G_yang <- ng_polyploid_grm(clean, ploidy = ploidy, method = "yang")
cat(sprintf("\nadditive GRM: VanRaden diag ~ %.2f | Yang diag ~ %.2f | markers used %d\n",
            mean(diag(G_vr)), mean(diag(G_yang)), attr(G_vr, "n_markers")))

# --- dominance GRM (digenic) -- distinct from additive; for clonal/heterosis crops ---
D_vr <- ng_polyploid_dominance_grm(clean, ploidy = ploidy, method = "vanraden")
cat(sprintf("dominance GRM diag ~ %.2f | differs from additive: %s\n",
            mean(diag(D_vr)), max(abs(D_vr - G_vr)) > 1e-6))

stopifnot(isSymmetric(unname(G_vr)), isSymmetric(unname(D_vr)),
          identical(attr(D_vr, "component"), "dominance"))
message("Polyploid QC + GRM example completed.")
