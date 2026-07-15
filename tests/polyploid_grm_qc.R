# Correct polyploid GRM (R/47) + ploidy-aware QC. Validates that ng_polyploid_grm is the
# allele-frequency-based VanRaden matrix generalized to ploidy (NOT the ploidy-midpoint shortcut),
# and that ng_polyploid_qc flags/cleans dosage data against 0..ploidy.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(47)

# --- 1. ploidy=2 GRM matches the textbook VanRaden formula G = WW'/(2 sum p(1-p)), W = M-2p ---
np <- 20L; nm <- 60L
M2 <- matrix(rbinom(np * nm, 2L, runif(nm, 0.15, 0.85)[col(matrix(0, np, nm))]), np, nm,
             dimnames = list(sprintf("S%02d", 1:np), sprintf("m%02d", 1:nm)))
p <- colMeans(M2) / 2
poly <- p > 0 & p < 1
Mk <- M2[, poly, drop = FALSE]; pk <- p[poly]
W <- sweep(Mk, 2, 2 * pk); G_manual <- tcrossprod(W) / (2 * sum(pk * (1 - pk)))
G <- ng_polyploid_grm(M2, ploidy = 2L)
stopifnot(max(abs(G - G_manual)) < 1e-8)   # exact VanRaden reduction at ploidy 2

# --- 2. it is NOT the ploidy-midpoint shortcut when allele freqs != 0.5 ---
G_old <- ng_poly4x_parent_relationship(M2, ploidy = 2L)
stopifnot(max(abs(G - G_old)) > 1e-6)      # genuinely different from the freq-0.5 version

# --- 3. any ploidy: correct scaling, sane matrix (symmetric, diagonal ~1 on average) ---
for (ploidy in c(4L, 6L)) {
  Mp <- matrix(rbinom(np * nm, ploidy, 0.5), np, nm, dimnames = list(rownames(M2), colnames(M2)))
  Gp <- ng_polyploid_grm(Mp, ploidy = ploidy)
  stopifnot(isSymmetric(unname(Gp)), identical(attr(Gp, "ploidy"), ploidy))
  stopifnot(abs(mean(diag(Gp)) - 1) < 0.25)          # relationship scale (~1 + F)
  stopifnot(abs(mean(Gp)) < 0.2)                      # centered: mean relationship ~ 0
}

# --- 4. missing-data imputation to expected dosage, and min_maf / monomorphic dropping ---
Mmiss <- M2; Mmiss[sample(length(Mmiss), 30)] <- NA
Gm <- ng_polyploid_grm(Mmiss, ploidy = 2L, impute_missing = TRUE)
stopifnot(all(is.finite(Gm)))
Mmono <- cbind(M2, mono = rep(2L, np))                # a monomorphic marker
Gmono <- ng_polyploid_grm(Mmono, ploidy = 2L, return_freq = TRUE)
stopifnot(attr(Gmono, "n_markers") == sum(poly))      # monomorphic marker excluded

# --- 4b. multi-method: Yang (GCTA) additive GRM is a valid alternative (diagonal ~1, centered) ---
Gy <- ng_polyploid_grm(M2, ploidy = 2L, method = "yang")
stopifnot(isSymmetric(unname(Gy)), abs(mean(diag(Gy)) - 1) < 0.3, abs(mean(Gy)) < 0.2)
stopifnot(identical(attr(Gy, "method"), "yang_polyploid"))
stopifnot(max(abs(Gy - G)) > 1e-8)          # Yang != VanRaden (different weighting)

# --- 4c. DOMINANCE GRM (digenic): symmetric, distinct from additive, both methods run ---
for (ploidy in c(2L, 4L)) {
  Mp <- matrix(rbinom(np * nm, ploidy, 0.5), np, nm, dimnames = list(rownames(M2), colnames(M2)))
  Gd <- ng_polyploid_dominance_grm(Mp, ploidy = ploidy, method = "vanraden")
  Ga <- ng_polyploid_grm(Mp, ploidy = ploidy, method = "vanraden")
  stopifnot(isSymmetric(unname(Gd)), identical(attr(Gd, "component"), "dominance"))
  stopifnot(abs(mean(diag(Gd)) - 1) < 0.3)         # scaled to relationship (~1)
  stopifnot(max(abs(Gd - Ga)) > 1e-6)              # dominance != additive relationship
  Gdy <- ng_polyploid_dominance_grm(Mp, ploidy = ploidy, method = "yang")
  stopifnot(isSymmetric(unname(Gdy)), identical(attr(Gdy, "method"), "yang_dominance_polyploid"))
}

# --- 5. QC: flags out-of-range, missingness, monomorphic, duplicates; returns a clean matrix ---
Mqc <- matrix(rbinom(np * nm, 4L, 0.5), np, nm, dimnames = list(sprintf("S%02d", 1:np), sprintf("m%02d", 1:nm)))
Mqc[1, 1] <- 9          # out of range for ploidy 4
Mqc[, 2] <- 4L          # monomorphic
Mqc[3, ] <- NA          # a mostly-missing sample
qc <- ng_polyploid_qc(Mqc, ploidy = 4L, max_missing_sample = 0.5, drop_monomorphic = TRUE)
stopifnot(qc$summary$n_out_of_range == 1L)
stopifnot(qc$marker_report$dropped[2])                       # monomorphic marker flagged
stopifnot(qc$sample_report$dropped[3])                       # missing sample flagged
stopifnot(ncol(qc$clean) < nm, nrow(qc$clean) < np)          # cleaned matrix is smaller
stopifnot(!qc$pass)                                          # out-of-range => pass is FALSE

# --- 6. GRM of the cleaned matrix is usable ---
Gclean <- ng_polyploid_grm(qc$clean, ploidy = 4L)
stopifnot(all(is.finite(Gclean)), nrow(Gclean) == nrow(qc$clean))

cat("polyploid GRM + QC test passed\n")
