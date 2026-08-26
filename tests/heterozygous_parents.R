helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# T1.5: ng_score_crosses() must refuse (or warn) on heterozygous parents because
# the DH/RIL recombination-variance kernel assumes inbred parents (dosage in
# {0, 2}). Before the fix, heterozygous dosages produced silently-wrong PMV.

# ---- ng_audit_inbred_dosage detects residual heterozygosity ------------------------------
geno_inbred <- rbind(
  P1 = c(0, 2, 0, 2, 0, 2, 0, 2, 0, 2),
  P2 = c(2, 0, 2, 0, 2, 0, 2, 0, 2, 0)
)
geno_outbred <- rbind(
  P1 = c(0, 2, 1, 2, 1, 2, 0, 2, 0, 2),   # 2/10 het markers
  P2 = c(2, 1, 2, 0, 2, 0, 2, 1, 2, 0)    # 2/10 het markers
)
colnames(geno_inbred)  <- paste0("M", 1:10)
colnames(geno_outbred) <- paste0("M", 1:10)

audit_in <- ng_audit_inbred_dosage(geno_inbred,  tolerance = 0.05, fraction_tolerance = 0.02)
audit_out <- ng_audit_inbred_dosage(geno_outbred, tolerance = 0.05, fraction_tolerance = 0.02)
stopifnot(length(audit_in$violators) == 0L)
stopifnot(length(audit_out$violators) == 2L)
stopifnot(setequal(audit_out$violators, c("P1", "P2")))

# ---- ng_score_crosses errors when assume_inbred = TRUE (default) -------------------------
mk <- data.frame(marker = paste0("M", 1:10), chr = 1L, pos_cm = seq(0, 90, length.out = 10))
beta <- setNames(rnorm(10, sd = 0.1), paste0("M", 1:10))
effects <- list(beta = beta, beta_var = setNames(rep(0.005, 10), paste0("M", 1:10)),
                reliability = 0.5, intercept = 0, marker_mean = colMeans(geno_outbred))

err <- tryCatch(
  ng_score_crosses(geno = geno_outbred, effects = effects, marker_map = mk,
                   ids = rownames(geno_outbred),
                   adjusted_pheno = setNames(rnorm(2L), rownames(geno_outbred)),
                   selection_prop = 0.5, use_cpp = FALSE),
  error = function(e) e
)
if (!inherits(err, "error")) {
  stop("Expected an error from ng_score_crosses with heterozygous parents under assume_inbred = TRUE")
}
if (!grepl("inbred", conditionMessage(err))) {
  stop("Error message should mention the inbred-parent assumption: ", conditionMessage(err))
}

# ---- RIL residual heterozygosity requires complete phase ---------------------------------
no_phase <- tryCatch(
  ng_score_crosses(geno = geno_outbred, effects = effects, marker_map = mk,
                   ids = rownames(geno_outbred),
                   adjusted_pheno = setNames(rnorm(2L), rownames(geno_outbred)),
                   selection_prop = 0.5, use_cpp = FALSE,
                   parent_type = "ril"),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("phased_haplotypes", no_phase, fixed = TRUE))
ph <- array(
  0, dim = c(nrow(geno_outbred), 2L, ncol(geno_outbred)),
  dimnames = list(rownames(geno_outbred), c("hap1", "hap2"), colnames(geno_outbred))
)
ph[, 1L, ] <- ifelse(geno_outbred == 2, 1, 0)
ph[, 2L, ] <- ifelse(geno_outbred >= 1, 1, 0)
res <- ng_score_crosses(
  geno = geno_outbred, effects = effects, marker_map = mk,
  ids = rownames(geno_outbred),
  adjusted_pheno = setNames(rnorm(2L), rownames(geno_outbred)),
  selection_prop = 0.5, use_cpp = FALSE,
  parent_type = "ril", phased_haplotypes = ph
)
stopifnot(nrow(res) == 1L)
stopifnot(is.finite(res$pmv))

# ---- Pure-inbred input must continue to work without warning -----------------------------
effects2 <- list(beta = beta, beta_var = setNames(rep(0.005, 10), paste0("M", 1:10)),
                 reliability = 0.5, intercept = 0, marker_mean = colMeans(geno_inbred))
res_in <- ng_score_crosses(geno = geno_inbred, effects = effects2, marker_map = mk,
                           ids = rownames(geno_inbred),
                           adjusted_pheno = setNames(rnorm(2L), rownames(geno_inbred)),
                           selection_prop = 0.5, use_cpp = FALSE)
stopifnot(nrow(res_in) == 1L)
stopifnot(is.finite(res_in$pmv))

cat("heterozygous_parents: 4/4 checks passed\n")
cat(sprintf("  max het fraction inbred=%.3f outbred=%.3f\n",
            audit_in$max_fraction, audit_out$max_fraction))
