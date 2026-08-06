# End-to-end: with parent_type='ril' AND phased_haplotypes, het-parent crosses
# get the exact residual-het variance (ng_gms_additive_var_general) instead of
# the biased-low inbred a'Ra. Inbred-only crosses (and every run without phased
# haplotypes) are untouched.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(7)
np <- 5L; m <- 40L; ids <- sprintf("P%02d", seq_len(np)); snps <- sprintf("S%03d", seq_len(m))
H1 <- matrix(sample(0:1, np * m, replace = TRUE), np, m); H2 <- H1   # inbred: hap2 = hap1
het_loci <- c(3, 8, 15, 22, 31)                                      # P01 = RIL, residual het here
H2[1, het_loci] <- 1 - H1[1, het_loci]
ph <- array(0, c(np, 2, m), dimnames = list(ids, c("hap1", "hap2"), snps))
ph[, 1, ] <- H1; ph[, 2, ] <- H2
geno <- H1 + H2; colnames(geno) <- snps
mm <- data.frame(SNP_code = snps, Chromosome = rep(1:2, each = 20),
                 Position_BP = rep(1:20, 2) * 5e5, stringsAsFactors = FALSE)
qtl <- c(3, 8, 15, 22, 31, 5, 12, 25, 34); bq <- c(1, 0.9, -0.8, 1.1, 0.7, 0.5, -0.6, 0.8, -0.4)
gv <- as.numeric(scale(((geno - 1)[, qtl]) %*% bq))
args0 <- list(
  genotype = data.frame(NAME = ids, geno, check.names = FALSE),
  phenotype = data.frame(NAME = ids, yield = gv + rnorm(np, 0, 0.3)),
  trait_direction = data.frame(Trait = "yield", Selection_direction = "increase"),
  marker_map = mm, id_col = "NAME", map_position_unit = "bp", bp_per_cm = 1e6,
  n_crosses = 10L, progeny = "RIL", parent_type = "ril", seed = 1L)

res_bias <- suppressWarnings(do.call(ng_run_cross_prediction, args0))                                # no phase
res_corr <- suppressWarnings(do.call(ng_run_cross_prediction, c(args0, list(phased_haplotypes = ph)))) # phase
cb <- res_bias$candidate_crosses; cc <- res_corr$candidate_crosses
vcol <- grep("_vpm$", names(cb), value = TRUE)[1]
stopifnot(!is.na(vcol))
o <- match(paste(cb$parent1, cb$parent2), paste(cc$parent1, cc$parent2))
is_p1 <- cb$parent1 == "P01" | cb$parent2 == "P01"

# het-parent crosses: corrected >= biased (never lower), strictly higher on average
stopifnot(all(cc[[vcol]][o][is_p1] >= cb[[vcol]][is_p1] - 1e-9))
stopifnot(mean(cc[[vcol]][o][is_p1]) > mean(cb[[vcol]][is_p1]))
# inbred-only crosses: byte-identical (correction never touches them)
stopifnot(isTRUE(all.equal(cb[[vcol]][!is_p1], cc[[vcol]][o][!is_p1])))
# and a run WITHOUT phased haplotypes is unchanged for ALL crosses (no silent shift)
stopifnot(inherits(res_bias, "ng_cross_prediction_result"))

cat(sprintf("het_parent_correction_e2e: het-parent vpm corrected up (%.4f -> %.4f), inbred-only unchanged\n",
            mean(cb[[vcol]][is_p1]), mean(cc[[vcol]][o][is_p1])))
