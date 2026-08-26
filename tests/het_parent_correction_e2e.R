# End-to-end: with parent_type='ril' AND phased_haplotypes, het-parent crosses
# get the exact residual-het variance (ng_gms_additive_var_general). A run
# without phase is blocked rather than returning the biased-low inbred a'Ra.
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

missing_phase <- tryCatch(
  do.call(ng_run_cross_prediction, args0),
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("phased_haplotypes", missing_phase, fixed = TRUE))
res_corr <- suppressWarnings(do.call(ng_run_cross_prediction, c(args0, list(phased_haplotypes = ph)))) # phase
cc <- res_corr$candidate_crosses
vcol <- grep("_vpm$", names(cc), value = TRUE)[1]
stopifnot(!is.na(vcol))
is_p1 <- cc$parent1 == "P01" | cc$parent2 == "P01"
stopifnot(any(is_p1), all(is.finite(cc[[vcol]])))

# Independently reconstruct the exact phased variance and the forbidden
# inbred-parent shortcut for every P01 cross. This proves the runner wires the
# exact kernel and that the missing-phase guard prevents a material downward bias.
markers <- colnames(res_corr$cleaned_data$genotype)
beta <- res_corr$marker_effects$yield$beta[markers]
R <- ng_recomb_decay_matrix(res_corr$cleaned_data$marker_map, target = "RIL")
hap_rows <- matrix(
  0, nrow = 2L * length(ids), ncol = length(markers),
  dimnames = list(as.vector(rbind(paste0(ids, "_HapA"), paste0(ids, "_HapB"))), markers)
)
for (k in seq_along(ids)) {
  hap_rows[2L * k - 1L, ] <- ph[ids[[k]], 1L, markers]
  hap_rows[2L * k, ] <- ph[ids[[k]], 2L, markers]
}
exact <- biased <- rep(NA_real_, nrow(cc))
for (i in which(is_p1)) {
  exact[[i]] <- ng_gms_additive_var_general(
    cc$parent1[[i]], cc$parent2[[i]], hap_rows, beta, R, target = "RIL"
  )[["VPM"]]
  i1 <- match(cc$parent1[[i]], ids); i2 <- match(cc$parent2[[i]], ids)
  a <- 0.5 * (geno[i1, markers] - geno[i2, markers]) * beta
  biased[[i]] <- as.numeric(crossprod(a, R %*% a))
}
stopifnot(max(abs(cc[[vcol]][is_p1] - exact[is_p1])) < 1e-8)
stopifnot(all(exact[is_p1] >= biased[is_p1] - 1e-9))
stopifnot(mean(exact[is_p1]) > mean(biased[is_p1]))

cat(sprintf("het_parent_correction_e2e: missing phase blocked; exact het-parent VPM exceeds biased a'Ra (%.4f -> %.4f)\n",
            mean(biased[is_p1]), mean(exact[is_p1])))
