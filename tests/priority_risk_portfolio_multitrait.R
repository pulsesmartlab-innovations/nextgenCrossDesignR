ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# Multi-trait cross-priority portfolio + risk: the axes are resolved on the SELECTION INDEX
# (see R/37). Same six columns and the same quadrant semantics as the single-trait path, so
# the reporting surface does not fork -- what differs is how level/upside/pev are derived.

# --- ng_multitrait_index_basis: raw-unit weights, orientation, shared basis ---
M <- cbind(yield = c(10, 12, 14, 16), protein = c(3, 2.5, 2, 1.5))
bs <- ng_multitrait_index_basis(M, c("maximize", "minimize"), c(yield = 0.6, protein = 0.4))
# w_k = coef_k * s_k / scale_k -- the minimize trait must carry a NEGATIVE raw-unit weight,
# so that a higher raw protein value lowers the index level.
stopifnot(bs$w[["yield"]] > 0, bs$w[["protein"]] < 0)
stopifnot(isTRUE(all.equal(unname(bs$w[["yield"]]), 0.6 / bs$scale[[1L]])))
stopifnot(isTRUE(all.equal(unname(bs$w[["protein"]]), -0.4 / bs$scale[[2L]])))
# level must be affine in w'm
lin <- as.numeric(M %*% bs$w)
stopifnot(isTRUE(all.equal(bs$level - mean(bs$level), lin - mean(lin))))
# coefficients are matched BY NAME when the matrix is named, not by position
bs_rev <- ng_multitrait_index_basis(M, c("maximize", "minimize"),
                                    c(protein = 0.4, yield = 0.6))
stopifnot(isTRUE(all.equal(bs$w, bs_rev$w)))
# supplying center/scale reuses the reference basis instead of re-deriving it: a SUBSET must
# land on the same axes as the full table (this is what keeps the selected plan comparable
# to the candidate pool).
sub <- M[2:3, , drop = FALSE]
bs_sub_own <- ng_multitrait_index_basis(sub, c("maximize", "minimize"), c(0.6, 0.4))
bs_sub_ref <- ng_multitrait_index_basis(sub, c("maximize", "minimize"), c(0.6, 0.4),
                                        center = bs$center, scale = bs$scale)
stopifnot(isTRUE(all.equal(bs_sub_ref$w, bs$w)))
stopifnot(isTRUE(all.equal(bs_sub_ref$level, bs$level[2:3])))
stopifnot(!isTRUE(all.equal(bs_sub_own$w, bs$w)))   # the bug this guards against
stopifnot(inherits(try(ng_multitrait_index_basis(M, "maximize", c(0.6, 0.4)),
                       silent = TRUE), "try-error"))     # direction length mismatch
cat("ng_multitrait_index_basis test passed\n")

# --- ng_multitrait_index_upside: exact w'Sw, off-diagonal included, vpm rescaling ---
cr <- data.frame(wf_var_a = c(4, 9), wf_var_b = c(1, 4), wf_cov_a_b = c(-1, 2))
w2 <- c(a = 2, b = 3)
up <- ng_multitrait_index_upside(cr, c("a", "b"), w2)
expect <- sqrt(4 * c(4, 9) + 9 * c(1, 4) + 2 * 2 * 3 * c(-1, 2))
stopifnot(isTRUE(all.equal(up, expect)))
# the off-diagonal must actually change the answer (a negative covariance SHRINKS the index
# spread -- the whole reason the exact within-family covariance is used)
diag_only <- sqrt(4 * c(4, 9) + 9 * c(1, 4))
stopifnot(up[[1L]] < diag_only[[1L]], up[[2L]] > diag_only[[2L]])
# reversed column naming (wf_cov_b_a) resolves to the same covariance
cr_rev <- data.frame(wf_var_a = c(4, 9), wf_var_b = c(1, 4), wf_cov_b_a = c(-1, 2))
stopifnot(isTRUE(all.equal(ng_multitrait_index_upside(cr_rev, c("a", "b"), w2), expect)))
# vpm rescaling: S is rescaled so its diagonal reproduces the authoritative per-trait vpm
vpm <- cbind(a = c(16, 9), b = c(1, 4))          # trait a's variance doubled in SD terms
up_rs <- ng_multitrait_index_upside(cr, c("a", "b"), w2, vpm = vpm)
r <- c(sqrt(16 / 4), 1)                            # per-row scale factors for cross 1
exp1 <- sqrt(4 * r[[1L]]^2 * 4 + 9 * 1 * 1 + 2 * 2 * 3 * r[[1L]] * 1 * (-1))
stopifnot(isTRUE(all.equal(up_rs[[1L]], exp1)))
stopifnot(isTRUE(all.equal(up_rs[[2L]], up[[2L]])))   # cross 2 unchanged (vpm == wf_var)
stopifnot(inherits(try(ng_multitrait_index_upside(cr, c("a", "zz"), w2), silent = TRUE),
                   "try-error"))                      # missing variance column is loud
stopifnot(all(ng_multitrait_index_upside(cr[0, , drop = FALSE], c("a", "b"), w2) == numeric(0)))
cat("ng_multitrait_index_upside test passed\n")

# --- attribution: shares sum to 1, risk driver is the argmax, negatives are meaningful ---
# Variance shares w_k(Sw)_k / w'Sw must reconstruct the reported upside exactly, or the
# attribution would be describing a different S from the one on the axis.
sw <- ng_mt_index_sw(cr, c("a", "b"), w2)
stopifnot(isTRUE(all.equal(sqrt(pmax(sw$total, 0)), up)))
vs <- ng_multitrait_variance_shares(cr, c("a", "b"), w2)
stopifnot(isTRUE(all.equal(unname(rowSums(vs)), c(1, 1))))
stopifnot(identical(colnames(vs), c("a", "b")))
# cross 1 carries a NEGATIVE covariance, so b's share must fall below its diagonal-only share
diag_share_b <- (9 * 1) / (4 * 4 + 9 * 1)
stopifnot(vs[1L, "b"] < diag_share_b)
# a trait antagonistic enough to REMOVE spread from the index shows a negative share -- that is
# the informative case ("protein narrows the index spread"), not an error
cr_anti <- data.frame(wf_var_a = 1, wf_var_b = 1, wf_cov_a_b = -0.9)
vs_anti <- ng_multitrait_variance_shares(cr_anti, c("a", "b"), c(a = 1, b = 0.5))
stopifnot(vs_anti[1L, "b"] < 0, isTRUE(all.equal(unname(rowSums(vs_anti)), 1)))

pv2 <- cbind(a = c(1, 4), b = c(2, 8))
ps <- ng_multitrait_pev_shares(pv2, w2, c("a", "b"))
stopifnot(isTRUE(all.equal(unname(rowSums(ps)), c(1, 1))))
stopifnot(all(ps >= 0), isTRUE(all.equal(unname(ps[1L, "a"]), 4 / 22)))  # no cross terms
drv <- ng_multitrait_risk_driver(ps, c("a", "b"))
stopifnot(identical(drv$trait, c("b", "b")))
stopifnot(isTRUE(all.equal(drv$share, c(18 / 22, 72 / 88))))
# a cross with no usable per-trait PEV has no driver, rather than a spurious first trait
drv_na <- ng_multitrait_risk_driver(matrix(c(NA_real_, 0.4, NA_real_, 0.6), 2, 2),
                                    c("a", "b"))
stopifnot(is.na(drv_na$trait[[1L]]), identical(drv_na$trait[[2L]], "b"))
cat("multi-trait attribution test passed\n")

# --- ng_multitrait_index_pev: block-diagonal weighted sum, NA propagates ---
pv <- cbind(a = c(1, 4), b = c(2, 8))
stopifnot(isTRUE(all.equal(ng_multitrait_index_pev(pv, w2), c(4 * 1 + 9 * 2, 4 * 4 + 9 * 8))))
stopifnot(is.na(ng_multitrait_index_pev(cbind(a = c(1, NA), b = c(2, 8)), w2)[[2L]]))
stopifnot(inherits(try(ng_multitrait_index_pev(pv, c(1, 2, 3)), silent = TRUE), "try-error"))
cat("ng_multitrait_index_pev test passed\n")

# --- e2e: a 2-trait run gets the SAME six columns the single-trait run gets ---
set.seed(11)
n <- 20L; mk <- 80L; gid <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * mk, 1, 0.5), n, mk,
             dimnames = list(gid, sprintf("M%03d", seq_len(mk))))
# antagonistic traits: shared QTL with opposite signs -> negative within-family covariance
b1 <- rnorm(mk, 0, 0.1)
b2 <- -0.8 * b1 + rnorm(mk, 0, 0.03)
y1 <- as.numeric(gm %*% b1) + rnorm(n, 0, 0.5)
y2 <- as.numeric(gm %*% b2) + rnorm(n, 0, 0.5)
genotype  <- data.frame(NAME = gid, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = gid, yield = y1, protein = y2, stringsAsFactors = FALSE)
runmm <- data.frame(SNP = colnames(gm), chr = rep(1:2, length.out = mk),
                    bp = rep(seq(0, 100, length.out = 40), 2)[seq_len(mk)] * 1e6)
dir2 <- data.frame(trait = c("yield", "protein"), column = c("yield", "protein"),
                   direction = c("increase", "increase"), weight = c(0.6, 0.4))
res <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = dir2,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 10L,
  max_crosses_per_parent = 4L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
want <- c("cross_level", "cross_upside", "cross_confidence", "risk_bin",
          "confidence_method", "portfolio_profile")
sel <- res$selected_crosses; cand <- res$candidate_crosses
stopifnot(all(want %in% names(sel)), all(want %in% names(cand)))
stopifnot(all(sel$cross_upside >= 0), any(is.finite(sel$cross_level)))
stopifnot(sel$confidence_method[[1L]] == "midparent_pev_index_partial")
stopifnot(diff(range(cand$cross_confidence, na.rm = TRUE)) > 0)
stopifnot(all(levels(sel$portfolio_profile) ==
                c("breakthrough", "workhorse", "long_shot", "deprioritize")))
d <- res$priority_risk_diagnostics
stopifnot(!is.null(d), identical(d$basis, "multi_trait_index"), d$n_traits == 2L,
          identical(d$upside_method, "exact_within_family_cov"),
          length(d$index_weights) == 2L)

# the exact within-family covariance rode the cross table and reproduces the reported vpm
stopifnot(all(c("wf_var_yield", "wf_var_protein", "wf_cov_yield_protein") %in% names(cand)))
stopifnot(max(abs(cand$wf_var_yield - cand$yield_vpm) / pmax(cand$yield_vpm, 1e-12)) < 1e-8)
stopifnot(mean(cand$wf_cov_yield_protein) < 0)      # traits really are antagonistic

# cross_upside is exactly sqrt(w'Sw), off-diagonal included
wv <- unlist(d$index_weights)
v_exact <- wv[["yield"]]^2 * cand$wf_var_yield + wv[["protein"]]^2 * cand$wf_var_protein +
  2 * wv[["yield"]] * wv[["protein"]] * cand$wf_cov_yield_protein
stopifnot(max(abs(cand$cross_upside - sqrt(pmax(v_exact, 0)))) < 1e-8)
# ...and the diagonal approximation it replaces is materially different here
v_diag <- wv[["yield"]]^2 * cand$wf_var_yield + wv[["protein"]]^2 * cand$wf_var_protein
stopifnot(median(sqrt(v_diag) / cand$cross_upside) > 1.2)

# cross_level is exactly affine in w'm (mid-parent GEBVs, raw trait units)
lin <- wv[["yield"]] * cand$yield_mean_gebv + wv[["protein"]] * cand$protein_mean_gebv
stopifnot(abs(stats::cor(cand$cross_level, lin) - 1) < 1e-10)

# REGRESSION GUARD: the plan and the candidate pool must be on the SAME axes. Each table
# derived its own IQR-based scale before the shared-basis fix, silently putting the selected
# crosses on a different index from the pool they were drawn out of.
key <- function(df) paste(pmin(df$parent1, df$parent2), pmax(df$parent1, df$parent2), sep = "|")
idx <- match(key(sel), key(cand))
stopifnot(!anyNA(idx))
stopifnot(max(abs(sel$cross_level  - cand$cross_level[idx]))  < 1e-8)
stopifnot(max(abs(sel$cross_upside - cand$cross_upside[idx])) < 1e-8)
cat("multi-trait portfolio + risk e2e test passed\n")

# --- portfolio_basis: rank indices are stamped as linearized, linear indices as exact ---
# `auto` promotes to `weighted` when trait weights are present, so multi_trait_score is a RANK
# index: the quadrant is defined by reinterpreting the trait weights as standardized-unit
# coefficients, and must be surfaced as indicative rather than as a decomposition of the merit.
stopifnot(identical(sel$portfolio_basis[[1L]], "linearized_rank_index"))
stopifnot(identical(d$portfolio_basis, "linearized_rank_index"))
stopifnot(is.character(d$portfolio_basis_note), nzchar(d$portfolio_basis_note))
stopifnot(grepl("economic_weight", d$portfolio_basis_note))   # tells the breeder the way out
stopifnot(grepl("linearized_rank_index", sel$cross_level_rule[[1L]]))

# economic weights promote `auto` to a SOLVED linear index -> the decomposition is exact
dir_econ <- data.frame(trait = c("yield", "protein"), column = c("yield", "protein"),
                       direction = c("increase", "increase"), economic_weight = c(3, 2))
res_e <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = dir_econ,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 10L,
  max_crosses_per_parent = 4L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
de <- res_e$priority_risk_diagnostics
stopifnot(identical(de$index_method, "economic_index"))
stopifnot(identical(de$portfolio_basis, "linear_index"))
stopifnot(is.na(de$portfolio_basis_note))
stopifnot(all(c("cross_level", "cross_upside") %in% names(res_e$selected_crosses)))
# basis classifier is driven by the method family, not by column presence
stopifnot(ng_multitrait_portfolio_basis("economic_index") == "linear_index",
          ng_multitrait_portfolio_basis("desired_gain")   == "linear_index",
          ng_multitrait_portfolio_basis("weighted")       == "linearized_rank_index",
          ng_multitrait_portfolio_basis("auto")           == "linearized_rank_index",
          ng_multitrait_portfolio_basis("threshold")      == "linearized_rank_index")
cat("multi-trait portfolio_basis test passed\n")

# --- e2e: a minimized trait must pull the index level DOWN ---
dir_min <- data.frame(trait = c("yield", "protein"), column = c("yield", "protein"),
                      direction = c("increase", "decrease"), weight = c(0.6, 0.4))
res_min <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = dir_min,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 10L,
  max_crosses_per_parent = 4L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
wm <- unlist(res_min$priority_risk_diagnostics$index_weights)
stopifnot(wm[["yield"]] > 0, wm[["protein"]] < 0)
cm <- res_min$candidate_crosses
stopifnot(stats::cor(cm$cross_level, cm$protein_mean_gebv) < 0)   # lower protein = better
stopifnot(all(cm$cross_upside >= 0))                              # variance stays positive
cat("multi-trait minimize-direction test passed\n")

# --- 3 traits: all pairwise covariance columns present and used ---
y3 <- as.numeric(gm %*% rnorm(mk, 0, 0.08)) + rnorm(n, 0, 0.5)
phen3 <- data.frame(NAME = gid, yield = y1, protein = y2, height = y3,
                    stringsAsFactors = FALSE)
dir3 <- data.frame(trait = c("yield", "protein", "height"),
                   column = c("yield", "protein", "height"),
                   direction = c("increase", "increase", "decrease"),
                   weight = c(0.5, 0.3, 0.2))
res3 <- ng_run_cross_prediction(
  phenotype = phen3, genotype = genotype, marker_map = runmm, trait_direction = dir3,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 10L,
  max_crosses_per_parent = 4L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
c3 <- res3$candidate_crosses
stopifnot(all(c("wf_cov_yield_protein", "wf_cov_yield_height", "wf_cov_protein_height")
              %in% names(c3)))
w3 <- unlist(res3$priority_risk_diagnostics$index_weights)
v3 <- w3[["yield"]]^2 * c3$wf_var_yield + w3[["protein"]]^2 * c3$wf_var_protein +
  w3[["height"]]^2 * c3$wf_var_height +
  2 * w3[["yield"]] * w3[["protein"]] * c3$wf_cov_yield_protein +
  2 * w3[["yield"]] * w3[["height"]] * c3$wf_cov_yield_height +
  2 * w3[["protein"]] * w3[["height"]] * c3$wf_cov_protein_height
stopifnot(max(abs(c3$cross_upside - sqrt(pmax(v3, 0)))) < 1e-8)
stopifnot(res3$priority_risk_diagnostics$n_traits == 3L)
cat("multi-trait 3-trait test passed\n")

# --- risk attribution end to end -----------------------------------------------------------
# There is deliberately NO weakest-link floor on cross_confidence: sum_k w_k^2 PEV_k is already
# the correct propagation of estimation error into the index (a lightly weighted trait with a
# large PEV still dominates the sum on its own), and the "unacceptable on a minor trait"
# question belongs to trait_checks / min_value-max_value, not to confidence. What IS provided is
# attribution: which trait carries the risk, and whether the index confidence is really an
# index quantity at all.
phen_n <- data.frame(NAME = gid, yield = y1, protein = y2, noise = rnorm(n),
                     stringsAsFactors = FALSE)
dir_n <- data.frame(trait = c("yield", "protein", "noise"),
                    column = c("yield", "protein", "noise"),
                    direction = c("increase", "increase", "increase"),
                    weight = c(0.5, 0.4, 0.1))
res_n <- ng_run_cross_prediction(
  phenotype = phen_n, genotype = genotype, marker_map = runmm, trait_direction = dir_n,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 8L,
  max_crosses_per_parent = 3L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
dn <- res_n$priority_risk_diagnostics; cn2 <- res_n$candidate_crosses
stopifnot(all(c("risk_driver_trait", "risk_driver_share") %in% names(cn2)))
stopifnot(all(c("risk_driver_trait", "risk_driver_share") %in% names(res_n$selected_crosses)))
stopifnot(all(cn2$risk_driver_share >= 0 & cn2$risk_driver_share <= 1, na.rm = TRUE))
stopifnot(all(cn2$risk_driver_trait %in% c("yield", "protein", "noise")))
# the driver is exactly the argmax of the per-trait PEV shares -- no floor, no reweighting
ps_n <- ng_multitrait_pev_shares(
  as.matrix(cn2[, paste0(c("yield", "protein", "noise"), "_midparent_pev")]),
  unlist(dn$index_weights), c("yield", "protein", "noise"))
stopifnot(identical(cn2$risk_driver_trait,
                    c("yield", "protein", "noise")[max.col(ps_n, ties.method = "first")]))
stopifnot(max(abs(cn2$risk_driver_share - apply(ps_n, 1L, max))) < 1e-12)
stopifnot(max(abs(rowSums(ps_n) - 1)) < 1e-8)

# per-trait index accounting
stopifnot(length(dn$index_traits) == 3L)
nm <- vapply(dn$index_traits, `[[`, character(1L), "trait")
stopifnot(identical(nm, c("yield", "protein", "noise")))
get1 <- function(trait, field) dn$index_traits[[match(trait, nm)]][[field]]
# shares are means over the candidate pool, so each set sums to ~1
stopifnot(abs(sum(vapply(dn$index_traits, `[[`, numeric(1L), "mean_variance_share")) - 1) < 1e-8)
stopifnot(abs(sum(vapply(dn$index_traits, `[[`, numeric(1L), "mean_pev_share")) - 1) < 1e-8)
# the pure-noise trait really is the worst-predicted one
stopifnot(get1("noise", "marker_effect_reliability") <
            get1("yield", "marker_effect_reliability"))

# Concentration guard. Per-trait PEVs are only comparable across traits when their residual
# variances are, and sigma_e^2 is estimated in sample -- a trait whose ridge lambda hits the
# grid floor interpolates its training rows and reports a near-zero PEV. When one trait then
# carries essentially all of the index PEV, risk_bin is a single-trait statement and the run
# must say so. Assert the invariant (note present exactly when concentration is extreme)
# rather than a particular trait, which depends on which lambda the CV happens to pick.
stopifnot(is.numeric(dn$pev_concentration), dn$pev_concentration <= 1 + 1e-8)
stopifnot(identical(is.na(dn$pev_concentration_note), !(dn$pev_concentration > 0.9)))
if (dn$pev_concentration > 0.9) {
  stopifnot(grepl("ridge_lambda", dn$pev_concentration_note),
            grepl(cn2$risk_driver_trait[[1L]], dn$pev_concentration_note, fixed = TRUE))
}
# the flag list is well-formed whatever it contains, and only ever names real traits
stopifnot(is.list(dn$risk_disproportionate_traits),
          all(unlist(dn$risk_disproportionate_traits) %in% nm))
cat("multi-trait risk attribution test passed\n")
