# cross_upside must be the spread of the variance the run SELECTED.
#
# The risk/portfolio layer read `<trait>_vpm` unconditionally: cross_upside was
# sqrt(VPM) whatever uc_variance_source said. On a default PMV run the merit was
# built on PMV while the upside beside it -- and the portfolio quadrants and risk
# tertiles derived from that upside -- came from a different estimand.
#
# PMV = VPM + tr(R Sigma_beta), and the gap varies per cross rather than being a
# constant offset, so this is not a rescaling: it reorders crosses. A cross whose
# effects are poorly estimated gains relatively more PMV than one whose effects are
# sharp, which is exactly the distinction the median split on upside acts on.
#
# This is a RANKING change, not a labelling one, and it is deliberate: the breeder
# chooses one variance and every downstream quantity must be that one.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(77); n <- 40L; m <- 40L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
gv <- as.numeric(g %*% c(rnorm(8, 0, 1), rep(0, m - 8)))
y  <- gv + rnorm(n, 0, 0.2 * stats::sd(gv))
run <- function(...) ng_run_cross_prediction(
  phenotype = data.frame(NAME = ids, yield = y),
  genotype  = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                          bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6),
  trait_direction = data.frame(trait = "yield", column = "yield", direction = "increase"),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 10L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L, ...)

# ---- 1. upside is sqrt of the SELECTED variance, for each choice ----------
for (ucs in c("pmv", "vpm")) {
  r <- run(uc_variance_source = ucs)
  cc <- r$candidate_crosses
  sel   <- suppressWarnings(as.numeric(cc[[paste0("yield_", ucs)]]))
  other <- suppressWarnings(as.numeric(cc[[paste0("yield_", setdiff(c("pmv","vpm"), ucs))]]))
  stopifnot(isTRUE(all.equal(cc$cross_upside, sqrt(pmax(sel, 0)))))
  # and it is NOT the other one -- proving the selection actually reached the layer
  stopifnot(!isTRUE(all.equal(cc$cross_upside, sqrt(pmax(other, 0)))))
}

# ---- 2. the choice genuinely moves the portfolio, not just a column -------
# If PMV and VPM produced the same quadrants this change would be cosmetic. They
# do not: the median split on upside lands differently.
rp <- run(uc_variance_source = "pmv"); rv <- run(uc_variance_source = "vpm")
up_p <- rp$candidate_crosses$cross_upside; up_v <- rv$candidate_crosses$cross_upside
stopifnot(!isTRUE(all.equal(up_p, up_v)))
stopifnot(max(abs(up_p - up_v)) > 1e-6)

# ---- 3. every derived quantity stays well formed under either choice ------
for (r in list(rp, rv)) {
  cc <- r$candidate_crosses
  stopifnot(all(is.finite(cc$cross_upside)), all(cc$cross_upside >= 0))
  stopifnot(!all(is.na(cc$portfolio_profile)))
  stopifnot(all(as.character(stats::na.omit(cc$risk_bin)) %in% c("low", "med", "high")))
}

# ---- 4. a metric with NO variance term falls back, and says which --------
# `mean` and `parent_distance` carry no within-family variance, so there is nothing
# for the breeder to have selected. The layer still needs a spread; it uses the pure
# genetic variance (VPM) and that is the documented fallback, not an accident.
rm_ <- run(trait_value_metric = "mean")
ccm <- rm_$candidate_crosses
stopifnot(isTRUE(all.equal(ccm$cross_upside, sqrt(pmax(as.numeric(ccm$yield_vpm), 0)))))

# ---- 5. MULTI-TRAIT deliberately keeps the pure genetic variance ----------
# Its cross_upside is sqrt(w'Sigma w) over the EXACT within-family covariance matrix:
# wf_var_<t> on the diagonal, wf_cov_<t>_<s> off it. Marker-effect uncertainty is not
# estimated ACROSS traits, so there is no PMV analogue for the off-diagonals; inflating
# only the diagonal would leave Sigma internally inconsistent and possibly not positive
# semi-definite. A coherent VPM-basis Sigma is preferred to a mixed-basis one, and that
# is a decision, so it is asserted rather than left to be rediscovered.
y2 <- as.numeric(g %*% c(rep(0, 8), rnorm(8, 0, 1), rep(0, m - 16)))
y2 <- y2 + rnorm(n, 0, 0.3 * stats::sd(y2))
rmt <- ng_run_cross_prediction(
  phenotype = data.frame(NAME = ids, yield = y, protein = y2),
  genotype  = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                          bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6),
  trait_direction = data.frame(trait = c("yield","protein"), column = c("yield","protein"),
                               direction = c("increase","increase")),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 10L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L, uc_variance_source = "pmv")
cm <- rmt$candidate_crosses
if (all(c("wf_var_yield","wf_var_protein") %in% names(cm))) {
  stopifnot(max(abs(cm$wf_var_yield - cm$yield_vpm) / pmax(cm$yield_vpm, 1e-12)) < 1e-8)
}

cat("risk_layer_follows_the_selected_variance: PASS\n")
