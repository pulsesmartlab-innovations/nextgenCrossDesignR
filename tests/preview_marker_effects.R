# The preview must predict exactly what the run will do.
#
# ng_preview_marker_effects() answers "should I launch this?" in seconds rather
# than hours: it runs QC and the marker-effect fits, reports per-trait
# cv_predictive_r2 and the verdict the gate will reach, and stops before any cross
# is scored. On a real 17-trait panel the run it previews takes ~23 hours.
#
# THE INVARIANT THAT MAKES IT WORTH HAVING. A preview that can disagree with the
# run is not a preview -- it is a second opinion, and the wrong one is worse than
# none. So this asserts the verdicts match trait for trait, and the numbers match
# BITWISE rather than approximately: the preview and the run share one code path
# and one seed, so anything short of identity means they have diverged.
#
# The seed matters more than it looks. The CV fold partition is drawn from the run
# seed, so the same data under a different seed gives a different
# cv_predictive_r2 -- which is how a trait sitting near the threshold can flip
# between a preview and its run. Sharing the path is what removes that risk.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

mk_cfg <- function(n, m, n_causal, noise_mult, metric = "usefulness",
                   thresh = 0.35, seed_data = 77L) {
  set.seed(seed_data)
  g <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
  ids <- sprintf("P%02d", seq_len(n))
  rownames(g) <- ids; colnames(g) <- sprintf("M%02d", seq_len(m))
  b <- c(rnorm(n_causal, sd = 1.5), rep(0, m - n_causal))
  gv <- as.numeric(g %*% b)
  list(genotype = data.frame(NAME = ids, g, check.names = FALSE, stringsAsFactors = FALSE),
       genotype_id_col = "NAME",
       phenotype = data.frame(NAME = ids, YIELD = gv + rnorm(n, sd = noise_mult * stats::sd(gv)),
                              stringsAsFactors = FALSE),
       phenotype_id_col = "NAME", traits_to_use = "YIELD",
       trait_direction = data.frame(Trait = "YIELD", Selection_direction = "increase",
                                    stringsAsFactors = FALSE),
       direction_trait_col = "Trait", direction_column_col = "Trait",
       direction_direction_col = "Selection_direction",
       marker_map = data.frame(SNP = colnames(g), chr = rep(1:2, each = m / 2),
                               cm = rep(seq(0, 70, length.out = m / 2), 2),
                               stringsAsFactors = FALSE),
       map_marker_col = "SNP", map_chr_col = "chr", map_pos_cm_col = "cm",
       map_position_unit = "cM", progeny = "DH", parent_type = "inbred",
       n_crosses = 6L, write_outputs = FALSE, write_figures = FALSE,
       run_posterior_prediction = FALSE, trait_value_metric = metric,
       min_cv_predictive_r2 = thresh, seed = 20260907L)
}

# ---- shape -----------------------------------------------------------------
cfg <- mk_cfg(40L, 60L, 10L, 1.2)
p <- ng_preview_marker_effects(cfg)
stopifnot(is.list(p), identical(p$schema, "ng_effect_preview.v1"))
stopifnot(all(c("effects", "verdicts", "advisories") %in% names(p)))
stopifnot(is.data.frame(p$effects), nrow(p$effects) == 1L)
stopifnot(all(c("trait", "cv_predictive_r2", "cv_available", "marker_effect_training_n")
              %in% names(p$effects)))
stopifnot(all(c("trait", "verdict", "reason", "remedy") %in% names(p$verdicts)))

# ---- it scores nothing -----------------------------------------------------
# The whole point is stopping before the expensive work, so the absence of a cross
# table is the observable -- not a stopwatch.
stopifnot(is.null(p$selected_crosses), is.null(p$candidate_crosses))

# ---- THE INVARIANT: preview == run, bitwise -------------------------------
r <- do.call(ng_run_cross_prediction, cfg)
stopifnot(identical(p$effects$cv_predictive_r2, r$effect_summary$cv_predictive_r2))
stopifnot(identical(as.integer(p$effects$marker_effect_training_n),
                    as.integer(r$effect_summary$marker_effect_training_n)))

# the projected verdict must be the decision the run actually took
proj <- p$verdicts$verdict[[1L]]
actual <- r$effect_summary$mean_source[[1L]]
stopifnot(identical(proj, if (identical(actual, "GEBV")) "gebv" else "fallback"))

# ---- it predicts a refusal without triggering one -------------------------
# A preview that errored on the configuration it is meant to warn you about would
# be useless: you would have to catch an exception to learn you should not launch.
bad <- mk_cfg(30L, 100L, 80L, 8.0)
pb <- ng_preview_marker_effects(bad)
stopifnot(identical(pb$verdicts$verdict[[1L]], "refuse"))
stopifnot(grepl("negative", pb$verdicts$reason[[1L]], fixed = TRUE))
stopifnot(nzchar(pb$verdicts$remedy[[1L]]))
# ...and the run it predicted does indeed refuse.
stopifnot(inherits(tryCatch(do.call(ng_run_cross_prediction, bad),
                            error = function(e) e), "error"))

# ---- the metric governs the verdict, in the preview too -------------------
pd <- ng_preview_marker_effects(mk_cfg(30L, 100L, 80L, 8.0, metric = "parent_distance"))
stopifnot(identical(pd$verdicts$verdict[[1L]], "not_applicable"))

cat("preview_marker_effects: PASS\n")
