# One control decides whether the marker effects are good enough for a trait.
#
# Reliability, in this package, was only ever a way to ask "is the marker-effect
# estimate for this trait any good?" -- nothing more. `cv_predictive_r2` answers exactly
# that question, and `min_cv_predictive_r2` is the bar a breeding programme sets for
# itself. `min_effect_reliability` was a second control for the same concept, gating a
# calibrated-reliability branch that nothing in the package can reach, kept as a hook for
# a PEV-based reliability that does not exist.
#
# Two adjacent controls for one concept, one of them inert, is the duplication this
# package has now removed twice: posterior_predictions$mean_source vs
# effect_summary$mean_source, and prediction_mode = "index_as_trait" vs value_kind. Both
# times the fix was to make ONE field answer the question. Reserving a hook for a
# quantity that may never arrive is not a third exception -- it is the same thing with a
# longer horizon.
#
# So min_effect_reliability is deprecated: it warns, it is ignored, and it is gone from
# the contract a frontend renders from.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(31); n <- 40L; m <- 40L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
gv <- as.numeric(g %*% c(rnorm(8, 0, 1), rep(0, m - 8)))
y  <- gv + rnorm(n, 0, 0.2 * stats::sd(gv))
run <- function(...) ng_run_cross_prediction(
  phenotype = data.frame(NAME = ids, yield = y),
  genotype = data.frame(NAME = ids, g, check.names = FALSE),
  marker_map = data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                          bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6),
  trait_direction = data.frame(trait = "yield", column = "yield", direction = "increase"),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 8L, write_outputs = FALSE, write_figures = FALSE,
  run_posterior_prediction = FALSE, seed = 5L, ...)

# ---- 1. supplying it warns, and names the control that replaced it -------
w <- NULL
r_dep <- withCallingHandlers(run(min_effect_reliability = 0.9),
  warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
hit <- grep("min_effect_reliability", w %||% character(), value = TRUE)
stopifnot(length(hit) > 0L)
stopifnot(grepl("deprecat", hit[[1L]], ignore.case = TRUE))
stopifnot(grepl("min_cv_predictive_r2", hit[[1L]], fixed = TRUE))

# ---- 2. and it is IGNORED, not merely announced --------------------------
# A deprecation that still changes the run is not a deprecation. 0.9 is a bar the
# fixture cannot clear, so if the parameter still governed anything the two runs would
# differ. They must not.
r_plain <- run()
stopifnot(identical(r_dep$effect_summary$mean_source[[1L]],
                    r_plain$effect_summary$mean_source[[1L]]))
stopifnot(isTRUE(all.equal(r_dep$candidate_crosses$yield_mean,
                           r_plain$candidate_crosses$yield_mean)))

# ---- 3. the live control still does govern the same run ------------------
# Proves the run is capable of responding to a threshold at all, so check 2 is evidence
# of the deprecation rather than of an inert fixture.
r_hi <- suppressWarnings(run(min_cv_predictive_r2 = 1.1))
stopifnot(!identical(r_hi$effect_summary$mean_source[[1L]],
                     r_plain$effect_summary$mean_source[[1L]]))

# ---- 4. it is gone from the contract a frontend renders from -------------
ctls <- ng_backend_controls()
ids_ctl <- vapply(ctls, function(c) c$id, character(1))
stopifnot(!("min_effect_reliability" %in% ids_ctl))
stopifnot("min_cv_predictive_r2" %in% ids_ctl)

# ---- 5. and the surviving control is described as the reliability decision
m <- Filter(function(c) identical(c$id, "min_cv_predictive_r2"), ctls)[[1L]]
note <- paste(unlist(m$note), collapse = " ")
stopifnot(grepl("marker", note, ignore.case = TRUE))
# It may NAME the retired control -- that is migration guidance -- but only as
# retired/replaced, never as somewhere to go and set something.
if (grepl("min_effect_reliability", note, fixed = TRUE)) {
  stopifnot(grepl("retired|replaced|deprecat", note, ignore.case = TRUE))
}
# and the rationale is recorded where a reader will meet it: the parents are
# phenotyped, so a PEV reliability for unphenotyped candidates is the wrong quantity
stopifnot(grepl("unphenotyped", note, ignore.case = TRUE))

cat("min_effect_reliability_retired: PASS\n")
