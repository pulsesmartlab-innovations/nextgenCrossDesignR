# The reliability verdict: one pure function, three tiers, two exemptions.
#
# ng_evaluate_marker_reliability() turns evidence the run already has -- per-trait
# cv_predictive_r2, whether the CV could run at all, and the resolved metric --
# into a verdict. It is deliberately a PURE function of that evidence so the gate
# inside the runner and the standalone preview cannot disagree: both call this.
#
# THE POLICY
#
#   cv unevaluable (NA: n < 10 and other degeneracies), or cv_predictive_r2 < 0
#       -> "refuse", IF the resolved metric's value contains a within-family
#          variance. Unevaluable means the fit has no cross-validatable basis;
#          negative means it predicts measurably worse than the trait mean.
#          Ranking crosses on a variance derived from either is not defensible.
#
#   0 <= cv_predictive_r2 < threshold
#       -> "fallback": the mean becomes the phenotypic mid-parent. This is what
#          already happens; the verdict exists so it can be announced.
#
#   cv_predictive_r2 >= threshold
#       -> "gebv".
#
# THE EXEMPTIONS. `parent_distance` uses only the GRM, and `mean` has no variance
# term at all -- its mid-parent falls back to phenotype when markers fail, so the
# delivered value carries no marker content. Both are honest non-genomic analyses,
# correctly labelled. What is refused is a *usefulness* whose variance comes from a
# model the engine has already rejected for the mean.
#
# Every verdict must carry a remedy a breeder can act on, which is why the
# unevaluable case is distinguished from the measured-and-failed case: only the
# first is fixed by supplying a training set.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

es <- function(...) {
  d <- data.frame(..., stringsAsFactors = FALSE)
  if (is.null(d$cv_available)) d$cv_available <- !is.na(d$cv_predictive_r2)
  if (is.null(d$cv_unavailable_reason))
    d$cv_unavailable_reason <- ifelse(d$cv_available, NA_character_, "n_lt_10")
  if (is.null(d$marker_effect_training_n)) d$marker_effect_training_n <- 40L
  d
}

ev <- function(summary, metric = "usefulness", thresh = 0.35) {
  ng_evaluate_marker_reliability(summary, trait_value_metric = metric,
                                 min_cv_predictive_r2 = thresh)
}

# ---- shape -----------------------------------------------------------------
v <- ev(es(trait = c("A", "B", "C"), cv_predictive_r2 = c(0.5, 0.2, -0.01)))
stopifnot(is.data.frame(v), nrow(v) == 3L)
stopifnot(all(c("trait", "cv_predictive_r2", "n", "verdict", "reason", "remedy") %in% names(v)))
stopifnot(all(v$verdict %in% c("gebv", "fallback", "refuse")))

# ---- the three tiers -------------------------------------------------------
stopifnot(identical(v$verdict, c("gebv", "fallback", "refuse")))

# negative R2 refuses: worse than predicting the trait mean
stopifnot(identical(ev(es(trait = "X", cv_predictive_r2 = -0.002))$verdict, "refuse"))
# exactly zero is not negative -- it is useless, not harmful; fallback, not refusal
stopifnot(identical(ev(es(trait = "X", cv_predictive_r2 = 0))$verdict, "fallback"))
# unevaluable refuses, and is NOT reported as a measurement failure
u <- ev(es(trait = "X", cv_predictive_r2 = NA_real_, cv_available = FALSE,
           cv_unavailable_reason = "n_lt_10", marker_effect_training_n = 8L))
stopifnot(identical(u$verdict, "refuse"))
stopifnot(grepl("cross-validation", u$reason, fixed = TRUE) || grepl("unevaluable", u$reason))
stopifnot(!grepl("worse", u$reason))          # that is the OTHER refusal's reason
stopifnot(grepl("training", u$remedy, fixed = TRUE))
stopifnot(identical(u$n, 8L))                 # the remedy must be able to name the count

# a measured failure names its number, not a training set as the only cure
m <- ev(es(trait = "X", cv_predictive_r2 = -0.5))
stopifnot(grepl("worse", m$reason, fixed = TRUE))

# ---- the threshold governs the gebv/fallback boundary ----------------------
stopifnot(identical(ev(es(trait = "X", cv_predictive_r2 = 0.30), thresh = 0.20)$verdict, "gebv"))
stopifnot(identical(ev(es(trait = "X", cv_predictive_r2 = 0.30), thresh = 0.35)$verdict, "fallback"))

# ---- exemptions ------------------------------------------------------------
# parent_distance: GRM only. Nothing here can be refused for marker quality.
for (r2 in list(-0.5, NA_real_, 0.9)) {
  d <- es(trait = "X", cv_predictive_r2 = as.numeric(r2),
          cv_available = !is.na(r2), cv_unavailable_reason = if (is.na(r2)) "n_lt_10" else NA_character_)
  stopifnot(!identical(ev(d, metric = "parent_distance")$verdict, "refuse"))
  stopifnot(!identical(ev(d, metric = "mean")$verdict, "refuse"))
}

# ...but the very same evidence refuses under every variance metric.
for (mt in c("usefulness", "pmv", "vpm")) {
  stopifnot(identical(ev(es(trait = "X", cv_predictive_r2 = -0.1), metric = mt)$verdict, "refuse"))
}

# ---- mixed traits: each judged on its own evidence -------------------------
mix <- ev(es(trait = c("YIELD", "HT", "DP"), cv_predictive_r2 = c(-0.008, 0.39, 0.21)),
          thresh = 0.35)
stopifnot(identical(mix$verdict, c("refuse", "gebv", "fallback")))
stopifnot(identical(mix$trait[mix$verdict == "refuse"], "YIELD"))

cat("reliability_gate_evaluator: PASS\n")
