# The marker-effect reliability verdict.
#
# One pure function of evidence the run already has -- per-trait cv_predictive_r2,
# whether the cross-validation could run at all, and the resolved metric. It is
# deliberately pure so the gate inside the runner and the standalone preview cannot
# disagree: both call this, and a preview that could diverge from the run it
# previews would be a lie about the run rather than a preview of it.
#
# THE POLICY. Three tiers, evaluated per trait:
#
#   unevaluable (cv_predictive_r2 is NA), or cv_predictive_r2 < 0
#       -> "refuse". Unevaluable means the fit has no cross-validatable basis at
#          all; negative means it predicts measurably worse than simply using the
#          trait mean. Ranking crosses on a within-family variance derived from
#          either is not defensible, and this package's stated posture is refusal
#          over silent substitution.
#
#   0 <= cv_predictive_r2 < min_cv_predictive_r2
#       -> "fallback". The cross mean becomes the phenotypic mid-parent. This is
#          what already happened; the verdict exists so it can be ANNOUNCED
#          instead of discovered afterwards in a spreadsheet.
#
#   cv_predictive_r2 >= min_cv_predictive_r2
#       -> "gebv".
#
# WHICH METRICS IT APPLIES TO. Exactly those whose delivered value contains a
# within-family variance -- and every within-family variance in this package is a
# quadratic form in beta-hat (VPM, PMV fast and full-posterior, the phased het
# kernel, both polyploid paths). Two metrics are exempt:
#
#   parent_distance   uses only the GRM; no beta-hat anywhere
#   mean              has no variance term, and its mid-parent falls back to the
#                     phenotype when markers fail, so the delivered value carries
#                     no marker content
#
# Both are honest non-genomic analyses, correctly labelled -- not loopholes. What
# is refused is a *usefulness* whose variance comes from a model the engine has
# already rejected for the mean. That incoherence is the reason this gate exists.
#
# WHY THE TWO REFUSAL REASONS ARE KEPT APART. "cross-validation could not run"
# (n < 10 and other degeneracies) and "predicts worse than the trait mean" are
# different situations with different remedies, and only the first is fixed by
# supplying a training set. Collapsing them into one NA is what made the remediable
# case invisible.

# Metrics whose delivered value contains a marker-derived within-family variance.
# `var_complex` is absent deliberately: it is rewritten to usefulness + pmv in
# ng_cp__build_ctx() before this is ever called, so only resolved tokens appear.
ng_variance_metrics <- function() c("usefulness", "pmv", "vpm")

ng_metric_uses_marker_variance <- function(trait_value_metric) {
  isTRUE(as.character(trait_value_metric)[[1L]] %in% ng_variance_metrics())
}

ng_evaluate_marker_reliability <- function(effect_summary,
                                           trait_value_metric = "usefulness",
                                           min_cv_predictive_r2 = 0.35) {
  es <- as.data.frame(effect_summary, stringsAsFactors = FALSE)
  if (!nrow(es)) {
    return(data.frame(trait = character(), cv_predictive_r2 = numeric(), n = integer(),
                      verdict = character(), reason = character(), remedy = character(),
                      stringsAsFactors = FALSE))
  }
  metric <- as.character(trait_value_metric)[[1L]]
  gated <- ng_metric_uses_marker_variance(metric)
  # parent_distance uses only the GRM: marker reliability is not merely tolerable
  # here, it is IRRELEVANT. Reporting a "phenotype fallback" for it would describe a
  # mean it does not compute. `mean` is different -- it has no variance term, but its
  # mid-parent genuinely does fall back, so it keeps the fallback verdict.
  if (identical(metric, "parent_distance")) {
    return(data.frame(
      trait = as.character(es$trait %||% seq_len(nrow(es))),
      cv_predictive_r2 = suppressWarnings(as.numeric(es$cv_predictive_r2 %||% NA_real_)),
      n = suppressWarnings(as.integer(es$marker_effect_training_n %||% NA_integer_)),
      verdict = "not_applicable",
      reason = "trait_value_metric = 'parent_distance' ranks on genomic relationship distance and uses no marker effects, so marker predictive ability does not bear on it",
      remedy = NA_character_, stringsAsFactors = FALSE))
  }
  thresh <- suppressWarnings(as.numeric(min_cv_predictive_r2)[[1L]])

  trait <- as.character(es$trait %||% seq_len(nrow(es)))
  r2 <- suppressWarnings(as.numeric(es$cv_predictive_r2 %||% rep(NA_real_, nrow(es))))
  n <- suppressWarnings(as.integer(es$marker_effect_training_n %||% rep(NA_integer_, nrow(es))))
  avail <- if (!is.null(es$cv_available)) as.logical(es$cv_available) else is.finite(r2)
  why <- as.character(es$cv_unavailable_reason %||% rep(NA_character_, nrow(es)))

  verdict <- character(nrow(es)); reason <- character(nrow(es)); remedy <- character(nrow(es))
  for (i in seq_len(nrow(es))) {
    unevaluable <- !isTRUE(avail[i]) || !is.finite(r2[i])
    if (unevaluable) {
      verdict[i] <- if (gated) "refuse" else "fallback"
      reason[i] <- sprintf(
        "cross-validation could not be computed (%s, n = %s), so cv_predictive_r2 is unevaluable and genomic means are impossible regardless of marker quality",
        if (is.na(why[i])) "unevaluable" else why[i], ifelse(is.na(n[i]), "unknown", n[i]))
      remedy[i] <- "supply a marker-effect training set (training_genotype / training_phenotype) to raise the number of genotyped+phenotyped records above the cross-validation minimum"
    } else if (r2[i] < 0) {
      verdict[i] <- if (gated) "refuse" else "fallback"
      reason[i] <- sprintf(
        "cv_predictive_r2 = %.4g is negative: the marker model predicts worse than the trait mean",
        r2[i])
      remedy[i] <- "supply a training set, or choose trait_value_metric = 'parent_distance' (marker-free) or 'mean' (phenotypic mid-parent)"
    } else if (r2[i] < thresh) {
      verdict[i] <- "fallback"
      reason[i] <- sprintf(
        "cv_predictive_r2 = %.4g is below min_cv_predictive_r2 = %.4g, so the cross mean falls back to the phenotypic mid-parent",
        r2[i], thresh)
      remedy[i] <- "supply a training set to sharpen the marker effects, or lower min_cv_predictive_r2 deliberately"
    } else {
      verdict[i] <- "gebv"
      reason[i] <- sprintf("cv_predictive_r2 = %.4g meets min_cv_predictive_r2 = %.4g", r2[i], thresh)
      remedy[i] <- NA_character_
    }
  }
  data.frame(trait = trait, cv_predictive_r2 = r2, n = n,
             verdict = verdict, reason = reason, remedy = remedy,
             stringsAsFactors = FALSE)
}
