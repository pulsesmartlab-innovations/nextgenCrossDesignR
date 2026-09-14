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

# Preview a run's marker effects without scoring a single cross.
#
# Answers "should I launch this?" in seconds where the run it previews can take
# many hours. It runs QC and the marker-effect fits, reports per-trait
# cv_predictive_r2 and the verdict the gate will reach, and stops.
#
# It routes through ng_cp__build_ctx() -> ng_cp__stage_qc() -> the SAME
# ng_cp__fit_trait_effects() the runner uses, and evaluates with the SAME
# ng_evaluate_marker_reliability(). That sharing is the whole point: a preview
# that fitted separately could drift from the run -- a different seed draws a
# different CV fold partition, and a trait near the threshold would flip between
# the two. A preview that can disagree with its run is a second opinion, and the
# wrong one is worse than none.
#
# It reports a refusal rather than raising one. A preview that errored on the very
# configuration it exists to warn about would force the caller to catch an
# exception in order to learn they should not launch.
ng_preview_marker_effects <- function(config) {
  config$run_posterior_prediction <- FALSE
  # ng_cp__build_ctx() expects a COMPLETE formal list -- it match.arg()s every enum,
  # and an absent one arrives as the full choices vector rather than a default. Back
  # -fill from ng_run_cross_prediction()'s own formals so the preview is configured
  # exactly as the run would be, rather than from a second set of defaults that
  # could drift.
  fm <- formals(ng_run_cross_prediction)
  for (nm in setdiff(names(fm), names(config))) {
    v <- fm[[nm]]
    config[nm] <- list(tryCatch(eval(v, envir = environment()), error = function(e) NULL))
  }
  ctx <- ng_cp__build_ctx(config)
  ctx <- ng_cp__stage_qc(ctx)
  if (identical(ctx$qc$status, "blocker")) {
    blockers <- ctx$qc$issues$message[ctx$qc$issues$severity == "blocker"]
    ng_stop("Input QC found blocker issues before the marker-effect preview: ",
            paste(utils::head(blockers, 5L), collapse = " | "))
  }
  # Input validity before model quality, exactly as the runner orders them.
  ng_validate_parent_dosage(ctx$geno, parent_type = ctx$parent_type,
                            ploidy = ctx$marker_ploidy %||% 2L,
                            phased_haplotypes = ctx$phased_haplotypes)

  training_set <- ng_run_cp_training_set(
    training_genotype = ctx$training_genotype, training_phenotype = ctx$training_phenotype,
    training_genotype_file = ctx$training_genotype_file,
    training_phenotype_file = ctx$training_phenotype_file,
    training_genotype_id_col = ctx$training_genotype_id_col %||% ctx$genotype_id_col,
    training_phenotype_id_col = ctx$training_phenotype_id_col %||% ctx$phenotype_id_col,
    parent_ids = ctx$ids, parent_markers = colnames(ctx$geno),
    trait_columns = ctx$trait_spec$column)

  effects <- ng_cp__fit_trait_effects(ctx$trait_spec, ctx$pheno, ctx$geno, ctx$ids,
                                      training_set, ctx$seed)
  thresh <- if (is.null(ctx$effect_gate_min_cv_predictive_r2)) ctx$min_cv_predictive_r2
            else ctx$effect_gate_min_cv_predictive_r2
  verdicts <- ng_evaluate_marker_reliability(effects,
                                             trait_value_metric = ctx$trait_value_metric,
                                             min_cv_predictive_r2 = thresh)
  adv <- ng_advisory_empty()
  for (k in seq_len(nrow(verdicts))) {
    if (verdicts$verdict[[k]] %in% c("gebv", "not_applicable")) next
    adv <- ng_advisory_add(adv,
      id = if (identical(verdicts$verdict[[k]], "refuse")) "would_refuse" else "would_fall_back",
      severity = "advice", stage = "preview", trait = verdicts$trait[[k]],
      message = sprintf("%s: %s. %s", verdicts$trait[[k]], verdicts$reason[[k]],
                        if (is.na(verdicts$remedy[[k]])) "" else verdicts$remedy[[k]]))
  }
  list(schema = "ng_effect_preview.v1",
       trait_value_metric = ctx$trait_value_metric,
       min_cv_predictive_r2 = thresh,
       effects = effects, verdicts = verdicts,
       advisories = adv, advisory_summary = ng_advisory_rollup(adv))
}
