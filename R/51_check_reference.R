# Check lines as REFERENCES, not filters. A check is a benchmark genotype (a released variety,
# a commercial check) that is never crossed: check x check would be a self, which
# ng_make_pairs() does not enumerate, so a check has no variance and can never be a row in the
# cross table. It contributes a scalar per trait, on the same mean_source the cross means use.
# See docs/design/2026-09-03-check-reference-lines-design.md.

# Align a separate check genotype matrix to the marker set already fixed by parent QC. A check
# matrix that silently disagrees on marker set or column order produces a wrong-but-plausible
# reference value, so every disagreement is an error rather than an intersection.
ng_align_check_geno <- function(check_geno, marker_names, ploidy = 2L) {
  marker_names <- as.character(marker_names)
  check_geno <- as.matrix(check_geno)
  if (is.null(colnames(check_geno))) {
    ng_stop("check_geno must have marker column names to verify alignment with the parent ",
            "marker set; an unnamed matrix cannot be checked for allele-order agreement")
  }
  check_geno <- ng_as_numeric_matrix(check_geno, "check_geno")
  if (is.null(rownames(check_geno))) ng_stop("check_geno must have check ids as rownames")
  miss <- setdiff(marker_names, colnames(check_geno))
  if (length(miss)) {
    ng_stop(sprintf(
      "check_geno is missing %d of %d markers used by the fitted effects (e.g. %s). The check ",
      length(miss), length(marker_names),
      paste(utils::head(miss, 5L), collapse = ", ")),
      "file must carry the same marker set, coding, and reference allele as the parent genotypes.")
  }
  out <- check_geno[, marker_names, drop = FALSE]
  ploidy <- as.integer(ploidy[[1L]])
  rng <- suppressWarnings(range(out, na.rm = TRUE))
  if (any(is.finite(rng)) && (min(rng) < 0 || max(rng) > ploidy)) {
    ng_stop(sprintf(
      "check_geno dosage values fall outside [0, %d] (observed %g..%g); the check file must use ",
      ploidy, rng[[1L]], rng[[2L]]),
      "the same allele coding as the parent genotypes")
  }
  out
}

# Resolve each check's value onto the SAME source the cross means used. ng_choose_mean_source()
# returns one of "GEBV", "GEBV_low_reliability", "GEBV_uncalibrated", "BLUP", "BLUE",
# "adjusted_pheno" -- per trait. A GEBV* source predicts from the check's markers with the
# trait's fitted effects; a phenotypic source reads the check's own record. A check with no
# record on a phenotypic source is NA (not evaluable) and must NEVER silently fall back to a
# GEBV, which would put the reference on a different scale from the axis it is drawn on.
ng_check_reference_value <- function(source, check_geno_aligned, effects, check_records = NULL) {
  source <- as.character(source)[[1L]]
  ids <- rownames(check_geno_aligned)
  if (startsWith(source, "GEBV")) {
    # Adapt effects structure for ng_predict_gebv: rename $effect to $beta if needed
    adapted_effects <- effects
    if (is.null(adapted_effects$beta) && !is.null(adapted_effects$effect)) {
      adapted_effects$beta <- adapted_effects$effect
    }
    if (is.null(adapted_effects$intercept)) {
      adapted_effects$intercept <- 0
    }
    return(stats::setNames(ng_predict_gebv(check_geno_aligned, adapted_effects), ids))
  }
  rec <- if (is.null(check_records)) NULL else check_records[[source]]
  if (is.null(rec)) return(stats::setNames(rep(NA_real_, length(ids)), ids))
  v <- suppressWarnings(as.numeric(rec[ids]))
  stats::setNames(v, ids)
}
