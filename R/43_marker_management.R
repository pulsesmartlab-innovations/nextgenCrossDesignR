# Marker steering and lethal-allele guarding.
#
# Two breeder-facing capabilities at the scoring layer (where the genotype matrix is
# available): (1) marker steering -- drive the frequency of a nominated major-gene allele
# up or down, or aim it at a target frequency; and (2) lethal-allele guarding -- avoid
# carrier x carrier matings at known deleterious recessive loci. Each produces per-cross
# columns the OCS optimizer can use as an extra objective term (marker steering) or as a
# hard candidate filter (lethal guarding). The idea of managing marker/genotype
# frequencies and recessive lethals inside mate selection follows Kinghorn (2011),
# Genet. Sel. Evol. 43:4.
#
# Genotype convention (matches the rest of the package): `geno` is a parents x markers
# dosage matrix with entries in 0..ploidy counting copies of the ALT allele; rownames are
# parent ids and colnames are marker ids. Everything is ploidy-aware, so it works for
# diploids and for the autotetraploid / subgenome model families.

# Specification of markers to steer toward a target. `direction` = "increase" drives the
# ALT-allele frequency up, "decrease" drives it down; alternatively give an explicit
# `target_freq` in [0, 1] to aim the expected progeny ALT frequency at that value.
ng_marker_target_spec <- function(marker,
                                  direction = "increase",
                                  target_freq = NA_real_,
                                  weight = 1) {
  marker <- as.character(marker)
  n <- length(marker)
  if (!n) ng_stop("marker_target_spec needs at least one marker")
  direction <- tolower(as.character(rep_len(direction, n)))
  bad <- !(direction %in% c("increase", "decrease"))
  if (any(bad)) ng_stop("direction must be 'increase' or 'decrease': ", paste(unique(direction[bad]), collapse = ", "))
  data.frame(
    marker = marker,
    direction = direction,
    target_freq = as.numeric(rep_len(target_freq, n)),
    weight = as.numeric(rep_len(weight, n)),
    stringsAsFactors = FALSE
  )
}

# Expected progeny ALT-allele frequency at a marker for a mid-parent cross, and a
# steering score (higher = better) summed over the target markers. For each target
# marker m and cross (p1, p2): freq = (dosage[p1,m] + dosage[p2,m]) / (2 * ploidy).
#   increase        : + weight * freq
#   decrease        : + weight * (1 - freq)
#   target_freq set : - weight * |freq - target_freq|  (closer is better)
# Returns a data.frame with marker_target_score plus one marker_freq_<m> column per
# target marker so the frontend can show per-locus progress.
ng_marker_target_scores <- function(geno, pairs, spec, ploidy = 2) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  ploidy <- as.numeric(ploidy)[[1L]]
  if (!is.finite(ploidy) || ploidy <= 0) ng_stop("ploidy must be a positive number")
  miss <- setdiff(spec$marker, colnames(geno))
  if (length(miss)) ng_stop("marker_target_spec markers not in geno: ", paste(miss, collapse = ", "))
  p1 <- match(as.character(pairs$parent1), rownames(geno))
  p2 <- match(as.character(pairs$parent2), rownames(geno))
  if (anyNA(p1) || anyNA(p2)) ng_stop("some cross parents are not rows of geno")
  n <- nrow(pairs)
  score <- numeric(n)
  freq_cols <- list()
  for (j in seq_len(nrow(spec))) {
    m <- spec$marker[[j]]
    d1 <- geno[p1, m]
    d2 <- geno[p2, m]
    freq <- (d1 + d2) / (2 * ploidy)
    freq[!is.finite(freq)] <- NA_real_
    w <- spec$weight[[j]]
    tgt <- spec$target_freq[[j]]
    contrib <- if (is.finite(tgt)) {
      -w * abs(freq - tgt)
    } else if (identical(spec$direction[[j]], "decrease")) {
      w * (1 - freq)
    } else {
      w * freq
    }
    contrib[!is.finite(contrib)] <- 0
    score <- score + contrib
    freq_cols[[paste0("marker_freq_", m)]] <- freq
  }
  out <- data.frame(marker_target_score = score, stringsAsFactors = FALSE)
  for (nm in names(freq_cols)) out[[nm]] <- freq_cols[[nm]]
  out
}

# Specification of lethal / strongly-deleterious recessive loci to manage. `risk_allele`
# names which allele is the deleterious one: "alt" (dosage counts risk copies directly)
# or "ref" (risk copies = ploidy - dosage).
ng_lethal_recessive_spec <- function(marker, risk_allele = "alt") {
  marker <- as.character(marker)
  n <- length(marker)
  if (!n) ng_stop("lethal_recessive_spec needs at least one marker")
  risk_allele <- tolower(as.character(rep_len(risk_allele, n)))
  bad <- !(risk_allele %in% c("alt", "ref"))
  if (any(bad)) ng_stop("risk_allele must be 'alt' or 'ref'")
  data.frame(marker = marker, risk_allele = risk_allele, stringsAsFactors = FALSE)
}

# Copies of the risk allele carried by each parent at each lethal locus (parents x loci).
ng_lethal_risk_dosage <- function(geno, spec, ploidy = 2) {
  geno <- ng_as_numeric_matrix(geno, "geno")
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  miss <- setdiff(spec$marker, colnames(geno))
  if (length(miss)) ng_stop("lethal_recessive_spec markers not in geno: ", paste(miss, collapse = ", "))
  d <- geno[, spec$marker, drop = FALSE]
  ref <- spec$risk_allele == "ref"
  if (any(ref)) d[, ref] <- as.numeric(ploidy) - d[, ref]
  colnames(d) <- spec$marker
  d
}

# Flag crosses that risk homozygous-affected progeny: both parents carry >= 1 copy of the
# risk allele at the same lethal locus (covers carrier x carrier, carrier x affected and
# affected x affected). Returns lethal_carrier_cross (logical) and lethal_risk_loci
# (count of shared at-risk loci) per cross.
ng_lethal_recessive_cross_risk <- function(geno, pairs, spec, ploidy = 2) {
  rd <- ng_lethal_risk_dosage(geno, spec, ploidy = ploidy)
  p1 <- match(as.character(pairs$parent1), rownames(rd))
  p2 <- match(as.character(pairs$parent2), rownames(rd))
  if (anyNA(p1) || anyNA(p2)) ng_stop("some cross parents are not rows of geno")
  carrier1 <- rd[p1, , drop = FALSE] >= 1
  carrier2 <- rd[p2, , drop = FALSE] >= 1
  both <- carrier1 & carrier2
  loci <- rowSums(both, na.rm = TRUE)
  data.frame(
    lethal_carrier_cross = loci > 0,
    lethal_risk_loci = as.integer(loci),
    stringsAsFactors = FALSE
  )
}

# Convenience: attach marker-management columns to a scored candidate table and, if
# requested, (a) blend the marker-target score into a new marker_adjusted_gain column and
# (b) drop carrier x carrier crosses so they never enter the mating plan. `scores` must
# carry parent1/parent2; `gain_col` is the base merit column to blend against. Returns the
# augmented (and optionally filtered) scores table.
ng_apply_marker_management <- function(scores,
                                       geno,
                                       marker_target_spec = NULL,
                                       lethal_spec = NULL,
                                       ploidy = 2,
                                       gain_col = "usefulness_pmv_gebv",
                                       lambda_marker = 0,
                                       drop_lethal_carrier_crosses = TRUE) {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(scores))) {
    ng_stop("scores must have parent1 and parent2 columns")
  }
  if (!is.null(marker_target_spec)) {
    mt <- ng_marker_target_scores(geno, scores, marker_target_spec, ploidy = ploidy)
    for (nm in names(mt)) scores[[nm]] <- mt[[nm]]
    if (is.finite(lambda_marker) && lambda_marker != 0 && gain_col %in% names(scores)) {
      # Blend on a comparable scale: rescale the marker score into the spread of the base
      # merit column so lambda_marker is an interpretable relative emphasis.
      g <- as.numeric(scores[[gain_col]])
      gsd <- stats::sd(g[is.finite(g)])
      msd <- stats::sd(scores$marker_target_score[is.finite(scores$marker_target_score)])
      scale <- if (is.finite(gsd) && is.finite(msd) && msd > 0) gsd / msd else 1
      scores$marker_adjusted_gain <- g + lambda_marker * scale * scores$marker_target_score
    }
  }
  if (!is.null(lethal_spec)) {
    lr <- ng_lethal_recessive_cross_risk(geno, scores, lethal_spec, ploidy = ploidy)
    scores$lethal_carrier_cross <- lr$lethal_carrier_cross
    scores$lethal_risk_loci <- lr$lethal_risk_loci
    if (isTRUE(drop_lethal_carrier_crosses)) {
      scores <- scores[!scores$lethal_carrier_cross, , drop = FALSE]
    }
  }
  scores
}
