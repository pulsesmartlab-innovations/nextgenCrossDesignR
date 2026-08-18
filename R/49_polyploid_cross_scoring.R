# Dominance-aware cross scoring for polyploid dosage data (optional; additive-only is the default).
#
# For CLONAL / outbred crops (cassava, sugarcane, potato) the value of a cross is the distribution
# of TOTAL genotypic value (additive + dominance) among its progeny clones. Predicting that from
# additive + dominance marker effects requires the PROGENY genotype distribution per locus:
#   * a gamete carries ploidy/2 alleles drawn without replacement from the parent's ploidy alleles
#     (random bivalent / polysomic segregation) -> # alt alleles ~ Hypergeometric;
#   * progeny dosage X = gamete_i + gamete_j.
# From X we get, per locus: E[X] (mid-parent dosage), Var(X) (additive segregation), and for the
# digenic dominance covariate H = X*(ploidy-X): E[H] (expected heterozygosity -> heterosis) and
# Var(H) (dominance segregation). The cross statistics are locus sums:
#   cross mean      = mid-parent breeding value + sum_k dom_k*(E[H_k] - hbar_k)   [heterosis]
#   within-fam var  = sum_k add_k^2 Var(X_k)  +  sum_k dom_k^2 Var(H_k)
#   usefulness      = mean + i*sqrt(var)
# The moments depend ONLY on the parental dosage pair, so a (ploidy+1)x(ploidy+1) table is
# precomputed once and looked up vectorized -- no per-cross simulation.
#
# VARIANCE MODEL -- "unlinked_phase_marginalized". The locus sums above set the between-locus
# term to zero (R = I), and that deserves a precise statement rather than the vague label
# "linkage-equilibrium approximation":
#   * It is UNBIASED, not an approximation. Autopolyploid parental phase -- which homologue
#     carries which allele -- is NOT identifiable from dosage. Averaged over the phase
#     configurations consistent with the observed dosages, the between-locus gamete covariance is
#     EXACTLY zero, even for completely linked loci (verified to 1e-16 by enumeration over all
#     phase configurations and all gametes). So sum_k a_k^2 Var(X_k) is the exact expectation of
#     the within-family variance given the information dosage actually carries.
#   * What it cannot do is DISCRIMINATE. For a fixed phase the covariance is not zero: for a
#     duplex x duplex pair of tightly linked loci it ranges over about [-1/3, +1/3]. Two crosses
#     with identical parental dosages but different phase therefore receive identical predictions.
#     Variance-based metrics (usefulness) consequently separate autopolyploid crosses less sharply
#     than the diploid path, where inbred parents make phase known and the exact recombination-
#     aware a'Ra kernel applies.
#   * Supplying phased polyploid haplotypes would make the exact computation possible; that is not
#     implemented. The allopolyploid subgenome path (R/18) IS recombination-aware when a per-
#     subgenome cM map is given, because disomic pairing makes phase tractable there.
# Digenic dominance only (trigenic / quadrigenic dominance components are not modelled).

# Gamete distribution: pmf of the number of alt alleles in a ploidy/2-allele gamete from a parent
# of dosage d, over 0..ploidy/2. Random chromosome segregation is Hypergeometric. DOUBLE REDUCTION
# (autopolyploids, ploidy >= 4) -- two sister-chromatid copies of one parental allele ending up in
# the same gamete -- is added at coefficient `dr`: with probability dr one gamete allele-pair is a
# duplicated (IBD) copy of a random parental allele (2 copies), the remaining ploidy/2 - 2 alleles
# drawn from the other ploidy - 1 alleles; with probability 1 - dr the gamete is pure Hypergeometric.
# dr is ignored for diploids (a 1-allele gamete has no pair to double-reduce).
ng_poly_gamete_pmf <- function(d, ploidy, dr = 0) {
  P <- as.integer(ploidy); g <- P %/% 2L
  base <- stats::dhyper(0:g, m = d, n = P - d, k = g)
  if (dr <= 0 || g < 2L) return(base)
  dr_pmf <- numeric(g + 1L)
  # A = alt (prob d/P): 2 duplicated alt + draw (g-2) from remaining (d-1 alt, P-d ref)
  if (d >= 1L) {
    hj <- stats::dhyper(0:(g - 2L), m = d - 1L, n = P - d, k = g - 2L)
    dr_pmf[(2L + (0:(g - 2L))) + 1L] <- dr_pmf[(2L + (0:(g - 2L))) + 1L] + (d / P) * hj
  }
  # A = ref (prob (P-d)/P): 0 alt from the dup + draw (g-2) from remaining (d alt, P-1-d ref)
  if (d <= P - 1L) {
    hj <- stats::dhyper(0:(g - 2L), m = d, n = P - 1L - d, k = g - 2L)
    dr_pmf[(0:(g - 2L)) + 1L] <- dr_pmf[(0:(g - 2L)) + 1L] + ((P - d) / P) * hj
  }
  (1 - dr) * base + dr * dr_pmf
}

# Precompute progeny-dosage moments for every parental dosage pair (di, dj) in 0..ploidy:
# mu = E[X], varX = Var(X), EH = E[X(ploidy-X)], varH = Var(X(ploidy-X)), covXH = Cov(X,H). Returns
# (P+1)x(P+1) matrices indexed [di+1, dj+1]. `double_reduction` sets the gamete DR coefficient.
ng_polyploid_progeny_moment_table <- function(ploidy, double_reduction = 0) {
  P <- as.integer(ploidy)
  gam <- lapply(0:P, ng_poly_gamete_pmf, ploidy = P, dr = double_reduction)   # gamete pmf per dosage
  x <- 0:P; Hx <- x * (P - x)
  mu <- varX <- EH <- varH <- covXH <- matrix(0, P + 1L, P + 1L)
  for (di in 0:P) for (dj in 0:P) {
    gi <- gam[[di + 1L]]; gj <- gam[[dj + 1L]]
    prog <- as.numeric(stats::convolve(gi, rev(gj), type = "open"))  # pmf of X = gi + gj over 0..P
    prog[prog < 0] <- 0; prog <- prog / sum(prog)
    m1 <- sum(x * prog); m2 <- sum(x^2 * prog)
    eh <- sum(Hx * prog); eh2 <- sum(Hx^2 * prog); exh <- sum(x * Hx * prog)
    mu[di + 1L, dj + 1L]   <- m1
    varX[di + 1L, dj + 1L] <- max(m2 - m1^2, 0)
    EH[di + 1L, dj + 1L]   <- eh
    varH[di + 1L, dj + 1L] <- max(eh2 - eh^2, 0)
    covXH[di + 1L, dj + 1L] <- exh - m1 * eh          # Cov(X, H) at the locus
  }
  list(mu = mu, varX = varX, EH = EH, varH = varH, covXH = covXH)
}

# Score candidate crosses on additive (+ optional dominance) genotypic value using marker effects
# from ng_polyploid_fit_effects. Returns parent1/parent2 + cross_mean (genotypic), add_var, dom_var,
# cross_var, cross_usefulness, mid_parent_bv, heterosis, pair_kinship, with a parent_kinship attribute.
ng_polyploid_score_crosses_dominance <- function(fit,
                                            dosage,
                                            pairs = NULL,
                                            selection_prop = 0.10,
                                            double_reduction = 0,
                                            grm_method = c("vanraden", "yang"),
                                            use_cpp = TRUE) {
  if (!inherits(fit, "ng_polyploid_effects")) ng_stop("fit must come from ng_polyploid_fit_effects")
  grm_method <- match.arg(grm_method)
  ploidy <- fit$ploidy
  M <- ng_polyploid_as_dosage_matrix(dosage, ploidy = ploidy, name = "dosage")
  miss <- setdiff(fit$markers, colnames(M))
  if (length(miss)) ng_stop("dosage is missing ", length(miss), " markers the model was fit on")
  M <- M[, fit$markers, drop = FALSE]
  ids <- rownames(M)
  storage.mode(M) <- "integer"

  mt <- ng_polyploid_progeny_moment_table(ploidy, double_reduction = double_reduction)
  ba <- fit$beta_add; bd <- fit$beta_dom; has_dom <- !is.null(bd)
  bd0 <- if (has_dom) bd else numeric(length(ba))
  cen_a <- ploidy * fit$allele_freq                     # additive centering
  hbar <- if (has_dom) fit$hbar else numeric(length(ba)) # dominance centering
  # The effects were fitted on the ORTHOGONAL dominance basis (R/48)
  #   D = (H - hbar) - b (X - ploidy p),   b = fit$b_orth (observed Cov(W,H)/Var(W); R/48)
  # so the progeny statistics must be taken on D, not on raw H. Everything needed follows from
  # the existing moment table:
  #   E[D]      = (E[H] - hbar) - b (E[X] - ploidy p)
  #   Var(D)    = Var(H) + b^2 Var(X) - 2 b Cov(X, H)
  #   Cov(X, D) = Cov(X, H) - b Var(X)
  # Scoring on raw H while the model was fitted on D would mis-state both the heterosis term and
  # the dominance variance whenever p != 0.5 -- i.e. at almost every marker.
  b_orth <- if (has_dom && !is.null(fit$b_orth)) as.numeric(fit$b_orth) else numeric(length(ba))
  intensity <- ng_selection_intensity(selection_prop)
  parent_kinship <- ng_polyploid_grm(M, ploidy = ploidy, method = grm_method)

  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = FALSE)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  p1 <- as.character(pairs$parent1); p2 <- as.character(pairs$parent2)
  unknown <- setdiff(unique(c(p1, p2)), ids)
  if (length(unknown)) ng_stop("pairs contains unknown parent IDs: ", paste(unknown, collapse = ", "))
  i1 <- match(p1, ids) - 1L; i2 <- match(p2, ids) - 1L  # 0-based parent indices

  if (isTRUE(use_cpp) && exists("ng_poly_dominance_scores_cpp", mode = "function", inherits = TRUE)) {
    # C++ accumulates the O(n_crosses x markers) loop over the moment table
    res <- ng_poly_dominance_scores_cpp(M, i1, i2, mt$mu, mt$varX, mt$EH, mt$varH, mt$covXH,
                                        ba, bd0, cen_a, hbar, b_orth, fit$intercept, has_dom)
    mid_bv <- res[, 1]; heterosis <- res[, 2]; add_var <- res[, 3]
    dom_var <- res[, 4]; cov_ad <- res[, 5]
  } else {
    n <- length(p1)
    mid_bv <- heterosis <- add_var <- dom_var <- cov_ad <- numeric(n)
    for (i in seq_len(n)) {
      di <- M[p1[i], ] + 1L; dj <- M[p2[i], ] + 1L      # 1-based index into the moment table
      ij <- cbind(di, dj)
      wx <- mt$mu[ij] - cen_a                           # E[X] - ploidy p
      mid_bv[i] <- fit$intercept + sum(ba * wx)
      add_var[i] <- sum(ba^2 * mt$varX[ij])
      if (has_dom) {
        ed  <- (mt$EH[ij] - hbar) - b_orth * wx                                  # E[D]
        vd  <- mt$varH[ij] + b_orth^2 * mt$varX[ij] - 2 * b_orth * mt$covXH[ij]  # Var(D)
        cxd <- mt$covXH[ij] - b_orth * mt$varX[ij]                               # Cov(X, D)
        heterosis[i] <- sum(bd * ed)                    # expected progeny dominance = heterosis
        dom_var[i] <- sum(bd^2 * pmax(vd, 0))
        cov_ad[i] <- 2 * sum(ba * bd * cxd)             # X and D share the draw
      }
    }
  }
  cross_mean <- mid_bv + heterosis
  cross_var <- pmax(add_var + dom_var + cov_ad, 0)
  out <- data.frame(parent1 = p1, parent2 = p2,
                    cross_mean = cross_mean, mid_parent_bv = mid_bv, heterosis = heterosis,
                    add_var = add_var, dom_var = dom_var, cross_var = cross_var,
                    cross_usefulness = cross_mean + intensity * sqrt(pmax(cross_var, 0)),
                    pair_kinship = ng_poly4x_pair_coancestry(parent_kinship, pairs),
                    stringsAsFactors = FALSE)
  attr(out, "parent_kinship") <- parent_kinship
  attr(out, "ploidy") <- as.integer(ploidy)
  attr(out, "has_dominance") <- !is.null(bd)
  # Travels with the numbers: an autopolyploid within-family variance is unbiased over unknown
  # parental phase but cannot resolve linkage-phase differences between crosses (see header).
  out$variance_model <- "unlinked_phase_marginalized"
  attr(out, "variance_model") <- "unlinked_phase_marginalized"
  attr(out, "dominance_model") <- if (is.null(bd)) NA_character_ else "digenic"
  attr(out, "double_reduction") <- as.numeric(double_reduction)
  out
}

# ---- EXACT phased autopolyploid within-family variance ---------------------------------------
#
# The dosage-only variance above is unbiased over unknown phase but cannot separate two crosses
# whose parents have identical dosages and different linkage phase. When the parental HOMOLOGUES
# are known (phased), the within-family additive variance is available exactly.
#
# Model: P homologues pair at random into P/2 bivalents; the gamete takes one chromatid per
# bivalent (P/2 alleles). Within a bivalent {i,j} the gamete carries h_i or h_j with probability
# 1/2 at a locus, and the SAME homologue at a second locus with probability (1 - r), so
#   Cov(g_k, g_l | bivalent{i,j}) = (1 - 2r_kl) (h_ik - h_jk)(h_il - h_jl) / 4.
# Averaging over the (P-1)!! pairings -- each unordered pair appears in (P-3)!! of them, i.e. a
# weight of 1/(P-1) on the sum over ALL pairs -- and using
#   sum_{i<j} (h_ik - h_jk)(h_il - h_jl) = P sum_i h_ik h_il - (sum_i h_ik)(sum_i h_il)
# gives, for marker effects b and the decay matrix R_kl = (1 - 2 r_kl),
#
#   Var(gamete value) = [ P * sum_i (b*h_i)' R (b*h_i)  -  (b*d)' R (b*d) ] / (4 (P - 1))
#
# with d = sum_i h_i the parental dosage. The cross variance is the sum over the two parents,
# because gametes from different parents are independent.
#
# Two properties worth knowing, both verified in tests/polyploid_phased_variance.R:
#   * under R = I (unlinked) it collapses EXACTLY to the hypergeometric d(P-d)/(4(P-1)), i.e. to
#     the dosage-only result -- phase changes the answer only through linkage, as it must;
#   * against directly simulated meiosis (random bivalent pairing + crossovers) it agrees to
#     Monte-Carlo error.
# Double reduction is NOT modelled on this path (random chromosome segregation only); use the
# dosage path with `double_reduction` if DR matters more than phase for your crop.
#
# Cost note: the per-parent term does not depend on the mate, so it is computed once per parent
# (P + 1 quadratic forms each) and every cross is then a sum of two precomputed numbers -- O(parents),
# not O(crosses).

# Split "<parent>_Hap<k>" rownames into parent id + homologue index. Mirrors the diploid
# "<parent>_HapA"/"_HapB" convention used by ng_gms_additive_var_general.
ng_poly_hap_parents <- function(haplotypes, ploidy) {
  rn <- rownames(haplotypes)
  if (is.null(rn)) ng_stop("phased_haplotypes must have rownames '<parent>_Hap<k>'")
  m <- regmatches(rn, regexec("^(.*)_Hap([0-9]+)$", rn))
  bad <- vapply(m, length, 0L) != 3L
  if (any(bad)) {
    ng_stop("phased_haplotypes rownames must look like '<parent>_Hap1'..'<parent>_Hap", ploidy,
            "'; offending: ", paste(utils::head(rn[bad], 3L), collapse = ", "))
  }
  parent <- vapply(m, `[`, character(1), 2L)
  idx <- as.integer(vapply(m, `[`, character(1), 3L))
  tab <- table(parent)
  wrong <- names(tab)[tab != ploidy]
  if (length(wrong)) {
    ng_stop("phased_haplotypes: ", length(wrong), " parent(s) do not have exactly ", ploidy,
            " homologues (e.g. ", paste(utils::head(wrong, 3L), collapse = ", "), ")")
  }
  if (any(idx < 1L | idx > ploidy)) ng_stop("phased_haplotypes homologue indices must be 1..", ploidy)
  list(parent = parent, index = idx)
}

# Per-parent gamete-value variance (the bracketed term above), summed over chromosomes.
# Returns a named numeric vector over parents.
ng_poly_phased_parent_var <- function(haplotypes, beta, marker_map, ploidy,
                                      recomb_model = c("haldane", "kosambi")) {
  recomb_model <- match.arg(recomb_model)
  P <- as.integer(ploidy)
  if (P < 2L || P %% 2L != 0L) ng_stop("phased polyploid variance needs an even ploidy >= 2")
  H <- as.matrix(haplotypes); storage.mode(H) <- "double"
  info <- ng_poly_hap_parents(H, P)
  mm <- ng_prepare_marker_map(marker_map, colnames(H), model = recomb_model)
  ord <- order(mm$chr_index, mm$pos_cm, mm$marker)
  mm <- mm[ord, , drop = FALSE]
  # Effects may arrive named by marker or positional in the haplotype column order. Resolve
  # explicitly and FAIL on an unmatched marker: silently coercing a missing effect to 0 yields a
  # plausible-looking variance (in the limit, exactly zero) instead of an error.
  if (is.null(names(beta))) {
    if (length(beta) != ncol(H)) {
      ng_stop("beta is unnamed, so it must have one value per haplotype column (",
              ncol(H), "); got ", length(beta))
    }
    beta <- stats::setNames(as.numeric(beta), colnames(H))
  }
  miss_b <- setdiff(mm$marker, names(beta))
  if (length(miss_b)) {
    ng_stop("beta is missing ", length(miss_b), " marker(s) present in the haplotypes, e.g. ",
            paste(utils::head(miss_b, 3L), collapse = ", "))
  }
  b <- as.numeric(beta[mm$marker])
  if (anyNA(b)) ng_stop("beta contains NA for ", sum(is.na(b)), " marker(s)")
  H <- H[, mm$marker, drop = FALSE]

  parents <- unique(info$parent)
  out <- stats::setNames(numeric(length(parents)), parents)
  for (cc in unique(mm$chr_index)) {
    idx <- which(mm$chr_index == cc)
    # R[k,l] = (1 - 2 r_kl) on this chromosome; cross-chromosome entries are 0 by construction,
    # so the quadratic form decomposes into independent per-chromosome blocks.
    Rc <- ng_recomb_decay_matrix(mm[idx, , drop = FALSE], model = recomb_model, target = "DH")
    bc <- b[idx]
    for (p in parents) {
      rows <- which(info$parent == p)
      Hp <- H[rows, idx, drop = FALSE]                       # P x m_chr homologues
      A <- sweep(Hp, 2L, bc, "*")                            # rows: b * h_i
      t1 <- sum(rowSums((A %*% Rc) * A))                     # sum_i (b h_i)' R (b h_i)
      bd <- bc * colSums(Hp)                                 # b * dosage
      t2 <- drop(bd %*% Rc %*% bd)
      out[[p]] <- out[[p]] + (P * t1 - t2)
    }
  }
  out / (4 * (P - 1))
}

# Exact within-family additive variance per cross, from phased parental homologues.
ng_poly_phased_within_family_var <- function(haplotypes, pairs, beta, marker_map, ploidy,
                                             recomb_model = c("haldane", "kosambi")) {
  recomb_model <- match.arg(recomb_model)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(pairs))) {
    ng_stop("pairs must contain parent1 and parent2")
  }
  pv <- ng_poly_phased_parent_var(haplotypes, beta, marker_map, ploidy,
                                  recomb_model = recomb_model)
  p1 <- as.character(pairs$parent1); p2 <- as.character(pairs$parent2)
  miss <- setdiff(unique(c(p1, p2)), names(pv))
  if (length(miss)) {
    ng_stop("phased_haplotypes is missing homologues for parent(s): ",
            paste(utils::head(miss, 5L), collapse = ", "))
  }
  unname(pv[p1] + pv[p2])
}
