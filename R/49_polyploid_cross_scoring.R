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
# precomputed once and looked up vectorized -- no per-cross simulation. (Linkage-equilibrium /
# unlinked-loci approximation for the segregation variance; digenic dominance is per-locus. Double
# reduction is not modelled here -- random chromosome segregation only.)

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
ng_poly_progeny_moment_table <- function(ploidy, double_reduction = 0) {
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
# from ng_fit_polyploid_effects. Returns parent1/parent2 + cross_mean (genotypic), add_var, dom_var,
# cross_var, cross_usefulness, mid_parent_bv, heterosis, pair_kinship, with a parent_K attribute.
ng_score_crosses_poly_dominance <- function(fit,
                                            dosage,
                                            pairs = NULL,
                                            selection_prop = 0.10,
                                            double_reduction = 0,
                                            grm_method = c("vanraden", "yang"),
                                            use_cpp = TRUE) {
  if (!inherits(fit, "ng_polyploid_effects")) ng_stop("fit must come from ng_fit_polyploid_effects")
  grm_method <- match.arg(grm_method)
  ploidy <- fit$ploidy
  M <- ng_poly4x_as_dosage_matrix(dosage, ploidy = ploidy, name = "dosage")
  miss <- setdiff(fit$markers, colnames(M))
  if (length(miss)) ng_stop("dosage is missing ", length(miss), " markers the model was fit on")
  M <- M[, fit$markers, drop = FALSE]
  ids <- rownames(M)
  storage.mode(M) <- "integer"

  mt <- ng_poly_progeny_moment_table(ploidy, double_reduction = double_reduction)
  ba <- fit$beta_add; bd <- fit$beta_dom; has_dom <- !is.null(bd)
  bd0 <- if (has_dom) bd else numeric(length(ba))
  cen_a <- ploidy * fit$allele_freq                     # additive centering
  hbar <- if (has_dom) fit$hbar else numeric(length(ba)) # dominance centering
  intensity <- ng_selection_intensity(selection_prop)
  parent_K <- ng_polyploid_grm(M, ploidy = ploidy, method = grm_method)

  if (is.null(pairs)) pairs <- ng_make_pairs(ids, include_self = FALSE)
  pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
  p1 <- as.character(pairs$parent1); p2 <- as.character(pairs$parent2)
  unknown <- setdiff(unique(c(p1, p2)), ids)
  if (length(unknown)) ng_stop("pairs contains unknown parent IDs: ", paste(unknown, collapse = ", "))
  i1 <- match(p1, ids) - 1L; i2 <- match(p2, ids) - 1L  # 0-based parent indices

  if (isTRUE(use_cpp) && exists("ng_poly_dominance_scores_cpp", mode = "function", inherits = TRUE)) {
    # C++ accumulates the O(n_crosses x markers) loop over the moment table
    res <- ng_poly_dominance_scores_cpp(M, i1, i2, mt$mu, mt$varX, mt$EH, mt$varH, mt$covXH,
                                        ba, bd0, cen_a, hbar, fit$intercept, has_dom)
    mid_bv <- res[, 1]; heterosis <- res[, 2]; add_var <- res[, 3]
    dom_var <- res[, 4]; cov_ad <- res[, 5]
  } else {
    n <- length(p1)
    mid_bv <- heterosis <- add_var <- dom_var <- cov_ad <- numeric(n)
    for (i in seq_len(n)) {
      di <- M[p1[i], ] + 1L; dj <- M[p2[i], ] + 1L      # 1-based index into the moment table
      ij <- cbind(di, dj)
      mid_bv[i] <- fit$intercept + sum(ba * (mt$mu[ij] - cen_a))
      add_var[i] <- sum(ba^2 * mt$varX[ij])
      if (has_dom) {
        heterosis[i] <- sum(bd * (mt$EH[ij] - hbar))    # expected progeny dominance = heterosis
        dom_var[i] <- sum(bd^2 * mt$varH[ij])
        cov_ad[i] <- 2 * sum(ba * bd * mt$covXH[ij])    # X and H share the draw
      }
    }
  }
  cross_mean <- mid_bv + heterosis
  cross_var <- pmax(add_var + dom_var + cov_ad, 0)
  out <- data.frame(parent1 = p1, parent2 = p2,
                    cross_mean = cross_mean, mid_parent_bv = mid_bv, heterosis = heterosis,
                    add_var = add_var, dom_var = dom_var, cross_var = cross_var,
                    cross_usefulness = cross_mean + intensity * sqrt(pmax(cross_var, 0)),
                    pair_kinship = ng_poly4x_pair_coancestry(parent_K, pairs),
                    stringsAsFactors = FALSE)
  attr(out, "parent_K") <- parent_K
  attr(out, "ploidy") <- as.integer(ploidy)
  attr(out, "has_dominance") <- !is.null(bd)
  out
}
