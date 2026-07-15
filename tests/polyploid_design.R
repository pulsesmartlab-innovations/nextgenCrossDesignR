# Any-ploidy polyploid mate design (R/46): ng_poly_score_crosses (analytic mean + kinship) and
# ng_design_crosses_poly (QC -> score -> native allocator with the full control suite).
# AlphaSimR-free: synthetic dosage 0..ploidy, so it is fast and deterministic.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(46)
for (ploidy in c(2L, 4L)) {
  np <- 16L; nm <- 30L
  ids <- sprintf("P%02d", seq_len(np))
  dosage <- matrix(rbinom(np * nm, ploidy, 0.5), np, nm, dimnames = list(ids, sprintf("M%02d", seq_len(nm))))
  effects <- rnorm(nm)

  # --- scorer contract ---
  scores <- ng_poly_score_crosses(dosage, effects, ploidy = ploidy)
  stopifnot(all(c("parent1", "parent2", "poly_mean", "pair_kinship") %in% names(scores)))
  stopifnot(nrow(scores) == choose(np, 2L))
  K <- attr(scores, "parent_K"); gebv <- attr(scores, "parent_gebv")
  stopifnot(is.matrix(K), nrow(K) == np, identical(attr(scores, "ploidy"), ploidy))
  # poly_mean equals the mid-parent GEBV exactly
  i1 <- scores$parent1[[1]]; i2 <- scores$parent2[[1]]
  stopifnot(abs(scores$poly_mean[[1]] - (gebv[[i1]] + gebv[[i2]]) / 2) < 1e-8)

  # --- one-call design with explicit effects ---
  plan <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, effects = effects,
                                 method = "greedy_local")
  stopifnot(nrow(plan) == 8L, identical(attr(plan, "ploidy"), ploidy))

  # --- effects estimated from phenotype ---
  pheno <- as.numeric(gebv) + rnorm(np, 0, 1); names(pheno) <- ids
  plan_ph <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, phenotype = pheno,
                                    method = "greedy_local")
  stopifnot(nrow(plan_ph) == 8L)

  # --- a forwarded breeder control (committed matings) is honoured (proves the reuse works) ---
  committed <- data.frame(parent1 = "P01", parent2 = "P02", stringsAsFactors = FALSE)
  plan_c <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, effects = effects,
                                   committed_crosses = committed, method = "greedy_local")
  key <- ng_group_pair_key(plan_c$parent1, plan_c$parent2)
  stopifnot(ng_group_pair_key("P01", "P02") %in% key)

  # --- a forwarded strategy dial is honoured ---
  plan_div <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, effects = effects,
                                     strategy = "diversity", method = "greedy_local")
  stopifnot(nrow(plan_div) == 8L, !is.null(attr(plan_div, "summary")$achieved_emphasis))

  # --- ploidy-aware QC runs by default and attaches a summary ---
  plan_qc <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, effects = effects,
                                    run_qc = TRUE, method = "greedy_local")
  stopifnot(nrow(plan_qc) == 8L, !is.null(attr(plan_qc, "qc")),
            attr(plan_qc, "qc")$ploidy == ploidy)

  # --- ADDITIVE path exposes within-family variance + usefulness (gain selectable) ---
  sc_add <- ng_poly_score_crosses(dosage, effects, ploidy = ploidy)
  stopifnot(all(c("poly_mean", "poly_var", "poly_usefulness") %in% names(sc_add)),
            all(sc_add$poly_var >= 0))
  plan_u <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, effects = effects,
                                   gain = "usefulness", method = "greedy_local")
  stopifnot(identical(attr(plan_u, "summary")$poly_gain_col, "poly_usefulness"))

  # --- optional dominance path (needs a phenotype) scores on genotypic value + heterosis ---
  pheno <- as.numeric(dosage %*% effects) + rnorm(np); names(pheno) <- ids
  plan_dom <- ng_design_crosses_poly(dosage, n_crosses = 8L, ploidy = ploidy, phenotype = pheno,
                                     dominance = TRUE, gain = "usefulness", method = "greedy_local")
  stopifnot(nrow(plan_dom) == 8L, isTRUE(attr(plan_dom, "summary")$dominance),
            "cross_usefulness" %in% names(plan_dom), "heterosis" %in% names(plan_dom))
  cat(sprintf("ploidy %d: score + design + controls + QC + dominance OK\n", ploidy))
}

cat("polyploid design (any-ploidy) test passed\n")
