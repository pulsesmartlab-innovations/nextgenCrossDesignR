# Additive + dominance marker-effect estimation (R/48). When dominance contributes to the trait,
# the genotypic-value prediction (additive + dominance) must track the TRUE total genetic value
# better than an additive-only breeding-value prediction -- the case that matters for clonal crops
# (cassava, potato) that select on genotypic value.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(48)
for (ploidy in c(2L, 4L)) {
  n <- 200L; m <- 120L
  ids <- sprintf("C%03d", seq_len(n))
  # allele freqs spread away from 0/1 so markers are informative
  freq <- runif(m, 0.2, 0.8)
  dosage <- matrix(rbinom(n * m, ploidy, rep(freq, each = n)), n, m,
                   dimnames = list(ids, sprintf("m%03d", seq_len(m))))
  p <- colMeans(dosage) / ploidy
  # true additive + dominance genetic values
  add_eff <- rnorm(m, 0, 1)
  dom_eff <- rnorm(m, 0, 1.2)                       # substantial dominance
  W <- sweep(dosage, 2, ploidy * p); H <- dosage * (ploidy - dosage); D <- sweep(H, 2, colMeans(H))
  true_bv <- as.numeric(W %*% add_eff)
  true_gv <- true_bv + as.numeric(D %*% dom_eff)    # total genotypic value
  pheno <- true_gv + rnorm(n, 0, sd(true_gv) * 0.5) # h2 ~ 0.8 on genotypic value
  names(pheno) <- ids

  train <- 1:150; test <- 151:200
  # additive-only vs additive+dominance
  fit_a <- ng_fit_polyploid_effects(dosage[train, ], pheno[train], ploidy = ploidy, model = "additive")
  fit_ad <- ng_fit_polyploid_effects(dosage[train, ], pheno[train], ploidy = ploidy, model = "additive_dominance")
  stopifnot(is.null(fit_a$beta_dom), !is.null(fit_ad$beta_dom))

  # predict genotypic value on held-out clones
  gv_a <- ng_predict_polyploid_value(fit_a, dosage[test, ], type = "genotypic")   # falls back to BV
  gv_ad <- ng_predict_polyploid_value(fit_ad, dosage[test, ], type = "genotypic")
  bv_ad <- ng_predict_polyploid_value(fit_ad, dosage[test, ], type = "breeding")

  acc_a <- cor(gv_a, true_gv[test])
  acc_ad <- cor(gv_ad, true_gv[test])
  # modelling dominance improves prediction of TOTAL genotypic value
  stopifnot(acc_ad > acc_a)
  # breeding value != genotypic value when dominance is present
  stopifnot(cor(bv_ad, gv_ad) < 0.999, sd(gv_ad - bv_ad) > 1e-6)
  cat(sprintf("ploidy %d: genotypic-value accuracy add-only %.3f -> add+dom %.3f\n", ploidy, acc_a, acc_ad))
}

cat("polyploid additive+dominance effects test passed\n")
