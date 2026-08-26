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
  fit_a <- ng_polyploid_fit_effects(dosage[train, ], pheno[train], ploidy = ploidy, model = "additive")
  fit_ad <- ng_polyploid_fit_effects(
    dosage[train, ], pheno[train], ploidy = ploidy,
    model = "additive_dominance", allow_experimental_dominance = TRUE
  )
  stopifnot(is.null(fit_a$beta_dom), !is.null(fit_ad$beta_dom))

  # predict genotypic value on held-out clones
  gv_a <- ng_polyploid_predict_value(fit_a, dosage[test, ], type = "genotypic")   # falls back to BV
  gv_ad <- ng_polyploid_predict_value(fit_ad, dosage[test, ], type = "genotypic")
  bv_ad <- ng_polyploid_predict_value(fit_ad, dosage[test, ], type = "breeding")

  acc_a <- cor(gv_a, true_gv[test])
  acc_ad <- cor(gv_ad, true_gv[test])
  # modelling dominance improves prediction of TOTAL genotypic value
  stopifnot(acc_ad > acc_a)
  # breeding value != genotypic value when dominance is present
  stopifnot(cor(bv_ad, gv_ad) < 0.999, sd(gv_ad - bv_ad) > 1e-6)
  cat(sprintf("ploidy %d: genotypic-value accuracy add-only %.3f -> add+dom %.3f\n", ploidy, acc_a, acc_ad))
}

cat("polyploid additive+dominance effects test passed\n")

# --- orthogonal dominance parameterization (0.19.x) -----------------------------------------
# H = d(ploidy-d) mean-centred is NOT orthogonal to the additive design: under HWE
# Cov(M,H) = ploidy(ploidy-1) p q (1-2p), zero only at p = 0.5. At low MAF the two columns reach
# r ~ 0.96, so ridge cannot separate additive from dominance and the split is decided by the
# penalty. The fitted design regresses H on W (observed coefficient) to remove exactly that.
set.seed(31)
Po <- 4L; no <- 60L; mo <- 50L
ido <- sprintf("O%02d", seq_len(no))
fo <- runif(mo, 0.05, 0.95)                       # skewed frequencies: the collinear regime
Mo <- matrix(rbinom(no * mo, Po, rep(fo, each = no)), no, mo,
             dimnames = list(ido, sprintf("W%03d", seq_len(mo))))
yo <- as.numeric(Mo %*% rnorm(mo, 0, .2) + (Mo * (Po - Mo)) %*% rnorm(mo, 0, .1)) + rnorm(no, 0, .5)
fo_fit <- ng_polyploid_fit_effects(
  Mo, yo, ploidy = Po, model = "additive_dominance", seed = 2L,
  allow_experimental_dominance = TRUE
)

stopifnot(!is.null(fo_fit$b_orth), length(fo_fit$b_orth) == mo, all(is.finite(fo_fit$b_orth)))
Wo <- sweep(Mo, 2L, Po * fo_fit$allele_freq, "-")
Ho <- Mo * (Po - Mo)
Do_plain <- sweep(Ho, 2L, colMeans(Ho), "-")                       # the old, collinear basis
Do_orth  <- Do_plain - sweep(Wo, 2L, fo_fit$b_orth, "*")           # the fitted basis
cor_plain <- vapply(seq_len(mo), function(k) suppressWarnings(stats::cor(Wo[, k], Do_plain[, k])), numeric(1))
cor_orth  <- vapply(seq_len(mo), function(k) suppressWarnings(stats::cor(Wo[, k], Do_orth[, k])), numeric(1))
# the problem is real...
stopifnot(max(abs(cor_plain), na.rm = TRUE) > 0.8)
# ...and the fitted basis is orthogonal to numerical precision, not merely "less correlated"
stopifnot(max(abs(cor_orth), na.rm = TRUE) < 1e-8)

# Under HWE the observed coefficient converges to the theoretical (ploidy-1)(1-2p); at ploidy 2
# that limit is the standard orthogonal dominance coding of Vitezica et al. (2013). Checked on the
# design directly -- fitting 20k samples would build a 20k x 20k dual and is beside the point.
set.seed(32)
for (Ph in c(2L, 4L)) for (ph in c(0.2, 0.5, 0.75)) {
  Mh <- rbinom(60000, Ph, ph)
  Wh <- Mh - Ph * (mean(Mh) / Ph)
  Hh <- Mh * (Ph - Mh)
  bh <- sum((Wh - mean(Wh)) * (Hh - mean(Hh))) / sum((Wh - mean(Wh))^2)
  stopifnot(abs(bh - (Ph - 1) * (1 - 2 * ph)) < 0.05)
}
# ploidy 2 HWE limit == published Vitezica orthogonal dominance coding (-2p^2, 2pq, -2q^2)
for (pv in c(0.1, 0.35, 0.6)) {
  qv <- 1 - pv; Mv <- 0:2; Hv <- Mv * (2 - Mv)
  Dv <- Hv - sum(dbinom(Mv, 2, pv) * Hv) - (2 - 1) * (1 - 2 * pv) * (Mv - 2 * pv)
  stopifnot(max(abs(Dv - c(-2 * pv^2, 2 * pv * qv, -2 * qv^2))) < 1e-12)
}

# predict() must use the SAME basis it fitted on: genotypic value = breeding value + D beta_dom
gvo <- ng_polyploid_predict_value(fo_fit, Mo, type = "genotypic")
bvo <- ng_polyploid_predict_value(fo_fit, Mo, type = "breeding")
stopifnot(max(abs((gvo - bvo) - as.numeric(Do_orth %*% fo_fit$beta_dom))) < 1e-9)

# additive-only fits carry no dominance machinery at all
fa_o <- ng_polyploid_fit_effects(Mo, yo, ploidy = Po, model = "additive", seed = 2L)
stopifnot(is.null(fa_o$beta_dom), is.null(fa_o$b_orth), is.null(fa_o$hbar))
cat("polyploid orthogonal dominance parameterization test passed\n")
