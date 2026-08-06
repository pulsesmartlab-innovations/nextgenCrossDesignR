# Exact within-cross additive variance for arbitrary phased parents
# (ng_gms_additive_var_general): generalizes a'Ra to residual-heterozygous
# (RIL) parents. See docs/design/residual-het-parent-variance.md.
#   DH  : Var = a'Ra + 1/2 b'[(1-r) o (1-2r) o (d1 d1' + d2 d2')] b
#   RIL : Var = a'R*a + 1/2 b'[R* o (d1 d1' + d2 d2')] b,  R* = (1-2r)/(1+2r)
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

hm <- function(a, b, cc, d, mk) {
  M <- rbind(P1_HapA = a, P1_HapB = b, P2_HapA = cc, P2_HapB = d); colnames(M) <- mk; M
}

## --- deterministic closed-form checks ------------------------------------
mm1 <- data.frame(marker = "m1", chr = 1L, chr_index = 1L, pos_cm = 0, stringsAsFactors = FALSE)
Rdh1 <- ng_recomb_decay_matrix(mm1, "haldane", "DH"); Rril1 <- ng_recomb_decay_matrix(mm1, "haldane", "RIL")
b1 <- c(m1 = 1)
# P1 = Aa (het) x P2 = aa -> 3/4 for both DH and RIL (single locus, r=0)
vhet_dh  <- ng_gms_additive_var_general("P1", "P2", hm(1, 0, 0, 0, "m1"), b1, Rdh1,  target = "DH")[["VPM"]]
vhet_ril <- ng_gms_additive_var_general("P1", "P2", hm(1, 0, 0, 0, "m1"), b1, Rril1, target = "RIL")[["VPM"]]
stopifnot(isTRUE(all.equal(vhet_dh, 0.75)), isTRUE(all.equal(vhet_ril, 0.75)))
# P1 = AA x P2 = aa (fully differing, inbred) -> 1
stopifnot(isTRUE(all.equal(ng_gms_additive_var_general("P1", "P2", hm(1, 1, 0, 0, "m1"), b1, Rdh1)[["VPM"]], 1)))
# monomorphic -> 0
stopifnot(isTRUE(all.equal(ng_gms_additive_var_general("P1", "P2", hm(0, 0, 0, 0, "m1"), b1, Rdh1)[["VPM"]], 0)))

## --- reduces exactly to a'Ra for inbred parents --------------------------
mm <- data.frame(marker = paste0("m", 1:5), chr = 1L, chr_index = 1L,
                 pos_cm = c(0, 15, 35, 60, 90), stringsAsFactors = FALSE)
beta <- c(m1 = 1, m2 = -0.7, m3 = 0.5, m4 = 0.8, m5 = -0.4)
Rdh <- ng_recomb_decay_matrix(mm, "haldane", "DH")
inbredH <- hm(c(0,0,0,0,0), c(0,0,0,0,0), c(1,1,1,1,1), c(1,1,1,1,1), mm$marker)
avec <- 0.5 * (c(0,0,0,0,0) - c(2,2,2,2,2)) * beta
ara <- as.numeric(t(avec) %*% Rdh %*% avec)
stopifnot(isTRUE(all.equal(ng_gms_additive_var_general("P1", "P2", inbredH, beta, Rdh)[["VPM"]], ara)))

## --- phase matters: coupling != repulsion (signed d_k d_l) ---------------
coup <- ng_gms_additive_var_general("P1", "P2", hm(c(1,1,0), c(0,0,0), c(0,0,1), c(0,0,1), paste0("m",1:3)),
          c(m1=1,m2=0.8,m3=-0.5), ng_recomb_decay_matrix(mm[1:3,], "haldane", "DH"))[["VPM"]]
repu <- ng_gms_additive_var_general("P1", "P2", hm(c(1,0,0), c(0,1,0), c(0,0,1), c(0,0,1), paste0("m",1:3)),
          c(m1=1,m2=0.8,m3=-0.5), ng_recomb_decay_matrix(mm[1:3,], "haldane", "DH"))[["VPM"]]
stopifnot(abs(coup - repu) > 0.2)                 # phase genuinely changes the variance

## --- PMV >= VPM with a PSD posterior; nonnegativity -----------------------
Sig <- diag(0.05, 5); dimnames(Sig) <- list(mm$marker, mm$marker)
mv <- ng_gms_additive_var_general("P1", "P2", hm(c(1,1,0,1,0), c(1,0,0,0,0), c(0,0,1,0,1), c(0,0,1,0,1), mm$marker),
        beta, Rdh, beta_cov = Sig)
stopifnot(mv[["VPM"]] >= 0, mv[["PMV"]] >= mv[["VPM"]] - 1e-9)

## --- Monte-Carlo guard (pure Haldane transmission) -----------------------
gam <- function(H, pos) { m <- ncol(H); g <- numeric(m); h <- sample.int(2L, 1L); g[1] <- H[h, 1]
  for (j in 2:m) { r <- 0.5 * (1 - exp(-2 * (pos[j] - pos[j-1]) / 100)); if (runif(1) < r) h <- 3L - h; g[j] <- H[h, j] }; g }
H1 <- rbind(c(1,1,0,1,0), c(1,0,0,0,0)); H2 <- rbind(c(0,0,1,0,1), c(0,0,1,0,1))
hmv <- hm(H1[1,], H1[2,], H2[1,], H2[2,], mm$marker)
set.seed(1)
vdh <- replicate(15000L, { f1 <- rbind(gam(H1, mm$pos_cm), gam(H2, mm$pos_cm)); sum(beta * (2 * gam(f1, mm$pos_cm))) })
mc_dh <- var(vdh)
cf_dh <- ng_gms_additive_var_general("P1", "P2", hmv, beta, Rdh, target = "DH")[["VPM"]]
stopifnot(abs(cf_dh - mc_dh) / mc_dh < 0.04)      # DH closed form within MC tolerance

cat(sprintf("gms_general_var: single-locus closed forms, inbred==a'Ra, phase-dependence, PMV, and DH MC (%.3f vs %.3f) passed\n",
            cf_dh, mc_dh))
