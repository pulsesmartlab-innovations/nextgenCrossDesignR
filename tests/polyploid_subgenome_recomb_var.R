## Recombination-aware within-family variance for the disomic-subgenome path.
##
## Proves the quantitative-genetics invariants of ng_polyploid_subgenome_score_crosses:
##   1. Backward compatibility: no map => linkage-equilibrium formula, unchanged.
##   2. Closed form: for two coupling/repulsion markers d cM apart the kernel returns
##      b1^2 + b2^2 +/- 2 b1 b2 exp(-2 d / 100)  (Haldane DH), matching theory exactly.
##   3. R = I reduction: markers on separate chromosomes => recomb-aware == LE.
##   4. Block-diagonal exactness: total recomb-aware variance == sum of the diploid
##      kernel applied independently per subgenome (disomic homoeologues never pair).
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

## ---- 1. Backward compatibility: no map => linkage_equilibrium, value unchanged ----
geno_by_subgenome <- list(
  A = matrix(c(2, 0, 2, 0,
               0, 2, 0, 2,
               2, 2, 0, 0), nrow = 3, byrow = TRUE,
             dimnames = list(c("P1", "P2", "P3"), c("A_m1", "A_m2", "A_m3", "A_m4"))),
  B = matrix(c(0, 2, 0, 2,
               2, 0, 2, 0,
               2, 0, 0, 2), nrow = 3, byrow = TRUE,
             dimnames = list(c("P1", "P2", "P3"), c("B_m1", "B_m2", "B_m3", "B_m4")))
)
effects_by_subgenome <- list(
  A = c(A_m1 = 0.5, A_m2 = -0.25, A_m3 = 0.4, A_m4 = 0.1),
  B = c(B_m1 = 0.25, B_m2 = 0.75, B_m3 = -0.3, B_m4 = 0.2)
)
pairs <- data.frame(parent1 = c("P1", "P1", "P2"),
                    parent2 = c("P2", "P3", "P3"), stringsAsFactors = FALSE)

le <- ng_polyploid_subgenome_score_crosses(geno_by_subgenome, effects_by_subgenome,
                                           candidate_pairs = pairs)
stopifnot(identical(unique(le$poly_variance_model), "linkage_equilibrium"))
stopifnot(identical(attr(le, "variance_model"), "linkage_equilibrium"))
# hand LE variance for cross 1 (P1 vs P2)
le_manual <- 0
for (sg in names(geno_by_subgenome)) {
  d <- geno_by_subgenome[[sg]]["P1", ] - geno_by_subgenome[[sg]]["P2", ]
  le_manual <- le_manual + sum((effects_by_subgenome[[sg]] * d)^2) / 4
}
stopifnot(isTRUE(all.equal(le$poly_var[[1]], le_manual, tolerance = 1e-12)))

## ---- 2. Closed form: coupling vs repulsion on one chromosome, two markers 5 cM apart ----
d_cm <- 5
decay <- exp(-2 * d_cm / 100)
# coupling: P1 = (2,2), P2 = (0,0); repulsion: P1 = (2,0), P2 = (0,2)
g_coupling <- list(S = matrix(c(2, 2, 0, 0), 2, 2, byrow = TRUE,
                              dimnames = list(c("P1", "P2"), c("m1", "m2"))))
g_repulsion <- list(S = matrix(c(2, 0, 0, 2), 2, 2, byrow = TRUE,
                               dimnames = list(c("P1", "P2"), c("m1", "m2"))))
eff <- list(S = c(m1 = 1, m2 = 1))
map1 <- list(S = data.frame(marker = c("m1", "m2"), chr = "1", pos_cm = c(0, d_cm),
                            stringsAsFactors = FALSE))
p1 <- data.frame(parent1 = "P1", parent2 = "P2", stringsAsFactors = FALSE)

sc_c <- ng_polyploid_subgenome_score_crosses(g_coupling, eff, candidate_pairs = p1,
                                             map_by_subgenome = map1)
sc_r <- ng_polyploid_subgenome_score_crosses(g_repulsion, eff, candidate_pairs = p1,
                                             map_by_subgenome = map1)
stopifnot(identical(unique(sc_c$poly_variance_model), "recombination_aware"))
theory_coupling  <- 1 + 1 + 2 * 1 * 1 * decay
theory_repulsion <- 1 + 1 - 2 * 1 * 1 * decay
stopifnot(isTRUE(all.equal(sc_c$poly_var[[1]], theory_coupling,  tolerance = 1e-8)))
stopifnot(isTRUE(all.equal(sc_r$poly_var[[1]], theory_repulsion, tolerance = 1e-8)))
# coupling raises variance above LE (=2), repulsion lowers it: linkage is signed
stopifnot(sc_c$poly_var[[1]] > 2, sc_r$poly_var[[1]] < 2)

## ---- 3. R = I reduction: each marker on its own chromosome => recomb-aware == LE ----
map_sep <- list(S = data.frame(marker = c("m1", "m2"), chr = c("1", "2"),
                               pos_cm = c(0, 0), stringsAsFactors = FALSE))
sc_sep <- ng_polyploid_subgenome_score_crosses(g_coupling, eff, candidate_pairs = p1,
                                               map_by_subgenome = map_sep)
stopifnot(isTRUE(all.equal(sc_sep$poly_var[[1]], 2, tolerance = 1e-10)))

## ---- 4. Block-diagonal exactness: total == sum of per-subgenome kernel calls ----
map_by_subgenome <- list(
  A = data.frame(marker = colnames(geno_by_subgenome$A), chr = c("A1", "A1", "A2", "A2"),
                 pos_cm = c(0, 8, 0, 12), stringsAsFactors = FALSE),
  B = data.frame(marker = colnames(geno_by_subgenome$B), chr = c("B1", "B1", "B1", "B2"),
                 pos_cm = c(0, 5, 15, 0), stringsAsFactors = FALSE)
)
sc_ra <- ng_polyploid_subgenome_score_crosses(geno_by_subgenome, effects_by_subgenome,
                                              candidate_pairs = pairs,
                                              map_by_subgenome = map_by_subgenome)
per_sg_sum <- numeric(nrow(pairs))
for (sg in names(geno_by_subgenome)) {
  mp <- ng_prepare_marker_map(map_by_subgenome[[sg]], colnames(geno_by_subgenome[[sg]]))
  v <- ng_dh_recomb_variance_pairs(geno = geno_by_subgenome[[sg]],
         beta = effects_by_subgenome[[sg]], beta_var = rep(0, ncol(geno_by_subgenome[[sg]])),
         marker_map = mp, ids = rownames(geno_by_subgenome[[sg]]), pairs = pairs, target = "DH")
  per_sg_sum <- per_sg_sum + as.numeric(v$vpm)
}
stopifnot(isTRUE(all.equal(sc_ra$poly_var, per_sg_sum, tolerance = 1e-10)))
# mean (gain) is variance-model-independent
stopifnot(isTRUE(all.equal(sc_ra$poly_gain, le$poly_gain, tolerance = 1e-12)))

## ---- validation: map names must match subgenomes ----
expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = function(e) e)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), ignore.case = TRUE))
}
expect_error(
  ng_polyploid_subgenome_score_crosses(geno_by_subgenome, effects_by_subgenome,
    candidate_pairs = pairs, map_by_subgenome = list(A = map_by_subgenome$A)),
  "map_by_subgenome names"
)

cat("polyploid_subgenome_recomb_var.R: PASS\n")
