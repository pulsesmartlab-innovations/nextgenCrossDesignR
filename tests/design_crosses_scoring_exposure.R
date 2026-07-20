# The documented pipeline entry point ng_design_crosses() must expose the two scientific
# scoring choices that ng_score_crosses() supports and ng_run_cross_prediction() already
# surfaces: the breeding target (DH vs RIL) and the map function (Haldane vs Kosambi).
# Both drive the recombination-aware within-family variance (vpm); defaults must
# stay DH / Haldane (backward compatible).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(414)
n <- 24L; m <- 150L
ids <- sprintf("P%02d", seq_len(n))
geno <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
               dimnames = list(ids, sprintf("M%03d", seq_len(m))))
beta <- rnorm(m, 0, 0.1)
y <- as.numeric(geno %*% beta) + rnorm(n); names(y) <- ids
# spread markers along 5 chromosomes so recombination distance (hence DH/RIL and
# Haldane/Kosambi) genuinely affects the within-family variance.
mm <- data.frame(marker = colnames(geno), chr = rep(1:5, length.out = m),
                 pos_cm = rep(seq(0, 120, length.out = ceiling(m / 5)), 5)[seq_len(m)],
                 stringsAsFactors = FALSE)
base <- list(geno = geno, y = y, marker_map = mm, ids = ids, n_crosses = 8L, use_cpp = FALSE)

d0    <- do.call(ng_design_crosses, base)
d_dh  <- do.call(ng_design_crosses, c(base, list(target = "DH")))
d_ril <- do.call(ng_design_crosses, c(base, list(target = "RIL")))
d_hal <- do.call(ng_design_crosses, c(base, list(recomb_model = "haldane")))
d_kos <- do.call(ng_design_crosses, c(base, list(recomb_model = "kosambi")))

# default is DH + Haldane (backward compatible)
stopifnot(isTRUE(all.equal(d0$scores$vpm, d_dh$scores$vpm)))
stopifnot(isTRUE(all.equal(d0$scores$vpm, d_hal$scores$vpm)))
# target and recomb_model reach the scorer and change the recombination variance
stopifnot(!isTRUE(all.equal(d_dh$scores$vpm, d_ril$scores$vpm)))
stopifnot(!isTRUE(all.equal(d_hal$scores$vpm, d_kos$scores$vpm)))

cat("ng_design_crosses scoring-exposure test passed\n")
