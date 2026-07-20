# window_cm caps the cM distance over which markers contribute to the within-family recombination
# variance (a'Ra). It is honored on the dense/banded (Kosambi, RIL) paths but was silently IGNORED
# on the fast Haldane-DH closed-form recursion, which sums the whole chromosome and cannot apply a
# cutoff. Fix: use the recursion only for window_cm = Inf; route finite windows to the dense/banded
# path (which honors the window). This test also checks the dense path matches the recursion at a
# window wide enough to cover every chromosome (consistency).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(99)
n <- 16L; m <- 90L; ids <- sprintf("P%02d", seq_len(n)); mk <- sprintf("S%03d", seq_len(m))
G <- matrix(2L * rbinom(n * m, 1, 0.5), n, m, dimnames = list(ids, mk)); storage.mode(G) <- "double"
y <- setNames(as.numeric(G %*% rnorm(m, 0, 0.2)) + rnorm(n), ids)
# 3 chromosomes, 30 markers each spread over 0..60 cM (2 cM spacing)
map <- data.frame(marker = mk, chr = rep(1:3, each = m / 3),
                  pos_cm = rep(seq(0, 60, length.out = m / 3), 3))
eff <- ng_fit_ridge_effects(G, y, ids)
pairs <- as.data.frame(t(utils::combn(ids, 2)), stringsAsFactors = FALSE)
names(pairs) <- c("parent1", "parent2")

score <- function(w) ng_score_crosses(G, eff, map, ids, pairs, target = "DH",
                                       recomb_model = "haldane", assume_inbred = TRUE,
                                       use_cpp = FALSE, window_cm = w)$vpm

v_inf   <- score(Inf)     # closed-form chromosome recursion
v_full  <- score(1000)    # finite window wider than any chromosome -> should equal the recursion
v_small <- score(5)       # real 5 cM window -> must change the variance

# A finite window that covers every chromosome must reproduce the full-chromosome recursion.
stopifnot(isTRUE(all.equal(v_inf, v_full, tolerance = 1e-6)))
# A real cutoff must actually change the recombination variance (previously it was ignored).
stopifnot(!isTRUE(all.equal(v_inf, v_small, tolerance = 1e-8)))
# Sanity: variances stay finite and non-negative.
stopifnot(all(is.finite(v_small)), all(v_small >= -1e-8))

cat("scoring_window_cm.R: PASS  (Inf vs full-window match; 5cM window changes variance)\n")
