# ng_add_p_superior_progeny_multitrait(cross_trait_cov = ...) must build Sigma_c from the EXACT
# within-family cross-trait covariance (matched by parent pair) and produce valid probabilities
# that differ from the population-correlation proxy.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(303)
n <- 8L; m <- 60L; nchr <- 3L
ids <- sprintf("P%02d", seq_len(n)); mk <- sprintf("M%03d", seq_len(m))
geno <- matrix(2L * rbinom(n * m, 1, 0.5), n, m, dimnames = list(ids, mk)); storage.mode(geno) <- "double"
map <- data.frame(marker = mk, chr = rep(seq_len(nchr), each = m / nchr),
                  pos_cm = rep(seq(0, 100, length.out = m / nchr), nchr), stringsAsFactors = FALSE)
b1 <- rnorm(m); b2 <- 0.5 * b1 + rnorm(m)
betas <- cbind(t1 = b1, t2 = b2); rownames(betas) <- mk

pairs <- as.data.frame(t(utils::combn(ids, 2)), stringsAsFactors = FALSE); names(pairs) <- c("parent1", "parent2")
ctc <- ng_cross_trait_within_family_cov(geno, betas, map, pairs = pairs, target = "DH")

# Cross means: mid-parent GEBV per trait.
gebv1 <- as.numeric(geno %*% b1); gebv2 <- as.numeric(geno %*% b2); names(gebv1) <- names(gebv2) <- ids
p1 <- pairs$parent1; p2 <- pairs$parent2
scores <- data.frame(parent1 = p1, parent2 = p2,
                     m_t1 = 0.5 * (gebv1[p1] + gebv1[p2]),
                     m_t2 = 0.5 * (gebv2[p1] + gebv2[p2]),
                     v_t1 = ctc$wf_var_t1, v_t2 = ctc$wf_var_t2,
                     stringsAsFactors = FALSE, row.names = NULL)
trait_specs <- data.frame(trait = c("t1", "t2"), mean_col = c("m_t1", "m_t2"),
                          var_col = c("v_t1", "v_t2"), stringsAsFactors = FALSE)
tau_lower <- c(stats::median(scores$m_t1), -Inf)     # t1 above median
tau_upper <- c(Inf, stats::median(scores$m_t2))      # t2 below median

ex <- ng_add_p_superior_progeny_multitrait(scores, trait_specs, tau_lower, tau_upper,
        k_progeny = 50L, cross_trait_cov = ctc, out_col = "p_exact")
px <- ng_add_p_superior_progeny_multitrait(scores, trait_specs, tau_lower, tau_upper,
        k_progeny = 50L, G_hat = diag(2), out_col = "p_proxy")   # zero cross-trait correlation

stopifnot(all(is.finite(ex$p_exact)), all(ex$p_exact >= 0 & ex$p_exact <= 1))
stopifnot(all(is.finite(px$p_proxy)), all(px$p_proxy >= 0 & px$p_proxy <= 1))
stopifnot(isTRUE(attr(ex, "p_superior_progeny_mt")$exact_cross_trait_cov))
# The exact within-family covariance (real cross-trait correlation) must move the probability vs
# the zero-correlation proxy.
stopifnot(max(abs(ex$p_exact - px$p_proxy)) > 1e-4)
# Row order independence: shuffle scores; matching is by parent pair, so results must be unchanged.
perm <- sample(nrow(scores))
ex2 <- ng_add_p_superior_progeny_multitrait(scores[perm, ], trait_specs, tau_lower, tau_upper,
         k_progeny = 50L, cross_trait_cov = ctc, out_col = "p_exact")
stopifnot(isTRUE(all.equal(ex$p_exact[perm], ex2$p_exact, tolerance = 1e-10)))

cat("multitrait_exact_cov_integration.R: PASS\n")
