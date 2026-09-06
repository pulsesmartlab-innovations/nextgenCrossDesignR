ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# I2 — ng_posterior_multitrait_cross_predict()'s threshold path (tau_lower_vec/tau_upper_vec,
# do_threshold = TRUE) calls ng_p_superior_progeny_multitrait() once per (posterior draw, cross
# pair), which calls mvtnorm::pmvnorm(). At dimension >= 3 traits pmvnorm() switches from a
# deterministic closed form (<= 2 traits) to the randomised GenzBretz lattice rule, which both
# DRAWS FROM and ADVANCES the ambient .Random.seed. Left unguarded (the pre-fix state of
# R/32_posterior_multitrait.R), a >= 3-trait posterior threshold call would perturb whatever the
# caller's ambient RNG stream happens to be -- this test fails RED against that unfixed code and
# passes GREEN once the call is scoped under its own fixed internal seed (R/32's
# NG_POSTERIOR_MT_THRESHOLD_PMVNORM_SEED via ng_with_rng_seed(), R/00_utils.R). Mirrors
# tests/check_reference_invariant.R's I1 case for the sibling check-reference bug fixed in
# 6734c5f, but tests the RNG-isolation property directly rather than via a downstream optimizer
# decision (this function has no optimizer step of its own to observe).

make_case <- function(n_traits) {
  set.seed(909)
  n <- 16L; m <- 40L
  ids <- paste0("L", seq_len(n))
  markers <- paste0("M", seq_len(m))
  geno <- matrix(2L * rbinom(n * m, 1L, 0.4), n, m, dimnames = list(ids, markers))
  mk <- data.frame(marker = markers, chr = rep(1:2, length.out = m),
                   pos_cm = rep(seq(0, 60, length.out = m / 2L), 2L),
                   stringsAsFactors = FALSE)
  trait_names <- paste0("t", seq_len(n_traits))
  Y <- matrix(NA_real_, n, n_traits, dimnames = list(ids, trait_names))
  for (j in seq_len(n_traits)) {
    Y[, j] <- as.numeric(geno %*% rnorm(m, sd = 0.05) + rnorm(n, sd = 0.4))
  }
  traits <- ng_multitrait_spec(trait = trait_names,
                               direction = rep("maximize", n_traits),
                               economic_weight = rep(1, n_traits))
  pair_ids <- ids[seq_len(8L)]
  pairs <- ng_make_pairs(pair_ids)
  G <- diag(n_traits); dimnames(G) <- list(trait_names, trait_names)
  list(geno = geno[pair_ids, , drop = FALSE], Y = Y[pair_ids, , drop = FALSE],
       traits = traits, marker_map = mk, ids = pair_ids, pairs = pairs, G = G,
       n_traits = n_traits)
}

run_threshold <- function(case, n_draws = 12L) {
  ng_posterior_multitrait_cross_predict(
    geno = case$geno, Y = case$Y, traits = case$traits, marker_map = case$marker_map,
    ids = case$ids, pairs = case$pairs,
    n_draws = n_draws, posterior_method = "closed_form",
    genetic_covariance_method = "beta_posterior",
    genetic_covariance = case$G,
    index_method = "economic_index", value_mode = "mean",
    use_cpp = FALSE, seed = 7L,
    tau_lower_vec = rep(0, case$n_traits), tau_upper_vec = rep(Inf, case$n_traits),
    threshold_k_progeny = 20L
  )
}

# ---- Test 1: at 2 traits (deterministic pmvnorm), the call must not perturb the ambient stream --
case2 <- make_case(2L)
set.seed(12345L)
ref_draw_2 <- { set.seed(12345L); runif(1) }
set.seed(12345L)
invisible(run_threshold(case2))
after_draw_2 <- runif(1)
stopifnot(identical(ref_draw_2, after_draw_2))
cat("Test 1 (2 traits, no RNG leak -- sanity control): OK\n")

# ---- Test 2 (RED against unfixed code): at >= 3 traits, the ambient stream must ALSO be
# untouched. This is exactly the scenario pmvnorm's randomised lattice rule can perturb, and the
# one the 2-trait case above cannot detect.
case3 <- make_case(3L)
set.seed(12345L)
ref_draw_3 <- { set.seed(12345L); runif(1) }
set.seed(12345L)
invisible(run_threshold(case3))
after_draw_3 <- runif(1)
if (!identical(ref_draw_3, after_draw_3)) {
  stop(sprintf(
    "RNG LEAK: ng_posterior_multitrait_cross_predict()'s threshold path perturbed the ambient RNG stream at n_traits=3 (ref=%.15f, after=%.15f). Every pmvnorm() call inside must be scoped under a fixed seed (see R/32_posterior_multitrait.R, mirroring the I1 fix in R/33_threshold_probability_multitrait.R).",
    ref_draw_3, after_draw_3))
}
cat("Test 2 (3 traits, no RNG leak): OK\n")

case4 <- make_case(4L)
set.seed(999L)
ref_draw_4 <- { set.seed(999L); runif(1) }
set.seed(999L)
invisible(run_threshold(case4))
after_draw_4 <- runif(1)
stopifnot(identical(ref_draw_4, after_draw_4))
cat("Test 3 (4 traits, no RNG leak): OK\n")

# ---- Test 4: the function's OWN result must also be reproducible run to run at >= 3 traits,
# independent of the ambient RNG state at call time (not just non-perturbing of it).
set.seed(111L)
out_a <- run_threshold(case3)
set.seed(222L)  # a DIFFERENT ambient seed
out_b <- run_threshold(case3)
stopifnot(isTRUE(all.equal(
  out_a$p_superior_progeny_mt_post_mean, out_b$p_superior_progeny_mt_post_mean,
  tolerance = 0
)))
cat("Test 4 (3-trait result independent of ambient RNG state): OK\n")

# ---- Test 5 (hazard guard): the fix must not collapse the posterior sample to a constant. A
# naive wrap that resets the SAME seed before every single draw/pair call (rather than once
# around the whole per-draw loop) is the wrong scope to reach for; verify draws still differ from
# each other by requiring non-trivial CI width across the posterior draws, the same style of
# guard tests/threshold_probability_multitrait.R's Test 9 already uses for the 2-trait case.
out_c <- run_threshold(case3, n_draws = 25L)
ci_widths <- out_c$p_superior_progeny_mt_post_upper - out_c$p_superior_progeny_mt_post_lower
stopifnot(all(is.finite(ci_widths)))
if (max(ci_widths, na.rm = TRUE) < 1e-3) {
  stop("posterior CI widths are all near zero at n_traits=3 -- the RNG scoping fix may have collapsed the posterior sample to a constant instead of preserving draw-to-draw variation")
}
cat(sprintf("Test 5 (posterior draws remain distinct, max CI width=%.4f): OK\n", max(ci_widths, na.rm = TRUE)))

cat("posterior_multitrait_rng_scoping.R: 5/5 checks passed\n")
