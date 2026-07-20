test_that("every parent pair is scored exactly once", {
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  expect_equal(nrow(s), choose(length(p$ids), 2L))
  expect_false(any(duplicated(
    paste(pmin(s$parent1, s$parent2), pmax(s$parent1, s$parent2))
  )))
  expect_true(all(s$parent1 != s$parent2))
})

test_that("cross scores are finite and variances are non-negative", {
  s <- ng_test_scores()
  expect_true(all(is.finite(s$cross_mean)))
  expect_true(all(is.finite(s$usefulness_pmv_gebv)))
  expect_true(all(s$parent_distance >= -1e-8))
  expect_true(all(s$vpm >= -1e-8))
  expect_true(all(s$pmv >= -1e-8))
})

test_that("PMV is never below VPM", {
  # PMV adds Var(beta) to the diagonal of the quadratic form, so it can only
  # widen the predicted family variance relative to the VPM point estimate.
  s <- ng_test_scores()
  expect_true(all(s$pmv + 1e-12 >= s$vpm))
})

test_that("recombination-aware variance is not the linkage-free variance", {
  # parent_distance ignores linkage; vpm applies the Haldane recursion.
  # If these ever coincide the recombination model has stopped doing anything.
  s <- ng_test_scores()
  expect_false(isTRUE(all.equal(s$parent_distance, s$vpm)))
})

test_that("the C++ kernel and the pure-R path agree", {
  skip_if_not(nzchar(Sys.getenv("NG_TEST_CPP", unset = "")),
              "Set NG_TEST_CPP=1 to compile and cross-check the C++ kernel.")
  p <- ng_test_panel()
  fit <- ng_fit_ridge_effects(p$geno, p$y, ids = p$ids, kfold = 3L, seed = 1L)
  args <- list(p$geno, fit, marker_map = p$marker_map, ids = p$ids,
               adjusted_pheno = p$y, selection_prop = 0.1)
  r_scores <- do.call(ng_score_crosses, c(args, list(use_cpp = FALSE)))
  cpp_scores <- do.call(ng_score_crosses, c(args, list(use_cpp = TRUE)))
  expect_equal(r_scores$vpm, cpp_scores$vpm, tolerance = 1e-8)
  expect_equal(r_scores$pmv, cpp_scores$pmv, tolerance = 1e-8)
})
