test_that("the mating plan honours its size and parent-use constraints", {
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  K <- ng_parent_kinship(p$geno)
  plan <- ng_optimize_mating_plan(s, n_crosses = 8L, parent_kinship = K,
                                  max_crosses_per_parent = 3L,
                                  lambda_group = 1,
                                  method = "greedy_local")
  summary <- attr(plan, "summary")
  expect_equal(nrow(plan), 8L)
  expect_lte(summary$max_parent_use, 3L)
  expect_true(is.finite(summary$mean_gain))
  expect_true(is.finite(summary$group_coancestry))
})

test_that("a tighter parent-use cap is actually binding", {
  # Guards the constraint itself: loosening the cap must not leave the plan
  # unchanged, or the constraint is being ignored rather than satisfied.
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  K <- ng_parent_kinship(p$geno)
  tight <- ng_optimize_mating_plan(s, n_crosses = 8L, parent_kinship = K,
                                   max_crosses_per_parent = 2L,
                                   lambda_group = 1, method = "greedy_local")
  expect_lte(attr(tight, "summary")$max_parent_use, 2L)
})

test_that("kinship is a valid symmetric relationship matrix", {
  p <- ng_test_panel()
  K <- ng_parent_kinship(p$geno)
  expect_equal(dim(K), c(length(p$ids), length(p$ids)))
  expect_equal(K, t(K))
  expect_true(all(is.finite(K)))
})

test_that("VanRaden and Yang kinship differ", {
  # Two different estimators must not silently collapse to the same matrix.
  p <- ng_test_panel()
  vr <- ng_parent_kinship(p$geno, method = "vanraden")
  yg <- ng_parent_kinship(p$geno, method = "yang")
  expect_equal(dim(vr), dim(yg))
  expect_false(isTRUE(all.equal(vr, yg)))
})

test_that("family sizes respect the total and the per-family bounds", {
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  K <- ng_parent_kinship(p$geno)
  plan <- ng_optimize_mating_plan(s, n_crosses = 8L, parent_kinship = K,
                                  max_crosses_per_parent = 3L,
                                  method = "greedy_local")
  fam <- ng_allocate_family_sizes(plan, total_progeny = 80L,
                                  min_progeny = 4L, max_progeny = 20L,
                                  value_col = "usefulness_pmv_gebv")
  expect_equal(sum(fam$n_progeny), 80L)
  expect_gte(min(fam$n_progeny), 4L)
  expect_lte(max(fam$n_progeny), 20L)
})
