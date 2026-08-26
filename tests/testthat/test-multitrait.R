ng_test_qg_covariances <- function(traits) {
  traits <- as.character(traits)
  p <- length(traits)
  P <- diag(1, p)
  G <- diag(0.6, p)
  dimnames(P) <- dimnames(G) <- list(traits, traits)
  list(P = P, G = G)
}

test_that("a multi-trait spec round-trips its traits and directions", {
  spec <- ng_multitrait_spec(
    trait = c("yield", "disease"),
    direction = c("maximize", "minimize"),
    economic_weight = c(2, 1)
  )
  expect_equal(spec$trait, c("yield", "disease"))
  expect_equal(spec$direction, c("maximize", "minimize"))
})

test_that("the multi-trait index scores every cross", {
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  s$yield <- s$cross_mean
  s$disease <- -s$cross_mean + rnorm(nrow(s), sd = 0.1)
  spec <- ng_multitrait_spec(
    trait = c("yield", "disease"),
    direction = c("maximize", "minimize"),
    economic_weight = c(2, 1)
  )
  qg <- ng_test_qg_covariances(spec$trait)
  scored <- ng_add_multitrait_score(
    s, spec, method = "economic_index",
    phenotypic_covariance = qg$P,
    genetic_covariance = qg$G
  )
  expect_equal(nrow(scored), nrow(s))
  expect_true(all(is.finite(scored$multi_trait_score)))
})

test_that("direction is honoured: a minimize trait is penalised, not rewarded", {
  # The clearest way this can silently break is a sign error, which would make
  # the index reward the trait the breeder asked to reduce.
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  s$only_trait <- s$cross_mean
  qg <- ng_test_qg_covariances("only_trait")
  up <- ng_add_multitrait_score(
    s, ng_multitrait_spec(trait = "only_trait", direction = "maximize",
                          economic_weight = 1),
    method = "economic_index",
    phenotypic_covariance = qg$P,
    genetic_covariance = qg$G)
  down <- ng_add_multitrait_score(
    s, ng_multitrait_spec(trait = "only_trait", direction = "minimize",
                          economic_weight = 1),
    method = "economic_index",
    phenotypic_covariance = qg$P,
    genetic_covariance = qg$G)
  # Flipping the direction must flip how the score ranks the same crosses.
  expect_lt(cor(up$multi_trait_score, down$multi_trait_score), 0)
})

test_that("the economic index reports which covariance it solved with", {
  p <- ng_test_panel()
  s <- ng_test_scores(p)
  s$yield <- s$cross_mean
  s$disease <- -s$cross_mean + rnorm(nrow(s), sd = 0.1)
  spec <- ng_multitrait_spec(trait = c("yield", "disease"),
                             direction = c("maximize", "minimize"),
                             economic_weight = c(2, 1))
  qg <- ng_test_qg_covariances(spec$trait)
  scored <- ng_add_multitrait_score(
    s, spec, method = "economic_index",
    phenotypic_covariance = qg$P,
    genetic_covariance = qg$G
  )
  meta <- attr(scored, "multi_trait")
  expect_true(nzchar(meta$economic_index_cov_source))
  expect_true(nzchar(meta$economic_index_cov_solve_form))
})
