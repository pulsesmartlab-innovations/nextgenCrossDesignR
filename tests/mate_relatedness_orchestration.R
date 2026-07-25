## Preventive orchestration of the per-cross relatedness penalty.
## The unified `mate_relatedness` control drives EXACTLY ONE of the two collinear lambdas,
## so stacking is impossible by construction; setting BOTH raw lambdas (or the unified
## control AND a raw lambda) is a HARD error.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = function(e) e)
  stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), ignore.case = TRUE))
}

set.seed(1)
ids <- sprintf("P%02d", 1:8)
pairs <- ng_make_pairs(ids, include_self = FALSE)
scores <- data.frame(parent1 = pairs$parent1, parent2 = pairs$parent2,
                     gain = rnorm(nrow(pairs), 10, 2),
                     pair_kinship = pmax(0, rnorm(nrow(pairs), 0.1, 0.1)),
                     stringsAsFactors = FALSE)
scores$expected_progeny_inbreeding <- pmax(0, scores$pair_kinship / 2)
K <- diag(length(ids)); dimnames(K) <- list(ids, ids)
plan <- function(...) ng_optimize_mating_plan(scores, n_crosses = 5, gain_col = "gain",
                                              parent_kinship = K, method = "greedy_local", ...)

## 1. unified 'avoid_inbreeding' == raw lambda_progeny_inbreeding (same single lambda)
a <- plan(mate_relatedness = "avoid_inbreeding", mate_relatedness_weight = 0.5)
b <- plan(lambda_progeny_inbreeding = 0.5)
stopifnot(identical(paste(a$parent1, a$parent2), paste(b$parent1, b$parent2)))

## 2. unified 'favor_complementarity' == raw lambda_mating
c1 <- plan(mate_relatedness = "favor_complementarity", mate_relatedness_weight = 0.5)
d1 <- plan(lambda_mating = 0.5)
stopifnot(identical(paste(c1$parent1, c1$parent2), paste(d1$parent1, d1$parent2)))

## 3. off (default) with a single raw lambda still works (advanced escape hatch)
stopifnot(nrow(plan(lambda_mating = 0.02)) == 5)
stopifnot(nrow(plan()) == 5)   # off + no raw lambdas: no relatedness penalty

## 4. HARD error: both raw lambdas > 0 (the collinear stacking that used to only warn)
expect_error(plan(lambda_mating = 0.02, lambda_progeny_inbreeding = 0.1),
             "both penalize parent-pair relatedness")

## 5. HARD error: unified control AND a raw lambda at once
expect_error(plan(mate_relatedness = "avoid_inbreeding", mate_relatedness_weight = 0.1,
                  lambda_mating = 0.02),
             "EITHER mate_relatedness OR the raw")

cat("mate_relatedness_orchestration.R: PASS\n")
