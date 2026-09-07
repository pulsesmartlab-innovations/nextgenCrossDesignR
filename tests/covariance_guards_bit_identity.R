helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 0.29.0 -- A VALID P/G PAIR MUST BE NUMERICALLY UNTOUCHED BY THE NEW GUARDS.
#
# The four covariance guards added in 0.29.0 refuse invalid input; they must not
# perturb a single digit of a run whose matrices satisfy all of them. The
# constants below were produced by a PRISTINE 0.28.0 tree (git archive HEAD of
# commit a42dda6, unpacked to /tmp and loaded with the same helper), as IEEE-754
# hex literals so they round-trip exactly.
#
# If a later change makes this file fail, the guards have stopped being free.

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }

# Comparing against constants frozen on ONE machine, from any machine.
#
# These references were captured on macOS. R's linear algebra runs through the
# platform's BLAS/LAPACK -- Accelerate there, OpenBLAS on the Linux CI runner --
# and the two are NOT bit-reproducible with each other: eigen() and solve()
# legitimately differ in the last bits. identical() against a frozen constant
# therefore asserts CROSS-PLATFORM bit reproducibility, which LAPACK cannot
# provide and which this file never meant to claim. (It failed exactly that way
# on the Linux harness while passing on macOS.)
#
# What this file does claim is that the 0.29.0 guards did not alter results for
# a valid P/G pair. That is a before-vs-after claim, so compare numerically at a
# tolerance far tighter than any real algorithmic change could hide beneath, and
# print the observed difference so genuine drift is visible in the log even on a
# pass. Non-numeric results (the crossing plan) stay on identical(): character
# comparison is exact on every platform.
ng_same_num <- function(actual, ref, label, tol = 1e-12) {
  a <- as.numeric(actual); r <- as.numeric(ref)
  if (length(a) != length(r))
    stop(label, ": length ", length(a), " vs reference ", length(r))
  rel <- max(abs(a - r) / pmax(abs(r), 1))
  cat(sprintf("  %-30s max rel diff = %.3e\n", label, rel))
  if (!isTRUE(rel <= tol))
    stop(label, ": max relative difference ", format(rel, digits = 6),
         " exceeds tolerance ", tol)
  invisible(TRUE)
}

# ---- reference values from the pristine 0.28.0 tree ------------------------
# IEEE-754 hex float literals: exact, not rounded decimal. as.numeric() parses them.
EI_SCORE <- as.numeric(c("-0x1.ff5a9e9544232p-3", "-0x1.a49e4a7861d9ap-1", "-0x1.7424152793fcep-7", "0x1.9e9685de5bd16p-2", "0x1.aa75e416cb098p-3", "0x1.4e727075486a6p-1", "0x1.0d9ccaa3f0ea6p-2", "-0x1.821f57d79850ep-2", "-0x1.2a5ce70cedc28p-1", "-0x1.d9df7861eb329p-7", "-0x1.0cbe86b946b91p+0", "0x1.8b5a1a81ebf9p-2", "0x1.0fd08671ed78ep-1", "-0x1.a12794f830604p-3", "0x1.fd5991ab8a9c7p-4", "0x1.92815ea83b178p-4", "-0x1.0ed615d1b767ap-3", "0x1.cafe2a53cde7cp-2", "-0x1.868505f4c9d5bp-1", "-0x1.de0e8e3884e54p-2", "0x1.938f3015c441p-2", "0x1.02229d90bd62fp+0", "-0x1.273e17cb7a1p-6", "0x1.625ff76559d4fp-2", "0x1.a4f1bc905ad94p-4", "0x1.ab002d5d95db7p-3", "-0x1.02fe81d0938a6p-1", "-0x1.8191bb6d6b46fp-1", "-0x1.9a73fb50f873p-1", "0x1.1621bd8eac0fcp-2", "0x1.30e265b6fa4f5p+0", "0x1.cc692c08741adp-2", "-0x1.85171a2208d03p-1", "-0x1.aaa9803728f9dp-4", "0x1.96910b5bb324bp-3", "-0x1.3eff8a9bee855p-2", "-0x1.c2510189b9b29p-8", "-0x1.54c1a60d5402p-5", "-0x1.913594a61ea93p-1", "0x1.72b2b50f011dp-1", "0x1.da3d06bf85b56p-1", "-0x1.670b787174504p-1", "-0x1.1926d99762bdbp-2", "0x1.0a22c2238b27p-3", "-0x1.f340b8517f3f8p-4", "0x1.c0976ae60ca9ep-2", "-0x1.f7fc27ecd4718p-3", "-0x1.653f2f505bc83p-2", "0x1.afed06eda516fp-1", "0x1.9dfb74799bc07p-3", "0x1.c18c5b341bc0dp-3", "-0x1.b29723a202957p-8", "0x1.a914eae9ba3a3p-3", "-0x1.286ddfd45adc8p-3", "-0x1.5842c43c9ca41p-2", "0x1.3b4ddd596f1f1p-1", "0x1.d85d0a93fc48dp-1", "0x1.fe19c90420c3cp-2", "0x1.5a979c29ebcc1p-5", "-0x1.12a37da7638adp-1"))

DG_SCORE <- as.numeric(c("-0x1.393fff1a26d61p-4", "-0x1.724105224fd89p-1", "0x1.1c22e349ffdfep-4", "0x1.77c7b896ffc8ap-4", "-0x1.cd88523c33d53p-6", "0x1.f4a4d35e15277p-1", "0x1.77087211d42e8p-6", "0x1.6fe78470e6f35p-4", "-0x1.aec29e7a04f1ap-5", "0x1.358122ae3b7dfp-7", "-0x1.e8b68fe2a0125p-2", "0x1.25d8a6dc0dbc7p-1", "0x1.2a69b2c72fc43p-1", "0x1.b2f1edb5326dfp-2", "-0x1.5fcee1a065593p-3", "0x1.37be961dbd5afp-1", "0x1.22a02b9c43c37p-2", "0x1.5adca8f6a58d2p-1", "-0x1.e1b879082550ap-2", "-0x1.150474c02fd56p+0", "0x1.8d1f1ba0aaddcp-2", "0x1.03c1fcc800aacp+0", "-0x1.6dc5ac419fc67p-2", "0x1.cc351d160d75ep-2", "0x1.673f5cc2b035p-4", "0x1.3d2d942e9a4fp-1", "-0x1.ffad5459c04fdp-2", "-0x1.12a9b789ac0bcp+0", "-0x1.1b8beb84b8c6ap+0", "0x1.5927dae144f9fp-2", "0x1.bb57df38fbf57p-1", "0x1.1dcf853b4222fp-1", "-0x1.bae0c9ad30632p-2", "0x1.816a17bc09fc3p-3", "0x1.23ea3bf71cd72p-2", "0x1.80a5dd62187a7p-4", "-0x1.66b1dcd106e24p-3", "-0x1.565ce62bdda86p-2", "-0x1.d2ab2a0a7545ep-1", "0x1.acca43526c328p-3", "0x1.93c3b343c0d3cp-1", "-0x1.69b7e54c74eadp-1", "0x1.9abe562396acap-2", "0x1.1953ab8ff1b3fp-3", "-0x1.0cf5275283cffp-5", "0x1.91ead9d7ee278p-1", "-0x1.a3aa4bdf79dd8p-2", "-0x1.5dce5dde59bcdp-1", "0x1.44a5a47e0b015p-1", "0x1.0fa6b238dc287p-1", "0x1.15cab59ad558p-3", "-0x1.8627ff414ec85p-2", "-0x1.04377ace43f2p-3", "-0x1.02726de3c5254p-4", "-0x1.57d38bc06c0b7p-3", "0x1.951cc47bb778bp-2", "0x1.2f171b52af48ap-1", "0x1.1a178c6c07292p-3", "-0x1.053ba79ef43efp-2", "-0x1.de745ed30df56p-2"))

EI_COEF <- as.numeric(c("0x1.570196611ee4bp-2", "0x1.93fe55bda1f83p-3", "0x1.deff3ec0101f3p-2"))
DG_COEF <- as.numeric(c("0x1.0a66ff39ef4c7p-1", "0x1.f1433869ee6cep-3", "0x1.e520caae54613p-3"))
EI_SD <- as.numeric("0x1.bd53f314a3deep-1")
DG_SD <- as.numeric("0x1.0a4864c42da4dp+0")
PLAN <- c("P2|P14", "P2|P5", "P3|P7", "P4|P7")

# ---- identical fixture to the one the 0.28.0 tree was run on ---------------
set.seed(20290101)
n <- 60L
ids <- paste0("P", seq_len(20L))
pairs <- t(utils::combn(ids, 2L))[seq_len(n), , drop = FALSE]
yield   <- stats::rnorm(n, 10, 2.0)
protein <- stats::rnorm(n, 13, 1.5) - 0.3 * yield
lodging <- stats::rnorm(n, 20, 4.0) + 0.4 * yield
scores <- data.frame(parent1 = pairs[, 1L], parent2 = pairs[, 2L],
                     yield = yield, protein = protein, lodging = lodging,
                     stringsAsFactors = FALSE)
traits <- ng_multitrait_spec(
  trait = c("yield", "protein", "lodging"),
  direction = c("maximize", "maximize", "minimize"),
  economic_weight = c(2, 1, 1), desired_change = c(4, 2, 3))
tn <- traits$trait
G <- matrix(c( 4.0, 1.2, -1.5,
               1.2, 2.0, -0.8,
              -1.5, -0.8, 9.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))
P <- matrix(c(10.0, 2.0, -2.5,
               2.0, 5.0, -1.2,
              -2.5, -1.2, 20.0), 3L, 3L, byrow = TRUE, dimnames = list(tn, tn))

# Every guard is satisfied, and visibly so.
stopifnot(min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) > 0)
stopifnot(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values) > 0)
stopifnot(min(eigen(P - G, symmetric = TRUE, only.values = TRUE)$values) > 0)
stopifnot(max(abs(P - t(P))) == 0, max(abs(G - t(G))) == 0)

ei <- ng_add_multitrait_score(scores, traits, method = "economic_index",
                              phenotypic_covariance = P, genetic_covariance = G)
dg <- ng_add_multitrait_score(scores, traits, method = "desired_gain",
                              phenotypic_covariance = P, genetic_covariance = G)
m_ei <- attr(ei, "multi_trait")
m_dg <- attr(dg, "multi_trait")

ng_same_num(ei$multi_trait_score, EI_SCORE, "EI_SCORE")
ng_same_num(dg$multi_trait_score, DG_SCORE, "DG_SCORE")
ok("economic_index and desired_gain scores are unchanged from 0.28.0")

ng_same_num(m_ei$economic_index_coefficients, EI_COEF, "EI_COEF")
ng_same_num(m_dg$desired_gain_coefficients, DG_COEF, "DG_COEF")
ng_same_num(m_ei$economic_index_index_sd, EI_SD, "EI_SD")
ng_same_num(m_dg$desired_gain_index_sd, DG_SD, "DG_SD")
ok("index coefficients and index SDs are unchanged from 0.28.0")

plan_scores <- scores
plan_scores$pair_kinship <- 0
K <- diag(length(ids)); dimnames(K) <- list(ids, ids)
plan <- ng_optimize_multitrait_mating_plan(
  scores = plan_scores, traits = traits, n_crosses = 4L, parent_kinship = K,
  multitrait_method = "economic_index", optimizer_method = "greedy_local",
  max_crosses_per_parent = 2L, phenotypic_covariance = P, genetic_covariance = G)
stopifnot(identical(paste(plan$parent1, plan$parent2, sep = "|"), PLAN))
ok("the selected crossing plan is the same set of crosses, in the same order")

cat(sprintf("covariance_guards_bit_identity.R: PASS (%d checks)\n", checks))
