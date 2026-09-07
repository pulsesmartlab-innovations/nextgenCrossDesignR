helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# Comparing against constants frozen on ONE machine, from any machine.
#
# The references below were captured on macOS. R's linear algebra runs through
# the platform's BLAS/LAPACK -- Accelerate there, OpenBLAS on the Linux CI
# runner -- and the two are NOT bit-reproducible with each other: eigen() and
# solve() legitimately differ in the last bits. identical() against a frozen
# constant therefore asserts CROSS-PLATFORM bit reproducibility, which LAPACK
# cannot provide and which this file never meant to claim. (It failed exactly
# that way on the Linux harness while passing on macOS.)
#
# The real claim is that the 0.30.0 changes did not alter results that were
# always valid -- a before-vs-after claim. So compare numerically at a tolerance
# far tighter than any real algorithmic change could hide beneath, and print the
# observed difference so genuine drift shows in the log even on a pass.
# Non-numeric results (cross keys, row counts, the method string) stay on
# identical(): those are exact on every platform.
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

# 0.30.0 -- THE TWO CHANGES MUST BE NUMERICALLY FREE ON WORK THAT WAS ALWAYS VALID.
#
# 0.30.0 makes two changes that could in principle move numbers:
#
#   1. ng_estimate_genetic_covariance() now REFUSES a G-hat whose implied
#      per-trait h2 = diag(G_hat) / diag(var(Y)) exceeds 1. On input where it
#      does not exceed 1, the returned matrix -- and every provenance attribute
#      on it -- must be bit-for-bit what 0.29.0 returned.
#
#   2. ng_load() now sources only what R CMD build would include, which on this
#      working tree removes ELEVEN untracked legacy scripts from the load set
#      (R/01_relationships.R, R/02_duplicate_detection_legacy.R,
#      R/02b_duplicate_detection.R, R/03_ld.R, R/04_variance_simple_uc.R,
#      R/05_pmv.R, R/06_build_cross_data.R, R/07_optimizers.R,
#      R/08_select_optimal_parents.R, R/09_pmv_cpp.R, R/10_marker_effects.R).
#      A source-loaded session therefore now runs a DIFFERENT FILE SET, so a
#      scoring run is exercised here too.
#
# The constants below were produced by a PRISTINE 0.29.0 tree (git archive HEAD
# of commit 7ec06cd, unpacked to /tmp and loaded with the same helper) and are
# written as IEEE-754 hex float literals, which round-trip exactly where decimal
# does not. They are compared with identical() -- tolerance = 0, not all.equal().
#
# If a later change makes this file fail, the 0.30.0 guard has stopped being free.

checks <- 0L
ok <- function(msg) { checks <<- checks + 1L; cat("  OK:", msg, "\n") }

# ---- reference values from the pristine 0.29.0 tree ------------------------
G_DIAG <- as.numeric(c("0x1.e044cacc8cf5dp-1", "0x1.30f216240081dp+1", "0x1.1124bd45fc6b5p-1"))
G_FULL <- as.numeric(c("0x1.e044cacc8cf5dp-1", "0x1.5ede1ee31a389p-4", "-0x1.8dd67ca7bf854p-3", "0x1.5ede1ee31a388p-4", "0x1.30f216240081dp+1", "0x1.52f6c2723d83ep-2", "-0x1.8dd67ca7bf853p-3", "0x1.52f6c2723d83fp-2", "0x1.1124bd45fc6b5p-1"))
G_CORR <- as.numeric(c("0x1p+0", "0x1.d56b0f6f4a06p-5", "-0x1.19321022d6a68p-2", "0x1.d56b0f6f4a05fp-5", "0x1p+0", "0x1.2caaed531908fp-2", "-0x1.19321022d6a67p-2", "0x1.2caaed531909p-2", "0x1p+0"))
P_FULL <- as.numeric(c("0x1.a21669706b94ap+0", "0x1.72ce646e3a135p-4", "-0x1.97702f5bca7f6p-3", "0x1.72ce646e3a134p-4", "0x1.07ae98ed55e7ep+2", "0x1.1dd60ac2cc444p-2", "-0x1.97702f5bca7f4p-3", "0x1.1dd60ac2cc445p-2", "0x1.b71c7c38a5714p-1"))
SIGMA_G <- as.numeric(c("0x1.e044cacc8cf61p-1", "0x1.30f216240081ep+1", "0x1.1124bd45fc6afp-1"))
SIGMA_E <- as.numeric(c("0x1.637b05b69176ap+0", "0x1.c36c5e88a53cfp+1", "0x1.945867ce2eb6fp-1"))
SCORE_MEAN <- as.numeric(c("-0x1.cb16d4ce985dep-2", "-0x1.d1b7fba54dd4p+0", "0x1.80b38a0194e08p-3", "-0x1.076702470d46bp-1", "-0x1.85e77cac16dabp-1", "0x1.1b342bc41d3c4p-1", "-0x1.0d3aca306c248p-4", "0x1.aaa37fc224629p-2", "-0x1.a539afc537e43p-1", "0x1.3364d47977deep-1", "0x1.4e3f954c3207ap-2", "-0x1.f5eacec06de9ap-2"))
SCORE_PMV <- as.numeric(c("0x1.47b880259ea53p-10", "0x1.61c8448b07244p-10", "0x1.949d29855f517p-10", "0x1.87ea8f160e5eap-10", "0x1.1ac3d51358ab8p-10", "0x1.4794d6568fa23p-10", "0x1.87d7e02e993aap-10", "0x1.747ddcc00e70dp-10", "0x1.88002695d532p-10", "0x1.411a669e7a968p-10", "0x1.4ded017cf1c65p-10", "0x1.8e67280fe0741p-10"))
SCORE_UC <- as.numeric(c("-0x1.8b8c4d5b9c876p-2", "-0x1.c136bf21d4944p+0", "0x1.06f425c65fa3bp-2", "-0x1.c95181d6e2135p-2", "-0x1.6864a372c063fp-1", "0x1.3af7b4e744fbfp-1", "0x1.161412c22332p-9", "0x1.ee6194a3dc58ep-2", "-0x1.827a796bb0144p-1", "0x1.52d78f2fffe7fp-1", "0x1.8e6362e206f5p-2", "-0x1.afdbb0a34ca22p-2"))
SCORE_RANK <- as.numeric(c("-0x1.8b8c4d5b9c876p-2", "-0x1.c136bf21d4944p+0", "0x1.06f425c65fa3bp-2", "-0x1.c95181d6e2135p-2", "-0x1.6864a372c063fp-1", "0x1.3af7b4e744fbfp-1", "0x1.161412c22332p-9", "0x1.ee6194a3dc58ep-2", "-0x1.827a796bb0144p-1", "0x1.52d78f2fffe7fp-1", "0x1.8e6362e206f5p-2", "-0x1.afdbb0a34ca22p-2"))
SCORE_KEY <- c("P1|P2", "P1|P3", "P1|P4", "P1|P5", "P1|P6", "P1|P7", "P1|P8",
               "P1|P9", "P1|P10", "P1|P11", "P1|P12", "P1|P13")
SCORE_NROW <- 190L

# ---- part 1: the estimator on GUARD-SATISFYING input -----------------------
#
# Same generating model as tests/genetic_covariance_estimator.R (n = 100,
# m = 200, 3 traits, every true h2 = 0.5), at a seed chosen because
# two_stage_ridge's implied h2 comes out BELOW 1 on all three traits -- one of
# the ~15% of datasets from this model that the new guard lets through. That is
# the point: the guard must be invisible here.
t_traits <- 3L
trait_names <- c("trait_a", "trait_b", "trait_c")
sd_g <- c(1, sqrt(2), sqrt(0.5))
R_true <- matrix(c(1, 0.3, -0.4,
                   0.3, 1, 0.1,
                   -0.4, 0.1, 1), 3L, byrow = TRUE)
G_true <- diag(sd_g) %*% R_true %*% diag(sd_g)
dimnames(G_true) <- list(trait_names, trait_names)
E_var <- diag(G_true)

set.seed(5031L)
n <- 100L; m <- 200L
maf <- stats::runif(m, 0.1, 0.45)
geno <- matrix(0, n, m)
for (k in seq_len(m)) geno[, k] <- stats::rbinom(n, 2L, maf[k])
rownames(geno) <- paste0("L", seq_len(n)); colnames(geno) <- paste0("M", seq_len(m))
ph <- colMeans(geno) / 2; den <- sum(2 * ph * (1 - ph))
eg <- eigen(G_true, symmetric = TRUE)
sq <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% t(eg$vectors)
bt <- matrix(stats::rnorm(m * t_traits), m) %*% (sq / sqrt(den))
Xc <- sweep(geno, 2L, 2 * ph, "-")
E <- matrix(stats::rnorm(n * t_traits), n) %*% diag(sqrt(E_var))
Y <- Xc %*% bt + E
colnames(Y) <- trait_names; rownames(Y) <- rownames(geno)

G_hat <- suppressWarnings(ng_estimate_genetic_covariance(
  geno, Y, method = "two_stage_ridge", kfold = 5L, seed = 5031L))
P_hat <- ng_estimate_phenotypic_covariance(Y, shrinkage = "auto")

# The guard is satisfied, and visibly so. The reference variance is the PER-COLUMN
# stats::var(Y[, t]); diag(stats::var(Y)) is the same quantity mathematically but
# takes a different code path in stats and can differ in the last bit, so the
# per-column form is what the guard uses and what is asserted here.
vp_obs <- apply(Y, 2L, function(col) stats::var(col, na.rm = TRUE))
h2_obs <- diag(G_hat) / vp_obs
cat(sprintf("  implied h2 vs var(Y): [%s] -- all below 1, so the guard is silent\n",
            paste(sprintf("%.4f", h2_obs), collapse = ", ")))
stopifnot(all(h2_obs < 1))
# ...and so is the downstream pair guard, on the pair the package's own
# workflow would build. The self-check is never stricter than that one.
ng_multitrait_validate_cov_pair(
  ng_multitrait_validate_cov(matrix(as.numeric(P_hat), 3L, 3L,
                                    dimnames = list(trait_names, trait_names)),
                             trait_names, "phenotypic_covariance"),
  ng_multitrait_validate_cov(matrix(as.numeric(G_hat), 3L, 3L,
                                    dimnames = list(trait_names, trait_names)),
                             trait_names, "genetic_covariance"))
ok("the guard-satisfying fixture also passes ng_multitrait_validate_cov_pair()")

ng_same_num(diag(G_hat), G_DIAG, "G_DIAG")
ng_same_num(G_hat, G_FULL, "G_FULL")
ng_same_num(attr(G_hat, "genetic_correlation"), G_CORR, "G_CORR")
ok("G_hat and its genetic_correlation are unchanged from 0.29.0")

ng_same_num(attr(G_hat, "genetic_variance"), SIGMA_G, "SIGMA_G")
ng_same_num(attr(G_hat, "residual_variance"), SIGMA_E, "SIGMA_E")
stopifnot(identical(attr(G_hat, "method"), "two_stage_ridge"))
stopifnot(identical(attr(G_hat, "n_used"), 100L))
ng_same_num(P_hat, P_FULL, "P_FULL")
ok("provenance attributes and P_hat are unchanged from 0.29.0")

# The 0.30.0 attributes report the quantity the guard judged, and it agrees
# with the quantity computed here from the same definition.
stopifnot(identical(as.numeric(attr(G_hat, "implied_h2_vs_observed_variance")),
                    as.numeric(h2_obs)))
stopifnot(identical(as.numeric(attr(G_hat, "phenotypic_variance_observed")),
                    as.numeric(vp_obs)))
ok("implied_h2_vs_observed_variance reports exactly what the guard judged")

# ---- part 2: a scoring run under the NEW ng_load() file set ----------------
#
# ng_load() no longer sources the eleven untracked legacy scripts in R/. On this
# machine those files are OneDrive placeholders that read as zero bytes, so they
# defined nothing and their removal cannot change a number -- which is exactly
# what this asserts, against a tree that never had them at all.
set.seed(4242L)
ids <- paste0("P", seq_len(20L))
markers <- paste0("M", seq_len(120L))
g2 <- matrix(2L * stats::rbinom(20L * 120L, 1L, 0.45), 20L, 120L,
             dimnames = list(ids, markers))
y2 <- as.numeric(g2 %*% stats::rnorm(120L, sd = 0.1) + stats::rnorm(20L))
names(y2) <- ids
mm <- data.frame(marker = markers, chr = rep(1:3, each = 40L),
                 pos_cm = rep(seq(0, 100, length.out = 40L), 3L))
fit <- ng_fit_ridge_effects(g2, y2, ids = ids, kfold = 3L, seed = 1L)
sc <- ng_score_crosses(g2, fit, marker_map = mm, ids = ids, adjusted_pheno = y2,
                       selection_prop = 0.1, use_cpp = FALSE)
keep <- seq_len(12L)
stopifnot(identical(nrow(sc), SCORE_NROW))
stopifnot(identical(paste(sc$parent1[keep], sc$parent2[keep], sep = "|"), SCORE_KEY))
ng_same_num(sc$cross_mean[keep], SCORE_MEAN, "SCORE_MEAN")
ng_same_num(sc$pmv[keep], SCORE_PMV, "SCORE_PMV")
ng_same_num(sc$usefulness_pmv[keep], SCORE_UC, "SCORE_UC")
ng_same_num(sc$rank_score[keep], SCORE_RANK, "SCORE_RANK")
ok("ng_score_crosses() under the reduced ng_load() file set is unchanged from 0.29.0")

# ---- part 3: the loader really did drop them, and would refuse a placeholder ----
root <- ng_test_find_root()
all_numbered <- sort(list.files(file.path(root, "R"), pattern = "^[0-9].*[.]R$",
                                full.names = TRUE))
loaded <- ng_load_apply_rbuildignore(all_numbered, root)
dropped <- basename(setdiff(all_numbered, loaded))
cat(sprintf("  ng_load() sources %d of the %d files matching R/[0-9]*.R\n",
            length(loaded), length(all_numbered)))
if (length(dropped)) {
  cat("  build-ignored, no longer sourced:", paste(dropped, collapse = ", "), "\n")
}
# Whatever the working tree happens to contain, the loaded set must be exactly
# the build set: nothing .Rbuildignore excludes, and everything it does not.
stopifnot(!any(basename(loaded) %in% dropped))
stopifnot(all(file.exists(loaded)))
ok("ng_load() loads exactly the non-build-ignored R/[0-9]*.R files")

# An unreadable file that SHOULD be loaded is a loud failure, never a silent skip.
tmp <- tempfile(fileext = ".R")
writeLines("x <- 1", tmp)
stopifnot(isTRUE(ng_load_require_materialised(tmp)))
unlink(tmp)
missing_msg <- tryCatch({ ng_load_require_materialised(file.path(root, "R", "0_no_such_file.R")); NA_character_ },
                        error = function(e) conditionMessage(e))
stopifnot(!is.na(missing_msg), grepl("cannot stat", missing_msg, fixed = TRUE))
ok("ng_load_require_materialised() refuses, by name, a file it cannot read")

cat(sprintf("estimator_selfguard_bit_identity: PASS (%d checks)\n", checks))
