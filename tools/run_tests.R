# Single entry point for the whole test suite, used locally and by CI.
#
# The package has two tiers of tests and they check different things:
#
#   1. tests/testthat/  - the fast gate. Ships with the package and runs under
#      R CMD check against the INSTALLED package via library(). It answers
#      "does the artifact we publish actually work, through its public API?"
#
#   2. tests/*.R        - the deep harness. Repository-only; loads the source
#      tree with ng_load() so it can reach internals, and covers the science
#      (recombination kernels, optimizers, calibration, polyploid families,
#      external baselines). It answers "is the method still correct?"
#
# Neither subsumes the other, so CI runs both and so should you. Run this before
# opening a pull request; CI runs the same script, so a green run here is the
# same evidence the repository gate uses.
#
# Usage:
#   Rscript tools/run_tests.R              # both tiers
#   Rscript tools/run_tests.R --fast       # testthat only
#   Rscript tools/run_tests.R --harness    # the tests/*.R scripts only
#
# Environment:
#   NG_TEST_CPP=1   also cross-check the C++ kernel against the pure-R path
#   NGCD_SKIP_CPP=1 skip compiling the kernel in the harness (faster)

args <- commandArgs(trailingOnly = TRUE)
run_fast <- !("--harness" %in% args)
run_harness <- !("--fast" %in% args)

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("Run this from the repository root (no DESCRIPTION here).", call. = FALSE)
}

failures <- character(0)
started <- Sys.time()

if (run_fast) {
  cat("\n== fast gate: tests/testthat (installed package, public API) ==\n")
  if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("testthat is not installed. install.packages(\"testthat\")", call. = FALSE)
  }
  if (!requireNamespace("nextgenCrossDesign", quietly = TRUE)) {
    stop("nextgenCrossDesign is not installed. The fast gate tests the installed\n",
         "package, so build and install it first:\n",
         "  R CMD INSTALL .", call. = FALSE)
  }
  res <- testthat::test_dir(file.path(root, "tests", "testthat"),
                            package = "nextgenCrossDesign",
                            stop_on_failure = FALSE,
                            reporter = "summary")
  df <- as.data.frame(res)
  n_fail <- sum(df$failed) + sum(df$error)
  if (n_fail > 0) failures <- c(failures, sprintf("testthat: %d failing", n_fail))
}

if (run_harness) {
  cat("\n== deep harness: tests/*.R (source tree, internals + science) ==\n")
  scripts <- setdiff(
    list.files(file.path(root, "tests"), pattern = "[.]R$", full.names = TRUE),
    file.path(root, "tests", c("helper_load.R", "testthat.R"))
  )
  scripts <- scripts[!grepl("/testthat/", scripts, fixed = TRUE)]
  rscript <- file.path(R.home("bin"), "Rscript")
  for (s in scripts) {
    nm <- basename(s)
    t0 <- Sys.time()
    output <- suppressWarnings(system2(rscript, shQuote(s), stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
    if (identical(status, 0L)) {
      cat(sprintf("  PASS  %-52s %4ds\n", nm, secs))
    } else {
      cat(sprintf("  FAIL  %-52s %4ds\n", nm, secs))
      if (length(output)) cat(paste0("        ", output, collapse = "\n"), "\n")
      failures <- c(failures, paste0("harness: ", nm))
    }
  }
}

elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "mins")), 1)
cat(sprintf("\n== finished in %s min ==\n", elapsed))
if (length(failures)) {
  cat("FAILURES:\n"); cat(paste0("  - ", failures, collapse = "\n"), "\n")
  quit(status = 1L)
}
cat("all green\n")
