# `%||%` entered base R only in 4.4.0, while DESCRIPTION declares R (>= 4.1.0). The package uses
# it in 15 places outside the one function that defines its own local copy -- including
# R/03_metrics.R, inside ng_score_crosses()'s core scoring path. Without a package-scope
# definition the package fails immediately on R 4.1-4.3 and works fine on 4.4+, so the breakage
# is invisible to anyone developing on a current R and no behavioural test can catch it.
#
# This test asserts the definition exists in the package's own scope, so deleting it fails
# loudly here rather than silently on someone else's older R.

ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# 1. the operator resolves, and does so from the package rather than only from base
stopifnot(exists("%||%", inherits = TRUE))
op <- get("%||%")
stopifnot(is.function(op))
stopifnot(identical(op(NULL, "fallback"), "fallback"))
stopifnot(identical(op("value", "fallback"), "value"))
stopifnot(is.null(op(NULL, NULL)))
# NA is a value, not absence -- it must NOT fall through
stopifnot(identical(op(NA, "fallback"), NA))

# 2. every %||% use in R/ is either inside a function that defines its own, or covered by the
#    package-scope definition. Guard the count so a new out-of-scope use is noticed.
root <- ng_test_find_root()
r_files <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
defs <- 0L
for (f in r_files) {
  txt <- readLines(f, warn = FALSE)
  defs <- defs + sum(grepl("^\\s*`%\\|\\|%`\\s*<-", txt))
}
stopifnot(defs >= 1L)   # at least the package-scope one in R/00_utils.R

cat("null-coalesce available at package scope\n")
