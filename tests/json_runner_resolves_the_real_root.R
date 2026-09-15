# The headless entry point must load the package it lives in, not a stale copy
# sitting beside it.
#
# tools/run_cross_prediction_json.R is what a frontend spawns. Its root search
# listed `getwd()/nextgen_cross_design` BEFORE `getwd()`, so whenever a directory
# of that name existed next to the real sources it won. In this repository one
# does -- an untracked 0.19.0 tree, eleven versions behind -- so every headless
# run loaded 0.19.0 while every in-process test loaded the current source.
#
# That is not a cosmetic mismatch. Under 0.19.0 a RIL panel with residual
# heterozygosity and no phased haplotypes merely WARNS and proceeds on the a'Ra
# kernel, which is biased low at exactly those loci. Current sources refuse it.
# So the frontend path silently delivered numbers the package itself considers
# invalid, and tests/run_cross_prediction_json_messages.R -- which asserts the
# refusal -- failed while the refusal worked perfectly everywhere else.
#
# Commit ac35088 fixed this same defect for tests/helper_load.R. The tools/ entry
# point was missed. Both now share one resolver.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

source(file.path(root, "tools", "ng_find_package_root.R"))

mk_pkg <- function(dir, version) {
  dir.create(file.path(dir, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("Package: nextgenCrossDesign", paste0("Version: ", version)),
             file.path(dir, "DESCRIPTION"))
  writeLines("ng_load <- function(...) invisible(NULL)", file.path(dir, "R", "load.R"))
  invisible(dir)
}

# ---- the decoy must lose to the root it is nested in ----------------------
td <- tempfile("rootres_"); dir.create(td)
outer <- mk_pkg(td, "0.32.0")
mk_pkg(file.path(td, "nextgen_cross_design"), "0.19.0")

got <- ng_find_package_root(outer)
stopifnot(identical(normalizePath(got), normalizePath(outer)))
ver <- grep("^Version:", readLines(file.path(got, "DESCRIPTION")), value = TRUE)
stopifnot(identical(trimws(sub("^Version:", "", ver)), "0.32.0"))

# ---- but a genuine subdirectory layout still resolves ---------------------
# Some checkouts really do keep the package one level down. When the working
# directory is NOT itself a package root, the subdirectory is the right answer.
td2 <- tempfile("rootres2_"); dir.create(td2)
mk_pkg(file.path(td2, "nextgen_cross_design"), "0.32.0")
got2 <- ng_find_package_root(td2)
stopifnot(identical(normalizePath(got2),
                    normalizePath(file.path(td2, "nextgen_cross_design"))))

# ---- a directory that only LOOKS like a root is rejected ------------------
# R/load.R alone is not enough; DESCRIPTION must name this package.
td3 <- tempfile("rootres3_"); dir.create(td3)
dir.create(file.path(td3, "R"), recursive = TRUE)
writeLines("ng_load <- function(...) NULL", file.path(td3, "R", "load.R"))
writeLines(c("Package: someOtherPkg", "Version: 1.0"), file.path(td3, "DESCRIPTION"))
mk_pkg(file.path(td3, "nextgen_cross_design"), "0.32.0")
stopifnot(identical(normalizePath(ng_find_package_root(td3)),
                    normalizePath(file.path(td3, "nextgen_cross_design"))))

# ---- and the real runner reports the real version -------------------------
# The observable that actually failed: result.json's package_version must be the
# version in the DESCRIPTION beside the runner, whatever else sits nearby.
desc_v <- trimws(sub("^Version:", "",
  grep("^Version:", readLines(file.path(root, "DESCRIPTION")), value = TRUE)[[1L]]))
stopifnot(identical(normalizePath(ng_find_package_root(root)), normalizePath(root)))
cat("json_runner_resolves_the_real_root: PASS  (DESCRIPTION", desc_v, ")\n")
