root_candidates <- unique(normalizePath(c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
), winslash = "/", mustWork = FALSE))
root_hits <- root_candidates[file.exists(file.path(root_candidates, "DESCRIPTION")) &
  file.exists(file.path(root_candidates, "R"))]
if (!length(root_hits)) stop("Could not locate nextgenCrossDesign package root", call. = FALSE)
root <- root_hits[[1]]

r_exe <- file.path(R.home("bin"), "R.exe")
if (!file.exists(r_exe)) r_exe <- file.path(R.home("bin"), "R")

stage_parent <- tempfile("ngcd_r_pkg_")
dir.create(stage_parent, recursive = TRUE)
stage <- file.path(stage_parent, "nextgenCrossDesign")
dir.create(stage)

copy_item <- function(name) {
  src <- file.path(root, name)
  if (!file.exists(src)) return(invisible(FALSE))
  dst <- file.path(stage, name)
  if (dir.exists(src)) {
    dir.create(dst, recursive = TRUE)
    file.copy(list.files(src, all.files = TRUE, no.. = TRUE, full.names = TRUE), dst, recursive = TRUE)
  } else {
    file.copy(src, dst)
  }
  invisible(TRUE)
}

for (name in c("DESCRIPTION", "LICENSE", "NAMESPACE", ".Rbuildignore", "README.md", "R", "src", "inst", "man")) {
  copy_item(name)
}

old <- setwd(stage_parent)
on.exit(setwd(old), add = TRUE)

message("Staged R package at: ", stage)
build_status <- system2(r_exe, c("CMD", "build", "nextgenCrossDesign"))
if (!identical(build_status, 0L)) stop("R CMD build failed with status ", build_status, call. = FALSE)

tarball <- list.files(stage_parent, pattern = "^nextgenCrossDesign_.*[.]tar[.]gz$", full.names = TRUE)
if (!length(tarball)) stop("R CMD build did not produce a source tarball", call. = FALSE)
check_status <- system2(r_exe, c("CMD", "check", "--no-manual", "--no-build-vignettes", tarball[[1L]]))
if (!identical(check_status, 0L)) stop("R CMD check failed with status ", check_status, call. = FALSE)

message("R backend package check passed: ", tarball[[1L]])
