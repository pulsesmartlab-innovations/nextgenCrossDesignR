# One resolver for "where is the package I am part of?".
#
# Every entry point that loads the source tree has to answer this, and two of them
# answered it differently. tests/helper_load.R tried getwd() before
# getwd()/nextgen_cross_design and verified DESCRIPTION; tools/run_cross_prediction_json.R
# tried the subdirectory FIRST and verified nothing. A working tree with an
# untracked nextgen_cross_design/ beside the real sources therefore ran the
# headless frontend path against that copy -- in this repository, a version eleven
# releases behind, whose RIL het-parent handling merely warns where current sources
# refuse. In-process tests loaded the current source and passed; the headless path
# loaded the old one and delivered numbers the package considers invalid.
#
# THE RULE: a directory that is ITSELF the package wins over any copy nested inside
# it. Nesting means the inner one is a vendored or stale duplicate; if the outer
# directory is a package root, it is the one being worked on. Only when the
# starting directory is not a package root does a subdirectory become the answer.
#
# A candidate qualifies only if it has R/load.R AND a DESCRIPTION naming this
# package -- R/load.R alone matches any project that happens to use the filename.
ng_find_package_root <- function(start = getwd()) {
  candidates <- c(
    Sys.getenv("NG_REPO_ROOT", unset = NA_character_),
    start,
    file.path(start, "nextgen_cross_design"),
    file.path(start, ".."),
    file.path(start, "..", "nextgen_cross_design")
  )
  candidates <- unique(normalizePath(
    candidates[!is.na(candidates) & nzchar(candidates)],
    winslash = "/", mustWork = FALSE))

  is_root <- vapply(candidates, function(path) {
    desc <- file.path(path, "DESCRIPTION")
    if (!file.exists(file.path(path, "R", "load.R")) || !file.exists(desc)) return(FALSE)
    any(grepl("^Package:\\s*nextgenCrossDesign\\s*$",
              readLines(desc, n = 20L, warn = FALSE)))
  }, logical(1L))

  roots <- candidates[is_root]
  if (!length(roots)) {
    stop("Could not locate nextgenCrossDesign root containing DESCRIPTION and R/load.R",
         call. = FALSE)
  }
  roots[[1L]]
}
