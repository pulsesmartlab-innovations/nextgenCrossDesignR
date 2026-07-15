ng_test_find_root <- function() {
  candidates <- c(
    Sys.getenv("NG_REPO_ROOT", unset = NA_character_),
    getwd(),
    file.path(getwd(), "nextgen_cross_design"),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", "nextgen_cross_design")
  )
  candidates <- unique(normalizePath(
    candidates[!is.na(candidates) & nzchar(candidates)],
    winslash = "/",
    mustWork = FALSE
  ))

  is_package_root <- vapply(candidates, function(path) {
    desc <- file.path(path, "DESCRIPTION")
    if (!file.exists(file.path(path, "R", "load.R")) || !file.exists(desc)) {
      return(FALSE)
    }
    any(grepl("^Package:\\s*nextgenCrossDesign\\s*$", readLines(desc, n = 20L, warn = FALSE)))
  }, logical(1L))

  roots <- candidates[is_package_root]
  if (!length(roots)) {
    stop("Could not locate nextgenCrossDesign root containing DESCRIPTION and R/load.R", call. = FALSE)
  }
  roots[[1L]]
}

root <- ng_test_find_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))

use_cpp <- isTRUE(get0("ng_test_use_cpp", ifnotfound = FALSE, inherits = FALSE))
ng_load(root, use_cpp = use_cpp, verbose = FALSE)
