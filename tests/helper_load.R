# Shared with tools/run_cross_prediction_json.R -- the two used to resolve the root
# differently, and the headless path lost to a stale nextgen_cross_design/ copy.
local({
  cands <- file.path(c(".", "..", "../..", "nextgen_cross_design",
                       "../nextgen_cross_design"), "tools", "ng_find_package_root.R")
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot locate tools/ng_find_package_root.R", call. = FALSE)
  source(hit[[1L]])
})
ng_test_find_root <- function() ng_find_package_root(getwd())

root <- ng_test_find_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))

use_cpp <- isTRUE(get0("ng_test_use_cpp", ifnotfound = FALSE, inherits = FALSE))
ng_load(root, use_cpp = use_cpp, verbose = FALSE)
