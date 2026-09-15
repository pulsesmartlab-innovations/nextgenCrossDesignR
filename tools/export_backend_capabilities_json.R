# Resolved through the one shared resolver. This used to probe
# getwd()/nextgen_cross_design BEFORE getwd(), so an untracked stale copy of the
# package beside the real sources won -- the defect that had two harness tests
# silently validating a 0.19.0 package while appearing green.
local({
  cands <- file.path(c(".", "..", "../..", "nextgen_cross_design",
                       "../nextgen_cross_design"), "tools", "ng_find_package_root.R")
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot locate tools/ng_find_package_root.R", call. = FALSE)
  source(hit[[1L]], local = FALSE)
})
root_candidates <- ng_find_package_root(getwd())
root <- root_candidates
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

out <- env_chr("NG_BACKEND_CAPABILITIES_OUT", file.path(root, "results", "backend_capabilities.json"))
path <- ng_write_backend_capability_registry_json(output_path = out)
cat("Wrote backend capability registry: ", path, "\n", sep = "")
