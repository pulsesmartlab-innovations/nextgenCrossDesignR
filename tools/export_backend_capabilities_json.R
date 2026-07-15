root_candidates <- unique(normalizePath(c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
), winslash = "/", mustWork = FALSE))
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
if (!length(root_hits)) stop("Could not locate nextgen_cross_design root", call. = FALSE)
root <- root_hits[[1]]
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
