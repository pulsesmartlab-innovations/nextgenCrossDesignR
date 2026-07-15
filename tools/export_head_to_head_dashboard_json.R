find_project_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(normalizePath(c(
    file.path(start, "nextgen_cross_design"),
    file.path(dirname(start), "nextgen_cross_design"),
    start,
    dirname(start),
    file.path(dirname(dirname(start)), "nextgen_cross_design"),
    dirname(dirname(start))
  ), winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "load.R"))) return(candidate)
  }
  stop("Could not locate nextgen_cross_design root containing R/load.R", call. = FALSE)
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

root <- find_project_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

prefix <- env_chr("NG_HEAD_TO_HEAD_PREFIX", "head_to_head_benchmark")
output_dir <- env_chr("NG_HEAD_TO_HEAD_OUTPUT_DIR", file.path(dirname(root), "results"))
out_path <- env_chr(
  "NG_HEAD_TO_HEAD_DASHBOARD_JSON_OUT",
  file.path(output_dir, paste0(prefix, "_dashboard.json"))
)

cat("Exporting head-to-head dashboard JSON\n")
cat("  prefix: ", prefix, "\n", sep = "")
cat("  output_dir: ", output_dir, "\n", sep = "")
cat("  out: ", out_path, "\n", sep = "")

written <- ng_write_head_to_head_dashboard_json(
  output_path = out_path,
  results_dir = output_dir,
  prefix = prefix
)

cat("Wrote head-to-head dashboard JSON: ", written, "\n", sep = "")
