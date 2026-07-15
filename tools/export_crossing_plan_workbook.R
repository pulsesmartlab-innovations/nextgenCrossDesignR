root_candidates <- unique(normalizePath(c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
), winslash = "/", mustWork = FALSE))
root_hits <- root_candidates[file.exists(file.path(root_candidates, "DESCRIPTION")) &
  file.exists(file.path(root_candidates, "R", "load.R"))]
if (!length(root_hits)) stop("Could not locate nextgenCrossDesign package root", call. = FALSE)
root <- root_hits[[1L]]
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

env_chr <- function(name, default = "") {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) default else value
}

env_int <- function(name, default) {
  value <- env_chr(name, "")
  if (!nzchar(value)) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (!is.finite(parsed) || parsed <= 0L) stop(name, " must be a positive integer", call. = FALSE)
  parsed
}

candidates_path <- env_chr("NG_CROSSING_PLAN_CANDIDATES")
if (!nzchar(candidates_path)) {
  stop("NG_CROSSING_PLAN_CANDIDATES must point to a candidate-line CSV", call. = FALSE)
}
criteria_path <- env_chr("NG_CROSSING_PLAN_CRITERIA")
out_path <- env_chr("NG_CROSSING_PLAN_OUT", file.path(root, "results", "crossing_plan.xlsx"))

candidates <- utils::read.csv(candidates_path, stringsAsFactors = FALSE, check.names = FALSE)
criteria <- if (nzchar(criteria_path)) {
  utils::read.csv(criteria_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  NULL
}

path <- ng_write_crossing_plan_workbook(
  output_path = out_path,
  candidates = candidates,
  criteria = criteria,
  top_n = env_int("NG_CROSSING_PLAN_TOP_N", 5L),
  cross_n = env_int("NG_CROSSING_PLAN_CROSS_N", 10L),
  title = env_chr("NG_CROSSING_PLAN_TITLE", NULL)
)
cat(path, "\n", sep = "")
