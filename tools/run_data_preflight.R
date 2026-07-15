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

env_chr <- function(name, default = "") {
  value <- Sys.getenv(name, unset = default)
  if (!nzchar(value)) default else value
}

env_num <- function(name, default = NA_real_) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) default else value
}

env_logical <- function(name, default = FALSE) {
  value <- tolower(trimws(Sys.getenv(name, unset = if (isTRUE(default)) "true" else "false")))
  value %in% c("1", "true", "yes", "y")
}

read_table <- function(path) {
  if (!nzchar(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

out <- env_chr("NG_PREFLIGHT_OUT", file.path(root, "results", "data_preflight.json"))
path <- ng_write_data_preflight_json(
  output_path = out,
  geno = read_table(env_chr("NG_PREFLIGHT_GENO")),
  phenotype = read_table(env_chr("NG_PREFLIGHT_PHENOTYPE")),
  candidate_pairs = read_table(env_chr("NG_PREFLIGHT_CANDIDATE_PAIRS")),
  trait_spec = read_table(env_chr("NG_PREFLIGHT_TRAIT_SPEC")),
  marker_map = read_table(env_chr("NG_PREFLIGHT_MARKER_MAP")),
  ploidy = env_num("NG_PREFLIGHT_PLOIDY", NA_real_),
  putative_duplicate_check = env_logical("NG_PREFLIGHT_PUTATIVE_DUPLICATES", FALSE),
  duplicate_threshold = env_num("NG_PREFLIGHT_DUPLICATE_THRESHOLD", 0.995),
  duplicate_maf_min = env_num("NG_PREFLIGHT_DUPLICATE_MAF_MIN", 0),
  duplicate_max_missing_prop = env_num("NG_PREFLIGHT_DUPLICATE_MAX_MISSING_PROP", 1),
  duplicate_min_compared_markers = as.integer(env_num("NG_PREFLIGHT_DUPLICATE_MIN_COMPARED_MARKERS", 100)),
  putative_duplicate_action = env_chr("NG_PREFLIGHT_DUPLICATE_ACTION", "report"),
  putative_duplicate_return_similarity = env_logical("NG_PREFLIGHT_DUPLICATE_RETURN_SIMILARITY", FALSE)
)
cat("Wrote data preflight JSON: ", path, "\n", sep = "")
