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

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (!is.finite(out)) as.integer(default) else out
}

env_chr <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

split_csv <- function(x) {
  x <- unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE)
  x <- trimws(x)
  x[nzchar(x) & !is.na(x)]
}

root <- find_project_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

prefix <- env_chr("NG_MULTITRAIT_VALIDATION_PREFIX", "multitrait_validation")
output_dir <- env_chr("NG_MULTITRAIT_VALIDATION_OUTPUT_DIR", file.path(dirname(root), "results"))
n_parents <- env_int("NG_MULTITRAIT_VALIDATION_PARENTS", 20L)
n_crosses <- env_int("NG_MULTITRAIT_VALIDATION_CROSSES", 6L)
seed <- env_int("NG_MULTITRAIT_VALIDATION_SEED", 1L)
allocator <- tolower(env_chr("NG_MULTITRAIT_VALIDATION_ALLOCATOR", "topn"))
methods <- split_csv(env_chr(
  "NG_MULTITRAIT_VALIDATION_METHODS",
  paste(ng_multitrait_validation_default_methods(), collapse = ",")
))

cat("Running multi-trait validation\n")
cat("  prefix: ", prefix, "\n", sep = "")
cat("  output_dir: ", output_dir, "\n", sep = "")
cat("  parents: ", n_parents, "\n", sep = "")
cat("  crosses: ", n_crosses, "\n", sep = "")
cat("  seed: ", seed, "\n", sep = "")
cat("  allocator: ", allocator, "\n", sep = "")
cat("  methods: ", paste(methods, collapse = ","), "\n", sep = "")

result <- ng_run_multitrait_validation(
  n_parents = n_parents,
  n_crosses = n_crosses,
  seed = seed,
  methods = methods,
  allocator = allocator,
  output_dir = output_dir,
  prefix = prefix
)

print(result$summary)
cat("Wrote multi-trait validation outputs with prefix: ", prefix, "\n", sep = "")
