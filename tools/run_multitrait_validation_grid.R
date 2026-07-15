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

env_int_vec <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(split_csv(value)))
  out <- out[is.finite(out)]
  if (!length(out)) as.integer(default) else out
}

root <- find_project_root()
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

prefix <- env_chr("NG_MULTITRAIT_GRID_PREFIX", "multitrait_grid")
output_dir <- env_chr("NG_MULTITRAIT_GRID_OUTPUT_DIR", file.path(dirname(root), "results"))
parent_sizes <- env_int_vec("NG_MULTITRAIT_GRID_PARENT_SIZES", c(20L, 40L, 60L, 80L))
reps <- env_int("NG_MULTITRAIT_GRID_REPS", 3L)
n_crosses <- env_int("NG_MULTITRAIT_GRID_CROSSES", 6L)
seed <- env_int("NG_MULTITRAIT_GRID_SEED", 1L)
allocator <- tolower(env_chr("NG_MULTITRAIT_GRID_ALLOCATOR", "topn"))
ocs_lambda_group <- suppressWarnings(as.numeric(env_chr("NG_MULTITRAIT_GRID_OCS_LAMBDA_GROUP", "0.05")))
if (!is.finite(ocs_lambda_group) || ocs_lambda_group < 0) ocs_lambda_group <- 0.05
ocs_lambda_mating <- suppressWarnings(as.numeric(env_chr("NG_MULTITRAIT_GRID_OCS_LAMBDA_MATING", "0")))
if (!is.finite(ocs_lambda_mating) || ocs_lambda_mating < 0) ocs_lambda_mating <- 0
methods <- split_csv(env_chr(
  "NG_MULTITRAIT_GRID_METHODS",
  paste(ng_multitrait_validation_default_methods(), collapse = ",")
))

cat("Running replicated multi-trait validation grid\n")
cat("  prefix: ", prefix, "\n", sep = "")
cat("  output_dir: ", output_dir, "\n", sep = "")
cat("  parent_sizes: ", paste(parent_sizes, collapse = ","), "\n", sep = "")
cat("  reps: ", reps, "\n", sep = "")
cat("  crosses: ", n_crosses, "\n", sep = "")
cat("  seed: ", seed, "\n", sep = "")
cat("  allocator: ", allocator, "\n", sep = "")
cat("  ocs_lambda_group: ", ocs_lambda_group, "\n", sep = "")
cat("  ocs_lambda_mating: ", ocs_lambda_mating, "\n", sep = "")
cat("  methods: ", paste(methods, collapse = ","), "\n", sep = "")

result <- ng_run_multitrait_validation_grid(
  parent_sizes = parent_sizes,
  reps = reps,
  n_crosses = n_crosses,
  seed = seed,
  methods = methods,
  allocator = allocator,
  ocs_lambda_group = ocs_lambda_group,
  ocs_lambda_mating = ocs_lambda_mating,
  output_dir = output_dir,
  prefix = prefix
)

print(result$winner_summary)
cat("Wrote replicated multi-trait validation grid outputs with prefix: ", prefix, "\n", sep = "")
