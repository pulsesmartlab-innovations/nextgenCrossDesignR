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

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (!is.finite(out)) as.integer(default) else out
}

env_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) as.numeric(default) else out
}

env_nullable_int <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(NULL)
  out <- suppressWarnings(as.integer(value))
  if (!is.finite(out)) NULL else out
}

env_nullable_num <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(NULL)
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) NULL else out
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

prefix <- env_chr("NG_MULTITRAIT_CROP_GRID_PREFIX", "multitrait_crop_grid")
output_dir <- env_chr("NG_MULTITRAIT_CROP_GRID_OUTPUT_DIR", file.path(dirname(root), "results"))
scenarios <- split_csv(env_chr(
  "NG_MULTITRAIT_CROP_GRID_SCENARIOS",
  "compact_selfing,maize_like,cassava_diploid"
))
parent_sizes <- env_int_vec("NG_MULTITRAIT_CROP_GRID_PARENT_SIZES", c(12L, 20L))
reps <- env_int("NG_MULTITRAIT_CROP_GRID_REPS", 1L)
n_crosses <- env_int("NG_MULTITRAIT_CROP_GRID_CROSSES", 4L)
realized_progeny <- env_int("NG_MULTITRAIT_CROP_GRID_REALIZED_PROGENY", 8L)
seed <- env_int("NG_MULTITRAIT_CROP_GRID_SEED", 1L)
allocator <- tolower(env_chr("NG_MULTITRAIT_CROP_GRID_ALLOCATOR", "topn"))
prediction_noise <- env_num("NG_MULTITRAIT_CROP_GRID_PREDICTION_NOISE", 0.20)
alphasimr_threads <- env_int("NG_MULTITRAIT_CROP_GRID_ALPHASIMR_THREADS", 1L)
ocs_lambda_group <- env_num("NG_MULTITRAIT_CROP_GRID_OCS_LAMBDA_GROUP", 0.05)
if (!is.finite(ocs_lambda_group) || ocs_lambda_group < 0) ocs_lambda_group <- 0.05
ocs_lambda_mating <- env_num("NG_MULTITRAIT_CROP_GRID_OCS_LAMBDA_MATING", 0)
if (!is.finite(ocs_lambda_mating) || ocs_lambda_mating < 0) ocs_lambda_mating <- 0
methods <- split_csv(env_chr(
  "NG_MULTITRAIT_CROP_GRID_METHODS",
  paste(ng_multitrait_validation_default_methods(), collapse = ",")
))

n_founders <- env_nullable_int("NG_MULTITRAIT_CROP_GRID_N_FOUNDERS")
n_chr <- env_nullable_int("NG_MULTITRAIT_CROP_GRID_N_CHR")
seg_sites <- env_nullable_int("NG_MULTITRAIT_CROP_GRID_SEG_SITES")
snp_per_chr <- env_nullable_int("NG_MULTITRAIT_CROP_GRID_SNP_PER_CHR")
qtl_per_chr <- env_nullable_int("NG_MULTITRAIT_CROP_GRID_QTL_PER_CHR")
genome_length_m <- env_nullable_num("NG_MULTITRAIT_CROP_GRID_GENOME_LENGTH_M")

cat("Running AlphaSimR multi-trait crop validation grid\n")
cat("  prefix: ", prefix, "\n", sep = "")
cat("  output_dir: ", output_dir, "\n", sep = "")
cat("  scenarios: ", paste(scenarios, collapse = ","), "\n", sep = "")
cat("  parent_sizes: ", paste(parent_sizes, collapse = ","), "\n", sep = "")
cat("  reps: ", reps, "\n", sep = "")
cat("  crosses: ", n_crosses, "\n", sep = "")
cat("  realized_progeny: ", realized_progeny, "\n", sep = "")
cat("  seed: ", seed, "\n", sep = "")
cat("  allocator: ", allocator, "\n", sep = "")
cat("  methods: ", paste(methods, collapse = ","), "\n", sep = "")
cat("  prediction_noise: ", prediction_noise, "\n", sep = "")
cat("  ocs_lambda_group: ", ocs_lambda_group, "\n", sep = "")
cat("  ocs_lambda_mating: ", ocs_lambda_mating, "\n", sep = "")

result <- ng_run_multitrait_crop_validation_grid(
  scenarios = scenarios,
  parent_sizes = parent_sizes,
  reps = reps,
  n_crosses = n_crosses,
  realized_progeny = realized_progeny,
  seed = seed,
  methods = methods,
  allocator = allocator,
  n_founders = n_founders,
  n_chr = n_chr,
  seg_sites = seg_sites,
  snp_per_chr = snp_per_chr,
  qtl_per_chr = qtl_per_chr,
  genome_length_m = genome_length_m,
  prediction_noise = prediction_noise,
  alphasimr_threads = alphasimr_threads,
  ocs_lambda_group = ocs_lambda_group,
  ocs_lambda_mating = ocs_lambda_mating,
  output_dir = output_dir,
  prefix = prefix
)

print(result$winner_summary)
cat("Wrote AlphaSimR multi-trait crop validation outputs with prefix: ", prefix, "\n", sep = "")
