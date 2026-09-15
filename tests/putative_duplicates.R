# Repository root FIRST. Checking `nextgen_cross_design/` ahead of getwd() made
# these tests load a STALE 0.19.0 copy of the package that sits in the working
# tree under exactly that name -- so they validated a package eleven versions old
# while appearing to cover the current one. The ones that failed were the lucky
# case; the ones that passed gave false assurance. ng_load() already resolves in
# this order; only these hand-rolled preambles inverted it.
root_candidates <- c(
  getwd(),
  file.path(".."),
  file.path(getwd(), "nextgen_cross_design"),
  file.path("..", "nextgen_cross_design")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

geno <- matrix(
  c(
    0, 1, 2, 0, 2, 1,
    0, 1, 2, 0, 2, 1,
    2, 2, 0, 1, 0, 1,
    0, 0, 1, 2, 1, 2
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(c("P01", "P01_copy", "P03", "P04"), paste0("M", 1:6))
)

dup <- ng_detect_putative_duplicates(
  geno,
  duplicate_threshold = 0.995,
  maf_min = 0,
  max_missing_prop = 1,
  min_compared_markers = 4
)

stopifnot(inherits(dup, "ng_putative_duplicates"))
stopifnot(nrow(dup$pairs) == 1L)
stopifnot(identical(dup$pairs$parent1[[1]], "P01"))
stopifnot(identical(dup$pairs$parent2[[1]], "P01_copy"))
stopifnot(dup$pairs$similarity[[1]] >= 0.995)
stopifnot(dup$summary$value[dup$summary$metric == "duplicate_pairs"] == 1L)
stopifnot(length(dup$clusters) == 1L)

tmp <- tempfile("ng_putative_duplicates_", fileext = ".png")
ng_plot_putative_duplicates(dup, output_path = tmp, max_ids = 20)
stopifnot(file.exists(tmp))
stopifnot(file.info(tmp)$size > 1000)

clean <- ng_detect_putative_duplicates(
  geno[c("P01", "P03", "P04"), ],
  duplicate_threshold = 0.995,
  min_compared_markers = 4
)
stopifnot(nrow(clean$pairs) == 0L)
tmp_clean <- tempfile("ng_putative_duplicates_clean_", fileext = ".png")
ng_plot_putative_duplicates(clean, output_path = tmp_clean)
stopifnot(file.exists(tmp_clean))
stopifnot(file.info(tmp_clean)$size > 1000)

cat("putative duplicate tests passed\n")
