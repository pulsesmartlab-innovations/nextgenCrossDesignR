ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

if (!exists("ng_dh_recomb_pairs_cpp", mode = "function", inherits = TRUE)) {
  message("C++ kernel unavailable; skipping consistency check.")
  quit(status = 0)
}

set.seed(22)
n <- 12
m <- 80
ids <- paste0("P", seq_len(n))
markers <- paste0("M", seq_len(m))
geno <- matrix(rbinom(n * m, 2, 0.45), n, m, dimnames = list(ids, markers))
beta <- stats::rnorm(m, sd = 0.05)
beta_var <- stats::runif(m, min = 0, max = 0.002)
marker_map <- data.frame(
  marker = markers,
  chr = rep(seq_len(4), each = m / 4),
  pos_cm = rep(seq(0, 100, length.out = m / 4), 4)
)
pairs <- ng_make_pairs(ids)
marker_map <- ng_prepare_marker_map(marker_map, markers)
sorted <- ng_sort_by_map(geno, beta, beta_var, marker_map)

cpp <- ng_dh_recomb_variance_pairs(sorted$geno, sorted$effects, sorted$beta_var,
                                   sorted$marker_map, ids, pairs, use_cpp = TRUE)
r <- ng_dh_recomb_variance_pairs(sorted$geno, sorted$effects, sorted$beta_var,
                                 sorted$marker_map, ids, pairs, use_cpp = FALSE)
stopifnot(max(abs(cpp$vpm - r$vpm)) < 1e-8)
stopifnot(max(abs(cpp$pmv - r$pmv)) < 1e-8)
message("C++ and R recombination kernels agree.")
