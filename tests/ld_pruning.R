ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(20260508)
geno <- matrix(
  sample(c(0, 1, 2, NA), 36L * 10L, replace = TRUE, prob = c(0.30, 0.35, 0.30, 0.05)),
  nrow = 36L
)
geno[, 2] <- geno[, 1]
geno[, 5] <- geno[, 4]
geno[, 9] <- 0
geno[1:2, 9] <- 1
rownames(geno) <- paste0("P", seq_len(nrow(geno)))
colnames(geno) <- paste0("M", seq_len(ncol(geno)))
marker_map <- data.frame(
  marker = colnames(geno), chr = 1L, pos_cm = seq(0, by = 1, length.out = ncol(geno)),
  stringsAsFactors = FALSE
)

r_pruned <- ng_ld_prune_markers(
  geno,
  window = 3L,
  r2_threshold = 0.8,
  maf_threshold = 0.03,
  backend = "r",
  marker_map = marker_map
)

stopifnot(is.character(r_pruned$keep_markers))
stopifnot(!"M9" %in% r_pruned$keep_markers)
stopifnot(sum(c("M1", "M2") %in% r_pruned$keep_markers) == 1L)
stopifnot(sum(c("M4", "M5") %in% r_pruned$keep_markers) == 1L)
stopifnot(r_pruned$report$markers_before == ncol(geno))
stopifnot(r_pruned$report$markers_after == length(r_pruned$keep_markers))
stopifnot(r_pruned$report$low_maf_removed >= 1L)

if (exists("ng_ld_prune_graph_cpp", mode = "function", inherits = TRUE)) {
  cpp_pruned <- ng_ld_prune_markers(
    geno,
    window = 3L,
    r2_threshold = 0.8,
    maf_threshold = 0.03,
    backend = "cpp",
    marker_map = marker_map
  )
  stopifnot(identical(cpp_pruned$keep_markers, r_pruned$keep_markers))
  stopifnot(identical(cpp_pruned$report$backend, "cpp"))
}

filtered <- ng_ld_prune_geno(
  geno,
  window = 3L,
  r2_threshold = 0.8,
  maf_threshold = 0.03,
  backend = "r",
  marker_map = marker_map
)
stopifnot(identical(colnames(filtered), r_pruned$keep_markers))
stopifnot(is.data.frame(attr(filtered, "ld_pruning_report")))

y <- stats::setNames(stats::rnorm(nrow(geno)), rownames(geno))
geno_design <- geno
geno_design[!is.finite(geno_design)] <- 0
geno_design <- 2L * (geno_design > 0)
dimnames(geno_design) <- dimnames(geno)
geno_design[, "M9"] <- 0L
plan <- suppressWarnings(ng_design_crosses(
  geno = geno_design,
  y = y,
  marker_map = marker_map,
  n_crosses = 4L,
  max_crosses_per_parent = 2L,
  ld_pruning = TRUE,
  ld_window = 3L,
  ld_r2_threshold = 0.8,
  ld_maf_threshold = 0.03,
  ld_backend = "r",
  use_cpp = FALSE,
  # This fixture tests graph pruning, not residual-heterozygous RIL variance.
  # Use the explicit inbred declaration after the pruning contract is checked.
  parent_type = "inbred"
))
stopifnot(is.data.frame(plan$ld_pruning_report))
stopifnot(plan$ld_pruning_report$markers_after == length(plan$effects$beta))
stopifnot(!"M9" %in% names(plan$effects$beta))

cat("LD pruning tests passed\n")
