root_candidates <- c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

scored <- data.frame(
  parent1 = rep(sprintf("P%02d", seq_len(10L)), each = 2L),
  parent2 = rep(sprintf("Q%02d", seq_len(10L)), times = 2L),
  pair_kinship = seq(-0.55, 0.65, length.out = 20L),
  multi_trait_score = seq(1.05, -0.35, length.out = 20L),
  multi_trait_threshold_violation = 0,
  stringsAsFactors = FALSE
)
selected <- ng_rank_cross_priority(scored[seq_len(8L), , drop = FALSE],
                                   breaks = c(0.25, 0.50, 0.75, 1))

png_path <- tempfile("priority_score_vs_kinship_", fileext = ".png")
written <- ng_plot_priority_score_vs_kinship(
  scored = scored,
  selected = selected,
  output_path = png_path
)
stopifnot(identical(normalizePath(png_path, winslash = "/", mustWork = TRUE), written))
stopifnot(file.exists(written))
stopifnot(file.info(written)$size > 1000)

pdf_path <- tempfile("priority_score_vs_kinship_", fileext = ".pdf")
written_pdf <- ng_plot_priority_score_vs_kinship(
  scored = scored,
  selected = selected,
  output_path = pdf_path
)
stopifnot(file.exists(written_pdf))
stopifnot(file.info(written_pdf)$size > 1000)

current_device_path <- tempfile("priority_score_vs_kinship_current_device_", fileext = ".png")
grDevices::png(current_device_path, width = 900, height = 560, res = 120)
result <- tryCatch({
  ng_plot_priority_score_vs_kinship(scored = scored, selected = selected)
}, finally = {
  grDevices::dev.off()
})
stopifnot(is.null(result))
stopifnot(file.exists(current_device_path))
stopifnot(file.info(current_device_path)$size > 1000)

cat("cross priority plot tests passed\n")
