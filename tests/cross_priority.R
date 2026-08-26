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

scores <- data.frame(
  parent1 = sprintf("P%03d", 1:100),
  parent2 = sprintf("Q%03d", 1:100),
  multi_trait_score = seq(100, 1),
  pair_kinship = seq(-0.5, 0.5, length.out = 100),
  multi_trait_threshold_violation = 0,
  stringsAsFactors = FALSE
)

ranked <- ng_rank_cross_priority(scores)
stopifnot(nrow(ranked) == 100L)
stopifnot(all(c("priority_tier", "priority_rank", "priority_index") %in% names(ranked)))
counts <- table(ranked$priority_tier)
stopifnot(identical(as.integer(counts[c("highly_priority", "priority", "medium_priority", "low_priority")]),
                    c(10L, 25L, 35L, 30L)))
stopifnot(identical(as.character(ranked$priority_tier[[1]]), "highly_priority"))
stopifnot(identical(ranked$priority_rank, seq_len(100L)))

ties <- data.frame(
  parent1 = c("A", "C"),
  parent2 = c("B", "D"),
  multi_trait_score = c(1, 1),
  pair_kinship = c(0.45, -0.25),
  stringsAsFactors = FALSE
)
ranked_ties <- ng_rank_cross_priority(ties, kinship_weight = 0.5)
stopifnot(identical(ranked_ties$parent1[[1]], "C"))

custom <- ng_rank_cross_priority(scores, breaks = c(0.20, 0.50, 0.80, 1))
custom_counts <- table(custom$priority_tier)
stopifnot(identical(as.integer(custom_counts[c("highly_priority", "priority", "medium_priority", "low_priority")]),
                    c(20L, 30L, 30L, 20L)))

set.seed(7101)
parents <- sprintf("D%03d", seq_len(50L))
geno <- matrix(2L * stats::rbinom(50L * 18L, size = 1L, prob = 0.5), nrow = 50L,
               dimnames = list(parents, sprintf("m%02d", seq_len(18L))))
y <- stats::setNames(rowSums(geno[, 1:4, drop = FALSE]) + stats::rnorm(50L, sd = 0.1), parents)
marker_map <- data.frame(marker = colnames(geno), chr = 1L,
                         pos_cm = seq(0, 85, length.out = ncol(geno)))
design <- suppressWarnings(ng_design_crosses(
  geno = geno,
  y = y,
  marker_map = marker_map,
  use_cpp = FALSE,
  parent_type = "inbred"
))
stopifnot(nrow(design$plan) == 100L)

cat("cross priority tests passed\n")
