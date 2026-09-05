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

# --- INTEGRATION: check_violation as a 4th, OFF-by-default priority component ----------------
scores_chk <- scores
scores_chk$check_violation <- rep(c(0L, 1L), length.out = 100L)

# check_weight = 0 (the default) -> tiers/ranks IDENTICAL to a run with no check_violation column
# at all: an existing run's tiers must never move just because the column now exists.
ranked_no_col <- ng_rank_cross_priority(scores)
ranked_col_off <- ng_rank_cross_priority(scores_chk)
stopifnot(identical(ranked_no_col$priority_rank, ranked_col_off$priority_rank))
stopifnot(identical(as.character(ranked_no_col$priority_tier), as.character(ranked_col_off$priority_tier)))
stopifnot(identical(ranked_no_col$priority_index, ranked_col_off$priority_index))
# reported (for inspection) even at weight 0, but contributes nothing to priority_index/rank/tier
stopifnot("priority_check_component" %in% names(ranked_col_off))
stopifnot(identical(ranked_col_off$priority_index, ranked_no_col$priority_index))

# a positive check_weight must NEVER veto (row count / evaluability unchanged), but must shift an
# otherwise-identical pair so the one failing checks ranks lower.
pair <- data.frame(
  parent1 = c("X1", "X2"), parent2 = c("Y1", "Y2"),
  multi_trait_score = c(5, 5),           # tied score
  pair_kinship = c(0, 0),                # tied kinship
  check_violation = c(0L, 2L),           # X1/Y1 passes both checks; X2/Y2 fails both
  stringsAsFactors = FALSE
)
ranked_pair <- ng_rank_cross_priority(pair, check_weight = 1, kinship_weight = 0)
stopifnot(nrow(ranked_pair) == 2L)        # never a veto -- both rows survive
stopifnot(identical(ranked_pair$parent1[[1L]], "X1"))    # 0 violations ranks ABOVE 2 violations
stopifnot(ranked_pair$priority_index[ranked_pair$parent1 == "X1"] >
          ranked_pair$priority_index[ranked_pair$parent1 == "X2"])
cat("check_violation priority component: check_weight=0 reproduces the no-column run; a positive\n",
    "check_weight demotes (never vetoes) an otherwise-identical cross that fails its checks\n")

cat("cross priority tests passed\n")
