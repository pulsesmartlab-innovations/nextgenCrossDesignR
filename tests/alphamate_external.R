helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

ids <- c("A", "B", "C", "D")
criterion <- c(A = 5, B = 4, C = 4, D = 3)
pairs <- ng_make_pairs(ids, include_self = FALSE)
p1 <- match(pairs$parent1, ids)
p2 <- match(pairs$parent2, ids)
scores <- data.frame(
  parent1 = pairs$parent1,
  parent2 = pairs$parent2,
  cross_mean = 0.5 * (criterion[p1] + criterion[p2]),
  pair_kinship = 0,
  stringsAsFactors = FALSE
)
scores$parent_distance <- seq_len(nrow(scores))
parent_kinship <- matrix(
  c(
    1.0, 0.9, 0.0, 0.0,
    0.9, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.2,
    0.0, 0.0, 0.2, 1.0
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(ids, ids)
)

recovered <- ng_parent_criterion_from_cross_mean(scores, "cross_mean")
stopifnot(all.equal(recovered[ids], criterion[ids], tolerance = 1e-8))

bad_plan <- data.frame(
  Parent1 = c("A", "Z"),
  Parent2 = c("A", "B"),
  MatingsCount = c(1L, 1L),
  stringsAsFactors = FALSE
)
bad_match <- tryCatch(
  ng_match_alphamate_plan(scores, bad_plan, n_crosses = 2L),
  error = function(e) e
)
stopifnot(inherits(bad_match, "error"))
stopifnot(grepl("self-cross", conditionMessage(bad_match), fixed = TRUE))
stopifnot(grepl("A x A", conditionMessage(bad_match), fixed = TRUE))
stopifnot(grepl("Z x B", conditionMessage(bad_match), fixed = TRUE))

exe <- normalizePath("external/AlphaMate/binaries/AlphaMate.exe", mustWork = FALSE)
if (!file.exists(exe)) {
  message("AlphaMate binary unavailable; parser/criterion tests passed")
  quit(save = "no", status = 0)
}

plan <- ng_select_alphamate(
  scores = scores,
  criterion_col = "cross_mean",
  n_crosses = 4L,
  parent_kinship = parent_kinship,
  executable = exe,
  runtime_path = "C:/Python/Lib/site-packages/torch/lib",
  target_degree = 45,
  max_contributions = 4L,
  number_of_parents = 4L,
  evol_iterations = 200L,
  evol_stop = 100L,
  workdir = tempfile("alphamate_test_"),
  keep_files = TRUE
)

stopifnot(nrow(plan) == 4L)
stopifnot(all(c("parent1", "parent2", "cross_mean") %in% names(plan)))
stopifnot(all(plan$parent1 %in% ids))
stopifnot(all(plan$parent2 %in% ids))
summary <- attr(plan, "summary")
stopifnot(is.list(summary))
stopifnot(identical(summary$lambda_parent_use_mode, "alphamate"))
stopifnot(isTRUE(grepl("AlphaMate.exe", summary$alphamate_executable, fixed = TRUE)))
stopifnot(identical(summary$alphamate_mode, "ModeOptTarget1"))
stopifnot(identical(summary$alphamate_exit_code, 0L))

long_ids <- paste0("ng_crop_aware_policy_ocs10_lps2_R1_C0_P", seq_along(ids))
names(long_ids) <- ids
long_scores <- scores
long_scores$parent1 <- unname(long_ids[long_scores$parent1])
long_scores$parent2 <- unname(long_ids[long_scores$parent2])
long_parent_K <- parent_kinship
rownames(long_parent_K) <- colnames(long_parent_K) <- unname(long_ids[ids])

long_plan <- ng_select_alphamate(
  scores = long_scores,
  criterion_col = "cross_mean",
  n_crosses = 4L,
  parent_kinship = long_parent_K,
  executable = exe,
  runtime_path = "C:/Python/Lib/site-packages/torch/lib",
  target_degree = 60,
  max_contributions = 4L,
  number_of_parents = 4L,
  evol_iterations = 200L,
  evol_stop = 100L,
  workdir = tempfile("alphamate_long_id_test_"),
  keep_files = FALSE
)

stopifnot(nrow(long_plan) == 4L)
stopifnot(all(long_plan$parent1 %in% unname(long_ids)))
stopifnot(all(long_plan$parent2 %in% unname(long_ids)))

cat("alphamate external tests passed\n")
