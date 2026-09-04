ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

ref <- list(active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                                stringsAsFactors = FALSE),
            values = list(yield = c(CHK_A = 6)), source = c(yield = "GEBV"))

# SINGLE TRAIT: the y axis IS the trait mean, so the line is exact
stopifnot(ng_check_line_value(ref, multi_trait_meta = NULL, trait = "yield") == 6)

# LINEAR INDEX: apply the same coefficients to the check's per-trait values
ref2 <- list(active = data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_A"),
                                 reject_if = "below", stringsAsFactors = FALSE),
             values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
             source = c(yield = "GEBV", protein = "GEBV"))
meta_lin <- list(method = "weighted", weights = c(yield = 0.5, protein = 0.5))
stopifnot(abs(ng_check_line_value(ref2, meta_lin) - 8) < 1e-8)

# RANK-BASED INDEX: a check has no rank, so there is NO valid line
meta_rank <- list(method = "rank_threshold")
stopifnot(is.na(ng_check_line_value(ref2, meta_rank)))

# an unevaluable check yields no line rather than a misplaced one
ref3 <- ref; ref3$values$yield <- c(CHK_A = NA_real_)
stopifnot(is.na(ng_check_line_value(ref3, NULL, trait = "yield")))

# the plot accepts a line and writes a file without error
scored <- data.frame(parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
                     multi_trait_score = c(9, 4), pair_kinship = c(0.1, 0.2),
                     priority_tier = c("priority", "low_priority"), stringsAsFactors = FALSE)
p <- file.path(tempdir(), "chk-plot.png")
ng_plot_priority_score_vs_kinship(scored, output_path = p, check_line = 6,
                                  check_label = "CHK_A")
stopifnot(file.exists(p), file.size(p) > 0)

# multi-trait: one panel per trait with a check, each with its own line
scored_mt <- data.frame(
  parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
  pair_kinship = c(0.1, 0.2),
  yield_mean = c(9, 4), yield_check_value = 6, yield_check_ok = c(TRUE, FALSE),
  protein_mean = c(11, 13), protein_check_value = 10, protein_check_ok = c(TRUE, TRUE),
  stringsAsFactors = FALSE)
ref_mt <- list(active = data.frame(trait = c("yield", "protein"), check = "CHK_A",
                                   reject_if = "below", stringsAsFactors = FALSE),
               values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
               source = c(yield = "GEBV", protein = "GEBV"))
p2 <- file.path(tempdir(), "chk-panels.png")
ng_plot_check_panels(scored_mt, ref_mt, output_path = p2)
stopifnot(file.exists(p2), file.size(p2) > 0)

# a run with no checks draws nothing and returns NULL rather than erroring
stopifnot(is.null(ng_plot_check_panels(scored_mt, NULL, output_path = p2)))

cat("check plot ok\n")
