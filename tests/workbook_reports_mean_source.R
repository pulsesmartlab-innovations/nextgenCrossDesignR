# The workbook must state which basis the cross mean came from.
#
# Regression test for the 2026-09-11 finding. `ng_choose_mean_source()` picks
# between GEBV and a phenotype fallback (BLUP / BLUE / adjusted_pheno), and that
# choice decides what `<trait>_mean` -- and therefore the ranking and the check
# comparison -- actually IS. It was recorded nowhere a reader would look:
#
#   * `ng_cpw_scoring_method()` took no arguments at all and emitted a fixed
#     description of the engine, identical for every run;
#   * `trait_mean_source` was assembled in the runner but did not reach
#     result.json;
#   * the only surviving record was inside `posterior_predictions`, which is
#     absent entirely when run_posterior_prediction = FALSE.
#
# So a 17-trait study ran to completion, was delivered, and was read by a
# breeder before anyone could tell that every trait had silently fallen back to
# the phenotype mean. The number was right; nothing said what it was.
#
# A spreadsheet outlives the session that explains it. If the basis is not on
# the face of the workbook, it is not reported.
# openxlsx is attached BEFORE .Rlib is prepended: the repo carries a
# Windows-built .Rlib whose Rcpp has no macOS shared object, and openxlsx pulls
# Rcpp in. Resolving it from the normal user library first avoids that shadow.
suppressPackageStartupMessages(library(openxlsx))
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

n <- 6L
crosses <- data.frame(
  parent1 = sprintf("P%02d", seq_len(n)),
  parent2 = sprintf("Q%02d", seq_len(n)),
  multi_trait_score = seq(0.9, 0.5, length.out = n),
  pair_kinship = seq(-0.2, 0.1, length.out = n),
  YIELD_mean = seq(103, 97, length.out = n),
  YIELD_mean_gebv = seq(97.0, 97.3, length.out = n),
  stringsAsFactors = FALSE
)
crosses <- ng_rank_cross_priority(crosses, breaks = c(0.25, 0.5, 0.75, 1),
                                  kinship_weight = 0.1)
trait_directions <- data.frame(trait = "YIELD", column = "YIELD_mean",
                               direction = "maximize", weight = 1,
                               stringsAsFactors = FALSE)

tables <- ng_cross_priority_workbook_tables(
  crosses = crosses,
  scored = crosses,
  trait_directions = trait_directions,
  n_crosses_requested = n,
  block_size = 3L,
  include_trait_gebv = TRUE,
  trait_mean_source = list(YIELD = "adjusted_pheno")
)

sm <- tables[["Scoring_Method"]]
stopifnot(is.data.frame(sm), nrow(sm) > 0)
blob <- paste(unlist(sm), collapse = " | ")

# THE ASSERTION: the sheet must name the basis actually used in THIS run.
stopifnot(grepl("adjusted_pheno", blob, fixed = TRUE))
# And it must name the trait it applies to, since a multi-trait run can mix bases.
stopifnot(grepl("YIELD", blob, fixed = TRUE))

# A run that used GEBV must say so instead -- i.e. the row is driven by the
# argument and is not boilerplate that happens to contain the word.
tables_gebv <- ng_cross_priority_workbook_tables(
  crosses = crosses, scored = crosses, trait_directions = trait_directions,
  n_crosses_requested = n, block_size = 3L, include_trait_gebv = TRUE,
  trait_mean_source = list(YIELD = "GEBV"))
blob_gebv <- paste(unlist(tables_gebv[["Scoring_Method"]]), collapse = " | ")
stopifnot(grepl("GEBV", blob_gebv, fixed = TRUE))
stopifnot(!grepl("adjusted_pheno", blob_gebv, fixed = TRUE))

# And through ng_write_cross_priority_workbook() -- the path a real run takes.
#
# This is not redundant with the assertions above. The first version of this fix
# added `trait_mean_source` as a formal of the writer but never FORWARDED it to
# ng_cross_priority_workbook_tables(), so every assertion above passed while a
# real run still produced a workbook with no basis row. Testing the tables
# function alone cannot see a broken hand-off between the two.
out <- tempfile(fileext = ".xlsx")
ng_write_cross_priority_workbook(
  output_path = out,
  crosses = crosses,
  scored = crosses,
  trait_directions = trait_directions,
  n_crosses_requested = n,
  block_size = 3L,
  include_trait_gebv = TRUE,
  trait_mean_source = list(YIELD = "adjusted_pheno")
)
stopifnot(file.exists(out))
sm_written <- openxlsx::read.xlsx(out, sheet = "Scoring_Method", startRow = 2)
blob_written <- paste(unlist(sm_written), collapse = " | ")
stopifnot(grepl("adjusted_pheno", blob_written, fixed = TRUE))
unlink(out)

message("workbook_reports_mean_source passed")
