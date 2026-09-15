# The workbook must not print a mid-parent on one basis next to a check
# comparison computed on another.
#
# Regression test for the 2026-09-11 finding. `ng_cpw_gebv_cols()` adds
# `<trait>_mean_gebv` (relabelled `<trait>_mid_parent_gebv`) when
# include_trait_gebv = TRUE, but NOTHING ever added the plain `<trait>_mean` --
# the BLUP/selection-basis mid-parent that `<trait>_vs_check` is computed from.
#
# The delivered 2026-09 barley workbooks therefore showed, for every one of 17
# traits, a mid-parent that could not reproduce the vs_check sitting beside it.
# For YIELD: mid_parent_gebv 97.013 against check 99.640 -- BELOW it -- while
# reporting vs_check = +3.840 and check_ok = TRUE. A breeder reading that sheet
# reasonably concludes the mid-parents are wrong.
#
# The invariant pinned here is the one that failed 17/17 and that no amount of
# reading column NAMES would have caught:
#
#     <trait>_mid_parent - <trait>_check_value == <trait>_vs_check
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

n <- 8L
mid_parent <- seq(103.5, 96.5, length.out = n)   # BLUP basis: what selection used
check_value <- 99.64
crosses <- data.frame(
  parent1 = sprintf("P%02d", seq_len(n)),
  parent2 = sprintf("Q%02d", seq_len(n)),
  multi_trait_score = seq(0.90, 0.55, length.out = n),
  pair_kinship = seq(-0.30, 0.10, length.out = n),
  YIELD_mean = mid_parent,
  # GEBV basis: deliberately on a DIFFERENT scale and nearly flat, exactly as in
  # the barley run, so a test that silently swapped the two bases would fail.
  YIELD_mean_gebv = seq(97.0, 97.3, length.out = n),
  YIELD_check_id = "ND GENESIS",
  YIELD_check_value = check_value,
  YIELD_vs_check = mid_parent - check_value,
  YIELD_check_ok = (mid_parent - check_value) > 0,
  stringsAsFactors = FALSE
)
crosses <- ng_rank_cross_priority(crosses, breaks = c(0.25, 0.50, 0.75, 1),
                                  kinship_weight = 0.10)

trait_directions <- data.frame(trait = "YIELD", column = "YIELD_mean",
                               direction = "maximize", weight = 1,
                               stringsAsFactors = FALSE)

tables <- ng_cross_priority_workbook_tables(
  crosses = crosses,
  scored = crosses,
  trait_directions = trait_directions,
  n_crosses_requested = n,
  block_size = 4L,
  include_trait_gebv = TRUE
)

sel <- tables[["Selected_All"]]
stopifnot(is.data.frame(sel), nrow(sel) == n)

# 1. The GEBV mid-parent is opt-in and was requested, so it must be present.
stopifnot("YIELD_mid_parent_gebv" %in% names(sel))

# 2. THE ASSERTION: the basis the check comparison used must be on the sheet too.
stopifnot("YIELD_mid_parent" %in% names(sel))

# 3. And it must actually reconcile -- the identity that failed 17/17.
lhs <- as.numeric(sel[["YIELD_mid_parent"]]) - as.numeric(sel[["YIELD_check_value"]])
stopifnot(isTRUE(all.equal(lhs, as.numeric(sel[["YIELD_vs_check"]]), tolerance = 1e-8)))

# 4. The two bases must remain distinct columns: printing the GEBV value under
#    the plain name would satisfy (2) and (3) is what would catch it.
stopifnot(!isTRUE(all.equal(as.numeric(sel[["YIELD_mid_parent"]]),
                            as.numeric(sel[["YIELD_mid_parent_gebv"]]))))

# 5. The same invariant for a MINIMIZE trait. `<trait>_vs_check` is signed so
#    positive always means favourable, so the subtraction flips with direction:
#      increase: mid_parent - check      decrease: check - mid_parent
#    Verified against all 17 delivered barley traits (8/8 and 9/9 respectively).
#    This assertion passes on the writer as it stands -- the writer copies
#    `<trait>_mean` and does not recompute vs_check, so there is no failing state
#    to observe first. It is kept deliberately, as a characterization guard: the
#    convention is undocumented in the codebase, and a future change that
#    recomputed vs_check with one fixed sign would silently corrupt the 9
#    minimize traits while leaving the 8 maximize traits correct. That asymmetry
#    is very hard to notice by eye, and it is the exact mistake made while
#    repairing the delivered files.
crosses$DISEASE_mean <- seq(2.0, 6.0, length.out = n)
crosses$DISEASE_mean_gebv <- seq(3.9, 4.1, length.out = n)
crosses$DISEASE_check_value <- 4.5
crosses$DISEASE_vs_check <- 4.5 - crosses$DISEASE_mean      # decrease: check - mid_parent

tables2 <- ng_cross_priority_workbook_tables(
  crosses = crosses,
  scored = crosses,
  trait_directions = rbind(trait_directions,
                           data.frame(trait = "DISEASE", column = "DISEASE_mean",
                                      direction = "minimize", weight = 1,
                                      stringsAsFactors = FALSE)),
  n_crosses_requested = n,
  block_size = 4L,
  include_trait_gebv = TRUE
)
sel2 <- tables2[["Selected_All"]]
stopifnot("DISEASE_mid_parent" %in% names(sel2))
stopifnot(isTRUE(all.equal(
  as.numeric(sel2[["DISEASE_check_value"]]) - as.numeric(sel2[["DISEASE_mid_parent"]]),
  as.numeric(sel2[["DISEASE_vs_check"]]), tolerance = 1e-8)))

message("workbook_mid_parent_basis passed")
