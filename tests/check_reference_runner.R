ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11)
n_p <- 12L; n_m <- 40L
ids <- paste0("P", seq_len(n_p))
markers <- paste0("m", seq_len(n_m))
gm <- matrix(2L * rbinom(n_p * n_m, 1, 0.4), nrow = n_p, dimnames = list(ids, markers))
genotype <- data.frame(id = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(marker = markers, chr = rep(1:4, length.out = n_m),
                         bp = rep(seq(0, 100, length.out = 10), 4)[seq_len(n_m)] * 1e6,
                         stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)
pheno <- data.frame(id = ids, yield = rnorm(n_p, 10, 2), stringsAsFactors = FALSE)

# the check is NOT among the parents -- the case the old veto rejected outright
chk <- matrix(2L * rbinom(2 * n_m, 1, 0.4), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), markers))

args <- list(genotype = genotype, phenotype = pheno, marker_map = marker_map,
             map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
             bp_per_cm = 1e6, id_col = "id",
             trait_direction = direction, n_crosses = 5L,
             write_outputs = FALSE, write_figures = FALSE, seed = 5L)

# The runner's default ridge fit never stamps a calibrated reliability (R/02_effects.R sets
# reliability_is_calibrated = FALSE), so ng_choose_mean_source() always resolves to the
# phenotypic "adjusted_pheno" source here (it is passed the raw phenotype unconditionally) --
# never a GEBV variant. check_records is exactly the mechanism the interface names for this
# case ("only consulted when the run's source for that trait is not a GEBV variant").
chk_records <- list(yield = list(adjusted_pheno = c(CHK_A = 11.5, CHK_B = 9.25)))

res <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_records = chk_records,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))

ref <- res$trait_check_reference
stopifnot(is.list(ref), nrow(ref$active) == 1L, ref$active$check == "CHK_A")
stopifnot(is.finite(ref$values$yield[["CHK_A"]]))
stopifnot(abs(ref$values$yield[["CHK_A"]] - 11.5) < 1e-8)
stopifnot(identical(ref$source[["yield"]], "adjusted_pheno"))
stopifnot(nzchar(ref$source[["yield"]]))

ct <- res$candidate_crosses
stopifnot(all(c("yield_check_value", "yield_vs_check", "yield_check_ok",
                "yield_p_beat_check", "checks_all_ok") %in% names(ct)))
# the check is never a parent of any candidate cross
stopifnot(!any(c(ct$parent1, ct$parent2) %in% c("CHK_A", "CHK_B")))

# a check id colliding with a parent id is rejected
err <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = matrix(0, nrow = 1, ncol = n_m,
                      dimnames = list("P1", markers)),
  check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "P1", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("also a candidate parent", err))

# a check named in trait_checks but absent from check_geno is rejected
err2 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "NOPE", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err2), grepl("NOPE", err2))

# check_progeny_size is REQUIRED with trait_checks -- never defaulted
err3 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err3), grepl("check_progeny_size", err3))
# and it must be a sensible count
err4 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 0L,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err4), grepl("check_progeny_size", err4))

# the size the user gave is what drives the probability, and is stamped on the result
stopifnot(identical(ref$progeny_size, 200L))

cat("runner wiring ok\n")

# --- Task 10: check_pheno -> real check values on a phenotype-source run ----
# check_pheno makes the check values REAL rather than NA on a phenotype-source run.
# check_pheno's id column must be the run's own id_col ("id" in `args`, above) --
# it is read the same way the phenotype table itself is, not by a separate id keyword.
cp <- data.frame(id = c("CHK_A", "CHK_B"), yield = c(12.5, 8.0), stringsAsFactors = FALSE)
res_cp <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(identical(res_cp$trait_check_reference$source[["yield"]], "adjusted_pheno"))
stopifnot(abs(res_cp$trait_check_reference$values$yield[["CHK_A"]] - 12.5) < 1e-8)
ct_cp <- res_cp$candidate_crosses
stopifnot(all(is.finite(ct_cp$yield_check_value)))
stopifnot(any(is.finite(ct_cp$yield_p_beat_check)))

# precedence: an explicit check_records value WINS over check_pheno for the same trait, even
# when both are supplied for it
res_prec <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp,
  check_records = list(yield = list(adjusted_pheno = c(CHK_A = 99))),
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(abs(res_prec$trait_check_reference$values$yield[["CHK_A"]] - 99) < 1e-8)

cat("task 10 runner ok\n")

# --- Task 10 fix round 1: check_pheno's id column resolves like every other table's ----------
# id_col OMITTED entirely: check_pheno keyed by NAME (the column real project phenotype files
# use, e.g. test_data/rich_pheno_new.csv) must still auto-detect and produce real check values,
# the same way phenotype/genotype auto-detect their own id column when id_col is not given.
args_auto <- args
args_auto$id_col <- NULL
cp_name <- data.frame(NAME = c("CHK_A", "CHK_B"), yield = c(15.0, 7.0), stringsAsFactors = FALSE)
res_auto <- do.call(ng_run_cross_prediction, c(args_auto, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp_name,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(identical(res_auto$trait_check_reference$source[["yield"]], "adjusted_pheno"))
stopifnot(abs(res_auto$trait_check_reference$values$yield[["CHK_A"]] - 15.0) < 1e-8)
ct_auto <- res_auto$candidate_crosses
stopifnot(all(is.finite(ct_auto$yield_check_value)))

# id_col GIVEN explicitly but absent from check_pheno: must still fail with the existing,
# clear message -- auto-detection must never mask a genuinely wrong explicit column name.
cp_wrong_col <- data.frame(NAME = c("CHK_A", "CHK_B"), yield = c(1, 2), stringsAsFactors = FALSE)
err_idcol <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp_wrong_col,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err_idcol), grepl("check_pheno", err_idcol),
          grepl("missing its id column", err_idcol), grepl("id", err_idcol))

cat("task 10 fix round 1 ok\n")

# --- Task 10 fix round 2: check_pheno resolution must be LAZY -------------------------------
# id_col omitted; the only trait ("yield") is supplied entirely via explicit check_records
# (which wins over check_pheno per the documented precedence), and check_pheno is present but
# keyed by a column OUTSIDE the auto-detect candidate list (accession_code, not id/name/...).
# Because check_pheno is never actually consulted for this trait, attaching it must NOT turn a
# working run into a failing one -- eager id-column resolution (fix round 1's first attempt)
# would fail here even though nothing ever reads check_pheno.
cp_unresolvable <- data.frame(accession_code = c("CHK_A", "CHK_B"), yield = c(1, 2),
                              stringsAsFactors = FALSE)
res_lazy <- do.call(ng_run_cross_prediction, c(args_auto, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp_unresolvable,
  check_records = list(yield = list(adjusted_pheno = c(CHK_A = 77))),
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(abs(res_lazy$trait_check_reference$values$yield[["CHK_A"]] - 77) < 1e-8)

cat("task 10 fix round 2 ok\n")

# --- I3: a GEBV-sourced trait must never trigger check_pheno id-column resolution -----------
# The stage_index loop gated check_pheno consultation on `is.null(recs) && !is.null(check_pheno)`
# alone, which does not test the SOURCE: for a GEBV-sourced trait, ng_check_records_from_pheno()
# itself returns NULL immediately (no check_records supplied --> recs stays NULL too), so the
# branch still ran and called ng_run_cp_id_col() on check_pheno -- hard-erroring on an
# unrecognised id column even though nothing in this trait's evaluation needs check_pheno at all.
# This codebase's own ridge fit never calibrates reliability (R/02_effects.R sets
# reliability_is_calibrated = FALSE unconditionally), so a GEBV* mean_source is otherwise
# unreachable end-to-end via the normal per-trait fitting path -- force it by monkey-patching the
# index stage (restored immediately after, success or failure) to override trait_mean_source
# right before the check-reference block consumes it. check_pheno is keyed by a column outside
# ng_run_cp_id_col()'s auto-detect candidate list, with id_col omitted so auto-detection (not an
# explicit column) is what would fire if check_pheno were consulted.
orig_index_stage <- ng_cp_pipeline$index
ng_cp_pipeline$index <- function(ctx) {
  ctx$trait_mean_source$yield <- "GEBV_uncalibrated"
  orig_index_stage(ctx)
}
cp_gebv_unresolvable <- data.frame(accession_code = c("CHK_A", "CHK_B"), yield = c(1, 2),
                                   stringsAsFactors = FALSE)
res_gebv <- tryCatch(
  do.call(ng_run_cross_prediction, c(args_auto, list(
    check_geno = chk, check_progeny_size = 200L, check_pheno = cp_gebv_unresolvable,
    trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) e)
ng_cp_pipeline$index <- orig_index_stage   # restore before any stopifnot can abort the script

stopifnot(!inherits(res_gebv, "error"))
stopifnot(identical(res_gebv$trait_check_reference$source[["yield"]], "GEBV_uncalibrated"))
# the check's value was predicted from its OWN markers (check_geno), never from check_pheno --
# proof that check_pheno was correctly never consulted for this trait.
stopifnot(is.finite(res_gebv$trait_check_reference$values$yield[["CHK_A"]]))

cat("I3 gebv check_pheno gate ok\n")

# --- I2: priority_check_weight was plumbed onto ng_run_cross_prediction() but never threaded to
# the ng_rank_cross_priority() call inside ng_cp__stage_rank() -- so a run could never let
# check_violation influence which tier a cross lands in, even though ng_rank_cross_priority()
# itself supports check_weight (tests/cross_priority.R). Prove the fix moves an ACTUAL cross's
# tier THROUGH THE RUNNER, not only via the exported function.
#
# Phase 1: a plain run (no synthetic override, `args` has no trait_checks/check_geno so
# cross_table carries no check_violation column at all) identifies which cross is the pool's
# best-ranked by score/kinship/threshold alone.
res_i2_base <- do.call(ng_run_cross_prediction, c(args, list(priority_check_weight = 0)))
sel_base <- res_i2_base$selected_crosses
stopifnot(!("check_violation" %in% names(sel_base)))
top_idx <- which.min(sel_base$priority_rank)
top_p1 <- sel_base$parent1[[top_idx]]; top_p2 <- sel_base$parent2[[top_idx]]
stopifnot(identical(as.character(sel_base$priority_tier[[top_idx]]), "highly_priority"))

# Phase 2: monkey-patch the index stage to attach a synthetic check_violation column that is 0
# everywhere except an EXTREME value on that exact SAME cross -- deterministic regardless of the
# biological objective's internal tie-breaking, and it does not affect scored_crosses$multi_trait_
# score, so the allocator selects the identical set of crosses either way (checks never
# veto/reselect).
orig_index_stage3 <- ng_cp_pipeline$index
ng_cp_pipeline$index <- function(ctx) {
  ctx <- orig_index_stage3(ctx)
  ct <- ctx$cross_table
  is_top <- ct$parent1 == top_p1 & ct$parent2 == top_p2
  stopifnot(sum(is_top) == 1L)
  ct$check_violation <- 0L
  ct$check_violation[is_top] <- 1000L
  ctx$cross_table <- ct
  ctx
}
res_i2_a <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(priority_check_weight = 0))),
                     error = function(e) e)
res_i2_b <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(priority_check_weight = 1e6))),
                     error = function(e) e)
ng_cp_pipeline$index <- orig_index_stage3   # restore before any stopifnot can abort the script

stopifnot(!inherits(res_i2_a, "error"), !inherits(res_i2_b, "error"))
sel_a <- res_i2_a$selected_crosses
sel_b <- res_i2_b$selected_crosses
stopifnot(nrow(sel_a) == nrow(sel_b))
# same crosses selected regardless of priority_check_weight -- checks never veto or reselect
key_a <- paste(sel_a$parent1, sel_a$parent2)
key_b <- paste(sel_b$parent1, sel_b$parent2)
stopifnot(setequal(key_a, key_b))

top_key <- paste(top_p1, top_p2)
tier_a <- as.character(sel_a$priority_tier[key_a == top_key])
tier_b <- as.character(sel_b$priority_tier[key_b == top_key])
# priority_check_weight = 0 (even with the synthetic column present) reproduces the no-check
# ranking exactly -- the top cross is still "highly_priority". At a check_weight overwhelming
# every other component, that SAME cross's (synthetic) maximal violation count swamps its top
# score and drops it to the worst tier -- proof the weight reaches ranking through the runner.
stopifnot(identical(tier_a, "highly_priority"))
stopifnot(!identical(tier_a, tier_b))
stopifnot(identical(tier_b, "low_priority"))
cat(sprintf("I2 ok: priority_check_weight reaches the runner -- top cross moved %s -> %s\n",
            tier_a, tier_b))
