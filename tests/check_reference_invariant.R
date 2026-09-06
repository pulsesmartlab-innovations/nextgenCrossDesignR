ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(202)
n_p <- 14L; n_m <- 50L
ids <- paste0("P", seq_len(n_p))
markers <- paste0("m", seq_len(n_m))
# Homozygous-only dosages (0/2). The brief's draft used rbinom(n, 2, p), which also yields 1s
# (heterozygous calls) and trips the inbred-line QC gate under the default
# parent_type = "inbred" (see check_reference_runner.R for the same fix in Task 6).
gm <- matrix(2L * rbinom(n_p * n_m, 1, 0.35), nrow = n_p, dimnames = list(ids, markers))
genotype <- data.frame(id = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(marker = markers, chr = rep(1:4, length.out = n_m),
                         bp = rep(seq(0, 100, length.out = ceiling(n_m / 4)), 4)[seq_len(n_m)] * 1e6,
                         stringsAsFactors = FALSE)
pheno <- data.frame(id = ids, yield = rnorm(n_p, 10, 2), protein = rnorm(n_p, 12, 1),
                    stringsAsFactors = FALSE)
direction <- data.frame(trait = c("yield", "protein"), column = c("yield", "protein"),
                        direction = c("increase", "increase"), stringsAsFactors = FALSE)

# The checks are NOT among the parents -- the case the pre-Task-6 veto rejected outright.
chk <- matrix(2L * rbinom(2 * n_m, 1, 0.35), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), markers))

# The runner's default ridge fit never stamps a calibrated reliability (R/02_effects.R sets
# reliability_is_calibrated = FALSE), so ng_choose_mean_source() always resolves to the
# phenotypic "adjusted_pheno" source here, never a GEBV variant (verified in Task 6). Supplying
# adjusted_pheno check_records exercises a real, non-NA code path instead of leaving the check
# columns at NA.
chk_records <- list(
  yield = list(adjusted_pheno = c(CHK_A = 11.0, CHK_B = 9.5)),
  protein = list(adjusted_pheno = c(CHK_A = 12.5, CHK_B = 11.8))
)

trait_checks <- data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_B"),
                           stringsAsFactors = FALSE)

args <- list(genotype = genotype, phenotype = pheno, marker_map = marker_map,
             map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
             bp_per_cm = 1e6, id_col = "id",
             trait_direction = direction, n_crosses = 6L,
             write_outputs = FALSE, write_figures = FALSE, seed = 5L)

without <- do.call(ng_run_cross_prediction, args)
with_ck <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_records = chk_records,
  trait_checks = trait_checks)))

# I4: the joint p_beat_all_checks caveat (Monte Carlo draw count/error, diagonal PEV
# approximation, clipped to the marginal minimum) must actually reach the run's returned
# trait_check_reference, not live only in an attribute nothing downstream reads.
note <- with_ck$trait_check_reference$p_beat_all_checks_note
stopifnot(is.character(note), length(note) == 1L, nzchar(note))
stopifnot(is.null(without$trait_check_reference))

# every column the run produced WITHOUT checks must be untouched BY checks
a <- without$candidate_crosses
b <- with_ck$candidate_crosses
stopifnot(nrow(a) == nrow(b))
n_compared <- 0L
for (nm in names(a)) {
  stopifnot(nm %in% names(b))
  stopifnot(isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 0)))
  n_compared <- n_compared + 1L
}

# the allocation is identical too -- checks must not shift which crosses are selected or how
# much of the budget each gets. The runner returns no "crossing_plan" field (that name does not
# exist on ng_run_cross_prediction()'s output); the real allocation outputs are
# selected_crosses (a data.frame) and plan_summary (a diagnostics list).
sa <- without$selected_crosses
sb <- with_ck$selected_crosses
stopifnot(nrow(sa) == nrow(sb))
n_selected_compared <- 0L
for (nm in names(sa)) {
  # selected_crosses also gains the 12 check columns (it is derived from the same
  # check-annotated scored table as candidate_crosses), so this must be column-wise, not a
  # whole-object all.equal() -- a whole-object comparison would fail spuriously on the added
  # columns even when the allocation itself is untouched.
  stopifnot(nm %in% names(sb))
  stopifnot(isTRUE(all.equal(sa[[nm]], sb[[nm]], tolerance = 0)))
  n_selected_compared <- n_selected_compared + 1L
}

# plan_summary carries no per-cross check columns, so it CAN be compared whole-object.
stopifnot(isTRUE(all.equal(without$plan_summary, with_ck$plan_summary, tolerance = 0)))

# and the only difference is ADDED columns -- exactly the check-reference columns the
# interface promises: per trait, <key>_check_id / _check_value / _vs_check / _check_ok /
# _p_beat_check; globally, checks_all_ok, p_beat_all_checks (because two traits are checked
# here), and check_violation (the INTEGRATION dimensionless wrong-side count -- exists ONLY on
# check runs; priority_check_component / priority_rule stay identical between the two runs
# because ng_rank_cross_priority()'s check component is gated on check_weight > 0, which neither
# run here opts into).
added <- setdiff(names(b), names(a))
stopifnot(length(added) > 0L)
allowed_suffix <- c("_check_id", "_check_value", "_vs_check", "_check_ok", "_p_beat_check")
allowed_exact <- c("checks_all_ok", "p_beat_all_checks", "check_violation")
is_allowed <- vapply(added, function(nm) {
  nm %in% allowed_exact ||
    any(vapply(allowed_suffix, function(s) grepl(paste0(s, "$"), nm), logical(1)))
}, logical(1))
stopifnot(all(is_allowed))

cat(sprintf(
  "invariant ok: %d candidate_crosses columns + %d selected_crosses columns + plan_summary identical, checks add columns and change nothing else\n",
  n_compared, n_selected_compared))

# --- I1: the SAME invariant at >= 3 checked traits, where mvtnorm::pmvnorm() (used by the joint
# p_beat_all_checks path, R/33/R/51) switches from a deterministic closed form (2 traits) to the
# randomised GenzBretz lattice rule -- which both draws from and advances the ambient
# .Random.seed. Left unguarded, a with-checks run's check-probability computation shifts the
# global RNG stream before optimizer = "evolution" (the capability registry's recommended
# option, R/25) runs with its default evol_seed = NULL (R/40_evolutionary_optimizer.R), i.e. off
# whatever the ambient stream happens to be -- so a with-checks run and a without-checks run can
# select DIFFERENT crossing plans even though checks are supposed to be a reference only. The
# 2-trait case above cannot see this bug (pmvnorm is deterministic at 2 dimensions), which is
# exactly why this second, 3-trait scenario exists. A large-enough, genuinely combinatorial
# allocation problem (20 parents, 15 crosses, max_crosses_per_parent constraint) is needed for
# the evolutionary optimizer's outcome to be sensitive to the ambient RNG stream at all -- a
# tiny problem's memetic warm-start + elitism converges to the same optimum regardless.
set.seed(303)
n_p3 <- 20L; n_m3 <- 50L
ids3 <- paste0("Q", seq_len(n_p3))
markers3 <- paste0("m", seq_len(n_m3))
gm3 <- matrix(2L * rbinom(n_p3 * n_m3, 1, 0.35), nrow = n_p3, dimnames = list(ids3, markers3))
genotype3 <- data.frame(id = ids3, gm3, check.names = FALSE, stringsAsFactors = FALSE)
marker_map3 <- data.frame(marker = markers3, chr = rep(1:4, length.out = n_m3),
                          bp = rep(seq(0, 100, length.out = ceiling(n_m3 / 4)), 4)[seq_len(n_m3)] * 1e6,
                          stringsAsFactors = FALSE)
pheno3 <- data.frame(id = ids3, yield = rnorm(n_p3, 10, 2), protein = rnorm(n_p3, 12, 1),
                     oil = rnorm(n_p3, 5, 0.5), stringsAsFactors = FALSE)
direction3 <- data.frame(trait = c("yield", "protein", "oil"),
                         column = c("yield", "protein", "oil"),
                         direction = c("increase", "increase", "increase"),
                         stringsAsFactors = FALSE)
chk3 <- matrix(2L * rbinom(2 * n_m3, 1, 0.35), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), markers3))
chk_records3 <- list(
  yield = list(adjusted_pheno = c(CHK_A = 11.0, CHK_B = 9.5)),
  protein = list(adjusted_pheno = c(CHK_A = 12.5, CHK_B = 11.8)),
  oil = list(adjusted_pheno = c(CHK_A = 5.2, CHK_B = 4.8))
)
trait_checks3 <- data.frame(trait = c("yield", "protein", "oil"),
                            check = c("CHK_A", "CHK_B", "CHK_A"),
                            stringsAsFactors = FALSE)

args3 <- list(genotype = genotype3, phenotype = pheno3, marker_map = marker_map3,
             map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
             bp_per_cm = 1e6, id_col = "id",
             trait_direction = direction3, n_crosses = 15L,
             max_crosses_per_parent = 3L,
             write_outputs = FALSE, write_figures = FALSE, seed = 5L,
             optimizer = "evolution", evol_solutions = 30L, evol_iterations = 40L,
             evol_stop = 10L)   # evol_seed left at its default NULL -- the failure scenario

without3 <- do.call(ng_run_cross_prediction, args3)
with_ck3 <- do.call(ng_run_cross_prediction, c(args3, list(
  check_geno = chk3, check_progeny_size = 200L, check_records = chk_records3,
  trait_checks = trait_checks3)))

sa3 <- without3$selected_crosses
sb3 <- with_ck3$selected_crosses
stopifnot(nrow(sa3) == nrow(sb3))
key_a3 <- paste(sa3$parent1, sa3$parent2)
key_b3 <- paste(sb3$parent1, sb3$parent2)
# the SAME set of crosses, in the SAME order -- checks must not perturb which plan the
# evolutionary optimizer converges to, even at >= 3 checked traits.
stopifnot(identical(key_a3, key_b3))
common3 <- intersect(names(sa3), names(sb3))
n3_compared <- 0L
for (nm in common3) {
  stopifnot(isTRUE(all.equal(sa3[[nm]], sb3[[nm]], tolerance = 0)))
  n3_compared <- n3_compared + 1L
}
stopifnot(isTRUE(all.equal(without3$plan_summary, with_ck3$plan_summary, tolerance = 0)))

cat(sprintf(
  "I1 3-trait invariant ok: %d shared selected_crosses columns identical under optimizer = 'evolution' with evol_seed = NULL, at 3 checked traits (pmvnorm dimension >= 3)\n",
  n3_compared))
