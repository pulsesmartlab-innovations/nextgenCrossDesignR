# Real-data validation of the check-reference feature.
#
# Not a package test: test_data/ is ~1MB of real breeding data that is deliberately not a
# fixture, and tools/ is .Rbuildignore'd. Run it by hand:
#
#     Rscript tools/validate_check_reference_real_data.R
#
# What it checks, on 160 barley lines x 3191 markers with the breeder's own check varieties:
#   1. The four named checks are removed from the candidate parent set and never appear as a
#      parent of any candidate cross.
#   2. Check values resolve to REAL numbers on the run's own mean source, not NA.
#   3. Flags fall on the correct side for a DECREASE trait -- B_glucan and Protein are both
#      "lower is better", which is the mirrored sgn = -1 path that the synthetic unit tests
#      exercise least.
#   4. The identity invariant holds on real data: a run with checks is numerically identical to
#      the same run without them, differing only by added reference columns.
#   5. P(beat check) behaves sensibly against a deliberately non-uniform check spread.

ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

DATA <- "test_data"
CHECKS <- c("AAC_SYNERGY", "CONLON", "ND_GENESIS", "PINNACLE")
# Both are decrease traits: the breeder wants the mid-parent BELOW the check.
CHECK_FOR <- c(B_glucan = "AAC_SYNERGY",   # stringent: 167.5 against a population mean of 218
               Protein  = "ND_GENESIS")     # stringent: 10.73 against a population mean of 12.19
PROGENY_PER_FAMILY <- 200L                  # the breeder's own figure; there is no default
N_CROSSES <- 20L

say <- function(...) cat(sprintf(...), "\n", sep = "")
ok  <- function(cond, what) {
  if (isTRUE(cond)) say("  PASS  %s", what)
  else { say("  FAIL  %s", what); stop("validation failed: ", what, call. = FALSE) }
}

# ---- load ------------------------------------------------------------------
geno_raw <- utils::read.csv(file.path(DATA, "rich_geno.csv"), check.names = FALSE,
                            stringsAsFactors = FALSE)
pheno_raw <- utils::read.csv(file.path(DATA, "rich_pheno_new.csv"), check.names = FALSE,
                             stringsAsFactors = FALSE)
map_raw <- utils::read.csv(file.path(DATA, "geno_map.csv"), stringsAsFactors = FALSE)
dir_raw <- utils::read.csv(file.path(DATA, "selection_direction.csv"), stringsAsFactors = FALSE)

say("loaded: %d genotyped lines x %d markers; %d phenotyped; %d mapped markers",
    nrow(geno_raw), ncol(geno_raw) - 1L, nrow(pheno_raw), nrow(map_raw))

# Keep only markers present in BOTH the genotype table and the map. Two of the 3191 are
# unmapped; the alignment gate would surface that rather than intersect silently, so do it here
# deliberately and say so.
gm <- setdiff(names(geno_raw), "NAME")
shared <- intersect(gm, map_raw$SNP_code)
if (length(shared) < length(gm)) {
  say("note: dropping %d genotype marker(s) absent from the map: %s",
      length(gm) - length(shared), paste(setdiff(gm, shared), collapse = ", "))
}

geno_all <- as.matrix(geno_raw[, shared, drop = FALSE])
storage.mode(geno_all) <- "numeric"
rownames(geno_all) <- geno_raw$NAME

map <- data.frame(marker = shared,
                  chr = map_raw$Chromosome[match(shared, map_raw$SNP_code)],
                  bp = map_raw$Position_BP[match(shared, map_raw$SNP_code)],
                  stringsAsFactors = FALSE)

traits <- names(CHECK_FOR)
dir_df <- data.frame(
  trait = traits, column = traits,
  direction = dir_raw$Selection_direction[match(traits, dir_raw$Trait)],
  stringsAsFactors = FALSE)
ok(all(dir_df$direction == "decrease"),
   "both validation traits are DECREASE traits (the mirrored path)")

# ---- split the checks OUT of the candidate parents --------------------------
# This is the change the redesign forces. Under the old veto these four had to BE candidate
# parents; now they must not be, or the runner rejects them by design.
is_check <- rownames(geno_all) %in% CHECKS
ok(sum(is_check) == length(CHECKS), "all four named checks found in the genotype table")

# The candidate parents go in as a DATA FRAME with an explicit id column: passing a bare matrix
# loses the rownames, because ng_run_cp_canonical_id_table() coerces with as.data.frame().
# check_geno stays a matrix -- ng_align_check_geno() wants rownames and reads them directly.
parent_geno <- data.frame(NAME = rownames(geno_all)[!is_check],
                          geno_all[!is_check, , drop = FALSE],
                          check.names = FALSE, stringsAsFactors = FALSE)
check_geno  <- geno_all[is_check, , drop = FALSE]
parent_pheno <- pheno_raw[!(pheno_raw$NAME %in% CHECKS), , drop = FALSE]
check_pheno  <- pheno_raw[  pheno_raw$NAME %in% CHECKS,  , drop = FALSE]
rownames(parent_pheno) <- parent_pheno$NAME

say("split: %d candidate parents, %d check lines", nrow(parent_geno), nrow(check_geno))

# The run's mean source is resolved per trait; supply the checks' own values on the phenotypic
# sources so the comparison is apples-to-apples whichever one the run picks. Task 10 added
# `check_pheno` (a table shaped like the phenotype file), which is the friendlier route; this
# harness uses the lower-level check_records so it also exercises the escape hatch, and because
# check_records wins by documented precedence it pins the value regardless of source resolution.
check_records <- lapply(traits, function(tr) {
  v <- stats::setNames(suppressWarnings(as.numeric(check_pheno[[tr]])), check_pheno$NAME)
  list(adjusted_pheno = v, BLUE = v, BLUP = v)
})
names(check_records) <- traits

trait_checks <- data.frame(trait = traits, check = unname(CHECK_FOR[traits]),
                           stringsAsFactors = FALSE)

# parent_type = "ril", not the default "inbred". These lines carry ~0.8% heterozygous markers
# overall (worst line 3.3%), which trips the inbred gate: that gate is inbred_marker_fraction
# = 0.02, NOT inbred_tolerance = 0.05, so 9 of 154 parents exceed it. Residual heterozygosity at
# a few loci is exactly what RILs retain after finite selfing, and the QC message itself names
# "ril" as the remedy. This is a property of the germplasm, not a workaround: declaring the
# material honestly is what lets the a'Ra kernel model it correctly.
base_args <- list(
  genotype = parent_geno, phenotype = parent_pheno, marker_map = map,
  trait_direction = dir_df, id_col = "NAME", bp_per_cm = 1e6,
  parent_type = "ril", n_crosses = N_CROSSES, seed = 11L)

# ---- run, with and without ---------------------------------------------------
say("running WITHOUT checks ...")
without <- do.call(ng_run_cross_prediction, base_args)
say("running WITH checks ...")
with_ck <- do.call(ng_run_cross_prediction, c(base_args, list(
  check_geno = check_geno, check_records = check_records,
  check_progeny_size = PROGENY_PER_FAMILY, trait_checks = trait_checks)))

ct <- with_ck$candidate_crosses
ref <- with_ck$trait_check_reference

# ---- 1. checks are never crossed --------------------------------------------
ok(!any(c(ct$parent1, ct$parent2) %in% CHECKS),
   "no check line appears as a parent of any candidate cross")
# parent_geno is a data frame keyed by its NAME column, so check that column -- rownames() here
# would be 1..n and the assertion would pass vacuously.
ok(!any(CHECKS %in% parent_geno$NAME),
   "no check line is in the candidate parent table")
ok(nrow(parent_geno) == nrow(geno_all) - length(CHECKS),
   "the candidate parent table is exactly the genotype table minus the checks")

# ---- 2. values resolve to real numbers --------------------------------------
for (tr in traits) {
  v <- ref$values[[tr]][[CHECK_FOR[[tr]]]]
  ok(is.finite(v), sprintf("%s: check %s resolved to a real value (%.3f) on source '%s'",
                           tr, CHECK_FOR[[tr]], v, ref$source[[tr]]))
}
ok(identical(ref$progeny_size, PROGENY_PER_FAMILY),
   "the progeny-per-family figure used is the one supplied, and is stamped on the result")

# ---- 3. direction is correct for a DECREASE trait ----------------------------
for (tr in traits) {
  key <- ng_run_cp_clean_trait_name(tr)
  mu  <- ct[[paste0(key, "_mean")]]
  tau <- ct[[paste0(key, "_check_value")]]
  okc <- ct[[paste0(key, "_check_ok")]]
  vs  <- ct[[paste0(key, "_vs_check")]]
  # decrease: good means BELOW the check, so ok <=> mu <= tau, and the margin is check - mean
  ok(all((okc == (mu <= tau))[is.finite(mu) & is.finite(tau)]),
     sprintf("%s: check_ok is TRUE exactly where the mid-parent is BELOW the check", tr))
  ok(all(abs(vs - (tau - mu)) < 1e-8, na.rm = TRUE),
     sprintf("%s: vs_check is direction-aware (check - mean), so positive means better", tr))
  say("        %s: %d of %d crosses on the wrong side (mid-parent above %s = %.2f)",
      tr, sum(okc %in% FALSE), length(okc), CHECK_FOR[[tr]], tau[[1L]])
}

# ---- 4. the identity invariant, on real data ---------------------------------
a <- without$candidate_crosses; b <- ct
ok(nrow(a) == nrow(b), "candidate cross count unchanged by adding checks")
n_same <- 0L
for (nm in names(a)) {
  stopifnot(nm %in% names(b))
  if (!isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 0))) {
    stop("INVARIANT BROKEN on real data, column: ", nm, call. = FALSE)
  }
  n_same <- n_same + 1L
}
ok(TRUE, sprintf("all %d pre-existing candidate columns identical at tolerance = 0", n_same))
for (nm in names(without$selected_crosses)) {
  if (!isTRUE(all.equal(without$selected_crosses[[nm]], with_ck$selected_crosses[[nm]],
                        tolerance = 0))) {
    stop("INVARIANT BROKEN on real data, selected_crosses column: ", nm, call. = FALSE)
  }
}
ok(isTRUE(all.equal(without$plan_summary, with_ck$plan_summary, tolerance = 0)),
   "the allocation plan is unchanged by adding checks")

# ---- 5. P(beat check) is sensible --------------------------------------------
for (tr in traits) {
  key <- ng_run_cp_clean_trait_name(tr)
  p <- ct[[paste0(key, "_p_beat_check")]]
  okc <- ct[[paste0(key, "_check_ok")]]
  ok(all(p >= 0 & p <= 1, na.rm = TRUE), sprintf("%s: P(beat check) lies in [0, 1]", tr))
  if (any(okc %in% FALSE) && any(is.finite(p[okc %in% FALSE]))) {
    say("        %s: of the %d below-par crosses, best tail probability = %.3f",
        tr, sum(okc %in% FALSE), max(p[okc %in% FALSE], na.rm = TRUE))
  }
}
if ("p_beat_all_checks" %in% names(ct)) {
  pj <- ct$p_beat_all_checks
  ok(all(pj >= 0 & pj <= 1, na.rm = TRUE), "P(beat all checks) lies in [0, 1]")
  for (tr in traits) {
    key <- ng_run_cp_clean_trait_name(tr)
    ok(all(pj <= ct[[paste0(key, "_p_beat_check")]] + 1e-8, na.rm = TRUE),
       sprintf("P(beat all checks) never exceeds the %s marginal", tr))
  }
}

say("")
say("ALL REAL-DATA CHECKS PASSED")
