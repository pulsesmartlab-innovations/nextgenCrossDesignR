ng_test_use_cpp <- TRUE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- synthetic inputs (self-contained; no external golden fixture) -------
# 40 lines x 240 markers over 6 chromosomes. The staged==monolith equivalence
# below holds for ANY valid input, and the monolith baseline is computed at
# test time (no pre-saved golden / no session-specific path), so this test
# runs anywhere (CI, clean checkout) not just where a fixture was captured.
set.seed(20240804L)
n_ind <- 40L; n_chr <- 6L; per_chr <- 40L; m <- n_chr * per_chr
ids  <- sprintf("L%03d", seq_len(n_ind))
snps <- sprintf("SNP_%03d", seq_len(m))
.Gmat <- matrix(sample(c(0L, 2L), n_ind * m, replace = TRUE), n_ind, m,
                dimnames = list(NULL, snps))
genotype <- data.frame(NAME = ids, .Gmat, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(
  SNP_code    = snps,
  Chromosome  = rep(seq_len(n_chr), each = per_chr),
  Position_BP = rep(seq_len(per_chr) * 5e5, n_chr),
  stringsAsFactors = FALSE
)
.qtl <- sort(sample(m, 25L)); .beta <- rnorm(length(.qtl))
.gv  <- as.numeric(scale((.Gmat[, .qtl] - 1) %*% .beta))
I <- list(
  pheno     = data.frame(NAME = ids, yield = .gv + rnorm(n_ind),
                         disease = -0.6 * .gv + rnorm(n_ind), stringsAsFactors = FALSE),
  dir       = data.frame(Trait = c("yield", "disease"),
                         Selection_direction = c("increase", "decrease"), stringsAsFactors = FALSE),
  pheno_idx = data.frame(NAME = ids, sel_index = .gv + rnorm(n_ind), stringsAsFactors = FALSE),
  dir_idx   = data.frame(Trait = "sel_index", Selection_direction = "increase", stringsAsFactors = FALSE)
)
base <- list(
  genotype = genotype, marker_map = marker_map, id_col = "NAME",
  map_position_unit = "bp", bp_per_cm = 1e6, n_crosses = 10L,
  progeny = "DH", optimizer = "mip_linear", seed = 1L
)

# stabilizer MUST also drop nested qc$generated_at (wall-clock, volatile)
stab <- function(x) {
  for (k in c("generated_at", "elapsed", "output_files")) x[[k]] <- NULL
  if (!is.null(x$qc)) x$qc$generated_at <- NULL
  x
}

# ng_cp__build_ctx()/ng_run_stage() expect a COMPLETE config (every formal of
# ng_run_cross_prediction resolved -- match.arg fields tolerate NULL/missing
# via match.arg's own "NULL -> choices[1]" rule, but plain-default fields do
# not). Build that complete config the same way the real function does:
# reflect its formals onto a throwaway function whose body is the exact
# `mget(names(formals()))` idiom R/39's driver uses, so unresolved args get
# R's normal lazy-default evaluation instead of being hand-transcribed here.
ng_test_full_config <- function(overrides) {
  f <- function() NULL
  formals(f) <- formals(ng_run_cross_prediction)
  body(f) <- quote(mget(names(formals())))
  do.call(f, overrides)
}

config_single <- ng_test_full_config(c(base, list(
  phenotype = I$pheno[, c("NAME", "yield")],
  trait_direction = I$dir[I$dir$Trait == "yield", ]
)))
config_multi <- ng_test_full_config(c(base, list(
  phenotype = I$pheno,
  trait_direction = I$dir,
  multi_trait_method = "auto"
)))
config_index <- ng_test_full_config(c(base, list(
  phenotype = I$pheno_idx,
  trait_direction = I$dir_idx,
  prediction_mode = "index_as_trait",
  index_col = "sel_index"
)))
config_ldprune <- ng_test_full_config(c(base, list(
  phenotype = I$pheno[, c("NAME", "yield")],
  trait_direction = I$dir[I$dir$Trait == "yield", ],
  ld_pruning = TRUE,
  ld_window = 50,
  ld_r2_threshold = 0.9
)))

drive <- function(config, rd = tempfile("pipe_")) {
  dir.create(rd)
  out <- NULL
  for (s in ng_cp_stage_order()) out <- ng_run_stage(s, rd, config)
  list(result = attr(out, "result"), run_dir = rd)
}

# --- 1. equivalence: staged driving == monolith golden -----------------
configs <- list(single = config_single, multi = config_multi, index = config_index, ldprune = config_ldprune)
run_dirs <- list()
for (nm in names(configs)) {
  driven <- drive(configs[[nm]])
  run_dirs[[nm]] <- driven$run_dir
  # baseline = the one-shot monolith on the SAME config, computed here (a
  # fixed seed makes it deterministic), so staged == monolith is verified
  # without any pre-saved golden.
  gold <- do.call(ng_run_cross_prediction, configs[[nm]])
  ok <- identical(stab(driven$result), stab(gold))
  if (!ok) {
    cat("MISMATCH for config:", nm, "\n")
    print(all.equal(stab(driven$result), stab(gold)))
  }
  stopifnot(ok)
  cat(sprintf("%-8s identical=TRUE\n", nm))
}

# --- 1b. equivalence with check lines configured -------------------------
# The loop above never exercises trait_checks/check_geno/check_records, so it cannot see a
# defect in how ng_run_stage()'s saveRDS()/readRDS() round trip carries the check-reference
# ctx fields across the stage boundary: check_geno (a matrix -- rownames carry the check
# identities), check_records (a nested list), trait_checks (a data.frame) all come from the
# ORIGINAL config and must still be intact when ng_cp__stage_index (R/39) reads them back out
# of the persisted predict.rds ctx; ctc (the exact within-family cross-trait covariance) is
# COMPUTED in stage_predict and must itself survive that same round trip into stage_index.
# Two checked traits (yield, disease, each with its own distinct check line) so p_beat_all_checks
# (the joint Monte Carlo over BOTH checks at once, R/51 + R/33) and ctc both actually participate
# -- a single checked trait would leave p_beat_all_checks absent entirely and never touch ctc.
set.seed(20240805L)
chk_ids <- c("CHK_A", "CHK_B")
chk_geno <- matrix(sample(c(0L, 2L), length(chk_ids) * m, replace = TRUE),
                   nrow = length(chk_ids), dimnames = list(chk_ids, snps))
trait_checks_df <- data.frame(trait = c("yield", "disease"), check = c("CHK_A", "CHK_B"),
                              stringsAsFactors = FALSE)
# adjusted_pheno records: the default ridge fit here never stamps a calibrated reliability, so
# ng_choose_mean_source() resolves both traits onto the phenotypic "adjusted_pheno" source (same
# reasoning as tests/check_reference_invariant.R) -- these are the records actually consulted.
chk_records <- list(
  yield   = list(adjusted_pheno = c(CHK_A = 0.3, CHK_B = -0.2)),
  disease = list(adjusted_pheno = c(CHK_A = -0.1, CHK_B = 0.4))
)
config_checks <- ng_test_full_config(c(base, list(
  phenotype = I$pheno,
  trait_direction = I$dir,
  multi_trait_method = "auto",
  trait_checks = trait_checks_df,
  check_geno = chk_geno,
  check_progeny_size = 200L,
  check_records = chk_records
)))

cat("checks-configured: driving staged pipeline...\n")
driven_checks <- drive(config_checks)
cat("checks-configured: computing one-shot gold...\n")
gold_checks <- do.call(ng_run_cross_prediction, config_checks)

ok_checks <- identical(stab(driven_checks$result), stab(gold_checks))
if (!ok_checks) {
  cat("MISMATCH for config: checks\n")
  print(all.equal(stab(driven_checks$result), stab(gold_checks)))
}
stopifnot(ok_checks)
cat("checks   identical=TRUE\n")

# Named, specific assertions on the check-bearing pieces -- so a future regression here names
# exactly which field diverged rather than just "not identical" (the requirement's ask).
cc_driven <- driven_checks$result$candidate_crosses
cc_gold   <- gold_checks$candidate_crosses
check_cols <- as.vector(outer(c("yield", "disease"),
                              c("_check_value", "_vs_check", "_check_ok", "_p_beat_check"),
                              paste0))
for (cn in c(check_cols, "checks_all_ok", "check_violation", "p_beat_all_checks")) {
  stopifnot(cn %in% names(cc_driven))
  stopifnot(cn %in% names(cc_gold))
  stopifnot(isTRUE(all.equal(cc_driven[[cn]], cc_gold[[cn]], tolerance = 0)))
}

tcr_driven <- driven_checks$result$trait_check_reference
tcr_gold   <- gold_checks$trait_check_reference
stopifnot(!is.null(tcr_driven), !is.null(tcr_gold))
stopifnot(isTRUE(all.equal(tcr_driven$active, tcr_gold$active, tolerance = 0)))
stopifnot(isTRUE(all.equal(tcr_driven$values, tcr_gold$values, tolerance = 0)))
stopifnot(identical(tcr_driven$source, tcr_gold$source))
stopifnot(identical(tcr_driven$progeny_size, tcr_gold$progeny_size))

cat(sprintf(
  "checks-configured equivalence ok: %d per-trait/global check columns + trait_check_reference (active/values/source/progeny_size) identical across staged and one-shot\n",
  length(check_cols) + 3L))

# --- 2. QC blocker gate --------------------------------------------------
# A genotype with a duplicated parent ID makes ng_preflight_input_tables
# flag a blocker-severity issue. ng_cp__stage_qc (R/39) does NOT throw on a
# blocker -- it sets ctx$qc (status="blocker", issues populated) and returns
# early (clean/align skipped), exactly where the monolith's ng_stop used to
# fire. ng_run_stage("qc", ...) therefore returns without error, persists
# qc.rds + qc.json + the manifest with status="blocker", and the downstream
# gate in ng_run_stage (any non-qc stage refuses to run when the persisted
# qc status is "blocker") is what actually stops the pipeline -- when
# ng_run_stage("predict", ...) is called next, it must error.
# The one-shot monolith (ng_run_cross_prediction) still throws on this same
# input -- see the "monolith still throws on blocker" check further below.
dup_geno <- base$genotype
dup_geno[2, "NAME"] <- dup_geno[1, "NAME"]
config_blocker <- ng_test_full_config(modifyList(base, list(
  genotype = dup_geno,
  phenotype = I$pheno[, c("NAME", "yield")],
  trait_direction = I$dir[I$dir$Trait == "yield", ]
)))
# Confirm ng_preflight_input_tables independently flags this input a blocker
# (the "or ng_preflight flags a blocker" branch of the requirement).
qc_direct <- ng_preflight_input_tables(
  geno = ng_run_cp_canonical_id_table(dup_geno, "NAME", "genotype"),
  phenotype = NULL, trait_spec = NULL, marker_map = NULL, ploidy = 2L
)
stopifnot(identical(qc_direct$status, "blocker"))
stopifnot(any(qc_direct$issues$severity == "blocker"))

rd_blocker <- tempfile("pipe_blocker_")
dir.create(rd_blocker)
qc_out <- ng_run_stage("qc", rd_blocker, config_blocker)
stopifnot(identical(qc_out$status, "blocker"))

qc_json_path <- file.path(rd_blocker, "artifacts", "qc.json")
stopifnot(file.exists(qc_json_path))
qc_json <- jsonlite::fromJSON(qc_json_path, simplifyVector = FALSE)
stopifnot(identical(qc_json$status, "blocker"))
stopifnot(length(qc_json$issues) > 0L || (is.data.frame(qc_json$issues) && nrow(qc_json$issues) > 0L))

predict_err <- tryCatch({ ng_run_stage("predict", rd_blocker, config_blocker); NULL }, error = function(e) e)
stopifnot(!is.null(predict_err))
cat("QC blocker gate: qc stage returned status=blocker with qc.json persisted; predict stage errored as expected\n")

# The one-shot monolith must still throw on this same blocker input, with
# the same message the staged qc stage would have raised pre-fix (verbatim
# original ng_stop text, reinstated in the driver's stage loop).
monolith_err <- tryCatch({ do.call(ng_run_cross_prediction, config_blocker); NULL }, error = function(e) e)
stopifnot(!is.null(monolith_err))
stopifnot(grepl("Input QC found blocker issues before cross prediction", conditionMessage(monolith_err), fixed = TRUE))
cat("Monolith blocker gate: ng_run_cross_prediction still throws on blocker input\n")

# --- 3. compute-once probe -----------------------------------------------
fit_log <- tempfile("fit_log_")
writeLines(character(0), fit_log)
invisible(trace(ng_fit_ridge_effects, tracer = quote(cat("FIT\n", file = fit_log, append = TRUE)), print = FALSE))
count_fits <- function() length(readLines(fit_log))

rd_once <- tempfile("pipe_once_")
dir.create(rd_once)
.probe <- ng_run_stage("qc", rd_once, config_single)
stopifnot(count_fits() == 0L)
.probe <- ng_run_stage("predict", rd_once, config_single)
n_after_predict <- count_fits()
stopifnot(n_after_predict > 0L)
.probe <- ng_run_stage("index", rd_once, config_single)
stopifnot(count_fits() == n_after_predict)
.probe <- ng_run_stage("allocate", rd_once, config_single)
stopifnot(count_fits() == n_after_predict)
.probe <- ng_run_stage("rank", rd_once, config_single)
stopifnot(count_fits() == n_after_predict)
rm(.probe)
untrace(ng_fit_ridge_effects)
cat("compute-once probe: fit fired only during predict (n_fits=", n_after_predict, ")\n", sep = "")

# --- 4. per-stage JSON presence + required keys --------------------------
art_dir <- file.path(run_dirs$single, "artifacts")
qc_json <- jsonlite::fromJSON(file.path(art_dir, "qc.json"), simplifyVector = FALSE)
stopifnot(all(c("status", "issues", "putative_duplicates") %in% names(qc_json)))
stopifnot("pairs" %in% names(qc_json$putative_duplicates))

predict_json <- jsonlite::fromJSON(file.path(art_dir, "predict.json"), simplifyVector = FALSE)
stopifnot(all(c("effect_summary", "ld_pruning_report", "n_candidates") %in% names(predict_json)))

index_json <- jsonlite::fromJSON(file.path(art_dir, "index.json"), simplifyVector = FALSE)
stopifnot(all(c("multi_trait_score", "objective") %in% names(index_json)))
stopifnot("method" %in% names(index_json$objective))

allocate_json <- jsonlite::fromJSON(file.path(art_dir, "allocate.json"), simplifyVector = FALSE)
stopifnot(all(c("plan_summary", "parent_use") %in% names(allocate_json)))

manifest_json <- jsonlite::fromJSON(file.path(art_dir, "manifest.json"), simplifyVector = FALSE)
stopifnot(identical(manifest_json$schema, "ng_stage_manifest.v1"))
stopifnot(all(ng_cp_stage_order() %in% names(manifest_json$stages)))
for (s in ng_cp_stage_order()) stopifnot(identical(manifest_json$stages[[s]]$status, if (identical(s, "qc")) "pass" else "done"))

stopifnot(!file.exists(file.path(art_dir, "rank.json")))
stopifnot(file.exists(file.path(art_dir, "rank.rds")))

cat("JSON artifact presence + required keys: OK\n")
cat("ALL STAGED PIPELINE TESTS PASSED\n")
