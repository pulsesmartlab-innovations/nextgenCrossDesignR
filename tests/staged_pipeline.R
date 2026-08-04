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
