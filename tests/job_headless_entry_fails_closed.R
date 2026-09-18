# A job that cannot run must END, and a job that runs must be READABLE.
#
# Two defects in the headless entry, both of which leave a breeder with a job they cannot act
# on and no supported way to clear it:
#
#   1. CONFIG VALIDATION LEFT A JOB QUEUED FOREVER. ng_job_create() ran before the
#      unknown-config-key stop(), which sat outside the tryCatch that marks a job failed. Such
#      a job stayed "queued" -- a state ng_job_status() never reinterprets (there is no
#      heartbeat to go stale) and ng_job_prune() protects unconditionally. ng_job_mark() is not
#      exported, so nothing could ever collect it. It would sit in the Jobs tab claiming it was
#      about to start.
#
#   2. job_dir AND batch_output_root WERE SILENTLY DECOUPLED. Workers write
#      <output_root>/<trait>/status.json; ng_job_status() scans <job_dir>. The entry defaulted
#      batch_output_root to dirname(result_path)/batch, so a caller who passed only job_dir got
#      a job whose traits were permanently invisible -- the record went "running" then
#      "finished" with nothing under it, forever.
#
# Driven through the real Rscript entry point, not by calling the functions: both defects are
# in the ORDER of that script, and an in-process test would not see the order at all.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])
root <- ng_test_find_root()

tmp <- file.path(tempdir(), paste0("ngcd_fc_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, recursive = TRUE); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

set.seed(11); n <- 16L; m <- 16L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
write.csv(data.frame(NAME = ids, g, check.names = FALSE), file.path(tmp, "geno.csv"), row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2)), file.path(tmp, "pheno.csv"), row.names = FALSE)
write.csv(data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                     Position_cM = rep(seq(0, 45, length.out = m / 4), times = 4)),
          file.path(tmp, "map.csv"), row.names = FALSE)
write.csv(data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
          file.path(tmp, "dir.csv"), row.names = FALSE)

base_cfg <- list(workflow = "batch", shared_dir = file.path(tmp, "shared"), batch_workers = 1L,
  genotype_file = file.path(tmp, "geno.csv"), genotype_id_col = "NAME",
  phenotype_file = file.path(tmp, "pheno.csv"), phenotype_id_col = "NAME",
  direction_file = file.path(tmp, "dir.csv"), direction_trait_col = "Trait",
  direction_column_col = "Trait", direction_direction_col = "Selection_direction",
  map_file = file.path(tmp, "map.csv"), map_marker_col = "SNP_code",
  map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  run_posterior_prediction = FALSE, min_cv_predictive_r2 = -1, effect_gate = "off", seed = 1L)

run_headless <- function(tag, cfg) {
  cfg_path <- file.path(tmp, paste0(tag, "_config.json"))
  res_path <- file.path(tmp, paste0(tag, "_result.json"))
  jsonlite::write_json(cfg, cfg_path, auto_unbox = TRUE)
  st <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                 c(shQuote(file.path(root, "tools", "run_cross_prediction_json.R")),
                   shQuote(cfg_path), shQuote(res_path)),
                 stdout = NULL, stderr = NULL))
  as.integer(st)
}

# --- 1. an unknown config key ends the job FAILED, not queued ---------------------------------
bad_dir <- file.path(tmp, "jobs", "job_badkey")
bad <- base_cfg; bad$job_dir <- bad_dir; bad$batch_output_root <- bad_dir
bad$not_a_runner_argument <- 42L
status <- run_headless("badkey", bad)
stopifnot(!identical(status, 0L))               # the process still reports failure
stopifnot(file.exists(file.path(bad_dir, "job.json")))
sb <- ng_job_status(bad_dir)
# The whole point: NOT "queued". A queued job is never reinterpreted and never pruned, so it
# would be unclearable through any supported route.
stopifnot(identical(sb$state, "failed"))
rec <- jsonlite::fromJSON(file.path(bad_dir, "job.json"), simplifyVector = TRUE)
stopifnot(grepl("not_a_runner_argument", as.character(rec$error_message), fixed = TRUE))

# --- 2. batch_output_root follows job_dir when it is not given ---------------------------------
# The config below names job_dir and NOTHING about where outputs go. Before the fix, the
# traits landed beside result.json and ng_job_status() -- which scans job_dir -- saw none.
solo_dir <- file.path(tmp, "jobs", "job_solo")
solo <- base_cfg; solo$job_dir <- solo_dir
stopifnot(is.null(solo$batch_output_root))
stopifnot(identical(run_headless("solo", solo), 0L))
s <- ng_job_status(solo_dir)
stopifnot(identical(s$state, "finished"))
stopifnot(identical(as.integer(s$n_traits), 2L))
stopifnot(identical(as.integer(s$n_done), 2L))
stopifnot(identical(sort(s$traits$trait), c("A", "B")))
# The deliverables are inside the job directory, which is the whole contract.
stopifnot(file.exists(file.path(solo_dir, "A", "result.json")))
stopifnot(file.exists(file.path(solo_dir, "manifest.json")))

# --- 3. a MISMATCH is refused, loudly, rather than producing an invisible job ------------------
mm_dir <- file.path(tmp, "jobs", "job_mismatch")
mm <- base_cfg; mm$job_dir <- mm_dir; mm$batch_output_root <- file.path(tmp, "elsewhere")
stopifnot(!identical(run_headless("mismatch", mm), 0L))
sm <- ng_job_status(mm_dir)
stopifnot(identical(sm$state, "failed"))
mrec <- jsonlite::fromJSON(file.path(mm_dir, "job.json"), simplifyVector = TRUE)
stopifnot(grepl("same directory as", as.character(mrec$error_message), fixed = TRUE))
# And it refused BEFORE doing the work, so nothing was written to the wrong place.
stopifnot(!dir.exists(file.path(tmp, "elsewhere")))

cat("PASS: job_headless_entry_fails_closed\n")
