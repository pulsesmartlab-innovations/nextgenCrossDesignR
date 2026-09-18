# The detached process must be self-describing from its FIRST act.
#
# The frontend will spawn one Rscript with wait = FALSE and then look for job.json. If the
# backend cannot load, or the config is bad, no record is ever written -- and "nothing
# happened" is the worst possible feedback for a run a breeder expects to take hours. So the
# headless entry writes the job record before it does anything expensive, and marks the job
# failed rather than vanishing.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])
root <- ng_test_find_root()

tmp <- file.path(tempdir(), paste0("ngcd_hw_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, recursive = TRUE); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

set.seed(9); n <- 16L; m <- 16L
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

job_dir <- file.path(tmp, "jobs", "job_1")
cfg <- list(workflow = "batch", job_dir = job_dir,
            batch_output_root = job_dir, shared_dir = file.path(tmp, "shared"),
            batch_workers = 1L,
            genotype_file = file.path(tmp, "geno.csv"), genotype_id_col = "NAME",
            phenotype_file = file.path(tmp, "pheno.csv"), phenotype_id_col = "NAME",
            direction_file = file.path(tmp, "dir.csv"), direction_trait_col = "Trait",
            direction_column_col = "Trait", direction_direction_col = "Selection_direction",
            map_file = file.path(tmp, "map.csv"), map_marker_col = "SNP_code",
            map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
            map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
            run_posterior_prediction = FALSE, min_cv_predictive_r2 = -1, seed = 1L)
cfg_path <- file.path(tmp, "config.json"); res_path <- file.path(tmp, "result.json")
jsonlite::write_json(cfg, cfg_path, auto_unbox = TRUE)

status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(file.path(root, "tools", "run_cross_prediction_json.R")),
                    shQuote(cfg_path), shQuote(res_path)),
                  stdout = NULL, stderr = NULL)
stopifnot(identical(as.integer(status), 0L))

# The job record exists and says the batch completed.
stopifnot(file.exists(file.path(job_dir, "job.json")))
s <- ng_job_status(job_dir)
stopifnot(identical(s$state, "finished"))
stopifnot(identical(as.integer(s$n_done), 2L))
stopifnot(identical(sort(s$traits$trait), c("A", "B")))

# --- a job that cannot run still leaves a record ---------------------------------------------
bad_dir <- file.path(tmp, "jobs", "job_bad")
bad <- cfg; bad$job_dir <- bad_dir; bad$batch_output_root <- bad_dir
bad$genotype_file <- file.path(tmp, "does_not_exist.csv")
bad_cfg <- file.path(tmp, "bad.json"); bad_res <- file.path(tmp, "bad_result.json")
jsonlite::write_json(bad, bad_cfg, auto_unbox = TRUE)
invisible(system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(file.path(root, "tools", "run_cross_prediction_json.R")),
                    shQuote(bad_cfg), shQuote(bad_res)),
                  stdout = NULL, stderr = NULL))
stopifnot(file.exists(file.path(bad_dir, "job.json")))
sb <- ng_job_status(bad_dir)
stopifnot(identical(sb$state, "failed"))

cat("PASS: job_workflow_runs_headlessly\n")
