# The headless JSON contract (tools/run_cross_prediction_json.R) must surface the
# backend's user-facing messages so the frontend can DISPLAY them:
#   - a residual-het RIL without phase -> a structured blocker
#   - a DH/inbred parent carrying het -> a structured
#     result$status = "error" with result$error$message, NOT an uncaught crash.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])
if (!requireNamespace("jsonlite", quietly = TRUE)) { cat("jsonlite unavailable; skipping\n"); quit(save = "no") }

runner <- file.path(root, "tools", "run_cross_prediction_json.R")  # `root` set by helper_load
stopifnot(file.exists(runner))

td <- tempfile("jsonmsg_"); dir.create(td)
set.seed(1); n <- 20L; m <- 60L; ids <- sprintf("L%03d", seq_len(n)); snps <- sprintf("S%03d", seq_len(m))
G <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), n, m); G[1, 1:3] <- 1L   # line 1 het 3/60 = 5%
colnames(G) <- snps
write.csv(data.frame(NAME = ids, G, check.names = FALSE), file.path(td, "geno.csv"), row.names = FALSE)
qtl <- sort(sample(m, 12)); gv <- as.numeric(scale((G[, qtl] - 1) %*% rnorm(12)))
write.csv(data.frame(NAME = ids, yield = gv + rnorm(n)), file.path(td, "pheno.csv"), row.names = FALSE)
write.csv(data.frame(SNP_code = snps, Chromosome = rep(1:3, each = 20), Position_BP = rep(1:20, 3) * 5e5),
          file.path(td, "map.csv"), row.names = FALSE)
write.csv(data.frame(Trait = "yield", Selection_direction = "increase"), file.path(td, "dir.csv"), row.names = FALSE)
base <- list(genotype_file = file.path(td, "geno.csv"), phenotype_file = file.path(td, "pheno.csv"),
             map_file = file.path(td, "map.csv"), direction_file = file.path(td, "dir.csv"),
             id_col = "NAME", map_position_unit = "bp", bp_per_cm = 1e6, n_crosses = 5L, progeny = "RIL", seed = 1L)
jsonlite::write_json(c(base, list(parent_type = "ril")), file.path(td, "cfg_ril.json"), auto_unbox = TRUE)
jsonlite::write_json(c(base, list(parent_type = "dh")),  file.path(td, "cfg_dh.json"),  auto_unbox = TRUE)

run_json <- function(cfg, out) {
  system2("Rscript", c(shQuote(runner), shQuote(cfg), shQuote(out)),
          stdout = FALSE, stderr = FALSE, env = "NGCD_SKIP_CPP=1")
  jsonlite::fromJSON(out, simplifyVector = TRUE, simplifyDataFrame = FALSE)
}

# RIL + het without phase: structured error, because the inbred a'Ra kernel is biased.
r_ril <- run_json(file.path(td, "cfg_ril.json"), file.path(td, "res_ril.json"))
stopifnot(identical(r_ril$status, "error"))
stopifnot(!is.null(r_ril$error_message), grepl("phased_haplotypes", r_ril$error_message, fixed = TRUE))

# DH + het: structured error result (not a crash) with the blocker message
r_dh <- run_json(file.path(td, "cfg_dh.json"), file.path(td, "res_dh.json"))
stopifnot(identical(r_dh$status, "error"))
stopifnot(!is.null(r_dh$error_message), grepl("BLOCKED", r_dh$error_message))

# --- staged pipeline path: same messages must surface per-stage ----------
run_stage <- function(cfg_file, stage, run_dir) {
  cfg <- jsonlite::fromJSON(cfg_file, simplifyVector = TRUE, simplifyDataFrame = TRUE)
  cfg$workflow <- "stage"; cfg$stage <- stage
  cf <- file.path(run_dir, paste0("cfg_", stage, ".json"))
  out <- file.path(run_dir, paste0("stage_", stage, ".json"))
  jsonlite::write_json(cfg, cf, auto_unbox = TRUE)
  system2("Rscript", c(shQuote(runner), shQuote(cf), shQuote(out)),
          stdout = FALSE, stderr = FALSE, env = "NGCD_SKIP_CPP=1")
  jsonlite::fromJSON(out, simplifyVector = TRUE, simplifyDataFrame = FALSE)
}
# RIL: the missing-phase blocker surfaces in the predict stage.
rd <- tempfile("stg_ril_"); dir.create(rd)
invisible(run_stage(file.path(td, "cfg_ril.json"), "qc", rd))
sp <- run_stage(file.path(td, "cfg_ril.json"), "predict", rd)
stopifnot(identical(sp$status, "error"), isTRUE(sp$error),
          !is.null(sp$error_message), grepl("phased_haplotypes", sp$error_message, fixed = TRUE))
# DH: the blocker surfaces as a structured error in the predict stage
rd2 <- tempfile("stg_dh_"); dir.create(rd2)
invisible(run_stage(file.path(td, "cfg_dh.json"), "qc", rd2))
sd <- run_stage(file.path(td, "cfg_dh.json"), "predict", rd2)
stopifnot(identical(sd$status, "error"), isTRUE(sd$error),
          !is.null(sd$error_message), grepl("BLOCKED", sd$error_message))

cat("run_cross_prediction_json_messages: phase and fixed-line blockers surfaced in BOTH monolith and staged JSON contracts\n")
