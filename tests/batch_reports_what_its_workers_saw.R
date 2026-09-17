# A batch must surface its workers' warnings, and must survive one job failing.
#
# TWO FAILURE MODES, ONE TEST, because they share a cause: a job runs in another process,
# and anything it does not deliberately carry home is lost.
#
# 1. WARNINGS. R/52_advisories.R:20 states the rule -- a warning() raised in a worker never
#    reaches the parent's calling handlers, so it never reaches the JSON envelope, and the
#    package calls a silent advisory channel "the exact defect class this work exists to
#    remove". The gate advisory that a trait fell back to the phenotypic mid-parent is
#    precisely the sort of thing a breeder must not lose: it changes how the plan is read.
#    So jobs return conditions as data and the parent re-emits them.
#
# 2. FAILURE. One bad trait out of seventeen must not cost the other sixteen. A batch may
#    run for hours; an all-or-nothing batch would be worse than the manual approach it
#    replaces, because at least the manual runs kept what had already finished.
#
# Both are asserted with the trait NAME attached. An unattributed warning in a 17-job batch
# tells a breeder nothing they can act on.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(929)
n <- 22L; m <- 28L
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n)); rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))
strong <- { set.seed(1L); b <- rep(0, m); b[1:6] <- rnorm(6, sd = 2)
            as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.1) }
# A trait the markers cannot predict -> the gate falls back and warns.
noise  <- { set.seed(2L); rnorm(n) }
# A trait with no usable data at all -> its job fails, and only its job.
empty  <- rep(NA_real_, n)

pheno <- data.frame(NAME = ids, GOOD = strong, WEAK = noise, EMPTY = empty,
                    stringsAsFactors = FALSE)
gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
map <- data.frame(SNP_code = colnames(geno), Chromosome = rep(1:4, length.out = m),
                  Position_cM = rep(seq(0, 100, length.out = m / 4), times = 4),
                  stringsAsFactors = FALSE)
dir <- data.frame(Trait = c("GOOD", "WEAK", "EMPTY"), Selection_direction = "increase",
                  stringsAsFactors = FALSE)

cfg <- list(
  genotype = gdf, genotype_id_col = "NAME", phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = dir, direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction", marker_map = map,
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 4L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  trait_value_metric = "mean", effect_gate = "on", seed = 20260917L)

root <- file.path(tempdir(), "ngcd_batch_conditions")
unlink(root, recursive = TRUE)
seen <- character(0)
b <- withCallingHandlers(
  ng_run_cross_prediction_batch(cfg, output_root = root, batch_workers = 2L),
  warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })

by_id <- function(id) Filter(function(j) identical(j$id, id), b$jobs)[[1L]]

# --- failure isolation ------------------------------------------------------------------
stopifnot(identical(by_id("EMPTY")$status, "error"))
stopifnot(is.character(by_id("EMPTY")$error_message), nzchar(by_id("EMPTY")$error_message))
# The other two completed and wrote their deliverables. That is the whole point.
stopifnot(identical(by_id("GOOD")$status, "ok"), identical(by_id("WEAK")$status, "ok"))
stopifnot(file.exists(file.path(root, "GOOD", "result.json")))
stopifnot(file.exists(file.path(root, "WEAK", "result.json")))
stopifnot(!file.exists(file.path(root, "EMPTY", "result.json")))

# --- the manifest records the failure rather than hiding it -----------------------------
mf <- jsonlite::fromJSON(b$manifest_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
statuses <- vapply(mf$jobs, function(j) j$status, character(1))
stopifnot(identical(sort(statuses), sort(c("ok", "ok", "error"))))

# --- worker warnings reached the parent, attributed ---------------------------------------
# WEAK cannot be predicted from these markers, so the gate falls back and says so. That
# warning was raised inside a daemon; it is here only because the job carried it home.
stopifnot(length(seen) > 0L)
stopifnot(any(grepl("WEAK", seen, fixed = TRUE)))
stopifnot(any(grepl("\\[WEAK\\]", seen)))

cat("PASS: batch_reports_what_its_workers_saw\n")
