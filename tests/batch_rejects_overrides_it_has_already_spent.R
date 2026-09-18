# A job may not override a setting the batch has already spent.
#
# The batch runs quality control and the shared marker-effect setup ONCE and hands the
# result to every job. A job that changed the duplicate threshold, the LD pruning window or
# the GRM method would therefore be scored against inputs it did not ask for.
#
# There are only three possible behaviours and two of them are wrong. Silently ignoring the
# override returns numbers that quietly contradict the request. Silently honouring it means
# the shared artefact was never shared, and the batch is 17 separate runs wearing a batch's
# name. So: refuse, name the key, and say what to do instead. The limitation is real and it
# is stated rather than discovered.
#
# Checked BEFORE any job starts. A batch may run for hours; learning at hour six that job
# twelve was inadmissible is not an acceptable way to find out.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# The list must name real settings. A typo, or a formal renamed by a later release, would
# silently stop protecting the key it was written for -- the guard would still "pass"
# while admitting exactly the override it exists to reject. The same check covers
# ng_cp__batch_artifact_keys (R/54), the artefact-content sibling of this override-guard
# list -- see the comment on ng_shared_artifact_key() in R/55 for why the two lists exist
# and are deliberately different.
formal_names <- names(formals(ng_run_cross_prediction))
unknown <- setdiff(ng_cp__batch_shared_keys, formal_names)
if (length(unknown)) {
  stop("ng_cp__batch_shared_keys names settings that are not runner arguments (renamed or ",
       "removed?): ", paste(unknown, collapse = ", "))
}
unknown_artifact <- setdiff(ng_cp__batch_artifact_keys, formal_names)
if (length(unknown_artifact)) {
  stop("ng_cp__batch_artifact_keys names settings that are not runner arguments (renamed or ",
       "removed?): ", paste(unknown_artifact, collapse = ", "))
}
# And it must actually cover the things the batch computes once.
for (k in c("duplicate_threshold", "ld_window", "grm_method", "training_genotype", "genotype")) {
  stopifnot(k %in% ng_cp__batch_shared_keys)
}
# traits_to_use must NOT be in the override guard -- selecting the trait is the entire point
# of a job.
stopifnot(!("traits_to_use" %in% ng_cp__batch_shared_keys))
# Nor may the override guard swallow ordinary per-job settings, which the breeder was
# promised they could vary.
for (k in c("n_crosses", "min_cv_predictive_r2", "trait_value_metric", "selection_prop")) {
  stopifnot(!(k %in% ng_cp__batch_shared_keys))
}
# But the ARTEFACT key must cover exactly the settings that feed trait_spec at QC time --
# traits_to_use included, precisely because it is overridable per job yet still bakes into
# the persisted trait_spec (ng_run_cp_trait_spec()/ng_run_cp_index_spec(), R/39).
for (k in c("traits_to_use", "trait_weights", "prediction_mode", "index_col", "index_direction")) {
  stopifnot(k %in% ng_cp__batch_artifact_keys)
}
# The artefact set is the shared-guard set PLUS those five, not a disjoint list.
stopifnot(identical(sort(ng_cp__batch_artifact_keys),
                    sort(unique(c(ng_cp__batch_shared_keys, "traits_to_use", "trait_weights",
                                 "prediction_mode", "index_col", "index_direction")))))

set.seed(1313)
n <- 20L; m <- 24L
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n)); rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))
trait_of <- function(s) { set.seed(s); b <- rep(0, m); b[1:5] <- rnorm(5, sd = 1.5)
  as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.3) }
pheno <- data.frame(NAME = ids, A = trait_of(1L), B = trait_of(2L), stringsAsFactors = FALSE)
gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
map <- data.frame(SNP_code = colnames(geno), Chromosome = rep(1:4, length.out = m),
                  Position_cM = rep(seq(0, 90, length.out = m / 4), times = 4),
                  stringsAsFactors = FALSE)
dir <- data.frame(Trait = c("A", "B"), Selection_direction = "increase", stringsAsFactors = FALSE)
cfg <- list(
  genotype = gdf, genotype_id_col = "NAME", phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = dir, direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction", marker_map = map,
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 4L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -Inf, seed = 20260917L)
root <- file.path(tempdir(), "ngcd_batch_reject"); unlink(root, recursive = TRUE)

bad <- list(A = list(id = "A", traits = "A"),
            B = list(id = "B", traits = "B", overrides = list(duplicate_threshold = 0.9)))
err <- tryCatch({ ng_run_cross_prediction_batch(cfg, jobs = bad, output_root = root,
                                                batch_workers = 1L); NULL },
                error = function(e) conditionMessage(e))
stopifnot(!is.null(err))
stopifnot(grepl("duplicate_threshold", err, fixed = TRUE))   # names the offending key
stopifnot(grepl("'B'", err, fixed = TRUE))                   # and the job it came from
stopifnot(grepl("own batch", err, fixed = TRUE))             # and the way forward
# Refused before doing the work, not after.
stopifnot(!file.exists(file.path(root, "A", "result.json")))

# --- a legitimate per-trait override is honoured ----------------------------------------
# "Some settings vary per trait" is the requirement this batch was built for: a breeder may
# hold a different acceptability bar for a well-phenotyped trait than for a sparse one.
good <- list(A = list(id = "A", traits = "A", overrides = list(n_crosses = 3L)),
             B = list(id = "B", traits = "B", overrides = list(n_crosses = 6L)))
unlink(root, recursive = TRUE)
b <- suppressWarnings(ng_run_cross_prediction_batch(cfg, jobs = good, output_root = root,
                                                    batch_workers = 1L))
for (j in b$jobs) if (!identical(j$status, "ok")) stop(j$id, ": ", j$error_message)
nsel <- function(id) nrow(jsonlite::fromJSON(file.path(root, id, "result.json"),
                                             simplifyDataFrame = TRUE)$selected_crosses)
stopifnot(identical(nsel("A"), 3L), identical(nsel("B"), 6L))

unlink(root, recursive = TRUE)
cat("PASS: batch_rejects_overrides_it_has_already_spent\n")
