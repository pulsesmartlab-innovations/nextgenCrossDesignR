# A batch must produce EXACTLY what running each trait on its own produces.
#
# This is the claim the whole batch design rests on. If a trait's numbers depended on
# whether it was scored alongside others, the batch would be a different analysis wearing
# the same name -- and the breeder, who ran the 17 by hand once already, would have no way
# to tell which set was right.
#
# It is credible only because two things were built that way on purpose: the candidate
# parent set is fixed at QC (a trait's missing phenotypes move only its TRAINING set), and
# the ridge-lambda fold partition is keyed on the run seed rather than the trait's row
# position. This test is what turns that reasoning into a fact that stays true.
#
# Compared on the full candidate table, not on the selected plan. A plan-only check would
# pass even if a variance moved in the twelfth decimal and merely reordered ties.
# helper_load.R resolves the package root and prepends the project library ONLY when it was
# built for this platform -- a blunt .libPaths(".Rlib") prepend silently shadows working
# packages with ones built for another OS, and a mirai daemon inherits the parent's library
# paths, so it fails inside a worker with an opaque "shared object not found".
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(707)
n <- 24L; m <- 32L
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n)); rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))
trait_of <- function(s) { set.seed(s); b <- rep(0, m); b[1:5] <- rnorm(5, sd = 1.5)
  as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.3) }
traits <- c("YIELD", "PROTEIN", "HEIGHT")
pheno <- data.frame(NAME = ids, YIELD = trait_of(1L), PROTEIN = trait_of(2L),
                    HEIGHT = trait_of(3L), stringsAsFactors = FALSE)
gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
map <- data.frame(SNP_code = colnames(geno), Chromosome = rep(1:4, length.out = m),
                  Position_cM = rep(seq(0, 100, length.out = m / 4), times = 4),
                  stringsAsFactors = FALSE)
dir <- data.frame(Trait = traits, Selection_direction = "increase", stringsAsFactors = FALSE)

base_config <- list(
  genotype = gdf, genotype_id_col = "NAME",
  phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = dir, direction_trait_col = "Trait",
  direction_column_col = "Trait", direction_direction_col = "Selection_direction",
  marker_map = map, map_marker_col = "SNP_code", map_chr_col = "Chromosome",
  map_pos_cm_col = "Position_cM", map_position_unit = "cM",
  progeny = "RIL", parent_type = "ril", n_crosses = 5L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -Inf, seed = 20260917L)

# --- each trait on its own, exactly as a breeder runs it today --------------------------
alone <- lapply(traits, function(tr) {
  suppressWarnings(do.call(ng_run_cross_prediction,
                           c(base_config, list(traits_to_use = tr))))
})
names(alone) <- traits

# --- the same three traits as one batch -------------------------------------------------
root <- file.path(tempdir(), paste0("ngcd_batch_", as.integer(runif(1, 1, 1e6))))
b <- suppressWarnings(ng_run_cross_prediction_batch(base_config, output_root = root,
                                                    batch_workers = 2L))
stopifnot(inherits(b, "ng_cross_prediction_batch"))
stopifnot(identical(length(b$jobs), 3L))
for (j in b$jobs) {
  if (!identical(j$status, "ok")) stop("job ", j$id, " failed: ", j$error_message)
}

# --- the comparison ---------------------------------------------------------------------
# Read each job back from the file the breeder is actually handed, not from an in-memory
# object -- the result.json IS the deliverable, so that is what has to match.
# Both sides go through the SAME writer before being compared. result.json is emitted at
# digits = 8, so comparing an in-memory solo result against a parsed batch file would fail
# on rounding alone and tell us nothing about whether the analyses agree.
num_cols <- function(df) df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
solo_dir <- file.path(root, "_solo"); dir.create(solo_dir, recursive = TRUE, showWarnings = FALSE)
read_env <- function(p) jsonlite::fromJSON(p, simplifyVector = TRUE, simplifyDataFrame = TRUE)

for (tr in traits) {
  path <- file.path(root, tr, "result.json")
  stopifnot(file.exists(path))
  got <- read_env(path)
  stopifnot(isTRUE(got$ok), identical(got$schema, "ng_run_result.v1"))

  solo <- read_env(ng_write_result_json(
    ng_run_result_envelope(alone[[tr]]), file.path(solo_dir, paste0(tr, ".json"))))
  a <- solo$candidate_crosses
  z <- got$candidate_crosses
  a <- a[order(a$parent1, a$parent2), , drop = FALSE]
  z <- z[order(z$parent1, z$parent2), , drop = FALSE]
  stopifnot(identical(nrow(a), nrow(z)))
  shared <- intersect(names(num_cols(a)), names(num_cols(z)))
  stopifnot(length(shared) >= 4L)
  for (cl in shared) {
    av <- as.numeric(a[[cl]]); zv <- as.numeric(z[[cl]])
    if (!isTRUE(all.equal(av, zv, tolerance = 0))) {
      stop("batch and solo disagree for trait ", tr, " on column ", cl)
    }
  }
  # The plan itself, and the model-quality evidence beside it.
  sa <- solo$selected_crosses; sz <- got$selected_crosses
  stopifnot(identical(paste(sa$parent1, sa$parent2), paste(sz$parent1, sz$parent2)))
  stopifnot(isTRUE(all.equal(as.numeric(solo$effect_summary$cv_predictive_r2),
                             as.numeric(got$effect_summary$cv_predictive_r2), tolerance = 0)))
  # The cross-mean basis is the decision the marker evidence drove; if a batch could change
  # it, a plan could silently switch between GEBV and phenotypic mid-parent.
  stopifnot(identical(as.character(solo$effect_summary$mean_source),
                      as.character(got$effect_summary$mean_source)))
}

# --- the manifest has to explain the run ------------------------------------------------
mf <- jsonlite::fromJSON(b$manifest_path, simplifyVector = TRUE)
stopifnot(identical(mf$schema, "ng_batch_manifest.v1"))
stopifnot(identical(as.integer(mf$n_jobs), 3L))
stopifnot(is.character(mf$worker_basis), nzchar(mf$worker_basis))

unlink(root, recursive = TRUE)
cat("PASS: batch_matches_individual_runs\n")
