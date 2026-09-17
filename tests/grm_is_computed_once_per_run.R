# The genomic relationship matrix must be built ONCE per run, not once per trait.
#
# ng_parent_kinship() is O(n^2 * m) in parents and markers, and it depends only on the
# genotype matrix and the GRM method -- nothing about a trait enters it. Yet a run
# rebuilt it inside ng_score_crosses() for every trait (R/03_metrics.R:223) and then
# once more in the allocate stage (R/39:1951). On a 17-trait panel that is 18 builds of
# one identical matrix.
#
# ng_score_crosses() already anticipated this: it has both `pairs` and `parent_kinship`
# formals, and the comment above the call site says they exist so a caller looping over
# traits does not recompute them. The runner simply never passed either. The cost is
# pure waste -- the matrices are bit-identical by construction, so nothing downstream
# can tell the difference.
#
# This matters most for the batch runner, where the whole point is to share
# trait-independent work across jobs; but it is worth fixing on its own, because every
# ordinary multi-trait run pays it today.
#
# Counted by rebinding the function in .GlobalEnv -- R/load.R sys.source()s every file
# there, so package functions resolve it lexically through the global environment and a
# counting wrapper is seen by every caller.
# helper_load.R resolves the package root and prepends the project library ONLY when it was
# built for this platform -- a blunt .libPaths(".Rlib") prepend silently shadows working
# packages with ones built for another OS, and a mirai daemon inherits the parent's library
# paths, so it fails inside a worker with an opaque "shared object not found".
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(404)
n <- 26L; m <- 36L
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))

trait_of <- function(seed) {
  set.seed(seed)
  b <- rep(0, m); b[1:6] <- rnorm(6, sd = 1.5)
  as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.3)
}
pheno <- data.frame(NAME = ids, YIELD = trait_of(1L), PROTEIN = trait_of(2L),
                    HEIGHT = trait_of(3L), stringsAsFactors = FALSE)
gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
map <- data.frame(SNP_code = colnames(geno),
                  Chromosome = rep(1:4, length.out = m),
                  Position_cM = rep(seq(0, 120, length.out = m / 4), times = 4),
                  stringsAsFactors = FALSE)
dir <- data.frame(Trait = c("YIELD", "PROTEIN", "HEIGHT"),
                  Selection_direction = "increase", stringsAsFactors = FALSE)

run_and_count <- function(traits) {
  orig <- ng_parent_kinship
  calls <- 0L
  assign("ng_parent_kinship",
         function(...) { calls <<- calls + 1L; orig(...) },
         envir = .GlobalEnv)
  on.exit(assign("ng_parent_kinship", orig, envir = .GlobalEnv), add = TRUE)
  res <- ng_run_cross_prediction(
    genotype = gdf, genotype_id_col = "NAME",
    phenotype = pheno, phenotype_id_col = "NAME", traits_to_use = traits,
    trait_direction = dir, direction_trait_col = "Trait",
    direction_column_col = "Trait", direction_direction_col = "Selection_direction",
    marker_map = map, map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_cm_col = "Position_cM", map_position_unit = "cM",
    progeny = "RIL", parent_type = "ril", n_crosses = 6L,
    write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
    seed = 20260917L)
  list(result = res, calls = calls)
}

three <- run_and_count(c("YIELD", "PROTEIN", "HEIGHT"))
one   <- run_and_count("YIELD")

# The GRM depends on genotype alone, so the trait COUNT must not change how many times
# it is built. Asserting equality (not just "few") is what makes this a real constraint:
# it fails the moment a per-trait build creeps back in.
stopifnot(identical(three$calls, one$calls))

# And that shared count is one. A second build would be a bit-identical matrix.
stopifnot(identical(three$calls, 1L))

# Reuse must not change any number. The kinship column is the one the allocate stage
# consumes, so it is the direct witness that the shared matrix reached the optimizer.
sc <- three$result$selected_crosses
stopifnot(is.data.frame(sc), nrow(sc) > 0L, "pair_kinship" %in% names(sc))
stopifnot(all(is.finite(sc$pair_kinship)))

cat("PASS: grm_is_computed_once_per_run\n")
