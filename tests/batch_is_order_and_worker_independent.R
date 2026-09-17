# A batch's results must not depend on HOW it was executed.
#
# Same jobs, different worker counts, different job order: the numbers must be identical.
# A breeder who reruns a batch on a smaller machine, or whose container is given fewer
# cores that day, must get the same crossing plan -- otherwise "the batch said so" means
# nothing.
#
# THIS TEST EXISTS BECAUSE IT CAUGHT A REAL DEFECT. mirai sets its daemons to
# L'Ecuyer-CMRG so concurrent workers draw independent streams. set.seed(s) does not mean
# the same thing under two generators, so the ridge-lambda CV fold partition -- seeded from
# the run seed -- came out different inside a daemon than in the parent. On this fixture the
# batch selected lambda = 8.86 where the same data run serially selected 1.62, and every
# marker-derived quantity moved with it: the variance, the usefulness, the index, the plan.
#
# Nothing errored. There was no warning. The only way to see it was to compare.
#
# The fix pins the daemons' RNG kind to the parent's, which is safe precisely because this
# package never relies on ambient randomness: every stochastic step is seeded explicitly
# from a derived seed. If a future mirai release changes its defaults again, this fails.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(818)
n <- 22L; m <- 28L
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

cfg <- list(
  genotype = gdf, genotype_id_col = "NAME", phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = dir, direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction", marker_map = map,
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 5L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -Inf, seed = 20260917L)

run_batch <- function(tag, jobs, workers) {
  root <- file.path(tempdir(), paste0("ngcd_wi_", tag))
  unlink(root, recursive = TRUE)
  b <- suppressWarnings(ng_run_cross_prediction_batch(cfg, jobs = jobs, output_root = root,
                                                      batch_workers = workers))
  for (j in b$jobs) if (!identical(j$status, "ok")) stop(tag, "/", j$id, ": ", j$error_message)
  out <- lapply(traits, function(tr)
    jsonlite::fromJSON(file.path(root, tr, "result.json"), simplifyDataFrame = TRUE))
  names(out) <- traits
  unlink(root, recursive = TRUE)
  out
}
as_jobs <- function(tr) { j <- lapply(tr, function(t) list(id = t, traits = t)); names(j) <- tr; j }

serial   <- run_batch("serial",   as_jobs(traits),      workers = 1L)
parallel <- run_batch("parallel", as_jobs(traits),      workers = 3L)
reversed <- run_batch("reversed", as_jobs(rev(traits)), workers = 2L)

compare <- function(x, y, what) {
  for (tr in traits) {
    a <- x[[tr]]$candidate_crosses; z <- y[[tr]]$candidate_crosses
    a <- a[order(a$parent1, a$parent2), , drop = FALSE]
    z <- z[order(z$parent1, z$parent2), , drop = FALSE]
    nm <- intersect(names(a)[vapply(a, is.numeric, logical(1))],
                    names(z)[vapply(z, is.numeric, logical(1))])
    stopifnot(length(nm) >= 4L)
    for (cl in nm) {
      if (!isTRUE(all.equal(as.numeric(a[[cl]]), as.numeric(z[[cl]]), tolerance = 0))) {
        stop(what, " changed ", tr, "/", cl)
      }
    }
    # The lambda is the quantity the RNG defect actually moved, so name it explicitly
    # rather than trusting it to fall out of the column sweep.
    stopifnot(isTRUE(all.equal(as.numeric(x[[tr]]$effect_summary$ridge_lambda),
                               as.numeric(y[[tr]]$effect_summary$ridge_lambda), tolerance = 0)))
    stopifnot(identical(paste(x[[tr]]$selected_crosses$parent1, x[[tr]]$selected_crosses$parent2),
                        paste(y[[tr]]$selected_crosses$parent1, y[[tr]]$selected_crosses$parent2)))
  }
}
compare(serial, parallel, "worker count")
compare(serial, reversed, "job order")

cat("PASS: batch_is_order_and_worker_independent\n")
