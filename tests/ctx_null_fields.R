ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- the setter preserves NULL where $<- would delete the key ---------------
ctx <- list(a = 1)
ctx <- ng_ctx_put(ctx, b = NULL, c = 3)
stopifnot(all(c("a", "b", "c") %in% names(ctx)))
stopifnot(is.null(ctx$b), ctx$c == 3)
# and the binding survives the list2env round trip a stage performs
f <- function(ctx) { list2env(ctx, environment()); is.null(b) }
stopifnot(isTRUE(f(ctx)))
# overwriting an existing key with NULL keeps it, too
ctx <- ng_ctx_put(ctx, c = NULL)
stopifnot("c" %in% names(ctx), is.null(ctx$c))

cat("ctx setter ok\n")

# --- prove the hazard is gone end to end: maximize NULL-valued ctx fields ---
# single trait, no training set, no posterior, no checks, no LD pruning, no
# lethal spec, write_outputs = TRUE -- the exact combination that exposed both
# live instances (trait_check_reference, ctc).
set.seed(404)
n_p <- 10L; n_m <- 30L
geno_mat <- matrix(rbinom(n_p * n_m, 1, 0.4) * 2L, nrow = n_p,
                    dimnames = list(paste0("P", seq_len(n_p)), paste0("m", seq_len(n_m))))
# NOTE (deviation from brief): ng_run_cp_canonical_id_table() always does
# as.data.frame(genotype) and then looks up id_col as a COLUMN -- a bare matrix with only
# rownames loses its IDs in that conversion and fails with "genotype is missing id_col: NAME".
# genotype must carry the id as an explicit column, matching every other runner test fixture
# (e.g. tests/check_reference_invariant.R).
geno <- data.frame(NAME = rownames(geno_mat), geno_mat, check.names = FALSE, stringsAsFactors = FALSE)
pheno <- data.frame(NAME = rownames(geno_mat), yield = rnorm(n_p, 10, 2), stringsAsFactors = FALSE)
map <- data.frame(marker = colnames(geno_mat), chr = 1L,
                  bp = seq_len(n_m) * 1e5, stringsAsFactors = FALSE)
dir_df <- data.frame(trait = "yield", column = "yield", direction = "increase",
                     stringsAsFactors = FALSE)
out_dir <- file.path(tempdir(), "ctx-null-run")

cat("running maximal-NULL end-to-end case ...\n")
res <- ng_run_cross_prediction(
  genotype = geno, phenotype = pheno, marker_map = map, trait_direction = dir_df,
  id_col = "NAME", bp_per_cm = 1e6, n_crosses = 3L,
  write_outputs = TRUE, write_figures = FALSE, output_dir = out_dir, seed = 7L)
stopifnot(is.list(res), nrow(res$candidate_crosses) > 0L)
stopifnot(length(list.files(out_dir, pattern = "[.]xlsx$")) >= 1L)

cat("ctx null-field end-to-end ok\n")
