# ng_run_cross_prediction() must let the user enlarge the marker-effect TRAINING set with
# extra individuals ("others") that are NOT candidate parents. The extra individuals augment
# the ridge fit only; they must be excluded from candidate crosses, allocation, and outputs.
# The real parents are exactly the main genotype/phenotype tables; training-only individuals
# come in via training_genotype / training_phenotype and never appear in selected_crosses.
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(4242)
n_par <- 12L; n_tr <- 48L; m <- 150L
markers <- sprintf("M%03d", seq_len(m))
par_ids  <- sprintf("PAR%02d", seq_len(n_par))
tr_ids   <- sprintf("OTH%02d", seq_len(n_tr))

# Parents are inbred lines (0/2); training-only individuals may be outbred (0/1/2).
freq <- runif(m, 0.2, 0.8)
G_par <- vapply(freq, function(p) 2L * rbinom(n_par, 1, p), integer(n_par))
G_tr  <- vapply(freq, function(p) rbinom(n_tr, 2, p), integer(n_tr))
colnames(G_par) <- colnames(G_tr) <- markers

# One heritable trait shared by parents and training individuals.
b <- rnorm(m) * (runif(m) < 0.15)
tv <- function(G) as.numeric(scale(G %*% b))
y_par <- 60 + 5 * (tv(G_par) + rnorm(n_par))
y_tr  <- 60 + 5 * (tv(G_tr)  + rnorm(n_tr))

genotype  <- data.frame(NAME = par_ids, G_par, check.names = FALSE)
phenotype <- data.frame(NAME = par_ids, yield = y_par)
train_geno  <- data.frame(NAME = tr_ids, G_tr, check.names = FALSE)
train_pheno <- data.frame(NAME = tr_ids, yield = y_tr)
mm  <- data.frame(SNP = markers, Chr = rep(1:6, length.out = m),
                  PosBP = rep(seq(0, 5e6, length.out = 25), 6)[seq_len(m)])
dir <- data.frame(trait = "yield", column = "yield", direction = "increase")

base_args <- list(
  phenotype = phenotype, genotype = genotype, marker_map = mm, trait_direction = dir,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "Chr", map_pos_col = "PosBP",
  map_pos_cm_divisor = 1e6, prediction_mode = "trait_by_trait",
  trait_value_metric = "var_complex", duplicate_action = "none",
  n_crosses = 10L, max_crosses_per_parent = 4L, optimizer = "greedy_local",
  assume_inbred = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 7L
)
run <- function(...) do.call(ng_run_cross_prediction, modifyList(base_args, list(...)))

## 1. With a training set, the run succeeds and NEVER crosses a training-only individual.
res <- run(training_genotype = train_geno, training_phenotype = train_pheno)
stopifnot(inherits(res, "ng_cross_prediction_result"))
sel_ids <- as.character(c(res$selected_crosses$parent1, res$selected_crosses$parent2))
stopifnot(all(sel_ids %in% par_ids))                 # only real parents are crossed
stopifnot(!any(tr_ids %in% sel_ids))                 # no training-only individual is crossed
cand_ids <- as.character(c(res$candidate_crosses$parent1, res$candidate_crosses$parent2))
stopifnot(all(cand_ids %in% par_ids))                # candidates are parents only

## 2. The result reports the training set explicitly (for the frontend to display).
aud <- res$input_match_audit
stopifnot(identical(as.integer(aud$training_only_count), n_tr))
stopifnot(identical(as.integer(aud$effect_training_n), n_par + n_tr))
stopifnot(setequal(aud$training_ids, tr_ids))
stopifnot(identical(as.integer(aud$matched_parent_count), n_par))   # parents unchanged
# per-trait effect summary records how many individuals trained the effects
stopifnot("marker_effect_training_n" %in% names(res$effect_summary))
stopifnot(res$effect_summary$marker_effect_training_n[[1L]] == n_par + n_tr)

## 3. A run WITHOUT the training set trains on parents only.
res0 <- run()
stopifnot(identical(as.integer(res0$input_match_audit$training_only_count), 0L))
stopifnot(identical(as.integer(res0$input_match_audit$effect_training_n), n_par))

## 4. Marker mismatch is a clear error, not a silent wrong fit.
bad_geno <- train_geno[, c("NAME", markers[-1])]        # drop one parent marker
err <- tryCatch(run(training_genotype = bad_geno, training_phenotype = train_pheno),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("marker", err, ignore.case = TRUE))

## 5. le does not use marker effects, so supplying training data warns (and is inert).
w <- tryCatch(
  withCallingHandlers(
    run(trait_value_metric = "le",
        training_genotype = train_geno, training_phenotype = train_pheno),
    warning = function(cnd) { message("caught: ", conditionMessage(cnd)); invokeRestart("muffleWarning") }
  ), error = function(e) e)
stopifnot(!inherits(w, "error"))                        # le + training still runs

cat("runner_training_set.R: PASS\n")
