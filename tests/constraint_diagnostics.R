# The run result must carry a structured `constraint_diagnostics` record so the frontend
# can surface visible run notes for what the breeder knobs actually did to the plan:
#   (2) breeder mating constraints  - plan shrink, min-unique relaxation, committed, min-use
#   (3) marker steering & lethal guarding - carrier x carrier drops, steering active/weight
#   (4) cost & logistics            - budget binding, cost/logistic emphasis
# These effects were previously only emitted to stderr (invisible to the out-of-process UI).
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11)
n <- 18L; m <- 90L
ids <- sprintf("P%02d", seq_len(n))
gm  <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
              dimnames = list(ids, sprintf("M%03d", seq_len(m))))
y   <- as.numeric(gm %*% rnorm(m, 0, 0.1)) + rnorm(n)
genotype  <- data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE)
runmm     <- data.frame(SNP_code = colnames(gm), Chromosome = rep(1:3, length.out = m),
                        Position_BP = rep(seq(0, 100, length.out = 30), 3)[seq_len(m)] * 1e6,
                        stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)
lethal <- ng_lethal_recessive_spec(marker = c("M001", "M002"), risk_allele = "alt")
mtar   <- ng_marker_target_spec(marker = c("M010", "M011"), direction = "increase", weight = 1)

run <- function(...) ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm,
  trait_direction = direction, id_col = "NAME",
  map_marker_col = "SNP_code", map_chr_col = "Chromosome",
  map_pos_col = "Position_BP", map_pos_cm_divisor = 1e6,
  trait_value_metric = "usefulness", n_crosses = 8L, max_crosses_per_parent = 3L,
  use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L, ...)

# --- 1. field present with lethal + marker steering active -----------------
res <- run(min_unique_parents = 10L, lethal_spec = lethal,
           marker_target_spec = mtar, lambda_marker = 0.5)
cd <- res$constraint_diagnostics
stopifnot(!is.null(cd), is.list(cd))
stopifnot(cd$allocation_method == "ocs")
stopifnot(cd$n_crosses_requested == 8L, is.finite(cd$n_crosses_delivered))

# (3) lethal guarding: active, 2 loci, dropped count = pre - post, all carriers removed
stopifnot(isTRUE(cd$lethal_active), cd$lethal_n_loci == 2L)
stopifnot(is.finite(cd$lethal_candidates_pre), cd$lethal_candidates_pre > 0L)
stopifnot(cd$lethal_dropped >= 0L, cd$lethal_dropped <= cd$lethal_candidates_pre)
stopifnot(isTRUE(cd$lethal_dropped_from_plan))
# no cross in the delivered plan may be a carrier x carrier mating
if (!is.null(res$candidate_crosses$lethal_carrier_cross))
  stopifnot(!any(isTRUE(res$candidate_crosses$lethal_carrier_cross)))

# (3) marker steering: active, blended into the criterion at the given weight
stopifnot(isTRUE(cd$marker_steering_active), cd$marker_n_loci == 2L)
stopifnot(isTRUE(all.equal(cd$marker_lambda, 0.5)), isTRUE(cd$marker_blended))

# --- 2. quiet run: no lethal/marker/cost knobs -> flags default off --------
res0 <- run()
cd0 <- res0$constraint_diagnostics
stopifnot(!is.null(cd0))
stopifnot(isFALSE(cd0$lethal_active), cd0$lethal_dropped == 0L)
stopifnot(isFALSE(cd0$marker_steering_active))
stopifnot(isFALSE(cd0$budget_active))
stopifnot(cd0$cost_emphasis == 0, cd0$logistic_emphasis == 0)

# --- 3. steering reported but NOT blended when weight is zero ---------------
resw0 <- run(marker_target_spec = mtar, lambda_marker = 0)
cdw0 <- resw0$constraint_diagnostics
stopifnot(isTRUE(cdw0$marker_steering_active), isFALSE(cdw0$marker_blended))

cat("constraint diagnostics test passed\n")
