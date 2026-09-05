ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11)
n_p <- 12L; n_m <- 40L
ids <- paste0("P", seq_len(n_p))
markers <- paste0("m", seq_len(n_m))
gm <- matrix(2L * rbinom(n_p * n_m, 1, 0.4), nrow = n_p, dimnames = list(ids, markers))
genotype <- data.frame(id = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(marker = markers, chr = rep(1:4, length.out = n_m),
                         bp = rep(seq(0, 100, length.out = 10), 4)[seq_len(n_m)] * 1e6,
                         stringsAsFactors = FALSE)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                        stringsAsFactors = FALSE)
pheno <- data.frame(id = ids, yield = rnorm(n_p, 10, 2), stringsAsFactors = FALSE)

# the check is NOT among the parents -- the case the old veto rejected outright
chk <- matrix(2L * rbinom(2 * n_m, 1, 0.4), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), markers))

args <- list(genotype = genotype, phenotype = pheno, marker_map = marker_map,
             map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
             bp_per_cm = 1e6, id_col = "id",
             trait_direction = direction, n_crosses = 5L,
             write_outputs = FALSE, write_figures = FALSE, seed = 5L)

# The runner's default ridge fit never stamps a calibrated reliability (R/02_effects.R sets
# reliability_is_calibrated = FALSE), so ng_choose_mean_source() always resolves to the
# phenotypic "adjusted_pheno" source here (it is passed the raw phenotype unconditionally) --
# never a GEBV variant. check_records is exactly the mechanism the interface names for this
# case ("only consulted when the run's source for that trait is not a GEBV variant").
chk_records <- list(yield = list(adjusted_pheno = c(CHK_A = 11.5, CHK_B = 9.25)))

res <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_records = chk_records,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))

ref <- res$trait_check_reference
stopifnot(is.list(ref), nrow(ref$active) == 1L, ref$active$check == "CHK_A")
stopifnot(is.finite(ref$values$yield[["CHK_A"]]))
stopifnot(abs(ref$values$yield[["CHK_A"]] - 11.5) < 1e-8)
stopifnot(identical(ref$source[["yield"]], "adjusted_pheno"))
stopifnot(nzchar(ref$source[["yield"]]))

ct <- res$candidate_crosses
stopifnot(all(c("yield_check_value", "yield_vs_check", "yield_check_ok",
                "yield_p_beat_check", "checks_all_ok") %in% names(ct)))
# the check is never a parent of any candidate cross
stopifnot(!any(c(ct$parent1, ct$parent2) %in% c("CHK_A", "CHK_B")))

# a check id colliding with a parent id is rejected
err <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = matrix(0, nrow = 1, ncol = n_m,
                      dimnames = list("P1", markers)),
  check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "P1", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("also a candidate parent", err))

# a check named in trait_checks but absent from check_geno is rejected
err2 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "NOPE", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err2), grepl("NOPE", err2))

# check_progeny_size is REQUIRED with trait_checks -- never defaulted
err3 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err3), grepl("check_progeny_size", err3))
# and it must be a sensible count
err4 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 0L,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err4), grepl("check_progeny_size", err4))

# the size the user gave is what drives the probability, and is stamped on the result
stopifnot(identical(ref$progeny_size, 200L))

cat("runner wiring ok\n")

# --- Task 10: check_pheno -> real check values on a phenotype-source run ----
# check_pheno makes the check values REAL rather than NA on a phenotype-source run.
# check_pheno's id column must be the run's own id_col ("id" in `args`, above) --
# it is read the same way the phenotype table itself is, not by a separate id keyword.
cp <- data.frame(id = c("CHK_A", "CHK_B"), yield = c(12.5, 8.0), stringsAsFactors = FALSE)
res_cp <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(identical(res_cp$trait_check_reference$source[["yield"]], "adjusted_pheno"))
stopifnot(abs(res_cp$trait_check_reference$values$yield[["CHK_A"]] - 12.5) < 1e-8)
ct_cp <- res_cp$candidate_crosses
stopifnot(all(is.finite(ct_cp$yield_check_value)))
stopifnot(any(is.finite(ct_cp$yield_p_beat_check)))

# precedence: an explicit check_records value WINS over check_pheno for the same trait, even
# when both are supplied for it
res_prec <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L, check_pheno = cp,
  check_records = list(yield = list(adjusted_pheno = c(CHK_A = 99))),
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))
stopifnot(abs(res_prec$trait_check_reference$values$yield[["CHK_A"]] - 99) < 1e-8)

cat("task 10 runner ok\n")
