ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---------------------------------------------------------------------------
# ng_check_line_value(): the plotted axis (multi_trait_score) is ALWAYS z-normalised by
# ng_add_multitrait_score() (R/19), even for a single trait. A check has no row in that
# computation, so it has no z of its own; it must be pushed through the SAME affine transform
# (oriented = sign * tau; z = (oriented - center) / scale) via the value_centers/value_scales/
# value_signs the runner retains in multi_trait_meta. Without that transform there is no honest
# line -- same discipline as the rank-based-index case, applied uniformly.
# ---------------------------------------------------------------------------

ref <- list(active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                                stringsAsFactors = FALSE),
            values = list(yield = c(CHK_A = 6)), source = c(yield = "GEBV"))

# SINGLE TRAIT, transform present (identity: center 0, scale 1, sign +1) -- the check maps
# straight through, so the line equals the raw check value in this degenerate case.
meta_single_identity <- list(value_centers = c(yield = 0), value_scales = c(yield = 1),
                             value_signs = c(yield = 1))
stopifnot(ng_check_line_value(ref, meta_single_identity, trait = "yield") == 6)

# SINGLE TRAIT, transform present and NON-trivial -- the check must land on the SAME
# standardized axis as the candidates, not in raw trait units.
meta_single_nontrivial <- list(value_centers = c(yield = 5), value_scales = c(yield = 2),
                               value_signs = c(yield = 1))
stopifnot(abs(ng_check_line_value(ref, meta_single_nontrivial, trait = "yield") - 0.5) < 1e-8)

# SINGLE TRAIT, NO transform available (e.g. multi_trait_meta = NULL, or a hand-built meta that
# never carried it) -- there is no honest way to place the check on a z-normalised axis, so this
# must be NA, not the raw check value. This is exactly the case that was silently wrong before:
# a real single-trait run's check line was drawn in raw units on a z-normalised axis and fell
# off the plotted range entirely.
stopifnot(is.na(ng_check_line_value(ref, multi_trait_meta = NULL, trait = "yield")))
stopifnot(is.na(ng_check_line_value(ref, list(method = "weighted"), trait = "yield")))

# LINEAR INDEX: apply the same coefficients to the check's per-trait Z values (identity
# transform here isolates the weighting logic from the standardization logic tested above)
ref2 <- list(active = data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_A"),
                                 reject_if = "below", stringsAsFactors = FALSE),
             values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
             source = c(yield = "GEBV", protein = "GEBV"))
meta_lin <- list(method = "weighted", weights = c(yield = 0.5, protein = 0.5),
                 value_centers = c(yield = 0, protein = 0),
                 value_scales = c(yield = 1, protein = 1),
                 value_signs = c(yield = 1, protein = 1))
stopifnot(abs(ng_check_line_value(ref2, meta_lin) - 8) < 1e-8)

# LINEAR INDEX, transform missing for one trait -- no honest combined line
meta_lin_partial <- meta_lin
meta_lin_partial$value_scales <- meta_lin_partial$value_scales["yield"]
stopifnot(is.na(ng_check_line_value(ref2, meta_lin_partial)))

# RANK-BASED INDEX: a check has no rank, so there is NO valid line, transform or not
meta_rank <- list(method = "rank_threshold", value_centers = c(yield = 0, protein = 0),
                  value_scales = c(yield = 1, protein = 1), value_signs = c(yield = 1, protein = 1))
stopifnot(is.na(ng_check_line_value(ref2, meta_rank)))

# an unevaluable check yields no line rather than a misplaced one
ref3 <- ref; ref3$values$yield <- c(CHK_A = NA_real_)
stopifnot(is.na(ng_check_line_value(ref3, meta_single_identity, trait = "yield")))

# ---------------------------------------------------------------------------
# ng_check_line_label(): names the check id, or names the blend when the traits contributing to
# the line use DIFFERENT checks (Finding 3 -- naming only the first check misrepresents a blend
# as a single check's value).
# ---------------------------------------------------------------------------
stopifnot(identical(ng_check_line_label(ref$active), "CHK_A"))
stopifnot(identical(ng_check_line_label(ref2$active), "CHK_A"))
ref_two_checks <- data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_B"),
                             reject_if = "below", stringsAsFactors = FALSE)
lbl <- ng_check_line_label(ref_two_checks)
stopifnot(grepl("CHK_A", lbl), grepl("CHK_B", lbl), !identical(lbl, "CHK_A"))

# ---------------------------------------------------------------------------
# Plot rendering: check_line/check_label are passed straight through as numeric/string, so
# these tests exercise the DRAWING code only (not ng_check_line_value's resolution logic).
# ---------------------------------------------------------------------------
scored <- data.frame(parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
                     multi_trait_score = c(9, 4), pair_kinship = c(0.1, 0.2),
                     priority_tier = c("priority", "low_priority"), stringsAsFactors = FALSE)
p <- file.path(tempdir(), "chk-plot.png")
ng_plot_priority_score_vs_kinship(scored, output_path = p, check_line = 6,
                                  check_label = "CHK_A")
stopifnot(file.exists(p), file.size(p) > 0)

# multi-trait: one panel per trait with a check, each with its own line
scored_mt <- data.frame(
  parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
  pair_kinship = c(0.1, 0.2),
  yield_mean = c(9, 4), yield_check_value = 6, yield_check_ok = c(TRUE, FALSE),
  protein_mean = c(11, 13), protein_check_value = 10, protein_check_ok = c(TRUE, TRUE),
  stringsAsFactors = FALSE)
ref_mt <- list(active = data.frame(trait = c("yield", "protein"), check = "CHK_A",
                                   reject_if = "below", stringsAsFactors = FALSE),
               values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
               source = c(yield = "GEBV", protein = "GEBV"))
p2 <- file.path(tempdir(), "chk-panels.png")
ng_plot_check_panels(scored_mt, ref_mt, output_path = p2)
stopifnot(file.exists(p2), file.size(p2) > 0)

# a run with no checks draws nothing and returns NULL rather than erroring
stopifnot(is.null(ng_plot_check_panels(scored_mt, NULL, output_path = p2)))

# ---------------------------------------------------------------------------
# REAL PIPELINE (Finding 2): the assertion that would actually have caught the axis-units bug.
# Every other check-reference test uses synthetic multi_trait_score values or write_figures =
# FALSE, so none of them can see whether the resolved line actually falls on the axis that gets
# plotted. This runs the FULL runner (single trait, default "auto" method -- the most common
# real call pattern) with a trait check attached, and requires the resolved line to fall inside
# the ACTUAL plotted data range (padded by one range-width, so this is a "not wildly off-axis"
# check, not a tight numerical equality). Confirmed RED against the pre-fix code (see task-9
# fix-round-1 report): the unfixed ng_check_line_value returned the check's raw 11.5 while the
# real run's multi_trait_score ranged -2.43..2.43 -- 11.5 is far outside even the padded range,
# and the assertion below failed for exactly that reason before the fix.
# ---------------------------------------------------------------------------
set.seed(11)
n_p <- 12L; n_m <- 40L
ids <- paste0("P", seq_len(n_p))
markers <- paste0("m", seq_len(n_m))
gm <- matrix(2L * rbinom(n_p * n_m, 1, 0.4), nrow = n_p, dimnames = list(ids, markers))
genotype <- data.frame(id = ids, gm, check.names = FALSE, stringsAsFactors = FALSE)
marker_map <- data.frame(marker = markers, chr = rep(1:4, length.out = n_m),
                         bp = rep(seq(0, 100, length.out = 10), 4)[seq_len(n_m)] * 1e6,
                         stringsAsFactors = FALSE)
rp_direction <- data.frame(trait = "yield", column = "yield", direction = "increase",
                          stringsAsFactors = FALSE)
rp_pheno <- data.frame(id = ids, yield = rnorm(n_p, 10, 2), stringsAsFactors = FALSE)
rp_chk <- matrix(2L * rbinom(2 * n_m, 1, 0.4), nrow = 2,
                 dimnames = list(c("CHK_A", "CHK_B"), markers))
rp_chk_records <- list(yield = list(adjusted_pheno = c(CHK_A = 11.5, CHK_B = 9.25)))
rp_res <- ng_run_cross_prediction(
  genotype = genotype, phenotype = rp_pheno, marker_map = marker_map,
  map_marker_col = "marker", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, id_col = "id",
  trait_direction = rp_direction, n_crosses = 5L,
  write_outputs = FALSE, write_figures = FALSE, seed = 5L,
  check_geno = rp_chk, check_progeny_size = 200L, check_records = rp_chk_records,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))

rp_mt <- attr(rp_res$candidate_crosses, "multi_trait")
rp_cl <- ng_check_line_value(rp_res$trait_check_reference, rp_mt)
stopifnot(is.finite(rp_cl))
rp_rng <- range(rp_res$candidate_crosses$multi_trait_score, na.rm = TRUE)
stopifnot(rp_cl >= rp_rng[1] - diff(rp_rng), rp_cl <= rp_rng[2] + diff(rp_rng))

cat("check plot ok\n")
