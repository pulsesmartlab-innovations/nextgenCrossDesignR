ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---------------------------------------------------------------------------
# ng_check_line_value() places a check on the axis that ACTUALLY exists, and that axis depends
# on the run's real scoring method family (multi_trait_meta$family):
#   - "rank_threshold" (auto) / "threshold" / "weighted_index" (weighted): the score is
#     z %*% weights where z = ng_rank_normalize() -- a RANK transform. The check is placed by
#     quantile lookup among the candidates' own oriented raw values (ng_check_rank_axis_z()),
#     which needs the candidate score table itself, not just metadata.
#   - "economic_index" / "desired_gain": the score is value_z %*% <SOLVED coefficients>. The
#     check maps through the retained affine value_z transform and combines with the solved
#     coefficients (not the raw input weights, which differ from the solved vector in general).
#   - anything else, or missing parameters: NA_real_ -- never a plausible-looking substitute.
# ---------------------------------------------------------------------------

# --- ng_check_rank_axis_z(): unit-test the quantile-placement primitive in isolation ---
# 7 candidates, increasing trait. Check value 5.5 sits strictly between candidates 5 and 6 of 7.
cand_x <- c(1, 2, 3, 4, 5, 6, 7)
z_check <- ng_check_rank_axis_z(cand_x, bigger_is_better = TRUE, tau_raw = 5.5)
# Independent hand computation (mirrors ng_rank_normalize(), not the package's own helper):
r <- rank(cand_x, ties.method = "average")
p <- (r - 0.5) / length(r)
out <- stats::qnorm(p)
mu <- mean(out); sdv <- stats::sd(out)
p_check_manual <- (sum(cand_x <= 5.5) + 0.5) / length(cand_x)  # 5 of 7 candidates are <= 5.5
z_check_manual <- (stats::qnorm(p_check_manual) - mu) / sdv
stopifnot(abs(z_check - z_check_manual) < 1e-10)
# A decreasing trait orients by -x throughout (candidates AND the check), mirroring
# ng_rank_normalize()'s bigger_is_better = FALSE. NOT a simple sign flip of the increasing-trait
# result: the "+0.5" continuity correction in p_check is asymmetric under negation whenever the
# check does not exactly tie a candidate (5 of 7 candidates are <= 5.5, but only 2 are >= 5.5,
# not the complementary 7-5=2 by coincidence here, so this is worth computing independently
# rather than assumed).
value_dec <- -cand_x
r_dec <- rank(value_dec, ties.method = "average")
p_dec <- (r_dec - 0.5) / length(r_dec)
out_dec <- stats::qnorm(p_dec)
mu_dec <- mean(out_dec); sd_dec <- stats::sd(out_dec)
p_check_dec_manual <- (sum(value_dec <= -5.5) + 0.5) / length(value_dec)
z_check_dec_manual <- (stats::qnorm(p_check_dec_manual) - mu_dec) / sd_dec
z_check_dec <- ng_check_rank_axis_z(cand_x, bigger_is_better = FALSE, tau_raw = 5.5)
stopifnot(abs(z_check_dec - z_check_dec_manual) < 1e-10)
# Degenerate inputs -> NA, never a fabricated number.
stopifnot(is.na(ng_check_rank_axis_z(c(1, NA), bigger_is_better = TRUE, tau_raw = 1)))
stopifnot(is.na(ng_check_rank_axis_z(cand_x, bigger_is_better = TRUE, tau_raw = NA_real_)))

# --- ng_check_line_value(): LINEAR families (economic_index / desired_gain) ---
# value_z's center/scale are computed from `_value` (mean +/- i*sigma under the default
# trait_value_metric = "usefulness"), which is NOT an affine function of the mean in general, so
# remapping the check's (mean-scale) tau through value_z is only exact when trait_value_metric ==
# "mean". Every call below therefore passes trait_value_metric = "mean" explicitly to exercise
# the linear-family math; the refusal-otherwise behaviour has its own tests further down.
ref <- list(active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                                stringsAsFactors = FALSE),
            values = list(yield = c(CHK_A = 6)), source = c(yield = "GEBV"))
meta_econ_single <- list(family = "economic_index",
                        economic_index_coefficients = c(yield = 1),
                        value_centers = c(yield = 5), value_scales = c(yield = 2),
                        value_signs = c(yield = 1))
stopifnot(abs(ng_check_line_value(ref, meta_econ_single, trait = "yield",
                                  trait_value_metric = "mean") - 0.5) < 1e-8)

ref2 <- list(active = data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_A"),
                                 reject_if = "below", stringsAsFactors = FALSE),
             values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
             source = c(yield = "GEBV", protein = "GEBV"))
# SOLVED coefficients differ from the raw/normalized weights on purpose, so this exercises that
# the function uses the coefficients, not $weights.
meta_econ_multi <- list(family = "economic_index",
                        economic_index_coefficients = c(yield = 0.7, protein = 0.3),
                        weights = c(yield = 0.5, protein = 0.5),
                        value_centers = c(yield = 0, protein = 0),
                        value_scales = c(yield = 1, protein = 1),
                        value_signs = c(yield = 1, protein = 1))
stopifnot(abs(ng_check_line_value(ref2, meta_econ_multi, trait_value_metric = "mean") -
             (6 * 0.7 + 10 * 0.3)) < 1e-8)

meta_desired <- list(family = "desired_gain",
                     desired_gain_coefficients = c(yield = 0.4, protein = 0.6),
                     value_centers = c(yield = 0, protein = 0),
                     value_scales = c(yield = 1, protein = 1),
                     value_signs = c(yield = 1, protein = 1))
stopifnot(abs(ng_check_line_value(ref2, meta_desired, trait_value_metric = "mean") -
             (6 * 0.4 + 10 * 0.6)) < 1e-8)

# Linear family but the solved coefficients are missing from the metadata -> NA.
meta_econ_no_coef <- list(family = "economic_index",
                          value_centers = c(yield = 0), value_scales = c(yield = 1),
                          value_signs = c(yield = 1))
stopifnot(is.na(ng_check_line_value(ref, meta_econ_no_coef, trait = "yield",
                                    trait_value_metric = "mean")))

# Linear family, coefficients present, BUT trait_value_metric is not "mean" (the default,
# "usefulness", included) -> NA. value_z is not an affine function of the mean under any other
# metric, so there is no exact remap and the function must refuse rather than draw a
# plausible-but-wrong line.
stopifnot(is.na(ng_check_line_value(ref, meta_econ_single, trait = "yield")))                       # metric omitted
stopifnot(is.na(ng_check_line_value(ref, meta_econ_single, trait = "yield",
                                    trait_value_metric = "usefulness")))
stopifnot(is.na(ng_check_line_value(ref, meta_econ_single, trait = "yield",
                                    trait_value_metric = "pmv")))

# Linear family, trait_value_metric = "mean", coefficients cover a 3-trait index, but only ONE
# trait is checked -> NA. Silently summing only the checked term is mathematically identical to
# imputing the other two traits' value_z at exactly their population median (value_z's own
# center), i.e. a plausible-looking substitute for a position the check does not actually have.
meta_econ_3trait <- list(family = "economic_index",
                        economic_index_coefficients = c(yield = 0.5, protein = 0.3, disease = 0.2),
                        value_centers = c(yield = 0, protein = 0, disease = 0),
                        value_scales = c(yield = 1, protein = 1, disease = 1),
                        value_signs = c(yield = 1, protein = 1, disease = 1))
stopifnot(is.na(ng_check_line_value(ref, meta_econ_3trait, trait = "yield",
                                    trait_value_metric = "mean")))
# ...but a check on EVERY index trait is fine (no traits omitted from the sum).
ref_3chk <- list(active = data.frame(trait = c("yield", "protein", "disease"),
                                     check = c("CHK_A", "CHK_A", "CHK_A"),
                                     reject_if = "below", stringsAsFactors = FALSE),
                 values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10), disease = c(CHK_A = 2)),
                 source = c(yield = "GEBV", protein = "GEBV", disease = "GEBV"))
stopifnot(abs(ng_check_line_value(ref_3chk, meta_econ_3trait, trait_value_metric = "mean") -
             (6 * 0.5 + 10 * 0.3 + 2 * 0.2)) < 1e-8)

# --- ng_check_line_value(): RANK families (auto / threshold / weighted) ---
# Small synthetic candidate table: one trait "yield". The check's tau is on the MEAN scale
# (ng_check_reference_value() resolves it onto whatever source produced the cross MEANS), so the
# rank placement must compare it against the candidates' "_mean" column, never "_value"
# (mean +/- i*sigma under the default metric) -- see the REAL PIPELINE section below for the
# worked example of what ranking against the wrong column actually does to the result.
rank_scores <- data.frame(yield_mean = c(8, 9, 10, 11, 12, 20, 21))
rank_traits <- data.frame(trait = "yield", column = "yield_value", direction = "maximize",
                          stringsAsFactors = FALSE)
meta_rank <- list(family = "rank_threshold", weights = c(yield = 1),
                  value_signs = c(yield = 1), traits = rank_traits)
cl_rank <- ng_check_line_value(ref, meta_rank, candidate_scores = rank_scores)
manual_rank <- ng_check_rank_axis_z(rank_scores$yield_mean, bigger_is_better = TRUE, tau_raw = 6)
stopifnot(abs(cl_rank - manual_rank) < 1e-10)

# Rank family but no candidate_scores supplied -> NA (this is exactly what a caller building
# multi_trait_meta by hand, or an older stored result, would hit).
stopifnot(is.na(ng_check_line_value(ref, meta_rank)))

# "threshold" and "weighted_index" families take the identical rank path.
meta_rank2 <- meta_rank; meta_rank2$family <- "threshold"
stopifnot(abs(ng_check_line_value(ref, meta_rank2, candidate_scores = rank_scores) - manual_rank) < 1e-10)
meta_rank3 <- meta_rank; meta_rank3$family <- "weighted_index"
stopifnot(abs(ng_check_line_value(ref, meta_rank3, candidate_scores = rank_scores) - manual_rank) < 1e-10)

# Any other family, or a meta with no family/method at all -> NA.
stopifnot(is.na(ng_check_line_value(ref, list(family = "other"), trait = "yield")))
stopifnot(is.na(ng_check_line_value(ref, NULL, trait = "yield")))
stopifnot(is.na(ng_check_line_value(ref, list())))

# an unevaluable check yields no line rather than a misplaced one
ref3 <- ref; ref3$values$yield <- c(CHK_A = NA_real_)
stopifnot(is.na(ng_check_line_value(ref3, meta_econ_single, trait = "yield",
                                    trait_value_metric = "mean")))
stopifnot(is.na(ng_check_line_value(ref3, meta_rank, candidate_scores = rank_scores)))

# ---------------------------------------------------------------------------
# ng_check_line_label(): names the check id, or names the blend when the traits contributing to
# the line use DIFFERENT checks (naming only the first check would misrepresent a blend as a
# single check's value).
# ---------------------------------------------------------------------------
stopifnot(identical(ng_check_line_label(ref$active), "CHK_A"))
stopifnot(identical(ng_check_line_label(ref2$active), "CHK_A"))
ref_two_checks <- data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_B"),
                             reject_if = "below", stringsAsFactors = FALSE)
lbl <- ng_check_line_label(ref_two_checks)
stopifnot(grepl("CHK_A", lbl), grepl("CHK_B", lbl), !identical(lbl, "CHK_A"))

# ---------------------------------------------------------------------------
# Plot rendering: a byte-identical PNG passed a file-existence check and hid the original
# off-axis bug (with-checks and no-checks PNGs came back byte-identical because the line was
# off the plotted range). Every rendering test below therefore asserts a BYTE DIFFERENCE against
# the same plot drawn with no line, not merely that a file was written.
# ---------------------------------------------------------------------------
scored <- data.frame(parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
                     multi_trait_score = c(9, 4), pair_kinship = c(0.1, 0.2),
                     priority_tier = c("priority", "low_priority"), stringsAsFactors = FALSE)
p_base <- file.path(tempdir(), "chk-plot-base.png")
ng_plot_priority_score_vs_kinship(scored, output_path = p_base)
p_line <- file.path(tempdir(), "chk-plot-line.png")
ng_plot_priority_score_vs_kinship(scored, output_path = p_line, check_line = 6,
                                  check_label = "CHK_A")
stopifnot(file.exists(p_line), file.size(p_line) > 0)
stopifnot(!identical(readBin(p_base, "raw", file.info(p_base)$size),
                     readBin(p_line, "raw", file.info(p_line)$size)))

# --- ng_plot_check_reference_line(): both orientations, byte-diff against a line-free baseline ---
draw_baseline <- function(path) {
  grDevices::png(path, width = 400, height = 300)
  on.exit(grDevices::dev.off())
  graphics::plot(1:10, 1:10, main = "orientation test")
}
draw_with_line <- function(path, mean_axis) {
  grDevices::png(path, width = 400, height = 300)
  on.exit(grDevices::dev.off())
  graphics::plot(1:10, 1:10, main = "orientation test")
  ng_plot_check_reference_line(5.5, label = "CHK_A", mean_axis = mean_axis)
}
p_orient_base <- file.path(tempdir(), "chk-orient-base.png")
p_orient_y <- file.path(tempdir(), "chk-orient-y.png")
p_orient_x <- file.path(tempdir(), "chk-orient-x.png")
draw_baseline(p_orient_base)
draw_with_line(p_orient_y, "y")
draw_with_line(p_orient_x, "x")
bytes_base <- readBin(p_orient_base, "raw", file.info(p_orient_base)$size)
bytes_y <- readBin(p_orient_y, "raw", file.info(p_orient_y)$size)
bytes_x <- readBin(p_orient_x, "raw", file.info(p_orient_x)$size)
stopifnot(!identical(bytes_base, bytes_y))   # horizontal line changes the plot
stopifnot(!identical(bytes_base, bytes_x))   # vertical line changes the plot
stopifnot(!identical(bytes_y, bytes_x))      # the two orientations are not the same drawing

# mean_axis = NULL (the default) draws nothing -- a plot with no mean-bearing axis gets no line,
# never a guessed default orientation.
p_orient_none <- file.path(tempdir(), "chk-orient-none.png")
grDevices::png(p_orient_none, width = 400, height = 300)
graphics::plot(1:10, 1:10, main = "orientation test")
ng_plot_check_reference_line(5.5, label = "CHK_A", mean_axis = NULL)
invisible(grDevices::dev.off())
stopifnot(identical(readBin(p_orient_none, "raw", file.info(p_orient_none)$size), bytes_base))

# multi-trait panels: one panel per trait with a check, each with its own line
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
scored_mt_no_check <- scored_mt
scored_mt_no_check$yield_check_value <- NA_real_
scored_mt_no_check$protein_check_value <- NA_real_
p2_base <- file.path(tempdir(), "chk-panels-base.png")
ng_plot_check_panels(scored_mt_no_check, ref_mt, output_path = p2_base)
stopifnot(!identical(readBin(p2, "raw", file.info(p2)$size),
                     readBin(p2_base, "raw", file.info(p2_base)$size)))

# a run with no checks draws nothing and returns NULL rather than erroring
stopifnot(is.null(ng_plot_check_panels(scored_mt, NULL, output_path = p2)))

# ---------------------------------------------------------------------------
# REAL PIPELINE: the resolved line on a DEFAULT ("auto") single-trait run must agree with an
# INDEPENDENTLY computed quantile placement to a tight tolerance -- not merely fall "in range".
# An earlier version of this test only checked range membership, which a genuinely wrong value_z
# -based placement (-0.327) satisfied while still disagreeing with the true rank placement
# (+0.688) by about 20% of the axis span. That gap is exactly what this test now catches.
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
stopifnot(identical(rp_mt$family, "rank_threshold"))  # confirms this hits the rank path, not value_z
rp_cl <- ng_check_line_value(rp_res$trait_check_reference, rp_mt,
                             candidate_scores = rp_res$candidate_crosses)
stopifnot(is.finite(rp_cl))

# Independent recomputation -- hand-rolled rank/qnorm arithmetic, NOT a call to the package's
# ng_check_rank_axis_z() helper, so this genuinely cross-checks the implementation rather than
# restating it. Ranks against yield_MEAN, not yield_value: the check's tau (11.5) is on the mean
# scale (ng_check_reference_value() resolves it onto whatever produced the cross means), while
# yield_value is mean + i*sigma under the run's default trait_value_metric = "usefulness" -- a
# DIFFERENT quantity. Ranking against yield_value was the C1 bug: it silently compared the check
# to the wrong column, which is exactly what this independent recomputation would have caught had
# it not itself restated the same wrong column (this is the fix to that self-check).
rp_x <- rp_res$candidate_crosses$yield_mean
rp_r <- rank(rp_x, ties.method = "average")
rp_p <- (rp_r - 0.5) / length(rp_r)
rp_out <- stats::qnorm(rp_p)
rp_mu <- mean(rp_out); rp_sd <- stats::sd(rp_out)
rp_p_check <- (sum(rp_x <= 11.5) + 0.5) / length(rp_x)
rp_indep <- (stats::qnorm(rp_p_check) - rp_mu) / rp_sd
stopifnot(abs(rp_cl - rp_indep) < 1e-8)

rp_rng <- range(rp_res$candidate_crosses$multi_trait_score, na.rm = TRUE)
stopifnot(rp_cl >= rp_rng[1] - diff(rp_rng), rp_cl <= rp_rng[2] + diff(rp_rng))

# The check line's position must agree with the Checks-sheet/column truth (yield_check_ok,
# R/51_check_reference.R), which compares the SAME yield_mean column to the SAME tau. rp_cl > 0
# means the check sits ABOVE the median of the candidates' rank-normal axis (few candidates
# reach that high), which must correspond to a MINORITY of candidates clearing the check.
frac_ok <- mean(rp_res$candidate_crosses$yield_check_ok)
stopifnot(identical(rp_cl > 0, frac_ok < 0.5))

cat("check plot ok\n")

# --- panel labels stay aligned when a checked trait is dropped ---------------
# ng_plot_check_panels() filters traits to those whose <key>_mean column exists in `scored`.
# traits, key and the check ids must be filtered together: indexing the UNFILTERED spec inside
# the loop mislabels every panel after a dropped trait. Unreachable through the runner (the
# attach step errors first on a missing mean column) but ng_plot_check_panels is exported, so a
# direct caller with a partial `scored` hits it.
spec_drop <- data.frame(trait = c("yield", "protein", "disease"),
                        check = c("CHK_Y", "CHK_P", "CHK_D"),
                        reject_if = "below",
                        column_key = c("yield", "protein", "disease"),
                        stringsAsFactors = FALSE)
scored_drop <- data.frame(
  pair_kinship = c(0.1, 0.2),
  yield_mean = c(9, 4), yield_check_value = 6, yield_check_ok = c(TRUE, FALSE),
  disease_mean = c(3, 5), disease_check_value = 4, disease_check_ok = c(TRUE, FALSE),
  stringsAsFactors = FALSE)                       # protein_mean deliberately absent

drop_seen <- character(0)
drop_env <- environment(ng_plot_check_panels)
drop_orig <- get("ng_plot_check_reference_line", envir = drop_env)
assign("ng_plot_check_reference_line",
       function(value, label = NULL, mean_axis = NULL, ...) {
         drop_seen <<- c(drop_seen, if (is.null(label)) NA_character_ else label)
         invisible(NULL)
       }, envir = drop_env)
ng_plot_check_panels(scored_drop, list(active = spec_drop),
                     output_path = file.path(tempdir(), "chk-drop.png"))
assign("ng_plot_check_reference_line", drop_orig, envir = drop_env)
# the middle trait is dropped, so the surviving panels must be labelled with THEIR OWN checks
stopifnot(identical(drop_seen, c("CHK_Y", "CHK_D")))

cat("check panel labels stay aligned across dropped traits\n")
