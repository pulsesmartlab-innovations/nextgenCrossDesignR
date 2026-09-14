# A run without a check line must not show check columns, a check sheet, or a check
# reference line on any plot.
#
# Every check artefact is conditional on the breeder having supplied one, and every
# such condition is a place where a NULL can leak through as a column full of NA, an
# empty sheet, or -- worst -- a dashed red line drawn at a value that came from
# nowhere. A reader cannot tell "no check was supplied" from "the check was not met"
# by looking at a blank cell, and a line on a plot is read as a threshold whether or
# not one exists.
#
# The test asserts BOTH directions. An absence-only test passes just as well when the
# check feature is broken outright, so the same fixture is run with a check to confirm
# the artefacts do appear when they should. Absence is only meaningful against a
# demonstrated presence.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(77L)
n <- 60L; m <- 40L; ids <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * m, 1, 0.5), n, m,
             dimnames = list(ids, sprintf("M%03d", seq_len(m))))
gv <- as.numeric(gm %*% c(rnorm(8L, 0, 1), rep(0, m - 8L)))
y  <- 10 + 2 * as.numeric(scale(gv + rnorm(n, 0, 0.2 * stats::sd(gv))))
chk <- matrix(2L * rbinom(2 * m, 1, 0.5), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), colnames(gm)))

base_args <- list(
  phenotype = data.frame(NAME = ids, yield = y, stringsAsFactors = FALSE),
  genotype  = data.frame(NAME = ids, gm, check.names = FALSE, stringsAsFactors = FALSE),
  marker_map = data.frame(SNP = colnames(gm), chr = rep(1:2, length.out = m),
                          bp = rep(seq(0, 100, length.out = m / 2), 2)[seq_len(m)] * 1e6),
  trait_direction = data.frame(trait = "yield", column = "yield", direction = "increase"),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  bp_per_cm = 1e6, n_crosses = 8L, write_outputs = TRUE, write_figures = TRUE,
  run_posterior_prediction = FALSE, seed = 5L)

CHECK_PAT <- "_check_id$|_check_value$|_vs_check$|_check_ok$|_p_beat_check$|checks_all_ok|p_beat_all_checks"

d_no <- tempfile("nochk_"); dir.create(d_no)
res_no <- do.call(ng_run_cross_prediction, c(base_args, list(output_dir = d_no)))

d_yes <- tempfile("chk_"); dir.create(d_yes)
res_yes <- do.call(ng_run_cross_prediction, c(base_args, list(
  output_dir = d_yes, check_geno = chk, check_progeny_size = 200L,
  check_records = list(yield = list(adjusted_pheno = c(CHK_A = 11.5, CHK_B = 9.25))),
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))

# ---- 1. the with-check run really does produce the artefacts --------------
# Establishes that the absences below mean something.
yes_cols <- grep(CHECK_PAT, names(res_yes$candidate_crosses), value = TRUE)
stopifnot(length(yes_cols) > 0L)
stopifnot(!is.null(res_yes$trait_check_reference))

# ---- 2. ...and the no-check run produces none of them --------------------
for (tbl in c("candidate_crosses", "selected_crosses")) {
  leaked <- grep(CHECK_PAT, names(res_no[[tbl]]), value = TRUE)
  if (length(leaked)) {
    stop("no check was supplied but ", tbl, " carries check column(s): ",
         paste(leaked, collapse = ", "), call. = FALSE)
  }
}
stopifnot(is.null(res_no$trait_check_reference))

# ---- 3. no check sheet, and no check column inside any sheet -------------
sheets_of <- function(dir) {
  wb <- file.path(dir, "crossing_plan.xlsx")
  if (!file.exists(wb)) return(NULL)
  openxlsx::getSheetNames(wb)
}
sh_yes <- sheets_of(d_yes); sh_no <- sheets_of(d_no)
if (!is.null(sh_yes)) {
  stopifnot(any(grepl("check", sh_yes, ignore.case = TRUE)))   # presence, first
  stopifnot(!any(grepl("check", sh_no, ignore.case = TRUE)))
  for (s in sh_no) {
    hdr <- tryCatch(names(openxlsx::read.xlsx(file.path(d_no, "crossing_plan.xlsx"),
                                              sheet = s, startRow = 2L)),
                    error = function(e) character())
    leaked <- grep(CHECK_PAT, hdr, value = TRUE)
    if (length(leaked)) {
      stop("no check was supplied but workbook sheet '", s, "' carries: ",
           paste(leaked, collapse = ", "), call. = FALSE)
    }
  }
}

# ---- 4. no check reference line is placed -------------------------------
mt_no <- attr(res_no$candidate_crosses, "multi_trait")
cl_no <- ng_check_line_value(res_no$trait_check_reference, mt_no,
                             candidate_scores = res_no$candidate_crosses,
                             trait_value_metric = "mean")
stopifnot(is.na(cl_no))

# ...and nothing is DRAWN. With no graphics device open, abline() errors, so a
# reference line that tried to draw would throw rather than silently appear. This
# is the observable: not "returns NULL", but "touched no device".
grDevices::graphics.off()
drew <- tryCatch({ ng_plot_check_reference_line(cl_no, label = "CHK_A", mean_axis = "y")
                   FALSE }, error = function(e) TRUE)
stopifnot(!drew)
# the same call with a real value DOES try to draw -- proving the probe works
would_draw <- tryCatch({ ng_plot_check_reference_line(11.5, label = "CHK_A", mean_axis = "y")
                         FALSE }, error = function(e) TRUE)
stopifnot(would_draw)

# ---- 5. no check figure on disk, and none announced in the index --------
pngs_no <- basename(list.files(d_no, pattern = "[.]png$", recursive = TRUE))
stopifnot(!any(grepl("check", pngs_no, ignore.case = TRUE)))
if (!is.null(sh_no) && "Figure_Index" %in% sh_no) {
  fi <- openxlsx::read.xlsx(file.path(d_no, "crossing_plan.xlsx"),
                            sheet = "Figure_Index", startRow = 2L)
  stopifnot(!any(grepl("check", unlist(lapply(fi, as.character)), ignore.case = TRUE)))
}

cat("no_check_means_no_check_artefacts: PASS\n")
