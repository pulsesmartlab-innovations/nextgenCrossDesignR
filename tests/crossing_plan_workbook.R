helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  message("openxlsx unavailable; skipping crossing plan workbook test")
  quit(status = 0)
}

candidates <- data.frame(
  line = paste0("NDP", 1:8),
  cross = c("N22P110-2", "N22P032-5", "N21P078", "N22P021-2", "N22P028-2", "N22P032-10", "N22P194-3", "N22P045-10"),
  seed_color = c("Yellow", "Yellow", "Green", "Yellow", "Yellow", "Yellow", "Green", "Green"),
  yield = c(87, 74, 77, 79, 74, 69, 72, 70),
  harvestability = c(2, 2, 2, 3, 2, 3, 3, 2),
  tsw = c(262, 282, 244, 253, 238, 277, 265, 209),
  stringsAsFactors = FALSE
)
criteria <- data.frame(
  trait = c("yield", "harvestability", "tsw"),
  label = c("Yield (bu/a)", "Harvestability", "TSW"),
  weight = c(0.5, 0.3, 0.2),
  direction = c("increase", "decrease", "increase"),
  meaning = c("Higher yield is better", "Lower lodging/harvestability score is better", "Higher seed weight is preferred"),
  stringsAsFactors = FALSE
)

tmp <- tempfile("ng_crossing_plan_", fileext = ".xlsx")
returned <- ng_write_crossing_plan_workbook(
  output_path = tmp,
  candidates = candidates,
  criteria = criteria,
  top_n = 5L,
  cross_n = 8L,
  title = "Top 5 Crossing Candidates with Crossing Plan - 2026"
)

stopifnot(file.exists(returned))
sheets <- openxlsx::getSheetNames(returned)
expected <- c("Top Parents", "Recommended Crosses", "Selection Index - All Lines", "Criteria", "Family Screening")
stopifnot(identical(sheets, expected))

top <- openxlsx::read.xlsx(returned, sheet = "Top Parents", startRow = 3)
stopifnot(nrow(top) == 5L)
stopifnot(!any(duplicated(top$Family.Root)))
stopifnot(all(top$Full.sib.exclusion.status == "Selected - unique family"))

crosses <- openxlsx::read.xlsx(returned, sheet = "Recommended Crosses", startRow = 3)
stopifnot(nrow(crosses) == 8L)
stopifnot(all(crosses$Full.Sib.Check == "Pass - different family roots"))
stopifnot(all(crosses$Parent.1 != crosses$Parent.2))
stopifnot(all(crosses$Tier %in% c("Tier 1", "Tier 2", "Tier 3")))

all_lines <- openxlsx::read.xlsx(returned, sheet = "Selection Index - All Lines", startRow = 3)
stopifnot(any(all_lines$Full.Sib.Excluded == "Yes"))
stopifnot(sum(all_lines$Selected.Top.5 == "Yes") == 5L)

family <- openxlsx::read.xlsx(returned, sheet = "Family Screening", startRow = 2)
stopifnot(any(family$Full.Sib.Excluded == "Yes"))
stopifnot(all(family$Within.Family.Rank >= 1L))

criteria_out <- openxlsx::read.xlsx(returned, sheet = "Criteria", startRow = 2)
stopifnot(any(criteria_out$Rule...Parameter == "Full-sib avoidance"))
stopifnot(any(criteria_out$Rule...Parameter == "Harvestability interpretation"))

entries <- utils::unzip(returned, list = TRUE)$Name
stopifnot(any(grepl("^xl/tables/table", entries)))
stopifnot(any(grepl("^xl/media/", entries)))

tmp_dir <- tempfile("ng_crossing_plan_cli_")
dir.create(tmp_dir, recursive = TRUE)
candidates_csv <- file.path(tmp_dir, "candidates.csv")
criteria_csv <- file.path(tmp_dir, "criteria.csv")
cli_out <- file.path(tmp_dir, "crossing_plan_cli.xlsx")
utils::write.csv(candidates, candidates_csv, row.names = FALSE)
utils::write.csv(criteria, criteria_csv, row.names = FALSE)
env <- c(
  NG_CROSSING_PLAN_CANDIDATES = normalizePath(candidates_csv, winslash = "/", mustWork = TRUE),
  NG_CROSSING_PLAN_CRITERIA = normalizePath(criteria_csv, winslash = "/", mustWork = TRUE),
  NG_CROSSING_PLAN_OUT = normalizePath(cli_out, winslash = "/", mustWork = FALSE),
  NG_CROSSING_PLAN_TOP_N = "4",
  NG_CROSSING_PLAN_CROSS_N = "5",
  NG_CROSSING_PLAN_TITLE = "CLI Crossing Plan"
)
old_env <- Sys.getenv(names(env), unset = NA_character_)
restore_env <- function() {
  for (name in names(old_env)) {
    if (is.na(old_env[[name]])) {
      Sys.unsetenv(name)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old_env[[name]], name)))
    }
  }
}
on.exit(restore_env(), add = TRUE)
do.call(Sys.setenv, as.list(env))
rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
cli_output <- system2(
  rscript,
  normalizePath(file.path(root, "tools", "export_crossing_plan_workbook.R"), winslash = "/", mustWork = TRUE),
  stdout = TRUE,
  stderr = TRUE
)
status <- attr(cli_output, "status")
if (is.null(status)) status <- 0L
if (!identical(as.integer(status), 0L)) print(cli_output)
stopifnot(length(status) == 1L, is.finite(status), status == 0L)
stopifnot(file.exists(cli_out))
cli_top <- openxlsx::read.xlsx(cli_out, sheet = "Top Parents", startRow = 3)
stopifnot(nrow(cli_top) == 4L)

cat("crossing plan workbook tests passed\n")
