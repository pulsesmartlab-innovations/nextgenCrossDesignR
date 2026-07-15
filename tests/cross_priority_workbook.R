root_candidates <- c(
  file.path(getwd(), "nextgen_cross_design"),
  getwd(),
  file.path("..", "nextgen_cross_design"),
  file.path("..")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

crosses <- data.frame(
  parent1 = sprintf("P%02d", seq_len(8L)),
  parent2 = sprintf("Q%02d", seq_len(8L)),
  multi_trait_score = seq(0.90, 0.55, length.out = 8L),
  pair_kinship = seq(-0.30, 0.10, length.out = 8L),
  multi_trait_threshold_violation = c(0, 0, 0.05, 0, 0.2, 0, 0, 0.4),
  pred_YIELD = seq(100, 86, length.out = 8L),
  pred_DISEASE = seq(2, 6, length.out = 8L),
  stringsAsFactors = FALSE
)
crosses <- ng_rank_cross_priority(
  crosses,
  breaks = c(0.25, 0.50, 0.75, 1),
  kinship_weight = 0.10
)

trait_directions <- data.frame(
  trait = c("YIELD", "DISEASE"),
  column = c("pred_YIELD", "pred_DISEASE"),
  direction = c("maximize", "minimize"),
  weight = c(1, 2),
  stringsAsFactors = FALSE
)

parent_use <- data.frame(
  parent = c("P01", "P02", "Q01"),
  crosses_selected = c(2L, 1L, 2L),
  stringsAsFactors = FALSE
)

duplicate_pairs <- data.frame(
  parent1 = "P07",
  parent2 = "Q07",
  similarity = 0.996,
  markers_compared = 1200L,
  stringsAsFactors = FALSE
)

tables <- ng_cross_priority_workbook_tables(
  crosses = crosses,
  scored = crosses,
  trait_directions = trait_directions,
  parent_use = parent_use,
  duplicate_pairs = duplicate_pairs,
  n_crosses_requested = 8L,
  block_size = 4L
)

stopifnot(all(c(
  "Dashboard", "Scoring_Method", "Trait_Directions", "Selected_All",
  "Highly_Priority", "Priority", "Medium_Priority", "Low_Priority",
  "Candidate_Crosses", "Parent_Use_QC", "Duplicate_QC"
) %in% names(tables)))

selected <- tables$Selected_All
stopifnot(nrow(selected) == 8L)
stopifnot(all(c(
  "Cross_ID", "Cross_Block", "Block_Position", "priority_tier",
  "Breeder_Rationale", "Breeder_Notes", "Final_Decision", "Crossing_Status",
  "top_favorable_traits", "top_risk_traits"
) %in% names(selected)))
stopifnot(!("Breeding.Rationale" %in% names(selected)))
stopifnot(!("Breeding_Rationale" %in% names(selected)))
stopifnot(all(selected$Breeder_Rationale == ""))
stopifnot(all(selected$Final_Decision == ""))

stopifnot(nrow(tables$Highly_Priority) == 2L)
stopifnot(nrow(tables$Priority) == 2L)
stopifnot(nrow(tables$Medium_Priority) == 2L)
stopifnot(nrow(tables$Low_Priority) == 2L)
stopifnot(nrow(tables$Duplicate_QC) == 1L)
stopifnot(any(grepl("editable", tables$Scoring_Method$Details, ignore.case = TRUE)))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  out_xlsx <- tempfile("cross_priority_workbook_", fileext = ".xlsx")
  written <- ng_write_cross_priority_workbook(
    output_path = out_xlsx,
    crosses = crosses,
    scored = crosses,
    trait_directions = trait_directions,
    parent_use = parent_use,
    duplicate_pairs = duplicate_pairs,
    n_crosses_requested = 8L,
    block_size = 4L
  )
  stopifnot(file.exists(written))
  xml_dir <- tempfile("cross_priority_workbook_xml_")
  dir.create(xml_dir)
  utils::unzip(written, files = "xl/workbook.xml", exdir = xml_dir)
  workbook_xml <- paste(readLines(file.path(xml_dir, "xl", "workbook.xml"), warn = FALSE), collapse = "\n")
  stopifnot(grepl("Highly_Priority", workbook_xml, fixed = TRUE))
  stopifnot(grepl("Low_Priority", workbook_xml, fixed = TRUE))
  all_xml <- paste(unlist(lapply(utils::unzip(written, list = TRUE)$Name, function(path) {
    if (!grepl("[.]xml$", path)) return("")
    con <- unz(written, path)
    on.exit(close(con), add = TRUE)
    paste(readLines(con, warn = FALSE), collapse = "\n")
  })), collapse = "\n")
  stopifnot(!grepl("Breeding.Rationale", all_xml, fixed = TRUE))
}

cat("cross priority workbook tests passed\n")
