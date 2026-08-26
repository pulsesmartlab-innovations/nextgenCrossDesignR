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

# --- portfolio + risk columns reach the workbook (0.19.0) -----------------------------------
# The workbook is what leaves the building. Before this, a breeder working from the spreadsheet
# saw priority tiers with no indication of which crosses were speculative.
set.seed(21)
n_p <- 12L; mk_p <- 40L; gid_p <- sprintf("Q%02d", seq_len(n_p))
gm_p <- matrix(2L * rbinom(n_p * mk_p, 1, 0.5), n_p, mk_p,
               dimnames = list(gid_p, sprintf("S%03d", seq_len(mk_p))))
bA <- rnorm(mk_p, 0, 0.1); bB <- -0.7 * bA + rnorm(mk_p, 0, 0.03)
phen_p <- data.frame(NAME = gid_p,
                     yield = as.numeric(gm_p %*% bA) + rnorm(n_p, 0, 0.5),
                     protein = as.numeric(gm_p %*% bB) + rnorm(n_p, 0, 0.5),
                     stringsAsFactors = FALSE)
res_p <- ng_run_cross_prediction(
  phenotype = phen_p,
  genotype = data.frame(NAME = gid_p, gm_p, check.names = FALSE, stringsAsFactors = FALSE),
  marker_map = data.frame(SNP = colnames(gm_p), chr = rep(1:2, length.out = mk_p),
                          bp = rep(seq(0, 100, length.out = 20), 2)[seq_len(mk_p)] * 1e6),
  trait_direction = data.frame(trait = c("yield", "protein"), column = c("yield", "protein"),
                               direction = c("increase", "decrease"), weight = c(0.6, 0.4)),
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 6L,
  max_crosses_per_parent = 3L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)

ti_p <- ng_cpw_trait_table(stats::setNames(c("increase", "decrease"), c("yield", "protein")),
                           res_p$selected_crosses)
sel_p <- ng_cpw_make_selected(res_p$selected_crosses, ti_p, NULL, NULL, 10L)
pf_cols <- c("portfolio_profile", "cross_level", "cross_upside", "risk_bin_within_candidate_pool",
             "cross_confidence_within_candidate_pool", "risk_driver_trait", "risk_driver_share",
             "portfolio_basis")
stopifnot(all(pf_cols %in% names(sel_p)))
stopifnot(nrow(sel_p) == nrow(res_p$selected_crosses))
# Values survive the transfer, not just the headers.
stopifnot(all(sel_p$risk_bin_within_candidate_pool %in% c("low", "med", "high")))
stopifnot(all(sel_p$portfolio_profile %in%
                c("breakthrough", "workhorse", "long_shot", "deprioritize")))
stopifnot(all(is.finite(sel_p$cross_upside)), all(sel_p$cross_upside >= 0))
stopifnot(all(sel_p$risk_driver_trait %in% c("yield", "protein")))
stopifnot(identical(sel_p$portfolio_basis[[1L]],
                    as.character(res_p$selected_crosses$portfolio_basis[[1L]])))
# Decision/free-text columns must stay rightmost so the breeder's writing area is unbroken.
for (cc in pf_cols) stopifnot(match(cc, names(sel_p)) < match("Breeder_Rationale", names(sel_p)))
stopifnot(match("Crossing_Status", names(sel_p)) > match("Breeder_Rationale", names(sel_p)))

# A run WITHOUT the annotation (pre-0.13 result, or an unresolvable index) must simply omit the
# columns rather than emit a sheet full of NA.
bare_p <- res_p$selected_crosses
bare_p[c("portfolio_profile", "cross_level", "cross_upside", "risk_bin", "cross_confidence",
         "risk_driver_trait", "risk_driver_share", "portfolio_basis")] <- NULL
sel_bare <- ng_cpw_make_selected(bare_p, ti_p, NULL, NULL, 10L)
stopifnot(!any(pf_cols %in% names(sel_bare)))
stopifnot("Breeder_Rationale" %in% names(sel_bare), nrow(sel_bare) == nrow(bare_p))

# An all-NA column is treated as absent (nothing to tell the breeder).
na_p <- res_p$selected_crosses
na_p$risk_driver_trait <- NA_character_
sel_na <- ng_cpw_make_selected(na_p, ti_p, NULL, NULL, 10L)
stopifnot(!("risk_driver_trait" %in% names(sel_na)), "portfolio_profile" %in% names(sel_na))
cat("workbook portfolio + risk column test passed\n")
