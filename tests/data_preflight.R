# Repository root FIRST. Checking `nextgen_cross_design/` ahead of getwd() made
# these tests load a STALE 0.19.0 copy of the package that sits in the working
# tree under exactly that name -- so they validated a package eleven versions old
# while appearing to cover the current one. The ones that failed were the lucky
# case; the ones that passed gave false assurance. ng_load() already resolves in
# this order; only these hand-rolled preambles inverted it.
root_candidates <- c(
  getwd(),
  file.path(".."),
  file.path(getwd(), "nextgen_cross_design"),
  file.path("..", "nextgen_cross_design")
)
root_hits <- root_candidates[file.exists(file.path(root_candidates, "R", "load.R"))]
stopifnot(length(root_hits) > 0L)
root <- normalizePath(root_hits[[1]], mustWork = TRUE)
source(file.path(root, "tools", "ng_project_libpath.R")); ng_prepend_project_lib(file.path(dirname(root), ".Rlib"))
source(file.path(root, "R", "load.R"))
ng_load(root, use_cpp = FALSE, verbose = FALSE)

geno <- data.frame(
  parent = c("P01", "P01", "P03"),
  M1 = c(0, 1, 2),
  M1 = c(0, 1, 2),
  M3 = c(0, 5, 2),
  check.names = FALSE
)
phenotype <- data.frame(
  parent = c("P01", "P02", "P02"),
  trait = c("yield", "yield", "yield"),
  value = c(10, 11, 12),
  stringsAsFactors = FALSE
)
candidate_pairs <- data.frame(
  parent1 = c("P01", "P02", "P02", "P03"),
  parent2 = c("P02", "P01", "P04", "P03"),
  stringsAsFactors = FALSE
)
trait_spec <- data.frame(
  trait = c("yield", "yield", "disease"),
  direction = c("maximize", "increase", "bad_direction"),
  min_value = c(NA, 12, 5),
  max_value = c(NA, 10, 2),
  stringsAsFactors = FALSE
)
marker_map <- data.frame(
  marker = c("M1", "M1", "M4"),
  chr = c(1, 1, 1),
  pos_cm = c(0, 1, 2),
  stringsAsFactors = FALSE
)
parent_kinship <- matrix(
  diag(3),
  nrow = 3,
  dimnames = list(c("P01", "P01", "P03"), c("P01", "P02", "P03"))
)

result <- ng_preflight_input_tables(
  geno = geno,
  phenotype = phenotype,
  candidate_pairs = candidate_pairs,
  trait_spec = trait_spec,
  marker_map = marker_map,
  parent_kinship = parent_kinship,
  ploidy = 2L,
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)

stopifnot(identical(result$schema_version, "ng_data_preflight.v1"))
stopifnot(identical(result$generated_at, "2026-05-07T00:00:00Z"))
stopifnot(identical(result$status, "blocker"))
stopifnot(is.data.frame(result$issues))
ids <- result$issues$id
stopifnot(all(c(
  "duplicate_parent_ids",
  "duplicate_marker_ids",
  "invalid_dosage_range",
  "duplicate_phenotype_rows",
  "genotype_phenotype_id_mismatch",
  "duplicate_candidate_crosses",
  "reciprocal_candidate_crosses",
  "unknown_candidate_parent",
  "self_cross",
  "duplicate_trait_names",
  "invalid_trait_direction",
  "invalid_trait_threshold",
  "duplicate_marker_map_entries",
  "missing_marker_map_entries",
  "parent_relationship_duplicate_ids",
  "parent_relationship_id_mismatch"
) %in% ids))
stopifnot(result$counts$blockers >= 10L)
stopifnot(result$tables$rows[result$tables$table == "geno"] == 3L)

clean <- ng_preflight_input_tables(
  geno = data.frame(parent = c("P01", "P02", "P03"), M1 = c(0, 1, 2), M2 = c(2, 1, 0)),
  phenotype = data.frame(parent = c("P01", "P02", "P03"), trait = c("yield", "yield", "yield"), value = c(10, 11, 12)),
  candidate_pairs = data.frame(parent1 = c("P01", "P02"), parent2 = c("P02", "P03")),
  trait_spec = data.frame(trait = c("yield", "disease"), direction = c("maximize", "minimize"), min_value = c(NA, NA), max_value = c(NA, NA)),
  marker_map = data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, 1)),
  ploidy = 2L,
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)
stopifnot(identical(clean$status, "pass"))
stopifnot(nrow(clean$issues) == 0L)

putative_dup <- ng_preflight_input_tables(
  geno = data.frame(parent = c("P01", "P02", "P03"), M1 = c(0, 0, 2), M2 = c(1, 1, 2), M3 = c(2, 2, 0), M4 = c(0, 0, 1)),
  phenotype = data.frame(parent = c("P01", "P02", "P03"), trait = c("yield", "yield", "yield"), value = c(10, 11, 12)),
  marker_map = data.frame(marker = c("M1", "M2", "M3", "M4"), chr = c(1, 1, 1, 1), pos_cm = 0:3),
  putative_duplicate_check = TRUE,
  duplicate_threshold = 0.995,
  duplicate_min_compared_markers = 4L,
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)
stopifnot("putative_duplicate_genotypes" %in% putative_dup$issues$id)
stopifnot(putative_dup$counts$warnings >= 1L)

remove_dup <- ng_preflight_input_tables(
  geno = data.frame(parent = c("P01", "P02", "P03"), M1 = c(0, 0, 2), M2 = c(1, 1, 2), M3 = c(2, 2, 0), M4 = c(0, 0, 1)),
  phenotype = data.frame(parent = c("P01", "P02", "P03"), trait = c("yield", "yield", "yield"), value = c(10, 11, 12)),
  candidate_pairs = data.frame(parent1 = c("P01", "P01", "P02"), parent2 = c("P02", "P03", "P03")),
  marker_map = data.frame(marker = c("M1", "M2", "M3", "M4"), chr = c(1, 1, 1, 1), pos_cm = 0:3),
  parent_kinship = matrix(
    diag(3),
    nrow = 3,
    dimnames = list(c("P01", "P02", "P03"), c("P01", "P02", "P03"))
  ),
  putative_duplicate_check = TRUE,
  putative_duplicate_action = "remove",
  duplicate_threshold = 0.995,
  duplicate_min_compared_markers = 4L,
  generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
)
stopifnot("putative_duplicate_genotypes_removed" %in% remove_dup$issues$id)
stopifnot(!is.null(remove_dup$cleaned_tables$geno))
stopifnot(!is.null(remove_dup$cleaning$putative_duplicates$removed_parents))
stopifnot(identical(remove_dup$cleaning$putative_duplicates$removed_parents$removed_parent, "P02"))
stopifnot(identical(remove_dup$cleaning$putative_duplicates$removed_parents$kept_parent, "P01"))
stopifnot(nrow(remove_dup$cleaned_tables$geno) == 2L)
stopifnot(identical(remove_dup$cleaned_tables$geno$parent, c("P01", "P03")))
stopifnot(identical(remove_dup$cleaned_tables$phenotype$parent, c("P01", "P03")))
stopifnot(nrow(remove_dup$cleaned_tables$candidate_pairs) == 1L)
stopifnot(identical(remove_dup$cleaned_tables$candidate_pairs$parent1, "P01"))
stopifnot(identical(remove_dup$cleaned_tables$candidate_pairs$parent2, "P03"))
stopifnot(identical(rownames(remove_dup$cleaned_tables$parent_kinship), c("P01", "P03")))
stopifnot(identical(colnames(remove_dup$cleaned_tables$parent_kinship), c("P01", "P03")))
stopifnot(remove_dup$cleaning$putative_duplicates$rows_removed$rows_removed[
  remove_dup$cleaning$putative_duplicates$rows_removed$table == "geno"
] == 1L)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  tmp <- tempfile("ng_data_preflight_")
  dir.create(tmp, recursive = TRUE)
  write.csv(geno, file.path(tmp, "geno.csv"), row.names = FALSE, quote = FALSE)
  write.csv(phenotype, file.path(tmp, "phenotype.csv"), row.names = FALSE, quote = FALSE)
  write.csv(candidate_pairs, file.path(tmp, "candidate_pairs.csv"), row.names = FALSE, quote = FALSE)
  write.csv(trait_spec, file.path(tmp, "trait_spec.csv"), row.names = FALSE, quote = FALSE)
  write.csv(marker_map, file.path(tmp, "marker_map.csv"), row.names = FALSE, quote = FALSE)
  cli_out <- file.path(tmp, "preflight.json")
  env <- c(
    NG_PREFLIGHT_GENO = file.path(tmp, "geno.csv"),
    NG_PREFLIGHT_PHENOTYPE = file.path(tmp, "phenotype.csv"),
    NG_PREFLIGHT_CANDIDATE_PAIRS = file.path(tmp, "candidate_pairs.csv"),
    NG_PREFLIGHT_TRAIT_SPEC = file.path(tmp, "trait_spec.csv"),
    NG_PREFLIGHT_MARKER_MAP = file.path(tmp, "marker_map.csv"),
    NG_PREFLIGHT_PLOIDY = "2",
    NG_PREFLIGHT_OUT = cli_out
  )
  old_env <- Sys.getenv(names(env), unset = NA_character_)
  restore_env <- function() {
    for (name in names(old_env)) {
      if (is.na(old_env[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, as.list(stats::setNames(old_env[[name]], name)))
    }
  }
  on.exit(restore_env(), add = TRUE)
  do.call(Sys.setenv, as.list(env))
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
  cli_output <- system2(
    rscript,
    normalizePath(file.path(root, "tools", "run_data_preflight.R"), winslash = "/", mustWork = FALSE),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(cli_output, "status")
  if (is.null(status)) status <- 0L
  if (!identical(as.integer(status), 0L)) print(cli_output)
  stopifnot(status == 0L)
  stopifnot(file.exists(cli_out))
  payload <- jsonlite::fromJSON(cli_out, simplifyVector = FALSE)
  stopifnot(identical(payload$schema_version, "ng_data_preflight.v1"))
  stopifnot(identical(payload$status, "blocker"))

  remove_json <- file.path(tmp, "preflight_remove.json")
  ng_write_data_preflight_json(
    output_path = remove_json,
    geno = data.frame(parent = c("P01", "P02", "P03"), M1 = c(0, 0, 2), M2 = c(1, 1, 2), M3 = c(2, 2, 0), M4 = c(0, 0, 1)),
    phenotype = data.frame(parent = c("P01", "P02", "P03"), trait = c("yield", "yield", "yield"), value = c(10, 11, 12)),
    marker_map = data.frame(marker = c("M1", "M2", "M3", "M4"), chr = c(1, 1, 1, 1), pos_cm = 0:3),
    putative_duplicate_check = TRUE,
    putative_duplicate_action = "remove",
    duplicate_threshold = 0.995,
    duplicate_min_compared_markers = 4L,
    generated_at = as.POSIXct("2026-05-07 00:00:00", tz = "UTC")
  )
  remove_payload <- jsonlite::fromJSON(remove_json, simplifyVector = FALSE)
  stopifnot(is.null(remove_payload[["cleaned_tables"]]))
  stopifnot(!is.null(remove_payload$cleaned_tables_summary))
  stopifnot(!is.null(remove_payload$cleaning$putative_duplicates$removed_parents))
  stopifnot(identical(remove_payload$cleaning$putative_duplicates$removed_parents[[1]]$removed_parent, "P02"))
}

cat("data preflight tests passed\n")
