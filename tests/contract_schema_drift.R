# Repository-only check: the frontend's parameter menu must stay in sync with
# the entrypoint it documents.
#
# docs/frontend/contracts/config_schema.json is what a frontend developer builds
# a run config from, and tools/run_cross_prediction_json.R rejects unknown keys.
# So a parameter documented there but absent from ng_run_cross_prediction()'s
# formals is a config the frontend can legally build and the backend will refuse
# at runtime. This lives here rather than in tests/testthat/ because docs/ is
# .Rbuildignore'd and so is not present in the installed package.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

schema_path <- file.path(root, "docs", "frontend", "contracts", "config_schema.json")
stopifnot(file.exists(schema_path))
schema <- jsonlite::fromJSON(schema_path, simplifyVector = FALSE)

stopifnot(identical(schema$schema, "ng_run_config.v1"))
stopifnot(identical(schema$entrypoint, "ng_run_cross_prediction"))

# The runner named by the contract must exist in this repository.
stopifnot(file.exists(file.path(root, schema$runner)))

declared <- unlist(lapply(schema$groups, function(g)
  vapply(g$params, function(p) p$name, character(1))), use.names = FALSE)
types <- unlist(lapply(schema$groups, function(g)
  vapply(g$params, function(p) if (is.null(p$type)) "" else p$type, character(1))),
  use.names = FALSE)

# "meta" entries document behaviour (e.g. optional-package gating) rather than
# naming an argument, so they are not expected to be formals.
declared <- declared[types != "meta"]
stopifnot(!anyDuplicated(declared))

formals_names <- names(formals(ng_run_cross_prediction))
undocumented_ok <- c(
  # In-memory alternatives to the *_file keys: a JSON config cannot carry a
  # matrix or data frame, so the contract deliberately exposes paths only.
  "genotype", "marker_map", "phenotype", "trait_direction",
  "training_genotype", "training_phenotype",
  # check_geno / check_records are also in-memory matrices/lists (a check genotype matrix
  # and per-source check phenotypic records); check_progeny_size IS documented in the schema.
  "check_geno", "check_records"
)

not_a_formal <- setdiff(declared, formals_names)
if (length(not_a_formal)) {
  stop("config_schema.json documents parameter(s) that ng_run_cross_prediction() ",
       "does not accept, so the headless runner would reject them: ",
       paste(not_a_formal, collapse = ", "), call. = FALSE)
}

missing_from_schema <- setdiff(formals_names, c(declared, undocumented_ok))
if (length(missing_from_schema)) {
  message("NOTE: ng_run_cross_prediction() argument(s) absent from ",
          "config_schema.json (frontend cannot discover them): ",
          paste(missing_from_schema, collapse = ", "))
}

cat("contract_schema_drift:",
    length(declared), "documented parameters, all valid arguments;",
    length(missing_from_schema), "argument(s) undocumented\n")
