# Every EXPORTED object must have a documentation alias.
#
# R CMD check enforces this ("checking for missing documentation entries"), and it is
# the ONLY tier that does. Neither tests/testthat nor tools/run_tests.R loads the man/
# directory at all, so both stayed green while CI went red on two Ubuntu legs -- the
# 0.35.0 export of ng_preview_marker_effects added a NAMESPACE line and no alias.
#
# That is a whole class of failure the local tiers are blind to, and it costs a full
# CI round trip to learn about. This test closes it: it reads the same two files
# R CMD check reads and reaches the same verdict in milliseconds.
#
# It checks alias coverage only, not prose quality. An alias is what makes `?fn` work
# and what check counts; whether the paragraph is any good is a review question, not
# something a parser can settle.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- ng_test_find_root()

# Parse NAMESPACE rather than asking the loaded namespace: R/load.R source()s every
# file into one environment, so getNamespaceExports() is unavailable here and ls()
# would return internals too. The file is the same declaration R CMD check reads.
ns <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
exported <- sub("^export\\(([^)]*)\\).*$", "\\1", grep("^export\\(", ns, value = TRUE))
exported <- gsub("[\"'`]", "", trimws(exported))
stopifnot(length(exported) > 100L)

rd <- list.files(file.path(root, "man"), pattern = "[.]Rd$", full.names = TRUE)
stopifnot(length(rd) > 0L)
aliased <- unlist(lapply(rd, function(f) {
  ln <- grep("^\\\\alias\\{", readLines(f, warn = FALSE), value = TRUE)
  gsub("[\"'`]", "", trimws(sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1", ln)))
}), use.names = FALSE)

undocumented <- setdiff(exported, aliased)
if (length(undocumented)) {
  stop(sprintf(
    "exported but not documented -- R CMD check will fail with 'Undocumented code objects': %s\nAdd \\alias{<name>} to man/nextgenCrossDesign-api.Rd (and say what it is for in \\details).",
    paste(undocumented, collapse = ", ")), call. = FALSE)
}

# The converse is not an error -- an alias may document a concept or an internal that
# is deliberately unexported -- but an alias naming nothing at all is a typo that
# would silently stop covering the export it was meant to cover.
api <- file.path(root, "man", "nextgenCrossDesign-api.Rd")
api_alias <- setdiff(gsub("[\"'`]", "", trimws(sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1",
  grep("^\\\\alias\\{", readLines(api, warn = FALSE), value = TRUE)))),
  "nextgenCrossDesign-api")
stray <- setdiff(api_alias, exported)
stopifnot(length(stray) == 0L)

cat("PASS: every_export_is_documented\n")
