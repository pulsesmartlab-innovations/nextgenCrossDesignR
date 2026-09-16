# Every entry point a frontend is expected to call must be EXPORTED.
#
# The workbench reaches this package as an installed dependency, so it can only call
# what NAMESPACE exports. An internal entry point leaves it two bad options: `:::`,
# which R CMD check flags and which couples a frontend to an internal carrying no
# stability promise, or reimplementing the thing.
#
# ng_preview_marker_effects() was exactly that. It exists so a breeder can ask "should
# I launch this?" in seconds rather than hours -- it runs QC and the marker fits, shares
# the fit path and seed with the run it previews, and stops before scoring a single
# cross. On a 17-trait panel it previews a run that takes hours. And the frontend could
# not call it, because it was never exported, which was discovered only when the
# frontend went to use it.
#
# The other two entry points were exported, so nothing announced the inconsistency.
# This test is what makes the set explicit: adding a frontend-facing function without
# exporting it now fails here rather than at integration time.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# The public surface a frontend depends on, each with the reason it is public.
frontend_entry_points <- c(
  ng_run_cross_prediction =
    "the run itself; tools/run_cross_prediction_json.R is the headless wrapper around it",
  ng_write_backend_capability_registry_json =
    "the capability registry a frontend renders its controls from",
  ng_preview_marker_effects =
    "answers 'should I launch this?' before the expensive work, sharing the fit path and seed with the run"
)

ns <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
exported <- sub("^export\\(([^)]*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))

missing <- setdiff(names(frontend_entry_points), exported)
if (length(missing)) {
  stop("frontend-facing entry point(s) not exported, so a frontend can only reach them ",
       "through ::: -- ",
       paste(sprintf("%s (%s)", missing, frontend_entry_points[missing]), collapse = "; "),
       call. = FALSE)
}

# ...and they must actually be callable from the namespace, not merely named in
# NAMESPACE. A stale export line naming a function that no longer exists would pass a
# text check and fail at load.
for (nm in names(frontend_entry_points)) {
  stopifnot(exists(nm, envir = asNamespace("nextgenCrossDesign"), inherits = FALSE) ||
            exists(nm, envir = .GlobalEnv, inherits = FALSE))
}

cat("frontend_entry_points_are_exported: PASS  (", length(frontend_entry_points),
    "entry points )\n")
