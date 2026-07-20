# The backend must stay free of any Shiny dependency: the workbench frontend is a
# separate package (nextgenCrossWorkbench) that drives this engine out-of-process
# via tools/run_cross_prediction_json.R, so Shiny must never leak in here.
source(file.path("tests", "helper_load.R"))

stopifnot(!exists("ng_shiny_app_dir", mode = "function", inherits = TRUE))
stopifnot(!exists("ng_run_shiny_app", mode = "function", inherits = TRUE))
stopifnot(!exists("ng_run_app", mode = "function", inherits = TRUE))

desc <- read.dcf(file.path(root, "DESCRIPTION"))[1, ]
deps <- paste(desc[c("Depends", "Imports", "Suggests", "Enhances")], collapse = "\n")
stopifnot(!grepl("\\bshiny\\b", deps, ignore.case = TRUE))

backend_app_dir <- file.path(root, "inst", "shiny", "nextgenCrossDesign")
stopifnot(!dir.exists(backend_app_dir))

cat("shiny boundary tests passed\n")
