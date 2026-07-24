## Progeny system: two genetic variance models (DH, RIL). The plural spellings
## DHs / RILs are accepted as backend aliases but are no longer separate dropdown
## choices (they resolve to the same DH / RIL recombination-variance kernels).
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

## ---- aliases still map to the two kernels ----
stopifnot(identical(ng_run_cp_target("DH"),   "DH"))
stopifnot(identical(ng_run_cp_target("DHs"),  "DH"))   # plural alias
stopifnot(identical(ng_run_cp_target("RIL"),  "RIL"))
stopifnot(identical(ng_run_cp_target("RILs"), "RIL"))  # plural alias
stopifnot(identical(ng_run_cp_target("doubled_haploids"), "DH"))

## ---- the dropdown (registry) now offers exactly the two models ----
ctl <- ng_backend_capability_registry()$controls
prog <- ctl[[which(vapply(ctl, function(x) identical(x$id, "progeny"), logical(1)))]]
choice_values <- vapply(prog$choices, function(ch) ch$value, character(1L))
stopifnot(setequal(choice_values, c("DH", "RIL")))

cat("progeny_alias.R: PASS\n")
