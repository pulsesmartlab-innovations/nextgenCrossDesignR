# A remedy the user cannot look up is not a remedy.
#
# The reliability gate refuses a run and tells the breeder what would fix it:
# supply training_genotype / training_phenotype, or switch trait_value_metric, or
# lower min_cv_predictive_r2 deliberately. Those instructions are the entire value
# of refusing rather than silently substituting -- the plan's words were "a training
# set is the only way across that line. Nothing says so."
#
# So every parameter the gate names must be findable in a surface that SHIPS with
# the installed package. docs/ is .Rbuildignore'd: a user who installs
# nextgenCrossDesign never sees docs/BACKEND_USER_GUIDE.md or the frontend contract.
# Only man/ and vignettes/ reach them. Documenting a remedy solely in docs/ leaves
# the error message pointing at something the recipient cannot read.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# ---- gather every parameter the gate's remedies actually name --------------
es <- data.frame(
  trait = c("A", "B", "C"),
  cv_predictive_r2 = c(NA_real_, -0.05, 0.10),
  marker_effect_training_n = c(6L, 40L, 40L),
  cv_available = c(FALSE, TRUE, TRUE),
  cv_unavailable_reason = c("n_lt_10", NA, NA),
  stringsAsFactors = FALSE)
v <- ng_evaluate_marker_reliability(es, trait_value_metric = "usefulness",
                                    min_cv_predictive_r2 = 0.35)
stopifnot(nrow(v) == 3L, any(v$verdict == "refuse"))
remedies <- paste(stats::na.omit(v$remedy), collapse = " ")
stopifnot(nzchar(remedies))

formals_names <- names(formals(ng_run_cross_prediction))
tokens <- unique(regmatches(remedies, gregexpr("[A-Za-z_][A-Za-z0-9_]*", remedies))[[1L]])
named <- intersect(tokens, formals_names)
if (!length(named)) stop("the gate's remedies name no parameter at all -- ",
                         "a refusal with no actionable instruction", call. = FALSE)

# ---- each must appear in a surface the installed package carries -----------
shipped <- character()
for (d in c("man", "vignettes")) {
  fs <- list.files(file.path(root, d), full.names = TRUE, recursive = TRUE)
  for (fp in fs) shipped <- c(shipped, readLines(fp, warn = FALSE))
}
shipped <- paste(shipped, collapse = "\n")

missing <- named[!vapply(named, function(nm) grepl(nm, shipped, fixed = TRUE), logical(1))]
if (length(missing)) {
  stop("the reliability gate tells users to reach for parameter(s) that no shipped ",
       "documentation (man/, vignettes/) mentions: ", paste(missing, collapse = ", "),
       " -- docs/ is .Rbuildignore'd and does not reach an installed user",
       call. = FALSE)
}

# ---- and the gate itself must be documented, not just its parameters ------
# Knowing min_cv_predictive_r2 exists does not tell a breeder that a run can be
# REFUSED outright, nor that two metrics are exempt from the refusal.
for (claim in c("min_cv_predictive_r2", "ng_preview_marker_effects", "parent_distance")) {
  if (!grepl(claim, shipped, fixed = TRUE)) {
    stop("shipped documentation never mentions ", claim,
         ", so the gate's behaviour is undiscoverable from an installed package",
         call. = FALSE)
  }
}

cat("gate_remedies_are_documented: PASS  (", paste(named, collapse = ", "), ")\n")
