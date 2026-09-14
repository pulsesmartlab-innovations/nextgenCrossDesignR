# The vignette is documentation nothing executes (`eval = FALSE`), so every claim
# it makes about the engine is unchecked by construction. It had drifted in four
# places at once:
#
#   * it named "var_complex" the default; the default is "usefulness";
#   * it described trait_value_metric = "pmv"/"vpm" as usefulness criteria, when
#     both rank on the RAW within-family variance with no cross mean in the score;
#   * it listed uc_variance_source = "le" -> parent_distance as a working setting;
#     the engine hard-errors on that pair at config time;
#   * it said var_complex falls back to "recombination variance or parent_distance";
#     its candidate list is c(pmv, vpm) and has never included parent_distance.
#
# A breeder configuring a run from that page would have chosen a metric that does
# something other than what they read. For a package meant for worldwide use, a
# doc that misdescribes the engine is a defect, not a cosmetic issue.
#
# This test cannot check prose. It checks the mechanical claims -- the token
# vocabulary and the stated default -- against formals(), which is where the four
# defects actually lived.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

vig <- readLines(file.path("vignettes", "nextgenCrossDesign.Rmd"), warn = FALSE)
txt <- paste(vig, collapse = "\n")
f <- formals(ng_run_cross_prediction)
allowed_metric <- eval(f$trait_value_metric)
allowed_source <- eval(f$uc_variance_source)

# ---- 1. every token the vignette shows is one the engine accepts -----------
grab <- function(param) {
  m <- regmatches(txt, gregexpr(sprintf('%s\\s*=\\s*"[^"]*"', param), txt))[[1L]]
  unique(sub('.*"([^"]*)".*', "\\1", m))
}
bad <- setdiff(grab("trait_value_metric"), allowed_metric)
if (length(bad)) stop("vignette shows trait_value_metric values the engine rejects: ",
                      paste(bad, collapse = ", "))
bad <- setdiff(grab("uc_variance_source"), allowed_source)
if (length(bad)) stop("vignette shows uc_variance_source values the engine rejects: ",
                      paste(bad, collapse = ", "))

# ---- 2. the stated default is the real default ----------------------------
# match.arg takes the first element, so formals() is the authority.
default_metric <- allowed_metric[[1L]]
default_source <- allowed_source[[1L]]
dl <- grep("DEFAULT", vig, value = TRUE)
if (!length(dl)) stop("vignette no longer states which trait_value_metric is the default")
if (!any(grepl(sprintf('"%s"', default_metric), dl, fixed = TRUE))) {
  stop("vignette names a default other than the real one (", default_metric, ")")
}

# ---- 3. the refused pair is documented AS refused, not as a setting --------
# usefulness + parent_distance is rejected in ng_cp__build_ctx(). The vignette must
# not present it as a choice; it previously did, under the "le" spelling.
i <- grep("uc_variance_source\\s*=\\s*\"(parent_distance|le)\"", vig)
if (!length(i)) stop("vignette no longer mentions the parent_distance/le variance source at all")
ctx <- paste(vig[i], collapse = " ")
if (!grepl("REFUSED|refused", ctx)) {
  stop("vignette lists uc_variance_source = parent_distance/le without marking it refused")
}

# ---- 4. it does not claim a fallback var_complex does not have -------------
if (grepl("falling back to recombination variance or `parent_distance`", txt, fixed = TRUE)) {
  stop("vignette claims var_complex falls back to parent_distance; its candidates are c(pmv, vpm)")
}

# ---- 5. and the claim itself is still true of the code --------------------
# Pin the fallback chain the corrected prose now describes, so the doc and the
# function cannot drift apart again from the OTHER side.
st <- data.frame(pmv = 1, vpm = 2)
stopifnot(identical(ng_run_cp_var_complex_col(st, "fast"), "pmv"))
stopifnot(identical(ng_run_cp_var_complex_col(data.frame(vpm = 2), "fast"), "vpm"))
stopifnot(inherits(tryCatch(ng_run_cp_var_complex_col(data.frame(parent_distance = 3), "fast"),
                            error = function(e) e), "error"))

# ---- 6. pmv/vpm really do rank on the raw variance ------------------------
sc <- data.frame(cross_mean_blend = c(100, 200), pmv = c(9, 1), vpm = c(4, 1))
stopifnot(identical(ng_run_cp_trait_value(sc, "maximize", "pmv"), c(9, 1)))
stopifnot(identical(ng_run_cp_trait_value(sc, "maximize", "vpm"), c(4, 1)))
# ...whereas usefulness puts the mean back in
u <- ng_run_cp_trait_value(sc, "maximize", "usefulness", uc_variance_source = "pmv")
stopifnot(all(u > c(100, 200)))

cat("vignette_matches_the_engine: PASS\n")
