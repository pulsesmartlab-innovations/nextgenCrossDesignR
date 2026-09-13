# The advisory channel's shape, and the fit-level facts it is built from.
#
# Phase 0 of the 0.32.0 transparency work. Nothing consumes these yet; this test
# pins the contract BEFORE anything depends on it, so the shape cannot drift into
# the runner and the frontend simultaneously.
#
# Two things are being established.
#
# 1. A post-effects advisory channel with the SAME column shape as qc$issues
#    (R/26_data_preflight.R:5-31). qc$issues is the right shape but the wrong
#    lifetime: it is finalised in the qc stage, long before marker effects exist,
#    so it cannot carry "this trait's markers do not predict" or "a training set
#    would help". Matching its columns means any consumer that can already render
#    qc issues renders advisories with the same code.
#
#    The severity domain is deliberately one level richer -- advice / warning /
#    blocker -- because most of what this channel carries is genuinely advice, and
#    calling advice a warning devalues the warnings.
#
# 2. Whether the cross-validation behind cv_predictive_r2 could run AT ALL, and if
#    not, why. ng_ridge_cv_predict() returns NULL when n < 10 (R/02_effects.R:164),
#    so cv_predictive_r2 becomes NA, so the GEBV branch of ng_choose_mean_source()
#    cannot fire, so the phenotype mid-parent is GUARANTEED regardless of marker
#    quality. Today that is indistinguishable from a model that was measured and
#    failed. The distinction matters because only the first is fixable by supplying
#    a training set.
#
# The empty-data-frame JSON round trip is asserted because a zero-row data.frame
# serialises to {} rather than [] under some settings, which would break a
# consumer on precisely the happy path where there is nothing to report.
suppressPackageStartupMessages(library(jsonlite))
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

ADVISORY_COLS <- c("id", "severity", "stage", "trait", "table", "field", "message", "count")

# ---- 1. constructors -------------------------------------------------------
empty <- ng_advisory_empty()
stopifnot(is.data.frame(empty), nrow(empty) == 0L)
stopifnot(identical(names(empty), ADVISORY_COLS))

one <- ng_advisory(id = "marker_effects_weak", severity = "warning", stage = "predict",
                   message = "YIELD: cv_predictive_r2 = -0.002", trait = "YIELD")
stopifnot(is.data.frame(one), nrow(one) == 1L)
stopifnot(identical(names(one), ADVISORY_COLS))
stopifnot(identical(one$id, "marker_effects_weak"), identical(one$trait, "YIELD"))
stopifnot(identical(one$count, 1L))

acc <- ng_advisory_add(empty, id = "training_set_would_help", severity = "advice",
                       stage = "predict", message = "supply training_genotype/training_phenotype")
acc <- rbind(acc, one)
stopifnot(nrow(acc) == 2L, identical(names(acc), ADVISORY_COLS))

# Severity domain is closed. A typo must not silently become a new severity.
stopifnot(all(acc$severity %in% c("advice", "warning", "blocker")))
stopifnot(inherits(tryCatch(ng_advisory(id = "x", severity = "urgent", stage = "predict",
                                        message = "m"), error = function(e) e), "error"))

# ---- 2. roll-up ------------------------------------------------------------
roll <- ng_advisory_rollup(acc)
stopifnot(identical(roll$schema, "ng_advisories.v1"))
stopifnot(identical(roll$counts$advice, 1L), identical(roll$counts$warnings, 1L))
stopifnot(identical(roll$counts$blockers, 0L), identical(roll$counts$issues, 2L))
stopifnot(identical(roll$status, "warning"))
stopifnot(identical(ng_advisory_rollup(empty)$status, "pass"))
stopifnot(identical(ng_advisory_rollup(
  ng_advisory_add(empty, id = "b", severity = "blocker", stage = "predict",
                  message = "m"))$status, "blocker"))

# ---- 3. JSON round trip, including the empty case --------------------------
rt <- function(d) fromJSON(toJSON(d, dataframe = "rows", na = "null"))
stopifnot(nrow(rt(acc)) == 2L)
stopifnot(identical(unname(unlist(toJSON(empty, dataframe = "rows"))), character(0)) ||
          identical(as.character(toJSON(empty, dataframe = "rows")), "[]"))

# ---- 4. emission happens, and only for warning/blocker ---------------------
# Advice is recorded but must not raise -- a run full of advice is not a run in
# trouble, and warning-fatigue is how real warnings get ignored.
w <- character(0)
withCallingHandlers(ng_emit_advisories(acc),
                    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
stopifnot(length(w) == 1L)
stopifnot(grepl("marker_effects_weak", w[[1L]], fixed = TRUE))

w2 <- character(0)
withCallingHandlers(ng_emit_advisories(empty),
                    warning = function(x) { w2 <<- c(w2, conditionMessage(x)); invokeRestart("muffleWarning") })
stopifnot(length(w2) == 0L)

# ---- 5. the fit reports whether CV could run, and why not ------------------
set.seed(7)
mk_fit <- function(n) {
  m <- 24L
  g <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
  rownames(g) <- paste0("P", seq_len(n)); colnames(g) <- paste0("M", seq_len(m))
  y <- setNames(as.numeric(g %*% c(rnorm(6), rep(0, m - 6))) + rnorm(n, sd = 0.2), rownames(g))
  ng_fit_ridge_effects(g, y, ids = rownames(g))
}

big <- mk_fit(30L)
stopifnot(isTRUE(big$cv_available))
stopifnot(is.na(big$cv_unavailable_reason))
stopifnot(identical(big$n_train, 30L))
stopifnot(is.finite(big$cv_predictive_r2))

# n < 10: CV cannot run, so cv_predictive_r2 is NA -- and the REASON is recorded,
# because this is the case a training set fixes.
small <- mk_fit(8L)
stopifnot(isFALSE(small$cv_available))
stopifnot(identical(small$cv_unavailable_reason, "n_lt_10"))
stopifnot(identical(small$n_train, 8L))
stopifnot(is.na(small$cv_predictive_r2))

cat("advisory_schema_contract: PASS\n")
