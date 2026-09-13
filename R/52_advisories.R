# Post-effects advisories: things the run learns only after marker effects exist.
#
# WHY A SECOND CHANNEL. `qc$issues` (R/26_data_preflight.R:5-31) is the right SHAPE
# -- a severity-graded, per-row, renderable table -- but the wrong LIFETIME. It is
# finalised in the qc stage (R/26:558-573), long before any marker effect is fitted,
# so it structurally cannot carry "this trait's markers do not predict out of
# sample", "cross-validation could not run at all", or "a training set would fix
# this". Those are the facts a breeder most needs and the ones that arrive latest.
#
# The columns deliberately MATCH `qc$issues` (plus `stage` and `trait`) so that any
# consumer able to render QC issues renders advisories with the same code.
#
# SEVERITY is one level richer than qc's blocker/warning:
#
#   advice   recorded, never raised. Most of this channel is genuinely advice
#            ("supplying a training set would raise cv_predictive_r2"), and calling
#            advice a warning is how real warnings come to be ignored.
#   warning  recorded AND raised as an R warning.
#   blocker  recorded AND raised; the caller decides whether to stop.
#
# WHERE THESE MAY BE EMITTED. Only from the PARENT process. `run_trait_job` executes
# under mclapply/parLapply (R/39:435-470), and a warning() raised in a forked child
# or a PSOCK worker never reaches the parent's calling handlers -- so it would never
# be collected into the JSON envelope's `warnings` array
# (tools/run_cross_prediction_json.R:110-114) and the advisory would be silent. A
# silent advisory channel is the exact defect class this work exists to remove.
# Build advisories from per-trait results AFTER they are collected, then call
# ng_emit_advisories() once, in the parent.
#
# The upside of that constraint: a parent-process warning() already reaches the
# frontend today through the existing `warnings` array, so the human-readable text
# is delivered with no frontend change at all.

ng_advisory_severities <- function() c("advice", "warning", "blocker")

ng_advisory_empty <- function() {
  data.frame(id = character(), severity = character(), stage = character(),
             trait = character(), table = character(), field = character(),
             message = character(), count = integer(),
             stringsAsFactors = FALSE)
}

ng_advisory <- function(id, severity, stage, message,
                        trait = NA_character_, table = NA_character_,
                        field = NA_character_, count = 1L) {
  sev <- as.character(severity)[[1L]]
  if (!sev %in% ng_advisory_severities()) {
    ng_stop("advisory severity must be one of: ",
            paste(ng_advisory_severities(), collapse = ", "), " (got '", sev, "')")
  }
  data.frame(id = as.character(id)[[1L]], severity = sev,
             stage = as.character(stage)[[1L]],
             trait = as.character(trait)[[1L]], table = as.character(table)[[1L]],
             field = as.character(field)[[1L]], message = as.character(message)[[1L]],
             count = as.integer(count)[[1L]], stringsAsFactors = FALSE)
}

ng_advisory_add <- function(advisories, id, severity, stage, message,
                            trait = NA_character_, table = NA_character_,
                            field = NA_character_, count = 1L) {
  rbind(advisories,
        ng_advisory(id = id, severity = severity, stage = stage, message = message,
                    trait = trait, table = table, field = field, count = count))
}

# Mirrors ng_preflight roll-up (R/26:558-573) so the two summaries read alike.
# NOTE the deliberate separation: this status never feeds the QC blocker gate. An
# advisory must not be able to flip a run's QC verdict -- the two channels answer
# different questions ("is the input usable?" vs "is the model worth trusting?").
ng_advisory_rollup <- function(advisories) {
  advisories <- if (is.null(advisories)) ng_advisory_empty() else advisories
  sev <- as.character(advisories$severity)
  list(
    schema = "ng_advisories.v1",
    counts = list(advice   = sum(sev == "advice"),
                  warnings = sum(sev == "warning"),
                  blockers = sum(sev == "blocker"),
                  issues   = nrow(advisories)),
    status = if (any(sev == "blocker")) "blocker"
             else if (any(sev == "warning")) "warning" else "pass"
  )
}

# Raise the non-advice rows as real R warnings. PARENT PROCESS ONLY -- see the
# header. Returns its input invisibly so it can be used inline.
ng_emit_advisories <- function(advisories) {
  advisories <- if (is.null(advisories)) ng_advisory_empty() else advisories
  loud <- which(as.character(advisories$severity) %in% c("warning", "blocker"))
  for (i in loud) {
    warning(sprintf("[%s] %s", advisories$id[[i]], advisories$message[[i]]),
            call. = FALSE)
  }
  invisible(advisories)
}
