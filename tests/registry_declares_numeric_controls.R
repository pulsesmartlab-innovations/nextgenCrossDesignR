# The capability registry must be able to describe a NUMBER.
#
# The registry is the single source of truth a frontend renders from -- helpers.R in the
# workbench says so, and drops any choice the backend marks experimental/guarded "so the
# backend can retract a capability without a frontend edit". That governance worked, but
# only for enums: there was one constructor, enum(), and every one of the 22 controls
# was type "enum". A numeric parameter could not be expressed at all.
#
# So every numeric knob was invisible to the contract, and the frontend had no choice but
# to hardcode them. The consequence was the worst available case:
#
#   min_effect_reliability   hardcoded in the UI as a numericInput(0.35)  -- and INERT.
#                            It gates the calibrated-reliability branch that nothing in
#                            this package can reach, so the dial governs no run.
#   min_cv_predictive_r2     absent from the frontend entirely -- and it is the threshold
#                            that decides whether cross means are genomic at all.
#
# A UI cannot be blamed for surfacing the wrong dial when the contract could not name the
# right one. This makes numeric capability expressible, declares the two that govern the
# reliability gate, and marks the inert one so a frontend can retract it.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

ctls <- ng_backend_controls()
by_id <- function(id) { h <- Filter(function(c) identical(c$id, id), ctls)
                        if (length(h)) h[[1L]] else NULL }
f <- formals(ng_run_cross_prediction)
eng_default <- function(nm, is_enum) {
  d <- f[[nm]]; v <- if (is.call(d)) eval(d) else d
  if (is_enum && length(v) > 1L) v[[1L]] else v
}

# ---- 1. the registry can express a number at all --------------------------
types <- unique(vapply(ctls, function(c) c$type %||% "?", character(1)))
if (!"number" %in% types) {
  stop("the registry still emits only ", paste(types, collapse = "/"),
       " -- no numeric parameter can be declared, so a frontend must hardcode every one",
       call. = FALSE)
}

# ---- 2. the knob that governs mean selection is declared ------------------
m <- by_id("min_cv_predictive_r2")
if (is.null(m)) stop("min_cv_predictive_r2 is not in the registry, so a frontend ",
                     "cannot discover the threshold that decides GEBV vs phenotypic ",
                     "cross means", call. = FALSE)
stopifnot(identical(m$type, "number"))
stopifnot(isTRUE(all.equal(m$default, eng_default("min_cv_predictive_r2", FALSE))))
stopifnot(is.finite(m$min), is.finite(m$max), m$min < m$max)

# ---- 3. the gate switch is declared, with the engine's own choices --------
g <- by_id("effect_gate")
if (is.null(g)) stop("effect_gate is not in the registry", call. = FALSE)
stopifnot(identical(g$type, "enum"))
vals <- vapply(g$choices, function(x) x$value, character(1))
stopifnot(setequal(vals, c("on", "off")))
stopifnot(identical(g$default, "on"))

# ---- 4. the INERT knob is declared AND marked, not silently omitted ------
# Omitting it would leave the frontend's hardcoded copy in place with nothing to say
# it governs nothing. Declaring it guarded is what lets a UI retract it.
r <- by_id("min_effect_reliability")
if (is.null(r)) stop("min_effect_reliability is absent: the frontend hardcodes this ",
                     "dial today, and nothing in the contract says it is inert",
                     call. = FALSE)
stopifnot(identical(r$status, "guarded"))
stopifnot(!is.null(r$note), nzchar(r$note))
stopifnot(grepl("reserved|inert|no effect|governs no", r$note, ignore.case = TRUE))

# ---- 5. every declared default RESOLVES to the engine's default ----------
# The registry speaks the breeder vocabulary ("reliable_family_variance") and the engine
# speaks tokens ("pmv"); ng_normalize_metric_token() is the bridge, so the invariant is a
# round-trip, not string equality. That is the stronger claim: it proves the value the
# contract advertises is one the engine actually accepts, in either vocabulary.
for (c in ctls) {
  if (!(c$id %in% names(f))) next
  is_enum <- identical(c$type, "enum")
  ed <- eng_default(c$id, is_enum)
  if (is.null(ed) || is.symbol(ed)) next
  same <- isTRUE(all.equal(as.character(c$default), as.character(ed))) ||
          isTRUE(all.equal(as.character(ng_normalize_metric_token(as.character(c$default))),
                           as.character(ed)))
  if (!same) {
    stop("registry says ", c$id, " defaults to '", c$default,
         "' which does not resolve to the engine's '", ed, "'", call. = FALSE)
  }
}

# ---- 5b. and every enum CHOICE the registry offers is one the engine takes -
# A contract may not advertise a value the backend refuses.
for (c in ctls) {
  if (!identical(c$type, "enum") || !(c$id %in% names(f))) next
  allowed <- eng_default(c$id, FALSE)
  allowed <- if (is.call(f[[c$id]]) || length(eval_allowed <- f[[c$id]]) <= 1L) NULL
             else as.character(eval_allowed)
  if (is.null(allowed)) next
  for (ch_i in c$choices) {
    v <- ch_i$value
    if (!(v %in% allowed) && !(ng_normalize_metric_token(v) %in% allowed)) {
      stop("registry offers ", c$id, " = '", v, "', which the engine does not accept",
           call. = FALSE)
    }
  }
}

# ---- 6. it survives the JSON bridge, which is what the frontend reads ----
# jsonlite silently re-types across the bridge; a numeric control that arrives as a
# string or an array is not renderable.
tmp <- tempfile(fileext = ".json")
invisible(ng_write_backend_capability_registry_json(output_path = tmp))
j <- jsonlite::fromJSON(tmp, simplifyVector = FALSE)
stopifnot(length(j$controls) == length(ctls))
jm <- Filter(function(c) identical(c$id, "min_cv_predictive_r2"), j$controls)[[1L]]
stopifnot(identical(jm$type, "number"))
stopifnot(is.numeric(jm$default), length(jm$default) == 1L)
stopifnot(is.numeric(jm$min), is.numeric(jm$max))

# ---- 7. the live threshold is described on the scale breeders read ------
# cv_predictive_r2 is an out-of-fold R2 against PHENOTYPE. Breeders reason in accuracy,
# and 0.35 R2 is r ~ 0.59 -- a demanding bar, not a modest one. A contract that states
# the number without the scale hands the breeder the acceptability decision and invites
# them to make it backwards. Being against phenotype rather than breeding value, it is
# also attenuated by heritability, so it UNDERSTATES GEBV accuracy: a breeder with
# h2 = 0.3 seeing R2 = 0.20 may reject markers that predict breeding value well.
mnote <- paste(unlist(m$note), collapse = " ")
stopifnot(nzchar(mnote))
stopifnot(grepl("0.59", mnote, fixed = TRUE))          # the correlation equivalent
stopifnot(grepl("phenotype", mnote, ignore.case = TRUE))
stopifnot(grepl("breeding value", mnote, ignore.case = TRUE))
# the threshold gates the MEAN; the variance stays marker-derived in every tier
stopifnot(grepl("variance", mnote, ignore.case = TRUE))
# one scalar covers every trait, and acceptability genuinely differs by trait
stopifnot(grepl("every trait|all traits|one threshold", mnote, ignore.case = TRUE))

# ---- 8. and it does NOT inherit the word "reliability" -------------------
# In quantitative genetics reliability is r^2(GEBV, TBV). This is not that, and the
# package reports the calibrated quantity as missing precisely because a phenotype CV
# statistic is not it. Two adjacent controls, one called reliability, would be read as
# two settings for one thing.
stopifnot(!grepl("reliability", as.character(m$label), ignore.case = TRUE))

# ---- 9. the guarded knob says what it is NOT -----------------------------
stopifnot(grepl("min_cv_predictive_r2", paste(unlist(r$note), collapse = " "), fixed = TRUE))

# ---- 10. the override reads as a decision, not as a diagnostics mode -----
# A breeder switching the gate off for a real run is not doing diagnostics; a label
# saying so misdescribes the choice they are making.
off <- Filter(function(x) identical(x$value, "off"), g$choices)[[1L]]
stopifnot(!grepl("diagnostics only", as.character(off$label), ignore.case = TRUE))
stopifnot(grepl("proceed|without the gate|override", as.character(off$label), ignore.case = TRUE))

cat("registry_declares_numeric_controls: PASS  (", length(ctls), "controls,",
    paste(types, collapse = "/"), ")\n")
