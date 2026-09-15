# The frontend contract must describe the engine, not a wish.
#
# tests/contract_schema_drift.R checks config_schema.json's parameter NAMES
# against formals(). It never checked `allowed` or `default`, and three defects
# lived in exactly that blind spot:
#
#   * trait_value_metric.default was "var_complex"; the real default is
#     "usefulness" -- so a UI prefilling from the contract mislabels the run;
#   * uc_variance_source.allowed offered "le", which ng_cp__build_ctx() hard-errors
#     on. The contract documented a config the backend refuses, so a UI could only
#     learn that by shipping it to a user;
#   * min_cv_predictive_r2 -- the knob that actually governs whether cross means
#     are genomic -- was absent entirely, while min_effect_reliability was present.
#     min_effect_reliability gates a branch nothing in this package can reach
#     (nothing sets reliability_is_calibrated = TRUE; see R/02_effects.R). The UI
#     could therefore expose only the threshold that governs nothing.
#
# A contract is a promise to someone who cannot read the source. Getting it wrong
# is not a documentation nit; it is the same defect class as an unlabelled number.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

schema <- jsonlite::fromJSON(file.path(root, "docs", "frontend", "contracts",
                                       "config_schema.json"), simplifyVector = FALSE)
params <- unlist(lapply(schema$groups, function(g) g$params), recursive = FALSE)
names(params) <- vapply(params, function(p) p$name, character(1))
f <- formals(ng_run_cross_prediction)

# The match.arg idiom -- a formal whose default IS its choice vector -- is the only
# place formals() carries an allowed set. A formal defaulting to a CALL (effect_gate
# = getOption(...), so a user can set it globally) carries only the resolved value,
# not the choices, so it must be validated by round-trip instead. Treating the
# resolved value as the allowed set would wrongly reject "off".
engine_allowed <- function(nm) {
  d <- f[[nm]]
  if (!is.call(d) && length(d) > 1L) as.character(d) else NULL
}
# match.arg's "first element is the default" rule applies to ENUMS only. A numeric
# formal may legitimately default to a vector (priority_breaks = c(0.1, .35, .7, 1)),
# where the whole vector is the default and taking element 1 would be nonsense.
engine_default <- function(nm, is_enum) {
  d <- f[[nm]]
  v <- if (is.call(d)) eval(d) else d
  if (is_enum && length(v) > 1L) v[[1L]] else v
}

# ---- 1. no `allowed` value the engine rejects -----------------------------
for (nm in names(params)) {
  p <- params[[nm]]
  if (is.null(p$allowed) || !(nm %in% names(f))) next
  eng <- engine_allowed(nm)
  if (is.null(eng)) next
  bad <- setdiff(unlist(p$allowed), eng)
  if (length(bad)) stop("config_schema.json offers ", nm, " value(s) the engine rejects: ",
                        paste(bad, collapse = ", "), call. = FALSE)
}

# ---- 2. every stated `default` is the engine's default --------------------
for (nm in names(params)) {
  p <- params[[nm]]
  if (is.null(p$default) || !(nm %in% names(f))) next
  eng <- engine_default(nm, identical(p$type, "enum"))
  if (is.null(eng) || is.symbol(eng)) next
  if (!isTRUE(all.equal(as.character(p$default), as.character(eng)))) {
    stop("config_schema.json says ", nm, " defaults to '", p$default,
         "'; the engine defaults to '", eng, "'", call. = FALSE)
  }
}

# ---- 3. usefulness + parent_distance is not offered as a choice -----------
# It is refused at config time, so a contract that lists it is documenting a
# config the backend will not run.
ucs <- unlist(params[["uc_variance_source"]]$allowed)
if (any(c("le", "parent_distance") %in% ucs)) {
  stop("config_schema.json still offers uc_variance_source = parent_distance/le, ",
       "which ng_cp__build_ctx() refuses when trait_value_metric = 'usefulness'",
       call. = FALSE)
}

# ---- 1b. and no engine capability MISSING from the contract ---------------
# The reverse direction, and the one that was absent. Check 1 fails only when the schema
# offers something the engine rejects; nothing failed when the ENGINE gained a value the
# schema never learned about. That is how multi_trait_method = "rank_sum" shipped in the
# backend while the contract still listed five methods: a frontend that renders from the
# contract simply could not offer it, and no test said so.
#
# The authority differs per parameter. Most enums carry their choices in the formal's
# default vector, but multi_trait_method's formal is just "auto" -- its real authority is
# ng_multitrait_methods(). Checking formals alone would miss exactly the case that
# motivated this.
engine_choices <- function(nm) {
  if (identical(nm, "multi_trait_method")) return(ng_multitrait_methods())
  d <- f[[nm]]
  if (!is.call(d) && length(d) > 1L) as.character(d) else NULL
}

# Deliberate narrowings, each with the reason it is not a gap. Same idea as
# `undocumented_ok` in tests/contract_schema_drift.R: an exception must be named and
# justified, never inferred from a structural quirk.
narrowed_on_purpose <- list(
  uc_variance_source = paste(
    "the formal accepts parent_distance/le, but ng_cp__build_ctx() hard-errors on them",
    "in the only context where this parameter applies (trait_value_metric =",
    "'usefulness'), so the contract must not advertise a run the backend refuses")
)

for (nm in names(params)) {
  p <- params[[nm]]
  if (is.null(p$allowed) || !(nm %in% names(f))) next
  eng <- engine_choices(nm)
  if (is.null(eng)) next
  missing_from_schema <- setdiff(eng, unlist(p$allowed))
  # a value the schema states in the breeder vocabulary still counts as declared
  if (length(missing_from_schema)) {
    declared <- vapply(unlist(p$allowed),
                       function(v) as.character(ng_normalize_metric_token(v)), character(1))
    missing_from_schema <- setdiff(missing_from_schema, declared)
  }
  if (length(missing_from_schema) && !is.null(narrowed_on_purpose[[nm]])) next
  if (length(missing_from_schema)) {
    stop("the engine accepts ", nm, " value(s) that config_schema.json never lists: ",
         paste(missing_from_schema, collapse = ", "),
         " -- a frontend rendering from the contract cannot offer them. Add them, or ",
         "name the narrowing in narrowed_on_purpose with its reason.", call. = FALSE)
  }
}

# ---- 3b. enums whose choices formals() cannot carry, validated by round-trip --
# ng_cp__build_ctx() match.arg()s these, so an accepted value survives and a
# rejected one errors. This is what keeps the contract honest for effect_gate,
# whose formal default is a getOption() call rather than a choice vector.
base <- lapply(f, function(x) if (is.call(x) || is.name(x)) tryCatch(eval(x), error = function(e) NULL) else x)
base <- base[!vapply(base, function(x) identical(x, quote(expr = )), logical(1))]
set.seed(2); nn <- 12L; mm <- 16L
gid <- sprintf("P%02d", seq_len(nn))
gm <- matrix(sample(c(0L, 2L), nn * mm, TRUE), nrow = nn,
             dimnames = list(gid, sprintf("M%02d", seq_len(mm))))
base$genotype <- data.frame(NAME = gid, gm, check.names = FALSE)
base$genotype_id_col <- "NAME"
base$phenotype <- data.frame(NAME = gid, YIELD = as.numeric(gm %*% rnorm(mm)))
base$phenotype_id_col <- "NAME"; base$traits_to_use <- "YIELD"
base$trait_direction <- data.frame(Trait = "YIELD", Selection_direction = "increase")
base$direction_trait_col <- "Trait"; base$direction_column_col <- "Trait"
base$direction_direction_col <- "Selection_direction"

for (nm in names(params)) {
  p <- params[[nm]]
  if (is.null(p$allowed) || !(nm %in% names(f))) next
  if (!is.null(engine_allowed(nm))) next          # already covered by check 1
  if (!is.call(f[[nm]])) next                     # not the getOption idiom
  for (v in unlist(p$allowed)) {
    cfg <- base; cfg[[nm]] <- v
    e <- tryCatch({ ng_cp__build_ctx(cfg); NULL }, error = function(e) e)
    if (!is.null(e) && grepl("should be one of|'arg' should", conditionMessage(e))) {
      stop("config_schema.json offers ", nm, " = '", v, "', which the engine rejects: ",
           conditionMessage(e), call. = FALSE)
    }
  }
}

# ---- 4. the knob that governs mean selection is discoverable --------------
# A UI that cannot find min_cv_predictive_r2 cannot expose the decision that
# determines whether a delivered cross mean is genomic or phenotypic.
for (nm in c("min_cv_predictive_r2", "effect_gate")) {
  if (!(nm %in% names(params))) {
    stop("config_schema.json omits ", nm, ", so the frontend cannot surface it",
         call. = FALSE)
  }
}

# ---- 5. and a knob that governs nothing is labelled as such ---------------
mer <- params[["min_effect_reliability"]]
if (!is.null(mer)) {
  d <- paste(unlist(mer$desc), collapse = " ")
  if (!grepl("reserved|not active|no effect|inert", d, ignore.case = TRUE)) {
    stop("config_schema.json exposes min_effect_reliability without saying it is ",
         "reserved; it gates a branch nothing in this package can reach", call. = FALSE)
  }
}

cat("contract_declares_what_the_engine_does: PASS\n")
