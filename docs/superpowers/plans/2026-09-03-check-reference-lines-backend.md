# Check Reference Lines — Backend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make check lines reference-only benchmarks — a check genotype supplies a per-trait
reference value that is reported and plotted, and never filters, penalises, or reorders a
single cross.

**Architecture:** Check genotypes arrive in a matrix separate from the candidate parents. They
are aligned to the marker set already fixed by parent QC, scored with the same per-trait marker
effects, and resolved onto the same `mean_source` the cross means use. The resulting scalar per
trait is attached as reference columns and drawn as a horizontal line. Checks never enter QC,
LD, GRM, or `ng_make_pairs()`, which is what makes the "reference only" claim testable.

**Tech Stack:** R (base + Matrix), openxlsx (workbook), grDevices/graphics (plots), plain-R
`stopifnot` test scripts under `tests/`, testthat under `tests/testthat/`.

**Spec:** `docs/design/2026-09-03-check-reference-lines-design.md`

## Global Constraints

- Package: `nextgenCrossDesign`, version `0.22.0` → **`0.23.0`** (breaking).
- Working branch: `feat/check-reference-lines` (already created; the spec commit is `b837878`).
- **No `Co-Authored-By: Claude` trailer on any commit in this project.** The user is the sole
  contributor of record.
- OS-agnostic: no shell-outs, no platform-specific paths, `file.path()` everywhere.
- **Removed public API (breaking):** `check_basis`, `exclude_threshold_violators`, and the
  parent-only check validation at `R/39_cross_prediction_runner.R:1215`.
- **The invariant:** a run with checks must produce cross predictions byte-identical to the
  same run without checks. Task 7 enforces it; every other task must not break it.
- Ties are never violations (`>=`, not `>`).
- **`%||%` is NOT available.** It is defined only as a local inside one function body
  (`R/39_cross_prediction_runner.R:1568`) and is neither exported nor imported by the package.
  Use explicit `if (is.null(x)) ... else ...`; do not define a local copy in a new file.
- **`ng_predict_gebv(geno, effects)` reads `effects$beta`** (named numeric over marker ids)
  **and `effects$intercept`** (`R/02_effects.R:188`). There is no `$effect` field anywhere in
  this package — any test fixture must use `beta`/`intercept`.
- Direction vocabulary: `reject_if = "below"` means the breeder wants the mid-parent **above**
  the check (an `increase` trait). `reject_if = "above"` means they want it **below**
  (a `decrease` trait).
- Test style: plain R scripts under `tests/`, sourcing `tests/helper_load.R`, asserting with
  `stopifnot()`. Run one with `Rscript tests/<name>.R` from the repo root.

---

### Task 1: Align a check genotype matrix to the parents' marker set

A check matrix from a separate file will not automatically share the parents' marker set,
column order, or allele coding. `ng_predict_gebv()` is a dot product, so a misaligned matrix
yields a wrong-but-plausible reference value with no warning. This task builds the gate.

**Files:**
- Create: `R/51_check_reference.R`
- Test: `tests/check_reference.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `ng_align_check_geno(check_geno, marker_names, ploidy = 2L)` → numeric matrix with
  columns exactly `marker_names` in that order, rownames preserved. Errors on any missing
  marker.

- [ ] **Step 1: Write the failing test**

Create `tests/check_reference.R`:

```r
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- Task 1: check-genotype alignment ---------------------------------------
markers <- c("m1", "m2", "m3")
chk <- matrix(c(0, 2, 2,
                2, 0, 0), nrow = 2, byrow = TRUE,
              dimnames = list(c("CHK_A", "CHK_B"), c("m3", "m1", "m2")))

# columns are reordered to the parents' marker order, rownames preserved
al <- ng_align_check_geno(chk, markers)
stopifnot(identical(colnames(al), markers))
stopifnot(identical(rownames(al), c("CHK_A", "CHK_B")))
stopifnot(al["CHK_A", "m1"] == 2, al["CHK_A", "m3"] == 0)

# extra markers in the check file are dropped, not an error
chk_extra <- cbind(chk, m9 = c(1, 1))
al2 <- ng_align_check_geno(chk_extra, markers)
stopifnot(identical(colnames(al2), markers))

# a MISSING marker is a hard error naming the count
err <- tryCatch(ng_align_check_geno(chk[, c("m1", "m2"), drop = FALSE], markers),
                error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("m3", err), grepl("1 of 3", err))

# unnamed columns are a hard error (cannot verify alignment)
bare <- matrix(0, nrow = 1, ncol = 3, dimnames = list("CHK_A", NULL))
err2 <- tryCatch(ng_align_check_geno(bare, markers), error = function(e) conditionMessage(e))
stopifnot(is.character(err2), grepl("column names", err2))

# dosages outside [0, ploidy] are a hard error (wrong coding / flipped file)
badcode <- matrix(c(0, 5, 1), nrow = 1, dimnames = list("CHK_A", markers))
err3 <- tryCatch(ng_align_check_geno(badcode, markers), error = function(e) conditionMessage(e))
stopifnot(is.character(err3), grepl("dosage", err3))

cat("task 1 ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference.R`
Expected: FAIL with `could not find function "ng_align_check_geno"`.

- [ ] **Step 3: Write minimal implementation**

Create `R/51_check_reference.R`:

```r
# Check lines as REFERENCES, not filters. A check is a benchmark genotype (a released variety,
# a commercial check) that is never crossed: check x check would be a self, which
# ng_make_pairs() does not enumerate, so a check has no variance and can never be a row in the
# cross table. It contributes a scalar per trait, on the same mean_source the cross means use.
# See docs/design/2026-09-03-check-reference-lines-design.md.

# Align a separate check genotype matrix to the marker set already fixed by parent QC. A check
# matrix that silently disagrees on marker set or column order produces a wrong-but-plausible
# reference value, so every disagreement is an error rather than an intersection.
ng_align_check_geno <- function(check_geno, marker_names, ploidy = 2L) {
  check_geno <- ng_as_numeric_matrix(check_geno, "check_geno")
  marker_names <- as.character(marker_names)
  if (is.null(colnames(check_geno))) {
    ng_stop("check_geno must have marker column names to verify alignment with the parent ",
            "marker set; an unnamed matrix cannot be checked for allele-order agreement")
  }
  if (is.null(rownames(check_geno))) ng_stop("check_geno must have check ids as rownames")
  miss <- setdiff(marker_names, colnames(check_geno))
  if (length(miss)) {
    ng_stop(sprintf(
      "check_geno is missing %d of %d markers used by the fitted effects (e.g. %s). The check ",
      length(miss), length(marker_names),
      paste(utils::head(miss, 5L), collapse = ", ")),
      "file must carry the same marker set, coding, and reference allele as the parent genotypes.")
  }
  out <- check_geno[, marker_names, drop = FALSE]
  ploidy <- as.integer(ploidy[[1L]])
  rng <- suppressWarnings(range(out, na.rm = TRUE))
  if (any(is.finite(rng)) && (min(rng) < 0 || max(rng) > ploidy)) {
    ng_stop(sprintf(
      "check_geno dosage values fall outside [0, %d] (observed %g..%g); the check file must use ",
      ploidy, rng[[1L]], rng[[2L]]),
      "the same allele coding as the parent genotypes")
  }
  out
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript tests/check_reference.R`
Expected: `task 1 ok`

- [ ] **Step 5: Commit**

```bash
git add R/51_check_reference.R tests/check_reference.R
git commit -m "feat(checks): align a separate check genotype matrix to the parent marker set"
```

---

### Task 2: Resolve the check's value onto the run's mean source

`ng_choose_mean_source()` decides **per trait** whether cross means are GEBV or a phenotypic
fallback. The check must land on that same source or the reference line is on a different scale
from the axis it is drawn on.

**Files:**
- Modify: `R/51_check_reference.R`
- Test: `tests/check_reference.R`

**Interfaces:**
- Consumes: `ng_align_check_geno()` from Task 1.
- **`ng_predict_gebv(geno, effects)` reads `effects$beta` (a named numeric over marker ids) and
  `effects$intercept`** (`R/02_effects.R:188`). Test fixtures must use those field names — there
  is no `$effect` field anywhere in this package.
- **`%||%` is NOT a package helper.** It exists only as a local inside one function body
  (`R/39_cross_prediction_runner.R:1568`) and is neither exported nor imported. Use explicit
  `if (is.null(x))` checks; do not define a local copy.
- Produces: `ng_check_reference_value(source, check_geno_aligned, effects, check_records = NULL)`
  → named numeric vector over `rownames(check_geno_aligned)`. Returns `NA_real_` for any check
  with no record on a non-GEBV source. `check_records` is a named list with any of
  `BLUP`, `BLUE`, `adjusted_pheno`, each a named numeric vector keyed by check id.

- [ ] **Step 1: Write the failing test**

Append to `tests/check_reference.R` (before the final `cat`, then move the `cat` to the end —
simplest is to append this block and add a new `cat("task 2 ok\n")` at the end of it):

```r
# --- Task 2: value resolution follows the run's mean source -----------------
eff <- list(beta = c(m1 = 1, m2 = 0.5, m3 = -2), intercept = 0)
al <- ng_align_check_geno(chk, markers)   # CHK_A = (m1=2, m2=2, m3=0); CHK_B = (0, 0, 2)

# GEBV source -> predicted from the check's own markers with the SAME effects
gv <- ng_check_reference_value("GEBV", al, eff)
stopifnot(is.numeric(gv), identical(names(gv), c("CHK_A", "CHK_B")))
stopifnot(abs(gv[["CHK_A"]] - (2 * 1 + 2 * 0.5 + 0 * -2)) < 1e-8)
stopifnot(abs(gv[["CHK_B"]] - (0 * 1 + 0 * 0.5 + 2 * -2)) < 1e-8)

# the low-reliability / uncalibrated GEBV variants are still GEBV
stopifnot(abs(ng_check_reference_value("GEBV_uncalibrated", al, eff)[["CHK_A"]] - gv[["CHK_A"]]) < 1e-8)
stopifnot(abs(ng_check_reference_value("GEBV_low_reliability", al, eff)[["CHK_A"]] - gv[["CHK_A"]]) < 1e-8)

# a phenotypic source reads the supplied record, NOT the markers
bv <- ng_check_reference_value("BLUE", al, eff,
                               check_records = list(BLUE = c(CHK_A = 11.5, CHK_B = 9.25)))
stopifnot(abs(bv[["CHK_A"]] - 11.5) < 1e-8, abs(bv[["CHK_B"]] - 9.25) < 1e-8)

# NOT-EVALUABLE: run source is BLUE, the check has no BLUE record -> NA, never a GEBV fallback
nv <- ng_check_reference_value("BLUE", al, eff, check_records = list(BLUE = c(CHK_A = 11.5)))
stopifnot(abs(nv[["CHK_A"]] - 11.5) < 1e-8, is.na(nv[["CHK_B"]]))

# no records at all on a phenotypic source -> all NA (still not a GEBV fallback)
allna <- ng_check_reference_value("adjusted_pheno", al, eff)
stopifnot(all(is.na(allna)))

cat("task 2 ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference.R`
Expected: FAIL with `could not find function "ng_check_reference_value"`.

- [ ] **Step 3: Write minimal implementation**

Append to `R/51_check_reference.R`:

```r
# Resolve each check's value onto the SAME source the cross means used. ng_choose_mean_source()
# returns one of "GEBV", "GEBV_low_reliability", "GEBV_uncalibrated", "BLUP", "BLUE",
# "adjusted_pheno" -- per trait. A GEBV* source predicts from the check's markers with the
# trait's fitted effects; a phenotypic source reads the check's own record. A check with no
# record on a phenotypic source is NA (not evaluable) and must NEVER silently fall back to a
# GEBV, which would put the reference on a different scale from the axis it is drawn on.
ng_check_reference_value <- function(source, check_geno_aligned, effects, check_records = NULL) {
  source <- as.character(source)[[1L]]
  ids <- rownames(check_geno_aligned)
  if (startsWith(source, "GEBV")) {
    return(stats::setNames(ng_predict_gebv(check_geno_aligned, effects), ids))
  }
  rec <- if (is.null(check_records)) NULL else check_records[[source]]
  if (is.null(rec)) return(stats::setNames(rep(NA_real_, length(ids)), ids))
  v <- suppressWarnings(as.numeric(rec[ids]))
  stats::setNames(v, ids)
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript tests/check_reference.R`
Expected: `task 1 ok` then `task 2 ok`

- [ ] **Step 5: Commit**

```bash
git add R/51_check_reference.R tests/check_reference.R
git commit -m "feat(checks): resolve check values onto the run's per-trait mean source"
```

---

### Task 3: Attach reference columns (never subset rows)

This replaces the veto. `ng_apply_trait_checks()` flagged **and dropped**; the replacement only
ever adds columns.

**Files:**
- Modify: `R/51_check_reference.R`
- Test: `tests/check_reference.R`

**Interfaces:**
- Consumes: nothing from Tasks 1-2 at call time (values are passed in).
- Produces: `ng_attach_check_reference(scores, spec, trait_values, check_values, k_progeny, mean_suffix = "_mean", sd_suffix = "_pmv_used")`.
  **`k_progeny` has NO default and is required** — progeny per family is a breeding-program
  fact the user supplies, never a number this package invents. It scales P(beat check)
  directly, so a silent default would drive a reported probability from a made-up input.
  → the same data frame with `nrow()` unchanged, plus per trait `<t>_check_id`,
  `<t>_check_value`, `<t>_vs_check`, `<t>_check_ok`, `<t>_p_beat_check`, plus `checks_all_ok`,
  and `attr(, "check_reference_diagnostics")`.
  `spec` is the data frame from `ng_trait_check_spec()` (columns `trait`, `check`, `reject_if`),
  optionally with a `column_key` column. **`column_key` is what column names are built from;
  `trait` is what gets reported.** The runner sets the key to the sanitised trait name because
  the cross table names its columns that way (`ng_run_cp_clean_trait_name()`,
  `R/39_cross_prediction_runner.R:545`: `make.names()` then dots to underscores). When the
  column is absent the fallback is `trait`.
  `trait_values` is unused here — the cross's own `<t>_mean` / `<t>_pmv_used` columns are read
  from `scores`. `check_values` is a named list keyed by trait, each a named numeric over check
  ids.

- [ ] **Step 1: Write the failing test**

Append to `tests/check_reference.R`:

```r
# --- Task 3: reference columns, no row loss ---------------------------------
scores <- data.frame(
  parent1     = c("P1", "P1", "P2"),
  parent2     = c("P2", "P3", "P3"),
  yield_mean  = c(10, 4, 7),      # check at 6 -> above, below, above
  yield_pmv_used = c(4, 4, 0),    # sd = 2, 2, 0
  matur_mean  = c(70, 80, 75),    # check at 75, DECREASE -> below is good
  matur_pmv_used = c(1, 1, 1),
  stringsAsFactors = FALSE)

spec <- ng_trait_check_spec(trait = c("yield", "matur"), check = c("CHK_A", "CHK_B"),
                            trait_direction = c(yield = "increase", matur = "decrease"))
cv <- list(yield = c(CHK_A = 6), matur = c(CHK_B = 75))

out <- ng_attach_check_reference(scores, spec, trait_values = NULL, check_values = cv,
                                 k_progeny = 50L)

# ROW COUNT IS UNCHANGED -- this is the whole point
stopifnot(nrow(out) == 3L)
stopifnot(identical(out$parent1, scores$parent1))

# check id and value carried through
stopifnot(all(out$yield_check_id == "CHK_A"), all(out$yield_check_value == 6))

# direction-aware margin: POSITIVE ALWAYS MEANS BETTER, both directions
stopifnot(abs(out$yield_vs_check - c(4, -2, 1)) < 1e-8)      # increase: mean - check
stopifnot(abs(out$matur_vs_check - c(5, -5, 0)) < 1e-8)      # decrease: check - mean

# ok flags; a tie is NOT a violation
stopifnot(identical(out$yield_check_ok, c(TRUE, FALSE, TRUE)))
stopifnot(identical(out$matur_check_ok, c(TRUE, FALSE, TRUE)))   # row 3 is an exact tie

# combined flag
stopifnot(identical(out$checks_all_ok, c(TRUE, FALSE, TRUE)))

# P(beat check), increase trait: 1 - pnorm((tau-mu)/sd)^k
expect1 <- 1 - stats::pnorm((6 - 10) / 2)^50
stopifnot(abs(out$yield_p_beat_check[[1L]] - expect1) < 1e-10)
# a BELOW-check cross can still have a high tail probability -- the reason this column exists
expect2 <- 1 - stats::pnorm((6 - 4) / 2)^50
stopifnot(abs(out$yield_p_beat_check[[2L]] - expect2) < 1e-10, out$yield_p_beat_check[[2L]] > 0.5)
# sd == 0 degenerates to the indicator
stopifnot(out$yield_p_beat_check[[3L]] == 1)

# P(beat check), DECREASE trait: mirrored, 1 - pnorm((mu-tau)/sd)^k
expect_d <- 1 - stats::pnorm((70 - 75) / 1)^50
stopifnot(abs(out$matur_p_beat_check[[1L]] - expect_d) < 1e-10)

# a NOT-EVALUABLE check value yields NA columns and never a FALSE flag
out_na <- ng_attach_check_reference(scores, spec, NULL,
                                    list(yield = c(CHK_A = NA_real_), matur = c(CHK_B = 75)),
                                    k_progeny = 50L)
stopifnot(all(is.na(out_na$yield_check_value)), all(is.na(out_na$yield_check_ok)))
stopifnot(identical(out_na$checks_all_ok, c(TRUE, FALSE, TRUE)))   # driven by matur alone

d <- attr(out, "check_reference_diagnostics")
stopifnot(is.list(d), d$n_wrong_side$yield == 1L, d$n_wrong_side$matur == 1L)
stopifnot(d$n_not_evaluable == 0L)

# A trait whose name is not a valid R name: the cross table sanitises it for column names
# (make.names + dots to underscores) while the spec keeps what the breeder typed. Columns are
# looked up and written by the KEY; diagnostics and reporting use the RAW name.
scores_sp <- data.frame(parent1 = "P1", parent2 = "P2",
                        `Days_to_flower_mean` = 70, `Days_to_flower_pmv_used` = 1,
                        check.names = FALSE, stringsAsFactors = FALSE)
spec_sp <- ng_trait_check_spec("Days to flower", "CHK_B",
                               trait_direction = c(`Days to flower` = "decrease"))
spec_sp$column_key <- "Days_to_flower"
out_sp <- ng_attach_check_reference(scores_sp, spec_sp, NULL,
                                    list(`Days to flower` = c(CHK_B = 75)), k_progeny = 50L)
stopifnot("Days_to_flower_check_value" %in% names(out_sp))   # written by key
stopifnot(out_sp$Days_to_flower_check_ok)                     # 70 < 75, decrease -> good
stopifnot(names(attr(out_sp, "check_reference_diagnostics")$n_wrong_side) == "Days to flower")

cat("task 3 ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference.R`
Expected: FAIL with `could not find function "ng_attach_check_reference"`.

- [ ] **Step 3: Write minimal implementation**

Append to `R/51_check_reference.R`:

```r
# Attach per-trait check reference columns to a scored cross table. This function NEVER changes
# nrow(scores) and NEVER reorders it: a check informs the breeder, it does not decide for them.
# (The 0.14.0 predecessor, ng_apply_trait_checks(), dropped rows -- that is the behaviour this
# replaces.)
ng_attach_check_reference <- function(scores, spec, trait_values = NULL, check_values,
                                      k_progeny,
                                      mean_suffix = "_mean", sd_suffix = "_pmv_used") {
  scores <- as.data.frame(scores, stringsAsFactors = FALSE, check.names = FALSE)
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  n <- nrow(scores)
  k_progeny <- as.integer(k_progeny[[1L]])
  ok_mat <- matrix(NA, nrow = n, ncol = nrow(spec), dimnames = list(NULL, spec$trait))
  n_not_evaluable <- 0L
  n_wrong <- list()
  # The cross table names its per-trait columns with the SANITISED trait name
  # (ng_run_cp_clean_trait_name: make.names + dots to underscores), while the spec carries the
  # raw name the breeder typed. Look columns up by the key, report by the raw name.
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  for (k in seq_len(nrow(spec))) {
    tr <- spec$trait[[k]]; ck <- spec$check[[k]]; kk <- key[[k]]
    # reject_if == "below" means the breeder wants the mid-parent ABOVE the check (increase
    # trait); sgn flips the comparison so that a positive margin always means "better".
    sgn <- if (identical(spec$reject_if[[k]], "below")) 1 else -1
    mcol <- paste0(kk, mean_suffix); scol <- paste0(kk, sd_suffix)
    if (!(mcol %in% names(scores))) ng_stop("scores missing mean column for check trait: ", mcol)
    mu <- suppressWarnings(as.numeric(scores[[mcol]]))
    v <- if (scol %in% names(scores)) suppressWarnings(as.numeric(scores[[scol]])) else rep(NA_real_, n)
    tau <- suppressWarnings(as.numeric(check_values[[tr]][[ck]]))
    if (!length(tau)) tau <- NA_real_
    margin <- sgn * (mu - tau)
    ok <- margin >= 0                       # ties are not violations
    scores[[paste0(kk, "_check_id")]] <- ck
    scores[[paste0(kk, "_check_value")]] <- tau
    scores[[paste0(kk, "_vs_check")]] <- margin
    scores[[paste0(kk, "_check_ok")]] <- ok
    p <- rep(NA_real_, n)
    usable <- is.finite(mu) & is.finite(v) & v >= 0 & is.finite(tau)
    if (any(usable)) {
      # sgn folds the decrease case into the same closed form: negating both mu and tau turns
      # P(at least one of k progeny >= tau) into P(at least one <= tau).
      p[usable] <- ng_p_superior_progeny(sgn * mu[usable], sqrt(v[usable]),
                                         sgn * tau, k_progeny)
    }
    scores[[paste0(kk, "_p_beat_check")]] <- p
    ok_mat[, k] <- ok
    n_not_evaluable <- n_not_evaluable + sum(is.na(ok))
    n_wrong[[tr]] <- sum(ok %in% FALSE)
  }
  # NA never counts as a failure: an unevaluable check is reported, not held against a cross.
  scores$checks_all_ok <- !apply(ok_mat, 1L, function(r) any(r %in% FALSE))
  attr(scores, "check_reference_diagnostics") <- list(
    active = spec, n_wrong_side = n_wrong,
    n_not_evaluable = as.integer(n_not_evaluable), n_candidates = n)
  rownames(scores) <- NULL
  scores
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript tests/check_reference.R`
Expected: `task 3 ok`

- [ ] **Step 5: Commit**

```bash
git add R/51_check_reference.R tests/check_reference.R
git commit -m "feat(checks): attach reference columns and P(beat check) without dropping rows"
```

---

### Task 4: Multi-trait P(beat every check at once)

`ng_p_superior_progeny_multitrait()` already takes `tau_lower` and `tau_upper`. Increase-traits
fill `tau_lower`, decrease-traits fill `tau_upper`.

**Files:**
- Modify: `R/51_check_reference.R`
- Test: `tests/check_reference.R`

**Interfaces:**
- Consumes: `ng_attach_check_reference()` output.
- Produces: `ng_check_tau_bounds(spec, check_values)` → `list(tau_lower = <named numeric>, tau_upper = <named numeric>)`, both over `spec$trait`, `-Inf`/`Inf` where that direction does not apply.
  `ng_attach_joint_check_probability(scores, spec, check_values, k_progeny, mean_suffix = "_mean", sd_suffix = "_pmv_used", cross_trait_cov = NULL)`
  → `scores` with one added column `p_beat_all_checks`.

- [ ] **Step 1: Write the failing test**

Append to `tests/check_reference.R`:

```r
# --- Task 4: tau bounds for the multi-trait joint probability ---------------
b <- ng_check_tau_bounds(spec, cv)
stopifnot(identical(names(b$tau_lower), c("yield", "matur")))
# increase trait -> a LOWER bound at the check; decrease trait -> an UPPER bound
stopifnot(b$tau_lower[["yield"]] == 6, is.infinite(b$tau_upper[["yield"]]), b$tau_upper[["yield"]] > 0)
stopifnot(b$tau_upper[["matur"]] == 75, is.infinite(b$tau_lower[["matur"]]), b$tau_lower[["matur"]] < 0)

# an unevaluable check widens to an unbounded side rather than excluding the trait
b_na <- ng_check_tau_bounds(spec, list(yield = c(CHK_A = NA_real_), matur = c(CHK_B = 75)))
stopifnot(is.infinite(b_na$tau_lower[["yield"]]), b_na$tau_lower[["yield"]] < 0)

# the joint column: P(a progeny beats EVERY check at once)
j <- ng_attach_joint_check_probability(scores, spec, cv, k_progeny = 50L)
stopifnot("p_beat_all_checks" %in% names(j), nrow(j) == 3L)
stopifnot(all(j$p_beat_all_checks >= 0 & j$p_beat_all_checks <= 1, na.rm = TRUE))
# beating both checks can never be more likely than beating either one alone
single <- ng_attach_check_reference(scores, spec, NULL, cv, k_progeny = 50L)
stopifnot(all(j$p_beat_all_checks <= single$yield_p_beat_check + 1e-8))
stopifnot(all(j$p_beat_all_checks <= single$matur_p_beat_check + 1e-8))

cat("task 4 ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference.R`
Expected: FAIL with `could not find function "ng_check_tau_bounds"`.

- [ ] **Step 3: Write minimal implementation**

Append to `R/51_check_reference.R`:

```r
# Translate a check spec into the (tau_lower, tau_upper) pair that
# ng_p_superior_progeny_multitrait() consumes: an increase trait bounds the progeny from below
# at its check, a decrease trait bounds it from above. A check with no evaluable value leaves
# that trait unbounded rather than dropping it, so the joint probability stays defined.
ng_check_tau_bounds <- function(spec, check_values) {
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  tr <- spec$trait
  lo <- stats::setNames(rep(-Inf, length(tr)), tr)
  hi <- stats::setNames(rep(Inf, length(tr)), tr)
  for (k in seq_along(tr)) {
    tau <- suppressWarnings(as.numeric(check_values[[tr[[k]]]][[spec$check[[k]]]]))
    if (!length(tau) || !is.finite(tau)) next
    if (identical(spec$reject_if[[k]], "below")) lo[[tr[[k]]]] <- tau else hi[[tr[[k]]]] <- tau
  }
  list(tau_lower = lo, tau_upper = hi)
}

# P(a progeny beats EVERY check at once). Thin wrapper over the existing multi-trait threshold
# machinery: the check values ARE the tau bounds, so no new probability model is introduced.
ng_attach_joint_check_probability <- function(scores, spec, check_values, k_progeny,
                                              mean_suffix = "_mean", sd_suffix = "_pmv_used",
                                              cross_trait_cov = NULL) {
  spec <- as.data.frame(spec, stringsAsFactors = FALSE)
  b <- ng_check_tau_bounds(spec, check_values)
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  trait_specs <- data.frame(
    trait = spec$trait,
    mean_col = paste0(key, mean_suffix),
    var_col = paste0(key, sd_suffix),
    stringsAsFactors = FALSE)
  ng_add_p_superior_progeny_multitrait(
    scores, trait_specs,
    tau_lower = as.numeric(b$tau_lower[spec$trait]),
    tau_upper = as.numeric(b$tau_upper[spec$trait]),
    k_progeny = k_progeny, cross_trait_cov = cross_trait_cov,
    out_col = "p_beat_all_checks")
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript tests/check_reference.R`
Expected: `task 4 ok`

- [ ] **Step 5: Commit**

```bash
git add R/51_check_reference.R tests/check_reference.R
git commit -m "feat(checks): tau bounds for the multi-trait joint beat-the-check probability"
```

---

### Task 5: Delete the veto

Remove the old behaviour before wiring the new one, so no code path can reach both.

**Files:**
- Delete: `R/44_trait_checks.R` — except `ng_trait_check_spec()`, which moves to
  `R/51_check_reference.R` unchanged.
- Delete: `tests/trait_checks.R`
- Modify: `R/39_cross_prediction_runner.R:1190-1250` (the `trait_checks` block)
- Modify: `NAMESPACE`, `man/nextgenCrossDesign-api.Rd`
- Modify: `vignettes/nextgenCrossDesign.Rmd:1746` — the API table lists
  `ng_apply_trait_checks()`, which this task deletes. This repo executes its vignette, so a
  stale reference fails `R CMD check` at Task 10. Replace it with
  `ng_attach_check_reference()` in that row.

**Interfaces:**
- Consumes: nothing.
- Produces: `ng_trait_check_spec(trait, check, direction = NA, basis = "gebv", trait_direction = NULL)`
  now lives in `R/51_check_reference.R`. **The `basis` argument is removed** — its signature
  becomes `ng_trait_check_spec(trait, check, direction = NA, trait_direction = NULL)` returning
  columns `trait`, `check`, `reject_if` (no `basis` column).

- [ ] **Step 1: Move the spec builder and drop `basis`**

Copy `ng_trait_check_spec()` from `R/44_trait_checks.R:8-38` into `R/51_check_reference.R`,
deleting the `basis` argument, the two `basis` lines, and the `basis` column:

```r
# Build a validated per-trait check spec (one check per trait). `direction` NA is filled from
# `trait_direction`: increase -> the breeder wants the mid-parent ABOVE the check, so a cross is
# on the wrong side when it is below. There is no `basis` argument: the check always follows the
# run's per-trait mean_source, which is what keeps the reference on the plotted scale.
ng_trait_check_spec <- function(trait, check, direction = NA, trait_direction = NULL) {
  trait <- as.character(trait); check <- as.character(check)
  n <- length(trait)
  if (!n) ng_stop("ng_trait_check_spec needs at least one trait")
  if (length(check) != n) ng_stop("check must be one per trait")
  direction <- tolower(as.character(rep_len(direction, n)))
  need <- is.na(direction) | !nzchar(direction) | direction == "na"
  if (any(need)) {
    if (is.null(trait_direction))
      ng_stop("direction is NA and no trait_direction supplied to resolve it for: ",
              paste(trait[need], collapse = ", "))
    td <- tolower(as.character(trait_direction[trait[need]]))
    if (anyNA(td)) ng_stop("trait_direction has no entry for: ",
                           paste(trait[need][is.na(td)], collapse = ", "))
    direction[need] <- ifelse(td == "increase", "below",
                       ifelse(td == "decrease", "above", NA_character_))
  }
  if (any(!direction %in% c("above", "below")))
    ng_stop("direction must resolve to 'above' or 'below'")
  dup <- unique(trait[duplicated(trait)])
  if (length(dup))
    ng_stop("ng_trait_check_spec: one check per trait; duplicate trait(s): ",
            paste(dup, collapse = ", "))
  data.frame(trait = trait, check = check, reject_if = direction, stringsAsFactors = FALSE)
}
```

- [ ] **Step 2: Delete the old files and the runner block**

```bash
git rm R/44_trait_checks.R tests/trait_checks.R
```

In `R/39_cross_prediction_runner.R`, delete the whole `if (!is.null(trait_checks)) { ... }`
block at lines 1190-1250 (from the `trait_check_diagnostics <- NULL` line through the closing
brace of the zero-survivors `ng_stop`), and remove the formals `check_basis` and
`exclude_threshold_violators` from the function signature. Leave `trait_checks` in the
signature — Task 6 rewires it.

- [ ] **Step 3: Drop the retired exports**

In `NAMESPACE`, remove `export(ng_apply_trait_checks)`. In
`man/nextgenCrossDesign-api.Rd`, remove the `ng_apply_trait_checks` alias.

- [ ] **Step 4: Keep the moved spec builder under test**

`tests/trait_checks.R` was the only coverage of `ng_trait_check_spec()` and it is gone. Append
its surviving assertions to `tests/check_reference.R`, minus everything about `basis`:

```r
# --- Task 5: the spec builder, moved and basis-free ------------------------
td <- c(yield = "increase", maturity = "decrease")
s <- ng_trait_check_spec(trait = c("yield", "maturity"), check = c("CkY", "CkM"),
                         trait_direction = td)
stopifnot(nrow(s) == 2L, !("basis" %in% names(s)))
stopifnot(s$reject_if[s$trait == "yield"] == "below")     # increase -> wrong side is below
stopifnot(s$reject_if[s$trait == "maturity"] == "above")  # decrease -> wrong side is above
s2 <- ng_trait_check_spec("protein", "CkP", direction = "above",
                          trait_direction = c(protein = "increase"))
stopifnot(s2$reject_if == "above")                        # explicit override wins
err_d <- tryCatch(ng_trait_check_spec("x", "C", direction = "sideways"),
                  error = function(e) conditionMessage(e))
stopifnot(is.character(err_d), grepl("direction", err_d))
err_dup <- tryCatch(ng_trait_check_spec(c("yield", "yield"), c("CkA", "CkB"),
                                        trait_direction = c(yield = "increase")),
                    error = function(e) conditionMessage(e))
stopifnot(is.character(err_dup), grepl("duplicate", err_dup))
# basis is GONE: passing it must be an unused-argument error, not silently ignored
err_b <- tryCatch(ng_trait_check_spec("yield", "CkY", basis = "phenotype",
                                      trait_direction = c(yield = "increase")),
                  error = function(e) conditionMessage(e))
stopifnot(is.character(err_b), grepl("unused argument", err_b))

cat("task 5 ok\n")
```

Run: `Rscript -e 'devtools::load_all("."); cat("loaded\n")'`
Expected: `loaded` with no error about a missing function.

Run: `Rscript tests/check_reference.R`
Expected: tasks 1-5 all pass.

- [ ] **Step 5: Commit**

```bash
git add -A R/ tests/ NAMESPACE man/
git commit -m "refactor(checks)!: remove the trait-check veto, keep the spec builder

BREAKING: exclude_threshold_violators and check_basis are gone. Checks no
longer filter candidate crosses; the reference path replaces them."
```

---

### Task 6: Wire the runner

**Files:**
- Modify: `R/39_cross_prediction_runner.R` (signature; new block where the veto was; the
  `utils::globalVariables()` list at :1631; the `ctx$trait_check_diagnostics` plumbing at
  :1279 / :1678 / :1719; JSON envelope)
- Modify: `docs/frontend/contracts/config_schema.json` (document `check_progeny_size`)
- Modify: `tests/contract_schema_drift.R` (`undocumented_ok`)
- Modify: `README.md` (lines ~170, ~181, ~335-336)
- Modify: `vignettes/nextgenCrossDesign.Rmd` (the `{r trait_checks}` section, ~line 793-820)
- Test: `tests/check_reference_runner.R`

**Interfaces:**
- Consumes: `ng_align_check_geno()`, `ng_check_reference_value()`, `ng_attach_check_reference()`,
  `ng_check_tau_bounds()`, `ng_trait_check_spec()`.
- Produces: new runner formals `check_geno = NULL`, `check_records = NULL`,
  `check_progeny_size = NULL` (**required whenever `trait_checks` is supplied** — hard error if
  missing or < 1; progeny per family is the breeder's figure, never a package default).
  Result gains `$trait_check_reference` =
  `list(active = <spec df>, values = <named list by trait>, source = <named chr by trait>, progeny_size = <int>, diagnostics = <list>)`.
- **`check_records` shape:** a list keyed by **trait**, each element a list keyed by **source**
  (`BLUP`, `BLUE`, `adjusted_pheno`), each of those a named numeric over check ids — e.g.
  `list(yield = list(BLUE = c(CHK_A = 11.5)))`. It is only consulted when the run's source for
  that trait is not a GEBV variant. `NULL` is normal: a GEBV run needs nothing here.

- [ ] **Step 1: Write the failing test**

Create `tests/check_reference_runner.R`:

```r
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(11)
n_p <- 12L; n_m <- 40L
geno <- matrix(rbinom(n_p * n_m, 2, 0.4), nrow = n_p,
               dimnames = list(paste0("P", seq_len(n_p)), paste0("m", seq_len(n_m))))
# the check is NOT among the parents -- the case the old veto rejected outright
chk <- matrix(rbinom(2 * n_m, 2, 0.4), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), colnames(geno)))
pheno <- data.frame(id = rownames(geno),
                    yield = rnorm(n_p, 10, 2),
                    stringsAsFactors = FALSE)

args <- list(geno = geno, pheno = pheno, id_col = "id",
             trait_direction = c(yield = "increase"), n_crosses = 5L)

res <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE))))

ref <- res$trait_check_reference
stopifnot(is.list(ref), nrow(ref$active) == 1L, ref$active$check == "CHK_A")
stopifnot(is.finite(ref$values$yield[["CHK_A"]]))
stopifnot(nzchar(ref$source[["yield"]]))

ct <- res$candidate_crosses
stopifnot(all(c("yield_check_value", "yield_vs_check", "yield_check_ok",
                "yield_p_beat_check", "checks_all_ok") %in% names(ct)))
# the check is never a parent of any candidate cross
stopifnot(!any(c(ct$parent1, ct$parent2) %in% c("CHK_A", "CHK_B")))

# a check id colliding with a parent id is rejected
err <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = matrix(0, nrow = 1, ncol = n_m,
                      dimnames = list("P1", colnames(geno))),
  check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "P1", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("also a candidate parent", err))

# a check named in trait_checks but absent from check_geno is rejected
err2 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L,
  trait_checks = data.frame(trait = "yield", check = "NOPE", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err2), grepl("NOPE", err2))

# check_progeny_size is REQUIRED with trait_checks -- never defaulted
err3 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err3), grepl("check_progeny_size", err3))
# and it must be a sensible count
err4 <- tryCatch(do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 0L,
  trait_checks = data.frame(trait = "yield", check = "CHK_A", stringsAsFactors = FALSE)))),
  error = function(e) conditionMessage(e))
stopifnot(is.character(err4), grepl("check_progeny_size", err4))

# the size the user gave is what drives the probability, and is stamped on the result
stopifnot(identical(ref$progeny_size, 200L))

cat("runner wiring ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference_runner.R`
Expected: FAIL with an unused-argument error for `check_geno`.

- [ ] **Step 3: Write the implementation**

Add to the `ng_run_cross_prediction()` signature, next to the existing `trait_checks`:

```r
                                   check_geno = NULL,
                                   check_records = NULL,
                                   check_progeny_size = NULL,
```

Replace the deleted veto block (same position — after the lethal guard, before
`ng_score_breeder_objective()`) with:

```r
  # Check lines as REFERENCES (Module: check reference). Attaches per-trait reference columns
  # and never touches the candidate set: no filtering, no penalty, no reordering. Checks are
  # kept out of geno entirely, so QC, LD, the GRM, and ng_make_pairs() are all untouched and a
  # run with checks is numerically identical to the same run without them.
  trait_check_reference <- NULL
  if (!is.null(trait_checks)) {
    if (is.null(check_geno)) {
      ng_stop("trait_checks needs check_geno: check lines are supplied in their own genotype ",
              "matrix, separate from the candidate parents, and are never crossed")
    }
    if (!identical(prediction_mode, "trait_by_trait")) {
      ng_stop("trait_checks requires prediction_mode = 'trait_by_trait' (checks are keyed by trait)")
    }
    if (!(inherits(trait_checks, "data.frame") && all(c("trait", "check") %in% names(trait_checks)))) {
      ng_stop("trait_checks must be a data.frame with trait + check columns")
    }
    # Progeny per family is the breeder's number, not ours. It scales P(beat check) directly,
    # so there is no default: a made-up family size would silently drive a reported probability.
    kp <- if (is.null(check_progeny_size)) NA_integer_ else
      suppressWarnings(as.integer(check_progeny_size[[1L]]))
    if (!length(kp) || is.na(kp) || kp < 1L) {
      ng_stop("check_progeny_size is required with trait_checks: give the number of progeny ",
              "you will raise per family. It sets P(beat check) -- the chance a cross throws a ",
              "line beating the check -- so it must be your program's figure, not a default.")
    }
    check_geno <- ng_align_check_geno(check_geno, colnames(geno), ploidy = marker_ploidy)
    clash <- intersect(rownames(check_geno), rownames(geno))
    if (length(clash)) {
      ng_stop("check id(s) also a candidate parent: ", paste(clash, collapse = ", "),
              ". A check is a benchmark, not breeding material; give it a distinct id or ",
              "remove it from the parent genotypes.")
    }
    tdir <- stats::setNames(direction_canonical$direction, direction_canonical$trait)
    tc_direction <- if (is.null(trait_checks$direction)) NA else trait_checks$direction
    tc_spec <- ng_trait_check_spec(trait_checks$trait, trait_checks$check,
                                   direction = tc_direction, trait_direction = tdir)
    # cross_table columns are named with the sanitised trait name; carry it as the lookup key
    # so a trait like "Days to flower" resolves to Days_to_flower_mean rather than erroring.
    tc_spec$column_key <- ng_run_cp_clean_trait_name(tc_spec$trait)
    missing_chk <- setdiff(tc_spec$check, rownames(check_geno))
    if (length(missing_chk)) {
      ng_stop("trait_checks names check line(s) absent from check_geno: ",
              paste(missing_chk, collapse = ", "))
    }
    check_values <- list(); check_source <- character(0)
    for (tr in unique(tc_spec$trait)) {
      if (!(tr %in% names(effects_list))) {
        ng_stop("trait_checks references a trait not present in trait_direction/effects: ", tr)
      }
      # The per-trait mean_source stamped on the scored table is the authority: the check must
      # land on whatever source produced the cross means for this trait, or the reference line
      # would sit on a different scale from the axis it is drawn on.
      src <- as.character(trait_mean_source[[tr]])
      check_source[[tr]] <- src
      check_values[[tr]] <- ng_check_reference_value(src, check_geno, effects_list[[tr]],
                                                     check_records = check_records[[tr]])
    }
    cross_table <- ng_attach_check_reference(cross_table, tc_spec, trait_values = NULL,
                                             check_values = check_values,
                                             k_progeny = kp)
    # Multi-trait only: P(a progeny beats every check at once). Skipped for a single check,
    # where p_beat_all_checks would just duplicate the per-trait column.
    if (nrow(tc_spec) > 1L) {
      cross_table <- ng_attach_joint_check_probability(
        cross_table, tc_spec, check_values, k_progeny = kp,
        cross_trait_cov = ctc)
    }
    trait_check_reference <- list(
      active = tc_spec, values = check_values, source = check_source,
      progeny_size = kp,
      diagnostics = attr(cross_table, "check_reference_diagnostics"))
    attr(cross_table, "check_reference_diagnostics") <- NULL
  }
```

In the per-trait loop at `R/39_cross_prediction_runner.R:1048-1065`, capture the source so the
block above can read it. Immediately after
`cross_table[[paste0(clean_trait, "_mean")]] <- scored_trait$cross_mean_blend`, add:

```r
    if (!exists("trait_mean_source", inherits = FALSE)) trait_mean_source <- list()
    trait_mean_source[[trait]] <- scored_trait$mean_source[[1L]]
```

and initialise `trait_mean_source <- list()` just before the `for (j in seq_along(trait_results))`
loop.

The joint-probability call above reads `ctc`, the exact within-family cross-trait covariance
table built at `R/39_cross_prediction_runner.R:1106`. It is currently assigned **inside** the
`if` block that computes it, so it does not exist when that branch is skipped. Hoist
`ctc <- NULL` to just before that `if` so the reference is always defined; when it is `NULL`
the multi-trait helper falls back to its `G_hat` proxy, which is the existing behaviour
elsewhere.

Finally, replace the old result field with the new one. `trait_check_diagnostics` is currently
plumbed through the ctx at three places — `ctx$trait_check_diagnostics <- ...` (:1279),
`trait_check_diagnostics <- ctx$trait_check_diagnostics` (:1678), and
`trait_check_diagnostics = trait_check_diagnostics` in the returned list (:1719). Rename all
three to `trait_check_reference`. Nothing should still reference the old name:

```bash
git grep -n "trait_check_diagnostics" -- R/   # must return nothing when you are done
```

- [ ] **Step 4: Keep R CMD check and the contract test clean**

The ctx fields are unpacked via `list2env()`, so `utils::globalVariables()` at
`R/39_cross_prediction_runner.R:1631` declares them for `R CMD check`. Remove the now-deleted
`"check_basis"` and `"exclude_threshold_violators"` entries and add the new ctx fields in
alphabetical position: `"check_geno"`, `"check_progeny_size"`, `"check_records"`,
`"trait_check_reference"`.

`tests/contract_schema_drift.R` compares `formals(ng_run_cross_prediction)` against
`docs/frontend/contracts/config_schema.json`. Only the *documented-but-not-a-formal* direction
is fatal, so removing `check_basis` is safe (it was never in the schema). Two additions keep
the NOTE honest:

- `check_progeny_size` is a plain number the frontend must be able to discover — add it to
  `config_schema.json` in the group that holds the other cross-filter parameters:

```json
{"name": "check_progeny_size", "type": "integer", "required": false,
 "desc": "Progeny per family you will raise. Required with trait_checks; sets P(beat check). No default -- it must be your program's figure."}
```

- `check_geno` and `check_records` are in-memory matrices/lists that a JSON config cannot
  carry, which is exactly what `undocumented_ok` in `tests/contract_schema_drift.R:36-40`
  exists for. Add them there beside `genotype` / `marker_map`.

- [ ] **Step 5: Re-document the user-facing API**

The veto is gone from the code but still described in two places, and one of them executes.

`vignettes/nextgenCrossDesign.Rmd` (~793-820) has a **live** `{r trait_checks}` chunk that
calls the runner with `check_basis = "gebv"` and a check of `"P001"` — a candidate parent,
which the new code rejects by design — and then reads `yield_check_violation`, `threshold_ok`,
and `r$trait_check_diagnostics`. Every one of those is removed or renamed. Rewrite the prose to
describe a reference rather than a filter, and replace the chunk with:

```r
# Check lines are benchmarks, never crossed: they come in their own genotype matrix.
chk <- geno[1:2, , drop = FALSE]
rownames(chk) <- c("CHECK_A", "CHECK_B")
tc <- data.frame(trait = "yield", check = "CHECK_A", stringsAsFactors = FALSE)
r <- ng_run_cross_prediction(phenotype = pheno, genotype = geno, marker_map = map,
  trait_direction = direction, id_col = "NAME", use_ocs = TRUE, n_crosses = 100,
  trait_checks = tc, check_geno = chk, check_progeny_size = 200,
  include_trait_gebv = TRUE)
# candidate_crosses gains per-trait reference columns -- nothing is dropped:
head(r$candidate_crosses[, c("parent1", "parent2", "yield_check_value",
                             "yield_vs_check", "yield_check_ok", "yield_p_beat_check")])
r$trait_check_reference$active          # which check is active for which trait
r$trait_check_reference$progeny_size    # the family size that produced P(beat check)
```

Note the chunk builds `chk` from rows of `geno` but renames them, so the ids do not collide
with candidate parents — a check id that IS a parent is now a hard error.

`README.md` describes `check_basis` and `exclude_threshold_violators` at ~170, ~181 and
~335-336. Rewrite those to the reference vocabulary: checks are supplied via `check_geno`,
require `check_progeny_size`, flag nothing out of the candidate pool, and add
`<trait>_vs_check` / `<trait>_p_beat_check` columns.

- [ ] **Step 6: Run test to verify it passes**

Run: `Rscript tests/check_reference_runner.R`
Expected: `runner wiring ok`

Then confirm the rewritten vignette chunk actually executes:

Run: `Rscript -e 'devtools::load_all("."); knitr::knit("vignettes/nextgenCrossDesign.Rmd", output = tempfile())'`
Expected: completes without error. This is the step that proves the vignette is not carrying a
removed argument.

- [ ] **Step 7: Commit**

```bash
git add R/39_cross_prediction_runner.R docs/frontend/contracts/config_schema.json \
        tests/contract_schema_drift.R tests/check_reference_runner.R \
        README.md vignettes/nextgenCrossDesign.Rmd
git commit -m "feat(checks): wire check_geno through the runner as a reference layer"
```

---

### Task 7: The invariant — checks change nothing

This is the test that makes "reference only" a fact rather than an intention.

**Files:**
- Test: `tests/check_reference_invariant.R`

**Interfaces:**
- Consumes: the runner from Task 6.
- Produces: nothing.

- [ ] **Step 1: Write the failing test**

Create `tests/check_reference_invariant.R`:

```r
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(202)
n_p <- 14L; n_m <- 50L
geno <- matrix(rbinom(n_p * n_m, 2, 0.35), nrow = n_p,
               dimnames = list(paste0("P", seq_len(n_p)), paste0("m", seq_len(n_m))))
chk <- matrix(rbinom(2 * n_m, 2, 0.35), nrow = 2,
              dimnames = list(c("CHK_A", "CHK_B"), colnames(geno)))
pheno <- data.frame(id = rownames(geno), yield = rnorm(n_p, 10, 2),
                    protein = rnorm(n_p, 12, 1), stringsAsFactors = FALSE)
args <- list(geno = geno, pheno = pheno, id_col = "id",
             trait_direction = c(yield = "increase", protein = "increase"),
             n_crosses = 6L)

without <- do.call(ng_run_cross_prediction, args)
with_ck <- do.call(ng_run_cross_prediction, c(args, list(
  check_geno = chk, check_progeny_size = 200L,
  trait_checks = data.frame(trait = c("yield", "protein"),
                            check = c("CHK_A", "CHK_B"), stringsAsFactors = FALSE))))

# every column the run produced WITHOUT checks must be untouched BY checks
a <- without$candidate_crosses
b <- with_ck$candidate_crosses
stopifnot(nrow(a) == nrow(b))
for (nm in names(a)) {
  stopifnot(nm %in% names(b))
  stopifnot(isTRUE(all.equal(a[[nm]], b[[nm]], tolerance = 0)))
}
# the plan is identical too -- allocation must not shift
stopifnot(isTRUE(all.equal(without$crossing_plan, with_ck$crossing_plan, tolerance = 0)))

# and the only difference is ADDED columns
added <- setdiff(names(b), names(a))
stopifnot(length(added) > 0L)
stopifnot(all(grepl("_check_|^checks_all_ok$|_vs_check$", added)))

cat("invariant ok: checks add columns and change nothing else\n")
```

- [ ] **Step 2: Run the test**

Run: `Rscript tests/check_reference_invariant.R`
Expected: PASS. If it fails, a check has leaked into QC, LD, the GRM, or pair enumeration —
find the leak; do not relax the tolerance.

- [ ] **Step 3: Commit**

```bash
git add tests/check_reference_invariant.R
git commit -m "test(checks): a run with checks must be numerically identical to one without"
```

---

### Task 8: Workbook — `Checks` sheet and per-trait columns

**Files:**
- Modify: `R/37_cross_priority_workbook.R` (`ng_cpw_make_selected` :174, `ng_cpw_candidate_table`, the sheet list at :348-360)
- Test: `tests/check_reference_workbook.R`

**Interfaces:**
- Consumes: `trait_check_reference` from Task 6.
- Produces: `ng_cpw_checks_sheet(trait_check_reference, crosses)` → data frame with columns
  `trait`, `check_id`, `direction`, `value`, `source`, `n_crosses_on_wrong_side`.
  `ng_build_cross_priority_workbook()` gains `trait_check_reference = NULL` and emits a
  `Checks` sheet when it is non-NULL.

- [ ] **Step 1: Write the failing test**

Create `tests/check_reference_workbook.R`:

```r
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

crosses <- data.frame(
  parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
  yield_check_id = "CHK_A", yield_check_value = 6,
  yield_vs_check = c(4, -2), yield_check_ok = c(TRUE, FALSE),
  yield_p_beat_check = c(0.98, 0.71), checks_all_ok = c(TRUE, FALSE),
  stringsAsFactors = FALSE)

ref <- list(
  active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                      stringsAsFactors = FALSE),
  values = list(yield = c(CHK_A = 6)),
  source = c(yield = "GEBV"),
  diagnostics = list(n_wrong_side = list(yield = 1L), n_not_evaluable = 0L, n_candidates = 2L))

sh <- ng_cpw_checks_sheet(ref, crosses)
stopifnot(nrow(sh) == 1L)
stopifnot(sh$trait == "yield", sh$check_id == "CHK_A")
stopifnot(sh$direction == "want above")          # reject_if 'below' reads as 'want above'
stopifnot(sh$value == 6, sh$source == "GEBV")
stopifnot(sh$n_crosses_on_wrong_side == "1 / 2")

# a decrease trait reads the other way round
ref2 <- ref; ref2$active$reject_if <- "above"
stopifnot(ng_cpw_checks_sheet(ref2, crosses)$direction == "want below")

cat("workbook checks sheet ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference_workbook.R`
Expected: FAIL with `could not find function "ng_cpw_checks_sheet"`.

- [ ] **Step 3: Write the implementation**

Add to `R/37_cross_priority_workbook.R`:

```r
# The `Checks` sheet: the reference itself, one row per trait. Breeders read the check first
# and the crosses second, so it gets its own sheet rather than being inferred from columns.
ng_cpw_checks_sheet <- function(trait_check_reference, crosses) {
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  if (!nrow(spec)) return(data.frame())
  n_total <- nrow(crosses)
  wrong <- trait_check_reference$diagnostics$n_wrong_side
  data.frame(
    trait = spec$trait,
    check_id = spec$check,
    direction = ifelse(spec$reject_if == "below", "want above", "want below"),
    value = vapply(seq_len(nrow(spec)), function(k) {
      suppressWarnings(as.numeric(trait_check_reference$values[[spec$trait[[k]]]][[spec$check[[k]]]]))
    }, numeric(1)),
    source = as.character(trait_check_reference$source[spec$trait]),
    n_crosses_on_wrong_side = vapply(spec$trait, function(tr) {
      w <- wrong[[tr]]
      sprintf("%d / %d", if (is.null(w)) NA_integer_ else as.integer(w), n_total)
    }, character(1)),
    stringsAsFactors = FALSE, row.names = NULL)
}
```

In `ng_build_cross_priority_workbook()`, add the formal `trait_check_reference = NULL` and,
immediately after the `Trait_Directions` entry in the `out <- list(...)` at :348, insert the
sheet conditionally:

```r
  if (!is.null(trait_check_reference)) {
    out$Checks <- ng_cpw_checks_sheet(trait_check_reference, selected)
  }
```

In `ng_cpw_make_selected()` and `ng_cpw_candidate_table()`, extend the retained-column vector
(`keep` at :259) so the check columns survive into the sheets:

```r
  keep <- c(keep, grep("_check_id$|_check_value$|_vs_check$|_check_ok$|_p_beat_check$",
                       names(crosses), value = TRUE),
            "checks_all_ok", "p_beat_all_checks")
  keep <- intersect(unique(keep), names(crosses))
```

- [ ] **Step 4: Wire the sheet into real runs**

`ng_run_cp_output_files()` (`R/39_cross_prediction_runner.R:609`) is the only call site of
`ng_write_cross_priority_workbook()` (:644). Without this step the sheet exists but no run
produces it. Add the formal:

```r
                                   include_trait_gebv = FALSE,
                                   trait_check_reference = NULL) {
```

pass it through to the workbook call:

```r
      include_trait_gebv = isTRUE(include_trait_gebv),
      trait_check_reference = trait_check_reference
```

thread it through `ng_write_cross_priority_workbook()` into
`ng_build_cross_priority_workbook()`, and at the `ng_run_cp_output_files(...)` call site in the
runner add `trait_check_reference = trait_check_reference`. The `NULL` default means a
no-checks run is byte-identical to today.

- [ ] **Step 5: Run test to verify it passes**

Run: `Rscript tests/check_reference_workbook.R`
Expected: `workbook checks sheet ok`

Run: `Rscript tests/cross_priority_workbook.R`
Expected: still passes (no-checks runs must produce the workbook they always did).

- [ ] **Step 6: Commit**

```bash
git add R/37_cross_priority_workbook.R R/39_cross_prediction_runner.R tests/check_reference_workbook.R
git commit -m "feat(checks): Checks sheet and per-trait reference columns in the workbook"
```

---

### Task 9: Plot — horizontal reference line, flag styling, opacity by P(beat check)

A check has no diversity coordinate, so it is a line spanning the full x-range, never a marker.

**Files:**
- Modify: `R/38_cross_priority_plots.R` (`ng_plot_priority_score_vs_kinship` :22)
- Test: `tests/check_reference_plot.R`

**Interfaces:**
- Consumes: `trait_check_reference` from Task 6.
- Produces: `ng_plot_priority_score_vs_kinship(..., check_line = NULL, check_label = NULL)`
  where `check_line` is a single numeric or `NULL`. Also
  `ng_check_line_value(trait_check_reference, multi_trait_meta, trait = NULL)` → single numeric
  or `NA_real_` when no valid line exists for the current index.

- [ ] **Step 1: Write the failing test**

Create `tests/check_reference_plot.R`:

```r
ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

ref <- list(active = data.frame(trait = "yield", check = "CHK_A", reject_if = "below",
                                stringsAsFactors = FALSE),
            values = list(yield = c(CHK_A = 6)), source = c(yield = "GEBV"))

# SINGLE TRAIT: the y axis IS the trait mean, so the line is exact
stopifnot(ng_check_line_value(ref, multi_trait_meta = NULL, trait = "yield") == 6)

# LINEAR INDEX: apply the same coefficients to the check's per-trait values
ref2 <- list(active = data.frame(trait = c("yield", "protein"), check = c("CHK_A", "CHK_A"),
                                 reject_if = "below", stringsAsFactors = FALSE),
             values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
             source = c(yield = "GEBV", protein = "GEBV"))
meta_lin <- list(method = "weighted", weights = c(yield = 0.5, protein = 0.5))
stopifnot(abs(ng_check_line_value(ref2, meta_lin) - 8) < 1e-8)

# RANK-BASED INDEX: a check has no rank, so there is NO valid line
meta_rank <- list(method = "rank_threshold")
stopifnot(is.na(ng_check_line_value(ref2, meta_rank)))

# an unevaluable check yields no line rather than a misplaced one
ref3 <- ref; ref3$values$yield <- c(CHK_A = NA_real_)
stopifnot(is.na(ng_check_line_value(ref3, NULL, trait = "yield")))

# the plot accepts a line and writes a file without error
scored <- data.frame(parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
                     multi_trait_score = c(9, 4), pair_kinship = c(0.1, 0.2),
                     priority_tier = c("priority", "low_priority"), stringsAsFactors = FALSE)
p <- file.path(tempdir(), "chk-plot.png")
ng_plot_priority_score_vs_kinship(scored, output_path = p, check_line = 6,
                                  check_label = "CHK_A")
stopifnot(file.exists(p), file.size(p) > 0)

cat("check plot ok\n")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript tests/check_reference_plot.R`
Expected: FAIL with `could not find function "ng_check_line_value"`.

- [ ] **Step 3: Write the implementation**

Add to `R/38_cross_priority_plots.R`:

```r
# The y value at which to draw the check reference line. Single-trait runs plot the trait mean
# itself, so the check's value is the line. A LINEAR index (weighted / economic) is a fixed
# combination of trait values, so the check's index value is exact. A RANK-based index is a
# function of the candidate distribution, and a check has no rank because it is not a cross --
# there is no honest line, so we return NA and the caller says so rather than drawing a number
# that looks authoritative and is not.
ng_check_line_value <- function(trait_check_reference, multi_trait_meta = NULL, trait = NULL) {
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  if (!nrow(spec)) return(NA_real_)
  val <- function(tr) {
    ck <- spec$check[[match(tr, spec$trait)]]
    suppressWarnings(as.numeric(trait_check_reference$values[[tr]][[ck]]))
  }
  if (!is.null(trait) || nrow(spec) == 1L) {
    tr <- if (is.null(trait)) spec$trait[[1L]] else trait
    if (!(tr %in% spec$trait)) return(NA_real_)
    v <- val(tr)
    return(if (length(v) && is.finite(v)) v else NA_real_)
  }
  method <- as.character(if (is.null(multi_trait_meta$method)) "" else multi_trait_meta$method)
  if (!(method %in% c("weighted", "economic_index", "desired_gain"))) return(NA_real_)
  w <- multi_trait_meta$weights
  if (is.null(w) || !length(w)) return(NA_real_)
  tr <- intersect(spec$trait, names(w))
  if (!length(tr)) return(NA_real_)
  v <- vapply(tr, val, numeric(1))
  if (any(!is.finite(v))) return(NA_real_)
  sum(v * as.numeric(w[tr]))
}
```

In `ng_plot_priority_score_vs_kinship()`, add the formals `check_line = NULL`,
`check_label = NULL` and, immediately after the plot region is drawn (after the existing
`plot()` call, before the legend), add:

```r
  if (!is.null(check_line) && is.finite(check_line)) {
    graphics::abline(h = check_line, lty = 2, lwd = 2, col = "#B00020")
    graphics::mtext(sprintf("check: %s", if (is.null(check_label)) "reference" else check_label),
                    side = 4, at = check_line, las = 1, cex = 0.7, col = "#B00020")
  }
```

- [ ] **Step 4: Add the per-trait panel for multi-trait runs**

One y axis cannot carry N check lines honestly, so multi-trait checks get small multiples — one
facet per trait, each on its own scale with its own line.

Append to `tests/check_reference_plot.R` (before its final `cat`):

```r
# multi-trait: one panel per trait with a check, each with its own line
scored_mt <- data.frame(
  parent1 = c("P1", "P1"), parent2 = c("P2", "P3"),
  pair_kinship = c(0.1, 0.2),
  yield_mean = c(9, 4), yield_check_value = 6, yield_check_ok = c(TRUE, FALSE),
  protein_mean = c(11, 13), protein_check_value = 10, protein_check_ok = c(TRUE, TRUE),
  stringsAsFactors = FALSE)
ref_mt <- list(active = data.frame(trait = c("yield", "protein"), check = "CHK_A",
                                   reject_if = "below", stringsAsFactors = FALSE),
               values = list(yield = c(CHK_A = 6), protein = c(CHK_A = 10)),
               source = c(yield = "GEBV", protein = "GEBV"))
p2 <- file.path(tempdir(), "chk-panels.png")
ng_plot_check_panels(scored_mt, ref_mt, output_path = p2)
stopifnot(file.exists(p2), file.size(p2) > 0)

# a run with no checks draws nothing and returns NULL rather than erroring
stopifnot(is.null(ng_plot_check_panels(scored_mt, NULL, output_path = p2)))
```

Add to `R/38_cross_priority_plots.R`:

```r
# Per-trait check panels: one facet per trait that has a check, y = that trait's mid-parent
# mean, x = pair kinship, with the trait's own check line. This is where multi-trait checks
# live, because a single index axis cannot carry several check lines on different scales.
ng_plot_check_panels <- function(scored, trait_check_reference, output_path = NULL,
                                 kinship_col = "pair_kinship", width = 10, height = 4,
                                 res = 150) {
  if (is.null(trait_check_reference)) return(invisible(NULL))
  spec <- as.data.frame(trait_check_reference$active, stringsAsFactors = FALSE)
  key <- if ("column_key" %in% names(spec)) as.character(spec$column_key) else as.character(spec$trait)
  keep <- paste0(key, "_mean") %in% names(scored)
  traits <- spec$trait[keep]; key <- key[keep]
  if (!length(traits)) return(invisible(NULL))
  if (!is.null(output_path)) {
    grDevices::png(output_path, width = width, height = height, units = "in", res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  op <- graphics::par(mfrow = c(1L, length(traits)), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  for (i in seq_along(traits)) {
    tr <- traits[[i]]; kk <- key[[i]]
    y <- suppressWarnings(as.numeric(scored[[paste0(kk, "_mean")]]))
    x <- suppressWarnings(as.numeric(scored[[kinship_col]]))
    ok <- scored[[paste0(kk, "_check_ok")]]
    graphics::plot(x, y, pch = 19, col = ifelse(ok %in% FALSE, "#BBBBBB", "#1F4E78"),
                   xlab = "Pair kinship", ylab = paste(tr, "mid-parent"), main = tr)
    tau <- suppressWarnings(as.numeric(scored[[paste0(kk, "_check_value")]][[1L]]))
    if (length(tau) && is.finite(tau)) {
      graphics::abline(h = tau, lty = 2, lwd = 2, col = "#B00020")
    }
  }
  invisible(output_path)
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `Rscript tests/check_reference_plot.R`
Expected: `check plot ok`

Run: `Rscript tests/cross_priority_plots.R`
Expected: still passes.

- [ ] **Step 6: Wire the line into real runs**

`ng_run_cp_output_files()` (`R/39_cross_prediction_runner.R:609`, formal added in Task 8) is the
only call site of `ng_plot_priority_score_vs_kinship()` (:630). Resolve the line there and pass
it, and emit the per-trait panel when more than one trait has a check:

```r
    chk_line <- if (is.null(trait_check_reference)) NA_real_ else
      ng_check_line_value(trait_check_reference, multi_trait_meta)
    ng_plot_priority_score_vs_kinship(
      scored = candidate_crosses,
      selected = selected_crosses,
      output_path = plot_path,
      check_line = if (is.finite(chk_line)) chk_line else NULL,
      check_label = if (is.finite(chk_line)) trait_check_reference$active$check[[1L]] else NULL
    )
    if (!is.null(trait_check_reference) && nrow(trait_check_reference$active) > 1L) {
      panel_path <- file.path(output_dir, "check_panels.png")
      ng_plot_check_panels(candidate_crosses, trait_check_reference, output_path = panel_path)
      files$check_panels_png <- normalizePath(panel_path, winslash = "/", mustWork = TRUE)
      figures <- rbind(figures, data.frame(figure = "check_panels",
                                           path = files$check_panels_png,
                                           stringsAsFactors = FALSE))
    }
```

`multi_trait_meta` is not currently available inside `ng_run_cp_output_files()`. Add it as a
formal (`multi_trait_meta = NULL`) and pass the runner's `ctx$multi_trait_meta` at the call
site — it is captured at `R/39_cross_prediction_runner.R:1391` region as
`ctx$multi_trait_meta <- attr(scored_crosses, "multi_trait")`.

- [ ] **Step 7: Commit**

```bash
git add R/38_cross_priority_plots.R R/39_cross_prediction_runner.R tests/check_reference_plot.R
git commit -m "feat(checks): check reference line and per-trait check panels"
```

---

### Task 10: Docs, exports, version bump, full check

**Files:**
- Modify: `NAMESPACE`, `man/nextgenCrossDesign-api.Rd`, `NEWS.md`, `DESCRIPTION`
- Modify: `docs/BACKEND_USER_GUIDE.md`

- [ ] **Step 1: Export the new API**

Add to `NAMESPACE`:

```
export(ng_align_check_geno)
export(ng_attach_check_reference)
export(ng_attach_joint_check_probability)
export(ng_check_line_value)
export(ng_check_reference_value)
export(ng_check_tau_bounds)
export(ng_plot_check_panels)
```

Add matching `\alias{}` entries to `man/nextgenCrossDesign-api.Rd`.

- [ ] **Step 2: Write the NEWS entry**

Prepend to `NEWS.md`:

```markdown
# nextgenCrossDesign 0.23.0

## Breaking

* Check lines are now **references, not filters**. `exclude_threshold_violators` and
  `check_basis` are removed, and checks no longer drop candidate crosses. A check is supplied
  in its own `check_geno` matrix, separate from the candidate parents, and is never crossed
  (check x check would be a self, which the package does not enumerate).
* `ng_apply_trait_checks()` is removed; use `ng_attach_check_reference()`.
* `ng_trait_check_spec()` loses its `basis` argument: the check always follows the run's
  per-trait `mean_source`, which is what keeps the reference on the same scale as the cross
  means it is compared against.

## New

* `check_geno`, `check_records`, `check_progeny_size` on `ng_run_cross_prediction()`.
  `check_progeny_size` is **required** whenever `trait_checks` is supplied and has no default:
  progeny per family scales P(beat check) directly, so it must be the breeding program's own
  figure rather than a number the package invents.
* Per-trait reference columns on the candidate table: `<trait>_check_value`, `<trait>_vs_check`
  (direction-aware: positive always means better), `<trait>_check_ok`,
  `<trait>_p_beat_check`, plus `checks_all_ok` and, for multi-trait runs,
  `p_beat_all_checks` — the probability that a progeny beats every check at once.
* `Checks` sheet in the cross-priority workbook.
* A horizontal check reference line on the score-versus-kinship plot, and
  `ng_plot_check_panels()` for per-trait small multiples in multi-trait runs.

## Bug fixes

* `ng_p_superior_progeny()` compared the wrong threshold for zero-variance crosses when `tau`
  varied per cross. `mu` was subset to the zero-variance rows while `tau` was not, so each such
  cross was scored against whichever threshold sat at position 1 rather than its own — silently
  returning a wrong probability at exactly the crosses whose answer is unambiguous. Scalar
  `tau`, which every in-package caller passes, was unaffected.
```

- [ ] **Step 3: Bump the version**

In `DESCRIPTION`, set `Version: 0.23.0`.

- [ ] **Step 4: Run the full suite and R CMD check**

```bash
Rscript -e 'for (f in list.files("tests", "\\.R$", full.names = TRUE)) { cat("--", f, "\n"); source(f) }'
R CMD build . && R CMD check --no-manual nextgenCrossDesign_0.23.0.tar.gz
```

Expected: every test script prints its `ok` line; `R CMD check` reports 0 errors, 0 warnings.

- [ ] **Step 5: Commit**

```bash
git add NAMESPACE man/ NEWS.md DESCRIPTION docs/
git commit -m "release: nextgenCrossDesign 0.23.0 — checks as references, not filters"
```
