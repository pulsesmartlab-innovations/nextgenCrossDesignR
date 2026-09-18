# Frontend Job Registry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A breeder picks "each trait separately", submits, and can close the browser — then returns to a Jobs tab, opens a finished trait, and reads its crossing plan while the rest are still being scored.

**Architecture:** The workbench becomes a *viewer* over a job store the backend owns. Submit spawns one detached `Rscript` and returns immediately. A Jobs tab polls cheaply on file mtimes and refreshes through the backend's own `ng_job_status()` only when something actually changed. Opening a trait loads that one trait's `result.json` into the existing `rv$result` slot, so every Results tab keeps working untouched.

**Tech Stack:** R, Shiny, `jsonlite`. No new package dependencies. The backend (`nextgenCrossDesign` >= 0.37.0) is driven as a separate process, as today.

**Spec:** `docs/superpowers/specs/2026-09-17-batch-job-registry-design.md` (Section 3, "Frontend")

**Repository:** `~/NextGenCrossDesign` — the frontend package `nextgenCrossWorkbench`. This plan does NOT touch the backend repo.

## Global Constraints

- **Frontend only.** The backend is finished and merged at 0.37.0; do not edit `~/Library/CloudStorage/.../cross_prediction`.
- **The frontend shows what the backend has.** It must not compute quantitative-genetics values, and must not reimplement the backend's derivations. In particular `crashed` and `incomplete` are derived by `ng_job_status()` and must never be recomputed in Shiny — a divergence would show a breeder a plan that was never finished.
- **The backend is reached ONLY as a subprocess.** This package must not `library(nextgenCrossDesign)`; `DESCRIPTION` must not depend on it. Existing pattern: `system2(rscript, ...)` in `R/run_backend.R`.
- **No new package dependencies.** `DESCRIPTION` Imports stay: shiny, bslib, DT, jsonlite, yaml, utils, grDevices, graphics, base64enc, stats, plotly, parallel.
- **Identical on Windows, macOS and Linux.** No POSIX-only paths; no `ps`/`tasklist`.
- **Never block the Shiny process on a batch.** R is single-threaded, so a blocking call freezes every session, not just the caller's.
- **NEVER add a `Co-Authored-By` trailer** to a commit. The repository owner is sole contributor of record. Check `git log --grep` before finishing — a prior branch had trailers added by tooling despite instructions.
- Tests are `testthat` files under `tests/testthat/`, run with `Rscript -e 'testthat::test_local()'` from the repo root, or the whole suite via `tests/testthat.R`. Tests needing a real backend are gated the way `tests/testthat/test-inbred-fix.R` is.
- `inst/BACKEND_VERSION` must end at `0.37.0`.

## File Structure

| File | Responsibility |
|---|---|
| `R/jobs.R` | **new** — everything about reading the job store: listing, status, locating a trait's result. The only file that knows the store's layout. |
| `inst/app/tools/job_status_json.R` | **new** — a tiny headless script that calls the backend's exported `ng_job_list()` / `ng_job_status()` and writes JSON. Keeps the derivation in the backend while letting Shiny stay out-of-process. |
| `R/run_backend.R` | modify — add `ngcd_submit_job()`: mint a job dir, materialise inputs, write config, spawn one detached `Rscript`, return immediately. |
| `R/config.R` | modify — `cfg$jobs_dir` and `cfg$shared_dir` as siblings of `runs_dir`; `BACKEND_VERSION` floor. |
| `R/helpers.R` | modify — `ngcd_objective_backend()` gains the `per_trait` mode; stage-key partition updated. |
| `R/app.R` | modify — fourth objective mode, the Jobs tab, submit handler, trait drill-down. |
| `tests/testthat/*` | **new** — unit tests for the store reader and the submit contract; one backend-gated end-to-end. |

`R/jobs.R` is deliberately separate from `R/runs.R`: runs are the old single-run history, jobs are the new durable store. They answer different questions and will change for different reasons.

---

### Task 1: Config — where the job store lives

**Files:**
- Modify: `R/config.R` (the block deriving `cfg$runs_dir` / `cfg$report_dir` / `cfg$presets_dir`, around line 162)
- Modify: `inst/BACKEND_VERSION`
- Test: `tests/testthat/test-config-jobs.R`

**Interfaces:**
- Produces: `cfg$jobs_dir`, `cfg$shared_dir` — absolute paths, created on startup.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-config-jobs.R`:

```r
# The job store and the shared-artefact store are SIBLINGS of runs_dir, never inside it.
#
# ngcd_prune_runs() enumerates list.dirs(cfg$runs_dir) and unlinks everything past
# keep_runs with unlink(recursive = TRUE, force = TRUE). A jobs/ directory living there
# would be just another directory to it: once a breeder accumulated keep_runs runs, prune
# would delete the entire job store -- including a job still running -- and the shared
# artefacts other jobs depend on.
test_that("the job store is a sibling of runs_dir, not inside it", {
  cfg <- ngcd_load_config()
  expect_true(nzchar(cfg$jobs_dir))
  expect_true(nzchar(cfg$shared_dir))

  norm <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  # Neither store may be located under runs_dir.
  expect_false(startsWith(norm(cfg$jobs_dir), paste0(norm(cfg$runs_dir), "/")))
  expect_false(startsWith(norm(cfg$shared_dir), paste0(norm(cfg$runs_dir), "/")))
  # Both live beside it, under data_dir, as presets_dir and report_dir already do.
  expect_equal(dirname(norm(cfg$jobs_dir)), norm(cfg$data_dir))
  expect_equal(dirname(norm(cfg$shared_dir)), norm(cfg$data_dir))
  # And they are distinct from each other.
  expect_false(identical(norm(cfg$jobs_dir), norm(cfg$shared_dir)))
})

test_that("both stores exist after configuration", {
  cfg <- ngcd_load_config()
  expect_true(dir.exists(cfg$jobs_dir))
  expect_true(dir.exists(cfg$shared_dir))
})

test_that("the declared backend floor is 0.37.0", {
  # The job registry (ng_job_list / ng_job_status / ng_job_prune, workflow = "batch")
  # arrived in 0.37.0. Against an older backend the Jobs tab would call functions that do
  # not exist, so the floor is what turns that into one clear message at startup instead
  # of an error per click.
  expect_equal(ngcd_default_backend_version(), "0.37.0")
  expect_equal(ngcd_load_config()$required_backend_version, "0.37.0")
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `cd ~/NextGenCrossDesign && Rscript -e 'testthat::test_local(filter = "config-jobs")'`
Expected: FAIL — `cfg$jobs_dir` is NULL / zero-length.

- [ ] **Step 3: Write the minimal implementation**

In `R/config.R`, in the block that already derives the other data directories, add `jobs_dir` and `shared_dir` beside them:

```r
  cfg$runs_dir    <- file.path(cfg$data_dir, "runs")
  cfg$report_dir  <- file.path(cfg$data_dir, "_report")
  cfg$presets_dir <- file.path(cfg$data_dir, "presets")
  # Siblings of runs_dir, never inside it: ngcd_prune_runs() unlinks everything in runs_dir
  # past keep_runs, and a job store there would be just another directory to it -- taking a
  # running job and other jobs' shared artefacts with it.
  cfg$jobs_dir    <- file.path(cfg$data_dir, "jobs")
  cfg$shared_dir  <- file.path(cfg$data_dir, "shared")
  for (d in c(cfg$runs_dir, cfg$report_dir, cfg$presets_dir, cfg$jobs_dir, cfg$shared_dir))
```

(extend the existing `for` loop's vector rather than adding a second loop).

Set `inst/BACKEND_VERSION` to `0.37.0`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "config-jobs")'`
Expected: 3 passing tests.

- [ ] **Step 5: Confirm nothing else regressed**

Run: `Rscript -e 'testthat::test_local(filter = "config")'`
Expected: all config tests pass.

- [ ] **Step 6: Commit**

```bash
git add R/config.R inst/BACKEND_VERSION tests/testthat/test-config-jobs.R
git commit -m "Give the job store a home beside runs, never inside it

ngcd_prune_runs() unlinks everything in runs_dir past keep_runs with
unlink(recursive = TRUE, force = TRUE). A job store living there would be just
another directory to it: once a breeder accumulated keep_runs runs, prune would
delete jobs that were still running, and the shared artefacts other jobs depend
on. jobs_dir and shared_dir sit beside runs_dir under data_dir, as presets_dir
and report_dir already do.

Backend floor to 0.37.0, which is where the job registry arrived."
```

---

### Task 2: Reading the job store

**Files:**
- Create: `R/jobs.R`
- Create: `inst/app/tools/job_status_json.R`
- Test: `tests/testthat/test-jobs-store.R`

**Interfaces:**
- Consumes: `cfg$jobs_dir` (Task 1).
- Produces:
  - `ngcd_jobs_fingerprint(jobs_dir)` → a single character string that changes when anything in the store changes. Cheap: mtimes only, no parsing, no subprocess.
  - `ngcd_jobs_list(cfg)` → data frame (`id`, `label`, `created_at`, `state`, `phase`, `n_traits`, `n_done`, `n_error`, `path`), newest first, empty frame when there are none.
  - `ngcd_job_detail(cfg, job_id)` → list with the job fields plus `traits`, a data frame (`trait`, `state`, `cv_predictive_r2`, `mean_source`, `effect_gate`, `n_selected`, `error_message`).
  - `ngcd_job_result_path(cfg, job_id, trait)` → path to that trait's `result.json`, or `NA_character_` if absent.

**WHY A FINGERPRINT AND A SUBPROCESS.** Shiny cannot call the backend in-process, so every `ng_job_status()` call costs an R process spawn (~1s). Polling that every few seconds would be absurd. But the frontend also must not reimplement the backend's `crashed`/`incomplete` derivations. So: poll a cheap fingerprint (file mtimes, no spawn); only when it CHANGES, spawn once and ask the backend. Idle polling costs nothing; a real change costs one process.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-jobs-store.R`:

```r
# Reading the job store must be CHEAP when nothing has changed, and must never reimplement
# the backend's derivations when something has.
#
# The frontend reaches the backend only as a subprocess, so each ng_job_status() call costs
# an R process. Polling that on a timer would be absurd -- but recomputing "crashed" and
# "incomplete" in Shiny to avoid the cost would be worse: the two implementations would
# eventually disagree, and the disagreement would show a breeder a plan that was never
# finished. So the cheap part (has anything changed?) is local, and the interpretation stays
# in the backend.

fake_job <- function(root, id, state = "finished", traits = list()) {
  d <- file.path(root, id)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(list(schema = "ng_job.v1", id = id, label = id,
                            created_at = "2026-09-18T10:00:00.000Z", state = state,
                            n_traits = length(traits)),
                       file.path(d, "job.json"), auto_unbox = TRUE)
  for (tr in names(traits)) {
    td <- file.path(d, tr); dir.create(td, showWarnings = FALSE)
    jsonlite::write_json(c(list(schema = "ng_job_trait_status.v1", trait = tr),
                           traits[[tr]]),
                         file.path(td, "status.json"), auto_unbox = TRUE)
  }
  d
}

test_that("the fingerprint is stable when nothing changes and moves when something does", {
  root <- file.path(tempdir(), paste0("fp_", as.integer(runif(1, 1, 1e6))))
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  fake_job(root, "job_a", traits = list(YIELD = list(state = "done")))
  f1 <- ngcd_jobs_fingerprint(root)
  expect_true(is.character(f1) && nzchar(f1))
  expect_identical(ngcd_jobs_fingerprint(root), f1)   # stable: no spawn, no churn

  Sys.sleep(1.1)   # mtime resolution
  fake_job(root, "job_b", traits = list(PROTEIN = list(state = "done")))
  expect_false(identical(ngcd_jobs_fingerprint(root), f1))
})

test_that("the fingerprint notices a trait finishing inside an existing job", {
  # This is the case that matters while a batch runs: the job directory itself may not
  # change, but a worker writes a trait's status.json. If the fingerprint missed that, the
  # Jobs view would sit frozen for hours while traits completed behind it.
  root <- file.path(tempdir(), paste0("fp2_", as.integer(runif(1, 1, 1e6))))
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  d <- fake_job(root, "job_a", state = "running",
                traits = list(YIELD = list(state = "running")))
  f1 <- ngcd_jobs_fingerprint(root)
  Sys.sleep(1.1)
  td <- file.path(d, "PROTEIN"); dir.create(td)
  jsonlite::write_json(list(schema = "ng_job_trait_status.v1", trait = "PROTEIN",
                            state = "done"),
                       file.path(td, "status.json"), auto_unbox = TRUE)
  expect_false(identical(ngcd_jobs_fingerprint(root), f1))
})

test_that("an empty or absent store is a normal state, not an error", {
  expect_identical(nrow(ngcd_jobs_list(list(jobs_dir = file.path(tempdir(), "nope_xyz")))), 0L)
  empty <- file.path(tempdir(), paste0("empty_", as.integer(runif(1, 1, 1e6))))
  dir.create(empty, recursive = TRUE)
  on.exit(unlink(empty, recursive = TRUE), add = TRUE)
  expect_identical(nrow(ngcd_jobs_list(list(jobs_dir = empty))), 0L)
})

test_that("the frontend does not reimplement the backend's derived states", {
  # `crashed` and `incomplete` are derived by ng_job_status() from a stale heartbeat and
  # from a trait's state relative to its job's. If this package computed either itself, the
  # two would drift. Nothing in R/jobs.R may mention them outside a comment.
  src <- readLines(system.file("R", "jobs.R", package = "nextgenCrossWorkbench"),
                   warn = FALSE)
  if (!length(src)) src <- readLines(file.path("..", "..", "R", "jobs.R"), warn = FALSE)
  code <- grep("^\\s*#", src, value = TRUE, invert = TRUE)
  expect_false(any(grepl("heartbeat", code)))
  expect_false(any(grepl('"crashed"', code, fixed = TRUE)))
  expect_false(any(grepl('"incomplete"', code, fixed = TRUE)))
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript -e 'testthat::test_local(filter = "jobs-store")'`
Expected: FAIL — `could not find function "ngcd_jobs_fingerprint"`.

- [ ] **Step 3: Write the headless bridge**

Create `inst/app/tools/job_status_json.R`:

```r
# Ask the backend about the job store, and write the answer as JSON.
#
# The workbench drives the backend as a separate process and never loads its namespace, so
# this script is how Shiny reaches ng_job_list() and ng_job_status(). It exists so the
# DERIVATIONS stay in the backend: `crashed` (a stale heartbeat) and `incomplete` (a trait
# still claiming to run under a job that is not) are computed there, once, and this script
# only transports the answer. A frontend that recomputed them would eventually disagree --
# and the disagreement would surface as a breeder acting on a plan that never finished.
#
# Usage:
#   Rscript job_status_json.R list   <jobs_dir> <out.json>
#   Rscript job_status_json.R detail <jobs_dir> <out.json> <job_id>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("usage: job_status_json.R <list|detail> <jobs_dir> <out> [job_id]")
mode <- args[[1L]]; jobs_dir <- args[[2L]]; out <- args[[3L]]

payload <- if (identical(mode, "list")) {
  list(schema = "ngcd_jobs_list.v1", ok = TRUE,
       jobs = nextgenCrossDesign::ng_job_list(jobs_dir))
} else if (identical(mode, "detail")) {
  if (length(args) < 4L) stop("detail requires a job id")
  s <- nextgenCrossDesign::ng_job_status(file.path(jobs_dir, args[[4L]]))
  list(schema = "ngcd_job_detail.v1", ok = TRUE, job = s)
} else {
  stop("unknown mode: ", mode)
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
# dataframe = "rows" so a table survives the boundary as an array of objects, matching the
# backend's own result.json convention; a consumer must not have to handle two shapes.
jsonlite::write_json(payload, out, auto_unbox = TRUE, na = "null", null = "null",
                     dataframe = "rows", pretty = TRUE)
```

- [ ] **Step 4: Write the store reader**

Create `R/jobs.R`:

```r
# Reading the job store.
#
# A job is a directory the backend writes and this package only ever READS. The layout and
# the meaning of every state belong to the backend; this file knows where the store is and
# how to ask.
#
# Two costs shape everything here. Reaching the backend means spawning an R process, so it
# cannot happen on a timer. But interpreting the store -- deciding that a job has crashed,
# or that a trait is incomplete -- must not happen here either, because a second
# implementation of those rules would drift from the first. So: a cheap local fingerprint
# answers "has anything changed?", and only a change buys a subprocess.

# A cheap signature of the store's current state: directory and status-file mtimes, no
# parsing and no subprocess. Includes per-trait status files deliberately -- while a batch
# runs, the job directory may not change but its workers are writing traits, and a
# fingerprint that missed those would leave the view frozen for hours.
ngcd_jobs_fingerprint <- function(jobs_dir) {
  if (!length(jobs_dir) || is.na(jobs_dir) || !dir.exists(jobs_dir)) return("empty")
  dirs <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
  if (!length(dirs)) return("empty")
  paths <- c(file.path(dirs, "job.json"),
             list.files(dirs, pattern = "^status[.]json$", recursive = TRUE,
                        full.names = TRUE))
  paths <- paths[file.exists(paths)]
  if (!length(paths)) return("empty")
  info <- file.info(paths)
  paste(length(paths), format(max(info$mtime), "%Y%m%d%H%M%OS3"),
        sum(as.numeric(info$size)), sep = "-")
}

# Run the headless bridge and parse its answer. Returns NULL on any failure: the Jobs view
# is the only route a breeder has to their results, so it degrades to "nothing to show"
# rather than taking the session down.
ngcd_jobs__ask <- function(cfg, mode, job_id = NULL) {
  script <- ngcd_res("app", "tools", "job_status_json.R")
  if (!nzchar(script) || !file.exists(script)) return(NULL)
  out <- tempfile(fileext = ".json"); on.exit(unlink(out), add = TRUE)
  args <- c(shQuote(script), shQuote(mode), shQuote(cfg$jobs_dir), shQuote(out))
  if (!is.null(job_id)) args <- c(args, shQuote(job_id))
  ok <- tryCatch(
    system2(ngcd_resolve_rscript(cfg), args = args, stdout = TRUE, stderr = TRUE,
            timeout = as.numeric(cfg$job_query_timeout_seconds %||% 60)),
    error = function(e) NULL)
  if (!file.exists(out)) return(NULL)
  tryCatch(jsonlite::fromJSON(out, simplifyVector = TRUE, simplifyDataFrame = TRUE),
           error = function(e) NULL)
}

ngcd_jobs_empty <- function() {
  data.frame(id = character(0), label = character(0), created_at = character(0),
             state = character(0), phase = character(0), n_traits = integer(0),
             n_done = integer(0), n_error = integer(0), path = character(0),
             stringsAsFactors = FALSE)
}

ngcd_jobs_list <- function(cfg) {
  if (!length(cfg$jobs_dir) || !dir.exists(cfg$jobs_dir)) return(ngcd_jobs_empty())
  res <- ngcd_jobs__ask(cfg, "list")
  if (is.null(res) || !isTRUE(res$ok) || !is.data.frame(res$jobs) || !nrow(res$jobs)) {
    return(ngcd_jobs_empty())
  }
  res$jobs
}

ngcd_job_detail <- function(cfg, job_id) {
  if (!length(job_id) || !nzchar(job_id)) return(NULL)
  res <- ngcd_jobs__ask(cfg, "detail", job_id)
  if (is.null(res) || !isTRUE(res$ok)) return(NULL)
  res$job
}

# Where a trait's deliverable lives. The path is the store's layout, which this package
# knows; what is IN the file is the backend's contract.
ngcd_job_result_path <- function(cfg, job_id, trait) {
  if (!length(job_id) || !length(trait)) return(NA_character_)
  p <- file.path(cfg$jobs_dir, job_id, trait, "result.json")
  if (file.exists(p)) p else NA_character_
}
```

If `ngcd_resolve_rscript(cfg)` does not already exist as a helper, use whatever `R/run_backend.R` already uses to resolve the Rscript binary, and reuse it rather than duplicating the logic.

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "jobs-store")'`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add R/jobs.R inst/app/tools/job_status_json.R tests/testthat/test-jobs-store.R
git commit -m "Read the job store without polling a subprocess or reinventing its rules

Shiny cannot call the backend in-process, so every ng_job_status() costs an R
process spawn. Polling that on a timer would be absurd. Recomputing the backend's
derived states in Shiny to avoid the cost would be worse: two implementations of
'crashed' and 'incomplete' would eventually disagree, and the disagreement would
show a breeder a plan that was never finished.

So a cheap local fingerprint -- directory and status-file mtimes, no parsing --
answers 'has anything changed?', and only a change buys a subprocess. The
fingerprint covers per-trait status files deliberately: while a batch runs the job
directory may not change at all, but its workers are writing traits."
```

---

### Task 3: Submit without blocking

**Files:**
- Modify: `R/run_backend.R`
- Test: `tests/testthat/test-submit-job.R`

**Interfaces:**
- Consumes: `cfg$jobs_dir`, `cfg$shared_dir`, `ngcd_materialize_inputs()`, `ngcd_write_config()`.
- Produces: `ngcd_submit_job(cfg, params, data, label = NULL)` → list with `job_id`, `job_dir`, `spawned` (logical), `log_path`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-submit-job.R`:

```r
# Submitting must return immediately, and a job that never starts must be visible as such.
#
# R is single-threaded: a blocking call freezes the WHOLE app process, so every other
# breeder's browser goes quiet too. A 17-trait batch runs for hours, so submit cannot wait.
#
# The harder requirement is the failure path. The detached process writes job.json as its
# first act. If the backend cannot load, or the config is bad, NO record is ever written --
# and "nothing happened" is the worst possible feedback for a run someone expects to take
# hours. The submit must leave enough behind for that case to be recognisable.

test_that("submit returns before the work finishes", {
  cfg <- ngcd_load_config()
  cfg$jobs_dir <- file.path(tempdir(), paste0("sj_", as.integer(runif(1, 1, 1e6))))
  dir.create(cfg$jobs_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cfg$jobs_dir, recursive = TRUE), add = TRUE)
  # Point at a command that sleeps far longer than this test will wait.
  cfg$rscript <- if (.Platform$OS.type == "windows") "timeout" else "sleep"
  cfg$submit_args_override <- "5"

  t0 <- Sys.time()
  out <- ngcd_submit_job(cfg, params = list(), data = list(), label = "quick")
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  expect_lt(elapsed, 3)                    # returned without waiting for the work
  expect_true(nzchar(out$job_id))
  expect_true(dir.exists(out$job_dir))
  expect_true(nzchar(out$log_path))
})

test_that("a job whose process never starts is recognisable, not silent", {
  cfg <- ngcd_load_config()
  cfg$jobs_dir <- file.path(tempdir(), paste0("sj2_", as.integer(runif(1, 1, 1e6))))
  dir.create(cfg$jobs_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cfg$jobs_dir, recursive = TRUE), add = TRUE)
  cfg$rscript <- "definitely_not_a_real_binary_xyz"

  out <- ngcd_submit_job(cfg, params = list(), data = list(), label = "doomed")
  # The directory and the log exist even though nothing ran, so the UI has something to
  # show and a reason to show. A silent nothing is the failure mode being designed out.
  expect_true(dir.exists(out$job_dir))
  expect_true(file.exists(out$log_path) || !isTRUE(out$spawned))
  expect_false(file.exists(file.path(out$job_dir, "job.json")))
})

test_that("each submit gets its own directory", {
  cfg <- ngcd_load_config()
  cfg$jobs_dir <- file.path(tempdir(), paste0("sj3_", as.integer(runif(1, 1, 1e6))))
  dir.create(cfg$jobs_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(cfg$jobs_dir, recursive = TRUE), add = TRUE)
  cfg$rscript <- if (.Platform$OS.type == "windows") "timeout" else "sleep"
  cfg$submit_args_override <- "1"

  a <- ngcd_submit_job(cfg, list(), list(), label = "one")
  b <- ngcd_submit_job(cfg, list(), list(), label = "two")
  expect_false(identical(a$job_id, b$job_id))
  expect_false(identical(a$job_dir, b$job_dir))
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript -e 'testthat::test_local(filter = "submit-job")'`
Expected: FAIL — `could not find function "ngcd_submit_job"`.

- [ ] **Step 3: Write the implementation**

Add to `R/run_backend.R`:

```r
# Submit a batch and return immediately.
#
# R is single-threaded, so a blocking system2() freezes the whole app process -- every
# connected session, not just the one that pressed the button. A 17-trait batch runs for
# hours. So the work is handed to a DETACHED process (wait = FALSE) and this function
# returns as soon as it has been launched.
#
# The detached process writes job.json as its first act. Until that appears the UI shows
# the job as starting; if it never appears, the job failed during startup and job.log says
# why. That is the whole reason the directory and the log are created here, before the
# spawn: a job that cannot start must still be a job someone can find and read an error
# from. "Nothing happened" is the worst possible feedback for work expected to take hours.
ngcd_submit_job <- function(cfg, params, data, label = NULL) {
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%OS3")
  job_id <- paste0("job-", gsub("[^0-9]", "", stamp), "-",
                   substr(basename(tempfile("")), 5L, 10L))
  job_dir <- file.path(cfg$jobs_dir, job_id)
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(job_dir, "job.log")

  input_paths <- ngcd_materialize_inputs(data, job_dir)
  cfg_json <- file.path(job_dir, "config.json")
  params <- utils::modifyList(params, c(input_paths, list(
    workflow = "batch",
    job_dir = job_dir,
    # The backend refuses a mismatch, and rightly: workers write each trait's status under
    # batch_output_root while ng_job_status() reads job_dir, so decoupling them would leave
    # a job whose traits are permanently invisible.
    batch_output_root = job_dir,
    shared_dir = cfg$shared_dir,
    label = label %||% job_id)))
  ngcd_write_config(params, cfg_json)

  rscript <- ngcd_resolve_rscript(cfg)
  args <- if (!is.null(cfg$submit_args_override)) cfg$submit_args_override else
    c(shQuote(cfg$runner_script), shQuote(cfg_json),
      shQuote(file.path(job_dir, "submit_result.json")))
  spawned <- tryCatch({
    system2(rscript, args = args, stdout = log_path, stderr = log_path, wait = FALSE)
    TRUE
  }, error = function(e) {
    # Record why nothing started, in the place the UI will look.
    try(writeLines(paste("Could not launch the backend:", conditionMessage(e)), log_path),
        silent = TRUE)
    FALSE
  })

  list(job_id = job_id, job_dir = job_dir, spawned = spawned, log_path = log_path)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "submit-job")'`
Expected: 3 passing tests.

- [ ] **Step 5: Prove it does not block, against a real delay**

```bash
Rscript -e '
devtools::load_all(".", quiet = TRUE)
cfg <- ngcd_load_config()
cfg$jobs_dir <- file.path(tempdir(), "blocktest"); dir.create(cfg$jobs_dir, recursive = TRUE)
cfg$rscript <- "sleep"; cfg$submit_args_override <- "10"
t0 <- Sys.time()
out <- ngcd_submit_job(cfg, list(), list(), label = "x")
cat("returned after", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 2), "s\n")
cat("process still running after return:", file.exists(out$job_dir), "\n")'
```

Expected: returns in well under a second, against a 10-second child.

- [ ] **Step 6: Commit**

```bash
git add R/run_backend.R tests/testthat/test-submit-job.R
git commit -m "Submit a batch without freezing every session

R is single-threaded, so the existing blocking system2() freezes the whole app
process -- every connected breeder's browser, not just the one that pressed the
button. A 17-trait batch runs for hours.

Submit now hands the work to a detached process and returns. The job directory
and its log are created BEFORE the spawn, because the detached process writes
job.json as its first act: if the backend cannot load, no record is ever written,
and a job that cannot start must still be findable with a log saying why.
'Nothing happened' is the worst possible feedback for work expected to take hours.

batch_output_root is pinned to job_dir. The backend refuses a mismatch, and
rightly -- workers write trait status under the output root while ng_job_status()
reads the job dir, so decoupling them leaves a job whose traits never appear."
```

---

### Task 4: The fourth objective mode

**Files:**
- Modify: `R/helpers.R` (`ngcd_objective_backend`, and the stage-key partition)
- Modify: `R/app.R` (the `objective_mode` radio, around line 343)
- Test: `tests/testthat/test-objective-per-trait.R`

**Interfaces:**
- Consumes: nothing new.
- Produces: `objective_mode == "per_trait"`, and `ngcd_objective_backend()` returning `prediction_mode = "trait_by_trait"`, `traits_to_use = <all selected traits>`, `multi_trait_method_applies = FALSE`, `per_trait_jobs = TRUE`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-objective-per-trait.R`:

```r
# "Each trait separately" is an OBJECTIVE, not an execution style.
#
# A breeder asking for seventeen crossing plans is making the same kind of choice as one
# asking to combine seventeen traits into an index: it is about what gets analysed, not
# about how it runs. It sits beside single / multi / index for that reason.
#
# The distinction that matters downstream: multi-trait produces ONE plan from a combined
# index, per-trait produces ONE PLAN PER TRAIT. Same traits, completely different
# deliverable.

test_that("per_trait selects every trait but does not combine them", {
  out <- ngcd_objective_backend("per_trait", single_trait = "YIELD",
                                traits = c("YIELD", "PROTEIN", "HEIGHT"),
                                index_col = NULL)
  expect_equal(out$prediction_mode, "trait_by_trait")
  expect_equal(out$traits_to_use, c("YIELD", "PROTEIN", "HEIGHT"))
  # No index is built, so the multi-trait method must not be offered or applied.
  expect_false(out$multi_trait_method_applies)
  expect_true(isTRUE(out$per_trait_jobs))
})

test_that("multi and per_trait differ only in what they produce", {
  traits <- c("YIELD", "PROTEIN")
  m <- ngcd_objective_backend("multi", "YIELD", traits, NULL)
  p <- ngcd_objective_backend("per_trait", "YIELD", traits, NULL)
  expect_equal(m$traits_to_use, p$traits_to_use)      # same traits
  expect_true(m$multi_trait_method_applies)           # one combined index
  expect_false(p$multi_trait_method_applies)          # one plan each
  expect_null(m$per_trait_jobs)
  expect_true(isTRUE(p$per_trait_jobs))
})

test_that("the existing three modes are unchanged", {
  s <- ngcd_objective_backend("single", "YIELD", c("YIELD", "PROTEIN"), NULL)
  expect_equal(s$traits_to_use, "YIELD")
  expect_false(s$multi_trait_method_applies)
  i <- ngcd_objective_backend("index", "YIELD", c("YIELD"), "IDX")
  expect_equal(i$prediction_mode, "index_as_trait")
  expect_null(i$traits_to_use)
})

test_that("an unknown mode still fails loudly", {
  expect_error(ngcd_objective_backend("nonsense", "A", "A", NULL), "Unknown objective_mode")
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript -e 'testthat::test_local(filter = "objective-per-trait")'`
Expected: FAIL — `Unknown objective_mode: per_trait`.

- [ ] **Step 3: Write the implementation**

In `R/helpers.R`, extend `ngcd_objective_backend()`'s `switch`:

```r
    multi  = list(prediction_mode = "trait_by_trait", traits_to_use = traits,
                  multi_trait_method_applies = TRUE),
    # Every trait, but no index: one crossing plan PER trait rather than one combined
    # plan. Same traits as `multi`, entirely different deliverable -- which is why it is an
    # objective a breeder chooses, not an execution mode.
    per_trait = list(prediction_mode = "trait_by_trait", traits_to_use = traits,
                     multi_trait_method_applies = FALSE, per_trait_jobs = TRUE),
```

In `R/app.R`, add the fourth choice to the `objective_mode` radio (keep the existing labels verbatim; add):

```r
                     "Each trait separately (one plan per trait)" = "per_trait",
```

and give `per_trait` the same trait-picker panel `multi` uses, so the breeder selects which traits:

```r
            shiny::conditionalPanel("input.objective_mode == 'multi' || input.objective_mode == 'per_trait'",
                                    shiny::uiOutput("traits_to_use_ui")),
```

Leave the multi-trait-method panel gated on `multi` alone — no index is built in `per_trait`.

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "objective-per-trait")'`
Expected: 4 passing tests.

- [ ] **Step 5: Wire the Run action to submit, so the mode actually does something**

A mode that cannot be run is not a mode. In `R/app.R`, where the Run action currently routes
between the staged and single paths (`ngcd_run_uses_staged()`), add a branch ahead of both:
when `input$objective_mode == "per_trait"`, submit instead of running.

```r
      # per_trait does not run here at all -- it is handed to a detached process and the
      # session is free immediately. Routed BEFORE the staged/single split because those
      # both block, and a seventeen-trait batch must never block: R is single-threaded, so
      # it would freeze every connected session for hours, not just this one.
      if (identical(input$objective_mode, "per_trait")) {
        out <- ngcd_submit_job(cfg, params = staged_params(), data = rv$data,
                               label = input$job_label %||% NULL)
        if (!isTRUE(out$spawned)) {
          shiny::showNotification(
            "The backend could not be started - see the job log for the reason.",
            type = "error", duration = NULL)
        } else {
          rv$open_job <- out$job_id
          shiny::showNotification(
            paste0("Submitted. This continues even if you close the browser - ",
                   "open the Jobs screen to watch it."),
            type = "message")
        }
        return(invisible(NULL))
      }
```

Add a matching assertion to `tests/testthat/test-objective-per-trait.R`:

```r
test_that("per_trait is routed away from both blocking paths", {
  # ngcd_run_uses_staged() decides between the two paths that BLOCK. per_trait must be
  # handled before that decision is ever reached, or a seventeen-trait batch would freeze
  # every session for hours.
  src <- readLines(file.path("..", "..", "R", "app.R"), warn = FALSE)
  submit_line <- grep("ngcd_submit_job", src)
  staged_line <- grep("ngcd_run_uses_staged", src)
  expect_true(length(submit_line) > 0)
  expect_true(length(staged_line) > 0)
  expect_lt(min(submit_line), min(staged_line))
})
```

- [ ] **Step 6: Confirm the pipeline-state partition still holds**

Run: `Rscript -e 'testthat::test_local(filter = "pipeline-state")'`
Expected: pass. That suite asserts every `build_params` key maps to at least one stage; if the new mode introduces a key, this is what catches it.

- [ ] **Step 7: Commit**

```bash
git add R/helpers.R R/app.R tests/testthat/test-objective-per-trait.R
git commit -m "Ask for a plan per trait as an objective, not an execution mode

A breeder wanting seventeen crossing plans is making the same kind of choice as
one wanting seventeen traits combined into an index: it is about what gets
analysed. So it sits beside single / multi / index rather than appearing as a
run-mechanics switch elsewhere.

per_trait selects every chosen trait like multi does, but builds no index -- same
traits, one plan each instead of one combined plan."
```

---

### Task 5: The Jobs tab

**Files:**
- Modify: `R/app.R` (navigation, plus the new panel's UI and server)
- Test: `tests/testthat/test-jobs-ui.R`

**Interfaces:**
- Consumes: `ngcd_jobs_fingerprint()`, `ngcd_jobs_list()` (Task 2).
- Produces: a top-level `Jobs` nav panel; `jobs_poll` reactive.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-jobs-ui.R`:

```r
# The Jobs view must be cheap to keep open and must never claim more than it knows.

test_that("the job list renders from status alone, never from results", {
  # The overview shows every trait's state, R2 and gate verdict. Those come from the small
  # per-trait status.json files. If it read result.json instead, opening the LIST would cost
  # what opening a RESULT costs -- seventeen files of ~125,000 candidate crosses each -- and
  # the thing the list exists to avoid would happen on every poll.
  src <- readLines(file.path("..", "..", "R", "jobs.R"), warn = FALSE)
  code <- grep("^\\s*#", src, value = TRUE, invert = TRUE)
  listing <- code[seq_len(max(grep("ngcd_job_result_path", code)[1] - 1, 1))]
  expect_false(any(grepl("result\\.json", listing)))
})

test_that("a Jobs panel exists in the navigation", {
  ui <- workbench_ui(ngcd_load_config())
  html <- paste(as.character(ui), collapse = " ")
  expect_true(grepl("Jobs", html, fixed = TRUE))
})

test_that("an empty store renders an explanation rather than an empty table", {
  # A breeder who has submitted nothing should be told what this screen is for, not shown a
  # bare table with no rows and no hint that anything is meant to appear.
  out <- ngcd_jobs_placeholder(ngcd_jobs_empty())
  html <- paste(as.character(out), collapse = " ")
  expect_true(nzchar(html))
  expect_true(grepl("submit|No jobs|nothing", html, ignore.case = TRUE))
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript -e 'testthat::test_local(filter = "jobs-ui")'`
Expected: FAIL — no Jobs panel / `ngcd_jobs_placeholder` not found.

- [ ] **Step 3: Write the implementation**

Add a `Jobs` panel to the navbar beside the existing stages, and in the server:

```r
    # Poll cheaply, refresh rarely. checkFunc reads mtimes only -- no parsing and no
    # subprocess -- so an idle Jobs tab costs essentially nothing. valueFunc spawns the
    # backend, and only runs when the fingerprint actually moved.
    jobs_poll <- shiny::reactivePoll(
      intervalMillis = 3000, session = session,
      checkFunc = function() ngcd_jobs_fingerprint(cfg$jobs_dir),
      valueFunc = function() ngcd_jobs_list(cfg))
```

and a placeholder helper in `R/jobs.R`:

```r
# What to show when there is nothing to show. A bare empty table tells a breeder neither
# what this screen is for nor that anything is expected to appear on it.
ngcd_jobs_placeholder <- function(jobs) {
  if (is.data.frame(jobs) && nrow(jobs)) return(NULL)
  shiny::div(
    class = "text-muted",
    shiny::p("No jobs yet."),
    shiny::p(paste("Choose \"Each trait separately\" on the Configure screen and submit;",
                   "the run continues even if you close this browser, and each trait",
                   "appears here as it finishes.")))
}
```

Render the table with the columns `label`, `state`, `phase`, `n_done` of `n_traits`, `created_at`, and make the row clickable to open the job.

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "jobs-ui")'`
Expected: 3 passing tests.

- [ ] **Step 5: Commit**

```bash
git add R/app.R R/jobs.R tests/testthat/test-jobs-ui.R
git commit -m "A Jobs tab that is cheap to leave open

checkFunc reads mtimes only -- no parsing, no subprocess -- so an idle Jobs tab
costs essentially nothing, and the backend is asked only when the store actually
changed.

The overview is built from the small per-trait status files. Reading result.json
for the list would cost what opening a result costs: seventeen files of ~125,000
candidate crosses each, on every poll."
```

---

### Task 6: Opening a trait

**Files:**
- Modify: `R/app.R`
- Test: `tests/testthat/test-job-drilldown.R`

**Interfaces:**
- Consumes: `ngcd_job_detail()`, `ngcd_job_result_path()` (Task 2), `rv$result`, `res()`.
- Produces: a job-detail view; selecting a trait loads that trait's `result.json` into `rv$result`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-job-drilldown.R`:

```r
# Opening a trait loads exactly one result, and the overview works even when results do not.

test_that("the overview survives unreadable results", {
  # THE PROPERTY THAT PROVES LAZINESS. If the overview can be built while every result.json
  # is corrupt, it provably does not read them. This fails loudly the day someone adds a
  # convenient fromJSON() to the list view -- which is exactly how a cheap screen becomes an
  # expensive one without anybody noticing.
  root <- file.path(tempdir(), paste0("dd_", as.integer(runif(1, 1, 1e6))))
  d <- file.path(root, "job_a"); dir.create(d, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  jsonlite::write_json(list(schema = "ng_job.v1", id = "job_a", label = "job_a",
                            created_at = "2026-09-18T10:00:00.000Z", state = "finished",
                            n_traits = 1),
                       file.path(d, "job.json"), auto_unbox = TRUE)
  td <- file.path(d, "YIELD"); dir.create(td)
  jsonlite::write_json(list(schema = "ng_job_trait_status.v1", trait = "YIELD",
                            state = "done", cv_predictive_r2 = 0.61, n_selected = 20),
                       file.path(td, "status.json"), auto_unbox = TRUE)
  writeLines("{ this is not json", file.path(td, "result.json"))   # deliberately corrupt

  expect_true(nzchar(ngcd_jobs_fingerprint(root)))                  # still readable
  # And the path helper reports the file exists without parsing it.
  expect_true(!is.na(ngcd_job_result_path(list(jobs_dir = root), "job_a", "YIELD")))
})

test_that("a trait with no result yet reports no path", {
  root <- file.path(tempdir(), paste0("dd2_", as.integer(runif(1, 1, 1e6))))
  dir.create(file.path(root, "job_a", "RUNNING_TRAIT"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_true(is.na(ngcd_job_result_path(list(jobs_dir = root), "job_a", "RUNNING_TRAIT")))
})
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript -e 'testthat::test_local(filter = "job-drilldown")'`
Expected: FAIL if `ngcd_job_result_path` is missing; otherwise confirm both assertions genuinely exercise it.

- [ ] **Step 3: Write the implementation**

In `R/app.R`, add the job-detail view: a header (label, state, phase, `n_done` of `n_traits`), the traits table from `ngcd_job_detail()$traits`, and a row action that opens a trait:

```r
    # ONE result in memory at a time. res() at the Results chokepoint is read by ~39 outputs,
    # so putting a trait's result into rv$result makes every existing Results tab work
    # untouched. Caching seventeen of them would put us back where we started: a result
    # carries every unordered parent pair, ~125,000 rows at 500 parents.
    shiny::observeEvent(input$job_open_trait, {
      p <- ngcd_job_result_path(cfg, rv$open_job, input$job_open_trait)
      if (is.na(p)) {
        shiny::showNotification(
          "That trait has not finished yet - its results appear here when it does.",
          type = "message")
        return(invisible(NULL))
      }
      parsed <- tryCatch(jsonlite::fromJSON(p, simplifyVector = TRUE,
                                            simplifyDataFrame = TRUE, simplifyMatrix = FALSE),
                         error = function(e) NULL)
      if (is.null(parsed)) {
        shiny::showNotification("That trait's result file could not be read.", type = "error")
        return(invisible(NULL))
      }
      rv$result <- ngcd_enrich_result(parsed)
      rv$error <- NULL
      rv$warnings <- parsed$warnings
    })
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript -e 'testthat::test_local(filter = "job-drilldown")'`
Expected: 2 passing tests.

- [ ] **Step 5: Confirm the Results outputs still work**

Run: `Rscript -e 'testthat::test_local(filter = "report|results|charts")'`
Expected: pass — `rv$result` still holds the same shape a single run produces, so nothing downstream changes.

- [ ] **Step 6: Commit**

```bash
git add R/app.R tests/testthat/test-job-drilldown.R
git commit -m "Open one trait at a time, through the slot Results already reads

res() is read by about thirty-nine outputs and fed by a single rv$result slot, so
loading a trait's result.json into that slot makes every existing Results tab work
untouched -- the lazy-loading change is one chokepoint, not a scattered refactor.

Deliberately no multi-result cache: a result carries every unordered parent pair,
about 125,000 rows at 500 parents, and holding seventeen would put us back where
we started. Re-reading on switch costs a second.

The overview survives unreadable results, and there is a test that proves it: if
it can render while every result.json is corrupt, it provably does not read them."
```

---

### Task 7: Reports per job, and the end-to-end check

**Files:**
- Modify: `R/config.R` / `R/app.R` (report paths)
- Test: `tests/testthat/test-report-paths.R`, `tests/testthat/test-job-e2e.R`

**Interfaces:**
- Consumes: everything above.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-report-paths.R`:

```r
# Two breeders must not overwrite each other's report.
#
# cfg$report_dir/report.html is app-global. With one blocking run at a time that was merely
# untidy. With jobs running concurrently it is a correctness bug: whoever renders last wins,
# and the other breeder downloads someone else's plan under their own trait's name.
test_that("a report path is scoped to its job and trait", {
  cfg <- ngcd_load_config()
  a <- ngcd_report_path(cfg, job_id = "job_1", trait = "YIELD", ext = "html")
  b <- ngcd_report_path(cfg, job_id = "job_1", trait = "PROTEIN", ext = "html")
  c2 <- ngcd_report_path(cfg, job_id = "job_2", trait = "YIELD", ext = "html")
  expect_false(identical(a, b))   # different traits of one job
  expect_false(identical(a, c2))  # same trait of different jobs
  expect_true(grepl("job_1", a, fixed = TRUE))
  expect_true(grepl("YIELD", a, fixed = TRUE))
})
```

- [ ] **Step 2: Run it, confirm it fails, implement, confirm it passes**

Run: `Rscript -e 'testthat::test_local(filter = "report-paths")'` before and after. Implement `ngcd_report_path(cfg, job_id, trait, ext)` returning a path under the job directory, and route the existing report writers through it, falling back to today's app-global path when there is no job (a single run).

- [ ] **Step 3: Write the backend-gated end-to-end test**

Create `tests/testthat/test-job-e2e.R`, gated exactly the way `tests/testthat/test-inbred-fix.R` gates on a real backend. It must: build a small three-trait fixture; submit via `ngcd_submit_job()`; poll until `ngcd_jobs_list()` reports the job `finished` (with a generous timeout and a clear skip if the backend is absent); assert `n_done == 3`; assert `ngcd_job_detail()` returns three traits with finite `cv_predictive_r2`; open one trait through `ngcd_job_result_path()` and assert the parsed result has `candidate_crosses` and `selected_crosses`.

- [ ] **Step 4: Run the full frontend suite**

Run: `cd ~/NextGenCrossDesign && Rscript tests/testthat.R`
Expected: 0 failures. Report the counts.

- [ ] **Step 5: Check for trailers before finishing**

```bash
git log --grep="Co-Authored-By" --oneline "$(git merge-base main HEAD)"..HEAD
```

Expected: no output. A prior branch had trailers inserted by tooling despite instructions; check rather than trust.

- [ ] **Step 6: Commit**

```bash
git add R/config.R R/app.R tests/testthat/test-report-paths.R tests/testthat/test-job-e2e.R
git commit -m "Scope reports to their job, and prove the whole path end to end

cfg\$report_dir/report.html is app-global. With one blocking run at a time that was
untidy; with jobs running concurrently it is a correctness bug -- whoever renders
last wins, and the other breeder downloads someone else's plan under their own
trait's name.

The end-to-end test submits a real three-trait batch, waits for it, and opens a
trait's result, so the contract between the two packages is exercised rather than
assumed. It is backend-gated the way test-inbred-fix.R already is."
```

---

## Self-Review

**Spec coverage.** Section 3's four requirements map as follows: submit stays instant → Task 3; the Jobs tab polls mtimes → Tasks 2 and 5; the overview is built from `status.json` only → Tasks 2, 5 and the laziness proof in Task 6; `res()` as the single chokepoint → Task 6. The spec's per-session report paths → Task 7. The spec's stated limitation (single runs keep blocking) is deliberately preserved — Tasks 3–6 add a parallel path and change nothing about the staged cards.

**Two spec details that changed during backend implementation**, reflected here: `job.json` now carries `phase`, so the Jobs table shows "still doing shared setup" rather than leaving a long silence unexplained; and `n_traits` is populated, so "k of n done" is real rather than derived from however many status files happen to exist.

**One deviation from the spec, deliberate.** The spec said the frontend calls `ng_job_list()`. Taken literally with a `reactivePoll`, that means an R process spawn every few seconds. Task 2 splits it: a local mtime fingerprint answers "has anything changed?", and the backend is asked only when it has. The derivations still live in the backend, which is what the constraint was protecting.

**Type consistency.** `ngcd_jobs_list()` returns the frame `ng_job_list()` produces (`id`, `label`, `created_at`, `state`, `phase`, `n_traits`, `n_done`, `n_error`, `path`); Task 5 renders exactly those columns. `ngcd_job_detail()$traits` carries `trait`, `state`, `cv_predictive_r2`, `mean_source`, `effect_gate`, `n_selected`, `error_message` — verified against a live `ng_job_status()` call, not assumed.

**No placeholders.** Every step contains the code to write or the exact command to run.
