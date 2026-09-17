# Batch Job Registry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make a multi-trait batch submittable and abandonable — a detached process runs it, workers report their own progress to disk, and a reader can tell at any moment which traits finished, which failed, and whether the job is still alive.

**Architecture:** A new `R/55_job_registry.R` owns the on-disk job record, the per-trait status files, the crash derivation and the content-addressed shared-artefact store. `R/54_batch_runner.R` keeps the batch itself and gains hooks into that registry: workers write status beside the results they already write, and the parent touches a heartbeat while it waits. Nothing about how a trait is scored changes.

**Tech Stack:** R (>= 4.1.0), base R plus `jsonlite` and `mirai` (both already Imports), `tools::md5sum` for file hashing. No new dependencies.

**Spec:** `docs/superpowers/specs/2026-09-17-batch-job-registry-design.md`

## Global Constraints

- **Backend only.** The frontend round follows; do not edit `~/NextGenCrossDesign`.
- **OS-agnostic.** Must behave identically on Windows, macOS and Linux. No fork-based parallelism, no `ps`/`tasklist` shelling for liveness, no POSIX-only paths.
- **No new package dependencies.** `DESCRIPTION` Imports stay: `stats, utils, grDevices, graphics, parallel, mirai, jsonlite, lpSolve, Rcpp, Matrix, mvtnorm, ggplot2, openxlsx`.
- **Every JSON write is atomic** — write `<path>.part`, then `file.rename()`. A reader polling a live directory must never parse a half-written file. Use the existing `ng_write_json_atomic(x, path)`.
- **Every JSON payload carries a `schema` field**, matching the existing `ng_run_result.v1` / `ng_batch_manifest.v1` convention.
- **Workers must never `warning()` or `message()`.** A condition raised in a mirai daemon never reaches the parent. Workers return conditions as data; only the parent signals. See `R/52_advisories.R:20-27`.
- **Every exported function needs an `\alias{}`** in `man/nextgenCrossDesign-api.Rd`. `tests/every_export_is_documented.R` fails the build otherwise, and `R CMD check` is the only tier that catches it.
- **Tests are plain `stopifnot()` scripts** in `tests/`, auto-discovered by `tools/run_tests.R`. They start with the `helper_load.R` preamble (never a bare `.libPaths(".Rlib")` prepend — the project `.Rlib` here is Windows-built and shadows working packages).
- **Never add `Co-Authored-By` trailers** to commits in this project.
- **Three tiers before pushing:** `tests/testthat` fast gate, `tools/run_tests.R` (~36 min local), and `R CMD check` — the last being the only tier that reads `man/`.
- Final version: `DESCRIPTION` `Version: 0.37.0`.

## File Structure

| File | Responsibility |
|---|---|
| `R/55_job_registry.R` | **new** — job record, per-trait status, crash derivation, job listing, retention, shared-artefact keying. Everything about *jobs as durable objects*. |
| `R/54_batch_runner.R` | modify — the batch calls into the registry: workers write status, the parent writes job state and touches the heartbeat, and the shared artefact is fetched from the store rather than written per job. |
| `tools/run_cross_prediction_json.R` | modify — `workflow = "batch"` accepts a `job_dir` and writes the job record, so a detached process is self-describing from its first act. |
| `man/nextgenCrossDesign-api.Rd`, `NAMESPACE`, `DESCRIPTION`, `NEWS.md` | modify — exports, aliases, version. |
| `tests/*.R` | **new** — one test file per task. |

The split is by responsibility: `R/54` answers "how is a batch computed", `R/55` answers "what is a job and how do we know its state". They change for different reasons.

---

### Task 1: Per-trait status files

**Files:**
- Create: `R/55_job_registry.R`
- Modify: `R/54_batch_runner.R` (`ng_cp__batch_run_one`, currently at `:128`)
- Test: `tests/job_trait_status_is_written.R`

**Interfaces:**
- Consumes: `ng_write_json_atomic(x, path)` from `R/54_batch_runner.R:197`.
- Produces: `ng_job_trait_status_write(job_dir, trait_id, state, fields = list())` returning the written path. `job_dir` is the per-JOB directory; the file lands at `<job_dir>/<trait_id>/status.json`. `state` is one of `"running"`, `"done"`, `"error"`.

- [ ] **Step 1: Write the failing test**

Create `tests/job_trait_status_is_written.R`:

```r
# A trait's state must be readable from disk WHILE the batch runs, and after a trait fails.
#
# The state cannot be inferred from result.json: that file appears only at the END, so a
# running trait would look identical to a queued one, and a FAILED trait writes no
# result.json at all -- so the one state a breeder most needs to see would be invisible.
#
# It is also deliberately SMALL. The overview table needs R2 and the gate verdict for every
# trait; reading those out of seventeen result.json files, each carrying ~125,000 candidate
# crosses, would make the list view as costly as the thing it exists to avoid opening.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

job_dir <- file.path(tempdir(), paste0("ngcd_ts_", as.integer(runif(1, 1, 1e6))))
dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(job_dir, recursive = TRUE), add = TRUE)

# --- entry: the trait announces itself before doing any work -----------------------------
p <- ng_job_trait_status_write(job_dir, "YIELD", "running")
stopifnot(file.exists(p))
st <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(st$schema, "ng_job_trait_status.v1"))
stopifnot(identical(st$trait, "YIELD"))
stopifnot(identical(st$state, "running"))
stopifnot(is.character(st$started_at), nzchar(st$started_at))

# --- exit: the summary a reader needs, without opening result.json -----------------------
ng_job_trait_status_write(job_dir, "YIELD", "done", list(
  cv_predictive_r2 = 0.61, mean_source = "gebv", effect_gate = "on", n_selected = 20L))
st <- jsonlite::fromJSON(file.path(job_dir, "YIELD", "status.json"), simplifyVector = TRUE)
stopifnot(identical(st$state, "done"))
stopifnot(isTRUE(all.equal(st$cv_predictive_r2, 0.61)))
stopifnot(identical(st$mean_source, "gebv"), identical(st$effect_gate, "on"))
stopifnot(identical(as.integer(st$n_selected), 20L))
stopifnot(is.character(st$finished_at), nzchar(st$finished_at))
# started_at must SURVIVE the second write -- elapsed time is otherwise unknowable, and a
# reader cannot tell a trait that took ten seconds from one that took ten hours.
stopifnot(is.character(st$started_at), nzchar(st$started_at))

# --- failure carries its reason ----------------------------------------------------------
ng_job_trait_status_write(job_dir, "EMPTY", "running")
ng_job_trait_status_write(job_dir, "EMPTY", "error", list(error_message = "no finite values"))
st <- jsonlite::fromJSON(file.path(job_dir, "EMPTY", "status.json"), simplifyVector = TRUE)
stopifnot(identical(st$state, "error"))
stopifnot(grepl("no finite values", st$error_message, fixed = TRUE))

# --- a bad state is refused, not written --------------------------------------------------
# "queued" and "finished" belong to the JOB, not a trait. Accepting them here would let two
# vocabularies drift into one file and make the reader guess which it is looking at.
bad <- tryCatch({ ng_job_trait_status_write(job_dir, "X", "finished"); NULL },
                error = function(e) conditionMessage(e))
stopifnot(!is.null(bad), grepl("finished", bad, fixed = TRUE))

cat("PASS: job_trait_status_is_written\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_trait_status_is_written.R`
Expected: FAIL with `could not find function "ng_job_trait_status_write"`

- [ ] **Step 3: Write the minimal implementation**

Create `R/55_job_registry.R`:

```r
# A job is a batch that outlives the process that asked for it.
#
# 0.36.0 made a 17-trait batch computable. This file makes one DURABLE: submitted, left
# running, and readable afterwards by anything that can see the directory. Everything here
# is about a job as an object on disk -- its record, its per-trait progress, whether it is
# still alive -- and nothing here knows how a cross is scored.
#
# The whole design rests on one asymmetry: a worker can report what it did, but it cannot
# report that it died. So state is written optimistically by whoever is alive, and the
# absence of an update is what the reader interprets. See ng_job_status().

ng_job_trait_states <- c("running", "done", "error")

# Record one trait's state beside the results it writes.
#
# Called twice per trait: once on entry so a reader can distinguish a trait being worked on
# from one not yet started, and once on exit with the summary. The summary is SMALL on
# purpose -- an overview of seventeen traits must not mean parsing seventeen result.json
# files of ~125,000 rows each.
ng_job_trait_status_write <- function(job_dir, trait_id, state, fields = list()) {
  if (!state %in% ng_job_trait_states) {
    ng_stop("trait state must be one of ", paste(ng_job_trait_states, collapse = ", "),
            "; got '", state, "'. Job-level states (queued, finished, crashed) belong in ",
            "job.json, not beside a trait.")
  }
  dir <- file.path(job_dir, trait_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, "status.json")
  now <- format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  # Preserve started_at across the second write: without it nothing can say how long a trait
  # took, and elapsed time is the first thing anyone asks of a job that is still running.
  prior <- if (file.exists(path)) {
    tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
  } else NULL
  rec <- c(list(
    schema = "ng_job_trait_status.v1",
    trait = trait_id,
    state = state,
    started_at = if (!is.null(prior$started_at)) prior$started_at else now,
    finished_at = if (identical(state, "running")) NULL else now
  ), fields)
  ng_write_json_atomic(rec, path)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_trait_status_is_written.R`
Expected: `PASS: job_trait_status_is_written`

- [ ] **Step 5: Wire it into the worker**

In `R/54_batch_runner.R`, inside `ng_cp__batch_run_one()`, immediately after the `tryCatch({` opening (before `shared <- get0(...)`), add:

```r
      ng_job_trait_status_write(output_root, job$id, "running")
```

Then replace the success `list(status = "ok", ...)` return so the status is persisted before the value is returned. After `es <- r$effect_summary` and before the `list(` that begins `status = "ok"`, insert:

```r
      summary_fields <- list(
        cv_predictive_r2 = if (!is.null(es) && "cv_predictive_r2" %in% names(es))
          as.numeric(es$cv_predictive_r2)[[1L]] else NA_real_,
        mean_source = if (!is.null(es) && "mean_source" %in% names(es))
          as.character(es$mean_source)[[1L]] else NA_character_,
        # The gate verdict needs no new computation -- effect_gate already sits on the
        # effect_summary row this function reads for cv_predictive_r2 (R/39:1432) -- but it
        # is a new FIELD, and it is what tells a reader whether a plan rests on markers the
        # engine would otherwise have refused.
        effect_gate = if (!is.null(es) && "effect_gate" %in% names(es))
          as.character(es$effect_gate)[[1L]] else NA_character_,
        n_selected = if (is.data.frame(r$selected_crosses)) nrow(r$selected_crosses) else NA_integer_)
      ng_job_trait_status_write(output_root, job$id, "done", summary_fields)
```

And in the `error = function(e)` handler of the same `tryCatch`, replace the body with:

```r
    }, error = function(e) {
      ng_job_trait_status_write(output_root, job$id, "error",
                                list(error_message = conditionMessage(e)))
      list(status = "error", traits = job$traits, error_message = conditionMessage(e))
    }),
```

- [ ] **Step 6: Verify a real batch now leaves status files behind**

Run: `Rscript tests/batch_reports_what_its_workers_saw.R`
Expected: `PASS` (unchanged — the existing batch tests must still pass)

Then confirm the new files appear:

```bash
Rscript -e '
helper <- c("tests/helper_load.R","helper_load.R"); source(helper[file.exists(helper)][[1L]])
set.seed(5); n <- 18L; m <- 20L
g <- matrix(sample(c(0L,2L), n*m, replace=TRUE), nrow=n); colnames(g) <- paste0("M",seq_len(m))
ids <- paste0("P",seq_len(n))
tr <- function(s){set.seed(s); b<-rep(0,m); b[1:4]<-rnorm(4,sd=1.5); as.numeric(scale(g%*%b))+rnorm(n,sd=.3)}
cfg <- list(genotype=data.frame(NAME=ids,g,check.names=FALSE), genotype_id_col="NAME",
  phenotype=data.frame(NAME=ids,A=tr(1),B=tr(2)), phenotype_id_col="NAME",
  trait_direction=data.frame(Trait=c("A","B"),Selection_direction="increase"),
  direction_trait_col="Trait", direction_column_col="Trait",
  direction_direction_col="Selection_direction",
  marker_map=data.frame(SNP_code=colnames(g),Chromosome=rep(1:4,length.out=m),
                        Position_cM=rep(seq(0,60,length.out=m/4),times=4)),
  map_marker_col="SNP_code", map_chr_col="Chromosome", map_pos_cm_col="Position_cM",
  map_position_unit="cM", progeny="RIL", parent_type="ril", n_crosses=3L,
  write_outputs=FALSE, write_figures=FALSE, run_posterior_prediction=FALSE,
  min_cv_predictive_r2=-1, seed=1L)
root <- file.path(tempdir(),"t1check"); unlink(root, recursive=TRUE)
invisible(suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root=root, batch_workers=1L)))
print(list.files(root, recursive=TRUE, pattern="status.json"))'
```

Expected: `"A/status.json" "B/status.json"`

- [ ] **Step 7: Commit**

```bash
git add R/55_job_registry.R R/54_batch_runner.R tests/job_trait_status_is_written.R
git commit -m "Workers report their own state to disk, not only through the socket

A trait's state could not be read while a batch ran, and a FAILED trait left no
trace at all -- it writes no result.json, so the one state a breeder most needs
to see was the one that was invisible.

Each worker now writes <trait>/status.json on entry and on exit. It is small on
purpose: an overview of seventeen traits must not mean parsing seventeen
result.json files of ~125,000 candidate crosses each.

The gate verdict is a new field, though not new computation -- effect_gate
already sits on the effect_summary row the worker reads for cv_predictive_r2."
```

---

### Task 2: The job record and heartbeat

**Files:**
- Modify: `R/55_job_registry.R`
- Modify: `R/54_batch_runner.R` (`ng_run_cross_prediction_batch` at `:211`, `ng_cp__batch_dispatch` at `:314`)
- Test: `tests/job_record_and_heartbeat.R`

**Interfaces:**
- Consumes: `ng_write_json_atomic()`, `ng_job_trait_status_write()` from Task 1.
- Produces:
  - `ng_job_create(job_dir, config = NULL, label = NULL, n_traits = NA_integer_)` → path to `job.json`, state `"queued"`.
  - `ng_job_mark(job_dir, state, fields = list())` → path; `state` in `c("queued","running","finished","failed")`.
  - `ng_job_heartbeat(job_dir)` → invisible path; touches `<job_dir>/heartbeat`.

- [ ] **Step 1: Write the failing test**

Create `tests/job_record_and_heartbeat.R`:

```r
# A job must say what it is before it starts, and prove it is alive while it runs.
#
# "Alive" cannot be asked of the process directly: checking a pid means `ps` on Unix and
# `tasklist` on Windows, and this package has a hard requirement to behave identically on
# both. So the running process leaves a heartbeat, and a reader interprets its staleness.
# The pid is recorded for a human debugging afterwards, never for the decision.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

job_dir <- file.path(tempdir(), paste0("ngcd_jr_", as.integer(runif(1, 1, 1e6))))
on.exit(unlink(job_dir, recursive = TRUE), add = TRUE)

# --- created before anything runs ---------------------------------------------------------
p <- ng_job_create(job_dir, config = list(seed = 7L), label = "17-trait run", n_traits = 17L)
stopifnot(file.exists(p))
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$schema, "ng_job.v1"))
stopifnot(identical(rec$state, "queued"))
stopifnot(identical(rec$label, "17-trait run"))
stopifnot(identical(as.integer(rec$n_traits), 17L))
stopifnot(is.character(rec$id), nzchar(rec$id))
stopifnot(is.character(rec$created_at), nzchar(rec$created_at))

# The config travels with the job, so a run is replayable from its own directory rather than
# from whatever the submitting session happened to still hold.
stopifnot(file.exists(file.path(job_dir, "config.json")))

# --- state transitions --------------------------------------------------------------------
ng_job_mark(job_dir, "running", list(pid = Sys.getpid()))
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$state, "running"))
stopifnot(identical(as.integer(rec$pid), Sys.getpid()))
# created_at must survive -- it is how a listing orders jobs.
stopifnot(is.character(rec$created_at), nzchar(rec$created_at))

ng_job_mark(job_dir, "finished")
rec <- jsonlite::fromJSON(p, simplifyVector = TRUE)
stopifnot(identical(rec$state, "finished"))
stopifnot(is.character(rec$finished_at), nzchar(rec$finished_at))

# --- a bad state is refused ----------------------------------------------------------------
# "crashed" is DERIVED, never written: the process that would have written it is the one
# that died. Accepting it here would invite code to write a state it cannot honestly know.
bad <- tryCatch({ ng_job_mark(job_dir, "crashed"); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(bad), grepl("derived", bad, fixed = TRUE))

# --- heartbeat ------------------------------------------------------------------------------
hb <- file.path(job_dir, "heartbeat")
ng_job_heartbeat(job_dir)
stopifnot(file.exists(hb))
t1 <- file.mtime(hb)
Sys.sleep(1.1)
ng_job_heartbeat(job_dir)
stopifnot(file.mtime(hb) > t1)   # it must actually advance, not merely exist

cat("PASS: job_record_and_heartbeat\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_record_and_heartbeat.R`
Expected: FAIL with `could not find function "ng_job_create"`

- [ ] **Step 3: Write the minimal implementation**

Append to `R/55_job_registry.R`:

```r
# A job may be written as any of these. `crashed` is absent on purpose -- see ng_job_status().
ng_job_states <- c("queued", "running", "finished", "failed")

ng_job__now <- function() format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")

# Create the job record. Written BEFORE the work starts, so a job that dies during startup
# is still a job someone can find and read an error out of, rather than an empty directory.
ng_job_create <- function(job_dir, config = NULL, label = NULL, n_traits = NA_integer_) {
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(config)) ng_write_json_atomic(config, file.path(job_dir, "config.json"))
  rec <- list(
    schema = "ng_job.v1",
    id = basename(job_dir),
    label = if (is.null(label)) basename(job_dir) else label,
    created_at = ng_job__now(),
    state = "queued",
    n_traits = if (is.na(n_traits)) NA_integer_ else as.integer(n_traits),
    pid = NA_integer_
  )
  ng_write_json_atomic(rec, file.path(job_dir, "job.json"))
}

# Move the job to a new state, preserving what was already recorded.
ng_job_mark <- function(job_dir, state, fields = list()) {
  if (!state %in% ng_job_states) {
    ng_stop("job state must be one of ", paste(ng_job_states, collapse = ", "),
            "; got '", state, "'. 'crashed' is derived by ng_job_status() from a stale ",
            "heartbeat and is never written -- the process that would write it is the one ",
            "that died.")
  }
  path <- file.path(job_dir, "job.json")
  rec <- if (file.exists(path)) {
    tryCatch(as.list(jsonlite::fromJSON(path, simplifyVector = TRUE)),
             error = function(e) list())
  } else list()
  rec$schema <- "ng_job.v1"
  if (is.null(rec$id)) rec$id <- basename(job_dir)
  if (is.null(rec$created_at)) rec$created_at <- ng_job__now()
  rec$state <- state
  if (state %in% c("finished", "failed")) rec$finished_at <- ng_job__now()
  for (nm in names(fields)) rec[[nm]] <- fields[[nm]]
  ng_write_json_atomic(rec, path)
}

# Prove the process is still alive. Only the mtime matters; the file's content never does.
ng_job_heartbeat <- function(job_dir) {
  path <- file.path(job_dir, "heartbeat")
  if (!file.exists(path)) {
    dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)
    file.create(path)
  } else {
    Sys.setFileTime(path, Sys.time())
  }
  invisible(path)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_record_and_heartbeat.R`
Expected: `PASS: job_record_and_heartbeat`

- [ ] **Step 5: Beat the heart while the batch waits**

In `R/54_batch_runner.R`, in `ng_cp__batch_dispatch()`, replace the single collection line `out <- m[]` (currently `:390`) with:

```r
  # Collect by polling rather than with a bare m[], for one reason: m[] blocks until every
  # job resolves, and a job that is working hard for two hours would be indistinguishable
  # from a job whose process died an hour ago. Touching an mtime is all this costs -- no file
  # content is rewritten, and progress itself still comes from the workers' own status files.
  if (!is.null(job_dir)) {
    while (any(mirai::unresolved(m))) {
      ng_job_heartbeat(job_dir)
      Sys.sleep(2)
    }
  }
  out <- m[]
```

Add `job_dir = NULL` to the `ng_cp__batch_dispatch()` formals, and pass it from
`ng_run_cross_prediction_batch()` at the existing call site (`:276`) as `job_dir = job_dir`.

- [ ] **Step 6: Give the batch a job directory to record into**

In `ng_run_cross_prediction_batch()`, add `job_dir = NULL` to the formals after
`memory_budget_bytes`. Immediately after the `dir.create(output_root, ...)` line, add:

```r
  # When a job directory is given, this batch IS a job: it says so before the expensive work
  # starts, and says how it ended. Without one the function behaves exactly as in 0.36.0, so
  # a direct caller pays nothing for machinery it is not using.
  if (!is.null(job_dir)) {
    if (!file.exists(file.path(job_dir, "job.json"))) ng_job_create(job_dir, config = NULL)
    ng_job_mark(job_dir, "running", list(pid = Sys.getpid()))
    ng_job_heartbeat(job_dir)
  }
```

And immediately before the `structure(list(manifest = manifest, ...))` return, add:

```r
  if (!is.null(job_dir)) {
    any_failed <- any(vapply(results, function(r) identical(r$status, "error"), logical(1)))
    # "finished" means the BATCH completed, not that every trait succeeded -- a failed trait
    # is recorded per trait, and marking the whole job failed would hide sixteen good plans
    # behind one bad one.
    ng_job_mark(job_dir, "finished", list(n_ok = sum(!vapply(
      results, function(r) identical(r$status, "error"), logical(1))),
      any_failed = any_failed))
  }
```

- [ ] **Step 7: Verify the existing batch tests still pass**

Run: `Rscript tests/batch_matches_individual_runs.R && Rscript tests/batch_is_order_and_worker_independent.R`
Expected: both `PASS`

- [ ] **Step 8: Commit**

```bash
git add R/55_job_registry.R R/54_batch_runner.R tests/job_record_and_heartbeat.R
git commit -m "A job says what it is before it starts, and proves it is alive while it runs

The job record is written BEFORE the work, so a batch that dies during startup is
still a job someone can find rather than an empty directory.

Liveness is a heartbeat, not a pid check: pid liveness means ps on Unix and
tasklist on Windows, and this package must behave identically on both. The pid is
recorded for a human debugging afterwards, never for the decision.

Collection polls unresolved() instead of a bare m[], because m[] blocks until
every job resolves -- a batch working hard for two hours would otherwise be
indistinguishable from one whose process died an hour ago. It touches an mtime and
nothing else; progress still comes from the workers' own status files.

'crashed' is deliberately not a writable state. The process that would write it is
the one that died."
```

---

### Task 3: Deriving job status, including crashed and incomplete

**Files:**
- Modify: `R/55_job_registry.R`
- Test: `tests/job_status_derives_crashed_and_incomplete.R`

**Interfaces:**
- Consumes: `ng_job_create()`, `ng_job_mark()`, `ng_job_heartbeat()`, `ng_job_trait_status_write()`.
- Produces: `ng_job_status(job_dir, stale_after_sec = 120)` → list with `id`, `label`, `created_at`, `state` (one of `queued`, `running`, `finished`, `failed`, `crashed`), `n_traits`, `n_done`, `n_error`, `n_running`, and `traits`, a data frame with columns `trait`, `state`, `cv_predictive_r2`, `mean_source`, `effect_gate`, `n_selected`, `error_message`.

- [ ] **Step 1: Write the failing test**

Create `tests/job_status_derives_crashed_and_incomplete.R`:

```r
# Job status is DERIVED, because the two states that matter most cannot be written.
#
# A worker writes "running" on entry. If it dies it cannot write its own epitaph. Normally
# the batch parent reconciles that trait to "error" -- but if the PARENT died too, nobody
# reconciles anything and the trait sits at "running" forever.
#
# So: a trait's state is only meaningful RELATIVE TO ITS JOB'S state. Under a crashed job,
# a trait still claiming "running" is reported as "incomplete". This rule lives here, in one
# tested place, precisely so a frontend never reimplements it and drifts.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

mk <- function(tag) {
  d <- file.path(tempdir(), paste0("ngcd_js_", tag, "_", as.integer(runif(1, 1, 1e6))))
  ng_job_create(d, config = NULL, label = tag, n_traits = 3L)
  d
}

# --- a healthy running job ----------------------------------------------------------------
d <- mk("live")
ng_job_mark(d, "running", list(pid = Sys.getpid()))
ng_job_heartbeat(d)
ng_job_trait_status_write(d, "A", "done", list(cv_predictive_r2 = 0.6, n_selected = 10L))
ng_job_trait_status_write(d, "B", "running")
s <- ng_job_status(d)
stopifnot(identical(s$state, "running"))
stopifnot(identical(as.integer(s$n_done), 1L))
stopifnot(identical(as.integer(s$n_running), 1L))
stopifnot(is.data.frame(s$traits), nrow(s$traits) == 2L)
stopifnot(all(c("trait", "state", "cv_predictive_r2") %in% names(s$traits)))
unlink(d, recursive = TRUE)

# --- a job whose process died --------------------------------------------------------------
d <- mk("dead")
ng_job_mark(d, "running", list(pid = 999999L))
ng_job_heartbeat(d)
Sys.setFileTime(file.path(d, "heartbeat"), Sys.time() - 3600)  # an hour of silence
ng_job_trait_status_write(d, "A", "done", list(cv_predictive_r2 = 0.6, n_selected = 10L))
ng_job_trait_status_write(d, "B", "running")                   # died mid-trait
s <- ng_job_status(d, stale_after_sec = 120)
stopifnot(identical(s$state, "crashed"))
# THE RULE: a trait cannot still be "running" under a job that is not.
b <- s$traits[s$traits$trait == "B", , drop = FALSE]
stopifnot(identical(as.character(b$state), "incomplete"))
# The finished trait keeps its result. A crash must not erase what completed.
a <- s$traits[s$traits$trait == "A", , drop = FALSE]
stopifnot(identical(as.character(a$state), "done"))
stopifnot(isTRUE(all.equal(as.numeric(a$cv_predictive_r2), 0.6)))
unlink(d, recursive = TRUE)

# --- a finished job is never reinterpreted ---------------------------------------------------
# Its heartbeat is stale by definition, because it stopped when the work did.
d <- mk("done")
ng_job_mark(d, "running", list(pid = Sys.getpid()))
ng_job_heartbeat(d)
Sys.setFileTime(file.path(d, "heartbeat"), Sys.time() - 86400)
ng_job_mark(d, "finished")
ng_job_trait_status_write(d, "A", "done", list(n_selected = 4L))
s <- ng_job_status(d)
stopifnot(identical(s$state, "finished"))
unlink(d, recursive = TRUE)

# --- a job that never started ----------------------------------------------------------------
d <- mk("queued")
s <- ng_job_status(d)
stopifnot(identical(s$state, "queued"))
stopifnot(identical(as.integer(s$n_done), 0L))
unlink(d, recursive = TRUE)

cat("PASS: job_status_derives_crashed_and_incomplete\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_status_derives_crashed_and_incomplete.R`
Expected: FAIL with `could not find function "ng_job_status"`

- [ ] **Step 3: Write the minimal implementation**

Append to `R/55_job_registry.R`:

```r
# Read a job's state, deriving what cannot be written.
#
# Two derivations live here and nowhere else:
#
#   crashed    -- the record says "running" but the heartbeat has gone quiet. Nothing writes
#                 this, because the writer would have to be the process that died.
#   incomplete -- a trait claiming "running" under a job that is not running. A worker that
#                 dies cannot write its own epitaph, and if the parent died too, nobody
#                 reconciles it.
#
# Keeping both in one function is the point. A frontend that recomputed them would eventually
# disagree with the backend about whether a job is alive, and the disagreement would surface
# as a breeder acting on a plan that was never finished.
ng_job_status <- function(job_dir, stale_after_sec = 120) {
  path <- file.path(job_dir, "job.json")
  if (!file.exists(path)) {
    ng_stop("no job record at ", job_dir,
            " -- a directory without job.json is not a job; if a process was launched for ",
            "it and never wrote one, it failed during startup and its log will say why.")
  }
  rec <- as.list(jsonlite::fromJSON(path, simplifyVector = TRUE))
  state <- rec$state %||% "queued"

  if (identical(state, "running")) {
    hb <- file.path(job_dir, "heartbeat")
    age <- if (file.exists(hb)) {
      as.numeric(difftime(Sys.time(), file.mtime(hb), units = "secs"))
    } else Inf
    if (is.finite(stale_after_sec) && age > stale_after_sec) state <- "crashed"
  }

  dirs <- list.dirs(job_dir, recursive = FALSE, full.names = TRUE)
  files <- file.path(dirs, "status.json")
  files <- files[file.exists(files)]
  live <- identical(state, "running")
  traits <- if (!length(files)) {
    data.frame(trait = character(0), state = character(0), cv_predictive_r2 = numeric(0),
               mean_source = character(0), effect_gate = character(0),
               n_selected = integer(0), error_message = character(0),
               stringsAsFactors = FALSE)
  } else {
    do.call(rbind, lapply(files, function(f) {
      st <- tryCatch(as.list(jsonlite::fromJSON(f, simplifyVector = TRUE)),
                     error = function(e) list())
      ts <- st$state %||% NA_character_
      # The rule: a trait's state is only meaningful relative to its job's.
      if (identical(ts, "running") && !live) ts <- "incomplete"
      chr <- function(x) if (is.null(x) || !length(x)) NA_character_ else as.character(x)[[1L]]
      num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[[1L]]
      data.frame(trait = chr(st$trait), state = ts, cv_predictive_r2 = num(st$cv_predictive_r2),
                 mean_source = chr(st$mean_source), effect_gate = chr(st$effect_gate),
                 n_selected = as.integer(num(st$n_selected)),
                 error_message = chr(st$error_message), stringsAsFactors = FALSE)
    }))
  }
  traits <- traits[order(traits$trait), , drop = FALSE]
  rownames(traits) <- NULL

  count <- function(w) sum(traits$state == w, na.rm = TRUE)
  list(schema = "ng_job_status.v1", id = rec$id %||% basename(job_dir),
       label = rec$label %||% basename(job_dir), created_at = rec$created_at,
       state = state, pid = rec$pid,
       n_traits = if (is.null(rec$n_traits)) nrow(traits) else as.integer(rec$n_traits),
       n_done = count("done"), n_error = count("error"),
       n_running = count("running"), n_incomplete = count("incomplete"),
       traits = traits, path = job_dir)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_status_derives_crashed_and_incomplete.R`
Expected: `PASS: job_status_derives_crashed_and_incomplete`

- [ ] **Step 5: Commit**

```bash
git add R/55_job_registry.R tests/job_status_derives_crashed_and_incomplete.R
git commit -m "Derive the two job states that cannot be written

A worker writes 'running' on entry; if it dies it cannot write its own epitaph.
Normally the parent reconciles that trait to an error -- but if the parent died
too, nobody reconciles anything and the trait sits at 'running' forever.

So a trait's state is only meaningful relative to its job's state: under a crashed
job, a trait still claiming 'running' is reported as 'incomplete'. Both that rule
and the stale-heartbeat crash derivation live in ng_job_status() and nowhere else,
so a frontend cannot reimplement them and drift -- the drift would surface as a
breeder acting on a plan that was never finished.

A crash does not erase what completed: finished traits keep their results."
```

---

### Task 4: Listing jobs

**Files:**
- Modify: `R/55_job_registry.R`
- Test: `tests/job_list_ignores_what_is_not_a_job.R`

**Interfaces:**
- Consumes: `ng_job_status()` from Task 3.
- Produces: `ng_job_list(jobs_dir, stale_after_sec = 120, limit = 100L)` → data frame with `id`, `label`, `created_at`, `state`, `n_traits`, `n_done`, `n_error`, `path`, newest first.

- [ ] **Step 1: Write the failing test**

Create `tests/job_list_ignores_what_is_not_a_job.R`:

```r
# A listing must survive the directory containing things that are not jobs.
#
# The jobs directory is a shared location on disk. Editors leave dotfiles, users drop
# folders in, an interrupted submit leaves a half-made directory with no job.json. A listing
# that errored on any of those would take the whole Jobs view down with it -- and the view
# is the only way a breeder reaches results at all.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_jl_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

for (nm in c("job_a", "job_b", "job_c")) {
  d <- file.path(root, nm)
  ng_job_create(d, config = NULL, label = paste("run", nm), n_traits = 2L)
  ng_job_mark(d, "finished")
  ng_job_trait_status_write(d, "A", "done", list(n_selected = 3L))
  Sys.sleep(1.1)   # distinct created_at, so ordering is testable rather than tied
}
# Things that are not jobs.
dir.create(file.path(root, "not_a_job"), showWarnings = FALSE)          # no job.json
dir.create(file.path(root, ".DS_Store_dir"), showWarnings = FALSE)
writeLines("x", file.path(root, "stray.txt"))                            # a file, not a dir

lst <- ng_job_list(root)
stopifnot(is.data.frame(lst))
stopifnot(identical(nrow(lst), 3L))                       # exactly the three real jobs
stopifnot(all(c("id", "label", "state", "n_done", "path") %in% names(lst)))
stopifnot(!("not_a_job" %in% lst$id), !("stray.txt" %in% lst$id))
# Newest first -- a breeder looking for what they just submitted should not have to scroll.
stopifnot(identical(lst$id, c("job_c", "job_b", "job_a")))
stopifnot(all(lst$state == "finished"))
stopifnot(all(as.integer(lst$n_done) == 1L))

# An empty or absent directory is a normal state on a fresh install, not an error.
stopifnot(identical(nrow(ng_job_list(file.path(root, "nope"))), 0L))
empty <- file.path(root, "empty"); dir.create(empty)
stopifnot(identical(nrow(ng_job_list(empty)), 0L))

cat("PASS: job_list_ignores_what_is_not_a_job\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_list_ignores_what_is_not_a_job.R`
Expected: FAIL with `could not find function "ng_job_list"`

- [ ] **Step 3: Write the minimal implementation**

Append to `R/55_job_registry.R`:

```r
# List the jobs in a directory, newest first.
#
# Tolerant by construction: the jobs directory is shared, and a stray folder, a dotfile or a
# half-made directory from an interrupted submit must not take the listing down. The listing
# is the only route a breeder has to their results, so it degrades by omitting what it cannot
# read rather than by failing.
ng_job_list <- function(jobs_dir, stale_after_sec = 120, limit = 100L) {
  empty <- data.frame(id = character(0), label = character(0), created_at = character(0),
                      state = character(0), n_traits = integer(0), n_done = integer(0),
                      n_error = integer(0), path = character(0), stringsAsFactors = FALSE)
  if (!length(jobs_dir) || is.na(jobs_dir) || !dir.exists(jobs_dir)) return(empty)
  dirs <- list.dirs(jobs_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[file.exists(file.path(dirs, "job.json"))]
  if (!length(dirs)) return(empty)
  dirs <- dirs[order(file.mtime(dirs), decreasing = TRUE)]
  if (is.finite(limit) && length(dirs) > limit) dirs <- dirs[seq_len(limit)]
  rows <- lapply(dirs, function(d) {
    s <- tryCatch(ng_job_status(d, stale_after_sec = stale_after_sec), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    data.frame(id = s$id, label = s$label,
               created_at = if (is.null(s$created_at)) NA_character_ else s$created_at,
               state = s$state, n_traits = as.integer(s$n_traits),
               n_done = as.integer(s$n_done), n_error = as.integer(s$n_error),
               path = d, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  out <- out[order(out$created_at, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_list_ignores_what_is_not_a_job.R`
Expected: `PASS: job_list_ignores_what_is_not_a_job`

- [ ] **Step 5: Commit**

```bash
git add R/55_job_registry.R tests/job_list_ignores_what_is_not_a_job.R
git commit -m "List jobs newest first, tolerating what is not a job

The jobs directory is shared: editors leave dotfiles, users drop folders in, an
interrupted submit leaves a directory with no job.json. The listing is the only
route a breeder has to their results, so it degrades by omitting what it cannot
read rather than by failing and taking the whole view down."
```

---

### Task 5: The shared artefact is shared across jobs

**Files:**
- Modify: `R/55_job_registry.R`
- Modify: `R/54_batch_runner.R` (`ng_run_cross_prediction_batch`, the `saveRDS` at `:262`)
- Test: `tests/shared_artifact_is_computed_once_across_jobs.R`

**Interfaces:**
- Consumes: `ng_cp__batch_shared_keys` (`R/54_batch_runner.R:70`), `ng_write_json_atomic()`.
- Produces:
  - `ng_shared_artifact_key(config)` → a single character hash.
  - `ng_shared_artifact_dir(shared_dir, key)` → path, created on demand.
  - `ng_shared_artifact_reference(shared_dir, key, job_id)` → records `job_id` in `meta.json`.

- [ ] **Step 1: Write the failing test**

Create `tests/shared_artifact_is_computed_once_across_jobs.R`:

```r
# The same data must be processed once, not once per job.
#
# Quality control, genotype cleaning, duplicate detection, LD PRUNING, the training-set
# alignment, the pair table and the GRM depend only on the genotypes, the map, and the
# settings that govern them -- never on which trait is being scored. 0.36.0 computes them
# once per batch; that is still once too often. A breeder who submits a second batch on the
# same data pays for all of it again, and LD pruning over a real marker panel is not cheap.
#
# The dangerous half of a cache is the KEY. A key that omitted an input a job depends on
# would silently reuse an artefact built from different data -- wrong numbers, no error.
# That is the same class of defect as the RNG-kind bug found in 0.36.0, which produced a
# different crossing plan from identical data and signalled nothing. So the safety
# assertions below matter more than the reuse assertion.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_sa_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(31); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
gfile <- file.path(root, "geno.csv"); pfile <- file.path(root, "pheno.csv")
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2)), pfile, row.names = FALSE)

cfg <- list(
  genotype_file = gfile, genotype_id_col = "NAME",
  phenotype_file = pfile, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L, ld_pruning = TRUE, grm_method = "vanraden")

k1 <- ng_shared_artifact_key(cfg)
stopifnot(is.character(k1), length(k1) == 1L, nzchar(k1))

# --- STABILITY: the same inputs give the same key -----------------------------------------
stopifnot(identical(ng_shared_artifact_key(cfg), k1))
# Settings that are NOT shared must not move it -- otherwise nothing would ever be reused,
# since per-trait settings are exactly what varies between jobs.
c2 <- cfg; c2$n_crosses <- 9L; c2$min_cv_predictive_r2 <- 0.5
stopifnot(identical(ng_shared_artifact_key(c2), k1))

# --- SAFETY: anything the artefact depends on must move it ---------------------------------
for (chg in list(
      function(x) { x$grm_method <- "yang"; x },
      function(x) { x$ld_pruning <- FALSE; x },
      function(x) { x$duplicate_threshold <- 0.97; x })) {
  stopifnot(!identical(ng_shared_artifact_key(chg(cfg)), k1))
}
# One byte of one input file must move it. Keying on a filename, or a name and a timestamp,
# would let edited data silently reuse an artefact built from the old data.
g2 <- g; g2[1, 1] <- if (g[1, 1] == 0L) 2L else 0L
write.csv(data.frame(NAME = ids, g2, check.names = FALSE), gfile, row.names = FALSE)
stopifnot(!identical(ng_shared_artifact_key(cfg), k1))
write.csv(data.frame(NAME = ids, g, check.names = FALSE), gfile, row.names = FALSE)
stopifnot(identical(ng_shared_artifact_key(cfg), k1))   # and restoring it restores the key

# --- SAFETY: an in-memory input must be hashed by CONTENT, not by shape --------------------
# Summarising with str() or dim() would be cheaper and would be a latent disaster: two
# different genotype matrices of the same shape would hash identically and silently reuse
# each other's artefact.
cm <- cfg
cm$marker_map$Position_cM[[1L]] <- cm$marker_map$Position_cM[[1L]] + 7   # same shape, new data
stopifnot(!identical(ng_shared_artifact_key(cm), k1))

# --- REUSE: a second batch on the same data does not redo the work -------------------------
shared_dir <- file.path(root, "shared")
qc_calls <- 0L
orig <- ng_cp__stage_qc
assign("ng_cp__stage_qc", function(...) { qc_calls <<- qc_calls + 1L; orig(...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__stage_qc", orig, envir = .GlobalEnv), add = TRUE)

b1 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = file.path(root, "out1"), batch_workers = 1L, shared_dir = shared_dir))
stopifnot(identical(qc_calls, 1L))

b2 <- suppressWarnings(ng_run_cross_prediction_batch(
  cfg, output_root = file.path(root, "out2"), batch_workers = 1L, shared_dir = shared_dir))
stopifnot(identical(qc_calls, 1L))   # the second batch ran NO quality control

# The second batch still produced real results, not an empty shell.
for (j in b2$jobs) stopifnot(identical(j$status, "ok"))
stopifnot(file.exists(file.path(root, "out2", "A", "result.json")))

# --- the artefact records who is using it, so retention can be safe ------------------------
meta <- jsonlite::fromJSON(file.path(shared_dir, k1, "meta.json"), simplifyVector = TRUE)
stopifnot(identical(meta$schema, "ng_shared_artifact.v1"))
stopifnot(length(meta$referenced_by) >= 1L)

cat("PASS: shared_artifact_is_computed_once_across_jobs\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/shared_artifact_is_computed_once_across_jobs.R`
Expected: FAIL with `could not find function "ng_shared_artifact_key"`

- [ ] **Step 3: Write the minimal implementation**

Append to `R/55_job_registry.R`:

```r
# The content key of the shared artefact for a configuration.
#
# Derived from ng_cp__batch_shared_keys -- the list 0.36.0 already refuses to let a job
# override, precisely because changing one of those settings would invalidate the shared
# work. That list therefore IS the definition of what the artefact depends on, and the
# existing test asserting it names only real runner arguments protects this cache too.
#
# File inputs are hashed by CONTENT, via tools::md5sum (base R, portable, no dependency).
# Keying on a path, or a path and a timestamp, would let edited data silently reuse an
# artefact built from the old data -- wrong numbers with nothing to signal them.
#
# In-memory inputs are hashed by SERIALISING them, never by summarising them. str() or
# dim() would be cheaper and would be a latent disaster: two different genotype matrices of
# the same shape would hash identically and silently reuse each other's artefact. A cache
# key must be a function of the content or it is not a cache key.
#
# Serialising a large matrix costs a temp file, paid once per batch against recomputing
# quality control. Callers submitting through the frontend materialise CSVs first anyway, so
# the path that matters takes the cheap md5sum-of-file branch.
ng_shared_artifact__digest <- function(x) {
  f <- tempfile(fileext = ".rds"); on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, compress = FALSE)   # uncompressed: deterministic bytes, and faster
  unname(tools::md5sum(f))
}

ng_shared_artifact_key <- function(config) {
  parts <- lapply(ng_cp__batch_shared_keys, function(k) {
    v <- config[[k]]
    if (is.null(v)) return(paste0(k, "=<null>"))
    if (is.character(v) && length(v) == 1L && !is.na(v) && file.exists(v)) {
      return(paste0(k, "=file:", unname(tools::md5sum(v))))
    }
    paste0(k, "=obj:", ng_shared_artifact__digest(v))
  })
  ng_shared_artifact__digest(paste(unlist(parts), collapse = "\n"))
}

ng_shared_artifact_dir <- function(shared_dir, key) {
  d <- file.path(shared_dir, key)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Record that a job depends on this artefact, so retention can refuse to delete it.
ng_shared_artifact_reference <- function(shared_dir, key, job_id) {
  d <- ng_shared_artifact_dir(shared_dir, key)
  path <- file.path(d, "meta.json")
  meta <- if (file.exists(path)) {
    tryCatch(as.list(jsonlite::fromJSON(path, simplifyVector = TRUE)), error = function(e) list())
  } else list()
  meta$schema <- "ng_shared_artifact.v1"
  meta$key <- key
  if (is.null(meta$created_at)) meta$created_at <- ng_job__now()
  meta$referenced_by <- unique(c(as.character(meta$referenced_by %||% character(0)),
                                 as.character(job_id)))
  ng_write_json_atomic(meta, path)
  invisible(path)
}
```

- [ ] **Step 4: Use the store in the batch**

In `R/54_batch_runner.R`, add `shared_dir = NULL` to `ng_run_cross_prediction_batch()`'s formals after `job_dir`. Replace the block that currently reads

```r
  shared_dir <- file.path(output_root, "_shared")
  dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)
  shared_path <- file.path(shared_dir, "shared.rds")
  saveRDS(ctx, shared_path)   # on disk, so a killed batch resumes without redoing QC
```

with:

```r
  # On disk, so a killed batch resumes without redoing QC -- and CONTENT-ADDRESSED when a
  # store is given, so a second batch on the same data reuses it rather than recomputing
  # quality control, LD pruning and the GRM it has already paid for.
  if (is.null(shared_dir)) shared_dir <- file.path(output_root, "_shared")
  key <- ng_shared_artifact_key(config)
  art_dir <- ng_shared_artifact_dir(shared_dir, key)
  shared_path <- file.path(art_dir, "shared.rds")
  saveRDS(ctx, shared_path)
  if (!is.null(job_dir)) {
    writeLines(key, file.path(job_dir, "shared_ref"))
    ng_shared_artifact_reference(shared_dir, key, basename(job_dir))
  }
```

Then hoist the reuse decision above the QC work. Replace

```r
  ctx <- ng_cp__build_ctx(config)
  ctx <- ng_cp__stage_qc(ctx)
```

with:

```r
  # Reuse before recomputing. The key covers every setting the batch spends plus the content
  # of each input file, so a hit means the artefact was built from exactly this data.
  reuse_dir <- if (is.null(shared_dir)) NULL else
    file.path(shared_dir, ng_shared_artifact_key(config))
  reuse_path <- if (is.null(reuse_dir)) NULL else file.path(reuse_dir, "shared.rds")
  if (!is.null(reuse_path) && file.exists(reuse_path)) {
    ctx <- readRDS(reuse_path)
  } else {
    ctx <- ng_cp__build_ctx(config)
    ctx <- ng_cp__stage_qc(ctx)
  }
```

Guard the two lines that follow so a reused context is not re-derived:

```r
  if (identical(ctx$qc$status, "blocker")) {
    ng_stop("QC blocker -- resolve before running a batch: ",
            paste(utils::head(ctx$qc$issues$message, 3L), collapse = "; "))
  }
  if (is.null(ctx$predict_prologue)) {
    ctx <- ng_ctx_put(ctx, predict_prologue = ng_cp__predict_prologue(ctx))
  }
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `Rscript tests/shared_artifact_is_computed_once_across_jobs.R`
Expected: `PASS: shared_artifact_is_computed_once_across_jobs`

- [ ] **Step 6: Verify nothing about the numbers changed**

Run: `Rscript tests/batch_matches_individual_runs.R`
Expected: `PASS` — a reused artefact must give the same answers as a freshly computed one.

- [ ] **Step 7: Commit**

```bash
git add R/55_job_registry.R R/54_batch_runner.R \
        tests/shared_artifact_is_computed_once_across_jobs.R
git commit -m "Process the same data once, across jobs and not merely within one

Quality control, cleaning, duplicate detection, LD pruning, the training-set
alignment, the pair table and the GRM depend only on the genotypes, the map and
the settings that govern them. 0.36.0 computed them once per batch; a breeder
submitting a second batch on the same data paid for all of it again.

The artefact is now content-addressed, keyed off ng_cp__batch_shared_keys -- the
list 0.36.0 already refuses to let a job override, precisely because changing one
invalidates the shared work. That list IS the definition of what the artefact
depends on, so the existing test asserting it names only real runner arguments
protects this cache too.

File inputs are hashed by content via tools::md5sum. Keying on a path, or a path
and a timestamp, would let edited data silently reuse an artefact built from the
old data -- the same shape of defect as the RNG-kind bug in 0.36.0, which produced
a different crossing plan from identical data and signalled nothing. The test
asserts the safety half harder than the reuse half."
```

---

### Task 6: Resume an unfinished job

**Files:**
- Modify: `R/54_batch_runner.R`
- Test: `tests/batch_resume_runs_only_what_is_unfinished.R`

**Interfaces:**
- Consumes: `ng_job_status()` (Task 3), the shared artefact store (Task 5).
- Produces: `ng_run_cross_prediction_batch(..., resume = FALSE)` — when `TRUE`, jobs whose trait status is already `done` are skipped.

- [ ] **Step 1: Write the failing test**

Create `tests/batch_resume_runs_only_what_is_unfinished.R`:

```r
# Resuming must finish a job, not restart it.
#
# Failure isolation (0.36.0) means one bad trait does not cost the other sixteen. But
# "the manifest records that it failed" is only half an answer; the other half is being able
# to finish. A batch may have run for hours, so re-running the traits that already succeeded
# would be its own kind of loss.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_rs_", as.integer(runif(1, 1, 1e6))))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

set.seed(77); n <- 18L; m <- 20L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
# EMPTY has no usable data, so its job fails -- and only its job.
pheno <- data.frame(NAME = ids, A = tr(1), B = tr(2), EMPTY = rep(NA_real_, n))
cfg <- list(
  genotype = data.frame(NAME = ids, g, check.names = FALSE), genotype_id_col = "NAME",
  phenotype = pheno, phenotype_id_col = "NAME",
  trait_direction = data.frame(Trait = c("A", "B", "EMPTY"), Selection_direction = "increase"),
  direction_trait_col = "Trait", direction_column_col = "Trait",
  direction_direction_col = "Selection_direction",
  marker_map = data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                          Position_cM = rep(seq(0, 60, length.out = m / 4), times = 4)),
  map_marker_col = "SNP_code", map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
  map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
  write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
  min_cv_predictive_r2 = -1, seed = 1L)

out <- file.path(root, "out")
b1 <- suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root = out, batch_workers = 1L))
st <- vapply(b1$jobs, function(j) j$status, character(1))
names(st) <- vapply(b1$jobs, function(j) j$id, character(1))
stopifnot(identical(unname(st[c("A", "B")]), c("ok", "ok")))
stopifnot(identical(unname(st["EMPTY"]), "error"))

# Record when A's result was written, so we can prove it is not rewritten.
a_mtime <- file.mtime(file.path(out, "A", "result.json"))
Sys.sleep(1.1)

# --- resume: only the unfinished trait runs -------------------------------------------------
ran <- character(0)
orig <- ng_cp__batch_run_one
assign("ng_cp__batch_run_one", function(job, ...) { ran <<- c(ran, job$id); orig(job, ...) },
       envir = .GlobalEnv)
on.exit(assign("ng_cp__batch_run_one", orig, envir = .GlobalEnv), add = TRUE)

b2 <- suppressWarnings(ng_run_cross_prediction_batch(cfg, output_root = out,
                                                     batch_workers = 1L, resume = TRUE))
stopifnot(identical(sort(ran), "EMPTY"))                 # A and B were not re-run
stopifnot(identical(file.mtime(file.path(out, "A", "result.json")), a_mtime))

# The manifest still describes ALL traits, not only the resumed one -- a resumed job is the
# same job, and a manifest listing one trait would misrepresent what the breeder has.
ids2 <- sort(vapply(b2$jobs, function(j) j$id, character(1)))
stopifnot(identical(ids2, c("A", "B", "EMPTY")))
carried <- Filter(function(j) identical(j$id, "A"), b2$jobs)[[1L]]
stopifnot(identical(carried$status, "ok"))

cat("PASS: batch_resume_runs_only_what_is_unfinished\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/batch_resume_runs_only_what_is_unfinished.R`
Expected: FAIL — `unused argument (resume = TRUE)`

- [ ] **Step 3: Write the minimal implementation**

In `R/54_batch_runner.R`, add `resume = FALSE` to `ng_run_cross_prediction_batch()`'s formals after `shared_dir`. Immediately before the `results <- ng_cp__batch_dispatch(...)` call, insert:

```r
  # Resuming finishes a job rather than restarting it. A batch may have run for hours, so
  # re-running the traits that already succeeded is its own kind of loss.
  carried <- list()
  run_jobs <- jobs
  if (isTRUE(resume)) {
    done_of <- function(id) {
      f <- file.path(output_root, id, "status.json")
      if (!file.exists(f)) return(FALSE)
      st <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
      identical(st$state, "done")
    }
    keep <- !vapply(names(jobs), done_of, logical(1))
    carried <- lapply(names(jobs)[!keep], function(id) {
      st <- jsonlite::fromJSON(file.path(output_root, id, "status.json"), simplifyVector = TRUE)
      list(id = id, traits = jobs[[id]]$traits, status = "ok",
           cv_predictive_r2 = as.numeric(st$cv_predictive_r2 %||% NA_real_),
           mean_source = as.character(st$mean_source %||% NA_character_),
           n_selected = as.integer(st$n_selected %||% NA_integer_),
           result_json = file.path(output_root, id, "result.json"),
           warnings = character(0), elapsed_sec = NA_real_, resumed = TRUE)
    })
    names(carried) <- names(jobs)[!keep]
    run_jobs <- jobs[keep]
  }

  results <- if (!length(run_jobs)) list() else
    ng_cp__batch_dispatch(run_jobs, min(workers, length(run_jobs)), shared_path, output_root,
                          generated_at, package_version, use_cpp = ctx$use_cpp,
                          job_dir = job_dir)
  # The manifest describes the whole job. A resumed batch that listed only the traits it
  # re-ran would misrepresent what the breeder actually has on disk.
  results <- c(results, carried)
  results <- results[order(match(names(results), names(jobs)))]
```

and delete the previous `results <- ng_cp__batch_dispatch(...)` line it replaces.

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/batch_resume_runs_only_what_is_unfinished.R`
Expected: `PASS: batch_resume_runs_only_what_is_unfinished`

- [ ] **Step 5: Commit**

```bash
git add R/54_batch_runner.R tests/batch_resume_runs_only_what_is_unfinished.R
git commit -m "Resume finishes a job rather than restarting it

0.36.0 isolated failures so one bad trait does not cost the other sixteen, but
'the manifest records that it failed' is only half an answer. The other half is
being able to finish: a batch may have run for hours, so re-running the traits
that already succeeded is its own kind of loss.

The manifest still describes every trait, not only the resumed ones -- a resumed
batch is the same job, and a manifest listing one trait would misrepresent what
the breeder has on disk."
```

---

### Task 7: Retention that cannot destroy a running job

**Files:**
- Modify: `R/55_job_registry.R`
- Test: `tests/job_prune_protects_running_and_referenced.R`

**Interfaces:**
- Consumes: `ng_job_status()`, `ng_job_list()`, the artefact `meta.json` from Task 5.
- Produces: `ng_job_prune(jobs_dir, shared_dir = NULL, keep = 20L, stale_after_sec = 120)` → invisible character vector of removed paths.

- [ ] **Step 1: Write the failing test**

Create `tests/job_prune_protects_running_and_referenced.R`:

```r
# Retention deletes with unlink(recursive = TRUE, force = TRUE). Everything it will not
# delete has to be asserted, because the cost of being wrong is a breeder's work.
#
# Two absolute constraints:
#   1. never delete a RUNNING job -- it would take hours of finished traits with it;
#   2. never delete a shared artefact any surviving job still references -- the artefact
#      outlives the batch that built it, which is the entire point of sharing it.
#
# The frontend's ngcd_prune_runs() is the cautionary example: it enumerates
# list.dirs(cfg$runs_dir) and unlinks everything past keep_runs, which is exactly why the
# jobs and shared stores are siblings of runs_dir rather than living inside it.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

root <- file.path(tempdir(), paste0("ngcd_pr_", as.integer(runif(1, 1, 1e6))))
jobs_dir <- file.path(root, "jobs"); shared_dir <- file.path(root, "shared")
dir.create(jobs_dir, recursive = TRUE); dir.create(shared_dir, recursive = TRUE)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

mk <- function(id, state, key) {
  d <- file.path(jobs_dir, id)
  ng_job_create(d, config = NULL, label = id, n_traits = 1L)
  ng_job_mark(d, state, list(pid = Sys.getpid()))
  if (identical(state, "running")) ng_job_heartbeat(d)
  writeLines(key, file.path(d, "shared_ref"))
  ng_shared_artifact_reference(shared_dir, key, id)
  saveRDS(list(1), file.path(shared_dir, key, "shared.rds"))
  Sys.sleep(0.3)
  d
}
old1 <- mk("old1", "finished", "KEYOLD")
old2 <- mk("old2", "finished", "KEYSHARED")
live <- mk("live", "running",  "KEYLIVE")     # oldest-but-one, and must survive anyway
new1 <- mk("new1", "finished", "KEYSHARED")   # shares KEYSHARED with old2

removed <- ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 1L)

# --- a running job is never removed, however old ---------------------------------------------
stopifnot(dir.exists(live))
stopifnot(!(live %in% removed))
# --- the newest survives, as retention intends ------------------------------------------------
stopifnot(dir.exists(new1))
# --- the genuinely old, finished ones go -------------------------------------------------------
stopifnot(!dir.exists(old1))

# --- artefacts: referenced ones stay, orphans go ------------------------------------------------
# KEYSHARED is still referenced by new1, so it must survive old2's removal. Deleting it would
# make a surviving job unresumable, and its results unexplainable.
stopifnot(dir.exists(file.path(shared_dir, "KEYSHARED")))
stopifnot(dir.exists(file.path(shared_dir, "KEYLIVE")))
# KEYOLD's only referent is gone.
stopifnot(!dir.exists(file.path(shared_dir, "KEYOLD")))

# --- nothing to do is not an error ----------------------------------------------------------------
stopifnot(identical(length(ng_job_prune(jobs_dir, shared_dir = shared_dir, keep = 100L)), 0L))

cat("PASS: job_prune_protects_running_and_referenced\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_prune_protects_running_and_referenced.R`
Expected: FAIL with `could not find function "ng_job_prune"`

- [ ] **Step 3: Write the minimal implementation**

Append to `R/55_job_registry.R`:

```r
# Remove old jobs, and any shared artefact nothing references any more.
#
# This function deletes with unlink(recursive = TRUE, force = TRUE), so its refusals matter
# more than its removals:
#
#   * a RUNNING job is never removed, however old. Deleting one would take hours of finished
#     traits with it, and the process writing into it would carry on writing into nothing.
#   * a shared artefact referenced by ANY surviving job is never removed. The artefact
#     outlives the batch that built it -- that is the point of sharing it -- so deletion is
#     driven by references, not by age.
#
# This is also why jobs_dir and shared_dir are siblings of the frontend's runs_dir rather
# than living inside it: ngcd_prune_runs() there enumerates list.dirs(runs_dir) and unlinks
# everything past keep_runs, and would happily take the whole job store with it.
ng_job_prune <- function(jobs_dir, shared_dir = NULL, keep = 20L, stale_after_sec = 120) {
  removed <- character(0)
  keep <- suppressWarnings(as.integer(keep))
  if (is.na(keep) || keep < 0L || !dir.exists(jobs_dir)) return(invisible(removed))

  lst <- ng_job_list(jobs_dir, stale_after_sec = stale_after_sec, limit = Inf)
  if (nrow(lst)) {
    protected <- lst$state %in% c("running", "queued")
    candidates <- lst[!protected, , drop = FALSE]
    if (nrow(candidates) > keep) {
      doomed <- candidates$path[seq.int(keep + 1L, nrow(candidates))]
      unlink(doomed, recursive = TRUE, force = TRUE)
      removed <- c(removed, doomed)
    }
  }

  if (!is.null(shared_dir) && dir.exists(shared_dir)) {
    survivors <- ng_job_list(jobs_dir, stale_after_sec = stale_after_sec, limit = Inf)
    refs <- unique(unlist(lapply(survivors$path, function(d) {
      f <- file.path(d, "shared_ref")
      if (file.exists(f)) trimws(readLines(f, warn = FALSE)) else character(0)
    })))
    for (d in list.dirs(shared_dir, recursive = FALSE, full.names = TRUE)) {
      if (!(basename(d) %in% refs)) {
        unlink(d, recursive = TRUE, force = TRUE)
        removed <- c(removed, d)
      }
    }
  }
  invisible(removed)
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_prune_protects_running_and_referenced.R`
Expected: `PASS: job_prune_protects_running_and_referenced`

- [ ] **Step 5: Commit**

```bash
git add R/55_job_registry.R tests/job_prune_protects_running_and_referenced.R
git commit -m "Retention that refuses to destroy work

Pruning deletes with unlink(recursive = TRUE, force = TRUE), so its refusals
matter more than its removals: a running job is never removed however old, and a
shared artefact referenced by any surviving job is never removed -- the artefact
outlives the batch that built it, which is the point of sharing it.

This is also why the job and shared stores are siblings of the frontend's
runs_dir rather than inside it. ngcd_prune_runs() enumerates list.dirs(runs_dir)
and unlinks everything past keep_runs; a job store living there would be just
another directory to it."
```

---

### Task 8: Public surface, headless entry, and release

**Files:**
- Modify: `NAMESPACE`, `man/nextgenCrossDesign-api.Rd`, `DESCRIPTION`, `NEWS.md`
- Modify: `tools/run_cross_prediction_json.R`
- Modify: `R/25_backend_capability_registry.R`
- Test: `tests/job_workflow_runs_headlessly.R`

**Interfaces:**
- Consumes: everything from Tasks 1–7.
- Produces: exports `ng_job_create`, `ng_job_status`, `ng_job_list`, `ng_job_prune`; `workflow = "batch"` accepting `job_dir`, `shared_dir`, `resume`.

- [ ] **Step 1: Write the failing test**

Create `tests/job_workflow_runs_headlessly.R`:

```r
# The detached process must be self-describing from its FIRST act.
#
# The frontend will spawn one Rscript with wait = FALSE and then look for job.json. If the
# backend cannot load, or the config is bad, no record is ever written -- and "nothing
# happened" is the worst possible feedback for a run a breeder expects to take hours. So the
# headless entry writes the job record before it does anything expensive, and marks the job
# failed rather than vanishing.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])
root <- ng_test_find_root()

tmp <- file.path(tempdir(), paste0("ngcd_hw_", as.integer(runif(1, 1, 1e6))))
dir.create(tmp, recursive = TRUE); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

set.seed(9); n <- 16L; m <- 16L
g <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
colnames(g) <- paste0("M", seq_len(m)); ids <- paste0("P", seq_len(n))
tr <- function(s) { set.seed(s); b <- rep(0, m); b[1:4] <- rnorm(4, sd = 1.5)
  as.numeric(scale(g %*% b)) + rnorm(n, sd = 0.3) }
write.csv(data.frame(NAME = ids, g, check.names = FALSE), file.path(tmp, "geno.csv"), row.names = FALSE)
write.csv(data.frame(NAME = ids, A = tr(1), B = tr(2)), file.path(tmp, "pheno.csv"), row.names = FALSE)
write.csv(data.frame(SNP_code = colnames(g), Chromosome = rep(1:4, length.out = m),
                     Position_cM = rep(seq(0, 45, length.out = m / 4), times = 4)),
          file.path(tmp, "map.csv"), row.names = FALSE)
write.csv(data.frame(Trait = c("A", "B"), Selection_direction = "increase"),
          file.path(tmp, "dir.csv"), row.names = FALSE)

job_dir <- file.path(tmp, "jobs", "job_1")
cfg <- list(workflow = "batch", job_dir = job_dir,
            batch_output_root = job_dir, shared_dir = file.path(tmp, "shared"),
            batch_workers = 1L,
            genotype_file = file.path(tmp, "geno.csv"), genotype_id_col = "NAME",
            phenotype_file = file.path(tmp, "pheno.csv"), phenotype_id_col = "NAME",
            direction_file = file.path(tmp, "dir.csv"), direction_trait_col = "Trait",
            direction_column_col = "Trait", direction_direction_col = "Selection_direction",
            map_file = file.path(tmp, "map.csv"), map_marker_col = "SNP_code",
            map_chr_col = "Chromosome", map_pos_cm_col = "Position_cM",
            map_position_unit = "cM", progeny = "RIL", parent_type = "ril", n_crosses = 3L,
            run_posterior_prediction = FALSE, min_cv_predictive_r2 = -1, seed = 1L)
cfg_path <- file.path(tmp, "config.json"); res_path <- file.path(tmp, "result.json")
jsonlite::write_json(cfg, cfg_path, auto_unbox = TRUE)

status <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(file.path(root, "tools", "run_cross_prediction_json.R")),
                    shQuote(cfg_path), shQuote(res_path)),
                  stdout = NULL, stderr = NULL)
stopifnot(identical(as.integer(status), 0L))

# The job record exists and says the batch completed.
stopifnot(file.exists(file.path(job_dir, "job.json")))
s <- ng_job_status(job_dir)
stopifnot(identical(s$state, "finished"))
stopifnot(identical(as.integer(s$n_done), 2L))
stopifnot(identical(sort(s$traits$trait), c("A", "B")))

# --- a job that cannot run still leaves a record ---------------------------------------------
bad_dir <- file.path(tmp, "jobs", "job_bad")
bad <- cfg; bad$job_dir <- bad_dir; bad$batch_output_root <- bad_dir
bad$genotype_file <- file.path(tmp, "does_not_exist.csv")
bad_cfg <- file.path(tmp, "bad.json"); bad_res <- file.path(tmp, "bad_result.json")
jsonlite::write_json(bad, bad_cfg, auto_unbox = TRUE)
invisible(system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(file.path(root, "tools", "run_cross_prediction_json.R")),
                    shQuote(bad_cfg), shQuote(bad_res)),
                  stdout = NULL, stderr = NULL))
stopifnot(file.exists(file.path(bad_dir, "job.json")))
sb <- ng_job_status(bad_dir)
stopifnot(identical(sb$state, "failed"))

cat("PASS: job_workflow_runs_headlessly\n")
```

- [ ] **Step 2: Run it to make sure it fails**

Run: `Rscript tests/job_workflow_runs_headlessly.R`
Expected: FAIL — no `job.json` is written (the control key is rejected as unknown).

- [ ] **Step 3: Thread the job keys through the headless entry**

In `tools/run_cross_prediction_json.R`, extend the control-key split (the block that already removes `workflow`, `stage`, `batch_jobs`, `batch_output_root`, `batch_workers`):

```r
batch_job_dir <- if (!is.null(cfg$job_dir)) as.character(cfg$job_dir) else NULL
batch_shared_dir <- if (!is.null(cfg$shared_dir)) as.character(cfg$shared_dir) else NULL
batch_resume <- isTRUE(cfg$resume)
cfg$job_dir <- NULL; cfg$shared_dir <- NULL; cfg$resume <- NULL
```

Immediately after that block — before any validation that could stop the script — add:

```r
# Write the job record FIRST. The frontend spawns this process detached and then looks for
# job.json; if the backend cannot load or the config is bad, a job that never wrote a record
# is indistinguishable from one that was never launched, and "nothing happened" is the worst
# possible feedback for a run expected to take hours.
if (!is.null(batch_job_dir)) {
  ng_job_create(batch_job_dir, config = cfg, label = basename(batch_job_dir))
}
```

In the `batched` dispatch branch, pass the new arguments:

```r
      ng_run_cross_prediction_batch(cfg, jobs = batch_jobs, output_root = root_dir,
                                    batch_workers = batch_workers,
                                    job_dir = batch_job_dir,
                                    shared_dir = batch_shared_dir,
                                    resume = batch_resume,
                                    generated_at = generated_at, package_version = version)
```

And in the `inherits(res, "ng_run_json_error")` branch, before `write_result(...)`, add:

```r
  if (!is.null(batch_job_dir) && file.exists(file.path(batch_job_dir, "job.json"))) {
    ng_job_mark(batch_job_dir, "failed", list(error_message = res$message))
  }
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `Rscript tests/job_workflow_runs_headlessly.R`
Expected: `PASS: job_workflow_runs_headlessly`

- [ ] **Step 5: Export and document the public surface**

In `NAMESPACE`, after `export(ng_run_cross_prediction_batch)`:

```r
export(ng_job_create)
export(ng_job_list)
export(ng_job_prune)
export(ng_job_status)
```

In `man/nextgenCrossDesign-api.Rd`, after `\alias{ng_run_cross_prediction_batch}`:

```
\alias{ng_job_create}
\alias{ng_job_list}
\alias{ng_job_prune}
\alias{ng_job_status}
```

And in the `\details` section, after the `ng_run_cross_prediction_batch()` paragraph:

```
\code{ng_job_create()}, \code{ng_job_status()}, \code{ng_job_list()} and
\code{ng_job_prune()} make a batch durable: a job is a directory that outlives
the process which ran it, holding the configuration, a record of the job, and one
small status file per trait written as that trait starts and finishes. A reader
can therefore see which traits finished, which failed, and which are still
running, without opening the results themselves. \code{ng_job_status()} derives
the two states nothing can write: \code{crashed}, when the record says running
but the heartbeat has gone quiet, and \code{incomplete}, for a trait still
claiming to run under a job that is not -- a worker that dies cannot write its
own epitaph. \code{ng_job_prune()} never removes a running job, nor a shared
artefact any surviving job still references.
```

- [ ] **Step 6: Verify the documentation guard passes**

Run: `Rscript tests/every_export_is_documented.R`
Expected: `PASS: every_export_is_documented`

- [ ] **Step 7: Declare the retention control in the registry**

In `R/25_backend_capability_registry.R`, add before the `batch_workers` control:

```r
    num("keep_jobs", "Finished jobs to keep", "run", 20, min = 1, max = 500, step = 1,
        note = paste0("How many finished batch jobs are retained before the oldest are ",
                      "removed. A running job is never removed however old it is, and a ",
                      "shared data artefact is never removed while any surviving job still ",
                      "references it -- so lowering this frees finished results, never work ",
                      "in progress."),
        ),
```

- [ ] **Step 8: Bump the version and write NEWS**

In `DESCRIPTION`, set `Version: 0.37.0`.

Prepend to `NEWS.md`:

```markdown
# nextgenCrossDesign 0.37.0

## A batch can be submitted and left running

0.36.0 made a 17-trait batch computable. This makes one durable. A job is a directory that
outlives the process which ran it: the configuration it was given, a record of the job, and
one small status file per trait, written as that trait starts and again as it finishes. A
reader can see which traits are done, which failed, and which are still running -- without
opening the results themselves, which for seventeen traits would mean parsing seventeen
files of ~125,000 candidate crosses each.

`ng_job_status()` derives the two states nothing can write. A worker writes "running" on
entry; if it dies, it cannot write its own epitaph, and if the batch parent died too, nobody
reconciles it. So `crashed` is derived from a stale heartbeat, and a trait still claiming to
run under a job that is not is reported as `incomplete`. Both rules live in one function, so
a frontend cannot reimplement them and drift -- the drift would surface as a breeder acting
on a plan that was never finished.

Liveness is heartbeat staleness rather than a pid check, because pid liveness means `ps` on
Unix and `tasklist` on Windows and this package must behave identically on both. The pid is
recorded for a human debugging afterwards, never for the decision.

## The same data is processed once, across jobs

Quality control, cleaning, duplicate detection, LD pruning, the training-set alignment, the
pair table and the GRM depend only on the genotypes, the map and the settings that govern
them. 0.36.0 computed them once per batch; a breeder submitting a second batch on the same
data paid for all of it again.

The artefact is now content-addressed, keyed off the list of settings a batch spends -- the
same list that refuses per-job overrides, because changing one invalidates the shared work.
Input files are hashed by content: keying on a path, or a path and a timestamp, would let
edited data silently reuse an artefact built from the old data.

## Resume

`resume = TRUE` re-runs only the traits that did not finish, without recomputing quality
control. Failure isolation said one bad trait does not cost the other sixteen; this is the
other half of that answer.

## Retention that refuses to destroy work

`ng_job_prune()` never removes a running job, however old, and never removes a shared
artefact any surviving job still references.
```

- [ ] **Step 9: Run all three tiers**

```bash
Rscript -e 'setwd("tests"); source("testthat.R")' || true   # see note below
Rscript tools/run_tests.R
```

For the fast gate, the package must be installed first — build from a `/tmp` git archive
(this repo is on OneDrive and an in-place `R CMD build` hangs):

```bash
rm -rf /tmp/ngcd37 && mkdir -p /tmp/ngcd37
git archive --format=tar "$(git write-tree)" | tar -x -C /tmp/ngcd37
export PATH="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools:$PATH"
cd /tmp/ngcd37 && R CMD build --no-build-vignettes --no-manual .
rm -rf /tmp/ngcdlib37 && mkdir -p /tmp/ngcdlib37
R CMD INSTALL -l /tmp/ngcdlib37 nextgenCrossDesign_0.37.0.tar.gz
cd - && (cd tests && R_LIBS=/tmp/ngcdlib37 Rscript testthat.R)
cd /tmp/ngcd37 && R CMD check --no-tests --no-examples --no-vignettes --no-manual \
  --no-build-vignettes nextgenCrossDesign_0.37.0.tar.gz
```

Expected: fast gate 0 failures; harness `all green`; `R CMD check` **0 errors, 0 notes**
(the two vignette WARNINGs are artifacts of `--no-build-vignettes` and do not appear in CI).

If `R CMD check` reports `no visible binding for global variable`, add the offending names to
the `utils::globalVariables()` declaration in `R/39_cross_prediction_runner.R` if they come
from `ctx`, or add a new declaration in `R/55_job_registry.R` if they are injected some other
way — follow the pattern and comment already in those files.

- [ ] **Step 10: Commit**

```bash
git add NAMESPACE man/nextgenCrossDesign-api.Rd DESCRIPTION NEWS.md \
        R/25_backend_capability_registry.R tools/run_cross_prediction_json.R \
        tests/job_workflow_runs_headlessly.R
git commit -m "Expose the job registry, and make a detached run self-describing (0.37.0)

The headless entry writes the job record before it does anything expensive. The
frontend will spawn one Rscript detached and then look for job.json; if the
backend cannot load or the config is bad, a job that never wrote a record is
indistinguishable from one that was never launched, and 'nothing happened' is the
worst possible feedback for a run expected to take hours. A job that cannot run is
marked failed rather than vanishing."
```

---

## Self-Review

**Spec coverage.** Section 1 (layout and status contract) → Tasks 1–3. Section 1b (shared
artefact across jobs, inputs stored once) → Task 5. Section 2 (backend changes: worker
status, heartbeat, crash derivation, exported job functions, resume) → Tasks 1, 2, 3, 6, 8.
Section 4 (crash handling, retention) → Tasks 3 and 7. Section 5 (verification) → the test in
each task, plus the three-tier run in Task 8, Step 9.

Section 3 (frontend) is deliberately absent — this plan is backend only, as the spec states.

**One spec item intentionally deferred:** the spec's job layout names `job.log` for the
detached process's stdout and stderr. That file is created by the *frontend* when it spawns
the process (`system2(..., stdout = logfile, stderr = logfile, wait = FALSE)`), so it belongs
to the frontend round; the backend never writes it. Task 8's failing-job test covers the
backend half of that contract — a failed job still leaves a readable record.

**Type consistency.** `ng_job_trait_status_write(job_dir, trait_id, state, fields)` is called
with `output_root` as its `job_dir` in Task 1 Step 5; that is correct, because a batch's
`output_root` and its job directory are the same directory when the batch is run as a job
(Task 8 passes `batch_output_root = job_dir`). `ng_job_status()` returns `state` and a
`traits` data frame in Task 3; Task 4's `ng_job_list()` and Task 7's `ng_job_prune()` both
consume exactly those fields. `ng_shared_artifact_key(config)` returns a single character in
Task 5 and is used as a directory name in both Task 5 and Task 7.

**No placeholders.** Every step contains the code to be written or the exact command to run.
