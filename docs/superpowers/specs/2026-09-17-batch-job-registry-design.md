# Batch job registry — design

**Date:** 2026-09-17
**Backend:** 0.37.0 (ships and merges first)
**Frontend:** the round after
**Depends on:** backend 0.36.0, which added `ng_run_cross_prediction_batch()`

## Context

Backend 0.36.0 made it possible to run many independent single-trait analyses as one
batch: seventeen traits, seventeen crossing plans, the trait-independent work done once.
It did not make that reachable from the workbench, and the way the workbench runs the
backend today makes it unreachable in principle.

`ngcd_run_backend()` (`R/run_backend.R:122`) is a blocking `system2()`, and so is the
per-stage `ngcd_run_stage()` (`R/run_backend.R:207`) — both call sites, not one. R is
single-threaded, so either freezes the **whole app process**, not the calling session. Every
other connected breeder's browser goes quiet for the duration. A batch that runs for hours
would turn seventeen short freezes into one very long one, which is worse than the manual
process it replaces.

The workbench also holds exactly one result. `rv$result` (`R/app.R:806`) is a single slot,
and a result is not small: `candidate_crosses` is every unordered parent pair, so 500
parents is ~125,000 rows, each carrying per-trait columns. Seventeen of those in a session
is not a storage strategy.

Intended outcome: a breeder submits a batch, closes the browser, comes back later, and
opens each trait's plan as they want it — with finished traits readable while the rest are
still running.

## Decisions

| Question | Decision |
|---|---|
| Execution model | Submit-and-detach. A job registry, not a blocking run. |
| Job ownership | A detached OS process. Survives session end and app restart; dies only with the host. |
| Landing view | Batch overview table first; drill into one trait. |
| Blast radius | Only batches become jobs. Single and index runs keep today's staged cards. |
| Progress reporting | Each worker writes its own small status file. |

The reference for the overall shape is PredictPro (`nexus.bigdata.ag.ndsu.edu/predictpro`):
a wizard ending in a Job Summary with **Submit**, a first-class **Jobs** tab, and results
pulled back later through a two-level selection (`Project` → `Model Result`), with large
tables paginated rather than materialised.

One structural difference shapes our design. In PredictPro a trait is a **column inside one
result** — the Best Models table lists traits as rows of a shared comparison. In our batch
each trait has **its own crossing plan**, with its own selected crosses. So "compare traits"
and "open a trait" are two different objects here, not two views of one table, and the
overview is a genuine summary rather than a slice.

## 1. Job store layout and status contract

```
<data_dir>/jobs/<job_id>/
├── job.json            ng_job.v1 — id, label, created_at, state, pid, n_traits
├── heartbeat           mtime, touched by the batch parent
├── job.log             stdout + stderr of the detached process
├── config.json         exactly what was submitted (replayable)
├── shared_ref          the content key of the shared artefact and inputs this job uses
├── manifest.json       ng_batch_manifest.v1, written at the end
├── YIELD/
│   ├── status.json     ng_job_trait_status.v1  — small, written by the worker
│   ├── result.json     the deliverable          — large, read only on drill-down
│   └── crossing_plan.xlsx
└── PROTEIN/ …
```

```
<data_dir>/shared/<content_key>/
├── inputs/             the materialised CSVs, stored ONCE and referenced by every job
├── shared.rds          QC + predict prologue — cleaned genotypes, LD-pruned markers,
│                       training-set alignment, pair table, GRM
└── meta.json           the key's inputs, and which jobs reference it
```

Job state: `queued → running → finished | failed | crashed`.
Trait state: `running → done | error`.

**`status.json` is separate from `result.json`, and small.** It carries trait, state,
timestamps, `cv_predictive_r2`, the cross-mean basis, the gate verdict, plan size and any
error message. State cannot be inferred from `result.json` instead, for three reasons:

- `result.json` appears only at the end, so a running trait would be indistinguishable
  from a queued one;
- a **failed** trait writes no `result.json` at all, so the state a breeder most needs to
  see would be the one that is invisible;
- the overview needs R² and the gate verdict for every trait, and reading those out of
  seventeen files of ~125,000 rows would make the list view as expensive as the thing it
  exists to avoid opening.

`crashed` is **derived, not written** — the process that would have written it is the one
that died.

Everything is schema-versioned to match `ng_run_result.v1` / `ng_batch_manifest.v1`, and
written through the atomic `.part`-then-rename path already in the package
(`ng_write_json_atomic`), because a frontend polling a live directory will otherwise
eventually parse a half-written file.

### Siblings of `runs_dir`, never inside it

`jobs/` and `shared/` sit beside `runs_dir` under `data_dir` — as `presets_dir` and
`report_dir` already do (`R/config.R:163`) — and **not** underneath it. This is not
tidiness. `ngcd_prune_runs()` (`R/run_backend.R:11`) enumerates `list.dirs(cfg$runs_dir)`
and calls `unlink(..., recursive = TRUE, force = TRUE)` on everything past the newest
`keep_runs` (default 20). A `jobs/` directory living there would be just another directory
to it: once a user accumulated twenty runs, prune would **delete the entire job tree and
the shared artefact store, including a running job**. `ngcd_run_index()`
(`R/runs.R:69`) would likewise list them as though they were runs.

So `cfg$jobs_dir` and `cfg$shared_dir` are new config entries alongside the existing three,
and job retention is a separate policy over `jobs_dir` with its own constraints.

## 1b. The shared artefact is shared ACROSS jobs, not just within one

Quality control, genotype cleaning, duplicate detection, **LD pruning**, the training-set
alignment, the pair table and the GRM depend only on the genotypes, the map, and the
settings that govern them. 0.36.0 computes them once per batch. That is still once too
often: a breeder who submits a second batch on the same data — a different metric, a
different acceptability bar, a re-run after fixing one trait — pays for all of it again,
and LD pruning over a real marker panel is not cheap.

So the artefact moves out of the job directory and becomes **content-addressed**:
`<data_dir>/shared/<content_key>/shared.rds`. A job records which key it used. If a job's
key already exists, it is reused and quality control never runs.

**The key is derived from `ng_cp__batch_shared_keys`.** That list already exists — it is
the set of settings a batch spends, which 0.36.0 refuses to let a job override, precisely
because changing one would invalidate the shared work. The same list therefore *is* the
definition of what the artefact depends on, and the existing test asserting the list names
only real runner arguments now protects the cache as well.

**Getting this key wrong is the dangerous failure.** A key that omits an input a job
actually depends on would silently reuse an artefact built from different data — wrong
numbers, no error, exactly the class of defect the RNG-kind bug in 0.36.0 turned out to be.
So the key covers every value in that list, plus a content digest of each input FILE via
`tools::md5sum` (no new dependency, portable).

Inputs supplied as in-memory data frames rather than file paths are **not cached** — they
fall back to computing fresh. Digesting a large in-memory genotype matrix would need either
a new dependency or a serialise-to-disk round trip, and the path that matters already
materialises CSVs to disk before submitting. Caching only what can be hashed cheaply and
honestly is better than hashing everything expensively or, worse, keying on something
weaker like a filename and a timestamp.

### The inputs themselves are stored once too

The materialised CSVs move into the same content-addressed directory rather than being
copied into each job. Ten batches on one genotype panel currently mean ten copies of that
panel on disk; referencing one copy means one. A job carries a `shared_ref`, not a private
duplicate of the data.

Retention gains a second absolute constraint: **never delete a shared artefact that any
surviving job references.** `meta.json` records the referencing job ids, so a shared
artefact — and the inputs beside it — outlives the job that created it and is collected
only when nothing points at it.

### What this does NOT eliminate, and why

On disk, the data is now stored once and processed once. **In memory it is still one copy
per concurrent worker**, and that is a consequence of the execution model rather than an
oversight.

mirai daemons are separate OS processes, which is exactly why they were chosen: identical
behaviour on Windows, macOS and Linux, and one worker's crash cannot take the batch with
it. Separate processes cannot share an R object's memory. Forked workers would share it
copy-on-write, but forking does not exist on Windows, and re-introducing a fork/PSOCK split
would give up the uniformity the batch was built for. Genuine shared memory (`bigmemory`,
`mmap`) would mean a new dependency and a platform-specific code path for a matrix the
scoring code reads through ordinary R subsetting.

What the design does instead is **bound** the duplication rather than pretend it away: the
shared context is read once per DAEMON (via `mirai::everywhere()`), never once per job, so a
daemon serving five traits reads it once; and the worker count is derived from the memory
actually available rather than from the core count, so the number of copies is chosen to
fit. That is the honest position — the cost is real, it is measured, and it is what sets
the concurrency.

## 2. Backend changes (0.37.0)

**Workers persist the summary they already produce.** `ng_cp__batch_run_one()` already
computes and returns trait, state, `cv_predictive_r2`, mean basis, plan size and error
message; today they travel back through the daemon socket and nowhere else. It also writes
them to `<trait>/status.json` — `running` on entry, `done`/`error` on exit.

One field is genuinely added: the **gate verdict** is not in the worker's current return
value. It is not new *computation* — `effect_gate` already sits on the `effect_summary` row
the worker reads for `cv_predictive_r2` (`R/39:1432`) — but it is a new field, and saying
"nothing new is invented" would have been an overstatement.

**A heartbeat, and the caveat it carries.** Workers report their own progress, so the parent
does not poll for that. But liveness needs someone alive to prove it, and
`mirai_map(...)[]` blocks until every job resolves. So the collection becomes a short loop
over `mirai::unresolved()` that touches the heartbeat and sleeps. This is polling, and it is
named rather than hidden — but it only touches an mtime, never rewrites content.

**Crash detection is heartbeat staleness, not pid liveness.** Pid checks differ per OS
(`ps` vs `tasklist`), and this package has a hard requirement to behave identically on
Windows, macOS and Linux. The pid is recorded for diagnostics; the decision rests on the
portable signal.

**The backend owns the job record.** Three new exported functions — `ng_job_create()`,
`ng_job_status()`, `ng_job_list()` — so the workbench calls them rather than hand-writing
JSON, per the standing rule that the frontend shows what the backend has and does not
manufacture what it lacks. `ng_job_status()` is where `crashed` and `incomplete` are
derived, so that rule lives in one tested place instead of being reimplemented in Shiny.

**Resume.** The shared artefact persists QC and the prologue, so `resume = TRUE`
re-runs only traits lacking a `done` status without recomputing quality control. This turns
failure isolation from "the manifest records that it failed" into "and here is how you
finish the job".

The manifest stays end-of-run. Progress is derived by the frontend from the status files.

## 3. Frontend

**Submit stays instant.** The Run tab gains a batch mode that materialises inputs
(`ngcd_materialize_inputs`), writes `config.json`, and spawns **one** detached `Rscript`
with `wait = FALSE`, stdout and stderr redirected to `job.log`. The detached process writes
`job.json` as its first act. Until it appears the UI shows *starting*; if it never appears
the job is *failed to start* and the log says why — a backend that cannot load produces no
job record at all, and "nothing happened" is the worst possible feedback.

**The Jobs tab is a poll over mtimes.** `reactivePoll` with a cheap `checkFunc` reading
directory mtimes and a `valueFunc` calling `ng_job_list()` — the same mtime-keyed-cache idea
`R/runs.R:81-85` already uses for the run browser. Newest first: label, created, state,
*k of n done*, elapsed.

**The overview is built from `status.json` only.** Seventeen small files; no `result.json`
parsed. That is what makes the list cheap enough to poll.

**`res()` is the leverage.** `res <- shiny::reactive(rv$result)` at `app.R:2652` is read by
**39 call sites**, and `rv$result` (`app.R:806`) is written at only three
(`app.R:2206`, `app.R:2613`, plus clears). Clicking a trait loads that trait's `result.json`
into `rv$result`, and every existing Results tab keeps working untouched. The lazy-loading
change really is one chokepoint rather than a scattered refactor — 39 readers, one writer to
change.

**Deliberately no multi-result cache.** Re-reading on switch costs a second or two and
bounds memory at one result. Caching seventeen puts us back where we started.

**One pre-existing bug comes into scope.** `cfg$report_dir/report.html` (`config.R:163`) is
app-global, so two concurrent users silently overwrite each other's report. Concurrent jobs
make that real rather than theoretical; reports move under the job/trait directory.

**Stated limitation:** this round does not fix single-run blocking. The staged cards keep
today's blocking `system2()`. Batch submit is non-blocking because it is detached, but a
single run still freezes every session.

## 4. Crash handling

| What dies | How it is seen |
|---|---|
| `Rscript` never starts | no `job.json` → *failed to start*, with `job.log` |
| Batch parent dies (OOM, container restart) | `job.json` says `running`, heartbeat stale → `crashed`; finished traits keep their results |
| A single worker dies | `mirai` returns a `miraiError`; the parent reconciles that trait to `error` |
| The app restarts | nothing happens to jobs; the tab repopulates from disk |
| Disk write interrupted | atomic rename means readers never see a partial file |

The hole is rows two and three together. A worker writes `running` on entry; if it dies it
**cannot write its own epitaph**. Normally the parent reconciles — but if the parent dies
too, traits sit at `running` forever.

**The rule: a trait's state is only meaningful relative to its job's state.** If the job is
`crashed`, any trait still claiming `running` is reported as `incomplete`. This is exactly
why the derivation belongs in `ng_job_status()` rather than in Shiny: one place, one rule,
tested once. A frontend reimplementing it would eventually disagree with the backend.

Pruning reuses the retention idea in `ngcd_prune_runs()` with two absolute constraints:
**never prune a running job**, and **never delete a shared artefact any surviving job
references**.

## 5. Verification

Backend, in the plain-`stopifnot` harness style:

- a worker writes `running` on entry and `done` on exit; a failed trait writes `error`, and
  the manifest agrees
- the heartbeat advances during a run; a stale heartbeat makes `ng_job_status()` report
  `crashed`
- a trait left at `running` under a crashed job reports `incomplete` — the reconciliation
  rule, asserted directly
- `resume = TRUE` re-runs only unfinished traits **and does not recompute quality control**,
  asserted against the shared artefact rather than inferred from timing
- schema contract for `ng_job.v1` and `ng_job_trait_status.v1`, guarded the way
  `test-contract.R` guards the capability registry
- `ng_job_list()` ordering, and that it ignores directories that are not jobs
- **a second batch on identical inputs reuses the shared artefact and does not re-run
  quality control or LD pruning** — asserted by counting the calls, the way
  `grm_is_computed_once_per_run.R` counts `ng_parent_kinship`
- **changing any one setting in `ng_cp__batch_shared_keys`, or one byte of an input file,
  produces a different key** — the safety half of the cache, and the half whose absence
  would silently return numbers computed from other data
- a shared artefact is not deleted while a job still references it
- **`ngcd_prune_runs()` never touches `jobs_dir` or `shared_dir`** — asserted by creating
  more than `keep_runs` runs alongside a job and confirming the job survives. This guards a
  destructive failure: the prune path calls `unlink(recursive = TRUE, force = TRUE)`.
- **the input CSVs exist once on disk no matter how many jobs use them**, asserted by
  counting files under the store rather than by inspecting any single job

Frontend, backend-gated as `test-inbred-fix.R` already is:

- **corrupt every `result.json`, then assert the overview still renders.** If the overview
  can be built while the large files are unreadable, it provably does not touch them. This
  turns "it is lazy" from a claim into a property, and fails loudly the day someone adds a
  convenient `fromJSON` to the list view.
- selecting a trait loads exactly one result into `rv$result`
- a job whose process never starts shows *failed to start* rather than hanging on *starting*
- the poll does not parse `result.json`

All three tiers before pushing: `tests/testthat` fast gate, `tools/run_tests.R`, and
`R CMD check` — the last being the only tier that reads `man/`, and the one that would catch
a new export without an alias.

## Deliberately not in scope

Converting single and index runs to jobs (the staged compute-once cards stay interactive).
Cross-trait comparison as a first-class artefact — which crosses recur across traits, where
plans agree or conflict — is a different question from "show me this trait's plan" and is
not currently computed. Multi-user identity and per-user job visibility: jobs are visible to
whoever can reach the deployment, as runs are today — and so, therefore, are shared
artefacts, which is the point of sharing them. Moving jobs off one machine (`crew`,
Slurm) — the job store is a directory contract, so that remains possible later without
changing the format.
