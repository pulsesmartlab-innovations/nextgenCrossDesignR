# v0.3.0 C++ porting — merit checks and decisions

This document records what was ported to C++, what was skipped, and the
measurements that grounded each decision. Approach: profile-or-measure
first, port only when the data supports it.

## Microbenchmark fixture

n = 60 parents, m = 1500 markers, 1770 candidate pairs (DH/RIL all-pair set).
Reproduce via `tools/cpp_speedup_bench.R`.

## Wave 1 (committed: `258e489`)

| Kernel | Status | Speedup | Equivalence vs R | Merit grounding |
|---|---|---:|---:|---|
| Full off-diagonal posterior PMV | C++ | 4.3× | 8.9e-16 | a-priori: dominant R loop |
| Banded Kosambi / RIL | C++ | 5.7× | 1.1e-15 | a-priori: dominant R loop |
| BCM closed-form posterior sampler | C++ | 1.7× | 1.2e-16 | a-priori: per-draw loop |

## Wave 2 (committed: this commit)

| Kernel | Status | Speedup | Equivalence vs R | Merit grounding |
|---|---|---:|---:|---|
| `ng_local_swap` greedy OCS optimizer | C++ | **116×** | objective gap = 0 | measured: 109% of `greedy_local` total time |

## Wave 2 skips (with data)

### BLAS-backed full-posterior kernel — SKIP

**Hypothesis going in:** plain C++ nested loops do per-pair O(m²) compute; BLAS
dgemv would vectorize / cache-block and give another 5–10× over current C++.

**Measurement:** R's `R %*% a` at m=1500 takes **24.8 ms per call**. At 1770
pairs × 2 multiplies/pair = 3540 BLAS-dgemv calls = **87.8 s of pure BLAS
time**. Current plain-loop C++ does the same workload in **17 s** — 5.2×
faster than the BLAS ceiling.

**Why C++ beats BLAS here:** two structural advantages that BLAS cannot
exploit on this workload:
- **Column-skip**: the C++ kernel skips columns where `a_k == 0` AND `d_k == 0`.
  For inbred DH/RIL parents (geno ∈ {0, 2}), about half of markers have
  matching dosage in both parents → `d_k = 0` → skip the entire column.
- **Fused two-multiply loop**: the kernel does both `Ra` and `Sd` per pair
  in one outer loop over k, fetching each R-column / R_had_Sigma-column once.
  BLAS would call dgemv twice, fetching each column twice.

**Verdict:** plain C++ is already faster than BLAS at this fixture due to
problem-specific sparsity. The earlier 4.3× speedup (commit `258e489`) is
essentially the speed ceiling. Reproduce via
`tools/_merit_check_blas_and_mcmc.R`.

### MCMC Gibbs sampler — DEFER

**Hypothesis going in:** per-iteration cost dominated by per-iter LAPACK
calls (`dsyrk` for `XX'`, `dpotrf` for Cholesky); R-level loop overhead
across iterations would be the only C++ port win.

**Measurement:** R `ng_sample_ridge_posterior_mcmc()` at n=60, m=1500 runs
100 burnin + 100 draws in **4.05 s** (~20 ms/iter). Projected to a realistic
run (500 burnin + 500 draws): **~20 s**. At full v0.0.x scale (n=400, m=5000)
this projects to ~10× larger work, so ~3–4 min per MCMC fit.

**Why a C++ port gains less here:** the per-iter cost decomposes as
`XX'` (LAPACK dsyrk) + `chol(A)` (LAPACK dpotrf) + BCM beta draw (already
in C++ via `ng_bcm_posterior_sampler_cpp`) + two `rgamma`/`rinvgamma` calls
(trivial). The two LAPACK calls already run at near-peak BLAS speed. C++
port would only avoid the per-iter R dispatch overhead — predicted gain
2–3×, bringing the projected 3–4 min to 1–2 min at full scale.

**Verdict:** real but modest gain (2–3×) on an opt-in code path most users
don't hit. Defer until users actually report MCMC wall time as a problem.
Reproduce via `tools/_merit_check_blas_and_mcmc.R`.

## Cumulative status

Hot paths covered by C++:

- Haldane DH chromosome recursion (`ng_dh_recomb_pairs_cpp`, v0.0.x)
- LD pruning graph (`ng_ld_prune_graph_cpp`, v0.0.x)
- Banded Kosambi / RIL decay kernel (`ng_dh_recomb_pairs_banded_cpp`, v0.3.0)
- Full off-diagonal posterior PMV (`ng_dh_recomb_pairs_full_posterior_cpp`, v0.3.0)
- BCM closed-form posterior sampler (`ng_bcm_posterior_sampler_cpp`, v0.3.0)
- Greedy local-swap OCS optimizer (`ng_local_swap_cpp`, v0.3.0)

Hot paths still in R (intentionally, after merit check):

- MCMC Gibbs sampler (modest gain, opt-in path) — deferred.
- BLAS-backed full-posterior PMV (current C++ beats BLAS via column-skip
  and fused-loop) — skipped permanently.

Hot paths in R (not yet measured):

- `ng_posterior_cross_predict` outer loop (runs `ng_score_crosses` per
  draw). Most per-iter cost is inside `ng_score_crosses`, which already
  uses the C++ recombination kernels. Outer loop's R overhead is probably
  small — would need a profile to confirm.
- `ng_posterior_multitrait_cross_predict` outer loop (similar structure).
- `ng_posterior_genetic_covariance` `beta_posterior` mode loop (similar).

Per the discipline above, defer these until microbenchmarks show they
matter.

## Reproducing the bench and merit checks

```bash
cd nextgen_cross_design
Rscript tools/cpp_speedup_bench.R                  # speedup table (4 ports)
Rscript tools/_merit_check_local_swap.R            # ng_local_swap verdict
Rscript tools/_merit_check_blas_and_mcmc.R         # BLAS-full-posterior + MCMC verdicts
Rscript tests/cpp_kernels_equivalence.R            # machine-precision equivalence
```
