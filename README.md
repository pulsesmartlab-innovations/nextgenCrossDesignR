# nextgenCrossDesign

**Genomic cross prediction and mate allocation for breeding programs.**

`nextgenCrossDesign` is an R package that predicts the genetic merit of every
candidate cross among a set of genotyped, phenotyped parents and then builds a
mating plan that trades genetic gain against diversity. It treats **cross
prediction** and **mate allocation** as two separate problems, so each can be
solved with the right method and audited independently.

<!-- badges -->
![version](https://img.shields.io/badge/version-0.21.0-blue)
![R](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)
![license](https://img.shields.io/badge/license-MIT-green)
![status](https://img.shields.io/badge/core%20path-validated%20(DH%2FRIL)-brightgreen)

> A guided Shiny front-end, **[nextgenCrossWorkbench](https://github.com/pulsesmartlab-innovations/NextGenCrossDesign)**,
> wraps this engine for point-and-click use. This repository is the backend engine
> and is fully usable on its own from R.

---

## What it does

The design is intentionally separated into four layers:

1. **Marker-effect estimation** — phenotype + genotype → additive (and optional
   dominance) marker effects.
2. **Cross-level metrics** — per-cross progeny **mean**, **within-family
   variance**, the **usefulness criterion** (mean + selection intensity × progeny
   SD), and prediction uncertainty.
3. **Candidate screening** — rank and shortlist crosses.
4. **Constrained mate allocation** — optimum-contribution-style selection that
   balances gain against group coancestry, with constraints on parent use,
   kinship, and group quotas.

The **within-family variance is recombination-aware**: for doubled-haploid / RIL
targets the C++ kernel computes the exact Haldane quadratic form

```text
Var(DH_ij) = a_ij' R a_ij ,   a_ijk = 0.5 (x_ik − x_jk) β_k ,
R_kl = 1 − 2 r_kl  (same chromosome, Haldane) , else 0
```

exactly, via a linear per-chromosome recursion rather than a dense marker×marker
matrix. The PMV extension propagates effect uncertainty with
`E(β_k²) = β_k² + Var(β_k)`.

## Highlights

- **Recombination-aware cross variance** (VPM / PMV) and the usefulness
  criterion, plus `var_simple`, expected-top-k, and kinship metrics.
- **Optimal mate allocation** — greedy, native evolutionary (GA), and MIP-based
  optimum-contribution allocation, with AlphaMate- and SimpleMating-style native
  components (no hard dependency on those tools).
- **Multi-trait selection** — direction-aware objectives with `auto`,
  `weighted`, `economic_index`, `desired_gain`, and soft/`strict` threshold
  methods.
- **Breeder decision controls** — a unified `mate_relatedness` dial (avoid
  inbreeding / favor complementarity), a **per-trait check-line veto** that
  screens candidates against a reference line before allocation, and a
  **portfolio & risk** decomposition (genetic level × within-family upside ×
  estimation risk) layered on the priority tiers — single-trait on the trait
  itself, multi-trait on the selection index.
- **Self-explaining plans** — every run returns `constraint_diagnostics` (what
  each allocation constraint did) and `priority_risk_diagnostics`, so a plan that
  came back smaller or different than requested explains itself.
- **Portfolio sizing** — sweep the number of crosses and recommend K by
  diminishing returns (elbow / kneedle), effective-population-size floor, or a
  coancestry budget.
- **Data QC + putative-duplicate removal** built into the preflight.
- **Polyploid model families** *(experimental)* — autotetraploid (tetrasomic,
  double reduction, dominance) and **allopolyploid disomic-subgenome** design
  with recombination-aware per-subgenome variance and VanRaden/Yang GRM.
- **Reproducible outputs** — priority-tiered cross lists, Excel workbooks, and
  JSON dashboard exports for downstream tools.

## Installation

```r
# from the source tarball attached to a release
# (github.com/pulsesmartlab-innovations/nextgenCrossDesignR/releases -> Assets):
install.packages("nextgenCrossDesign_0.21.0.tar.gz", repos = NULL, type = "source")

# or directly from GitHub:
# install.packages("remotes")
remotes::install_github("pulsesmartlab-innovations/nextgenCrossDesignR")
```

The package compiles a small C++ kernel via `Rcpp`, so a toolchain is needed at
install time (Rtools on Windows, Xcode CLT on macOS, `build-essential` on Linux).
Runtime imports: `Rcpp`, `lpSolve`, `Matrix`, `mvtnorm`, `ggplot2`, `openxlsx`,
`jsonlite` (plus base/recommended packages).

### Developer load (from a clone, without installing)

```r
source("R/load.R")
ng_load(use_cpp = TRUE)   # set NGCD_SKIP_CPP=1 or use_cpp = FALSE to skip the C++ build
```

`ng_load()` sources every `R/[0-9]*.R` file in numeric order and compiles
`src/ng_kernels.cpp`; the package falls back to pure R when the kernel is skipped.

## Quick start

The recommended one call runs the whole workflow — QC, duplicate handling,
marker-effect estimation, cross scoring, allocation, and priority ranking:

```r
library(nextgenCrossDesign)

result <- ng_run_cross_prediction(
  genotype  = geno,        # parents × markers (0/1/2 dosage)
  phenotype = pheno,       # parent IDs + trait(s)
  marker_map = map,        # marker, chromosome, position (bp or cM)
  n_crosses = 100
)

plan <- result$selected_crosses          # priority-tiered mating plan
ng_write_cross_priority_workbook(result) # reviewable Excel workbook
```

A copy-paste script is at
[`inst/examples/00_user_friendly_cross_prediction.R`](inst/examples/00_user_friendly_cross_prediction.R);
the lower-level, layer-by-layer workflow is
[`inst/examples/01_chronological_cross_prediction_workflow.R`](inst/examples/01_chronological_cross_prediction_workflow.R).

## Core capabilities

### Multi-trait selection

Build one selection objective from several traits with explicit directions
(increase some, decrease others):

```r
traits <- ng_multitrait_spec(
  trait     = c("yield", "lodging", "disease"),
  direction = c("maximize", "minimize", "minimize"),
  weight    = c(2, 1, 1)
)
plan <- ng_optimize_multitrait_mating_plan(scores, traits, n_crosses = 100)
```

Methods: `auto` (rank-normalized, equal weights unless supplied),
`economic_index` (Smith–Hazel; requires caller-supplied P and G),
`desired_gain` (Pesek–Baker; requires caller-supplied P and G), and `weighted`.
Thresholds are **soft by default** (missing a target is penalized, not
discarded); use `strict_thresholds = TRUE` for hard quality/market cutoffs. For
user-facing workflows prefer `ng_breeder_selection_objective()`.

### Breeder decision controls

Three controls on `ng_run_cross_prediction()` add real breeding decisions on top
of the objective, without touching how crosses are scored:

```r
result <- ng_run_cross_prediction(
  genotype = geno, phenotype = pheno, marker_map = map, n_crosses = 100,
  # 1. Relatedness intent as one dial (not two stacked lambdas):
  mate_relatedness = "avoid_inbreeding", mate_relatedness_weight = 20,
  # 2. Per-trait check-line veto — flag/remove crosses whose mid-parent for a
  #    trait is worse than a reference line, on GEBV or phenotype basis:
  trait_checks = data.frame(trait = "yield", check = "P001"),
  check_basis = "gebv", exclude_threshold_violators = FALSE
)
```

- **Unified mate-relatedness** — `mate_relatedness` names the intent
  (`avoid_inbreeding` / `favor_complementarity` / `off`) and sets the parent-pair
  penalty for you; setting it *and* the raw `lambda_mating` /
  `lambda_progeny_inbreeding` is a hard error, so two relatedness terms can never
  silently stack.
- **Per-trait check-line veto** — a candidate-eligibility filter applied *before*
  allocation (like the lethal guard), so the objective/index is never touched and
  it works with any multi-trait method. Flags by default; `exclude_threshold_violators = TRUE`
  drops violators. Not-evaluable crosses are counted, never silently passed.
- **Portfolio & risk profiling** — the runner attaches `cross_level`,
  `cross_upside`, a `risk_bin` from the mid-parent prediction-error variance, a
  `portfolio_profile` (`breakthrough` / `workhorse` / `long_shot` /
  `deprioritize`), and `portfolio_basis`. Single-trait resolves the axes on the
  trait itself (`cross_upside` = √VPM, the pure within-family genetic SD);
  multi-trait resolves them on the **selection index** (`cross_level` = w′m over
  mid-parent GEBVs, `cross_upside` = √(w′Sw) using the exact recombination-aware
  cross-trait covariance), and adds `risk_driver_trait` naming the component trait
  carrying most of the index uncertainty. Read `portfolio_basis` before
  interpreting the axes: for the rank-based index methods (`auto` / `weighted` /
  `threshold`) it is `linearized_rank_index`, meaning the quadrant is indicative
  rather than a decomposition of `multi_trait_score`. `cross_confidence` and
  `risk_bin` are within-run quantities and are not comparable across runs.

Every run also returns a `constraint_diagnostics` block (crosses requested vs
delivered, which caps bound the plan, lethal/marker/budget activity) so a plan
explains itself, plus `priority_risk_diagnostics` for the portfolio/risk layer.

### Portfolio size and diminishing returns

```r
curve <- ng_optimize_mating_plan_curve(scores, K_range = 3:20)
ng_plot_diminishing_returns(curve)   # needs ggplot2
```

Recommends K under `elbow_relative`, `elbow_kneedle`, `ne_target`, or
`coancestry_budget`, reporting total/marginal gain, group coancestry, effective
population size, and unique parent use. Opt-in decision support — it does not
change optimizer defaults.

### Data QC and putative-duplicate removal

`ng_preflight_input_tables(putative_duplicate_check = TRUE, putative_duplicate_action = "remove")`
detects near-identical genotype profiles, keeps one parent per duplicate cluster,
and returns cleaned tables plus audit metadata. Direct APIs
(`ng_detect_putative_duplicates()`, `ng_plot_putative_duplicates()`,
`ng_write_data_preflight_json()`) remain exported for reporting.

### Polyploid model families *(experimental)*

Two distinct, biologically-correct paths — pick by inheritance, not crop:

- **Autotetraploid (tetrasomic)** — `ng_poly4x_*` and `ng_polyploid_grm()` /
  `ng_polyploid_qc()` at ploidy 4/6, with allele-frequency GRM (VanRaden or
  Yang), digenic dominance GRM, additive+dominance effects, double reduction, and
  dominance-aware (heterosis + variance) cross scoring. For potato/cassava-like
  clonal crops. One-call: `ng_polyploid_design_crosses()`.
- **Allopolyploid disomic-subgenome** — `ng_polyploid_subgenome_*` for species
  whose subgenomes are inherited *diploidly* (each coded 0..2), e.g. true
  allopolyploids. Everything is subgenome-aware: per-subgenome QC, marker
  effects, VanRaden/Yang GRM, and **recombination-aware within-family variance**
  computed exactly per subgenome and summed — which is exact under disomic
  inheritance because homoeologues do not pair at meiosis. Falls back to a
  linkage-equilibrium approximation when no chromosome + cM map is supplied.

> Most polyploid crops (wheat, canola) are genotyped and analysed as diploid — use
> the standard path for those. The polyploid families above are for data you
> genuinely want modelled as autopolyploid or as separate diploid subgenomes.
> See [`VALIDATED_STATE.md`](VALIDATED_STATE.md) for what is and is not validated.

## Documentation

- **Vignette** — [`vignettes/nextgenCrossDesign.Rmd`](vignettes/nextgenCrossDesign.Rmd):
  an R-user tutorial starting from `ng_run_cross_prediction()`.
- **Backend user guide** — [`docs/BACKEND_USER_GUIDE.md`](docs/BACKEND_USER_GUIDE.md):
  recommended defaults by breeding scenario, input schemas, native components,
  validation runners, and report exports.
- **Validated state** — [`VALIDATED_STATE.md`](VALIDATED_STATE.md): the
  authoritative, short-form record of evidence-backed claims and the limits on
  PopVar/SimpleMating/AlphaMate/polyploid comparisons. Read this before describing
  a method as "better than" an external tool. [`BENCHMARK_NOTES.md`](BENCHMARK_NOTES.md)
  is the long-form history.
- **Statistical release gate** —
  [`docs/STATISTICAL_RELEASE_GATE.md`](docs/STATISTICAL_RELEASE_GATE.md): a
  reproducible installed-package audit of quantitative-genetic identities,
  compiled-kernel parity, LD graph pruning, relationship scales, formal
  multi-trait indices, probability metrics, and allocation constraints. This
  code gate is separate from the required historical forward-validation gate.
- **Capability registry** — `ng_backend_capability_registry()` enumerates the
  methods, breeding systems, and UI controls a front-end can discover.

## Testing and package build

Tests are plain `Rscript` files (each sources `tests/helper_load.R`); many use
`stopifnot()`, so a passing run exits silently with status 0:

```bash
Rscript tests/smoke_test.R
```

The package gate is the staged build wrapper (do **not** run `R CMD build` from
the repo root — non-package files confuse it):

```bash
Rscript tools/check_r_package.R   # staged R CMD build + R CMD check
```

Benchmark and validation runners live in `tools/` (AlphaSimR benchmarks,
head-to-head studies, parent-size grids, metric-merit studies). Common env
overrides: `NG_USE_CPP=1` (C++ kernel for 5K-marker runs),
`NG_ALPHASIMR_THREADS=1` (required for fair same-seed comparisons),
`NG_SHARED_SCORING=1` (shared score table across methods). See the runner headers
and `docs/BACKEND_USER_GUIDE.md` for the full menu.

The production-namespace statistical gate installs the package first and does
not source the test suite:

```bash
mkdir -p /tmp/ngcd_release_gate_lib
R CMD INSTALL --preclean --no-multiarch --library=/tmp/ngcd_release_gate_lib .
NGCD_RELEASE_LIB=/tmp/ngcd_release_gate_lib Rscript tools/run_statistical_release_gate.R
```

## Design note: DH/RIL target

The default target is **DH/RIL cross design**, not generic F1 progeny variance.
For inbred-line crossing the relevant segregation variance is the recombination
variance of the F1 haplotype contrast — the quadratic form shown above. Because
Haldane decay is exponential, each chromosome's quadratic form is computed
exactly with a linear recursion. The first correctness check for any metric is
not whether it wins a benchmark but whether predicted within-family variance has
the right **direction and scale** against realized simulated family variance.

## Versioning and changelog

Releases are tagged from `v0.4.0` through `v0.14.0`; later versions are on `main`
but not yet tagged, so install from a built tarball or from `main` for those. See
[GitHub Releases](https://github.com/pulsesmartlab-innovations/nextgenCrossDesignR/releases)
for full notes. Recent highlights:

- **0.21.0** — exact phased autopolyploid within-family variance (`phased_haplotypes` +
  `marker_map`, `variance_model = "phased_exact"`; validated against simulated meiosis and
  collapsing exactly to the dosage result when unlinked) and posterior-ON cross confidence
  (`confidence_method = "posterior_ci"` on the metric the run actually ranks on, plus
  `prob_top_tier`).
- **0.20.0** — polyploid quantitative-genetics audit: identified additive/dominance
  split (dominance orthogonalized against the additive design), correct
  frequency-centred GRM for `poly4x` coancestry, simulated-variance Monte-Carlo SE
  surfaced, and the autopolyploid `variance_model` stated
  (`uniform_phase_prior_expectation` — an expectation under a uniform prior over
  compatible phase, but it cannot separate crosses differing only in linkage phase).
- **0.19.0** — multi-trait portfolio & risk on the selection index: `cross_level` =
  w′m, `cross_upside` = √(w′Sw) from the exact cross-trait covariance, plus
  `portfolio_basis` and per-cross risk attribution.
- **0.14.0** — per-trait check-line veto (`trait_checks` / `check_basis` /
  `exclude_threshold_violators`): flag or remove candidate crosses whose
  mid-parent for a trait is worse than a reference line, before allocation.
- **0.13.0** — single-trait portfolio & risk decision layer (`cross_level`,
  `cross_upside`, `risk_bin`, `portfolio_profile`) plus `priority_risk_diagnostics`.
- **0.12.0** — `constraint_diagnostics` on every run (what each allocation
  constraint did).
- **0.11.0** — unified `mate_relatedness` control (hard error on stacking with the
  raw relatedness lambdas).
- **0.9.0** — recombination-aware within-family variance for the allopolyploid
  disomic-subgenome path (exact per-subgenome `a′R a`, summed) and VanRaden/Yang
  GRM for subgenomes.
- **0.7.0** — public-API parameter-vocabulary rename; core packages promoted to
  Imports.
- **0.6.0 / 0.5.0** — metric score-column and polyploid API renames.
- **0.4.0** — ploidy-general polyploid path (additive default, dominance
  optional) with allele-frequency GRM, dominance GRM, and ploidy-aware QC.

## Citation and license

Released under the **MIT License** (see [`LICENSE`](LICENSE)). If you use
`nextgenCrossDesign` in published work, please cite the package and the
usefulness-criterion / recombination-variance literature referenced in the
vignette (Meuwissen 1997; Lehermeier et al. 2017; VanRaden 2008).
