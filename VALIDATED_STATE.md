# Validated Software State

Last reviewed: 2026-08-25 for the installed-package statistical release gate.

## Statistical Release Gate (2026-08-25)

The production package was installed with its compiled kernels and evaluated
directly, without sourcing package test files. The reproducible audit and its
machine-readable results are in
[`docs/STATISTICAL_RELEASE_GATE.md`](docs/STATISTICAL_RELEASE_GATE.md) and
[`docs/STATISTICAL_RELEASE_GATE_RESULTS.csv`](docs/STATISTICAL_RELEASE_GATE_RESULTS.csv).

- **Mathematical/software quantitative-genetics gate: PASS (29/29).** This
  includes exact DH enumeration, the RIL-infinity identity, fast-versus-dense
  PMV, diagonal full-posterior reduction, 60 randomized graph-LD oracle cases,
  relationship/coancestry scales, optimizer hard constraints, Smith-Hazel and
  Pesek-Baker identities, polyploid model-domain checks, probability identities,
  and the installed one-call workflow. Graph LD pruning and fast PMV remain
  enabled and were explicitly validated.
- **R package gate: PASS.** The staged `R CMD build` and `R CMD check
  --no-manual --no-build-vignettes` completed with `Status: OK`.
- **Historical forward-validation gate: NOT RUN.** The code gate does not
  establish crop-specific prediction accuracy or realized genetic gain.
- **Unrestricted worldwide production release: HOLD** until leakage-free
  historical cross-by-progeny forward validation and an independent
  quantitative-genetics review pass. This hold is an evidence boundary, not a
  defect in graph LD pruning or fast PMV.

## Phased Autopolyploid Variance & Posterior Confidence Note (v0.21.0, 2026-08-18)

Closes two items that were carried as **Disallowed claims** in the 2026-08-16 note: the
autopolyploid within-family variance could not resolve linkage phase, and `cross_confidence`
had no posterior path (so it was a within-run ranking only, and P(top-N) did not exist).

Allowed claims:

- **Exact phased autopolyploid within-family additive variance.** With phased parental
  homologues and a cM map, for P homologues under random bivalent pairing:

      Var(gamete) = [ P * sum_i (b h_i)' R (b h_i) - (b d)' R (b d) ] / (4 (P - 1))

  with `d = sum_i h_i` the dosage and `R_kl = (1 - 2 r_kl)` the decay matrix the package
  already builds; the cross variance is the sum over the two parents. Derived from the bivalent
  model (within a bivalent the gamete keeps the same homologue across two loci with probability
  `1 - r`, and the pairing average puts weight `1/(P-1)` on the sum over all homologue pairs).
  Validated three ways, all in `tests/polyploid_phased_variance.R`:
  * against **directly simulated meiosis** (random bivalent pairing + crossovers along the
    chromosome): 0.43%-2.10% relative error at n = 15000-30000 draws, i.e. within Monte-Carlo
    error (~1%);
  * the **unlinked limit collapses EXACTLY to the dosage-only moment table** (max abs diff
    3.6e-15 through the function, 5.8e-15 through the scorer) -- so phase can change the answer
    only through LINKAGE, never on its own;
  * the **single-locus case equals the hypergeometric** `d(P-d)/(4(P-1))` exactly, for ploidy
    2/4/6 and every dosage.
  Exposed as `phased_haplotypes` + `marker_map` on `ng_polyploid_score_crosses()`, which then
  stamps `variance_model = "phased_exact"` instead of `"unlinked_phase_marginalized"`.
  Rownames follow the diploid convention generalized to ploidy: `<parent>_Hap1`..`<parent>_HapP`.
- **Phase changes the DECISION, not just the number.** Spearman between dosage-based and
  phased usefulness is 0.926 on a test panel -- the cross ORDER moves. This matters more than
  the posterior work below, because the autopolyploid variance feeds `cross_usefulness`, which
  feeds the index, which feeds allocation. Posterior confidence does not: the allocator
  optimizes `multi_trait_score` (or `marker_adjusted_gain`), and posterior output never reaches
  it in `ng_run_cross_prediction` -- `ng_optimize_robust_mating_plan` exists but is called only
  from tests and the vignette.
- **Cost is O(parents), not O(crosses).** The per-parent term does not depend on the mate, so it
  is computed once per parent (P + 1 quadratic forms) and each cross is a sum of two precomputed
  numbers.
- **Posterior-ON confidence summarizes the SELECTED metric (design-doc item F2).**
  `ng_posterior_cross_predict()` takes a `value_fun` and reports `ranked_value_post_sd`, so the
  runner summarizes the value it actually ranks on rather than the hardcoded
  `usefulness_pmv_gebv`. With posterior on, `cross_confidence` derives from the absolute
  posterior SD of that value on its native scale (never a CV) and
  `confidence_method = "posterior_ci"` -- with **no `_partial` suffix**, because unlike the
  mid-parent PEV fallback it covers the merit's variance term as well as its mean.
  `prob_top_tier` is populated from P(cross in top-N) at N = the plan size.
- **The F2 fix is measurable, not cosmetic.** Confidence computed under `trait_value_metric =
  "usefulness"` versus `"mean"` correlates at **Spearman -0.21** on the same data -- nearly
  opposite orderings. Under the previous behaviour both runs would have received the identical
  usefulness-based interval.

Disallowed claims:

- **The phased path covers ADDITIVE variance only.** `ng_polyploid_score_crosses_dominance()`
  has no phased path; its dominance (and additive) variance still comes from the dosage moment
  table. Do not describe dominance-aware polyploid scoring as phase-exact.
- **Double reduction is not modelled on the phased path** (random chromosome segregation only).
  Supplying `phased_haplotypes` together with `double_reduction > 0` is a hard error rather than
  a silent choice between them; if DR matters more than phase for a crop, use the dosage path.
- **Phase must be supplied, and is not inferred.** Without `phased_haplotypes` the autopolyploid
  variance remains `unlinked_phase_marginalized` -- unbiased over unknown phase but unable to
  separate crosses differing only in linkage phase. The package does not phase dosage data.
- **`prob_top_tier` is not the confidence axis.** It is a joint merit x uncertainty quantity
  ("is this cross genuinely top-N?"), reported as its own continuous column and never binned
  into `risk_bin` or relabelled as confidence.
- **No posterior path for MULTI-TRAIT runs.** Measured, not assumed. Two blockers:
  * *Definitional.* The per-draw index the design doc specifies is linear, but `auto` promotes
    to `weighted`, a RANK index, for which no linear per-draw value exists. Rank-normalizing
    within each draw gives every draw the same marginal distribution by construction, so the
    across-draw spread reflects rank REORDERING rather than magnitude uncertainty -- not a
    confidence signal.
  * *The cheap approximation does not work.* Holding the within-family covariance at its point
    estimate and varying only the mean per draw retains 98.5-99.3% of the posterior SD's
    MAGNITUDE but only 0.37-0.98 Spearman of its per-cross ORDERING (risk-bin agreement
    44%/53%/89% across three settings; three bins give 33% by chance). The mean term's
    uncertainty is largely shared across crosses, so the between-cross signal lives
    disproportionately in the variance term. Since `cross_confidence` is a within-run
    normalization and `risk_bin` is within-run tertiles, only the ordering matters -- so the
    approximation fails exactly where it would be used, and worst in small, low-h2 programs.
  * The full per-draw path is affordable at small scale and not at production scale: measured
    0.14 s per pass at 50 parents / 1k markers / T=2 (70 s for 500 draws), 1.3 s at 100 / 2k
    (11 min), and 40 s at 200 parents / 5k markers / T=3 (5.6 h).
  Multi-trait runs therefore keep the mid-parent PEV path and its existing disclosure.
- **No selection-superiority claim** for either feature. The phased variance is an accuracy and
  identifiability improvement; the posterior work is a reporting surface.

## Multi-Trait Portfolio & Polyploid QG Audit Note (v0.20.0, 2026-08-16)

Two efforts. (a) The cross-priority risk/portfolio layer, single-trait since v0.13.0, now
resolves its two axes on the SELECTION INDEX so multi-trait runs get the same six columns
(previously they got none). (b) A quantitative-genetics audit of the polyploid path, which had
not received the adversarial treatment the diploid path had; it validated most of the genetics
by direct simulation and found three defects.

Allowed claims:

- **Multi-trait portfolio axes.** `cross_level = w'm` over per-trait mid-parent GEBVs and
  `cross_upside = sqrt(w'Sw)`, with `S` the EXACT recombination-aware within-family cross-trait
  covariance (`ng_cross_trait_within_family_cov`, per-cross `a_t' R a_s`). The diagonal shortcut
  `sqrt(sum w^2 Var)` is not used and is not equivalent: on an antagonistic trait pair it
  overstated index SD by 1.85x at the median and up to 7x. `S` is rescaled so its diagonal
  reproduces the run's reported `vpm`, which makes `T == 1` collapse to exactly `sqrt(vpm)` and
  honours het-parent corrected variances. The index basis is resolved once on the candidate pool,
  so the selected plan and the candidate pool sit on the same axes.
- **Polyploid gamete transmission.** The hypergeometric gamete pmf and the classical double
  reduction model (with probability alpha the gamete is two IBD copies of one randomly chosen
  parental allele) reproduce directly simulated meiosis: max |empirical - analytic| = 0.0039
  against a Monte-Carlo SE of ~0.005, over ploidy 4 and 6, dr in {0, 0.25}, all parental dosages.
- **Double reduction preserves the mid-parent.** `E[gamete] = d/2` to 1.8e-15 for ploidy 2/4/6/8
  and every dr in [0, 1]. Double reduction redistributes gametes toward homozygosity but cannot
  shift expected dose, so predicted cross MEANS stay correct at any dr.
- **Polyploid cross moments.** The progeny-moment table reproduces simulated progeny moments, and
  `E[X]` is exactly the mid-parent dosage. End to end, the analytic `cross_mean` / `cross_var`
  reproduce the genotypic values of SIMULATED progeny to 0.001% (mean) and 1.1% (sd), the latter
  being the Monte-Carlo error at n = 3000.
- **Polyploid GRM.** `ng_polyploid_grm()` reduces EXACTLY to textbook VanRaden at ploidy 2
  (max abs diff 0) and keeps mean(diag) ~ 1 at ploidy 2/4/6 with unrelated off-diagonals ~ 0
  (-0.005 measured on 200 unrelated autotetraploids).
- **Identified additive/dominance split.** The dominance design is now the residual of
  `H = d(ploidy - d)` regressed on the additive design, using the OBSERVED per-marker coefficient
  (statistical parameterization; Alvarez-Castro & Carlborg 2007, Vitezica et al. 2013). Measured
  `cor(W, D)` falls from 0.983 to 0.000 in the fitted sample. Under HWE the coefficient converges
  to `(ploidy - 1)(1 - 2p)`, and at ploidy 2 that limit reproduces the published orthogonal
  dominance coding `(-2p^2, 2pq, -2q^2)` to 1e-16. The change is an exact reparameterization of
  the same model space (total genotypic value invariant to 4e-16).
- **Autopolyploid variance is UNBIASED over unknown phase.** Autopolyploid parental phase is not
  identifiable from dosage. Enumerating every phase configuration consistent with the observed
  dosages and every gamete, the between-locus gamete covariance averages to exactly zero
  (max 1.1e-16) even for COMPLETELY LINKED loci. The locus-sum variance
  (`variance_model = "unlinked_phase_marginalized"`) is therefore the exact expectation given the
  information dosage carries, not an approximation.
- **Allopolyploid subgenome variance** remains exact under disomic inheritance: the
  cross-subgenome block of `R` is identically zero, so the per-subgenome `a'Ra` block-diagonal
  decomposition is exact (recombination-aware when a per-subgenome cM map is supplied).
- **R / C++ parity** for the polyploid dominance kernel: all score columns agree to 4e-16.

Disallowed claims:

- **The multi-trait quadrant is NOT a decomposition of `multi_trait_score` for rank-based index
  methods.** `auto` / `weighted` / `threshold` combine rank-normalized traits, so no linear index
  exists in genetic units and the axes come from reinterpreting trait weights as standardized-unit
  coefficients (`portfolio_basis = "linearized_rank_index"`). Since `auto` promotes to `weighted`
  whenever trait weights are present, this is the common case, and the frontend MUST badge it.
  Only `economic_index` / `desired_gain` give `portfolio_basis = "linear_index"`. This is a
  deliberate override of the 2026-07-25 design's "omit the quadrant for non-linear methods" rule;
  the measured basis for it (`pearson(cross_level, multi_trait_score)` = 0.871 weighted vs 0.873
  economic_index) is recorded as Rev 7 of that design doc.
- **`cross_confidence` and `risk_bin` are within-run only** -- a min-max normalization and
  tertiles of the crosses on screen. A plan always contains roughly one third "high risk"
  regardless of how well it is estimated. Not comparable across runs; not an absolute statement
  about plan quality.
- **The index PEV is a block-diagonal approximation.** Traits are fitted by INDEPENDENT univariate
  ridges, so no cross-trait estimation-error covariance exists; `sum_k w_k^2 PEV_k` is a ranking
  input, not a calibrated interval. Per-trait PEVs are moreover only comparable across traits when
  their residual variances are, and sigma_e^2 is estimated in-sample: a trait whose ridge lambda
  lands on the grid floor can interpolate its training rows and report a near-zero PEV. When one
  trait then carries >90% of the index PEV the run says so (`pev_concentration_note`) and
  `risk_bin` must not be read as an index-wide statement.
- **The autopolyploid variance cannot DISCRIMINATE on linkage phase.** At fixed phase the
  between-locus covariance is real -- about [-1/3, +1/3] for a duplex x duplex tightly linked pair
  -- so two crosses with identical parental dosages but different phase receive identical
  predictions. Variance-based metrics therefore separate autopolyploid crosses less sharply than
  the diploid path, where inbred parents make phase known and the exact `a'Ra` kernel applies.
  Making it exact would need phased polyploid haplotypes, which are NOT implemented.
- **Polyploid dominance is DIGENIC only.** Trigenic and quadrigenic dominance components are not
  modelled.
- **The `poly4x` (AlphaSimR) path reports a SIMULATED variance**, with relative Monte-Carlo SE
  ~ sqrt(2/(n-1)): 29% at 25 progeny, 10% at 200. Below ~50 progeny the sampling error exceeds the
  real spread between many crosses and `poly4x_var` / `poly4x_usefulness` rankings are
  substantially noise. The default is now 200 and per-cross `poly4x_var_se` is reported; no
  precision claim beyond that SE.
- **No selection-superiority claim** for either effort. The portfolio layer is a reporting and
  decision surface, not a change to the ranked merit; the polyploid work is accuracy and
  identifiability, not demonstrated long-term gain.
- **Prior additive/dominance splits are not comparable** with post-0.20.0 ones. The
  reparameterization changes how ridge distributes effects between the components (that split was
  previously decided by the penalty rather than the data), so historical `add_var` / `dom_var`
  numbers should not be compared across the change.

Method note (why the audit was possible at all): `ng_load()` compiles `src/ng_kernels.cpp` via
`Rcpp::sourceCpp`, bypassing `RcppExports`, so the whole suite passed against a changed C++
signature while the generated exports were stale -- installation broke and only `R CMD check`
caught it. `tests/cpp_kernels_equivalence.R` now compares the exported signatures in the kernel
source against `R/RcppExports.R` and was verified to fire on a deliberately reverted signature.

## Residual-Heterozygous-Parent Variance & Message Contract Note (v0.17.1, 2026-08-07)

Adds the exact within-cross additive variance for arbitrary **phased** parents
(`ng_gms_additive_var_general`), generalizing the inbred `a'Ra` kernel to RIL
parents that carry residual heterozygosity; the `parent_type` heterozygosity
governance (`inbred`/`dh`/`ril`); the headless-JSON message contract
(warnings + structured blockers); and the staged JSON runner.

Allowed claims:

- The exact **DH** and **RIL-infinite** within-cross additive variance for phased
  parents (inbred OR residual-heterozygous) is correct for the package's Haldane
  (no-interference) recombination model: validated to Monte-Carlo precision
  against the assumption-free pure-Haldane transmission (DH -0.24%, RIL +0.49%;
  coupling / repulsion / coexisting-Delta,delta / both-het / random configs,
  max |err| 0.31%), and reduces byte-identically to `a'Ra` for inbred parents.
  Derived and adversarially QG-reviewed; see
  `docs/design/residual-het-parent-variance.md`.
- `parent_type` governs the heterozygosity audit independent of the progeny
  target: `dh` (zero-noise floor) and `inbred` block het as a data error; `ril`
  accepts residual het. Wired into `ng_score_crosses`: with phased haplotypes,
  het-parent crosses under `parent_type = "ril"` get the exact variance;
  inbred-parent crosses and every no-phase / DH / inbred path are byte-identical
  to prior behavior.
- The headless JSON runner surfaces backend warnings (`result$warnings[]`) and
  blockers (`ok = false` + `error_message`) instead of crashing; the staged JSON
  runner (`workflow = "stage"`) works and surfaces the same messages per stage.

Disallowed claims:

- No selection-superiority claim. The RAH finite-lineage tail metric was a
  **powered validated NULL** (it does not beat analytical usefulness or mean
  selection); the exact variance correction fixes *accuracy*, not proven
  long-term gain.
- Finite-generation RIL (F2:F3, F3:F4, ...) variance is NOT closed-form here;
  only DH and RIL-infinite are. Finite-RIL remains simulation-only (RAH module).
- The correction is exact only for the Haldane no-interference model and only
  with phased haplotypes; dosage-only input recovers only the diagonal part and
  stays biased low (with the honest advisory).

Frontend-surfacing governance: **experimental capabilities must NOT surface in
the frontend.** The RAH-PMV finite-lineage module (validated NULL) is parked on
its own branch and is NOT in the published package, so it cannot surface.
Anything the backend capability registry marks `status = "experimental"` or
`"guarded"` (e.g. complex-polyploid modules, finite-RIL analytics) must be gated
out of the UI -- the frontend is expected to gate on the backend capability
registry (`ng_backend_capabilities.v1`) for exactly this.

## External Convex Mate-Allocation Benchmark Note (2026-07-03)

REMOVED FROM PACKAGE: the external convex-optimization mate-allocation package (Endelman 2024,
Genetics iyae193) is NO LONGER integrated or depended on -- it never beat the native allocator on
realized gain in any tested setting, so it was dropped (no wrapper, no exports, no DESCRIPTION
dependency). The benchmark record below is retained ONLY as the evidence basis for the
native-allocator-default decision. Method labels coma_oma/coma_useful/coma_* refer to that removed
external comparator in the historical runs.


`tools/run_mate_allocation_study.R` (recurrent mode) benchmarks COMA (Endelman, Genetics 2024,
iyae193; installed via github jendelman/COMA) as an EXACT external mate-allocation baseline
against our native allocators, in OUTBRED populations WITH DOMINANCE (diploid AND autotetraploid,
runMacs founders for realistic LD) -- the setting where mid-parent heterosis, hence COMA's OMA,
is defined (inbred lines have dominance = 0, so OMA == OCS). Each method runs an INDEPENDENT
recurrent genomic-selection lineage from the SAME founders (allocate on GEBV via an augmented
additive+dominance ridge; evaluate on TRUE gv), so differences are attributable to the ALLOCATOR,
not parent selection. Methods: `coma_oma` (COMA exact convex OMA on predicted mean/GPMP),
`ng_useful` (our usefulness = mean + i*within-family SD), `ng_mean` (our mean control). Three
configs: A (8 cycles, meanDD 0.5, 10 reps), B (15 cycles, meanDD 0.5, 8 reps), C (8 cycles,
strong dominance meanDD 1.2, 10 reps); GEBV accuracy ~0.67-0.81 throughout.

Allowed claims (scoped to these configs):
- We are COMPETITIVE WITH the published convex-optimization reference. COMA's exact OMA does NOT
  outperform our native GREEDY usefulness/mean allocation on realized genetic gain in any config
  or ploidy -- final-cycle gains overlap within ~1-2 SE. Under STRONG dominance (config C) our
  simple `ng_mean` significantly LEADS in diploid (3.42 vs COMA 2.97, ~2 SE) and both native
  methods edge COMA in tetraploid. A simpler, OS-agnostic greedy allocator matches a
  second-order-cone convex solver here.
- Our USEFULNESS criterion robustly RETAINS MORE GENETIC VARIANCE than both COMA and mean
  allocation, at equal-or-lower coancestry, in every config and at nearly every cycle (e.g.
  config A final-cycle paired: `ng_useful` +0.037 var vs COMA diploid, +0.058 var and +0.18 gain
  vs COMA tetraploid; config B 15-cycle: `ng_useful` highest var AND lowest coancestry both
  ploidies). This is its designed diversity-preservation advantage, made empirical.
- COMA integration is CERTIFIED (R/45_coma_baseline.R): coancestry-scale K (diagonal ~1/ploidy,
  matching COMA's own potato vignette kinship ~0.27 at 4x), dF as an upper inbreeding-rate bound,
  per-mating contribution cap to yield n discrete crosses. Verified feasible/valid for ploidy 2 and 4.
- VARIABLE FAMILY SIZE (the realistic breeding setting) tested and COMA still not superior. Equal
  family size is a simulation convenience; real programs have unequal families (fertility, seed set,
  germination, crossability). A variable-family study (tools/run_mate_allocation_study.R,
  NG_MATEALLOC_FAMILY_MODE=variable; 6 reps x 6 cyc, diploid+4x) gives each method a fixed total
  progeny budget split across matings: COMA by its optimal contributions (ng_select_coma
  total_progeny=, its native mode -- unequal families), native by OCS selection + score-weighted
  family sizes (ng_allocate_family_sizes). Result: the NATIVE allocator is STRICTLY BETTER than COMA
  -- MORE gain at LOWER coancestry AND more retained variance -- in both ploidies (paired COMA-native:
  diploid -0.18..-0.29 gain at +0.014 coancestry; tetraploid -0.18..-0.24 gain at +0.038 coancestry).
  COMA runs near its ΔF ceiling (coancestry ~0.05), spending inbreeding for one-generation
  contribution-optimality, while the native allocator stays diverse (~0.015) and out-gains it over
  cycles -- diversity preservation beats spending it in recurrent selection. (A 2-rep pilot spuriously
  suggested COMA led; 6 reps reversed it -- pilots are noise.) So across BOTH equal-family and
  variable-family settings, COMA does not beat the native allocator on realized gain.
- USEFULNESS-AWARE OMA tested and NOT an improvement. A 2x2 study {mean, usefulness criterion} x
  {native greedy, COMA exact convex} (main 10 reps x 10 cyc; long 8 reps x 18 cyc) added `coma_useful`
  = our within-family-variance usefulness fed as the per-mating merit into COMA's exact convex ΔF
  solve -- a criterion COMA's own package does NOT offer (it optimizes only the cross mean). Result:
  `coma_useful` does NOT beat COMA's mean OMA (`coma_oma`) on realized gain -- behind in diploid
  (6.11 vs 6.24 at cyc18) and tetraploid (8.17 vs 8.28), while retaining more variance (diploid).
  Its diploid gain gap vs `coma_oma` shrinks monotonically as its variance gap grows (cyc6/12/18:
  gain -0.33/-0.29/-0.13, var +0.03/+0.08/+0.12) -- a directional hint of a crossover well beyond
  18 cycles, unconfirmed. The overall tetraploid top method was the NATIVE greedy usefulness allocator.
  So fusing our variance into COMA's OMA is a diversity lever, not a gain improvement, at <=18-cycle horizons.

Disallowed claims:
- Do NOT claim usefulness (or any method), or the usefulness-aware OMA, wins on GAIN: final-cycle
  gains are within SE and the usefulness-aware OMA is behind COMA's mean OMA on gain at every tested
  horizon. The retained-variance advantage did NOT convert into a gain lead even at 18 cycles.
  Usefulness's benefit at these horizons is VARIANCE/DIVERSITY RETENTION, not demonstrated gain.
- Do NOT claim COMA is inferior. Its continuous OMA is DISCRETIZED to a fixed number of
  equal-family crosses here (the operational "make N crosses" constraint), which does not exercise
  COMA's variable-family-size contribution allocation; and the augmented-ridge dominance
  estimation is weak, under-exercising COMA's heterosis-GPMP edge. The finding is PARITY under a
  fixed-N-crosses operating point, not COMA deficiency.
- Single genome architecture, single dF target, approximate (not exact) coancestry matching across
  allocators. Absolute gains are run-specific; the PARITY-on-gain and usefulness-VARIANCE-RETENTION
  patterns are the takeaways, both directional (8-10 reps).

## Metric Variants in Recurrent Selection Note (2026-07-03)

An 11-method recurrent RIL study (`tools/run_ril_breeding_program_benchmark.R`, 10 reps x
25 cycles, realistic config) raced ALL the cross-scoring metric variants in the multi-cycle
loop -- `uc_pmv_ocs`, `uc_vpm_ocs`, `pmv_ocs`, `vpm_ocs`, `var_simple_ocs`, `mean_ocs`,
`var_complex_ocs` -- alongside `var_complex_evolution`, `strategy_balanced`,
`simple_usefa_select`, and `package_user_api`. This is the RECURRENT counterpart to the
single-generation Metric-Merit Study Note below (0 skips, all 11 x 25 cycles).

Allowed claims (scoped to this config):
- SINGLE-GENERATION RANKING != RECURRENT REALIZED GAIN. `var_simple` -- the WEAKEST metric in
  the single-generation ranking study (0.443, last) -- is the SECOND-BEST recurrent performer
  (cycle-25 gain 11.17 +/- 0.25, ~4 SE above the usefulness cluster) and retains by far the
  most polymorphism (0.161 vs ~0.03-0.04). Mechanism: `var_simple` is a relationship-DISTANCE
  criterion that favours unrelated/complementary parents, preserving variance -- a poor
  single-generation merit predictor but a diversity-preserving one that pays off over cycles.
- The PMV/VPM/uc usefulness variants cluster tightly on recurrent gain (~9.6-10.2, overlapping
  in SE) -- they rank crosses nearly identically -- and all deplete diversity hard (He ~0.01,
  poly ~0.03), plateauing lower than diversity-preserving methods.
- `mean_ocs` and `simple_usefa_select` are lowest with the most depleted diversity.
- Prediction accuracy is similar across all metrics (~0.84-0.87): the differences are about
  variance management over cycles, not prediction quality.
- Practical guidance: for SHORT-TERM single-generation selection use PMV/VPM usefulness
  (`var_simple` is weak); for LONG-TERM recurrent selection use explicit diversity management
  on a good merit metric (the strategy dial, cycle-25 11.63) rather than relying on
  `var_simple`'s incidental diversity effect or a pure-usefulness metric that exhausts variance.

Disallowed claims:
- Single config, 10 reps -> directional. Do not claim `var_simple` is a good MERIT metric
  (it is not; its recurrent edge is diversity preservation, and the strategy dial dominates it).
  Absolute gains are run-specific; the single-gen-vs-recurrent DISSOCIATION and the relative
  orderings are the takeaways.

## Optimizer Axis in Recurrent Selection Note (2026-07-04)

`tools/run_method_optimizer_metric_study.R` (new-user/breeder recurrent RIL study, 13 arms x
10 reps x 25 cycles, driven through `ng_run_cross_prediction`, evolution the standard optimizer)
includes an explicit optimizer axis on a FIXED metric (`uc`/`pmv`): `uc_pmv_ocs` (evolution),
`uc_pmv_mip` (mip_contribution), `uc_pmv_greedy` (greedy_local).

Allowed claims (config-scoped, directional):
- On REALIZED recurrent gain, **evolution matches exact MIP and both clearly beat greedy** --
  cycle-25 gain evolution 7.20 +/- 0.30, MIP 7.07 +/- 0.44, greedy 6.78 +/- 0.30 (n=10). The
  evolution-vs-MIP paired difference (0.13 +/- 0.29; evolution ahead in 6/10 reps) is NOT
  significant.
- Evolution is the sound DEFAULT for real/recurrent/large programs: level with exact MIP on
  realized gain, warm-started from greedy (never worse), and it scales where MIP is
  size/time-guarded. This complements the single-shot objective bake-off
  (`tools/run_optimizer_benchmark.R`) where mip_contribution ranks first on ACHIEVED OBJECTIVE
  and evolution is a close second.

Disallowed claims:
- Do NOT claim evolution "beats" MIP on realized recurrent gain (the difference is within noise).
  The optimizer changes allocation only, not prediction accuracy. Absolute gains are run-specific.

## Metric-Merit Study Note (2026-07-03)

A controlled AlphaSimR study (`tools/run_metric_merit_study.R`, 8 conditions x 15 reps:
training size {120, 400} x heritability {0.25, 0.60} x architecture {oligogenic 5 QTL/chr,
polygenic 40 QTL/chr}; runMacs GENERIC for realistic marker-QTL LD; DH progeny) graded the
cross-scoring metrics against an UNBIASED true-effect yardstick -- every candidate cross is
re-scored with the true QTL effects (true usefulness = mid-parent true GV + i * true
within-family SD) so a metric cannot flatter itself. Measures: Spearman ranking accuracy and
top-K true merit (allocation-free, so they isolate metric quality, not allocator quality).

Allowed claims:
- The package's UC-VPM usefulness and SimpleMating's additive usefulness produce IDENTICAL
  cross rankings (max |diff| 1e-4 across all 120 runs) -- an independent cross-validation of
  the package's recombination-variance usefulness.
- UC-PMV and UC-VPM are near-equivalent in cross RANKING (overall accuracy 0.469 vs 0.471);
  PMV's posterior term shifts predicted magnitudes, not order, so its value is in calibration.
- Usefulness (UC-PMV/VPM) beats mean-only only when the within-family variance is meaningful
  AND well-estimated: clearest for oligogenic traits with adequate training + heritability
  (UC-PMV minus mean up to +0.058 +/- 0.006). For polygenic traits usefulness ~= mean; for
  oligogenic + poorly-estimated effects (small training, low h2) the noisy variance term can
  slightly reduce ranking accuracy (-0.020 +/- 0.013).
- var_simple is the weakest merit metric (overall 0.443); it is a diversity/QC proxy, not a
  cross-merit criterion.
- UC-PMV is a safe, near-best default across all conditions tested.

Disallowed claims:
- Do not claim any single metric is universally best; merit is condition-dependent (above).
- These are SINGLE-GENERATION ranking results on DH progeny under one demography (runMacs
  GENERIC); do not extend to multi-cycle realized gain without the recurrent study.

## RIL Recurrent-Study Note (2026-07-03)

A 10-rep x 25-cycle recurrent RIL study (`tools/run_ril_breeding_program_benchmark.R`,
realistic config: 100 founders -> 60 parents, 40 crosses, ~2500 SNP / 300 QTL, 500 training
lines) compared six methods: `var_complex_ocs` (MIP), `var_complex_evolution`, the balanced
strategy dial (`strategy_balanced`), `mean_ocs`, `simple_usefa_select` (SimpleMating, RIL),
and `package_user_api` -- the package driven through its top-level API (`ng_run_cross_prediction`)
exactly as an end user would, the same black-box way SimpleMating is invoked. All six ran all
25 cycles (0 skips). Cycle-25 gain (mean +/- SE, n=10): strategy 11.60 +/- 0.22 > OCS 10.10 +/-
0.19 = user-API 10.10 +/- 0.12 > evolution 9.89 +/- 0.22 > mean 9.49 +/- 0.24 > SimpleMating
9.26 +/- 0.23.

Allowed claims (scoped to this config):
- USER-API VALIDATION: the package driven via its top-level API (`package_user_api`)
  reproduces the internal `var_complex_ocs` path on long-term realized gain (both cycle-25
  10.10; mean absolute per-cycle difference 0.57 on a ~10 scale, converging by ~cycle 15).
  The small early-cycle gap is expected -- the user API trains effects on the candidate
  parents' own phenotypes and scores with `assume_inbred = FALSE`, the realistic user setting.
- The gain-diversity strategy dial (balanced) sustains higher LONG-TERM gain (11.60 +/- 0.22)
  than the fixed-light-lambda methods by preserving genetic variance; those plateau by ~cycle 15.
- SimpleMating (RIL) is competitive at a single generation but plateaus lowest in recurrent
  selection under its default culling.

Disallowed claims:
- The evolution optimizer and OCS/MIP are WITHIN NOISE on long-term realized gain (this run
  OCS 10.10 vs evolution 9.89; the earlier 5-method run had evolution 9.92 vs OCS 9.66 -- the
  ordering flips within ~1 SE). Do NOT claim evolution is superior for realized gain; its
  validated advantage is on achieved objective in the optimizer bake-off, not realized gain.
- Single config, 10 reps -> directional, not definitive. The OCS/SimpleMating baselines used
  fixed/untuned diversity settings; this is NOT a claim that the dial beats a delta-F-tuned
  OCS. Use `target_coancestry=` (constrained OCS) to compare methods at matched inbreeding.

## Optimizer Improvements Note (2026-07)

Added a native evolutionary (memetic genetic-algorithm) mate-allocation optimizer
`ng_evolutionary_mate_allocation()` / `ng_optimize_mating_plan(method = "evolution")`
that optimizes the same objective as the OCS/MIP path
(`ng_plan_objective_contribution`), warm-started from `greedy_local` with elitism (so
it is never worse than greedy) and a C++ `ng_local_swap` memetic hill-climb. Also
hardened the existing optimizers: the C++ `ng_local_swap_cpp` kernel now models the
`lambda_parent_use` term (no more slow R fallback) with deterministic tie-breaking;
the lpSolve MIP paths are size/time-guarded and fall back to `greedy_local`
(`summary$mip_fallback`) instead of hanging; and `greedy_local`/`repair_local` now
actively enforce `min_unique_parents` rather than only warning.

Allowed claim:
- On a fixed-input optimizer bake-off (`tools/run_optimizer_benchmark.R`), the
  evolutionary optimizer is feasible, reproducible under a fixed seed, and achieves
  an objective >= greedy/repair and competitive with the MIP path.

Disallowed claim:
- Do not claim the evolutionary optimizer is superior for realized genetic gain. The RIL
  Recurrent-Study Note (2026-07-03) shows it is WITHIN NOISE of OCS on long-term realized gain
  (the ordering flips between runs within ~1 SE). Its validated strength is achieved objective
  in the fixed-input optimizer bake-off, not realized gain in a breeding program.

## Current Package Artifact Note (v0.4.0)

Version 0.4.0 rebuilds the polyploid path correctly (ploidy-general, additive default,
dominance OPTIONAL) and removes the COMA integration (benchmarked, never beat the native
allocator; a third-party package the backend no longer uses). New/changed:
- Correct allele-frequency-based polyploid GRM `ng_polyploid_grm(method = "vanraden" | "yang")`
  (replaces the old ploidy-midpoint kinship) and digenic `ng_polyploid_dominance_grm()`;
  ploidy-aware `ng_polyploid_qc()`.
  > **CORRECTION (2026-08-16, v0.20.0).** "Replaces the old ploidy-midpoint kinship" OVERSTATED
  > what v0.4.0 shipped. The correct GRM replaced the midpoint form on the `ng_polyploid_design_*`
  > path only; the `poly4x` (AlphaSimR) path went on building `parent_kinship` from
  > `ng_poly4x_parent_relationship` and handing it to `ng_poly4x_ocs` as the coancestry the
  > diversity penalty acts on, until v0.20.0. The midpoint form centres dosage on `ploidy/2`, i.e.
  > assumes every allele frequency is 0.5: on 200 UNRELATED autotetraploids it reports a mean
  > off-diagonal of 1.4724 (range 1.203..1.718, diag 2.4724) where VanRaden gives -0.0050 (range
  > -0.148..0.145, diag 0.9952), and -- the part that matters for OCS, which consumes an ORDERING
  > -- the Spearman rank correlation between the two sets of pair values is only 0.565. Any
  > `poly4x` diversity/coancestry result produced before v0.20.0 was computed on that matrix.
- Additive+dominance marker effects `ng_fit_polyploid_effects()` and value prediction
  `ng_predict_polyploid_value(type = "genotypic" | "breeding")` for clonal-crop selection.
- Dominance-aware cross scoring `ng_score_crosses_poly_dominance()` (heterosis-inclusive mean +
  within-family additive+dominance variance, double reduction) via a progeny-moment table and a
  C++ kernel `ng_poly_dominance_scores_cpp`; validated against gamete simulation.
- Any-ploidy one-call `ng_design_crosses_poly(dominance =, gain =, double_reduction =, grm_method =)`.
- Examples 29-31, a vignette "Polyploid and dominance-aware design" section, and capability-registry
  entries `polyploid_grm_qc` / `polyploid_dominance_design`.
- REMOVED: `ng_select_coma` / `ng_coma_*` exports, the COMA/CVXR optional dependencies, and the
  COMA docs/examples. Evidence retained in the benchmark notes as the native-default rationale.

## Previous Package Artifact Note (v0.3.13)

Version 0.3.13 adds runnable examples for the practical workflow gaps raised
before independent user testing: duplicate QC reporting/removal, exact
phenotype/genotype/direction/map matching, base-pair marker-map positions,
cross-number sweeps, priority workbook and figure outputs, allocation-method
comparison, posterior robust mate allocation, and full-posterior PMV shortlist
reruns. The new examples are numbered
`14_qc_duplicate_removal_and_reporting.R` through
`20_full_posterior_pmv_shortlist.R`.

## Previous Package Artifact Note (v0.3.12)

Version 0.3.12 adds separate, runnable examples for each user-facing
`prediction_mode` and `multi_trait_method` path. The new examples are numbered
so users can work through them in order:
`08_prediction_mode_trait_by_trait.R`,
`09_prediction_mode_index_as_trait.R`, `10_multitrait_method_auto.R`,
`11_multitrait_method_weighted.R`, `12_multitrait_method_economic_index.R`,
and `13_multitrait_method_desired_gain.R`. Each script keeps the relevant
parameter block explicit and checks that the expected wrapper path was used.

## Previous Package Artifact Note (v0.3.11)

Version 0.3.11 exposes the PMV/posterior workflow arguments directly in
`ng_run_cross_prediction()`: `method_varPMV`, `ril_mode`,
`run_posterior_prediction`, `posterior_method`, `nIter`, `burnIn`, and
`use_parallel`. The arguments are wired into the actual workflow rather than
recorded as metadata only. `method_varPMV = "full_posterior"` requests
`ng_fit_ridge_effects(return_beta_cov_full = TRUE)` and passes the full beta
covariance into `ng_score_crosses()`. Trait values then use the
`dh_pmv_var_full_posterior` column through the `*_pmv_used` column. When
`run_posterior_prediction = TRUE`, the wrapper calls
`ng_fit_ridge_effects_posterior()` and `ng_posterior_cross_predict()` and
returns the per-trait tables in `result$posterior_predictions`.

The new example at `inst/examples/07_trait_value_metric_parameter_guide.R`
shows when to use each trait-value metric, the exposed PMV/posterior settings,
and a small posterior-prediction row.

## Previous Package Artifact Note (v0.3.10)

Version 0.3.10 promotes `lpSolve` from `Suggests` to `Imports` because
`mip_linear`, `mip_contribution`, and `optimizer = "auto"` MIP dispatch are
user-facing OCS optimizer paths. The optimizer guide at
`inst/examples/06_optimizer_parameter_guide.R` therefore expects all optimizer
rows to run after package installation rather than skipping MIP rows.

## Previous Package Artifact Note (v0.3.9)

Version 0.3.9 keeps the v0.3.8 wrapper API and adds a dedicated package OCS
example for `allocation_method = "ocs"` at `inst/examples/05_ocs_user_run.R`.
The example exposes the OCS-specific capacity, kinship, optimizer, and penalty
settings separately from native AlphaMate-style and external AlphaMate examples.

## Previous Package Artifact Note (v0.3.8)

Version 0.3.8 keeps the v0.3.7 user-facing wrapper API and adds a dedicated
example for `allocation_method = "alphamate_executable"` at
`inst/examples/04_alphamate_executable_user_run.R`. The example is guarded so it
prints setup instructions when `NG_ALPHAMATE_EXE` does not point to an installed
AlphaMate binary, and it exposes the executable-only parameters separately from
the native AlphaMate-style allocator.

## Previous Package Artifact Note (v0.3.7)

Version 0.3.7 keeps the v0.3.6 explicit input-matching contract and adds
user-facing method selection inside `ng_run_cross_prediction()`. The wrapper now
accepts `trait_value_metric = "var_complex"` as a package-native,
PopVar-inspired complex usefulness metric without requiring PopVar, and exposes
`allocation_method = "ocs"`, `"alphamate_style"`, or
`"alphamate_executable"`. The native AlphaMate-style path uses
`ng_alphamate_style_select()`; the executable path uses `ng_select_alphamate()`
and fails clearly when the configured AlphaMate executable is unavailable.

A synthetic source test covers `var_complex` direction-aware trait values,
native AlphaMate-style allocation, and the missing-executable error path. The
test is also runnable against the installed source tarball through
`NG_TEST_INSTALLED_LIB`.

The method-comparison evidence below remains the retained validation context;
v0.3.7 changed the user workflow surface and wrapper exposure, not the
underlying validated method-performance claims.

## Previous Package Artifact Note (v0.3.6)

Version 0.3.6 keeps `ng_run_cross_prediction()` as the user-facing workflow for
phenotype, genotype, map, and trait-direction files and adds an explicit
matching contract for independent user testing. The wrapper now accepts
separate phenotype and genotype ID columns, direction-file trait/phenotype
column/direction columns, marker-map marker/chromosome/base-pair position
columns, `map_position_unit = "bp"`, and `bp_per_cm`. It canonicalizes IDs
internally, runs package QC and putative duplicate removal, requires exact
phenotype/genotype parent matching after QC, requires exact marker-map/genotype
marker matching, reorders phenotype and marker-map tables to the genotype
order, keeps `pos_bp`, derives `pos_cm`, and returns `result$input_match_audit`.

A synthetic matching-contract test covers shuffled phenotype rows, different
phenotype/genotype ID column names, custom direction-file column names,
base-pair marker-map positions, exact marker order, missing marker-map errors,
and missing phenotype-trait errors. The prior v0.3.5 workflow test still covers
both `prediction_mode = "trait_by_trait"` and `prediction_mode =
"index_as_trait"` and verifies duplicate removal plus PNG figure output.

## Previous Package Artifact Note (v0.3.5)

Version 0.3.5 added `ng_run_cross_prediction()` as the user-facing workflow for
phenotype, genotype, map, and trait-direction files. The wrapper runs package
QC, putative duplicate-genotype handling, trait-by-trait marker-effect
estimation, cross prediction, multi-trait objective construction, OCS
allocation, priority ranking, and optional workbook/figure output.

## Previous Package Artifact Note (v0.3.4)

The current independent-testing artifact is
`dist/nextgenCrossDesign_0.3.4.tar.gz`. Version 0.3.4 adds package-embedded
putative duplicate genotype removal through `ng_preflight_input_tables()` with
`putative_duplicate_action = "remove"`, returns cleaned analysis tables in
`qc$cleaned_tables`, and exposes the feature through
`ng_backend_capability_registry()`.

The data-rich package-direct workflow was verified from the source tarball on
2026-06-20. It removed 2 duplicate parents inside package QC, scored 12,403
candidate crosses from 158 retained parents, selected 100 crosses, produced the
10/25/35/30 priority-tier split, wrote 8 figures, and wrote
`breeder_crossing_plan.xlsx`.

This note updates package readiness and workflow verification only. The older
method-comparison and AlphaMate validation evidence below remains the retained
statistical-validation context for the model families described there.

> **v0.2.0 STATUS — 6-rep follow-up (2026-05-24): doubled reps tighten
> variance but do not push per-band Welch p below α=0.05 against best
> AlphaMate. Directional headline holds at sign-test p = 0.063.**
>
> The 6-rep follow-up (7 sizes × 6 reps × 3 cycles × 5000 markers × 400
> training × up-to-10 methods) was run at
> `results/revalidation_v0_2_0_5k_6rep_*`. AlphaMate's exact binary
> crashed reproducibly at p=70 (rep 1 cycles 1-2 on two separate
> attempts, with and without `alphamate_opt60`). Removing AlphaMate
> let p=70 and p=80 complete cleanly, confirming AlphaMate as the
> failure mode. Final coverage: 7/7 sizes with internal methods + 5/7
> sizes (p=20–60) with all three AlphaMate target degrees.
>
> ## Headline at n=6 reps (5 bands with AlphaMate)
>
> | Parents | Frontier mean ± SD | Best AlphaMate (target, mean ± SD) | Δ %     | Welch p |
> |---:|---|---|---:|---:|
> | 20 | 7.293 ± 0.873 | `opt30` 6.808 ± 0.948 | **+7.13%**  | 0.38 |
> | 30 | 7.649 ± 1.028 | `opt45` 7.324 ± 0.860 | **+4.43%**  | 0.57 |
> | 40 | 8.263 ± 0.862 | `opt30` 7.728 ± 0.809 | **+6.93%**  | 0.29 |
> | 50 | 8.705 ± 0.804 | `opt30` 7.998 ± 0.760 | **+8.83%**  | 0.15 |
> | 60 | 9.155 ± 1.215 | `opt30` 8.809 ± 1.016 | **+3.93%**  | 0.61 |
> | 70 | 9.175 ± 1.318 | (AlphaMate crashed)   | —           | —    |
> | 80 | 9.495 ± 1.351 | (AlphaMate crashed)   | —           | —    |
>
> **Frontier wins 5/5 bands where AlphaMate ran.** Sign-test under
> H₀=0.5: p = 2 × 0.5⁵ = 0.063. No per-band Welch reaches α=0.05; the
> closest is p=50 (Δ +8.83%, Welch p=0.149). **Doubling reps from 3 to
> 6 did NOT push per-band p-values below α=0.05 against the best
> AlphaMate target.** The noise floor (SDs of 0.8–1.3 on means of 7–9)
> is too high relative to the 4–9% effect sizes; detection at α=0.05
> would need ~12+ reps (i.e., another doubling).
>
> ## Top-3 consistency at n=6 reps (7 bands, methods that ran at each)
>
> | Method | Top-3 finishes |
> |---|---:|
> | `ng_meta_router_ocs10_lps2` | **5/7** (down from 6/7 at 3 reps) |
> | `ng_pmv_blend_balanced_ocs10_lps2` | 4/7 (up from 1/7 at 3 reps) |
> | `ng_meta_portfolio_ocs10_lps2` | 4/7 (was 5/7 at 3 reps) |
> | `ng_frontier_policy_ocs10_lps2` | 4/7 (unchanged from 3 reps) |
> | `var_simple_ocs10_lps2` | 2/7 |
> | `ng_recomb_gebv_ocs10_lps2` | 2/7 (was 3/7) |
> | `var_simple_topn` | 0/7 (was 1/7) |
> | `alphamate_opt30 / opt45 / opt60` | 0/7 each (unchanged) |
>
> `ng_meta_router` is still the most consistent winner. `ng_pmv_blend_balanced`
> jumped from 1/7 to 4/7 — the 3-rep grid underestimated it. The frontier
> policy holds steady at 4/7.
>
> ## AlphaMate stability note
>
> The exact AlphaMate binary at `external/AlphaMate/binaries/AlphaMate.exe`
> (with `NG_ALPHAMATE_RUNTIME_PATH=C:/Python/Lib/site-packages/torch/lib`)
> crashed mid-cycle for p=70 on two separate 6-rep attempts. The 3-rep
> grid completed cleanly at the same scale, so this looks like
> intermittent memory pressure or AlphaMate I/O race condition that
> manifests when the per-scenario reps × cycle count is high. Worth a
> separate binary-level audit (memory usage, runtime DLL version,
> AlphaMate's stdout/stderr) before publishing any "vs AlphaMate at
> p ≥ 70" claim. The 3-rep × 3-cycle 5K refinement grid below has
> stable AlphaMate evidence for all 7 bands; treat it as the
> authoritative "vs AlphaMate" run.
>
> ## Older context: the 5K refinement (3 reps × 3 cycles × 5K markers)
>
> The metric correctness bugs are fixed (see `docs/V0_1_0_MIGRATION.md` and
> `docs/V0_2_0_RELEASE.md`). Three grids were run on v0.2.0 in order of
> increasing scale:
>
> 1. A **replicated grid** (4 parent sizes × 3 reps × 1 cycle × 2000
>    markers × 300 training) at `results/revalidation_v0_2_0_replicated_*`.
> 2. A **full grid** (7 parent sizes × 3 reps × 2 cycles × 2500 markers
>    × 400 training × 8 methods including `alphamate_opt45`) at
>    `results/revalidation_v0_2_0_full_*`.
> 3. A **5K-marker refinement** (7 parent sizes × 3 reps × **3 cycles**
>    × **5000 markers** × 400 training × 10 methods including
>    `alphamate_opt30`, `opt45`, `opt60`) at
>    `results/revalidation_v0_2_0_5k_*`. **This matches the exact v0.0.x
>    configuration** the original headline was built on.
>
> Analysis reproduced via `tools/analyze_revalidation_grid.R` (env var
> `NG_ANALYSIS_PREFIX` selects the grid).
>
> ## Headline (5K refinement — full v0.0.x scale)
>
> **The v0.0.x claim "ng_frontier_policy_ocs10_lps2 beat the best real
> AlphaMate target in 7/7 parent-size bands" IS reproduced** at the exact
> v0.0.x configuration on the corrected metrics. The frontier beats the
> best of {opt30, opt45, opt60} in 7/7 bands; effect sizes are smaller
> (+0.56% to +11.11%) than the 2-cycle 2500-marker estimate (+3.71% to
> +13.90%), but the directional claim is robust:
>
> | Parents | Frontier mean ± SD | Best AlphaMate (target, mean ± SD) | Δ % | Welch p |
> |---:|---|---|---:|---:|
> | 20 | 7.178 ± 0.509 | `opt30` 6.461 ± 0.428 | **+11.11%** | 0.14 |
> | 30 | 7.217 ± 0.569 | `opt30` 7.064 ± 1.094 | **+2.17%** | 0.84 |
> | 40 | 7.965 ± 0.695 | `opt45` 7.637 ± 0.863 | **+4.30%** | 0.64 |
> | 50 | 8.315 ± 0.783 | `opt30` 7.817 ± 0.767 | **+6.38%** | 0.48 |
> | 60 | 8.737 ± 1.012 | `opt30` 8.582 ± 1.040 | **+1.81%** | 0.86 |
> | 70 | 8.673 ± 0.903 | `opt30` 8.435 ± 0.608 | **+2.83%** | 0.73 |
> | 80 | 8.792 ± 0.576 | `opt30` 8.743 ± 0.297 | **+0.56%** | 0.90 |
>
> No per-band Welch p reaches α=0.05 at 3 reps × 3 cycles. The cross-band
> sign-test (7/7 wins under H₀=0.5) gives p = 2 × 0.5⁷ = **0.016** —
> the consistency of the directional pattern is statistically reliable
> even though no single-band gap is.
>
> ## Important correction to the previous (2-cycle / `opt45`-only) claim
>
> The earlier 2-cycle full grid (which tested only `alphamate_opt45`)
> reported AlphaMate as rank 8/8 in every band and frontier deltas of
> +3.7% to +13.9%. The 5K refinement shows that conclusion was partially
> wrong: **`opt30` is the best AlphaMate target in 6/7 bands** (with
> `opt45` winning only at p=40). With `opt30` in the picture, AlphaMate
> is much more competitive — its rank ranges from 4 (at p=80) to 8 (at
> p=20). The "AlphaMate is uniformly worst" framing from the prior grid
> was an artifact of testing only `opt45`.
>
> The previous head-to-head numbers used `opt45` as the AlphaMate
> reference (the only one in that grid). The headline numbers above use
> the best AlphaMate target per band, which is the correct v0.0.x
> comparison.
>
> ## Per-parent-size rankings at full v0.0.x scale (5K refinement)
>
> Top 5 + AlphaMate ranks by mean top-10 GV at cycle 3 (n=3 reps):
>
> | Parents | Rank 1 | Rank 2 | Rank 3 | Rank 4 | Best AlphaMate | Worst AlphaMate |
> |---:|---|---|---|---|---|---|
> | 20 | `ng_frontier_policy` 7.178 | `ng_recomb_gebv` 7.178 | `var_simple_topn` 7.117 | `ng_pmv_blend_balanced` 7.008 | `opt30` 6.461 (#8) | `opt60` 6.153 (#10) |
> | 30 | `ng_meta_router` 7.675 | `ng_recomb_gebv` 7.540 | `var_simple_ocs10_lps2` 7.531 | `ng_pmv_blend_balanced` 7.485 | `opt30` 7.064 (#8) | `opt60` 6.613 (#10) |
> | 40 | `ng_meta_portfolio` 8.051 | `ng_recomb_gebv` 7.992 | `ng_meta_router` 7.967 | `ng_frontier_policy` 7.965 | `opt45` 7.637 (#7) | `opt60` 6.725 (#10) |
> | 50 | `ng_meta_portfolio` 8.366 | `ng_frontier_policy` 8.315 | `ng_meta_router` 8.315 | `var_simple_ocs10_lps2` 8.256 | `opt30` 7.817 (#7) | `opt60` 7.336 (#10) |
> | 60 | `ng_meta_portfolio` 8.823 | `ng_meta_router` 8.823 | `ng_frontier_policy` 8.737 | `ng_recomb_gebv` 8.717 | `opt30` 8.582 (#6) | `opt60` 7.656 (#10) |
> | 70 | `ng_pmv_blend_balanced` 8.902 | `ng_meta_portfolio` 8.884 | `ng_meta_router` 8.884 | `ng_frontier_policy` 8.673 | `opt30` 8.435 (#8) | `opt60` 7.523 (#10) |
> | 80 | `ng_frontier_policy` 8.792 | `ng_meta_portfolio` 8.792 | `ng_meta_router` 8.792 | `opt30` 8.743 | `opt30` 8.743 (#4) | `opt60` 8.075 (#10) |
>
> **Frontier policy is rank 4/10 or better in 6/7 bands** (the exception
> is p=30 where it drops to rank 7). The best AlphaMate target is
> `opt30` in 6/7 bands, `opt45` only at p=40, and `opt60` is the worst
> AlphaMate target in **every band** — a robust finding that argues
> against using `opt60` as a default on AlphaSimR DH/RIL setups.
>
> ## Frontier vs best non-frontier internal method (5K refinement)
>
> | Parents | Frontier mean ± SD | Best non-frontier (mean ± SD) | Δ % | Welch p |
> |---:|---|---|---:|---:|
> | 20 | 7.178 ± 0.509 | `ng_recomb_gebv_ocs10_lps2` 7.178 ± 0.509 | +0.00% | 1.00 |
> | 30 | 7.217 ± 0.569 | `ng_meta_router_ocs10_lps2` 7.675 ± 0.848 | +6.35% | 0.49 |
> | 40 | 7.965 ± 0.695 | `ng_meta_portfolio_ocs10_lps2` 8.051 ± 0.729 | +1.07% | 0.89 |
> | 50 | 8.315 ± 0.783 | `ng_meta_portfolio_ocs10_lps2` 8.366 ± 0.716 | +0.62% | 0.94 |
> | 60 | 8.737 ± 1.012 | `ng_meta_portfolio_ocs10_lps2` 8.823 ± 0.879 | +0.98% | 0.92 |
> | 70 | 8.673 ± 0.903 | `ng_pmv_blend_balanced_ocs10_lps2` 8.902 ± 0.608 | +2.64% | 0.74 |
> | 80 | 8.792 ± 0.576 | `ng_meta_portfolio_ocs10_lps2` 8.792 ± 0.576 | +0.00% | 1.00 |
>
> The frontier is tied or within 1.1% of the best non-frontier method in
> 5/7 bands (p=20, 40, 50, 60, 80). At p=30 a 6.35% gap to
> `ng_meta_router` appears but is not significant (Welch p=0.49).
>
> ## Top-3 consistency across all 7 bands (5K refinement)
>
> | Method | Top-3 finishes (of 7) |
> |---|---:|
> | `ng_meta_router_ocs10_lps2` | **6** |
> | `ng_meta_portfolio_ocs10_lps2` | 5 |
> | `ng_frontier_policy_ocs10_lps2` | 4 |
> | `ng_recomb_gebv_ocs10_lps2` | 3 |
> | `var_simple_topn` | 1 |
> | `var_simple_ocs10_lps2` | 1 |
> | `ng_pmv_blend_balanced_ocs10_lps2` | 1 |
> | `alphamate_opt30` | **0** |
> | `alphamate_opt45` | **0** |
> | `alphamate_opt60` | **0** |
>
> ## Operational defaults on v0.2.0 (updated to reflect 5K refinement)
>
> 1. **`ng_meta_router_ocs10_lps2` is the most consistent performer**
>    (6/7 top-3 finishes) — promoted over `ng_meta_portfolio` and
>    `ng_frontier_policy` based on the 5K-refinement evidence.
> 2. **`ng_meta_portfolio_ocs10_lps2`** is the second-most-consistent
>    (5/7 top-3) — recommended fallback.
> 3. **`ng_frontier_policy_ocs10_lps2`** remains evidence-backed: rank
>    4/10 or better in 6/7 bands, ties for best at p=20 and p=80, beats
>    best AlphaMate in 7/7 bands (sign-test p=0.016). Continue to ship
>    it as a band-aware default for users who want explicit per-band
>    method selection.
> 4. **`ng_recomb_gebv_ocs10_lps2`** is a strong third choice (3/7
>    top-3, tied for #1 at p=20).
> 5. **AlphaMate target degree matters a lot.** `opt30` is the
>    best-of-AlphaMate in 6/7 bands and competitive (within 11% of
>    frontier) in all 7. `opt45` is mid-pack. **`opt60` is the worst
>    method in every band** — do not use without re-checking AlphaMate
>    tuning.
> 6. **`var_simple_topn`** is the always-available zero-tuning baseline
>    — beats best AlphaMate in 5/7 bands at 5K scale.
>
> ## Limits of the 5K refinement (final-scale)
>
> - **3 reps × 3 cycles = 9 realisations per cell** is still not enough
>   to detect per-band differences at α=0.05 (Welch p ≥ 0.14). The
>   cross-band sign-test (7/7 wins) is the only test that reaches
>   significance.
> - **Replicate variance is large.** SDs are 0.3–1.1 GV units on means of
>   7–9. Practical implication: a single-rep run cannot reliably rank
>   methods; always run ≥ 3 reps for any method-promotion decision.
> - **AlphaMate runtime configuration**: `NG_ALPHAMATE_RUNTIME_PATH=
>   C:/Python/Lib/site-packages/torch/lib`. The 8/10 / 9/10 / 10/10
>   ranks for `opt45` / `opt30` / `opt60` are reproducible on this
>   binary + runtime config but should be re-tested if the binary or
>   `libiomp5md.dll` version changes.
>
> ## Older 2-cycle / `opt45`-only grid (preserved for context, superseded)
>
> | Parents | Frontier mean ± SD | AlphaMate opt45 mean ± SD | Δ %     | Welch p |
> |---:|---|---|---:|---:|
> | 20 | 6.492 ± 1.259      | 6.260 ± 1.267             | **+3.71%** | 0.83 |
> | 30 | 6.659 ± 1.069      | 5.934 ± 1.161             | **+12.22%**| 0.47 |
> | 40 | 6.389 ± 0.807      | 6.060 ± 0.923             | **+5.44%** | 0.67 |
> | 50 | 7.069 ± 1.370      | 6.364 ± 0.637             | **+11.08%**| 0.48 |
> | 60 | 6.911 ± 0.759      | 6.067 ± 0.597             | **+13.90%**| 0.21 |
> | 70 | 6.763 ± 0.759      | 6.218 ± 0.891             | **+8.78%** | 0.47 |
> | 80 | 6.851 ± 0.718      | 6.178 ± 0.724             | **+10.90%**| 0.32 |
>
> **None of the per-size differences reach α=0.05 at 3 reps × 2 cycles
> (Welch p ≥ 0.21).** The headline is a *directional* claim across 7/7
> bands, not a per-band hypothesis test. Sign-test (frontier wins in 7
> of 7 bands under H₀ = 0.5) gives p = 2 × 0.5⁷ = 0.016 — the
> *consistency* of the directional pattern is statistically reliable
> even though no single band's gap is.
>
> **Multi-cycle compounding is the load-bearing variable.** On the
> 1-cycle replicated grid the frontier was rank 6/8 at p=20 and p=40
> (the AlphaMate-free comparison); on the 2-cycle full grid it is rank
> 4/8 or better in all 7 bands. The v0.0.x evidence used 3 cycles —
> reproducing the directional claim required restoring multi-cycle
> selection, which compounds the band-dispatcher's modest per-cycle
> advantage.
>
> **AlphaMate (opt45) ranks 8/8 in every single parent size** — beaten
> by every internal method including the simplest baseline
> `var_simple_topn`. This is a much bigger story than the headline: real
> AlphaMate, with its default tuning at target degree 45, is the worst
> performer of the 8 methods tested. The frontier-vs-AlphaMate Δ is
> driven as much by AlphaMate under-performing as by the frontier
> winning. Worth a separate investigation: whether AlphaMate's defaults
> on this AlphaSimR setup are mis-tuned, whether other target degrees
> (opt30, opt60) recover, and whether `NG_ALPHAMATE_RUNTIME_PATH`
> affects results.
>
> **Frontier policy vs the best non-frontier internal method:**
>
> | Parents | Frontier mean ± SD | Best non-frontier (mean ± SD)             | Δ %    | Welch p |
> |---:|---|---|---:|---:|
> | 20 | 6.492 ± 1.259      | `ng_pmv_blend_balanced_ocs10_lps2` 6.739 ± 1.423 | +3.80% | 0.83 |
> | 30 | 6.659 ± 1.069      | `ng_recomb_gebv_ocs10_lps2` 6.904 ± 1.161 | +3.68% | 0.80 |
> | 40 | 6.389 ± 0.807      | `ng_meta_portfolio_ocs10_lps2` 6.664 ± 0.916 | +4.31% | 0.72 |
> | 50 | 7.069 ± 1.370      | `ng_meta_portfolio_ocs10_lps2` 7.174 ± 1.195 | +1.47% | 0.93 |
> | 60 | 6.911 ± 0.759      | `ng_recomb_gebv_ocs10_lps2` 6.913 ± 0.826 | +0.03% | 1.00 |
> | 70 | 6.763 ± 0.759      | `var_simple_topn` 6.776 ± 0.695 | +0.19% | 0.98 |
> | 80 | 6.851 ± 0.718      | `ng_meta_portfolio_ocs10_lps2` 6.851 ± 0.718 | +0.00% | 1.00 |
>
> The frontier is competitive (within 4.3% of the best non-frontier
> method) in every band, and tied or within 0.2% at p ≥ 60.
>
> **Top-3 consistency across all 7 parent sizes:**
>
> | Method | Top-3 finishes (of 7) |
> |---|---:|
> | `ng_meta_portfolio_ocs10_lps2` | 5 |
> | `ng_frontier_policy_ocs10_lps2`| 5 |
> | `ng_recomb_gebv_ocs10_lps2`    | 4 |
> | `ng_meta_router_ocs10_lps2`    | 3 |
> | `var_simple_topn`              | 2 |
> | `var_simple_ocs10_lps2`        | 1 |
> | `ng_pmv_blend_balanced_ocs10_lps2` | 1 |
> | `alphamate_opt45`              | **0** |
>
> **Operational defaults on v0.2.0:**
>
> 1. **`ng_frontier_policy_ocs10_lps2` is restored as the recommended
>    default for multi-cycle DH/RIL programs.** It is top-3 in 5/7 bands,
>    rank 4 or better in all 7, and beats real AlphaMate (opt45) in 7/7
>    on the corrected metrics with multi-cycle compounding. The
>    directional pattern is robust (sign-test p = 0.016).
> 2. **`ng_meta_portfolio_ocs10_lps2` is the recommended fallback** —
>    tied with the frontier on top-3 consistency (5/7), simpler
>    dispatch logic, no band-table dependence.
> 3. **`ng_recomb_gebv_ocs10_lps2`** is a strong third choice
>    (4/7 top-3 finishes, identical to frontier at p=60, beats frontier
>    at p=30).
> 4. **`var_simple_topn` is the always-available baseline** — beats
>    AlphaMate (opt45) in 7/7 bands; not the best method but a useful
>    control.
> 5. **`alphamate_opt45` (real AlphaMate) is rank 8/8 in every band**
>    on these defaults — investigate AlphaMate tuning separately before
>    trusting the comparison as evidence about the original AlphaMate
>    algorithm.
>
> **Limits of this re-validation:**
> - **Marker density**: 2500 vs the 5000-marker v0.0.x evidence. The
>   v0.0.x headline used 5K markers; we recovered the *directional*
>   claim at half marker density, which suggests robustness, but
>   absolute effect sizes may differ at 5K.
> - **AlphaMate target degrees**: only `opt45` (middle); `opt30` /
>   `opt60` not tested. The v0.0.x claim was "vs best of opt30/45/60";
>   we tested only the middle. If `alphamate_opt45` is unusually weak,
>   the comparison may understate AlphaMate.
> - **Single AlphaMate config**: `NG_ALPHAMATE_RUNTIME_PATH=
>   C:/Python/Lib/site-packages/torch/lib` — verify the binary
>   produces output identical to AlphaMate documented behavior on a
>   stand-alone test before treating "AlphaMate rank 8/8" as evidence
>   about the algorithm.
> - **Statistical detectability**: 3 reps × 2 cycles gives 6
>   independent realisations per cell. The per-band Welch p-values are
>   ≥ 0.21, so individual-band claims are directional, not significant.
>   The cross-band sign-test (7/7 wins under H₀ = 0.5) reaches p =
>   0.016 — significant.
>
> **v0.2.0 native novelty (no re-validation needed — these are properties
> of the corrected math, not of any benchmark grid):**
> - Closed-form ridge posterior over β with credible intervals on
>   per-cross usefulness via `ng_posterior_cross_predict()`.
> - `P(superior progeny ≥ τ)` as a first-class per-cross metric.
> - Posterior rank stability (`posterior_topn_prob_N`) for OCS.
> - Full off-diagonal posterior PMV `a'Ra + d'(R⊙Σ_β)d` at native package
>   speed; matches genomicMateSelectR-style PMV to machine precision.
> - Posterior over the additive genetic covariance G across traits
>   coupled to the per-trait β posteriors via BCM
>   (`ng_posterior_genetic_covariance()`).
> - End-to-end multi-trait posterior cross prediction with index-level
>   credible intervals (`ng_posterior_multitrait_cross_predict()` — v0.3.0
>   preview).
> - Built-in genetic and phenotypic covariance estimators
>   (`ng_estimate_genetic_covariance()`, `ng_estimate_phenotypic_covariance()`).
>
> These differentiate v0.2.0 from PopVar, SimpleMating, AlphaMate, and
> genomicMateSelectR independently of any benchmark grid: no other tool
> reports posterior credible intervals, threshold-clearing probabilities,
> or rank-stability probabilities at the cross level.

This file records the current evidence-backed state of `nextgen_cross_design`.
It is a guardrail against losing benchmark context in ignored `results/`
artifacts or long-form notes.

## Incorporated Package Assets

The local `master` branch tracks the package source tree needed for a clean
checkout:

- `R/00_utils.R` through `R/22_head_to_head_benchmark.R`
- `R/load.R`
- `src/ng_kernels.cpp`
- package tests under `tests/`
- benchmark runners and summarizers under `tools/`
- package metadata and framework docs

Legacy root-level workspace material is intentionally ignored. This includes
root `R/`, `src/`, `tests/`, `tools/`, `results/`, `external/`, large local
data files, RStudio state, and old standalone scripts. The package code lives
under `nextgen_cross_design/`.

## Evidence-Backed Claims

### Real AlphaMate, Diploid DH/RIL Harness

The strongest current positive result is against the official AlphaMate
executable in the validated diploid DH/RIL AlphaSimR setting.

Evidence source: `BENCHMARK_NOTES.md`, section "Full 5K real-AlphaMate
parent-size validation".

Scope:

- parent sizes 20, 30, 40, 50, 60, 70, and 80
- 3 AlphaSimR replicates
- 3 breeding cycles
- 5K SNPs
- 400 effect-training individuals
- real AlphaMate target degrees `opt30`, `opt45`, and `opt60`

Allowed claim:

- `ng_frontier_policy_ocs10_lps2` beat the best real AlphaMate target across
  all seven parent sizes for average cycle top-10 genetic value in this DH/RIL
  validation grid.
- The average top-10 delta was positive versus both `var_simple_topn` and the
  best AlphaMate target per parent size.

Do not generalize this to every crop, ploidy system, generation scheme, or
trait architecture without a new validation grid.

### PopVar and SimpleMating

PopVar and SimpleMating are incorporated as explicit external-baseline families
and benchmark controls, not just informal references.

Implemented branches include:

- PopVar top-N criteria: `popvar_mu_topn`, `popvar_var_topn`,
  `popvar_uc_topn`, `popvar_musp_topn`
- SimpleMating top-N and constrained criteria: `simple_mpv_topn`,
  `simple_usefa_topn`, `simple_mpv_selectN`, `simple_usefa_selectN`
- adaptive OCS variants using PopVar and SimpleMating scores

Allowed claim:

- The framework is competitive with PopVar and SimpleMating controls and
  contains native formula-compatible recombination usefulness paths that match
  PopVar/SimpleMating-style usefulness in exact all-pair tiers.

Disallowed claim:

- Do not claim the current framework is generally better than PopVar or
  SimpleMating. The balanced external parent-size grid tied many cells, won
  some, and lost some, with the 80-parent shortlist tier favoring
  PopVar/SimpleMating-style comparators.

### Polyploid State

The package now has three separate polyploid-related tracks:

- `ng_poly4x_*`: true autotetraploid 4x model family using AlphaSimR 4x dosage,
  sampled progeny scoring, and policy modes `gain`, `diversity`, and `ocs`
- `ng_poly_model_select()`: model-family router for autotetraploid,
  wheat-like allopolyploid, and complex-polyploid guard cases
- `ng_poly_subgenome_*`: Phase 2A wheat-like disomic subgenome dosage scoring
  and policy selection

The polyploid policy wrappers (`ng_poly4x_ocs`, `ng_poly4x_policy`, `ng_poly_policy`)
now forward the full shared mate-selection control set (strategy dial /
`target_coancestry` / committed matings / group permission+quota / cost+logistics /
`method="evolution"`) to the generic allocator, so polyploid users get the same
controls as diploid (test: `tests/poly4x_strategy_constraints.R`). COMA is integrated
as an optional external arbitrary-ploidy mate-allocation baseline -- see the COMA
Mate-Allocation Study Note above for the outbred+dominance benchmark.

Allowed claim:

- Autotetraploid 4x and wheat-like disomic subgenome support exist as
  test-passing model families.
- The default true-4x grid now evaluates `gain`, `diversity`, and `ocs` policy
  modes end to end for potato/cassava-like autotetraploid scenarios. The
  initial `poly4x_grid` evidence was smoke-level only: one replicate, one
  cycle, parent sizes 20 and 40.
- The 2-replicate, 2-cycle `poly4x_grid_2rep2cycle_20260506` run supports a
  scenario-aware 4x conclusion: policy modes won the potato p20/p40 top10
  cells, while legacy `ng_poly4x_usefulness_topn` and `ng_poly4x_ocs` won the
  cassava p20/p40 top10 cells.
- `ng_poly4x_policy_select(mode = "auto")` now encodes that evidence-scoped
  scenario-aware selector. Explicit `gain`, `diversity`, and `ocs` modes remain
  available and unchanged.

Disallowed claim:

- Do not claim broad superiority for true polyploid breeding programs yet.
  Crop/polyploid AlphaMate smoke tests showed mixed results: wins for some
  diploidized approximations, but failures for cassava tetraploid stress and
  sugarcane polyploid stress. Complex polyploids remain guarded rather than
  solved.
- Do not promote one universal true-4x policy mode from the current evidence.
  The replicated 4x grid is mixed by scenario.

### Multi-Trait Selection

The package now includes a practical multi-trait objective layer for score
tables that already contain cross-level predictions or criteria.

Implemented behavior:

- `ng_multitrait_spec()` records trait columns, maximize/minimize directions,
  optional weights, optional lower/upper thresholds, desired-change inputs, and
  economic weights.
- `ng_add_multitrait_score()` creates `multi_trait_score` from rank-normalized
  trait values and records the resolved weights as diagnostics.
- `method = "auto"` supports breeders who know trait directions but not
  reliable economic weights by using equal normalized weights unless weights are
  supplied.
- `method = "economic_index"` supports breeders who know relative economic
  weights but not exact desired gains. It orients traits to beneficial
  direction, estimates the candidate-table trait covariance, solves stabilized
  covariance-aware economic-index coefficients, and records target and
  predicted response diagnostics.
- `method = "desired_gain"` now implements a covariance-aware economic
  desired-gain index. It orients traits so the beneficial direction is positive,
  estimates the trait covariance in the candidate cross table, solves stabilized
  desired-gain coefficients, and records the coefficients, economic weights,
  target response, and predicted response.
- Thresholds are soft penalties by default and can be made strict with
  `strict_thresholds = TRUE`.
- `ng_optimize_multitrait_mating_plan()` sends the multi-trait score through
  the existing constrained OCS optimizer.
- `ng_run_multitrait_validation()` and
  `tools/run_multitrait_validation.R` provide a deterministic smoke harness
  that compares `auto`, `weighted`, `economic_index`, `desired_gain`, and
  `threshold` against realized multi-trait outcomes in a synthetic trade-off
  scenario.
- `ng_run_multitrait_validation_grid()` and
  `tools/run_multitrait_validation_grid.R` repeat that comparison across
  parent sizes and replicates, then write candidate scores, config, and
  tie-aware winner summaries by realized metric.
- `ng_run_multitrait_crop_validation_grid()` and
  `tools/run_multitrait_crop_validation_grid.R` add the first AlphaSimR-backed
  crop multi-trait validation harness. It uses crop-genome metadata,
  crop-specific genetic covariance profiles, DH base parents, marker-derived
  parent relationships with configurable OCS group-coancestry penalties,
  realized DH progeny family means, scenario-order-stable seeds, score-table
  CSV output, and tie-aware winner summaries.
- `ng_run_head_to_head_benchmark()` and
  `tools/run_head_to_head_benchmark.R` add a formal multi-trait head-to-head
  benchmark harness. It compares NextGen multi-trait OCS candidates against
  PopVar-, SimpleMating-, and AlphaMate-style baselines on the same crop
  scenario, parent set, realized family means, trait directions, and crossing
  budget. The outputs include summary, selections, scores, method registry,
  pairwise candidate-versus-baseline comparisons, config, and tie-aware winner
  CSVs.
- The head-to-head method registry explicitly records `implementation`,
  `exact_external_status`, and `fallback_reason`. Current PopVar,
  SimpleMating, and AlphaMate entries in this multi-trait CI harness are
  style-proxy baselines, not exact external package/executable runs.

Allowed claim:

- Multi-trait positive/negative objective construction, covariance-aware
  economic weighted-index solving, economic desired-gain coefficient solving,
  OCS integration, a deterministic smoke validation harness, and a replicated
  synthetic validation grid are implemented and test-covered.
- The first AlphaSimR-backed crop multi-trait validation grid is implemented
  and test-covered at smoke scale. The 2026-05-06 extensive pass included the
  full R test suite, a reduced-genome 4-scenario x 3-parent-size x 2-replicate
  top-N grid, a metadata-scale 3-scenario grid, and a small OCS allocator check.
- A formal multi-trait head-to-head benchmark contract is implemented and
  test-covered. The first reduced smoke run wrote all expected comparison
  artifacts and preserved the style-proxy baseline caveats in machine-readable
  output fields.
- The first replicated formal multi-trait head-to-head style-proxy grid was run
  on 2026-05-07 across 3 scenarios, 3 parent sizes, and 2 replicates. It wrote
  108 summary rows, 324 selected crosses, 1284 scored candidate pairs, 1188
  candidate-vs-baseline comparison rows, 54 winner rows, 18 config rows, and a
  6-method registry. Results were mixed by metric: the AlphaMate-style weighted
  OCS proxy was the most frequent primary winner, while NextGen desired-gain
  and economic-index OCS were strongest on yield-oriented and diversity-use
  metrics.

Disallowed claim:

- Do not claim evidence-backed superiority for multi-trait selection yet. The
  current multi-trait runners are smoke/synthetic validation harnesses plus
  reduced-scale AlphaSimR crop grids and style-proxy head-to-head grids, not
  replicated production-scale crop-specific exact-external benchmarks.
- Do not describe `popvar_style_*`, `simplemate_style_*`, or
  `alphamate_style_*` rows from `run_head_to_head_benchmark.R` as exact
  PopVar, SimpleMating, or AlphaMate results. Exact external runs must use the
  package/executable-backed benchmark paths and should record package/binary
  availability.

## Current Practical Default

For the validated DH/RIL AlphaSimR evidence, use the frontier/crop-aware policy
layers rather than a single hard-coded method:

- `ng_frontier_policy_ocs*` for the current DH/RIL parent-count frontier
- `ng_crop_aware_policy_ocs*` when crop metadata are available
- `ng_poly4x_policy()` for true autotetraploid 4x experiments
- `ng_poly4x_policy_select(mode = "auto")` when a scenario-aware 4x policy
  dispatch is needed for potato/cassava-like validation runs
- `ng_poly_subgenome_score_crosses()` plus `ng_poly_policy()` for wheat-like
  allopolyploid disomic subgenome experiments
- `ng_multitrait_select_topn()` or `ng_optimize_multitrait_mating_plan()` when
  a breeding objective has multiple traits with mixed positive and negative
  selection directions
- `ng_run_multitrait_validation()` as the first smoke comparison before larger
  multi-trait validation grids
- `ng_run_multitrait_validation_grid()` for replicated synthetic checks across
  parent sizes before investing in crop-specific multi-trait AlphaSimR runs
- `ng_run_multitrait_crop_validation_grid()` for the first crop-specific
  AlphaSimR multi-trait smoke grids with realistic trait covariance profiles
- `ng_run_head_to_head_benchmark()` for formal same-input multi-trait
  candidate-versus-baseline comparisons, while treating style-proxy baseline
  rows as CI contracts rather than exact external-tool evidence

Keep PopVar, SimpleMating, AlphaMate, `var_simple`, recombination-GEBV,
PMV-balanced, meta-selector, and meta-portfolio branches in validation grids.

## Next Scientific Step

The next improvement should separate metric quality from allocation quality and
then validate by crop/ploidy scope:

1. Scale `run_head_to_head_benchmark.R` from reduced smoke settings to
   replicated crop-specific parent-size grids.
2. Add optional exact-external multi-trait jobs where PopVar, SimpleMating, or
   AlphaMate can be run on the same scalarized multi-trait score table.
3. Run the same allocation style over PopVar, SimpleMating, `var_simple`, exact
   recombination usefulness, calibrated PMV, and hybrid scores.
4. Run the same score through OCS, SimpleMating-style constrained selection,
   and top-N allocation.
5. Repeat parent-size validation for 20 through 80 parents.
6. Add a separate true-polyploid validation grid before making claims for
   autotetraploid or complex polyploid breeding programs.

Two further items follow directly from the 2026-08-16 audit, and both remove a stated
Disallowed claim rather than adding a feature for its own sake:

7. **Posterior-ON confidence + `prob_top_tier`** (the F2 item of the cross-priority design).
   Today `cross_confidence` is a within-run normalization of a posterior-OFF mid-parent PEV, so
   it can only RANK crosses inside one plan. A posterior path would make it a calibrated
   interval and support "probability this cross beats the check", retiring the within-run-only
   restriction. The multi-trait case additionally needs the per-draw INDEX value, and the
   block-diagonal index PEV is the approximation it would replace.
8. **Phased polyploid haplotypes -> exact autopolyploid within-family variance.** The current
   locus sum is unbiased over unknown phase but cannot separate two crosses differing only in
   linkage phase (covariance about [-1/3, +1/3] at fixed phase for a duplex x duplex tightly
   linked pair). Supplying phase would make the computation exact and let variance-based metrics
   discriminate autopolyploid crosses as sharply as the diploid path already does. This is the
   polyploid analogue of what `phased_haplotypes` already does for residual-heterozygous RIL
   parents (v0.17.1), so the precedent and the kernel shape both exist.
