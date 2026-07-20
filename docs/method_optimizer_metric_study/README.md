# Method × Optimizer × Metric study — replicate it and apply it

A breeder-practical, recurrent RIL breeding-program simulation that races
`nextgenCrossDesign`'s cross-selection **methods**, **metrics**, and **optimizers**
against **SimpleMating**, driving the package exactly as an end user would
(`ng_run_cross_prediction`).

- **Script:** `tools/run_method_optimizer_metric_study.R`
- **Raw output:** `results/method_optimizer_metric_full_raw.csv` (arm × rep × cycle)
- **Trajectory report:** `trajectories.html` (in this folder — open in a browser)
- **Cycle-25 summary:** `results_cycle25_summary.csv` (in this folder)

There are two ways to use this document:

1. **Replicate / adapt the benchmark** — run the simulation on your own scenario to
   compare methods before you trust one (§4).
2. **Apply the winner to a real program** — take your real genotypes + phenotypes and
   design crosses with the right method and parameters (§5). This is the part most
   breeders want.

---

## 1. How to read an "arm" name

Every competitor in the study is named **`<metric>_<allocation>`** — a *scoring metric*
(how good is one cross?) combined with an *allocation method* (which set of crosses do I
commit to?). Two exceptions: `strategy_balanced` names the diversity **dial** instead of
an allocation, and `random`/`simple_usefa_select` are the controls.

| Arm | Metric (score a cross by…) | Allocation | Optimizer |
|---|---|---|---|
| `mean_ocs` | `mean` — mid-parent breeding value only | OCS | evolution |
| `var_simple_ocs` | `var_simple` — parental relationship *distance* | OCS | evolution |
| `uc_vpm_ocs` | `uc`+`vpm` — usefulness on recombination variance | OCS | evolution |
| `uc_pmv_ocs` | `uc`+`pmv` — usefulness on posterior variance | OCS | evolution |
| `pmv_ocs` | `pmv` — posterior within-family variance | OCS | evolution |
| `vpm_ocs` | `vpm` — recombination within-family variance | OCS | evolution |
| `var_complex_ocs` | `var_complex` — PopVar-style usefulness (default) | OCS | evolution |
| `strategy_balanced` | `var_complex` usefulness | **strategy dial @ balanced** | evolution |
| `var_complex_alphamate` | `var_complex` usefulness | **AlphaMate-style** frontier | evolution |
| `uc_pmv_greedy` | `uc`+`pmv` | OCS | **greedy_local** |
| `uc_pmv_mip` | `uc`+`pmv` | OCS | **MIP** |
| `simple_usefa_select` | SimpleMating usefulness (external) | SimpleMating cull | — |
| `random` | none (control) | random | — |

To change what's raced, edit the **`arm_registry()`** table at the top of
`tools/run_method_optimizer_metric_study.R`. Each row is one call to the `arm()` helper:

```r
arm(name, engine, metric, uc_variance_source, allocation, optimizer, strategy)
# e.g. race var_complex usefulness under OCS + evolution:
arm("var_complex_ocs", "package", "var_complex", "pmv", "ocs", "evolution")
# e.g. the diversity dial on the same metric:
arm("strategy_diversity", "package", "var_complex", "pmv", "ocs", "evolution", strategy = "diversity")
```

`engine` is `"package"` (drive `ng_run_cross_prediction`), `"simplemating"` (external
comparator), or `"random"` (control).

---

## 2. The metrics — how to pick one

Set with **`trait_value_metric`** (and `uc_variance_source` when the metric is `uc`).

| `trait_value_metric` | What it scores | Use when |
|---|---|---|
| `mean` | mid-parent GEBV only (expected progeny mean) | highly polygenic trait, or small / low-h² training — the variance term adds noise |
| `uc` | **usefulness** = mean + i·SD; pick the SD source with `uc_variance_source = "pmv"` \| `"vpm"` \| `"var_simple"` | you want the *best progeny*, not the average one |
| `pmv` | usefulness on **posterior mean variance** (adds marker-effect uncertainty to recombination variance) | default-quality usefulness |
| `vpm` | usefulness on **recombination variance** only | faster usefulness; ranks almost like `pmv` |
| `var_complex` | **PopVar-style** complex usefulness — the package default | a safe default; best for oligogenic traits with adequate training |
| `var_simple` | parental **relationship distance** (a diversity proxy, *not* a merit predictor) | you are managing diversity — see the evidence note below |

`i` is the selection intensity implied by `selection_prop` (default 0.10).

**Evidence (from this study and the package's own studies):**
- **Single round / short term:** a usefulness metric (`var_complex`, `pmv`, `vpm`, `uc`)
  ranks crosses best. They are near-interchangeable. `mean` is competitive only when
  the trait is very polygenic or the training set is weak.
- **Long-term recurrent selection (recycling parents many cycles):** pure usefulness
  metrics **exhaust genetic variance and plateau**. Manage diversity explicitly instead
  — either the **`var_simple`** distance metric, or (better) keep a good merit metric and
  add a **coancestry penalty / strategy dial / target coancestry** (§3). In this study
  `var_simple_ocs` and `strategy_balanced` finished 1st and 2nd at cycle 25; SimpleMating
  and the pure-usefulness cluster finished lower with depleted diversity.

---

## 3. The allocation methods and the strategy dial

**`allocation_method`** decides *which* crosses to commit to, given the scores and how
related the parents are:

| `allocation_method` | What it does |
|---|---|
| `ocs` (default) | Optimal Contribution Selection: maximize gain **subject to controlling group coancestry** (relatedness build-up). Tune with `lambda_group`, `lambda_mating`, `target_coancestry`. |
| `alphamate_style` | native AlphaMate-style **gain–diversity frontier** selection to a target degree. |
| `alphamate_executable` | call the external AlphaMate binary (needs the executable). |

**The strategy dial** sets *where on the gain–diversity frontier to sit*, on top of a
merit metric. Set **`strategy`** to one of:

| `strategy` | Frontier degree | Meaning |
|---|---|---|
| `high_gain` | 15° | push gain, spend diversity |
| `balanced` | 45° | balance gain and diversity |
| `diversity` | 75° | protect diversity, give up some gain |

Or set a continuous **`diversity_emphasis`** (0–100), or a hard **`target_coancestry`**
ceiling. These are the levers a recurrent-selection breeder should reach for.

---

## 4. The optimizers — **what is used here**

The **optimizer** is the *engine that solves the OCS allocation problem* — it does not
change the metric or the objective, only how well/fast the mating plan is found. Set with
**`optimizer`**.

| `optimizer` | What it is | Speed | Quality |
|---|---|---|---|
| `auto` | dispatches to `mip_contribution` when a diversity penalty is active + `lpSolve` is present, else `mip_linear` / `greedy_local`; falls back to greedy on oversized/timed-out MIP | adaptive | best available |
| `mip_contribution` (alias **`mip`**, `ocs`) | exact-style MIP OCS with the coancestry penalty | slowest (size/time-guarded) | best objective |
| `mip_linear` (alias `lp`) | MIP without the coancestry penalty (use when `lambda_group = 0`) | slow | exact for gain-only |
| `greedy_local` (alias `egsi`) | fast greedy + local swap | fastest | slightly below MIP |
| `repair_local` | top-N + capacity repair + local swap | fast | near greedy |
| **`evolution`** (aliases `ga`, `de`, `memetic`) | native **memetic genetic algorithm**; warm-started from greedy with elitism, so never worse than greedy; tune with `evol_solutions` / `evol_iterations` / `evol_stop` | moderate | competitive with MIP |

**What THIS study uses:**
- **`evolution` is the standard optimizer on every allocation arm** (all the `*_ocs`
  arms, `strategy_balanced`, `var_complex_alphamate`). This is the deliberate choice so
  the metric comparison is not confounded by different optimizers.
- The study **also compares optimizers** on a single fixed metric (`uc/pmv`): the
  `uc_pmv_greedy` (greedy_local), `uc_pmv_mip` (MIP), and `uc_pmv_ocs` (evolution) arms.
- `random` and `simple_usefa_select` do not use the package optimizer.

**Result — which optimizer to use:** over the 25-cycle recurrent horizon,
**evolution ≈ MIP, and both clearly beat greedy** on *realized* gain — evolution
7.20 ± 0.30, MIP 7.07 ± 0.44, greedy 6.78 ± 0.30 at cycle 25 for `uc/pmv` (the
evolution–MIP paired difference, 0.13 ± 0.29, is **not** significant). Exact MIP wins the
one-shot objective bake-off, but on realized recurrent gain evolution matches it **and**
scales to large candidate sets where MIP is size/time-guarded — so **`evolution` is the
recommended default** for real/recurrent programs (never worse than greedy; use `auto`/`mip`
when you want the exact objective for a single small mating decision).

---

## 5. Replicate the benchmark

Run from the repository root (`Rscript` on the PATH; AlphaSimR installed):

```bash
# defaults: 10 reps × 25 cycles × 13 arms, evolution optimizer, parallel, C++ kernel
Rscript tools/run_method_optimizer_metric_study.R
```

Everything is configured by `NG_MOM_*` environment variables (defaults in parentheses).
Defaults are a tractable-but-realistic architecture; scale up to match your program.

| Env var | Default | Meaning |
|---|---|---|
| `NG_MOM_REPS` | 10 | independent replicates (run in parallel) |
| `NG_MOM_CYCLES` | 25 | recurrent breeding cycles per rep |
| `NG_MOM_WORKERS` | 0 | parallel workers (0 = cores − 2; forking, macOS/Linux) |
| `NG_MOM_SEED` | 20260703 | base RNG seed (fully reproducible) |
| `NG_MOM_USE_CPP` | 1 | use the C++ recombination-variance kernel |
| `NG_MOM_N_FOUNDER` | 80 | founder lines simulated per rep |
| `NG_MOM_N_PARENTS` | 30 | elite parents selected each cycle |
| `NG_MOM_N_CROSSES` | 40 | crosses made each cycle |
| `NG_MOM_N_CHR` | 5 | chromosomes |
| `NG_MOM_QTL_PER_CHR` | 40 | QTL per chromosome (× chr = total QTL) |
| `NG_MOM_SNP_PER_CHR` | 200 | SNP markers per chromosome |
| `NG_MOM_SEG_SITES` | 500 | segregating sites per chromosome (runMacs) |
| `NG_MOM_FOUNDER_H2` | 0.5 | founder trait heritability |
| `NG_MOM_F2_PER_CROSS` | 80 | F2 family size per cross |
| `NG_MOM_BULK_SIZE` | 5 | bulk size at F3/F4 |
| `NG_MOM_WITHIN_FAMILY_PROP` | 0.10 | within-family selection proportion (F2) |
| `NG_MOM_BETWEEN_FAMILY_PROP` | 0.20 | across-family selection proportion (PYT) |
| `NG_MOM_AYT_PROP` | 0.70 | AYT selection proportion |
| `NG_MOM_STAGE_H2` | 0.03,0.15,0.40,0.60,0.70 | F2..F6 stage heritabilities |
| `NG_MOM_SELECTION_PROP` | 0.10 | selection intensity for the usefulness `i` |
| `NG_MOM_METHOD_VARPMV` | fast | PMV mode (`fast` or `full_posterior`) |
| `NG_MOM_MIN_EFFECT_RELIABILITY` | 0.20 | reliability floor for scoring |
| `NG_MOM_RECOMBINATION_MODEL` | haldane | map function (`haldane` or `kosambi`) |
| `NG_MOM_MAX_USES` | 6 | max crosses per parent |
| `NG_MOM_LAMBDA_GROUP` | 0.05 | OCS group-coancestry penalty |
| `NG_MOM_LAMBDA_MATING` | 0.02 | OCS pairwise-relatedness penalty |
| `NG_MOM_EVOL_SOLUTIONS` | 100 | GA population size |
| `NG_MOM_EVOL_ITERATIONS` | 150 | GA generations |
| `NG_MOM_EVOL_STOP` | 30 | GA early-stop patience |
| `NG_MOM_INCLUDE_SIMPLEMATING` | 1 | include the external SimpleMating arm |
| `NG_MOM_SIMPLEMATING_GITHUB` | 1 | build SimpleMating from GitHub if missing |
| `NG_MOM_OUTPUT_DIR` | results/ | output directory |
| `NG_MOM_OUTPUT_PREFIX` | method_optimizer_metric | output filename prefix |

Example — a larger, wheat-scale run:

```bash
NG_MOM_N_FOUNDER=200 NG_MOM_N_PARENTS=60 NG_MOM_N_CROSSES=100 \
NG_MOM_N_CHR=12 NG_MOM_QTL_PER_CHR=60 NG_MOM_SNP_PER_CHR=500 NG_MOM_SEG_SITES=1644 \
Rscript tools/run_method_optimizer_metric_study.R
```

**Outputs:** `<prefix>_raw.csv` / `.rds` — one row per arm × rep × cycle with the six
breeder metrics: `genetic_gain_index`, `var_index_gv` (additive variance),
`expected_heterozygosity`, `prop_polymorphic_markers`, `mean_maf`,
`mean_parent_relationship`. The console prints the ranked cycle-25 summary.

---

## 6. Apply it to a **real breeding program**

For real crossing decisions you do **not** run the simulation — you call
`ng_run_cross_prediction()` directly on *your* data. The simulation above just tells you
*which settings to trust*.

### 6.1 Prepare four tables (CSV or in-memory data frames)

1. **Genotype** — one row per candidate parent, columns = markers, dosages **0 / 1 / 2**
   (allele count). First column is the line ID. RIL/DH parents should be near-homozygous
   (0/2); residual heterozygous calls are handled but should be rare.
   ```
   NAME, SNP0001, SNP0002, SNP0003, ...
   Line_A,     0,      2,      2, ...
   Line_B,     2,      0,      2, ...
   ```
2. **Phenotype** — one row per parent, the trait(s) you selected them on.
   ```
   NAME, yield
   Line_A, 61.2
   Line_B, 58.9
   ```
3. **Marker map** — marker, chromosome, position (bp or cM).
   ```
   SNP_code, Chromosome, Position_BP
   SNP0001, 1, 1250000
   ```
4. **Trait direction** — which way each trait should go.
   ```
   trait, column, direction
   yield, yield, increase
   ```

### 6.2 The call

```r
source("R/load.R"); ng_load()   # or library(nextgenCrossDesign) if installed

result <- ng_run_cross_prediction(
  # --- your data ---
  phenotype_file = "phenotype.csv",
  genotype_file  = "genotype.csv",
  map_file       = "marker_map.csv",
  direction_file = "trait_direction.csv",
  id_col            = "NAME",
  map_marker_col    = "SNP_code",
  map_chr_col       = "Chromosome",
  map_pos_col       = "Position_BP",
  map_position_unit = "bp",       # positions are base pairs
  bp_per_cm         = 1e6,        # ~1 cM per Mb; set to your genome's ratio

  # --- HOW to score crosses (the metric) ---
  prediction_mode    = "trait_by_trait",
  trait_value_metric = "var_complex",   # see §2; var_simple for long recurrent programs
  uc_variance_source = "pmv",           # only used when trait_value_metric = "usefulness"
  progeny            = "RIL",            # "RIL" or "DH" — must match your progeny type
  selection_prop     = 0.10,            # selection intensity for usefulness

  # --- HOW to allocate crosses (allocation + optimizer + diversity) ---
  allocation_method  = "ocs",
  optimizer          = "evolution",     # recommended (see §4)
  n_crosses          = 100,
  max_uses_per_parent= 6,
  lambda_group       = 0.05,            # diversity management; or:
  # strategy         = "balanced",      # the friendly dial, or
  # target_coancestry= 0.05,            # a hard relatedness ceiling

  # --- outputs ---
  write_outputs = TRUE, output_dir = "out", output_file = "crossing_plan.xlsx"
)

result$selected_crosses   # the recommended crossing plan (parent1, parent2, …)
result$candidate_crosses  # every candidate cross with its scores
result$plan_summary       # gain, group coancestry, diversity of the plan
```

### 6.3 Which settings to choose (decision guide)

| Your situation | `trait_value_metric` | diversity lever | `optimizer` |
|---|---|---|---|
| One round, want the best progeny | `var_complex` (or `uc`/`pmv`/`vpm`) | `lambda_group = 0.05` | `evolution` |
| Highly polygenic trait / weak training | `mean` | `lambda_group = 0.05` | `evolution` |
| **Recurrent program, recycling parents many cycles** | keep `var_complex` **and** add a dial | `strategy = "balanced"` or `target_coancestry` | `evolution` |
| Want maximum diversity retention | `var_simple` | `lambda_group ≥ 0.05` | `evolution` |
| Multiple traits | set `prediction_mode` + `multi_trait_method` (see BACKEND_USER_GUIDE.md) | as above | `evolution` |

**Rules of thumb from the evidence:**
- Use **`evolution`** as your optimizer — it matches exact MIP on realized gain, beats
  greedy, scales to large candidate sets, and is never worse than greedy (warm-started).
- For a **single cycle**, a usefulness metric (`var_complex`/`pmv`/`vpm`/`uc`) is best.
- For a **multi-cycle recurrent program**, do **not** chase pure usefulness — manage
  diversity explicitly (strategy dial / `target_coancestry` / `lambda_group`), or use
  `var_simple`, or you will deplete variance and gain will plateau.
- Match **`progeny`** to your real system (`RIL` vs `DH`) — it changes the within-family
  variance and therefore the usefulness ranking.

---

## 7. Files in this folder

| File | What it is |
|---|---|
| `README.md` | this guide |
| `results_cycle25_summary.csv` | cycle-25 mean per arm (gain, diversity, additive variance) |
| `trajectories.html` | interactive 25-cycle gain & diversity trajectories (open in a browser) |

See also `docs/BACKEND_USER_GUIDE.md` for the full package API and multi-trait selection.
