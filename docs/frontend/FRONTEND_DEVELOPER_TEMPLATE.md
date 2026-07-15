# Frontend Developer Template for nextgenCrossDesign

This document is the handoff brief for building a professional frontend around the `nextgenCrossDesign` R backend. The goal is a breeder-facing workbench that makes advanced cross prediction, mate allocation, validation, and reporting usable without exposing users to raw package internals.

The backend is the source of truth. The frontend should discover capability from `ng_backend_capability_registry()` and should run approved workflows through the existing R runner contract rather than reimplementing quantitative genetics logic in TypeScript.

> **Start here:** `docs/frontend/README.md` is the step-by-step build order and single index over this
> kit. Drive runs through the headless `tools/run_cross_prediction_json.R` (JSON config →
> `ng_run_result.v1`); build the config form from `docs/frontend/contracts/config_schema.json`; code
> result views against `docs/frontend/contracts/example_result.json`. Capabilities added since this
> brief was first written and to surface in the IA: the **marker-effect training set**
> (`training_genotype`/`training_phenotype` — extra effect-only individuals that are never candidate
> parents; shown via `input_match_audit$training_only_count`/`effect_training_n`), `grm_method`
> (vanraden/yang), on-runner marker steering (`marker_target_spec`/`lambda_marker`), and the **exact
> within-family cross-trait covariance** for the multi-trait threshold probability
> (`ng_cross_trait_within_family_cov` → `ng_add_p_superior_progeny_multitrait(cross_trait_cov=)`).

## Product Position

Build a **breeding decision workbench** for planning crosses, comparing allocation strategies, validating methods, and exporting breeder-facing plans.

The software should feel like an operational research tool: clear, calm, dense enough for repeated use, and explicit about evidence limits. Avoid a marketing landing-page style. The first screen should be a working command center with data readiness, workflow status, and the next useful action.

## Research References

- BreedBase positions itself as a browser-based breeding ecosystem covering field layouts, phenotyping, genotyping, genomic selection, and predictions: https://breedbase.org/
- BreedBase navigation shows useful domain groupings: Search, Manage, Analyze, Breeder Tools, field trials, genotyping, traits, kinship, mixed models, and selection index: https://breedbase.org/
- BMS Pro emphasizes a centralized and intuitive system where data flows between tasks, teams, and seasons, with modules for database queries, breeding activities, statistical analysis, genotyping, and field data capture: https://www.bmspro.io/
- BrAPI is the interoperability reference for plant phenotype and genotype data exchange, modular breeding data types, flexible search, field data collection, genotype visualization, and data transfer: https://brapi.org/
- Nielsen Norman Group usability heuristics should guide the interaction model: system status, breeder language, user control, consistency, error prevention, recognition over recall, efficiency, minimalism, error recovery, and focused help: https://www.nngroup.com/articles/ten-usability-heuristics/
- W3C accessibility principles should guide implementation: perceivable information, operable navigation, understandable content, and robust compatibility: https://www.w3.org/WAI/fundamentals/accessibility-principles/
- Scientific workflow systems such as Galaxy are relevant as mental models for saved histories, reusable workflows, datasets, job status, and reproducibility: https://galaxyproject.org/

MateSel was requested as a reference, but I could not locate a reliable indexed MateSel webpage or YouTube result from the available search environment. If a direct MateSel video or manual is available, review it before final UI polishing, especially for mate allocation controls and breeder terminology.

## Target Users

Primary users:

- Plant breeders deciding which crosses to make this cycle.
- Quantitative geneticists validating method choices, trait directions, and allocation settings.
- Breeding program leads reviewing evidence, risk, and crossing-plan exports.

Secondary users:

- Data managers preparing genotype, phenotype, trait, map, and candidate-cross files.
- Developers integrating the R backend into a lab-server deployment.

## UX Principles

1. **Start with readiness.** The UI should show whether the data are usable before exposing scoring or allocation actions.
2. **Use breeder language first.** Say "recommended crosses", "parent use", "coancestry", "trait direction", and "validation status" before raw function names.
3. **Keep backend evidence visible.** Every method card should show whether it is validated, supported, experimental, guarded, exact external, or style proxy.
4. **Prevent high-cost mistakes.** Do not allow a run to start when IDs are mismatched, duplicate crosses exist, trait directions are missing, or complex polyploid workflows are unsupported. Putative duplicate genotype profiles should be handled through backend QC and its cleaned-table output, not by frontend-only logic.
5. **Make expert controls discoverable.** Defaults should be safe, but advanced controls for OCS penalties, K sweeps, thresholds, and external tools must be accessible.
6. **Separate decision and evidence.** Results should show both selected crosses and why those crosses were selected.
7. **Design for repeated work.** Users will revisit data, adjust objectives, rerun workflows, compare outputs, and export plans.

## Information Architecture

Use a persistent workbench shell with these modules:

| Module | Purpose | Backend source |
| --- | --- | --- |
| Command Center | Current project status, active workflow, latest run, next action | Registry + run catalog |
| Data QC | Import files, inspect schemas, run preflight, fix blockers, review duplicate-removal audit | `ng_preflight_input_tables()`, `ng_write_data_preflight_json()` |
| Traits | Build multi-trait objective, directions, weights, thresholds, desired gains | `ng_multitrait_spec()`, `ng_add_multitrait_score()` |
| Scoring | Configure DH/RIL PMV, usefulness, PopVar-style, SimpleMating-style, reliability metrics | `ng_score_crosses()`, `ng_popvar_style_scores()` |
| Allocation | Choose OCS, family-size allocation, K sweep, parent-use and coancestry constraints | `ng_optimize_mating_plan()`, `ng_optimize_mating_plan_curve()`, `ng_allocate_family_sizes()` |
| Validation | Run smoke checks, crop grids, family calibration, benchmark workflows | scripts under `tools/` |
| Polyploid | Guarded 4x and allopolyploid workflows with explicit validation labels | `ng_poly4x_*`, `ng_poly_subgenome_*`, `ng_poly_model_select()` |
| External Tools | PopVar, SimpleMating, AlphaMate exact/proxy readiness | external integration registry |
| Reports | Dashboard JSON, HTML evidence report, Excel crossing plan, logs | report/export scripts |

## Page Requirements

### Command Center

Show:

- Data readiness status.
- Default breeding system and validation scope.
- Active objective method.
- Latest run state and artifact links.
- Recommended next action.
- Evidence caveats: validated DH/RIL, experimental polyploid, exact/proxy external status.

Primary actions:

- Start data preflight.
- Build objective.
- Queue workflow.
- Open latest report.

### Data QC

Required UI:

- File slots for genotype, phenotype, trait spec, marker map, candidate pairs, parent relationship matrix.
- Schema preview with column mapping.
- Issue list grouped by blocker/warning/info.
- "Run preflight" action that calls the backend preflight workflow.
- Putative duplicate genotype panel showing duplicate pairs, kept parents, removed parents, and cleaned-table row counts from `qc$cleaning$putative_duplicates`.
- Fix guidance written in breeder/data-manager language.

Block run launch when:

- Parent/sample IDs do not align.
- Marker IDs are duplicated or missing from map/effects.
- Candidate crosses contain exact or reciprocal duplicates.
- Trait names are missing or duplicated.
- Polyploid dosage values are outside declared range.

When users enable duplicate cleanup, launch downstream scoring only with the
backend-returned `qc$cleaned_tables`; do not reimplement duplicate removal in
the frontend.

### Traits

Required UI:

- Trait table with trait name, source column, direction, weight, optional min/max threshold, threshold strictness, desired change, and economic weight.
- Segmented method selector: Auto, Weighted, Economic Index, Desired Gain, Threshold.
- Inline explanation of soft vs strict thresholds.
- Preview card explaining how the backend will orient maximize/minimize traits.

Do not hide trait direction. A disease score, lodging score, plant height risk, or maturity risk trait must never default silently to "maximize".

### Scoring

Required UI:

- Method family cards for DH/RIL PMV/usefulness, PopVar-style, SimpleMating-style, AlphaMate-style, frontier policy, crop-aware policy.
- Candidate-cross score table with sortable columns.
- Gain/diversity tradeoff visual.
- Reliability and marker-density notes.
- Exact/proxy labels for external baselines.

### Allocation

Required UI:

- Controls for number of crosses K, max crosses per parent, minimum unique parents, group coancestry penalty, parent-use penalty, and allocation method.
- K sweep and diminishing-returns panel using `ng_optimize_mating_plan_curve()`.
- Plan summary: mean gain, total gain, group coancestry, effective population size estimate, unique parents, max parent use.
- Recommended crosses table with parent1, parent2, score, reason, parent-use impact, and family-size recommendation.

### Validation

Required UI:

- Workflow cards for head-to-head benchmark, multi-trait validation, crop validation grid, family calibration, AlphaSimR benchmark, polyploid benchmarks.
- Runtime estimates from registry.
- Queue/run/log status.
- Result summary with pass/fail and evidence scope.

### Polyploid

Required UI:

- Breeding-system selector with status badges: default, supported, stress screen, experimental, guarded.
- 4x and allopolyploid workflow panels.
- Guardrail warning for complex polyploids: do not allow unsupported runs without empirical generator/config.
- Evidence label on every result.


### Reports

Required UI:

- Run list with filters by scenario, crop, method, status, date.
- Open dashboard JSON summary.
- Open HTML visual report.
- Export Excel crossing plan.
- Download logs and manifest.
- Claim-boundary panel copied from validated state: do not overstate multi-trait or polyploid superiority.

## Backend Contracts

The UI should use the backend capability registry as its menu and capability source:

```r
source("R/load.R")
ng_load(use_cpp = TRUE)
ng_write_backend_capability_registry_json("results/backend_capabilities.json")
```

Existing frontend APIs already point in the right direction:

- `frontend/src/server/capabilities.ts`
- `frontend/src/server/r-runner.ts`
- `frontend/worker/r-worker.mjs`
- `frontend/src/lib/types.ts`

When adding a new backend workflow, update both R workflow maps:

- `frontend/src/server/r-runner.ts`
- `frontend/worker/r-worker.mjs`

## Visual Design Direction

Use a restrained professional palette:

- Background: cool gray/off-white.
- Primary action: deep green or blue-green.
- Warning/guarded status: amber.
- Error/blocker: red.
- Experimental status: violet or steel.
- Avoid one-note palettes and decorative gradients.

Layout:

- Persistent left navigation.
- Top status bar with project, data readiness, and active run.
- Main content uses full-width workbench sections, not nested cards inside cards.
- Tables should be dense, sortable, and readable.
- Cards are for repeated items such as methods, workflows, artifacts, and issue summaries.
- Keep border radius at 8px or less.

Use icons from `lucide-react` in the real Next.js frontend for navigation and command buttons. Prefer familiar icons for run, upload, download, warning, settings, file, table, chart, and workflow controls.

## State Model

Recommended frontend state:

```ts
type WorkbenchStage =
  | "data"
  | "traits"
  | "scoring"
  | "allocation"
  | "validation"
  | "reports";

type EvidenceStatus =
  | "validated"
  | "supported"
  | "style_proxy"
  | "exact_external"
  | "experimental"
  | "guarded"
  | "blocked";

interface WorkbenchRunDraft {
  breedingSystem: string;
  traitMethod: string;
  scoringFamily: string;
  allocationMethod: string;
  workflow: string;
  strictThresholds: boolean;
  kRange?: [number, number];
  maxCrossesPerParent?: number;
  lambdaGroup?: number;
  lambdaParentUse?: number;
}
```

## Acceptance Checklist

A frontend implementation is acceptable when:

- Data preflight blockers visibly stop workflow launch.
- Trait directions are required and visible.
- Soft vs strict thresholds are clear.
- Default DH/RIL path is easy to run.
- Polyploid workflows are labeled experimental or guarded.
- The run queue shows queued, running, completed, failed, and cancelled states.
- Reports show both recommended crosses and evidence caveats.
- All critical actions have loading, success, empty, warning, and error states.
- Navigation and forms work with keyboard controls.
- Color is never the only way status is conveyed.
- Tables remain usable on laptop-width screens.
- The UI does not claim broad superiority beyond `VALIDATED_STATE.md`.

## Working Prototype

The companion prototype lives in `frontend-template/`.

Open:

```text
frontend-template/index.html
```

It is dependency-free and is meant as a visual and interaction reference. The production implementation should translate it into the existing Next.js application using the current `frontend/` stack.
