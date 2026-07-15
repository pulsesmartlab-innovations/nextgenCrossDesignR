# v0.3.4 QC and Priority Workflow

Version 0.3.4 prepares `nextgenCrossDesign` for independent package testing on
the data-rich workflow.

## Main Change

Putative duplicate genotype handling is now embedded in package QC. Users do
not need to run a separate duplicate-detection script before analysis.

Use:

```r
qc <- ng_preflight_input_tables(
  geno = geno,
  phenotype = phenotype,
  trait_spec = trait_spec,
  marker_map = marker_map,
  ploidy = 2,
  putative_duplicate_check = TRUE,
  putative_duplicate_action = "remove",
  duplicate_threshold = 0.995
)

geno <- qc$cleaned_tables$geno
phenotype <- qc$cleaned_tables$phenotype
```

The package keeps one parent per putative duplicate cluster, removes redundant
parents from genotype and phenotype tables, filters candidate pairs and parent
relationship matrices when supplied, and returns an audit trail in:

```r
qc$cleaning$putative_duplicates
```

## Public API

The duplicate-QC and practical crossing workflow are exposed through the package
API:

- `ng_preflight_input_tables()`
- `ng_write_data_preflight_json()`
- `ng_detect_putative_duplicates()`
- `ng_plot_putative_duplicates()`
- `ng_breeder_selection_objective()`
- `ng_score_breeder_objective()`
- `ng_optimize_breeder_selection_plan()`
- `ng_rank_cross_priority()`
- `ng_cross_priority_summary()`
- `ng_plot_priority_score_vs_kinship()`
- `ng_write_cross_priority_workbook()`

The backend capability registry includes
`putative_duplicate_genotypes_removed` under Data QC so frontend and workflow
integrations can discover the feature from `ng_backend_capability_registry()`.

## Reproducible User Script

The user-facing script is:

```text
data_rich/run_user_step_by_step_nextgenCrossDesign.R
```

From the project root:

```powershell
Rscript data_rich/run_user_step_by_step_nextgenCrossDesign.R
```

The script installs the local source tarball:

```text
dist/nextgenCrossDesign_0.3.4.tar.gz
```

Then it loads genotype, phenotype, marker map, and trait-direction files; runs
package QC with duplicate removal; scores all candidate crosses for all traits;
selects a practical 100-cross plan; ranks the plan into priority tiers; writes
figures; and writes the breeder-facing workbook.

## Expected Data-Rich Output

The latest verified run on `data_rich` produced:

- 160 raw parents.
- 2 putative duplicate pairs.
- 2 duplicate parents removed by package QC.
- 158 parents retained for scoring.
- 12,403 candidate crosses scored.
- 100 recommended crosses.
- Priority tiers: 10 highly priority, 25 priority, 35 medium priority, 30 low
  priority.
- 8 figure files.
- `breeder_crossing_plan.xlsx`.

Key QC files:

- `qc/putative_duplicate_pairs.csv`
- `qc/putative_duplicate_removed_parents.csv`
- `qc/putative_duplicate_rows_removed.csv`
- `qc/data_preflight_cleaned_tables.csv`
- `qc/scoring_warnings.txt`

## Verification Commands

Targeted package checks used for this release:

```powershell
Rscript tests\data_preflight.R
Rscript tests\putative_duplicates.R
Rscript tests\backend_capability_registry.R
```

Full `R CMD check` was attempted, but this Windows environment could not
complete it reliably because the check shell could not consistently access
Rtools `g++`; the Rtools Bash path segfaulted before writing a useful check log.
The final tarball was instead verified by direct source install from `dist/`,
targeted package tests, installed API inspection, and the full data-rich
workflow run.
