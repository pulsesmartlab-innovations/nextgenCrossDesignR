$ErrorActionPreference = "Stop"

function Get-EnvValue($name, $default) {
  $value = [Environment]::GetEnvironmentVariable($name)
  if ([string]::IsNullOrWhiteSpace($value)) { return $default }
  return $value
}

$gridPrefix = Get-EnvValue "NG_GRID_PREFIX" "nextgen_parent_grid"
$parentSizes = (Get-EnvValue "NG_GRID_PARENT_SIZES" "20,30,40,50,60,70,80").Split(",") |
  ForEach-Object { [int]$_.Trim() }
$trainingMode = (Get-EnvValue "NG_GRID_EFFECT_TRAINING_MODE" "same").ToLowerInvariant()
$fixedEffectTrainingN = Get-EnvValue "NG_GRID_EFFECT_TRAINING_N" ""
$fixedTopCrosses = Get-EnvValue "NG_GRID_TOP_CROSSES" ""

$env:NG_USE_CPP = Get-EnvValue "NG_GRID_USE_CPP" "1"
$env:NG_ALPHASIMR_THREADS = Get-EnvValue "NG_GRID_ALPHASIMR_THREADS" "1"
$env:NG_REPS = Get-EnvValue "NG_GRID_REPS" "2"
$env:NG_CYCLES = Get-EnvValue "NG_GRID_CYCLES" "3"
$env:NG_N_FOUNDERS = Get-EnvValue "NG_GRID_N_FOUNDERS" "120"
$env:NG_N_CHR = Get-EnvValue "NG_GRID_N_CHR" "5"
$env:NG_GENOME_LENGTH_M = Get-EnvValue "NG_GRID_GENOME_LENGTH_M" "1"
$env:NG_SEG_SITES = Get-EnvValue "NG_GRID_SEG_SITES" "1200"
$env:NG_SNP_PER_CHR = Get-EnvValue "NG_GRID_SNP_PER_CHR" "1000"
$env:NG_QTL_PER_CHR = Get-EnvValue "NG_GRID_QTL_PER_CHR" "40"
$env:NG_PROGENY_PER_CROSS = Get-EnvValue "NG_GRID_PROGENY_PER_CROSS" "40"
$env:NG_PHENO_REPS = Get-EnvValue "NG_GRID_PHENO_REPS" "3"
$env:NG_PHENO_H2 = Get-EnvValue "NG_GRID_PHENO_H2" "0.5"
$env:NG_EFFECT_H2_PRIOR = Get-EnvValue "NG_GRID_EFFECT_H2_PRIOR" "0.5"
$env:NG_EFFECT_KFOLD = Get-EnvValue "NG_GRID_EFFECT_KFOLD" "5"
$env:NG_LAMBDA_GROUP = Get-EnvValue "NG_GRID_LAMBDA_GROUP" "1"
$env:NG_LAMBDA_MATING = Get-EnvValue "NG_GRID_LAMBDA_MATING" "0"
$env:NG_LAMBDA_PARENT_USE = Get-EnvValue "NG_GRID_LAMBDA_PARENT_USE" "2"
$env:NG_LAMBDA_PARENT_USE_MODE = Get-EnvValue "NG_GRID_LAMBDA_PARENT_USE_MODE" "adaptive"
$env:NG_OCS_ITER = Get-EnvValue "NG_GRID_OCS_ITER" "5"
$env:NG_LOCAL_ITER = Get-EnvValue "NG_GRID_LOCAL_ITER" "1000"
$env:NG_CALIBRATION_MIN_N = Get-EnvValue "NG_GRID_CALIBRATION_MIN_N" "20"
$env:NG_CALIBRATION_POOL = Get-EnvValue "NG_GRID_CALIBRATION_POOL" "global"
$env:NG_TOPK_PROP = Get-EnvValue "NG_GRID_TOPK_PROP" "0.10"
$env:NG_FAMILY_ALLOCATION_METHOD = Get-EnvValue "NG_GRID_FAMILY_ALLOCATION_METHOD" "marginal_topk"
$env:NG_SHARED_SCORING = Get-EnvValue "NG_GRID_SHARED_SCORING" "1"
$env:NG_EXTERNAL_SHORTLIST_N = Get-EnvValue "NG_GRID_EXTERNAL_SHORTLIST_N" ""
$env:NG_EXTERNAL_SHORTLIST_MULTIPLIER = Get-EnvValue "NG_GRID_EXTERNAL_SHORTLIST_MULTIPLIER" ""
$env:NG_EXTERNAL_SHORTLIST_SCORE_COL = Get-EnvValue "NG_GRID_EXTERNAL_SHORTLIST_SCORE_COL" "etk_dh_pmv_var_blend_cal,uc_dh_blend,var_simple,mpv"
$env:NG_SIMPLEMATING_MIN_CROSS = Get-EnvValue "NG_GRID_SIMPLEMATING_MIN_CROSS" "1"
$env:NG_SIMPLEMATING_MAX_SEARCH = Get-EnvValue "NG_GRID_SIMPLEMATING_MAX_SEARCH" "100000"
$env:NG_SIMPLEMATING_CULLING_K = Get-EnvValue "NG_GRID_SIMPLEMATING_CULLING_K" ""
$env:NG_SIMPLEMATING_THREADS = Get-EnvValue "NG_GRID_SIMPLEMATING_THREADS" "1"
$env:NG_METHODS = Get-EnvValue "NG_GRID_METHODS" "var_simple_topn,var_simple_ocs10_lps2,ng_recomb_gebv_ocs10_lps2,ng_pmv_blend_balanced_ocs10_lps2,ng_meta_portfolio_ocs10_lps2,ng_meta_selector_ocs10_lps2,ng_meta_router_ocs10_lps2,ng_frontier_policy_ocs10_lps2"

foreach ($nParents in $parentSizes) {
  if ([string]::IsNullOrWhiteSpace($fixedTopCrosses)) {
    $topCrosses = [Math]::Max(5, [Math]::Min(20, [int][Math]::Round($nParents / 4.0)))
  } else {
    $topCrosses = [int]$fixedTopCrosses
  }
  if ([string]::IsNullOrWhiteSpace($fixedEffectTrainingN)) {
    $effectTrainingN = switch ($trainingMode) {
      "min80" { [Math]::Max($nParents, 80) }
      "min160" { [Math]::Max($nParents, 160) }
      default { $nParents }
    }
  } else {
    $effectTrainingN = [Math]::Max($nParents, [int]$fixedEffectTrainingN)
  }
  $env:NG_N_PARENTS = [string]$nParents
  $env:NG_TOP_CROSSES = [string]$topCrosses
  $env:NG_EFFECT_TRAINING_N = [string]$effectTrainingN
  $env:NG_OCS_MAX_CROSSES_PER_PARENT = [string][Math]::Min(10, [Math]::Max(2, $topCrosses))
  $env:NG_OUTPUT_PREFIX = "${gridPrefix}_${nParents}p"
  $env:NG_TOTAL_PROGENY = [string]($topCrosses * [int]$env:NG_PROGENY_PER_CROSS)
  $env:NG_FAMILY_SELECTED_TOP_N = [string]$nParents

  Write-Host "Running parent-size scenario: parents=$nParents top_crosses=$topCrosses effect_training_n=$effectTrainingN"
  Rscript tools\run_alphasimr_benchmark.R
  if ($LASTEXITCODE -ne 0) {
    throw "AlphaSimR benchmark failed for parents=$nParents with exit code $LASTEXITCODE"
  }
}

$env:NG_GRID_PREFIX = $gridPrefix
$env:NG_GRID_PARENT_SIZES = ($parentSizes -join ",")
Rscript tools\summarize_parent_size_grid.R
if ($LASTEXITCODE -ne 0) {
  throw "Parent-size grid summary failed with exit code $LASTEXITCODE"
}
