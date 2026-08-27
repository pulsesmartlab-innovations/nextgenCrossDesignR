$ErrorActionPreference = "Stop"

function Get-EnvValue($name, $default) {
  $value = [Environment]::GetEnvironmentVariable($name)
  if ([string]::IsNullOrWhiteSpace($value)) { return $default }
  return $value
}

function Set-DefaultEnv($name, $value) {
  if ([string]::IsNullOrWhiteSpace([Environment]::GetEnvironmentVariable($name))) {
    [Environment]::SetEnvironmentVariable($name, $value, "Process")
  }
}

function Invoke-Rscript($script) {
  Write-Host "Running: Rscript $script"
  Rscript $script
  if ($LASTEXITCODE -ne 0) {
    throw "Rscript $script failed with exit code $LASTEXITCODE"
  }
}

function Invoke-PowershellScript($script) {
  Write-Host "Running: powershell -ExecutionPolicy Bypass -File $script"
  powershell -ExecutionPolicy Bypass -File $script
  if ($LASTEXITCODE -ne 0) {
    throw "$script failed with exit code $LASTEXITCODE"
  }
}

$phase = (Get-EnvValue "NG_VALIDATION_PHASE" "smoke").ToLowerInvariant()
$parentSizes = Get-EnvValue "NG_VALIDATION_PARENT_SIZES" "20,30,40,50,60,70,80"
$reps = Get-EnvValue "NG_VALIDATION_REPS" "3"
$cycles = Get-EnvValue "NG_VALIDATION_CYCLES" "3"
$effectTrainingN = Get-EnvValue "NG_VALIDATION_EFFECT_TRAINING_N" "400"
$useCpp = Get-EnvValue "NG_VALIDATION_USE_CPP" "1"

Set-DefaultEnv "NG_ALPHASIMR_THREADS" "1"
Set-DefaultEnv "NG_SHARED_SCORING" "1"

Write-Host "nextgen_cross_design validation phase: $phase"

if ($phase -in @("smoke", "all")) {
  Invoke-Rscript "nextgen_cross_design\tests\smoke_test.R"
  Invoke-Rscript "nextgen_cross_design\tests\cpp_consistency.R"
}

if ($phase -in @("family", "all")) {
  Set-DefaultEnv "NG_OUTPUT_PREFIX" "diagnostic_family_5k"
  Set-DefaultEnv "NG_CAL_PARENT_SIZES" $parentSizes
  Set-DefaultEnv "NG_CAL_REPS" $reps
  Set-DefaultEnv "NG_EFFECT_TRAINING_N" $effectTrainingN
  Set-DefaultEnv "NG_USE_CPP" $useCpp
  Set-DefaultEnv "NG_CAL_INCLUDE_EXTERNAL" "1"
  Set-DefaultEnv "NG_CAL_INCLUDE_GMS" "1"
  Invoke-PowershellScript "nextgen_cross_design\tools\run_family_calibration_grid.ps1"
}

if ($phase -in @("allocator", "all")) {
  Set-DefaultEnv "NG_EXTERNAL_GRID_PREFIX" "diagnostic_allocator_5k"
  Set-DefaultEnv "NG_EXTERNAL_GRID_PARENT_SIZES" $parentSizes
  Set-DefaultEnv "NG_EXTERNAL_GRID_REPS" $reps
  Set-DefaultEnv "NG_EXTERNAL_GRID_CYCLES" "1"
  Set-DefaultEnv "NG_EXTERNAL_GRID_EFFECT_TRAINING_N" $effectTrainingN
  Set-DefaultEnv "NG_EXTERNAL_GRID_USE_CPP" $useCpp
  Set-DefaultEnv "NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER" "2"
  Invoke-PowershellScript "nextgen_cross_design\tools\run_allocator_crosscheck_grid.ps1"
}

if ($phase -in @("grid", "all")) {
  Set-DefaultEnv "NG_GRID_PREFIX" "diagnostic_parent_grid_5k"
  Set-DefaultEnv "NG_GRID_PARENT_SIZES" $parentSizes
  Set-DefaultEnv "NG_GRID_REPS" $reps
  Set-DefaultEnv "NG_GRID_CYCLES" $cycles
  Set-DefaultEnv "NG_GRID_EFFECT_TRAINING_N" $effectTrainingN
  Set-DefaultEnv "NG_GRID_USE_CPP" $useCpp
  Set-DefaultEnv "NG_GRID_EXTERNAL_SHORTLIST_MULTIPLIER" "20"
  Set-DefaultEnv "NG_GRID_METHODS" "var_simple_topn,var_simple_ocs10_lps2,popvar_uc_ocs10_lps1,simple_usefa_ocs10_lps1,ng_recomb_gebv_ocs10_lps2,ng_pmv_blend_balanced_ocs10_lps2,ng_meta_portfolio_ocs10_lps2,ng_meta_selector_ocs10_lps2,ng_meta_router_ocs10_lps2,ng_frontier_policy_ocs10_lps2"
  Invoke-PowershellScript "nextgen_cross_design\tools\run_parent_size_grid.ps1"
}

if (!($phase -in @("smoke", "family", "allocator", "grid", "all"))) {
  throw "Unknown NG_VALIDATION_PHASE '$phase'. Use smoke, family, allocator, grid, or all."
}

Write-Host "Validation phase completed: $phase"
