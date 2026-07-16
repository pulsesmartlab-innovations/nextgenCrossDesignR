$ErrorActionPreference = "Stop"

function Set-DefaultEnv($name, $value) {
  if ([string]::IsNullOrWhiteSpace([Environment]::GetEnvironmentVariable($name))) {
    [Environment]::SetEnvironmentVariable($name, $value, "Process")
  }
}

Set-DefaultEnv "NG_VALIDATION_PHASE" "grid"
Set-DefaultEnv "NG_VALIDATION_USE_CPP" "1"
Set-DefaultEnv "NG_VALIDATION_REPS" "3"
Set-DefaultEnv "NG_VALIDATION_CYCLES" "3"
Set-DefaultEnv "NG_VALIDATION_PARENT_SIZES" "20,30,40,50,60,70,80"
Set-DefaultEnv "NG_VALIDATION_EFFECT_TRAINING_N" "400"
Set-DefaultEnv "NG_GRID_PREFIX" "frontier_policy_parent_grid_5k"
Set-DefaultEnv "NG_GRID_EXTERNAL_SHORTLIST_MULTIPLIER" "20"
Set-DefaultEnv "NG_ALPHASIMR_THREADS" "1"

powershell -ExecutionPolicy Bypass -File tools\run_framework_validation.ps1
if ($LASTEXITCODE -ne 0) {
  throw "frontier policy validation grid failed with exit code $LASTEXITCODE"
}
