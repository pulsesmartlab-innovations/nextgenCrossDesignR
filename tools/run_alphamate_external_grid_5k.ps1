$ErrorActionPreference = "Stop"

function Set-DefaultEnv($name, $value) {
  if ([string]::IsNullOrWhiteSpace([Environment]::GetEnvironmentVariable($name))) {
    [Environment]::SetEnvironmentVariable($name, $value, "Process")
  }
}

Set-DefaultEnv "NG_GRID_PREFIX" "alphamate_external_grid_5k"
Set-DefaultEnv "NG_GRID_PARENT_SIZES" "20,30,40,50,60,70,80"
Set-DefaultEnv "NG_GRID_USE_CPP" "1"
Set-DefaultEnv "NG_GRID_REPS" "3"
Set-DefaultEnv "NG_GRID_CYCLES" "3"
Set-DefaultEnv "NG_GRID_N_FOUNDERS" "120"
Set-DefaultEnv "NG_GRID_N_CHR" "5"
Set-DefaultEnv "NG_GRID_GENOME_LENGTH_M" "1"
Set-DefaultEnv "NG_GRID_SEG_SITES" "1200"
Set-DefaultEnv "NG_GRID_SNP_PER_CHR" "1000"
Set-DefaultEnv "NG_GRID_QTL_PER_CHR" "40"
Set-DefaultEnv "NG_GRID_PROGENY_PER_CROSS" "40"
Set-DefaultEnv "NG_GRID_PHENO_REPS" "3"
Set-DefaultEnv "NG_GRID_PHENO_H2" "0.5"
Set-DefaultEnv "NG_GRID_EFFECT_H2_PRIOR" "0.5"
Set-DefaultEnv "NG_GRID_EFFECT_KFOLD" "5"
Set-DefaultEnv "NG_GRID_EFFECT_TRAINING_N" "400"
Set-DefaultEnv "NG_GRID_METHODS" "var_simple_topn,alphamate_opt30,alphamate_opt45,alphamate_opt60,ng_frontier_policy_ocs10_lps2"
Set-DefaultEnv "NG_ALPHAMATE_RUNTIME_PATH" "C:\Python\Lib\site-packages\torch\lib"
Set-DefaultEnv "NG_ALPHAMATE_EVOL_ITERATIONS" "1000"
Set-DefaultEnv "NG_ALPHAMATE_EVOL_STOP" "200"
Set-DefaultEnv "NG_ALPHAMATE_THREADS" "1"
Set-DefaultEnv "NG_ALPHASIMR_THREADS" "1"
Set-DefaultEnv "NG_SHARED_SCORING" "1"

powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_parent_size_grid.ps1
if ($LASTEXITCODE -ne 0) {
  throw "AlphaMate external 5K grid failed with exit code $LASTEXITCODE"
}
