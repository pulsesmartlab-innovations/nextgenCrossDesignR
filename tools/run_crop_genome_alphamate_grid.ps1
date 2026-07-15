$ErrorActionPreference = "Stop"

function Set-DefaultEnv($name, $value) {
  if ([string]::IsNullOrWhiteSpace([Environment]::GetEnvironmentVariable($name))) {
    [Environment]::SetEnvironmentVariable($name, $value, "Process")
  }
}

Set-DefaultEnv "NG_CROP_GRID_PREFIX" "crop_genome_alphamate_grid"
Set-DefaultEnv "NG_CROP_GRID_SCENARIOS" "compact_selfing,maize_like,bread_wheat_hexaploid_approx,barley_like,field_pea_like,potato_tetraploid_stress,cassava_diploid,cassava_tetraploid_stress,sugarcane_polyploid_stress"
Set-DefaultEnv "NG_CROP_GRID_PARENT_SIZES" "20,40"
Set-DefaultEnv "NG_CROP_GRID_REPS" "1"
Set-DefaultEnv "NG_CROP_GRID_CYCLES" "1"
Set-DefaultEnv "NG_CROP_GRID_USE_CPP" "1"
Set-DefaultEnv "NG_CROP_GRID_METHODS" "var_simple_topn,alphamate_opt30,alphamate_opt45,alphamate_opt60,ng_frontier_policy_ocs10_lps2,ng_crop_aware_policy_ocs10_lps2"
Set-DefaultEnv "NG_CROP_GRID_EXTERNAL_SHORTLIST_MULTIPLIER" ""
Set-DefaultEnv "NG_CROP_GRID_ALPHASIMR_THREADS" "1"
Set-DefaultEnv "NG_CROP_GRID_SHARED_SCORING" "1"
Set-DefaultEnv "NG_ALPHAMATE_RUNTIME_PATH" "C:\Python\Lib\site-packages\torch\lib"
Set-DefaultEnv "NG_ALPHAMATE_EVOL_ITERATIONS" "300"
Set-DefaultEnv "NG_ALPHAMATE_EVOL_STOP" "120"
Set-DefaultEnv "NG_ALPHAMATE_THREADS" "1"

Rscript nextgen_cross_design\tools\run_crop_genome_scenarios.R
if ($LASTEXITCODE -ne 0) {
  throw "Crop-genome AlphaMate grid failed with exit code $LASTEXITCODE"
}
