$ErrorActionPreference = "Stop"

function Set-DefaultEnv($name, $value) {
  $current = [Environment]::GetEnvironmentVariable($name)
  if ([string]::IsNullOrWhiteSpace($current)) {
    Set-Item -Path "Env:$name" -Value $value
  }
}

Set-DefaultEnv "NG_GRID_PREFIX" "mean_source_parent_grid_5k_5rep_train400"
Set-DefaultEnv "NG_GRID_PARENT_SIZES" "20,30,40,50,60,70,80"
Set-DefaultEnv "NG_GRID_REPS" "5"
Set-DefaultEnv "NG_GRID_CYCLES" "1"
Set-DefaultEnv "NG_GRID_N_FOUNDERS" "120"
Set-DefaultEnv "NG_GRID_N_CHR" "5"
Set-DefaultEnv "NG_GRID_SEG_SITES" "1200"
Set-DefaultEnv "NG_GRID_SNP_PER_CHR" "1000"
Set-DefaultEnv "NG_GRID_QTL_PER_CHR" "40"
Set-DefaultEnv "NG_GRID_TOP_CROSSES" "10"
Set-DefaultEnv "NG_GRID_PROGENY_PER_CROSS" "40"
Set-DefaultEnv "NG_GRID_PHENO_REPS" "3"
Set-DefaultEnv "NG_GRID_PHENO_H2" "0.5"
Set-DefaultEnv "NG_GRID_EFFECT_TRAINING_N" "400"
Set-DefaultEnv "NG_GRID_EFFECT_H2_PRIOR" "0.5"
Set-DefaultEnv "NG_GRID_EFFECT_KFOLD" "5"
Set-DefaultEnv "NG_GRID_USE_CPP" "0"
Set-DefaultEnv "NG_GRID_SHARED_SCORING" "1"
Set-DefaultEnv "NG_GRID_CALIBRATION_POOL" "global"
Set-DefaultEnv "NG_GRID_CALIBRATION_MIN_N" "20"
Set-DefaultEnv "NG_GRID_EXTERNAL_SHORTLIST_MULTIPLIER" "5"
Set-DefaultEnv "NG_GRID_EXTERNAL_SHORTLIST_SCORE_COL" "etk_dh_pmv_scaled_var_blend_cal,etk_vpm_blend_cal,etk_dh_pmv_scaled_var_gebv_cal,etk_vpm_gebv_cal,etk_dh_pmv_scaled_var_adj_cal,etk_vpm_adj_cal,mid_parent_value,parent_distance"
Set-DefaultEnv "NG_GRID_METHODS" "var_simple_topn,popvar_uc_topn,simple_usefa_topn,ng_recomb_gebv_cal_topn,ng_recomb_adj_cal_topn,ng_recomb_blend_cal_topn,ng_pmv_gebv_cal_topn,ng_pmv_adj_cal_topn,ng_pmv_blend_cal_topn,popvar_uc_ocs10_lps1,simple_usefa_ocs10_lps1,ng_recomb_gebv_ocs10_lps1,ng_recomb_adj_ocs10_lps1,ng_recomb_blend_ocs10_lps1,ng_pmv_gebv_ocs10_lps1,ng_pmv_adj_ocs10_lps1,ng_pmv_blend_ocs10_lps1"

$repoRoot = Resolve-Path (Join-Path $PSScriptRoot "..\..")
Push-Location $repoRoot
try {
  & (Join-Path $PSScriptRoot "run_parent_size_grid.ps1")
} finally {
  Pop-Location
}
