$ErrorActionPreference = "Stop"

function Get-EnvValue($name, $default) {
  $value = [Environment]::GetEnvironmentVariable($name)
  if ([string]::IsNullOrWhiteSpace($value)) { return $default }
  return $value
}

# Exact package baselines are intentionally kept smaller by default.
# PopVar and SimpleMating usefulness become expensive for full 80-parent,
# 5K-marker all-pair diallels, so use this as a correctness/novelty screen.
$env:NG_USE_CPP = Get-EnvValue "NG_EXTERNAL_USE_CPP" "0"
$env:NG_REPS = Get-EnvValue "NG_EXTERNAL_REPS" "1"
$env:NG_CYCLES = Get-EnvValue "NG_EXTERNAL_CYCLES" "1"
$env:NG_N_PARENTS = Get-EnvValue "NG_EXTERNAL_N_PARENTS" "40"
$env:NG_N_FOUNDERS = Get-EnvValue "NG_EXTERNAL_N_FOUNDERS" "80"
$env:NG_N_CHR = Get-EnvValue "NG_EXTERNAL_N_CHR" "5"
$env:NG_SEG_SITES = Get-EnvValue "NG_EXTERNAL_SEG_SITES" "1200"
$env:NG_SNP_PER_CHR = Get-EnvValue "NG_EXTERNAL_SNP_PER_CHR" "1000"
$env:NG_QTL_PER_CHR = Get-EnvValue "NG_EXTERNAL_QTL_PER_CHR" "40"
$env:NG_TOP_CROSSES = Get-EnvValue "NG_EXTERNAL_TOP_CROSSES" "10"
$env:NG_PROGENY_PER_CROSS = Get-EnvValue "NG_EXTERNAL_PROGENY_PER_CROSS" "40"
$env:NG_PHENO_REPS = Get-EnvValue "NG_EXTERNAL_PHENO_REPS" "3"
$env:NG_PHENO_H2 = Get-EnvValue "NG_EXTERNAL_PHENO_H2" "0.5"
$env:NG_EFFECT_TRAINING_N = Get-EnvValue "NG_EXTERNAL_EFFECT_TRAINING_N" $env:NG_N_PARENTS
$env:NG_EFFECT_H2_PRIOR = Get-EnvValue "NG_EXTERNAL_EFFECT_H2_PRIOR" "0.5"
$env:NG_EFFECT_KFOLD = Get-EnvValue "NG_EXTERNAL_EFFECT_KFOLD" "0"
$env:NG_CALIBRATION_POOL = Get-EnvValue "NG_EXTERNAL_CALIBRATION_POOL" "global"
$env:NG_LAMBDA_PARENT_USE_MODE = Get-EnvValue "NG_EXTERNAL_LAMBDA_PARENT_USE_MODE" "adaptive"
$env:NG_LAMBDA_PARENT_USE = Get-EnvValue "NG_EXTERNAL_LAMBDA_PARENT_USE" "2"
$env:NG_OCS_MAX_CROSSES_PER_PARENT = Get-EnvValue "NG_EXTERNAL_OCS_MAX_CROSSES_PER_PARENT" "10"
$env:NG_METHODS = Get-EnvValue "NG_EXTERNAL_METHODS" "var_simple_topn,popvar_musp_topn,popvar_uc_topn,simple_mpv_topn,simple_usefa_topn,simple_usefa_select4,ng_ocs_mip10_lps1,ng_ocs_mip10_lps2"
$env:NG_OUTPUT_PREFIX = Get-EnvValue "NG_EXTERNAL_OUTPUT_PREFIX" "external_40p_5k_all_exact"

Rscript tools\run_alphasimr_benchmark.R
