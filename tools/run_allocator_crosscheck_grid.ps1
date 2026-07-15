$ErrorActionPreference = "Stop"

# Diagnostic grid to separate metric quality from mate-allocation quality.
# It compares top-N ranking, OCS allocation, and SimpleMating-style constrained
# selection for the same score families.
$env:NG_EXTERNAL_GRID_PREFIX = if ($env:NG_EXTERNAL_GRID_PREFIX) { $env:NG_EXTERNAL_GRID_PREFIX } else { "allocator_crosscheck_5k" }
$env:NG_EXTERNAL_GRID_PARENT_SIZES = if ($env:NG_EXTERNAL_GRID_PARENT_SIZES) { $env:NG_EXTERNAL_GRID_PARENT_SIZES } else { "50,60,70,80" }
$env:NG_EXTERNAL_GRID_REPS = if ($env:NG_EXTERNAL_GRID_REPS) { $env:NG_EXTERNAL_GRID_REPS } else { "2" }
$env:NG_EXTERNAL_GRID_CYCLES = if ($env:NG_EXTERNAL_GRID_CYCLES) { $env:NG_EXTERNAL_GRID_CYCLES } else { "1" }
$env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER = if ($env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER) { $env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER } else { "2" }
$env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL = if ($env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL) { $env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL } else { "etk_dh_pmv_var_blend_cal,uc_dh_blend,var_simple,mpv" }
$env:NG_EXTERNAL_GRID_USE_CPP = if ($env:NG_EXTERNAL_GRID_USE_CPP) { $env:NG_EXTERNAL_GRID_USE_CPP } else { "0" }
$env:NG_SHARED_SCORING = if ($env:NG_SHARED_SCORING) { $env:NG_SHARED_SCORING } else { "1" }
$env:NG_EXTERNAL_GRID_METHODS = if ($env:NG_EXTERNAL_GRID_METHODS) {
  $env:NG_EXTERNAL_GRID_METHODS
} else {
  "var_simple_topn,var_simple_ocs10_lps1,var_simple_select4,ng_pmv_blend_cal_topn,ng_ocs_mip10_lps1,ng_hybrid_select4,popvar_musp_topn,popvar_musp_ocs10_lps1,popvar_musp_select4,simple_usefa_topn,simple_usefa_ocs10_lps1,simple_usefa_select4"
}

powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_external_parent_size_grid.ps1
