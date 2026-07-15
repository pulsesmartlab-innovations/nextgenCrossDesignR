$ErrorActionPreference = "Stop"

# Larger parent-size screen using fast all-pair scoring followed by exact
# PopVar/SimpleMating rescoring on a shortlist. This is not the same as the
# exact all-pair external benchmark; it is the scalable validation design.
$env:NG_EXTERNAL_GRID_PREFIX = if ($env:NG_EXTERNAL_GRID_PREFIX) { $env:NG_EXTERNAL_GRID_PREFIX } else { "external_parent_grid_shortlist_5k" }
$env:NG_EXTERNAL_GRID_PARENT_SIZES = if ($env:NG_EXTERNAL_GRID_PARENT_SIZES) { $env:NG_EXTERNAL_GRID_PARENT_SIZES } else { "50,60,70,80" }
$env:NG_EXTERNAL_GRID_REPS = if ($env:NG_EXTERNAL_GRID_REPS) { $env:NG_EXTERNAL_GRID_REPS } else { "1" }
$env:NG_EXTERNAL_GRID_CYCLES = if ($env:NG_EXTERNAL_GRID_CYCLES) { $env:NG_EXTERNAL_GRID_CYCLES } else { "1" }
$env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER = if ($env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER) { $env:NG_EXTERNAL_GRID_SHORTLIST_MULTIPLIER } else { "2" }
$env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL = if ($env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL) { $env:NG_EXTERNAL_GRID_SHORTLIST_SCORE_COL } else { "etk_dh_pmv_var_blend_cal,uc_dh_blend,var_simple,mpv" }
$env:NG_EXTERNAL_GRID_USE_CPP = if ($env:NG_EXTERNAL_GRID_USE_CPP) { $env:NG_EXTERNAL_GRID_USE_CPP } else { "0" }

powershell -ExecutionPolicy Bypass -File nextgen_cross_design\tools\run_external_parent_size_grid.ps1
