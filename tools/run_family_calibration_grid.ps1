$ErrorActionPreference = "Stop"

# Family-level calibration benchmark. This diagnoses whether each score predicts
# realized family mean, variance, top-tail value, and maximum value before any
# mate-allocation optimization is applied.
$env:NG_OUTPUT_PREFIX = if ($env:NG_OUTPUT_PREFIX) { $env:NG_OUTPUT_PREFIX } else { "family_calibration_5k" }
$env:NG_CAL_PARENT_SIZES = if ($env:NG_CAL_PARENT_SIZES) { $env:NG_CAL_PARENT_SIZES } else { "80" }
$env:NG_CAL_REPS = if ($env:NG_CAL_REPS) { $env:NG_CAL_REPS } else { "5" }
$env:NG_CAL_N_FAMILIES = if ($env:NG_CAL_N_FAMILIES) { $env:NG_CAL_N_FAMILIES } else { "80" }
$env:NG_CAL_RANDOM_FAMILIES = if ($env:NG_CAL_RANDOM_FAMILIES) { $env:NG_CAL_RANDOM_FAMILIES } else { "40" }
$env:NG_CAL_TOP_PER_METRIC = if ($env:NG_CAL_TOP_PER_METRIC) { $env:NG_CAL_TOP_PER_METRIC } else { "5" }
$env:NG_CAL_BOTTOM_PER_METRIC = if ($env:NG_CAL_BOTTOM_PER_METRIC) { $env:NG_CAL_BOTTOM_PER_METRIC } else { "2" }
$env:NG_PROGENY_PER_CROSS = if ($env:NG_PROGENY_PER_CROSS) { $env:NG_PROGENY_PER_CROSS } else { "80" }
$env:NG_PHENO_REPS = if ($env:NG_PHENO_REPS) { $env:NG_PHENO_REPS } else { "3" }
$env:NG_CAL_INCLUDE_EXTERNAL = if ($env:NG_CAL_INCLUDE_EXTERNAL) { $env:NG_CAL_INCLUDE_EXTERNAL } else { "1" }
$env:NG_CAL_INCLUDE_GMS = if ($env:NG_CAL_INCLUDE_GMS) { $env:NG_CAL_INCLUDE_GMS } else { "1" }
$env:NG_CAL_GMS_MAX_MARKERS = if ($env:NG_CAL_GMS_MAX_MARKERS) { $env:NG_CAL_GMS_MAX_MARKERS } else { "5000" }
$env:NG_USE_CPP = if ($env:NG_USE_CPP) { $env:NG_USE_CPP } else { "0" }
$env:NG_ALPHASIMR_THREADS = if ($env:NG_ALPHASIMR_THREADS) { $env:NG_ALPHASIMR_THREADS } else { "1" }

Rscript nextgen_cross_design\tools\run_family_calibration_benchmark.R
