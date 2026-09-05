$ErrorActionPreference = "Stop"

$root = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$script = Join-Path $root "tools\run_poly4x_benchmark.R"

if (-not (Test-Path -LiteralPath $script -PathType Leaf)) {
  throw "Poly4x benchmark runner not found: $script"
}

if (-not (Get-Command "Rscript" -ErrorAction SilentlyContinue)) {
  throw "Rscript is required to run the poly4x grid, but it was not found on PATH."
}

function Get-ProcessEnv {
  param([Parameter(Mandatory = $true)][string]$Name)
  [Environment]::GetEnvironmentVariable($Name, "Process")
}

function Set-ProcessEnv {
  param(
    [Parameter(Mandatory = $true)][string]$Name,
    [AllowNull()][string]$Value
  )
  [Environment]::SetEnvironmentVariable($Name, $Value, "Process")
}

function Split-GridValue {
  param([Parameter(Mandatory = $true)][string]$Value)
  $Value.Split(",") | ForEach-Object { $_.Trim() } | Where-Object { $_ }
}

function Read-PositiveInt {
  param(
    [Parameter(Mandatory = $true)][string]$Name,
    [Parameter(Mandatory = $true)][string]$Value
  )
  $parsed = 0
  if (-not [int]::TryParse($Value, [ref]$parsed) -or $parsed -le 0) {
    throw "$Name must be a positive integer; got '$Value'."
  }
  $parsed
}

$envNames = @(
  "NG_POLY4X_SCENARIO",
  "NG_POLY4X_N_PARENTS",
  "NG_POLY4X_PREFIX",
  "NG_POLY4X_METHODS",
  "NG_POLY4X_POLICY_MODES",
  "NG_POLY4X_REPS",
  "NG_POLY4X_CYCLES",
  "NG_POLY4X_TOP_CROSSES",
  "NG_POLY4X_SCORE_PROGENY_PER_CROSS",
  "NG_POLY4X_REALIZED_PROGENY_PER_CROSS"
)

$savedEnv = @{}
foreach ($name in $envNames) {
  $value = Get-ProcessEnv $name
  $savedEnv[$name] = @{
    Exists = $null -ne $value
    Value = $value
  }
}

try {
  $scenarios = if (Get-ProcessEnv "NG_POLY4X_GRID_SCENARIOS") {
    Get-ProcessEnv "NG_POLY4X_GRID_SCENARIOS"
  } else {
    "potato_autotetraploid_4x,cassava_autotetraploid_4x"
  }
  $parentSizes = if (Get-ProcessEnv "NG_POLY4X_GRID_PARENT_SIZES") {
    Get-ProcessEnv "NG_POLY4X_GRID_PARENT_SIZES"
  } else {
    "20,40"
  }
  $prefix = if (Get-ProcessEnv "NG_POLY4X_GRID_PREFIX") {
    Get-ProcessEnv "NG_POLY4X_GRID_PREFIX"
  } else {
    "poly4x_grid"
  }

  $scenarioList = Split-GridValue $scenarios
  $parentList = Split-GridValue $parentSizes
  if (-not $scenarioList) {
    throw "NG_POLY4X_GRID_SCENARIOS did not contain any scenarios."
  }
  if (-not $parentList) {
    throw "NG_POLY4X_GRID_PARENT_SIZES did not contain any parent sizes."
  }

  $callerTopCrosses = Get-ProcessEnv "NG_POLY4X_TOP_CROSSES"
  $callerRealizedProgeny = Get-ProcessEnv "NG_POLY4X_REALIZED_PROGENY_PER_CROSS"

  foreach ($scenario in $scenarioList) {
    foreach ($parentSizeValue in $parentList) {
      $parentSize = Read-PositiveInt "NG_POLY4X_N_PARENTS" $parentSizeValue

      if (-not (Get-ProcessEnv "NG_POLY4X_METHODS")) {
        Set-ProcessEnv "NG_POLY4X_METHODS" "ng_poly4x_var_topn,ng_poly4x_usefulness_topn,ng_poly4x_ocs"
      }
      if (-not (Get-ProcessEnv "NG_POLY4X_POLICY_MODES")) {
        Set-ProcessEnv "NG_POLY4X_POLICY_MODES" "gain,diversity,ocs"
      }
      if (-not (Get-ProcessEnv "NG_POLY4X_REPS")) {
        Set-ProcessEnv "NG_POLY4X_REPS" "1"
      }
      if (-not (Get-ProcessEnv "NG_POLY4X_CYCLES")) {
        Set-ProcessEnv "NG_POLY4X_CYCLES" "1"
      }
      if (-not (Get-ProcessEnv "NG_POLY4X_SCORE_PROGENY_PER_CROSS")) {
        Set-ProcessEnv "NG_POLY4X_SCORE_PROGENY_PER_CROSS" "6"
      }
      if (-not (Get-ProcessEnv "NG_POLY4X_REALIZED_PROGENY_PER_CROSS")) {
        Set-ProcessEnv "NG_POLY4X_REALIZED_PROGENY_PER_CROSS" "8"
      }

      $realizedProgeny = Read-PositiveInt "NG_POLY4X_REALIZED_PROGENY_PER_CROSS" (Get-ProcessEnv "NG_POLY4X_REALIZED_PROGENY_PER_CROSS")
      $topCrossesValue = if ($callerTopCrosses) { $callerTopCrosses } else { "4" }
      $topCrosses = Read-PositiveInt "NG_POLY4X_TOP_CROSSES" $topCrossesValue
      $minTopCrosses = [int][Math]::Ceiling($parentSize / [double]$realizedProgeny)

      if ($parentSize -gt ($topCrosses * $realizedProgeny)) {
        if ($callerTopCrosses) {
          throw "NG_POLY4X_TOP_CROSSES=$topCrosses with NG_POLY4X_REALIZED_PROGENY_PER_CROSS=$realizedProgeny cannot cover parent size $parentSize. Set NG_POLY4X_TOP_CROSSES to at least $minTopCrosses or increase NG_POLY4X_REALIZED_PROGENY_PER_CROSS."
        }
        $topCrosses = $minTopCrosses
      }

      Set-ProcessEnv "NG_POLY4X_TOP_CROSSES" ([string]$topCrosses)
      if (-not $callerRealizedProgeny) {
        Set-ProcessEnv "NG_POLY4X_REALIZED_PROGENY_PER_CROSS" "8"
      }

      Set-ProcessEnv "NG_POLY4X_SCENARIO" $scenario
      Set-ProcessEnv "NG_POLY4X_N_PARENTS" ([string]$parentSize)
      Set-ProcessEnv "NG_POLY4X_PREFIX" "${prefix}_${scenario}_p${parentSize}"

      Write-Host "Running $scenario parents=$parentSize top_crosses=$topCrosses realized_progeny=$realizedProgeny prefix=$(Get-ProcessEnv "NG_POLY4X_PREFIX")"
      Push-Location $root
      try {
        Rscript --vanilla $script
      } finally {
        Pop-Location
      }
    }
  }
} finally {
  foreach ($name in $envNames) {
    if ($savedEnv[$name].Exists) {
      Set-ProcessEnv $name $savedEnv[$name].Value
    } else {
      Set-ProcessEnv $name $null
    }
  }
}
