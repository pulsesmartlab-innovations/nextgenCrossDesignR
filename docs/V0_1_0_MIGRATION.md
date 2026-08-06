# `nextgenCrossDesign` 0.1.0 — Migration & Statistical-Correctness Notes

This release is a coordinated correctness fix-up of the cross-prediction
metric and optimization layer. It is not API-compatible with 0.0.x because
several columns and helper functions encoded statistical bugs that cannot be
fixed in place without breaking the numerical contract.

The headline benchmark claim from `VALIDATED_STATE.md` (`ng_frontier_policy_ocs10_lps2`
beat AlphaMate target-degree controls across 7/7 parent-size bands) is
**invalidated** until re-run; see "Re-validation required" below.

## Correctness bugs fixed

| # | Where | Bug | Fix |
|---|-------|-----|-----|
| **T1.1** | `R/03_metrics.R`, `src/ng_kernels.cpp` | `recomb_model = "kosambi"` was silently ignored — both the R chromosome recursion and the C++ kernel hardcoded the Haldane decay. | Kosambi requests now route to a dense O(M²) path (`ng_dh_recomb_variance_pairs_dense`). The recursion remains Haldane-only by construction (relies on multiplicativity). |
| **T1.2** | `R/03_metrics.R` | `target = "RIL"` was silently treated as DH. | New `ng_progeny_decay(d, model, target)` implements the Haldane–Waddington (1931) equilibrium `R = 2r / (1 + 2r)` and threads `target` through scoring. RIL routes through the dense path. |
| **T1.3** | `R/02_effects.R` | `beta_var <- sigma_e2 / (colSums(X^2) + lambda)` — this is 1 / diag(X′X + λI), **not** diag((X′X + λI)⁻¹). Under LD between markers the formula systematically under-states posterior marker-effect uncertainty (mean 3–4× error on test data). PMV inherits the bias. | Dual-form Woodbury identity: `diag((X′X + λI)⁻¹)_k = (1/λ)·(1 − x_k′ A⁻¹ x_k)` where `A = XX′ + λI`. Reuses the existing `solve(A, X)` and removes the duplicate solve from `df_eff`. |
| **T1.4** | `R/03_metrics.R` | `pmv_scale = median(var_simple) / median(dh_pmv)` rescaled PMV so its median matched the `var_simple` relatedness-distance metric (different units). `gated_var = effect_rel·dh_pmv + (1 − effect_rel)·var_simple` and `hybrid_var = w·dh_pmv_scaled + (1 − w)·var_simple` mixed (genetic value)² with `2pq`-scaled IBS directly. `μ + i·√V` is dimensionally invalid on the resulting columns. | Removed `pmv_scale`, `dh_pmv_scaled_var`, `gated_var`, `hybrid_var`, and every `uc_dh_scaled` / `uc_gated` / `uc_hybrid` variant. `var_simple` is now reported only as a relatedness/diversity column (consumed by the OCS penalty layer), not as a variance. |
| **T1.5** | `R/03_metrics.R`, `R/00_utils.R` | The DH/RIL variance formula assumes both parents are inbred (dosage ∈ {0, 2}). Heterozygous-parent input produced silently-wrong PMV. | New `ng_audit_inbred_dosage()` detects residual heterozygosity. `ng_score_crosses()` gains `assume_inbred = TRUE` (default) — errors out on residual-het parents, since the DH/RIL `a'Ra` kernel omits their parental gametic (`p(1−p)`) segregation and biases vpm/pmv low. Pass `assume_inbred = FALSE` to downgrade to a warning and proceed with the approximate kernel. (No exact analytical DH/RIL variance for residual-het parents exists in the package; `ng_exact_gms_additive_var()` computes the gametic-MS / F1-offspring variance, a different quantity — zero for inbred parents — not the DH/RIL cross variance.) |
| **T1.6** | `R/19_multi_trait_selection.R` | `economic_index` and `desired_gain` solved `b = Sigma_hat⁻¹ a` where `Sigma_hat` was the **sample Pearson covariance of candidate cross predictions** — not Smith–Hazel `P`, not Pesek–Baker `G`. Also applied a ridge penalty twice (line 227 inflated the diagonal, line 228 added a second ridge). | New `ng_multitrait_index_covariance()` accepts user-supplied `phenotypic_covariance` and `genetic_covariance`. Routes correctly: `b = P⁻¹ G a` (Smith–Hazel), `b = G⁻¹ a` (Pesek–Baker), `b = Sigma_hat⁻¹ a` (legacy proxy, with `cov_source = "phenotypic_proxy"` recorded). Single ridge. Predicted response = `G b` when `G` is supplied. |

## Rigor improvements (no behavior change for legacy users)

| # | Area | Change |
|---|------|--------|
| **T2.1** | `dh_recomb_rel_var = dh_var * recomb_rel_scale` | Removed. The ad-hoc `recomb_rel_scale = min(4, 1 / max(rel, 0.20))` had no statistical basis; effect-uncertainty inflation is the job of PMV (now trustworthy after T1.3). |
| **T2.2** | Reliability gate | The `cross_mean_blend = effect_rel·MPV + (1 − effect_rel)·adjusted` continues to provide a smooth alternative to the hard `min_effect_reliability` cliff used by `cross_mean`. |
| **T2.3** | `ng_selection_intensity(p, n_progeny = NULL)` | New optional `n_progeny`. When supplied, returns the deterministic Blom-approximation expected mean of the top `k = round(n·p)` order statistics of N(0, 1), replacing the infinite-N `φ(z)/p` formula. For DH/RIL family sizes of 50–200 this corrects a 5–20 % overestimate. |
| **T2.4** | `ng_validate_metric_calibration()` | Added Fisher-z 95 % CI for Pearson r, Spearman ρ, and HC3 heteroscedasticity-consistent slope SE. RMSE now uses fitted residuals, not raw differences. |
| **T2.5** | `ng_fit_variance_calibrator()` | Default `weighting = "inv_pred_sq"` (WLS with weights ∝ 1 / pred²) reflects the χ²-noise model on realized within-family variance. Negative-slope fallback replaced with a non-negative-slope (NNLS-through-origin) re-fit instead of the `median(y) / median(x)` heuristic. |
| **T2.7** | `ng_add_multitrait_score(threshold_penalty_autoscale = TRUE)` | New default. Threshold penalty is auto-scaled to the IQR of the unpenalized weighted score so a 1-SD violation costs ~1 × IQR of the index. Disable with `threshold_penalty_autoscale = FALSE` for legacy behavior. |
| **T2.8** | `ng_alphamate_style_select()` | Plan summary now includes `alphamate_achieved_degree`, `alphamate_degree_gap`, and `style_proxy = TRUE`. The proxy scalarizes the gain–coancestry frontier rather than solving the target-degree constraint exactly; downstream reports must not equate this with real AlphaMate. |
| **T2.9** | `ng_simplemating_style_select()` | New canonical name; `ng_simplemating_select_crosses_native()` retained as a deprecated alias. Plan summary now carries `style_proxy = TRUE` and a `style_proxy_note`. |
| **T2.10** | `ng_optimize_balanced_usefulness()` | Now emits a `warning()` at every `min_unique_parents` relaxation step and records `balanced_min_unique_relaxation_attempts` in the plan summary, so a silently weakened constraint is visible to the caller. |

## Renamed / removed columns

Removed (no replacement; do not use):
- `dh_pmv_scaled_var`, `dh_pmv_scale`
- `dh_recomb_rel_var`, `dh_recomb_rel_scale`
- `gated_var`, `hybrid_var`
- `uc_dh_scaled`, `uc_gated`, `uc_hybrid` (and all `_gebv` / `_adj` / `_blend` variants)
- `uc_recomb_rel` (and variants)

Retained (now the canonical set):
- Variances (genetic-value² units): `dh_recomb_var` (VPM), `dh_pmv_var` (PMV with diagonal posterior β).
- Relatedness (distance units): `var_simple`, `pair_kinship` — use as OCS diversity penalty, NOT in `μ + i·√V`.
- Means: `cross_mean`, `cross_mean_gebv`, `cross_mean_adjusted_pheno`, `cross_mean_blend`.
- Usefulness: `uc_recomb`, `uc_dh` (× `_gebv`, `_adj`, `_blend`).
- Calibrated: `dh_recomb_var_cal`, `dh_pmv_var_cal`, `var_simple_cal`, and `etk_*_cal` derived from them.
- Default `rank_score` in `ng_score_crosses()` output: `uc_dh_gebv` (or `uc_dh_blend` below `min_effect_reliability`).

Renamed:
- Frontier-policy method `ng_pmv_scaled_blend_balanced_ocs10_lps2` → `ng_pmv_blend_balanced_ocs10_lps2`.
- `ng_simplemating_select_crosses_native()` → `ng_simplemating_style_select()` (old name deprecated, still works).

## Defaults that changed

- `ng_optimize_mating_plan(gain_col = ...)`: was `"uc_gated"`, now `"uc_dh_gebv"`.
- `ng_pareto_mate_allocation(gain_col = ...)`: was `"uc_gated"`, now `"uc_dh_gebv"`.
- `ng_design_crosses()` interior optimizer uses `"uc_dh_gebv"`.
- `ng_external_shortlist_indices(score_col = ...)`: was `"etk_hybrid_var_cal"`, now `"etk_dh_pmv_var_cal"`.

## Re-validation required

`VALIDATED_STATE.md` claims that depend on the changed metrics are invalidated
until re-run on the corrected code. In particular:

- The "ng_frontier_policy_ocs10_lps2 beat best AlphaMate target in 7/7 parent-size bands"
  claim was built on `etk_dh_pmv_scaled_var_blend_cal`, which no longer exists.
- All `ng_pmv_scaled_*` / `ng_recomb_rel_*` / `ng_hybrid_*` methods in
  `tools/run_alphasimr_benchmark.R` are stale and should be either deleted or
  redirected to the in-units equivalents before re-running.
- Per-family `realized_var` should be reported with `n_progeny` so the WLS
  calibrator weighting (T2.5) can be honored.

Recommended re-validation sequence:
1. `Rscript tests/ridge_beta_var.R` (sanity check the dual-form variance).
2. `Rscript tests/recomb_kosambi.R` + `Rscript tests/ril_variance.R`.
3. `Rscript tools/check_r_package.R` (full R CMD check).
4. `powershell -File tools/run_framework_validation.ps1` with `NG_VALIDATION_PHASE=family`
   to re-calibrate the variance metrics on a clean seed.
5. `powershell -File tools/run_framework_validation.ps1` with `NG_VALIDATION_PHASE=grid`
   to re-run the parent-size frontier and decide whether the prior headline
   claim still holds.

## Why these changes are scientifically necessary

Cross-prediction reliability ≠ "method X looks like it ranks well in one run."
A method shipped to breeders worldwide is asserted on the dimensional
correctness of its inputs and the validity of its underlying probabilistic
model. Under-stated marker-effect variance (T1.3) leads to under-stated PMV,
which makes high-`dh_pmv_var` crosses look more favorable than they actually
are; mixing relatedness distance with genetic-value variance (T1.4) makes
`μ + i·σ` not a usefulness in any honest sense; ignoring user-requested map
functions (T1.1) silently changes the answer for crops where Kosambi is
preferred (long-chromosome species); treating RIL as DH (T1.2) understates
recombination in self-pollinated programs. None of these would be visible in
a single AlphaSimR diagnostic run, which is why this audit was needed.
