\name{nextgenCrossDesign-api}
\alias{nextgenCrossDesign-api}
\alias{ng_evolutionary_mate_allocation}
\alias{ng_with_rng_seed}
\alias{ng_choose_frontier_point}
\alias{ng_select_by_strategy}
\alias{ng_select_by_target_coancestry}
\alias{ng_progeny_inbreeding_histogram}
\alias{ng_apply_marker_management}
\alias{ng_marker_target_spec}
\alias{ng_marker_target_scores}
\alias{ng_lethal_recessive_spec}
\alias{ng_lethal_recessive_cross_risk}
\alias{ng_trait_check_spec}
\alias{ng_align_check_geno}
\alias{ng_check_reference_value}
\alias{ng_attach_check_reference}
\alias{ng_check_tau_bounds}
\alias{ng_attach_joint_check_probability}
\alias{ng_check_records_from_pheno}
\alias{ng_plot_check_panels}
\alias{ng_plot_check_reference_line}
\alias{ng_check_line_value}
\alias{ng_frontier_export_payload}
\alias{ng_write_frontier_json}
\alias{ng_write_progeny_inbreeding_histogram_json}
\alias{ng_run_comparison_payload}
\alias{ng_write_run_comparison_json}
\alias{ng_add_balanced_usefulness_score}
\alias{ng_add_external_baseline_scores}
\alias{ng_add_gms_vpm_scores}
\alias{ng_add_multitrait_score}
\alias{ng_add_popvar_scores}
\alias{ng_add_simplemating_scores}
\alias{ng_polyploid_score_crosses}
\alias{ng_polyploid_design_crosses}
\alias{ng_polyploid_grm}
\alias{ng_polyploid_dominance_grm}
\alias{ng_polyploid_qc}
\alias{ng_polyploid_fit_effects}
\alias{ng_polyploid_predict_value}
\alias{ng_polyploid_score_crosses_dominance}
\alias{ng_polyploid_progeny_moment_table}
\alias{ng_polyploid_allele_freq}
\alias{ng_alphamate_style_select}
\alias{ng_allocate_family_sizes}
\alias{ng_alphamate_default_executable}
\alias{ng_apply_family_variance_calibrators}
\alias{ng_backend_capability_registry}
\alias{ng_breeder_selection_objective}
\alias{ng_compare_selected_families}
\alias{ng_crop_aware_policy_select}
\alias{ng_crop_genome_env}
\alias{ng_crop_genome_scenarios}
\alias{ng_crop_genome_select}
\alias{ng_cross_priority_summary}
\alias{ng_cross_priority_workbook_tables}
\alias{ng_design_crosses}
\alias{ng_detect_putative_duplicates}
\alias{ng_build_cross_trait_covariance}
\alias{ng_cross_trait_within_family_cov}
\alias{ng_estimate_genetic_covariance}
\alias{ng_estimate_phenotypic_covariance}
\alias{ng_posterior_genetic_covariance}
\alias{ng_posterior_multitrait_cross_predict}
\alias{ng_expected_excess_above_threshold}
\alias{ng_external_methods_needed}
\alias{ng_family_metric_registry}
\alias{ng_family_metric_summary}
\alias{ng_family_selection_threshold}
\alias{ng_fit_family_variance_calibrators}
\alias{ng_fit_ridge_effects}
\alias{ng_fit_ridge_effects_posterior}
\alias{ng_frontier_default_policy}
\alias{ng_frontier_policy_from_spec}
\alias{ng_frontier_policy_select}
\alias{ng_head_to_head_compare_summary}
\alias{ng_head_to_head_dashboard_payload}
\alias{ng_head_to_head_default_methods}
\alias{ng_ld_prune_geno}
\alias{ng_ld_prune_markers}
\alias{ng_load}
\alias{ng_load_cpp}
\alias{ng_make_pairs}
\alias{ng_multitrait_crop_validation_scenario}
\alias{ng_multitrait_policy_select}
\alias{ng_multitrait_select_topn}
\alias{ng_multitrait_spec}
\alias{ng_multitrait_validation_default_methods}
\alias{ng_multitrait_validation_scenario}
\alias{ng_optimize_balanced_usefulness}
\alias{ng_optimize_mating_plan}
\alias{ng_optimize_mating_plan_curve}
\alias{ng_plot_diminishing_returns}
\alias{ng_plot_priority_score_vs_kinship}
\alias{ng_plot_putative_duplicates}
\alias{ng_optimize_breeder_selection_plan}
\alias{ng_optimize_multitrait_mating_plan}
\alias{ng_optimize_robust_mating_plan}
\alias{ng_pareto_mate_allocation}
\alias{ng_parent_kinship}
\alias{ng_p_superior_progeny}
\alias{ng_p_superior_progeny_multitrait}
\alias{ng_add_p_superior_progeny}
\alias{ng_add_p_superior_progeny_multitrait}
\alias{ng_posterior_cross_predict}
\alias{ng_sample_ridge_posterior_bcm}
\alias{ng_sample_ridge_posterior_mcmc}
\alias{ng_polyploid_as_dosage_matrix}
\alias{ng_polyploid_model_select}
\alias{ng_polyploid_policy}
\alias{ng_polyploid_policy_modes}
\alias{ng_polyploid_subgenome_as_dosage_list}
\alias{ng_polyploid_subgenome_grm}
\alias{ng_polyploid_subgenome_score_crosses}
\alias{ng_predict_gebv}
\alias{ng_preflight_input_tables}
\alias{ng_prepare_marker_map}
\alias{ng_run_alphamate}
\alias{ng_run_cross_prediction}
\alias{ng_run_head_to_head_benchmark}
\alias{ng_run_multitrait_crop_validation_grid}
\alias{ng_run_multitrait_validation}
\alias{ng_run_multitrait_validation_grid}
\alias{ng_rank_cross_priority}
\alias{ng_score_crosses}
\alias{ng_score_breeder_objective}
\alias{ng_select_alphamate}
\alias{ng_simplemating_build_crosses}
\alias{ng_simplemating_past_thinning}
\alias{ng_simplemating_relate_thinning}
\alias{ng_simplemating_select_crosses_native}
\alias{ng_simplemating_style_select}
\alias{ng_validate_metric_calibration}
\alias{ng_popvar_style_scores}
\alias{ng_write_backend_capability_registry_json}
\alias{ng_write_crossing_plan_workbook}
\alias{ng_write_cross_priority_workbook}
\alias{ng_write_data_preflight_json}
\alias{ng_write_head_to_head_dashboard_json}
\alias{ng_write_head_to_head_visual_report}
\alias{ng_run_stage}
\alias{ng_poly_run_stage}
\title{Public NextGen Cross Design Backend API}
\description{
Public backend functions for cross scoring, mate allocation, multi-trait
objectives, validation, benchmark reporting, external baseline comparison,
data preflight, and scoped polyploid workflows.
}
\value{
Most functions return data frames, lists, or file paths containing cross scores,
mating plans, validation summaries, benchmark reports, or JSON-ready payloads.
}
\details{
This manual page documents the package's public API surface. Internal helpers
remain available within the namespace but are intentionally not exported.

Native external-tool-inspired components are integrated directly in the package
so the backend does not require PopVar, SimpleMating, or AlphaMate to operate.
Use \code{ng_run_cross_prediction(trait_value_metric = "var_complex")} for the
native complex usefulness metric, \code{ng_popvar_style_scores()} for
lower-level PopVar-style compatibility scoring,
\code{ng_simplemating_relate_thinning()}, \code{ng_simplemating_build_crosses()},
\code{ng_simplemating_past_thinning()}, and
\code{ng_simplemating_select_crosses_native()} for SimpleMating-style parent and
cross handling, and \code{ng_alphamate_style_select()} for AlphaMate-style
gain-diversity target allocation. Exact external calls remain optional through
the PopVar, SimpleMating, and AlphaMate integration wrappers, or through
\code{ng_run_cross_prediction(allocation_method = "alphamate_executable")} for
AlphaMate.

\code{ng_run_cross_prediction()} exposes the main trait-value and uncertainty
controls directly: \code{trait_value_metric}, \code{uc_variance_source},
\code{method_varPMV = "fast"} or \code{"full_posterior"}, \code{ril_mode =
"infinite"}, \code{run_posterior_prediction}, \code{posterior_method},
\code{n_iter}, \code{burn_in}, and \code{use_parallel}. When
\code{method_varPMV = "full_posterior"}, the wrapper fits and uses the full
marker-effect covariance in PMV scoring. When
\code{run_posterior_prediction = TRUE}, posterior cross-prediction tables are
returned in \code{result$posterior_predictions}.

\code{method_varPMV = "fast"} is the diagonal marker-effect-posterior PMV and
is retained as a valid approximation; \code{"full_posterior"} additionally
uses off-diagonal marker-effect covariance. A genetic marker map is required
for recombination-aware variance. Physical positions require an explicit
\code{bp_per_cm}; the package does not infer a universal physical-to-genetic
conversion.

Ridge \code{cv_predictive_r2} is out-of-fold phenotype predictive ability, not
breeding-value reliability or accuracy squared. Native fits therefore report
uncalibrated reliability as missing. Formal \code{economic_index} scoring
(Smith-Hazel, \eqn{b = P^{-1} G a}) requires both a phenotypic covariance
matrix P and an additive-genetic covariance matrix G. Formal
\code{desired_gain} scoring (Pesek-Baker, \eqn{b = G^{-1} d}) requires only G;
P never enters that solve and is used solely to standardise the reported
predicted response by the index standard deviation, which is returned as
\code{NA} when P is absent. Candidate-cross score covariance is not used as a
substitute for either matrix under any method.

A caller-supplied P or G carrying dimnames is REORDERED BY NAME to the trait
order and must name the trait set exactly; a missing, extra, misspelled or
one-sided label is an error. A matrix with no dimnames is interpreted
POSITIONALLY, in trait order. Send dimnames whenever column order cannot be
guaranteed -- a JSON matrix round trip drops them.

A caller-supplied P or G is validated and INVALID INPUT IS REFUSED, never
warned about, because an invalid covariance propagates into the index weights,
the cross ranking, the allocation and the joint probabilities alike. Every
tolerance is RELATIVE to the scale of the matrix being checked
(\eqn{10^{-8} \times \max(1, \mathrm{scale})}), because a covariance carries the
square of the trait's unit and its magnitude therefore says nothing about
validity. The checks are: finite entries; symmetry, judged against the largest
element magnitude; strictly positive diagonal variances; every implied
correlation \eqn{|M_{ts}| \le \sqrt{M_{tt} M_{ss}}} inside \eqn{[-1, 1]},
reported as a correlation; and positive semidefiniteness by eigenvalue.

When BOTH P and G are supplied, \eqn{P - G} must itself be positive
semidefinite. Since \eqn{P = G + R}, the residual covariance \eqn{R = P - G} is
a covariance matrix; this is the multivariate generalisation of
\eqn{0 \le h^2 \le 1}, and a violation means some trait contrast is assigned
more genetic than phenotypic variance. The per-trait special case
\eqn{\mathrm{diag}(G) \le \mathrm{diag}(P)} is reported first, naming the trait
and its implied \eqn{h^2}; otherwise the offending eigenvalue and the trait
contrast that realises it are reported. The check runs after by-name alignment,
so trait order cannot affect it. Note that
\code{ng_estimate_genetic_covariance(method = "two_stage_ridge")} frequently
produces a \code{G} that fails this check against
\code{ng_estimate_phenotypic_covariance()}, because its lambda-inversion
diagonal is not a calibrated genetic variance.

The matrix a formal index INVERTS must additionally be positive DEFINITE, not
merely semidefinite: P for \code{economic_index}, G for \code{desired_gain}.
Numerically singular input, or a condition number above \eqn{10^{8}}, is
refused rather than passed to the index ridge, which would return coefficients
determined by the regulariser instead of by the supplied matrix. The threshold
is the point at which the default ridge stops being a rounding correction; mild
ill-conditioning below it is accepted, and raising \code{ridge} does not lift
it. A merely semidefinite matrix is NOT refused for a method that does not
invert it -- a singular G is accepted under \code{economic_index}, where G
enters only as the forward projection \eqn{G a}.

\code{ng_estimate_genetic_covariance(method="auto")} requires the sommer
multivariate REML engine. The explicit \code{two_stage_ridge} option is a
warned diagnostic heuristic. Its \code{Y} argument must contain phenotypes or
BLUEs. BLUPs are already shrunk toward the mean
(\eqn{Var(BLUP) = \sigma^2_g - PEV}) and must be deregressed
(Garrick-Taylor-Dekkers) before use; GEBVs must not be used at all, being a
linear function of the same markers the GRM is built from. A guard warns when
the supplied \code{Y} shows essentially no residual variation against the GRM,
and the returned matrix carries provenance attributes
(\code{method}, \code{formal_variance_component_estimate},
\code{implied_heritability}, \code{shrunk_input_suspected}, and others).
The legacy \code{cross_mean_blend} column is retained for compatibility but is
now exactly the midpoint of the explicitly selected mean source; it is not an
uncalibrated phenotype/GEBV interpolation.
Superior-progeny probabilities assume independent Gaussian progeny draws. In
posterior prediction they use the actual per-draw GEBV family mean and a
within-family genetic variance; a custom ranking score is not treated as a
family mean.

Relationship units are explicit in scored and optimized outputs:
\code{pair_relationship} and \code{group_relationship} are on the VanRaden
additive-relationship scale, while \code{pair_kinship} and
\code{group_coancestry} are one half of those quantities. Expected immediate
progeny inbreeding and \code{target_coancestry} use coancestry units.
For \code{parent_type="ril"}, fully inbred parents use the standard RIL
kernel. Any residual-heterozygous parent requires complete 0/1 phased
haplotypes whose allele sums match the dosage matrix; the package does not
return the known downward-biased inbred approximation.

Use \code{ng_write_crossing_plan_workbook()} to create breeder-facing Excel
workbooks with top parents, recommended crosses, full-sib screening, editable
criteria, and an embedded selection-index chart.

Use \code{ng_ld_prune_markers()} and \code{ng_ld_prune_geno()} to remove
low-MAF and high-LD redundant markers before formal model fitting or cross
design. The pruning workflow has a native C++ backend for large dosage matrices
and an R fallback.

Use \code{ng_preflight_input_tables()} before scoring. Set
\code{putative_duplicate_check = TRUE} to run the package IBS 0-1 screen for
near-identical genotype profiles. Set
\code{putative_duplicate_action = "remove"} to keep one parent per duplicate
cluster and return cleaned input tables in \code{qc$cleaned_tables}; the audit
metadata is returned in \code{qc$cleaning$putative_duplicates}. Use
\code{ng_write_data_preflight_json()} for compact UI or batch logs, and use
\code{ng_detect_putative_duplicates()} plus \code{ng_plot_putative_duplicates()}
when duplicate-pair tables or duplicate heatmaps are needed for QC review.

Use \code{ng_rank_cross_priority()} after selecting a practical crossing list
such as 100 crosses to add \code{priority_rank}, \code{priority_index}, and
breeder-facing \code{priority_tier} values. The default tiers are
\code{highly_priority}, \code{priority}, \code{medium_priority}, and
\code{low_priority}; use \code{ng_cross_priority_summary()} to summarize the
number of selected crosses in each tier.

Use \code{ng_write_cross_priority_workbook()} to create a breeder-facing Excel
workbook from a selected crossing plan. The workbook reports objective evidence
such as priority tier, score, kinship, threshold violation, parent-use QC, and
duplicate QC, while leaving \code{Breeder_Rationale},
\code{Breeder_Notes}, \code{Final_Decision}, and \code{Crossing_Status}
editable because final rationale differs by breeding program.

Use \code{ng_plot_priority_score_vs_kinship()} to show all candidate crosses in
score-versus-relatedness space and overlay the selected crossing plan colored by
priority tier.

A check is a benchmark genotype (a released variety or commercial check) that is
never crossed and contributes one scalar reference value per trait; it is
reference-only and never filters or reorders candidate crosses. Supply
\code{ng_run_cross_prediction(trait_checks = ng_trait_check_spec(...),
check_geno = ..., check_progeny_size = ...)}, with the check's own genotype
matrix kept separate from the candidate parents. \code{check_progeny_size} is
required whenever \code{trait_checks} is supplied and has no default, because
progeny per family scales P(beat check) directly. Which input the check value
is read from -- predicted from \code{check_geno} versus read from
\code{check_pheno}/\code{check_records} -- is decided by the run's own
resolved \code{mean_source} per trait, not by the caller, so the reference is
always on the same scale as the cross means it is compared against. Use
\code{ng_check_records_from_pheno()} to convert a phenotype-shaped
\code{check_pheno} table into the \code{check_records} structure consumed by
\code{ng_check_reference_value()}, and \code{ng_align_check_geno()} to verify
\code{check_geno} carries the same marker set, order, and dosage coding as the
fitted effects. \code{ng_attach_check_reference()} adds per-trait
\code{<trait>_check_value}, \code{<trait>_vs_check} (direction-aware: positive
always means better), \code{<trait>_check_ok}, and
\code{<trait>_p_beat_check} columns plus \code{checks_all_ok}, without
changing \code{nrow(scores)} or its order. \code{ng_attach_joint_check_probability()}
(built on \code{ng_check_tau_bounds()}) adds \code{p_beat_all_checks}, the
probability a progeny beats every check at once. Where any checked trait's
posterior effect variance (PEV) is available, this is a fixed, seeded
150-draw Monte Carlo estimate (typical error up to about 0.03), computed
under a diagonal (independent-across-traits) approximation to the
cross-trait posterior effect covariance, and clipped to the row-wise
minimum of the per-trait marginals -- see \code{trait_check_reference$
p_beat_all_checks_note} on a completed run's return value for this same
caveat attached to the specific run. A completed run reports the
check in a \code{Checks} workbook sheet and, on the score-versus-kinship plot,
as a reference line on whichever axis carries the mean; use
\code{ng_plot_check_panels()} for per-trait small multiples in multi-trait
runs, \code{ng_check_line_value()} to compute where that line falls on a
scored table's axis for a given trait, and \code{ng_plot_check_reference_line()}
to draw it -- as a horizontal or vertical line, via \code{mean_axis = "y"} or
\code{"x"} -- on any plot of your own that carries the trait mean on one axis.

Use \code{ng_run_stage()} to drive the cross-prediction pipeline one stage at a
time (\code{qc}, \code{predict}, \code{index}, \code{allocate}, \code{rank}),
persisting each stage's accumulated context to a run directory so a caller such
as a front end can compute each stage once and resume from any completed stage.
Driving all stages in order is byte-identical to a single
\code{ng_run_cross_prediction()} call because both paths run the same internal
pipeline; QC blockers are returned as a stage status and refuse downstream
stages until resolved.
}
\keyword{methods}
