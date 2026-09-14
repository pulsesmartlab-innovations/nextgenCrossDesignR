# The gate fires before scoring, and says what to do about it.
#
# cv_predictive_r2 is known immediately after the ridge fit -- which is cheap --
# while the dense O(M^2) recombination kernel, the posterior draws and the
# cross-trait covariance all run afterwards and unconditionally. So a breeder whose
# markers cannot support genomic means pays the full cost to obtain a plan whose
# cross means carry no genomic content. That happened in production: a 17-trait
# study where every trait fell back to the phenotypic mid-parent, discovered only
# when a breeder queried the delivered spreadsheet.
#
# The gate refuses that run instead, before the expensive work, and names the fix.
#
# It applies exactly to metrics whose value contains a marker-derived within-family
# variance. parent_distance (GRM only) and mean (no variance term) are exempt and
# must still run on the identical data -- they are the designed escape hatches, and
# a gate that blocked them would leave a breeder with weak markers no way forward
# at all.
#
# The fallback tier must WARN rather than refuse, and the warning must be raised in
# the PARENT process: run_trait_job executes under mclapply/parLapply, and a
# warning() from a forked child or PSOCK worker never reaches the caller's handlers,
# so it would never appear in the JSON envelope. An advisory nobody receives is the
# defect this work exists to remove, so the test captures the warning rather than
# trusting that one was issued.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

mk <- function(n, m, n_causal, noise_mult, seed = 41) {
  set.seed(seed)
  g <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
  ids <- sprintf("P%02d", seq_len(n))
  rownames(g) <- ids; colnames(g) <- sprintf("M%02d", seq_len(m))
  b <- c(rnorm(n_causal, sd = 1.5), rep(0, m - n_causal))
  gv <- as.numeric(g %*% b)
  y <- gv + rnorm(n, sd = noise_mult * stats::sd(gv))
  list(
    geno = data.frame(NAME = ids, g, check.names = FALSE, stringsAsFactors = FALSE),
    pheno = data.frame(NAME = ids, YIELD = y, stringsAsFactors = FALSE),
    map = data.frame(SNP = colnames(g), chr = rep(1:2, each = m / 2),
                     cm = rep(seq(0, 70, length.out = m / 2), 2), stringsAsFactors = FALSE))
}

run_with <- function(fx, metric = "usefulness", extra = list()) {
  args <- c(list(
    genotype = fx$geno, genotype_id_col = "NAME",
    phenotype = fx$pheno, phenotype_id_col = "NAME", traits_to_use = "YIELD",
    trait_direction = data.frame(Trait = "YIELD", Selection_direction = "increase",
                                 stringsAsFactors = FALSE),
    direction_trait_col = "Trait", direction_column_col = "Trait",
    direction_direction_col = "Selection_direction",
    marker_map = fx$map, map_marker_col = "SNP", map_chr_col = "chr",
    map_pos_cm_col = "cm", map_position_unit = "cM",
    progeny = "DH", parent_type = "inbred", n_crosses = 6L,
    write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
    trait_value_metric = metric, seed = 20260907L), extra)
  do.call(ng_run_cross_prediction, args)
}

# ---- fixtures: one refusing, one falling back ------------------------------
# NOTE the fixtures are calibrated against the RUN seed (20260907), not the fit's
# default: the CV fold partition is drawn from the run seed, so the same data gives
# a different cv_predictive_r2 under a different seed. Calibrating against the wrong
# seed is how a "negative" fixture turns out to be +0.05.
bad  <- mk(30L, 100L, 80L, 8.0)   # cv_predictive_r2 ~ -0.12 under the run seed
weak <- mk(40L, 60L, 10L, 1.2)    # 0 < cv_predictive_r2 < 0.35

# ---- 1. a variance metric on unusable markers is REFUSED -------------------
e <- tryCatch({ run_with(bad, "usefulness"); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(e))
stopifnot(grepl("YIELD", e, fixed = TRUE))          # names the offending trait
stopifnot(grepl("cv_predictive_r2", e, fixed = TRUE))  # and the evidence
stopifnot(grepl("parent_distance", e, fixed = TRUE))   # and a way forward

# every variance metric, same verdict
for (mt in c("pmv", "vpm")) {
  stopifnot(inherits(tryCatch(run_with(bad, mt), error = function(e) e), "error"))
}
# including the deprecated alias, which resolves to usefulness + pmv
stopifnot(inherits(tryCatch(run_with(bad, "var_complex"), error = function(e) e), "error"))

# ---- 2. the exemptions still run on the SAME data --------------------------
for (mt in c("parent_distance", "mean")) {
  r <- run_with(bad, mt)
  stopifnot(is.list(r), nrow(r$selected_crosses) > 0L)
}

# ---- 3. the fallback tier warns in the PARENT and completes ----------------
w <- character(0)
r_weak <- withCallingHandlers(
  run_with(weak, "usefulness", list(min_cv_predictive_r2 = 0.99)),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })

stopifnot(is.list(r_weak), nrow(r_weak$selected_crosses) > 0L)
stopifnot(identical(r_weak$effect_summary$mean_source[[1L]], "adjusted_pheno"))
# The warning must have actually reached the caller -- not merely "been issued"
# somewhere inside a worker where nobody can hear it.
stopifnot(any(grepl("YIELD", w, fixed = TRUE)))
stopifnot(any(grepl("phenotypic mid-parent", w, fixed = TRUE)))
stopifnot(any(grepl("training", w, fixed = TRUE)))

# ---- 4. strong markers: no gate, no noise ---------------------------------
good <- mk(40L, 60L, 10L, 0.3)
w2 <- character(0)
r_good <- withCallingHandlers(run_with(good, "usefulness", list(min_cv_predictive_r2 = 0.35)),
  warning = function(x) { w2 <<- c(w2, conditionMessage(x)); invokeRestart("muffleWarning") })
stopifnot(identical(r_good$effect_summary$mean_source[[1L]], "GEBV"))
stopifnot(!any(grepl("phenotypic mid-parent", w2, fixed = TRUE)))

# ---- 5. the escape hatch for tooling --------------------------------------
# Benchmark harnesses legitimately drive weak-marker scenarios. They must be able
# to opt out explicitly -- but only explicitly.
r_off <- run_with(bad, "usefulness", list(effect_gate = "off"))
stopifnot(is.list(r_off), nrow(r_off$selected_crosses) > 0L)

# ---- 6. the verdicts are in the artifact ----------------------------------
stopifnot("advisories" %in% names(r_weak))
stopifnot(any(r_weak$advisories$id == "mean_source_phenotype_fallback"))
stopifnot(all(c("id", "severity", "stage", "trait", "message") %in% names(r_weak$advisories)))

cat("reliability_gate_runner: PASS\n")
