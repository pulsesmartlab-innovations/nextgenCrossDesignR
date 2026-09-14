# An invalid multi_trait_method must be rejected BEFORE any modelling.
#
# 2026-09-11/12: a 17-trait joint-modelling run failed after 2,437 minutes -- 40.6
# HOURS -- with
#
#   'arg' should be one of "auto", "economic_index", "desired_gain", "weighted", "threshold"
#
# The config had asked for multi_trait_method = "smith_hazel". That is not a
# parameter value; it is the INTERNAL source label R/19 attaches to the
# economic_index solve (b = P^-1 G a), and the package's own comments describe
# that method as "economic_index (Smith-Hazel, ...)". Writing the name the
# documentation uses is the obvious mistake to make.
#
# The cost was not the typo. It was WHERE the typo was caught:
#
#   * with run_posterior_prediction = FALSE, R/19's own check rejects it in ~2s;
#   * with posteriors ON, that path is not reached first, and the value survives
#     until match.arg() inside ng_posterior_multitrait_cross_predict() (R/32),
#     called from R/39 only AFTER every per-trait model, all cross scoring and
#     the full posterior draw set.
#
# So the fast validation existed and was simply bypassed by the expensive
# configuration -- the one where failing late costs the most. A parameter whose
# legal values are known before any data is touched must be checked before any
# data is touched.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

# Synthetic fixture: random genotypes against a random phenotype, so the marker
# model has no signal BY CONSTRUCTION (cv_predictive_r2 is negative or unevaluable).
# The reliability gate is therefore correct to refuse it -- this file tests scoring,
# outputs and plumbing, not marker quality, so it opts out explicitly.
options(ngcd.effect_gate = "off")

set.seed(3)
n <- 24L; m <- 30L
g <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
ids <- paste0("P", seq_len(n))
rownames(g) <- ids; colnames(g) <- paste0("M", seq_len(m))
gdf <- data.frame(NAME = ids, g, check.names = FALSE, stringsAsFactors = FALSE)
ph  <- data.frame(NAME = ids, A = rnorm(n), B = rnorm(n), stringsAsFactors = FALSE)
map <- data.frame(SNP = colnames(g), chr = rep(1:2, length.out = m),
                  cm = rep(seq(0, 60, length.out = m / 2), 2), stringsAsFactors = FALSE)
dir <- data.frame(Trait = c("A", "B"),
                  Selection_direction = c("increase", "decrease"),
                  stringsAsFactors = FALSE)

run <- function(method, posteriors) {
  ng_run_cross_prediction(
    genotype = gdf, genotype_id_col = "NAME",
    phenotype = ph, phenotype_id_col = "NAME", traits_to_use = c("A", "B"),
    trait_direction = dir, direction_trait_col = "Trait",
    direction_column_col = "Trait", direction_direction_col = "Selection_direction",
    marker_map = map, map_marker_col = "SNP", map_chr_col = "chr",
    map_pos_cm_col = "cm", map_position_unit = "cM",
    progeny = "RIL", parent_type = "ril", n_crosses = 5L,
    multi_trait_method = method,
    run_posterior_prediction = posteriors, posterior_method = "closed_form",
    write_outputs = FALSE, write_figures = FALSE, seed = 1L)
}

# THE ASSERTION. Posteriors ON is the configuration that cost 40 hours, so that is
# the one pinned here. Evidence of fail-fast is structural, not a stopwatch: the
# posterior path emits an unmistakable O(M^2) cost warning, so if that warning
# appears the expensive work ran before the value was checked.
warned <- FALSE
msg <- withCallingHandlers(
  tryCatch({ run("smith_hazel", TRUE); NA_character_ },
           error = function(e) conditionMessage(e)),
  warning = function(w) {
    if (grepl("O(M^2)", conditionMessage(w), fixed = TRUE)) warned <<- TRUE
    invokeRestart("muffleWarning")
  })

stopifnot(!is.na(msg))                    # it must error at all
stopifnot(!warned)                        # ...and BEFORE the posterior pipeline ran

# The message must be actionable: name the legal values, and name the one the
# caller almost certainly wanted. "'arg' should be one of ..." does neither.
stopifnot(grepl("multi_trait_method", msg, fixed = TRUE))
stopifnot(grepl("economic_index", msg, fixed = TRUE))
stopifnot(grepl("smith_hazel", msg, fixed = TRUE))

# A legal value must still be accepted on the same path. "auto" is used rather
# than "weighted" because weighted additionally requires per-trait weights, and
# this test is about method VALIDATION, not about weight validation.
ok <- run("auto", FALSE)
stopifnot(is.list(ok), isTRUE(nrow(ok$selected_crosses) > 0L))

message("multitrait_method_fails_fast passed")
