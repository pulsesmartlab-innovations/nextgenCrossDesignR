# Dead defensive branches are not free: they claim a case that cannot occur, and a
# reader trusts the claim.
#
# ng_cp__stage_rank() carried `uc_variance_source %in% c("parent_distance","le")`
# in its effect_based_x test, and `"le"` in the metric list. Both are unreachable:
# ng_cp__build_ctx() canonicalises "le" -> "parent_distance" and hard-errors on the
# usefulness + parent_distance pair before any stage runs. Reading that branch, one
# would reasonably conclude the combination is SUPPORTED downstream -- the exact
# opposite of what the engine does.
#
# Removing it is only safe while these two invariants hold. This test is what makes
# the removal safe: it fails the moment either invariant is weakened, rather than
# leaving a silently wrong risk annotation behind.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

base <- formals(ng_run_cross_prediction)
cfg0 <- lapply(base, function(x) if (is.call(x) || is.name(x)) eval(x) else x)
cfg0 <- cfg0[!vapply(cfg0, function(x) identical(x, quote(expr = )), logical(1))]
mk <- function(...) { o <- cfg0; for (nm in names(list(...))) o[[nm]] <- list(...)[[nm]]; o }

set.seed(3); n <- 12L; m <- 20L
ids <- sprintf("P%02d", seq_len(n))
g <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
cfg0$genotype <- data.frame(NAME = ids, g, check.names = FALSE)
cfg0$genotype_id_col <- "NAME"
cfg0$phenotype <- data.frame(NAME = ids, YIELD = as.numeric(g %*% rnorm(m)))
cfg0$phenotype_id_col <- "NAME"
cfg0$traits_to_use <- "YIELD"
cfg0$trait_direction <- data.frame(Trait = "YIELD", Selection_direction = "increase")
cfg0$direction_trait_col <- "Trait"; cfg0$direction_column_col <- "Trait"
cfg0$direction_direction_col <- "Selection_direction"
cfg0$n_crosses <- 4L; cfg0$write_outputs <- FALSE; cfg0$write_figures <- FALSE
cfg0$run_posterior_prediction <- FALSE

# ---- 1. "le" never reaches a stage: it is canonicalised at config time -----
c_le <- ng_cp__build_ctx(mk(trait_value_metric = "le", uc_variance_source = "le"))
stopifnot(identical(c_le$trait_value_metric, "parent_distance"))
stopifnot(identical(c_le$uc_variance_source, "parent_distance"))
# the ORIGINAL request is still reported -- canonicalising must not erase it
stopifnot(identical(c_le$uc_variance_source_input, "le"))

# ---- 2. usefulness + parent_distance never reaches a stage: it is refused --
e <- tryCatch(ng_cp__build_ctx(mk(trait_value_metric = "usefulness",
                                  uc_variance_source = "parent_distance")),
              error = function(e) e)
stopifnot(inherits(e, "error"))
stopifnot(grepl("not a trait variance", conditionMessage(e), fixed = TRUE))

# ---- 3. so every stage sees only resolved tokens --------------------------
for (mtr in c("usefulness", "pmv", "vpm", "parent_distance", "var_complex", "mean")) {
  ctx <- ng_cp__build_ctx(mk(trait_value_metric = mtr))
  stopifnot(!identical(ctx$trait_value_metric, "le"))
  stopifnot(!identical(ctx$uc_variance_source, "le"))
}

cat("resolved_tokens_reach_the_stages: PASS\n")
