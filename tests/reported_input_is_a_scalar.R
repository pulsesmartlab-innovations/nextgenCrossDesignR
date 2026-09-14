# A field that reports "what you asked for" must report one value.
#
# ng_cp__build_ctx() captures trait_value_metric_input / uc_variance_source_input BEFORE
# match.arg() resolves them. When the caller supplies nothing, the formal's default is
# still its whole choice VECTOR, so the captured "input" was
# c("pmv","vpm","parent_distance","le") -- four values, presented as the one the caller
# chose. Through the JSON bridge that becomes an array where every consumer expects a
# string.
#
# It is the release's own defect shape: a reported field that is not the operative
# value. The run used "pmv"; the field said "pmv, vpm, parent_distance, le".
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

base <- formals(ng_run_cross_prediction)
cfg <- lapply(base, function(x) if (is.call(x) || is.name(x)) eval(x) else x)
cfg <- cfg[!vapply(cfg, function(x) identical(x, quote(expr = )), logical(1))]
set.seed(9); n <- 12L; m <- 16L; ids <- sprintf("P%02d", seq_len(n))
g <- matrix(sample(c(0L, 2L), n * m, TRUE), n, m,
            dimnames = list(ids, sprintf("M%02d", seq_len(m))))
cfg$genotype <- data.frame(NAME = ids, g, check.names = FALSE)
cfg$genotype_id_col <- "NAME"
cfg$phenotype <- data.frame(NAME = ids, YIELD = as.numeric(g %*% rnorm(m)))
cfg$phenotype_id_col <- "NAME"; cfg$traits_to_use <- "YIELD"
cfg$trait_direction <- data.frame(Trait = "YIELD", Selection_direction = "increase")
cfg$direction_trait_col <- "Trait"; cfg$direction_column_col <- "Trait"
cfg$direction_direction_col <- "Selection_direction"
mk <- function(...) { o <- cfg; a <- list(...); for (nm in names(a)) o[[nm]] <- a[[nm]]; o }

# ---- nothing supplied: the field reports the effective default, as ONE value ----
ctx <- ng_cp__build_ctx(mk())
for (f in c("trait_value_metric_input", "uc_variance_source_input")) {
  v <- ctx[[f]]
  if (length(v) != 1L) {
    stop(f, " reports ", length(v), " values (", paste(v, collapse = ", "),
         ") when the caller supplied none; it must report the one value the run used",
         call. = FALSE)
  }
}
stopifnot(identical(ctx$uc_variance_source_input, "pmv"))
stopifnot(identical(ctx$trait_value_metric_input, "usefulness"))

# ---- supplied explicitly: the caller's own word survives verbatim -------------
# This is the whole point of the field and must not be lost to the scalarisation.
ctx2 <- ng_cp__build_ctx(mk(trait_value_metric = "var_complex"))
stopifnot(identical(ctx2$trait_value_metric_input, "var_complex"))
stopifnot(identical(ctx2$trait_value_metric, "usefulness"))     # resolved, not echoed

ctx3 <- ng_cp__build_ctx(mk(trait_value_metric = "le"))
stopifnot(identical(ctx3$trait_value_metric_input, "le"))       # deprecated spelling kept
stopifnot(identical(ctx3$trait_value_metric, "parent_distance"))

ctx4 <- ng_cp__build_ctx(mk(uc_variance_source = "vpm"))
stopifnot(identical(ctx4$uc_variance_source_input, "vpm"))

cat("reported_input_is_a_scalar: PASS\n")
