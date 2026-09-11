# effect_summary must carry the cross-mean basis, so it survives into result.json.
#
# `trait_mean_source` is assembled by the runner and returned at top level, but the
# JSON app builds result.json from an explicit whitelist of fields and does not pick
# it up -- so through the app the basis is invisible again, which is the whole defect
# this line of work exists to close. The app is owned elsewhere and is not edited here.
#
# `effect_summary` IS on that whitelist, and the basis belongs there on the merits:
# effect_summary already reports cv_predictive_r2 and ridge_lambda -- how well the
# marker model did -- and the basis is the decision taken ON that evidence. A reader
# who sees cv_predictive_r2 = -0.002 next to mean_source = "adjusted_pheno" can
# immediately see why the phenotype mean was used.
#
# Asserted against a REAL run rather than a helper, because the value has to survive
# assembly, and a unit test of a formatting function would not have caught the
# original defect either.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(11)
n <- 30L; m <- 40L
geno <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), nrow = n)
ids <- paste0("P", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- paste0("M", seq_len(m))
b <- rep(0, m); b[1:6] <- rnorm(6, sd = 1.5)
y <- as.numeric(scale(geno %*% b)) + rnorm(n, sd = 0.25)

pheno <- data.frame(NAME = ids, YIELD = y, stringsAsFactors = FALSE)
gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
map <- data.frame(SNP_code = colnames(geno),
                  Chromosome = rep(1:4, length.out = m),
                  Position_cM = rep(seq(0, 120, length.out = m / 4), times = 4),
                  stringsAsFactors = FALSE)
dir <- data.frame(Trait = "YIELD", Selection_direction = "increase",
                  stringsAsFactors = FALSE)

run <- function(thresh) {
  ng_run_cross_prediction(
    genotype = gdf, genotype_id_col = "NAME",
    phenotype = pheno, phenotype_id_col = "NAME", traits_to_use = "YIELD",
    trait_direction = dir, direction_trait_col = "Trait",
    direction_column_col = "Trait", direction_direction_col = "Selection_direction",
    marker_map = map, map_marker_col = "SNP_code", map_chr_col = "Chromosome",
    map_pos_cm_col = "Position_cM", map_position_unit = "cM",
    progeny = "RIL", parent_type = "ril", n_crosses = 8L,
    write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
    min_cv_predictive_r2 = thresh, seed = 20260907L)
}

# A bar the markers cannot clear -> phenotype basis, and effect_summary says so.
lo <- run(0.99)
es <- lo$effect_summary
stopifnot(is.data.frame(es), nrow(es) >= 1L)
stopifnot("mean_source" %in% names(es))
stopifnot(identical(as.character(es$mean_source[[1L]]), "adjusted_pheno"))

# The evidence the decision rests on must sit beside it.
stopifnot("cv_predictive_r2" %in% names(es))

# A bar they do clear -> GEBV basis, same column. The field tracks the decision
# rather than reporting a constant.
hi <- run(0.01)
stopifnot(identical(as.character(hi$effect_summary$mean_source[[1L]]), "GEBV"))

message("effect_summary_reports_mean_source passed")
