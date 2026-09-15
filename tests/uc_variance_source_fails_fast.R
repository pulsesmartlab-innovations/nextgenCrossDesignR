# An impossible metric/variance pairing must be rejected at config time.
#
# `uc_variance_source` passes match.arg with "parent_distance" / "le", and the
# runner then hard-errors on that combination under trait_value_metric =
# "usefulness" -- correctly, because a relationship distance is not a trait
# variance and mu + i*sqrt(distance) is not a quantity. But the error was raised
# during scoring: per trait, after marker effects had been fitted, and possibly
# inside a parallel worker where the message is mangled or lost.
#
# The legality of the pairing does not depend on the data. It is knowable before a
# single genotype is read, so it must be settled before a single genotype is read.
# This is the same lesson as multi_trait_method = "smith_hazel", which sat behind
# the entire posterior pipeline and cost 40.6 hours on a real run.
#
# The observable for "before any fitting" is that ng_cp__build_ctx() -- which
# touches no data -- throws on its own. If it does, no worker can have started.
.libPaths(c(normalizePath(".Rlib", mustWork = FALSE), .libPaths()))
source(file.path("R", "load.R"))
ng_load(use_cpp = FALSE, verbose = FALSE)

set.seed(37)
n <- 20L; m <- 30L
geno <- matrix(sample(c(0L, 2L), n * m, TRUE), nrow = n)
ids <- sprintf("P%02d", seq_len(n))
rownames(geno) <- ids; colnames(geno) <- sprintf("M%02d", seq_len(m))
y <- as.numeric(geno %*% c(rnorm(6), rep(0, m - 6))) + rnorm(n, sd = 0.3)

gdf <- data.frame(NAME = ids, geno, check.names = FALSE, stringsAsFactors = FALSE)
ph  <- data.frame(NAME = ids, YIELD = y, stringsAsFactors = FALSE)
map <- data.frame(SNP = colnames(geno), chr = rep(1:2, each = m / 2),
                  cm = rep(seq(0, 60, length.out = m / 2), 2), stringsAsFactors = FALSE)
dir <- data.frame(Trait = "YIELD", Selection_direction = "increase",
                  stringsAsFactors = FALSE)

cfg <- function(metric, ucvs) {
  list(genotype = gdf, genotype_id_col = "NAME",
       phenotype = ph, phenotype_id_col = "NAME", traits_to_use = "YIELD",
       trait_direction = dir, direction_trait_col = "Trait",
       direction_column_col = "Trait", direction_direction_col = "Selection_direction",
       marker_map = map, map_marker_col = "SNP", map_chr_col = "chr",
       map_pos_cm_col = "cm", map_position_unit = "cM",
       progeny = "DH", parent_type = "inbred", n_crosses = 5L,
       write_outputs = FALSE, write_figures = FALSE, run_posterior_prediction = FALSE,
       trait_value_metric = metric, uc_variance_source = ucvs, seed = 1L)
}

# ---- the impossible pairing is refused, and refused EARLY ------------------
for (bad in c("parent_distance", "le")) {
  e <- tryCatch({ do.call(ng_run_cross_prediction, cfg("usefulness", bad)); NULL },
                error = function(e) conditionMessage(e))
  stopifnot(!is.null(e))
  # The message must say what is wrong, not merely that something is.
  stopifnot(grepl("not a trait variance", e, fixed = TRUE))
}

# ---- legal pairings still run ---------------------------------------------
for (good in c("pmv", "vpm")) {
  r <- do.call(ng_run_cross_prediction, cfg("usefulness", good))
  stopifnot(is.list(r), nrow(r$selected_crosses) > 0L)
  stopifnot(identical(r$settings$uc_variance_source, good))
}

# ---- parent_distance as the METRIC is legal and unaffected ----------------
# It is the one marker-free choice in the package; refusing it here would remove
# the only escape hatch available when marker effects cannot be trusted.
r_pd <- do.call(ng_run_cross_prediction, cfg("parent_distance", "pmv"))
stopifnot(is.list(r_pd), nrow(r_pd$selected_crosses) > 0L)
stopifnot(identical(r_pd$settings$trait_value_metric, "parent_distance"))

# ---- the `le` alias resolves to its canonical name ------------------------
r_le <- do.call(ng_run_cross_prediction, cfg("le", "pmv"))
stopifnot(identical(r_le$settings$trait_value_metric, "parent_distance"))
stopifnot(identical(r_le$settings$trait_value_metric_input, "le"))

cat("uc_variance_source_fails_fast: PASS\n")
