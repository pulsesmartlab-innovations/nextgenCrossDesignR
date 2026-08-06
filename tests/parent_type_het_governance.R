# parent_type governs the residual-heterozygosity audit, independent of the
# progeny `target`. A DH / fully-fixed line is homozygous by construction, so a
# heterozygous locus beyond the QC tolerance is a data error -> BLOCK. RILs
# legitimately retain residual het -> proceed. Legacy `assume_inbred` is
# deprecated and reconciled (TRUE->inbred, FALSE->ril) with a one-time warning.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

## --- unit: the reconciler + policy ---------------------------------------
stopifnot(identical(ng_reconcile_parent_type(), "inbred"))           # default
stopifnot(identical(ng_reconcile_parent_type("dh"), "dh"))
stopifnot(identical(ng_reconcile_parent_type("ril"), "ril"))
stopifnot(ng_parent_type_blocks_het("inbred"), ng_parent_type_blocks_het("dh"))
stopifnot(!ng_parent_type_blocks_het("ril"))
# legacy assume_inbred is deprecated: warns and maps
wmap <- tryCatch({ ng_reconcile_parent_type("inbred", assume_inbred = FALSE); NULL },
                 warning = function(w) conditionMessage(w))
stopifnot(!is.null(wmap), grepl("deprecated", wmap))
stopifnot(identical(suppressWarnings(ng_reconcile_parent_type("inbred", assume_inbred = FALSE)), "ril"))
stopifnot(identical(suppressWarnings(ng_reconcile_parent_type("ril",    assume_inbred = TRUE)),  "inbred"))

## --- integration: end-to-end governance through ng_run_cross_prediction ---
set.seed(20240804L)
n <- 30L; m <- 240L; ids <- sprintf("L%03d", seq_len(n)); snps <- sprintf("S%03d", seq_len(m))
G <- matrix(sample(c(0L, 2L), n * m, replace = TRUE), n, m, dimnames = list(NULL, snps))
Ghet <- G; Ghet[1, 1:8] <- 1L                     # 8/120 = 6.7% het in one line (> 2% QC tol)
mm <- data.frame(SNP_code = snps, Chromosome = rep(1:6, each = 40),
                 Position_BP = rep(1:40, 6) * 5e5, stringsAsFactors = FALSE)
qtl <- sort(sample(m, 20)); bq <- rnorm(length(qtl))
gv <- as.numeric(scale((G[, qtl] - 1) %*% bq))
run <- function(geno, ...) do.call(ng_run_cross_prediction, c(list(
  genotype = data.frame(NAME = ids, geno, check.names = FALSE),
  phenotype = data.frame(NAME = ids, yield = gv + rnorm(n)),
  trait_direction = data.frame(Trait = "yield", Selection_direction = "increase"),
  marker_map = mm, id_col = "NAME", map_position_unit = "bp", bp_per_cm = 1e6,
  n_crosses = 8L, progeny = "DH", seed = 1L), list(...)))

# clean inbred data: default parent_type = "inbred" scores without error
stopifnot(inherits(run(G), "ng_cross_prediction_result"))
# DH / inbred declaration + real het -> hard blocker
err_dh <- tryCatch({ run(Ghet, parent_type = "dh"); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(err_dh), grepl("BLOCKED", err_dh))
err_default <- tryCatch({ run(Ghet); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(err_default), grepl("BLOCKED", err_default))     # default blocks too
# RIL declaration + het -> proceeds (with the biased-kernel warning)
stopifnot(inherits(suppressWarnings(run(Ghet, parent_type = "ril")), "ng_cross_prediction_result"))
# legacy assume_inbred=FALSE proceeds (mapped to ril) with a deprecation warning
stopifnot(inherits(suppressWarnings(run(Ghet, assume_inbred = FALSE)), "ng_cross_prediction_result"))

## --- DH strict floor (0.5%) -----------------------------------------------
# A DH line is 100% homozygous by construction, so 'dh' uses a strict 0.5%
# fraction floor: het ABOVE the floor blocks (contamination / mislabelled RIL),
# but genotyping-noise BELOW the floor passes so real DH panels are not
# false-blocked. Both are below the looser 2% 'inbred' tolerance.
Gabove <- G; Gabove[1, 1:2] <- 1L                  # 2/240 = 0.83% > 0.5% floor
err_dh_above <- tryCatch({ run(Gabove, parent_type = "dh"); NULL }, error = function(e) conditionMessage(e))
stopifnot(!is.null(err_dh_above), grepl("BLOCKED", err_dh_above))       # dh: above floor blocks
Gbelow <- G; Gbelow[1, 1] <- 1L                    # 1/240 = 0.42% < 0.5% floor (genotyping noise)
stopifnot(inherits(run(Gbelow, parent_type = "dh"), "ng_cross_prediction_result"))       # dh: below floor passes
stopifnot(inherits(run(Gabove, parent_type = "inbred"), "ng_cross_prediction_result"))   # inbred: 0.83% < 2%, tolerated
stopifnot(inherits(suppressWarnings(run(Gabove, parent_type = "ril")), "ng_cross_prediction_result"))

cat("parent_type_het_governance: DH 0.5% floor (above blocks, noise passes); inbred tolerates; RIL proceeds; legacy reconciled\n")
