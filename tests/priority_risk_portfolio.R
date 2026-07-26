ng_test_use_cpp <- FALSE
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- ng_midparent_pev: quadratic form, centered, per-cross variation ---
set.seed(1)
m <- 12L; ids <- sprintf("P%d", 1:5)
geno <- matrix(rbinom(length(ids) * m, 2, 0.5), length(ids), m, dimnames = list(ids, NULL))
B <- crossprod(matrix(rnorm(m * m), m)) / m         # a valid m x m covariance
pairs <- data.frame(parent1 = c("P1","P1","P3"), parent2 = c("P2","P3","P4"),
                    stringsAsFactors = FALSE)
pev <- ng_midparent_pev(geno, pairs, B)
# reference: 1/4 (xc1+xc2)' B (xc1+xc2) with xc = geno - colMeans
mm <- colMeans(geno); Xc <- sweep(geno, 2, mm, "-")
ref <- vapply(seq_len(nrow(pairs)), function(k) {
  s <- Xc[pairs$parent1[k], ] + Xc[pairs$parent2[k], ]
  0.25 * as.numeric(t(s) %*% B %*% s)
}, numeric(1))
stopifnot(isTRUE(all.equal(pev, ref)))
stopifnot(all(pev >= 0), diff(range(pev)) > 0)        # PSD form, varies by cross
stopifnot(is.na(ng_midparent_pev(geno, pairs, NULL)[1]))  # NULL Sigma -> NA
cat("ng_midparent_pev test passed\n")

# --- ng_cross_confidence: direction, bins, method suffix, NA path ---
pev2 <- c(0.01, 0.04, 0.09, 0.16, 0.25)          # spread = 0.1..0.5
cc <- ng_cross_confidence(pev2, effect_based_x = TRUE)
stopifnot(cc$confidence_method == "midparent_pev_partial")
# higher PEV (spread) => lower confidence
stopifnot(cc$cross_confidence[1] > cc$cross_confidence[5])
stopifnot(as.character(cc$risk_bin[1]) == "low", as.character(cc$risk_bin[5]) == "high")
stopifnot(is.ordered(cc$risk_bin))
cc2 <- ng_cross_confidence(pev2, effect_based_x = FALSE)
stopifnot(cc2$confidence_method == "midparent_pev")
# NULL/degenerate -> reliability, NA bins
cn <- ng_cross_confidence(rep(NA_real_, 4), effect_based_x = TRUE)
stopifnot(cn$confidence_method == "reliability", all(is.na(cn$risk_bin)))
# Degenerate-spread: distinct pev that all clamp to 0 -> reliability fallback
cd <- ng_cross_confidence(c(-1, -2, -3), effect_based_x = TRUE)
stopifnot(cd$confidence_method == "reliability", all(is.na(cd$cross_confidence)), all(is.na(cd$risk_bin)))
cat("ng_cross_confidence test passed\n")

# --- portfolio profile: median cuts, 4 quadrants ---
lvl <- c(10, 10, 1, 1); ups <- c(5, 0.1, 5, 0.1)     # (hi,hi)(hi,lo)(lo,hi)(lo,lo)
pf <- ng_cross_portfolio_profile(lvl, ups)
stopifnot(as.character(pf) == c("breakthrough","workhorse","long_shot","deprioritize"))

# --- annotator: cross_upside = sqrt(vpm), all columns present, level echoed ---
cr <- data.frame(parent1 = c("P1","P3"), parent2 = c("P2","P4"), stringsAsFactors = FALSE)
ann <- ng_annotate_cross_priority(cr, level = c(3, 1), vpm = c(4, 9),
                                  pev = c(0.02, 0.08), effect_based_x = TRUE)
stopifnot(isTRUE(all.equal(ann$cross_upside, c(2, 3))))     # sqrt(vpm)
stopifnot(isTRUE(all.equal(ann$cross_level, c(3, 1))))
stopifnot(all(c("cross_confidence","risk_bin","confidence_method","portfolio_profile")
              %in% names(ann)))
stopifnot(ann$confidence_method[1] == "midparent_pev_partial")
cat("ng_annotate_cross_priority test passed\n")

# --- summary: counts reconcile with rows ---
df <- data.frame(
  priority_tier = factor(c("a","a","b","b"), levels = c("a","b")),
  portfolio_profile = factor(c("breakthrough","workhorse","breakthrough","deprioritize"),
                             levels = c("breakthrough","workhorse","long_shot","deprioritize")),
  cross_level = c(4,3,2,1), cross_upside = c(2,1,2,1), cross_confidence = c(.9,.8,.7,.6),
  stringsAsFactors = FALSE)
sm <- ng_cross_portfolio_summary(df)
stopifnot(sum(sm$n) == nrow(df))
stopifnot(all(c("priority_tier","portfolio_profile","n","mean_level","mean_upside",
                "mean_confidence") %in% names(sm)))

# NA portfolio_profile must form its own "(unclassified)" group, not be dropped
df_na <- data.frame(
  priority_tier = c("a","a","b"),
  portfolio_profile = c("breakthrough", NA, "workhorse"),
  cross_level = c(4,3,2), cross_upside = c(2,1,2), cross_confidence = c(.9,.8,.7),
  stringsAsFactors = FALSE)
sm_na <- ng_cross_portfolio_summary(df_na)
stopifnot(sum(sm_na$n) == nrow(df_na))
stopifnot("(unclassified)" %in% sm_na$portfolio_profile)

# 0-row input -> 0-row data.frame with the six expected columns, not NULL
sm_empty <- ng_cross_portfolio_summary(df[0, , drop = FALSE])
stopifnot(is.data.frame(sm_empty), nrow(sm_empty) == 0L)
stopifnot(identical(names(sm_empty),
                    c("priority_tier","portfolio_profile","n","mean_level",
                      "mean_upside","mean_confidence")))

# All-NA cross_level within a group -> mean_level is NA, not NaN
df_allna <- data.frame(
  priority_tier = c("a","a"),
  portfolio_profile = c("breakthrough","breakthrough"),
  cross_level = c(NA_real_, NA_real_), cross_upside = c(1,2), cross_confidence = c(.5,.6),
  stringsAsFactors = FALSE)
sm_allna <- ng_cross_portfolio_summary(df_allna)
stopifnot(!any(is.nan(sm_allna$mean_level)))
stopifnot(all(is.na(sm_allna$mean_level)))
cat("ng_cross_portfolio_summary test passed\n")

# --- e2e: single-trait run gets portfolio + risk columns + diagnostics ---
set.seed(7)
n <- 16L; mk <- 60L; gid <- sprintf("P%02d", seq_len(n))
gm <- matrix(2L * rbinom(n * mk, 1, 0.5), n, mk, dimnames = list(gid, sprintf("M%03d", seq_len(mk))))
y  <- as.numeric(gm %*% rnorm(mk, 0, 0.1)) + rnorm(n)
genotype  <- data.frame(NAME = gid, gm, check.names = FALSE, stringsAsFactors = FALSE)
phenotype <- data.frame(NAME = gid, yield = y, stringsAsFactors = FALSE)
runmm <- data.frame(SNP = colnames(gm), chr = rep(1:2, length.out = mk),
                    bp = rep(seq(0, 100, length.out = 30), 2)[seq_len(mk)] * 1e6)
direction <- data.frame(trait = "yield", column = "yield", direction = "increase")
res <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "usefulness", n_crosses = 8L,
  max_crosses_per_parent = 3L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE,
  seed = 5L)
sc <- res$selected_crosses
stopifnot(all(c("cross_level","cross_upside","cross_confidence","risk_bin",
                "confidence_method","portfolio_profile") %in% names(sc)))
stopifnot(all(sc$cross_upside >= 0), diff(range(sc$cross_confidence, na.rm = TRUE)) > 0)
stopifnot(sc$confidence_method[[1L]] == "midparent_pev_partial")   # usefulness = effect-based X
stopifnot(!is.null(res$priority_risk_diagnostics))
stopifnot(all(c("cross_level","cross_upside","risk_bin","portfolio_profile") %in% names(res$candidate_crosses)))
# mean metric -> mean fully covers risk
res_m <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "mean", n_crosses = 8L,
  max_crosses_per_parent = 3L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L)
stopifnot(res_m$selected_crosses$confidence_method[[1L]] == "midparent_pev")
# le metric -> deprecated alias for parent_distance; must also get the mean-covers path
# (not midparent_pev_partial), since it maps to the same effect-free variance column.
res_le <- ng_run_cross_prediction(
  phenotype = phenotype, genotype = genotype, marker_map = runmm, trait_direction = direction,
  id_col = "NAME", map_marker_col = "SNP", map_chr_col = "chr", map_pos_col = "bp",
  map_pos_cm_divisor = 1e6, trait_value_metric = "le", n_crosses = 8L,
  max_crosses_per_parent = 3L, use_ocs = TRUE, write_outputs = FALSE, write_figures = FALSE, seed = 5L)
stopifnot(res_le$selected_crosses$confidence_method[[1L]] == "midparent_pev")
cat("priority risk portfolio e2e test passed\n")

# --- json export: priority_risk_diagnostics round-trip ---
env <- list(selected_crosses = res$selected_crosses,
            priority_risk_diagnostics = res$priority_risk_diagnostics)
sanitize <- function(x) { if (is.data.frame(x)) { for (j in seq_along(x))
  if (is.numeric(x[[j]])) x[[j]][!is.finite(x[[j]])] <- NA; return(x) }
  if (is.list(x)) return(lapply(x, sanitize)); if (is.numeric(x)) x[!is.finite(x)] <- NA; x }
tf <- tempfile(fileext = ".json")
jsonlite::write_json(sanitize(env), tf, auto_unbox = TRUE, na = "null", null = "null",
                     dataframe = "rows", digits = 8)
back <- jsonlite::read_json(tf, simplifyVector = TRUE)
stopifnot("portfolio_profile" %in% names(back$selected_crosses))
stopifnot(!is.null(back$priority_risk_diagnostics$confidence_method))
cat("priority risk portfolio json test passed\n")
