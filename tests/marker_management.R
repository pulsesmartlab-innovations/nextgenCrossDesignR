# Marker steering and lethal-allele guarding.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

set.seed(404)
parents <- sprintf("P%02d", 1:8)
markers <- sprintf("M%02d", 1:5)
geno <- matrix(sample(0:2, length(parents) * length(markers), replace = TRUE),
               nrow = length(parents), dimnames = list(parents, markers))

# Deterministic marker for steering: M01 dosage increasing across parents.
geno[, "M01"] <- rep(0:1, length.out = length(parents)) * 2  # 0 or 2
# Lethal locus M05: make P01, P02 carriers (dosage 1), others homozygous-safe (0).
geno[, "M05"] <- 0
geno[c("P01", "P02"), "M05"] <- 1

pairs <- data.frame(
  parent1 = c("P01", "P03", "P05", "P01"),
  parent2 = c("P02", "P04", "P06", "P03"),
  stringsAsFactors = FALSE)

# --- 1. marker steering: expected progeny ALT freq and directional score ---
spec_inc <- ng_marker_target_spec("M01", direction = "increase", weight = 1)
mt <- ng_marker_target_scores(geno, pairs, spec_inc, ploidy = 2)
# freq = (d1 + d2) / (2*ploidy); check first cross P01xP02
f_expected <- (geno["P01", "M01"] + geno["P02", "M01"]) / (2 * 2)
stopifnot(isTRUE(all.equal(mt$marker_freq_M01[[1]], unname(f_expected))))
# increasing direction: higher freq => higher score (monotone in freq)
stopifnot(isTRUE(all.equal(mt$marker_target_score, mt$marker_freq_M01)))

# decrease direction flips the score to (1 - freq)
spec_dec <- ng_marker_target_spec("M01", direction = "decrease", weight = 1)
mt_dec <- ng_marker_target_scores(geno, pairs, spec_dec, ploidy = 2)
stopifnot(isTRUE(all.equal(mt_dec$marker_target_score, 1 - mt$marker_freq_M01)))

# target_freq: score is negative distance to target
spec_tgt <- ng_marker_target_spec("M01", target_freq = 0.5, weight = 2)
mt_tgt <- ng_marker_target_scores(geno, pairs, spec_tgt, ploidy = 2)
stopifnot(isTRUE(all.equal(mt_tgt$marker_target_score, -2 * abs(mt$marker_freq_M01 - 0.5))))

# --- 2. lethal-recessive: only P01xP02 is carrier x carrier ---
lspec <- ng_lethal_recessive_spec("M05", risk_allele = "alt")
lr <- ng_lethal_recessive_cross_risk(geno, pairs, lspec, ploidy = 2)
stopifnot(identical(lr$lethal_carrier_cross, c(TRUE, FALSE, FALSE, FALSE)))
stopifnot(identical(lr$lethal_risk_loci, c(1L, 0L, 0L, 0L)))

# --- 3. convenience wrapper: columns added, carrier x carrier dropped ---
scores <- pairs
scores$uc_dh_gebv <- c(10, 9, 8, 7)
scores$pair_kinship <- 0
aug <- ng_apply_marker_management(scores, geno,
                                  marker_target_spec = spec_inc,
                                  lethal_spec = lspec, ploidy = 2,
                                  lambda_marker = 0.5,
                                  drop_lethal_carrier_crosses = TRUE)
stopifnot("marker_target_score" %in% names(aug))
stopifnot("marker_adjusted_gain" %in% names(aug))
# P01xP02 dropped (carrier x carrier); the other three remain
stopifnot(nrow(aug) == 3L)
stopifnot(!any(aug$parent1 == "P01" & aug$parent2 == "P02"))

cat("marker management test passed\n")
