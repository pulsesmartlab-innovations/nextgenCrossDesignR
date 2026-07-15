# Progeny inbreeding as a first-class metric + histogram.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

# --- 1. expected_progeny_inbreeding == parental coancestry (pair_kinship / 2) ---
# Hand-built relationship matrix so we know the answer exactly.
parents <- c("A", "B", "C", "D")
K <- matrix(c(
  1.00, 0.50, 0.10, 0.00,
  0.50, 1.00, 0.20, 0.05,
  0.10, 0.20, 1.00, 0.30,
  0.00, 0.05, 0.30, 1.00), 4, 4, byrow = TRUE,
  dimnames = list(parents, parents))
pairs <- data.frame(parent1 = c("A", "A", "C"), parent2 = c("B", "C", "D"),
                    stringsAsFactors = FALSE)
rv <- ng_pair_relationship_variance(pairs, K)
# ng_score_crosses computes expected_progeny_inbreeding = pair_kinship / 2; check the
# relationship the metric encodes directly against the hand matrix.
expected_f <- c(0.50, 0.10, 0.30) / 2
stopifnot(isTRUE(all.equal(rv$pair_kinship, c(0.50, 0.10, 0.30))))
stopifnot(isTRUE(all.equal(pmax(0, rv$pair_kinship / 2), expected_f)))

# --- 2. histogram bins sum to n, mean attribute matches ---
set.seed(7)
vals <- runif(200, 0, 0.4)
hgram <- ng_progeny_inbreeding_histogram(vals, breaks = 10L)
stopifnot(sum(hgram$count) == length(vals))
stopifnot(isTRUE(all.equal(sum(hgram$proportion), 1)))
stopifnot(isTRUE(all.equal(attr(hgram, "mean_progeny_inbreeding"), mean(vals))))
# data.frame input path (uses expected_progeny_inbreeding column)
df <- data.frame(expected_progeny_inbreeding = vals)
hgram2 <- ng_progeny_inbreeding_histogram(df, breaks = 10L)
stopifnot(sum(hgram2$count) == nrow(df))

# --- 3. lambda_progeny_inbreeding > 0 lowers mean progeny F in the chosen plan ---
set.seed(303)
np <- 20L
pp <- sprintf("P%02d", seq_len(np))
cmb <- t(utils::combn(np, 2L))
scores <- data.frame(parent1 = pp[cmb[, 1]], parent2 = pp[cmb[, 2]],
                     stringsAsFactors = FALSE)
scores$uc_dh_gebv <- rnorm(nrow(scores), 5, 1)
L <- matrix(rnorm(np * np, 0, 0.3), np, np)
G <- crossprod(L) / np
diag(G) <- diag(G) + 1
dimnames(G) <- list(pp, pp)
scores$pair_kinship <- G[cbind(match(scores$parent1, pp), match(scores$parent2, pp))]
scores$expected_progeny_inbreeding <- pmax(0, scores$pair_kinship / 2)

base <- ng_optimize_mating_plan(scores, 12L, parent_K = G, lambda_progeny_inbreeding = 0)
pen  <- ng_optimize_mating_plan(scores, 12L, parent_K = G, lambda_progeny_inbreeding = 50)
sb <- attr(base, "summary"); sp <- attr(pen, "summary")
stopifnot(is.finite(sb$mean_progeny_inbreeding), is.finite(sp$mean_progeny_inbreeding))
stopifnot(sp$mean_progeny_inbreeding <= sb$mean_progeny_inbreeding + 1e-8)
stopifnot(isTRUE(all.equal(sp$lambda_progeny_inbreeding, 50)))

# --- 4. no naming/semantic conflict with the pre-existing inbreeding vocabulary ---
# The package's only pre-existing "inbreed*" family is the parent-state assumption
# (assume_inbred / ng_audit_inbred_dosage): whether each PARENT is an inbred (homozygous)
# line, needed for the DH/RIL variance kernel. The new progeny-inbreeding family is a
# PAIR quantity: the coancestry between the two parents. They are different concepts on
# different objects, so they must never be the same name and must coexist independently.
set.seed(909)
ng <- 14L; mg <- 120L; gid <- sprintf("P%02d", seq_len(ng))
geno_g <- matrix(2L * rbinom(ng * mg, 1, 0.5), ng, mg,
                 dimnames = list(gid, sprintf("M%03d", seq_len(mg))))
yg <- as.numeric(geno_g %*% rnorm(mg, 0, 0.1)) + rnorm(ng); names(yg) <- gid
mmg <- data.frame(marker = colnames(geno_g), chr = rep(1:3, length.out = mg),
                  pos_cm = rep(seq(0, 100, length.out = 40), 3)[seq_len(mg)])
dg <- ng_design_crosses(geno_g, yg, marker_map = mmg, ids = gid, n_crosses = 8,
                        assume_inbred = TRUE, use_cpp = FALSE)

# (a) distinct kinds: assume_inbred is a parameter; expected_progeny_inbreeding is a column
stopifnot("assume_inbred" %in% names(formals(ng_design_crosses)))
stopifnot("expected_progeny_inbreeding" %in% names(dg$scores))
# (b) no bare `inbreeding` column silently created that could be confused
stopifnot(!("inbreeding" %in% names(dg$scores)))
# (c) the two are independent: the column is parent coancestry (pair_kinship/2), computed
#     regardless of the inbred assumption, and it varies across crosses (not a flag)
stopifnot(isTRUE(all.equal(dg$scores$expected_progeny_inbreeding,
                           pmax(0, dg$scores$pair_kinship / 2))))
stopifnot(diff(range(dg$scores$expected_progeny_inbreeding)) > 0)
# (d) the parameter-name families are disjoint stems: no ng_* API argument is named
#     exactly "inbreeding", and the parent-state stem "assume_inbred" is never reused for
#     the progeny quantity.
opt_args <- names(formals(ng_optimize_mating_plan))
stopifnot(!("inbreeding" %in% opt_args))            # no ambiguous bare name
stopifnot("lambda_progeny_inbreeding" %in% opt_args) # progeny penalty is fully qualified
stopifnot(!("assume_inbred" %in% opt_args))          # parent-state flag not leaked here

# --- 5. intended-use overlap with the pre-existing lambda_mating is handled honestly ---
# expected_progeny_inbreeding = max(0, pair_kinship/2), so lambda_progeny_inbreeding and
# lambda_mating act on the SAME parent-pair axis and add. Setting both must (a) flag the
# overlap in the summary and (b) warn once.
set.seed(112)
no <- 18L; po <- sprintf("P%02d", seq_len(no))
cbo <- t(utils::combn(no, 2L))
so <- data.frame(parent1 = po[cbo[, 1]], parent2 = po[cbo[, 2]], stringsAsFactors = FALSE)
so$uc_dh_gebv <- rnorm(nrow(so), 8, 2)
Lo <- matrix(rnorm(no * no, 0, 0.3), no, no); Go <- crossprod(Lo) / no
diag(Go) <- diag(Go) + 1; dimnames(Go) <- list(po, po)
so$pair_kinship <- Go[cbind(match(so$parent1, po), match(so$parent2, po))]
so$expected_progeny_inbreeding <- pmax(0, so$pair_kinship / 2)

options(ngcd.warned_relatedness_overlap = NULL)  # reset session-once guard
got_warn <- FALSE
both <- withCallingHandlers(
  ng_optimize_mating_plan(so, 10L, parent_K = Go, lambda_mating = 1,
                          lambda_progeny_inbreeding = 10),
  warning = function(w) { if (grepl("both penalize parent-pair", conditionMessage(w))) got_warn <<- TRUE; invokeRestart("muffleWarning") })
stopifnot(got_warn)
stopifnot(isTRUE(attr(both, "summary")$relatedness_penalty_overlap))
# using only one knob does NOT flag overlap
one <- ng_optimize_mating_plan(so, 10L, parent_K = Go, lambda_progeny_inbreeding = 10)
stopifnot(isFALSE(attr(one, "summary")$relatedness_penalty_overlap))

cat("progeny inbreeding test passed\n")
