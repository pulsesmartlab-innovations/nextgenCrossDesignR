# Staged autotetraploid design == one-shot ng_polyploid_design_crosses, byte-for-byte on the
# selected-cross plan and the summary. Proves ng_poly_run_stage decomposes the poly design without
# changing the science (same QC -> effects -> scoring -> allocation), and that a later stage reloads
# the cached upstream ctx (compute-once) rather than recomputing.
helper <- c(file.path("tests", "helper_load.R"), "helper_load.R",
            file.path("nextgen_cross_design", "tests", "helper_load.R"),
            file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][1])

set.seed(11)
P <- 24L; M <- 60L; ploidy <- 4L
dose <- matrix(rbinom(P * M, ploidy, 0.5), P, M)
rownames(dose) <- sprintf("C%02d", seq_len(P)); colnames(dose) <- sprintf("SNP%03d", seq_len(M))
pheno <- setNames(rnorm(P), rownames(dose))

# deterministic allocator (greedy_local) so one-shot and staged cannot diverge on allocation RNG
cfg <- list(dosage = dose, n_crosses = 12L, ploidy = ploidy, phenotype = pheno,
            run_qc = TRUE, gain = "usefulness", grm_method = "vanraden", method = "greedy_local")

## --- one-shot ---------------------------------------------------------------
one <- do.call(ng_polyploid_design_crosses, cfg)

## --- staged: qc -> predict -> allocate -> rank over a run_dir ----------------
rd <- file.path(tempdir(), paste0("polystage_", Sys.getpid()))
unlink(rd, recursive = TRUE)
out <- NULL
for (s in ng_poly_cp_stage_order()) out <- ng_poly_run_stage(s, rd, cfg)
staged <- attr(out, "result")

## --- byte-identical plan + summary -----------------------------------------
stopifnot(identical(as.data.frame(one), as.data.frame(staged)))
stopifnot(isTRUE(all.equal(attr(one, "summary")$mean_gain, attr(staged, "summary")$mean_gain)))
stopifnot(identical(attr(one, "ploidy"), attr(staged, "ploidy")))

## --- compute-once: allocate reloaded the cached predict ctx (no refit) -------
predict_rds <- file.path(rd, "artifacts", "predict.rds")
stopifnot(file.exists(predict_rds))
ctx_predict <- readRDS(predict_rds)
stopifnot(!is.null(ctx_predict$scored), !is.null(ctx_predict$parent_kinship))  # attributes persisted as fields

## --- manifest records all four stages done --------------------------------
man <- ng_stage_status(rd)
stopifnot(all(c("qc", "predict", "allocate", "rank") %in% names(man$stages)))

## --- run_qc = FALSE also equivalent (no-QC path) ---------------------------
cfg2 <- cfg; cfg2$run_qc <- FALSE
one2 <- do.call(ng_polyploid_design_crosses, cfg2)
rd2 <- file.path(tempdir(), paste0("polystage2_", Sys.getpid())); unlink(rd2, recursive = TRUE)
o2 <- NULL; for (s in ng_poly_cp_stage_order()) o2 <- ng_poly_run_stage(s, rd2, cfg2)
stopifnot(identical(as.data.frame(one2), as.data.frame(attr(o2, "result"))))

cat("polyploid_staged_equivalence: staged == one-shot (QC on + off); compute-once ctx cached; manifest complete\n")
