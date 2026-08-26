helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

ids <- paste0("P", 1:6)
K <- diag(length(ids))
rownames(K) <- colnames(K) <- ids
K["P1", "P2"] <- K["P2", "P1"] <- 0.7
K["P2", "P3"] <- K["P3", "P2"] <- 0.65
K["P4", "P5"] <- K["P5", "P4"] <- 0.8
criterion <- data.frame(
  Genotype = ids,
  Criterion = c(4, 8, 6, 7, 3, 5),
  stringsAsFactors = FALSE
)

keep <- ng_simplemating_relate_thinning(K, criterion, threshold = 0.5, max_per_cluster = 1L)
stopifnot(identical(sort(keep), c("P2", "P4", "P6")))
stopifnot(attr(keep, "summary")$n_kept == 3L)

crosses <- ng_simplemating_build_crosses(
  moms = c("P1", "P2", "P4"),
  dads = ids,
  keep = keep,
  criterion = criterion,
  parent_kinship = K,
  include_self = FALSE,
  max_pair_kinship = 0.75
)
stopifnot(nrow(crosses) > 0L)
stopifnot(all(c("parent1", "parent2", "criterion_mean", "Y", "pair_kinship") %in% names(crosses)))
stopifnot(!any(crosses$parent1 == crosses$parent2))
stopifnot(all(crosses$pair_kinship <= 0.75))

past <- data.frame(parent1 = "P2", parent2 = "P4", stringsAsFactors = FALSE)
filtered <- ng_simplemating_past_thinning(crosses, past)
stopifnot(attr(filtered, "n_removed") == sum(ng_unordered_pair_key(crosses$parent1, crosses$parent2) == ng_unordered_pair_key("P2", "P4")))
stopifnot(!any(ng_unordered_pair_key(filtered$parent1, filtered$parent2) == ng_unordered_pair_key("P2", "P4")))

scores <- filtered
scores$cross_mean <- scores$criterion_mean
scores$parent_distance <- seq_len(nrow(scores)) / 10
scores <- ng_popvar_style_scores(scores, selection_prop = 0.2, out_prefix = "popvar_style")
stopifnot(all(is.finite(scores$popvar_style_uc)))
stopifnot(all(scores$popvar_style_status == "native_proxy"))

plan_simple <- ng_simplemating_select_crosses_native(
  scores = scores,
  score_col = "Y",
  n_crosses = min(2L, nrow(scores)),
  parent_kinship = K,
  max_crosses_per_parent = 2L,
  culling_pairwise_k = 0.75,
  method = "greedy_local"
)
stopifnot(nrow(plan_simple) == min(2L, nrow(scores)))
stopifnot(identical(attr(plan_simple, "summary")$simplemating_style, "native_proxy"))

plan_alpha <- ng_alphamate_style_select(
  scores = scores,
  criterion_col = "cross_mean",
  n_crosses = min(2L, nrow(scores)),
  parent_kinship = K,
  mode = "ModeOptTarget1",
  target_degree = 45,
  max_contributions = 2L,
  method = "greedy_local"
)
stopifnot(nrow(plan_alpha) == min(2L, nrow(scores)))
stopifnot(identical(attr(plan_alpha, "summary")$alphamate_style, "native_proxy"))

fallback_scores <- data.frame(
  parent1 = c("P1", "P1"),
  parent2 = c("P2", "P3"),
  cross_mean = c(1, 2),
  parent_distance = c(0.2, 0.4),
  stringsAsFactors = FALSE
)
fallback_scores <- ng_apply_popvar_native_proxy(fallback_scores, tail_p = 0.1)
stopifnot(all(is.finite(fallback_scores$popvar_uc)))
fallback_scores <- ng_apply_simplemating_native_proxy(fallback_scores, prop_sel = 0.1)
stopifnot(all(is.finite(fallback_scores$simple_usefa)))

dummy_geno <- matrix(c(0, 2, 1, 2, 0, 1), nrow = 3, byrow = TRUE,
                     dimnames = list(c("P1", "P2", "P3"), c("M1", "M2")))
dummy_effects <- list(beta = c(M1 = 0.1, M2 = 0.2))
dummy_map <- data.frame(marker = c("M1", "M2"), chr = c(1, 1), pos_cm = c(0, 10))
native_popvar <- ng_add_popvar_scores(
  scores = data.frame(
    parent1 = c("P1", "P1"),
    parent2 = c("P2", "P3"),
    cross_mean = c(1.0, 1.5),
    parent_distance = c(0.2, 0.3),
    stringsAsFactors = FALSE
  ),
  geno = dummy_geno,
  effects = dummy_effects,
  marker_map = dummy_map,
  tail_p = 0.1,
  engine = "native"
)
stopifnot(all(native_popvar$popvar_status == "native_proxy"))
stopifnot(all(is.finite(native_popvar$popvar_uc)))

native_simple <- ng_add_simplemating_scores(
  scores = data.frame(
    parent1 = c("P1", "P1"),
    parent2 = c("P2", "P3"),
    cross_mean = c(1.0, 1.5),
    parent_distance = c(0.2, 0.3),
    stringsAsFactors = FALSE
  ),
  geno = dummy_geno,
  effects = dummy_effects,
  marker_map = dummy_map,
  prop_sel = 0.1,
  engine = "native"
)
stopifnot(all(native_simple$simple_status == "native_proxy"))
stopifnot(all(is.finite(native_simple$simple_usefa)))

native_external <- ng_add_external_baseline_scores(
  scores = data.frame(
    parent1 = c("P1", "P1"),
    parent2 = c("P2", "P3"),
    cross_mean = c(1.0, 1.5),
    parent_distance = c(0.2, 0.3),
    stringsAsFactors = FALSE
  ),
  geno = dummy_geno,
  effects = dummy_effects,
  marker_map = dummy_map,
  methods = c("popvar_uc_native_check", "simple_usefa_native_check"),
  selection_prop = 0.1,
  popvar_engine = "native",
  simplemating_engine = "native"
)
stopifnot(all(native_external$popvar_status == "native_proxy"))
stopifnot(all(native_external$simple_status == "native_proxy"))
stopifnot(all(is.finite(native_external$popvar_uc)))
stopifnot(all(is.finite(native_external$simple_usefa)))

target_ids <- paste0("T", 1:5)
target_K <- matrix(0.02, nrow = length(target_ids), ncol = length(target_ids),
                   dimnames = list(target_ids, target_ids))
diag(target_K) <- 1
target_K["T1", "T2"] <- target_K["T2", "T1"] <- 0.90
target_K["T1", "T3"] <- target_K["T3", "T1"] <- 0.80
target_K["T2", "T3"] <- target_K["T3", "T2"] <- 0.80
target_scores <- data.frame(
  parent1 = c("T1", "T1", "T2", "T4", "T4", "T5"),
  parent2 = c("T2", "T3", "T3", "T5", "T1", "T2"),
  cross_mean = c(10.0, 9.8, 9.7, 6.0, 6.2, 6.1),
  pair_kinship = c(0.90, 0.80, 0.80, 0.02, 0.02, 0.02) / 2,
  stringsAsFactors = FALSE
)
gain_target <- ng_alphamate_style_select(
  scores = target_scores,
  criterion_col = "cross_mean",
  n_crosses = 2L,
  parent_kinship = target_K,
  mode = "ModeOptTarget1",
  target_degree = 0,
  max_contributions = 2L,
  method = "greedy_local",
  local_iter = 100L
)
diverse_target <- ng_alphamate_style_select(
  scores = target_scores,
  criterion_col = "cross_mean",
  n_crosses = 2L,
  parent_kinship = target_K,
  mode = "ModeOptTarget1",
  target_degree = 100,
  max_contributions = 2L,
  method = "greedy_local",
  local_iter = 100L
)
gain_summary <- attr(gain_target, "summary")
diverse_summary <- attr(diverse_target, "summary")
stopifnot(isTRUE(gain_summary$mean_gain >= diverse_summary$mean_gain))
stopifnot(isTRUE(diverse_summary$group_coancestry < gain_summary$group_coancestry))
stopifnot(isTRUE(diverse_summary$alphamate_frontier_n > 1L))
stopifnot(identical(diverse_summary$alphamate_selection_rule, "target_frontier"))

cat("external compatibility tests passed\n")
