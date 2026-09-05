helper <- c(file.path("tests", "helper_load.R"), "helper_load.R", file.path("nextgen_cross_design", "tests", "helper_load.R"), file.path("..", "tests", "helper_load.R"))
source(helper[file.exists(helper)][[1L]])

diag <- data.frame(
  candidate_i = 1:3,
  family = c("var_simple", "recomb_gebv", "portfolio"),
  router_score = c(1.00, 0.75, 0.30),
  plan_z = c(-0.20, 1.10, 0.20),
  gain_z = c(0.00, 0.90, 0.10),
  stringsAsFactors = FALSE
)

guarded <- ng_meta_router_regret_guard(
  diag,
  selected_index = 1L,
  enabled = TRUE,
  plan_weight = 0.60,
  gain_weight = 0.40,
  min_advantage = 0.75,
  max_router_penalty = 0.40
)
stopifnot(isTRUE(guarded$override))
stopifnot(guarded$selected_index == 2L)
stopifnot(guarded$original_index == 1L)
stopifnot(guarded$selected_family == "recomb_gebv")
stopifnot(guarded$original_family == "var_simple")
stopifnot(is.finite(guarded$advantage))
stopifnot(is.finite(guarded$router_penalty))

blocked <- ng_meta_router_regret_guard(
  diag,
  selected_index = 1L,
  enabled = TRUE,
  plan_weight = 0.60,
  gain_weight = 0.40,
  min_advantage = 0.75,
  max_router_penalty = 0.10
)
stopifnot(!isTRUE(blocked$override))
stopifnot(blocked$selected_index == 1L)

disabled <- ng_meta_router_regret_guard(
  diag,
  selected_index = 1L,
  enabled = FALSE,
  plan_weight = 0.60,
  gain_weight = 0.40,
  min_advantage = 0.75,
  max_router_penalty = 0.40
)
stopifnot(!isTRUE(disabled$override))
stopifnot(disabled$selected_index == 1L)

cat("meta router regret guard smoke test passed\n")
