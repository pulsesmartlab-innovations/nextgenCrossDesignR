ng_meta_router_regret_guard <- function(diag,
                                        selected_index,
                                        enabled = TRUE,
                                        plan_weight = 0.60,
                                        gain_weight = 0.40,
                                        min_advantage = 0.75,
                                        max_router_penalty = 0.40) {
  n <- if (is.data.frame(diag)) nrow(diag) else 0L
  selected_index <- as.integer(selected_index[[1]])
  empty_result <- function(reason) {
    list(
      selected_index = selected_index,
      original_index = selected_index,
      override = FALSE,
      selected_family = NA_character_,
      original_family = NA_character_,
      guard_score = NA_real_,
      original_guard_score = NA_real_,
      advantage = NA_real_,
      router_penalty = NA_real_,
      reason = reason
    )
  }

  if (!isTRUE(enabled)) return(empty_result("disabled"))
  if (!n || !is.finite(selected_index) || selected_index < 1L || selected_index > n) {
    return(empty_result("invalid_selection"))
  }

  for (col in c("router_score", "plan_z", "gain_z")) {
    if (!(col %in% names(diag))) diag[[col]] <- NA_real_
  }
  if (!("family" %in% names(diag))) diag$family <- rep(NA_character_, n)

  plan_weight <- as.numeric(plan_weight[[1]])
  gain_weight <- as.numeric(gain_weight[[1]])
  if (!is.finite(plan_weight) || plan_weight < 0) plan_weight <- 0.60
  if (!is.finite(gain_weight) || gain_weight < 0) gain_weight <- 0.40
  weight_sum <- plan_weight + gain_weight
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    plan_weight <- 0.60
    gain_weight <- 0.40
    weight_sum <- 1
  }
  plan_weight <- plan_weight / weight_sum
  gain_weight <- gain_weight / weight_sum

  min_advantage <- as.numeric(min_advantage[[1]])
  if (!is.finite(min_advantage) || min_advantage < 0) min_advantage <- 0.75
  max_router_penalty <- as.numeric(max_router_penalty[[1]])
  if (!is.finite(max_router_penalty) || max_router_penalty < 0) max_router_penalty <- 0.40

  plan_z <- suppressWarnings(as.numeric(diag$plan_z))
  gain_z <- suppressWarnings(as.numeric(diag$gain_z))
  router_score <- suppressWarnings(as.numeric(diag$router_score))
  guard_score <- plan_weight * plan_z + gain_weight * gain_z
  guard_score[!is.finite(guard_score)] <- -Inf
  router_score[!is.finite(router_score)] <- -Inf

  original_guard <- guard_score[[selected_index]]
  original_router <- router_score[[selected_index]]
  if (!is.finite(original_guard) || !is.finite(original_router)) {
    out <- empty_result("invalid_original")
    out$selected_family <- as.character(diag$family[[selected_index]])
    out$original_family <- out$selected_family
    return(out)
  }

  advantage <- guard_score - original_guard
  router_penalty <- original_router - router_score
  eligible <- seq_len(n) != selected_index &
    is.finite(advantage) &
    is.finite(router_penalty) &
    advantage >= min_advantage &
    router_penalty <= max_router_penalty

  if (!any(eligible)) {
    other_idx <- seq_len(n) != selected_index
    best_advantage <- if (any(other_idx & is.finite(advantage))) {
      max(advantage[other_idx], na.rm = TRUE)
    } else {
      NA_real_
    }
    return(list(
      selected_index = selected_index,
      original_index = selected_index,
      override = FALSE,
      selected_family = as.character(diag$family[[selected_index]]),
      original_family = as.character(diag$family[[selected_index]]),
      guard_score = original_guard,
      original_guard_score = original_guard,
      advantage = best_advantage,
      router_penalty = NA_real_,
      reason = "no_guard_candidate"
    ))
  }

  eligible_idx <- which(eligible)
  best_rank <- order(
    -advantage[eligible_idx],
    router_penalty[eligible_idx],
    -router_score[eligible_idx]
  )
  guarded_index <- eligible_idx[[best_rank[[1]]]]

  list(
    selected_index = guarded_index,
    original_index = selected_index,
    override = TRUE,
    selected_family = as.character(diag$family[[guarded_index]]),
    original_family = as.character(diag$family[[selected_index]]),
    guard_score = guard_score[[guarded_index]],
    original_guard_score = original_guard,
    advantage = advantage[[guarded_index]],
    router_penalty = router_penalty[[guarded_index]],
    reason = "guard_override"
  )
}
