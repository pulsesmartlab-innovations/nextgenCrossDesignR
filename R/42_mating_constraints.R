# Breeder mating constraints, applied in the fast "fixing" (repair) style rather than as
# objective penalties (Kinghorn 2011 found repair ~100x faster than penalising, and it is
# what our greedy/repair/evolution paths already do). These operate on the integer
# `selected` cross-index vector every allocator produces, so they apply uniformly across
# greedy_local / repair_local / evolution / MIP.
#
#   * ng_filter_group_permission  - candidate pre-filter: drop crosses whose parents'
#     groups are not an allowed group x group pairing (mating-group legality).
#   * ng_enforce_min_use_if_used   - a parent used at all must be used >= min_use times,
#     else it is dropped to zero (anti-singleton batch economics). Wired to the
#     previously-unused min_crosses_per_parent argument.
#   * ng_enforce_group_quota       - cap the number of crosses per group x group pairing.
#   * ng_resolve_committed_crosses / ng_apply_committed_crosses - lock a set of
#     already-decided matings into the plan and optimise the rest around them.
#
# Committed crosses are treated as PROTECTED throughout: the min-use, quota, capacity
# and trim passes never drop them.

# Unordered group-pair key for two group labels (vectorised).
ng_group_pair_key <- function(g1, g2) {
  g1 <- as.character(g1)
  g2 <- as.character(g2)
  lo <- ifelse(g1 <= g2, g1, g2)
  hi <- ifelse(g1 <= g2, g2, g1)
  paste(lo, hi, sep = "||")
}

# Map each candidate cross to its parents' group labels. parent_group is a named vector
# (parent id -> group). The public permission filter rejects missing parent/group mappings.
ng_cross_groups <- function(scores, parent_group) {
  pg <- parent_group[as.character(scores$parent1)]
  qg <- parent_group[as.character(scores$parent2)]
  list(g1 = unname(pg), g2 = unname(qg))
}

ng_validate_parent_group_mapping <- function(scores, parent_group) {
  if (is.null(parent_group)) ng_stop("parent_group is required for group constraints")
  if (is.null(names(parent_group)) || anyDuplicated(names(parent_group))) {
    ng_stop("parent_group must be a uniquely named vector indexed by parent ID")
  }
  requested <- unique(c(as.character(scores$parent1), as.character(scores$parent2)))
  missing_parent <- setdiff(requested, names(parent_group))
  if (length(missing_parent)) {
    ng_stop("parent_group is missing candidate parent IDs: ",
            paste(utils::head(missing_parent, 6L), collapse = ", "))
  }
  mapped <- as.character(parent_group[requested])
  if (anyNA(mapped) || any(!nzchar(trimws(mapped)))) {
    ng_stop("every candidate parent must have a non-missing group for group constraints")
  }
  invisible(TRUE)
}

ng_filter_group_permission <- function(scores, parent_group, group_permission) {
  if (is.null(parent_group) && is.null(group_permission)) return(scores)
  if (is.null(group_permission)) return(scores)
  ng_validate_parent_group_mapping(scores, parent_group)
  gr <- ng_cross_groups(scores, parent_group)
  g1 <- gr$g1; g2 <- gr$g2
  gp <- as.matrix(group_permission)
  if (is.null(rownames(gp)) || is.null(colnames(gp))) {
    ng_stop("group_permission must have group row and column names")
  }
  used_groups <- unique(c(g1, g2))
  if (anyNA(used_groups) || any(!nzchar(trimws(used_groups)))) {
    ng_stop("every candidate parent must have a non-missing group when group_permission is supplied")
  }
  missing_group <- union(setdiff(used_groups, rownames(gp)), setdiff(used_groups, colnames(gp)))
  if (length(missing_group)) {
    ng_stop("group_permission does not cover groups: ", paste(missing_group, collapse = ", "))
  }
  allowed <- vapply(seq_len(nrow(scores)), function(i) {
    a <- g1[[i]]; b <- g2[[i]]
    isTRUE(as.logical(gp[a, b])) || isTRUE(as.logical(gp[b, a]))
  }, logical(1))
  scores[allowed, , drop = FALSE]
}

# Final, explicit audit after all repair stages. Protected/committed crosses may
# make a requested cap infeasible; that override must be reported, never hidden.
ng_audit_mating_constraints <- function(scores, selected, parents,
                                        max_crosses_per_parent,
                                        min_use = NA_real_,
                                        min_unique_parents = NULL,
                                        parent_group = NULL,
                                        group_quota = NULL,
                                        cost_col = NULL,
                                        budget = Inf) {
  selected <- as.integer(selected)
  counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
  capacity <- names(counts)[counts > max_crosses_per_parent]
  used <- names(counts)[counts > 0]
  min_use_bad <- if (is.finite(min_use) && min_use > 1) used[counts[used] < min_use] else character(0)
  min_unique_bad <- if (!is.null(min_unique_parents) && is.finite(min_unique_parents) &&
                        length(used) < min_unique_parents) {
    sprintf("%d<%d", length(used), as.integer(min_unique_parents))
  } else character(0)
  quota_bad <- character(0)
  if (!is.null(parent_group) && !is.null(group_quota) && length(group_quota)) {
    gr <- ng_cross_groups(scores, parent_group)
    tab <- table(ng_group_pair_key(gr$g1, gr$g2)[selected])
    audit_parts <- strsplit(names(group_quota), "||", fixed = TRUE)
    audit_keys <- vapply(audit_parts, function(x) ng_group_pair_key(x[[1L]], x[[2L]]), character(1L))
    names(group_quota) <- audit_keys
    quota_bad <- audit_keys[vapply(audit_keys, function(k) {
      observed <- if (is.na(tab[k])) 0 else as.numeric(tab[k])
      observed > as.numeric(group_quota[[k]])
    }, logical(1L))]
  }
  total_cost <- NA_real_
  over_budget <- FALSE
  if (!is.null(cost_col) && cost_col %in% names(scores)) {
    total_cost <- sum(as.numeric(scores[[cost_col]][selected]), na.rm = TRUE)
    over_budget <- is.finite(budget) && total_cost > budget + 1e-9
  }
  violations <- c(
    if (length(capacity)) paste0("parent_capacity:", paste(capacity, collapse = ",")),
    if (length(min_use_bad)) paste0("min_use:", paste(min_use_bad, collapse = ",")),
    if (length(min_unique_bad)) paste0("min_unique_parents:", min_unique_bad),
    if (length(quota_bad)) paste0("group_quota:", paste(quota_bad, collapse = ",")),
    if (over_budget) "budget"
  )
  list(
    all_hard_constraints_satisfied = !length(violations),
    violations = violations,
    capacity_violating_parents = capacity,
    min_use_violating_parents = min_use_bad,
    min_unique_parents_violation = min_unique_bad,
    quota_violations = quota_bad,
    total_cost = total_cost,
    over_budget = over_budget
  )
}

# Drop-based min-use-if-used: any used parent below the floor is banned and its
# (non-protected) crosses removed, then the plan is refilled without reintroducing a
# banned parent. Banning is monotone (a parent is never un-banned), so this terminates
# in at most #parents rounds. Dropping rather than force-promoting keeps the plan
# objective-consistent (we refill with the best remaining feasible crosses).
ng_enforce_min_use_if_used <- function(scores, selected, parents,
                                       max_crosses_per_parent, min_use, n_crosses,
                                       protected = integer(0)) {
  if (is.null(min_use) || !is.finite(min_use) || min_use <= 1) return(selected)
  min_use <- as.integer(min_use)
  p1 <- as.character(scores$parent1)
  p2 <- as.character(scores$parent2)
  banned <- character(0)
  protected <- intersect(protected, selected)
  guard <- 0L
  repeat {
    guard <- guard + 1L
    if (guard > length(parents) + 2L) break
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    used <- names(counts)[counts > 0]
    viol <- setdiff(used[counts[used] < min_use], banned)
    if (!length(viol)) break
    banned <- union(banned, viol)
    # protected crosses may pin a banned parent below the floor; we cannot drop those,
    # so leave them (a warning is emitted by the caller if the floor stays violated).
    drop <- selected[(p1[selected] %in% banned | p2[selected] %in% banned) &
                       !(selected %in% protected)]
    if (length(drop)) selected <- setdiff(selected, drop)
    blocked_idx <- match(banned, parents)
    blocked_idx <- blocked_idx[is.finite(blocked_idx)]
    selected <- ng_fill_best(scores, selected, parents, max_crosses_per_parent,
                             n_crosses, blocked_parents = blocked_idx)
  }
  selected
}

ng_enforce_group_quota <- function(scores, selected, parents, parent_group,
                                   group_quota, max_crosses_per_parent, n_crosses,
                                   protected = integer(0)) {
  if (is.null(group_quota) || !length(group_quota)) return(selected)
  ng_validate_parent_group_mapping(scores, parent_group)
  if (is.null(names(group_quota)) || any(!nzchar(names(group_quota))) || anyDuplicated(names(group_quota))) {
    ng_stop("group_quota must be a uniquely named vector keyed as 'group1||group2'")
  }
  key_parts <- strsplit(names(group_quota), "||", fixed = TRUE)
  if (any(lengths(key_parts) != 2L) || any(!nzchar(unlist(key_parts, use.names = FALSE)))) {
    ng_stop("group_quota names must use exactly 'group1||group2'")
  }
  known_groups <- unique(as.character(parent_group))
  key_groups <- unlist(key_parts, use.names = FALSE)
  unknown_groups <- setdiff(key_groups, known_groups)
  if (length(unknown_groups)) {
    ng_stop("group_quota names contain unknown groups: ", paste(unknown_groups, collapse = ", "))
  }
  canonical_keys <- vapply(key_parts, function(x) ng_group_pair_key(x[[1L]], x[[2L]]), character(1L))
  if (anyDuplicated(canonical_keys)) {
    ng_stop("group_quota contains duplicate unordered group-pair keys")
  }
  names(group_quota) <- canonical_keys
  quota_value <- suppressWarnings(as.numeric(group_quota))
  if (any(!is.finite(quota_value)) || any(quota_value < 0) ||
      any(abs(quota_value - round(quota_value)) > 1e-8)) {
    ng_stop("group_quota values must be finite non-negative integers")
  }
  group_quota[] <- as.integer(round(quota_value))
  gr <- ng_cross_groups(scores, parent_group)
  key <- ng_group_pair_key(gr$g1, gr$g2)
  quota_keys <- names(group_quota)
  protected <- intersect(protected, selected)
  guard <- 0L
  repeat {
    guard <- guard + 1L
    if (guard > length(group_quota) + 2L) break
    tab <- table(key[selected])
    over <- character(0)
    for (qk in quota_keys) {
      if (!is.na(tab[qk]) && tab[qk] > group_quota[[qk]]) over <- c(over, qk)
    }
    if (!length(over)) break
    # drop lowest-gain non-protected crosses in over-quota group pairs down to the cap
    for (qk in over) {
      in_pair <- selected[key[selected] == qk & !(selected %in% protected)]
      excess <- as.integer(tab[qk] - group_quota[[qk]])
      if (excess > 0 && length(in_pair)) {
        drop <- in_pair[order(scores$.linear_gain[in_pair])][seq_len(min(excess, length(in_pair)))]
        selected <- setdiff(selected, drop)
      }
    }
    # refill, but do not exceed any quota: block crosses in already-full pairs
    tab2 <- table(key[selected])
    full_keys <- quota_keys[vapply(quota_keys, function(qk) {
      !is.na(tab2[qk]) && tab2[qk] >= group_quota[[qk]]
    }, logical(1))]
    selected <- ng_fill_best(scores, selected, parents, max_crosses_per_parent,
                             n_crosses, blocked_pair_keys = full_keys, pair_keys = key)
  }
  selected
}

# Resolve committed matings (a data.frame with parent1/parent2, order-insensitive) to
# candidate row indices. Errors if a committed cross is not among the candidates.
ng_resolve_committed_crosses <- function(scores, committed_crosses) {
  if (is.null(committed_crosses) || !nrow(committed_crosses)) return(integer(0))
  cc <- as.data.frame(committed_crosses, stringsAsFactors = FALSE)
  if (!all(c("parent1", "parent2") %in% names(cc))) {
    ng_stop("committed_crosses must have parent1 and parent2 columns")
  }
  cand_key <- ng_group_pair_key(scores$parent1, scores$parent2)
  want_key <- ng_group_pair_key(cc$parent1, cc$parent2)
  idx <- match(want_key, cand_key)
  if (anyNA(idx)) {
    ng_stop("committed_crosses contains crosses that are not candidate crosses: ",
            paste(want_key[is.na(idx)], collapse = ", "))
  }
  unique(idx)
}

# Lock committed crosses into the plan and repair the rest around them: union them in,
# relieve any per-parent capacity overflow by dropping the lowest-gain NON-committed
# crosses, trim to n_crosses by dropping lowest-gain non-committed, then fill if short.
ng_apply_committed_crosses <- function(scores, selected, committed_idx, parents,
                                       max_crosses_per_parent, n_crosses) {
  if (!length(committed_idx)) return(selected)
  if (length(committed_idx) > n_crosses) {
    ng_stop("more committed_crosses than n_crosses")
  }
  selected <- union(committed_idx, selected)
  p1 <- as.character(scores$parent1)
  p2 <- as.character(scores$parent2)
  # capacity relief protecting committed
  guard <- 0L
  repeat {
    guard <- guard + 1L
    if (guard > 4L * n_crosses + 10L) break
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    over <- names(counts)[counts > max_crosses_per_parent]
    if (!length(over)) break
    droppable <- selected[(p1[selected] %in% over | p2[selected] %in% over) &
                            !(selected %in% committed_idx)]
    if (!length(droppable)) break  # only committed overflow capacity -> unavoidable
    drop <- droppable[which.min(scores$.linear_gain[droppable])]
    selected <- setdiff(selected, drop)
  }
  # trim to n_crosses, never dropping committed
  if (length(selected) > n_crosses) {
    extra <- setdiff(selected, committed_idx)
    extra <- extra[order(scores$.linear_gain[extra])]
    n_drop <- length(selected) - n_crosses
    selected <- setdiff(selected, head(extra, n_drop))
  }
  if (length(selected) < n_crosses) {
    selected <- ng_fill_best(scores, selected, parents, max_crosses_per_parent, n_crosses)
  }
  selected
}

# Hard budget cap on total cross cost. Drops the worst value-for-money (lowest
# .linear_gain per unit cost) NON-protected crosses until the plan is within budget, then
# refills empty slots with the best crosses that still fit the remaining budget and
# capacity. If the budget cannot support n_crosses, the shorter plan is returned (the
# caller warns) rather than silently overspending.
ng_enforce_budget <- function(scores, selected, cost_col, budget, parents,
                              max_crosses_per_parent, n_crosses, protected = integer(0)) {
  if (is.null(cost_col) || !(cost_col %in% names(scores)) || !is.finite(budget)) return(selected)
  cost <- as.numeric(scores[[cost_col]]); cost[!is.finite(cost)] <- 0
  protected <- intersect(protected, selected)
  guard <- 0L
  repeat {
    guard <- guard + 1L
    if (guard > 4L * n_crosses + 10L) break
    if (sum(cost[selected]) <= budget) break
    droppable <- setdiff(selected, protected)
    if (!length(droppable)) break  # only committed crosses exceed budget: unavoidable
    gpc <- scores$.linear_gain[droppable] / pmax(cost[droppable], 1e-9)
    selected <- setdiff(selected, droppable[which.min(gpc)])
  }
  if (length(selected) < n_crosses) {
    counts <- ng_parent_counts(scores[selected, , drop = FALSE], parents)
    remaining <- budget - sum(cost[selected])
    p1 <- as.character(scores$parent1); p2 <- as.character(scores$parent2)
    ord <- order(scores$.linear_gain, decreasing = TRUE)
    for (idx in ord) {
      if (length(selected) >= n_crosses) break
      if (idx %in% selected) next
      if (cost[idx] > remaining) next
      a <- p1[idx]; b <- p2[idx]
      if (a == b) {
        if (counts[a] + 2L > max_crosses_per_parent) next
      } else if (counts[a] + 1L > max_crosses_per_parent ||
                 counts[b] + 1L > max_crosses_per_parent) next
      selected <- c(selected, idx)
      counts[a] <- counts[a] + 1L; counts[b] <- counts[b] + 1L
      remaining <- remaining - cost[idx]
    }
  }
  selected
}
