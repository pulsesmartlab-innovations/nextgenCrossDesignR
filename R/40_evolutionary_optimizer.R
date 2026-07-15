# Native evolutionary (memetic genetic-algorithm) mate-allocation optimizer.
#
# AlphaMate's real solver is a differential-evolution metaheuristic; the package
# previously only had a lambda-grid scalarization proxy. This adds a genuine native
# metaheuristic that optimizes the SAME objective as the MIP/greedy paths
# (ng_plan_objective_contribution): maximize
#   sum(.linear_gain[sel]) - lambda_group*(c'Kc) - lambda_parent_use*sum(c^2)
# subject to exactly n_crosses, per-parent capacity, and min_unique_parents.
#
# Encoding: an individual is the length-n_crosses integer vector of selected candidate
# cross-row indices (the exact object every existing repair/objective function
# consumes -- no decode layer). Constraints are handled by REPAIR (always feasible),
# not penalty. It is MEMETIC: each generation the elite is hill-climbed with the
# C++ ng_local_swap kernel, and the population is warm-started from greedy_local, so
# with elitism the result is guaranteed >= the greedy solution.
ng_evolutionary_mate_allocation <- function(scores,
                                            n_crosses,
                                            parents,
                                            parent_K,
                                            max_crosses_per_parent,
                                            min_unique_parents = NULL,
                                            lambda_group = 0,
                                            lambda_parent_use = 0,
                                            local_iter = 200L,
                                            evol_solutions = 100L,
                                            evol_iterations = 200L,
                                            evol_stop = 40L,
                                            p_mutation = 0.3,
                                            warm_start = TRUE,
                                            seed = NULL) {
  run <- function() {
    n_pairs <- nrow(scores)
    g <- scores$.linear_gain
    gsd <- stats::sd(g[is.finite(g)]); if (!is.finite(gsd) || gsd <= 0) gsd <- 1
    # Precompute integer parent indices + kinship matrix ONCE so the per-individual
    # fitness is pure vectorized arithmetic (tabulate + a small quadratic form) with
    # no data.frame subsetting -- the fitness is evaluated population x generation
    # times and this is what dominates runtime.
    np <- length(parents)
    p1i <- match(as.character(scores$parent1), parents)
    p2i <- match(as.character(scores$parent2), parents)
    Kmat <- as.matrix(parent_K[parents, parents, drop = FALSE]); storage.mode(Kmat) <- "double"
    total_slots <- 2 * n_crosses
    softmax_sample <- function(pool, k) {
      k <- min(k, length(pool)); if (k <= 0) return(integer(0))
      w <- exp((g[pool] - max(g[pool])) / gsd)
      w[!is.finite(w) | w < 0] <- 0
      if (sum(w) <= 0) w <- rep(1, length(pool))
      pool[sample.int(length(pool), k, replace = FALSE, prob = w)]
    }
    repair <- function(idx) {
      idx <- unique(as.integer(idx))
      idx <- ng_repair_capacity(scores, idx, parents, max_crosses_per_parent, n_crosses)
      ng_enforce_min_unique_parents(scores, idx, parents, max_crosses_per_parent, min_unique_parents, n_crosses)
    }
    fitness <- function(idx) {
      cnt <- tabulate(c(p1i[idx], p2i[idx]), nbins = np)
      base <- sum(g[idx])
      if (lambda_group <= 0 && lambda_parent_use <= 0) return(base)
      cv <- cnt / total_slots
      pen <- 0
      if (lambda_group > 0) pen <- pen + lambda_group * as.numeric(crossprod(cv, Kmat %*% cv))
      if (lambda_parent_use > 0) pen <- pen + lambda_parent_use * sum(cv * cv)
      base - pen
    }
    local_improve <- function(idx) {
      if (local_iter <= 0) return(idx)
      ng_local_swap(scores, idx, parents, parent_K, max_crosses_per_parent, lambda_group, local_iter,
                    lambda_parent_use = lambda_parent_use)
    }
    ord <- order(g, decreasing = TRUE)

    # --- initial population ---
    N <- max(4L, as.integer(evol_solutions))
    pop <- vector("list", N)
    pop[[1]] <- if (isTRUE(warm_start)) {
      ng_greedy_local(scores, n_crosses, parents, parent_K, max_crosses_per_parent,
                      min_unique_parents, lambda_group, local_iter, lambda_parent_use = lambda_parent_use)
    } else repair(ord[seq_len(n_crosses)])
    pop[[2]] <- repair(ord[seq_len(n_crosses)])
    for (i in 3:N) pop[[i]] <- repair(softmax_sample(seq_len(n_pairs), n_crosses))
    fits <- vapply(pop, fitness, numeric(1))
    best <- pop[[which.max(fits)]]; best_fit <- max(fits)

    tournament <- function() {
      cand <- sample.int(N, min(3L, N))
      cand[which.max(fits[cand])]
    }
    crossover <- function(a, b) {
      inter <- intersect(a, b)
      if (length(inter) >= n_crosses) return(inter[seq_len(n_crosses)])
      sdiff <- setdiff(union(a, b), inter)
      fill <- softmax_sample(sdiff, n_crosses - length(inter))
      c(inter, fill)
    }
    mutate <- function(idx) {
      k <- sample.int(3L, 1L)
      off <- setdiff(seq_len(n_pairs), idx)
      if (!length(off)) return(idx)
      drop <- if (length(idx) >= k) sample(idx, k) else idx
      add <- softmax_sample(off, k)
      c(setdiff(idx, drop), add)
    }

    n_elite <- max(1L, round(0.1 * N))
    no_improve <- 0L
    for (gen in seq_len(max(1L, as.integer(evol_iterations)))) {
      ord_f <- order(fits, decreasing = TRUE)
      new_pop <- pop[ord_f[seq_len(n_elite)]]
      new_pop[[1]] <- local_improve(new_pop[[1]])          # memetic hill-climb on the best
      while (length(new_pop) < N) {
        child <- crossover(pop[[tournament()]], pop[[tournament()]])
        if (stats::runif(1) < p_mutation) child <- mutate(child)
        new_pop[[length(new_pop) + 1L]] <- repair(child)
      }
      pop <- new_pop
      fits <- vapply(pop, fitness, numeric(1))
      gi <- which.max(fits)
      if (fits[gi] > best_fit + 1e-10) { best_fit <- fits[gi]; best <- pop[[gi]]; no_improve <- 0L }
      else no_improve <- no_improve + 1L
      if (no_improve >= max(1L, as.integer(evol_stop))) break
    }
    best <- repair(local_improve(best))
    best
  }
  if (is.null(seed)) run() else ng_with_rng_seed(seed, run())
}
