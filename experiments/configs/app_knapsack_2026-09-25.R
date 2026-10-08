# Application (a): two heuristics for random 0/1 knapsack instances. FROZEN 2026-09-25 as a
# self-contained stand-in until the lab's lcdaGRASP (GRASP vs Reactive GRASP) data is
# available; the protocol is identical, only `gen` changes.
# Population P: n_items ~ {50, 100, 200}, weights ~ U(1, 100), values correlated with
# weights (value = weight + U(0, 50)), capacity = 0.5 * sum(weights).
# A = greedy by value/weight ratio; B = GRASP-like randomised greedy (RCL alpha = 0.3)
# with 20 restarts and first-improvement swap local search, seeded by `seed`.
# Loss = min(gap / 0.05, 1) with gap = (UB - value) / UB and UB the Dantzig LP bound: the
# relative gap truncated at 5 % (gaps here are ~0.1-1 %, so bounds [0, 1] on the raw gap would
# waste the scale; truncation is the PAR-k-type operational choice of the formulation, §7).
# Margin 0.02 in truncated units = 0.1 percentage points of gap. Cost = number
# of objective evaluations (A: 1; B: 20 restarts * (1 + local search moves)). Data are
# expensive to generate, so the runner caches each (dgp, rep) stream (cache_data = TRUE).
local({
  greedy <- function(w, v, cap, order) {
    take <- logical(length(w)); load <- 0
    for (i in order) if (load + w[i] <= cap) { take[i] <- TRUE; load <- load + w[i] }
    take
  }
  local_search <- function(take, w, v, cap) {
    evals <- 0L; improved <- TRUE
    while (improved) {
      improved <- FALSE
      ins <- which(!take); outs <- which(take); load <- sum(w[take])
      for (i in ins) {
        if (load + w[i] <= cap) { take[i] <- TRUE; load <- load + w[i]; evals <- evals + 1L; improved <- TRUE; break }
        for (j in outs) { evals <- evals + 1L
          if (load - w[j] + w[i] <= cap && v[i] > v[j]) { take[j] <- FALSE; take[i] <- TRUE; load <- load - w[j] + w[i]; improved <- TRUE; break } }
        if (improved) break
      }
    }
    list(take = take, evals = evals)
  }
  gen_instance <- function(seed) {
    n <- sample(c(50L, 100L, 200L), 1); w <- runif(n, 1, 100); v <- w + runif(n, 0, 50); cap <- 0.5 * sum(w)
    r <- v / w; o <- order(-r)
    cw <- cumsum(w[o]); k <- which(cw > cap)[1]
    ub <- sum(v[o][seq_len(k - 1)]) + (cap - if (k > 1) cw[k - 1] else 0) * r[o][k]
    # A: deterministic greedy
    ta <- greedy(w, v, cap, o); fa <- sum(v[ta]); ca <- 1L
    # B: randomised greedy with RCL + local search, seeded
    set.seed(seed); best <- -Inf; cb <- 0L
    for (rep in 1:20) {
      remaining <- seq_len(n); ord <- integer(0)
      while (length(remaining)) {
        rr <- r[remaining]; thr <- max(rr) - 0.3 * (max(rr) - min(rr))
        pick <- remaining[rr >= thr]; pick <- if (length(pick) > 1) sample(pick, 1) else pick
        ord <- c(ord, pick); remaining <- setdiff(remaining, pick)
      }
      tb <- greedy(w, v, cap, ord); ls <- local_search(tb, w, v, cap); cb <- cb + 1L + ls$evals
      best <- max(best, sum(v[ls$take]))
    }
    c(loss_a = min((ub - fa) / ub / 0.05, 1), loss_b = min((ub - best) / ub / 0.05, 1), cost = ca + cb)
  }
  list(
    name = "app_knapsack_2026-09-25", seed = 20260925L,
    alpha = 0.05, margin = 0.02, bounds = c(0, 1), n_max = 600L, R = 100L, cache_data = TRUE,
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "fixed:final"),
    dgp = data.frame(scenario = "greedy_vs_grasp"),
    mu_of = function(p) NA_real_,
    truth_of = function(p, margin) NA_character_,
    gen = function(p, n) {
      seeds <- sample.int(1e9, n)
      l <- t(vapply(seeds, gen_instance, numeric(3)))
      data.frame(instance = seq_len(n), seed = seeds, loss_a = l[, 1], loss_b = l[, 2], cost = l[, 3])
    }
  )
})
