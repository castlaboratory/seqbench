# Application (distribution logistics): two constructive heuristics for the capacitated
# vehicle routing problem (CVRP). FROZEN 2026-10-23, before any loss was computed.
#
# Question. Which constructive heuristic should a distribution planner deploy on random
# delivery problems like these: the sweep algorithm (Gillett & Miller 1974) or the
# Clarke-Wright savings algorithm (Clarke & Wright 1964), both followed by 2-opt inside
# every route? Margin delta = 0.01: one percentage point of the gap to a lower bound.
#
# Population P (one draw per instance, generator in the spirit of Uchoa et al. 2017):
#   customers n ~ U{50..150} in the unit square; depot central (0.5, 0.5), eccentric
#   (0, 0) or uniform, with probability 1/3 each; customer positions uniform, clustered
#   (3-8 Gaussian clusters, sd 0.05, truncated to the square) or half/half; demands unitary,
#   U{1..10} or U{5..100}; capacity Q = max(max demand, ceiling(r * mean demand)) with
#   r ~ U(4, 25) customers per route on average. Euclidean distances.
# A = sweep + 2-opt: customers sorted by polar angle around the depot from a random start
#   angle (the instance's seed), routes filled greedily in that order.
# B = parallel Clarke-Wright savings + 2-opt.
# Loss L = 1 - LB / z in [0, 1): relative gap of the solution cost z to the lower bound
#   LB = max(MST of depot and customers, 2 * sum(q_i d_0i) / Q) (every CVRP solution
#   spans all nodes and visits each customer once; the second is the radial bound). The
#   bounds [0, 1] hold a priori for any feasible solution: no truncation, no clipping.
# Cost = CPU seconds of the two heuristics.
#
# Operating characteristics are measured on streams resampled from a frozen pool of
# `pool_n` i.i.d. instances (experiments/run/make_pool.R). The resampling population is the
# pool's empirical distribution, whose mean `mu_of()` is known exactly, so coverage and false
# declarations are measurable; the pool in generation order is the "live" stream of fresh
# instances reported in the article.
local({
  dist_mat <- function(xy) as.matrix(stats::dist(xy))
  route_len <- function(r, D) { p <- c(1L, r, 1L); sum(D[cbind(p[-length(p)], p[-1])]) }
  two_opt <- function(r, D) {                        # first improvement, vectorised over j
    if (length(r) < 3L) return(r)
    repeat {
      p <- c(1L, r, 1L); k <- length(p); improved <- FALSE
      for (i in 1:(k - 3L)) {
        j <- (i + 2L):(k - 1L)
        delta <- D[p[i], p[j]] + D[p[i + 1L], p[j + 1L]] - D[p[i], p[i + 1L]] - D[cbind(p[j], p[j + 1L])]
        b <- which(delta < -1e-12)
        if (length(b)) { jj <- j[b[1]]; p[(i + 1L):jj] <- rev(p[(i + 1L):jj]); improved <- TRUE; break }
      }
      r <- p[2:(k - 1L)]
      if (!improved) return(r)
    }
  }
  sweep <- function(xy, q, Q, D, start) {
    ang <- atan2(xy[-1, 2] - xy[1, 2], xy[-1, 1] - xy[1, 1])
    ang <- (ang - start) %% (2 * pi)
    ord <- order(ang) + 1L                           # node indices (depot = 1)
    routes <- list(); cur <- integer(0); load <- 0
    for (v in ord) {
      if (load + q[v] > Q) { routes[[length(routes) + 1L]] <- cur; cur <- integer(0); load <- 0 }
      cur <- c(cur, v); load <- load + q[v]
    }
    routes[[length(routes) + 1L]] <- cur
    lapply(routes, two_opt, D = D)
  }
  savings <- function(q, Q, D) {
    n <- nrow(D); cust <- 2:n
    pr <- which(upper.tri(D[cust, cust]), arr.ind = TRUE)
    i <- cust[pr[, 1]]; j <- cust[pr[, 2]]
    s <- D[1, i] + D[1, j] - D[cbind(i, j)]
    o <- order(-s); i <- i[o]; j <- j[o]; s <- s[o]
    keep <- s > 0; i <- i[keep]; j <- j[keep]
    rid <- seq_len(n); routes <- as.list(seq_len(n)); load <- q
    for (k in seq_along(i)) {
      a <- i[k]; b <- j[k]; ra <- rid[a]; rb <- rid[b]
      if (ra == rb || load[ra] + load[rb] > Q) next
      A <- routes[[ra]]; B <- routes[[rb]]
      if (A[length(A)] == a && B[1] == b) new <- c(A, B)
      else if (A[1] == a && B[length(B)] == b) new <- c(B, A)
      else if (A[length(A)] == a && B[length(B)] == b) new <- c(A, rev(B))
      else if (A[1] == a && B[1] == b) new <- c(rev(A), B)
      else next                                      # a or b is interior to its route
      routes[[ra]] <- new; routes[rb] <- list(NULL); load[ra] <- load[ra] + load[rb]; rid[new] <- ra
    }
    lapply(Filter(Negate(is.null), routes[cust]), two_opt, D = D)
  }
  mst_len <- function(D) {                           # Prim, O(n^2)
    n <- nrow(D); inT <- c(TRUE, rep(FALSE, n - 1L)); best <- D[1, ]; tot <- 0
    for (k in 2:n) { best[inT] <- Inf; v <- which.min(best); tot <- tot + best[v]; inT[v] <- TRUE; best <- pmin(best, D[v, ]) }
    tot
  }
  check_solution <- function(routes, n, q, Q) {
    v <- unlist(routes)
    stopifnot(length(v) == n - 1L, setequal(v, 2:n), all(vapply(routes, function(r) sum(q[r]), 0) <= Q))
  }
  instance <- function(seed) {
    set.seed(seed)
    n <- sample(50:150, 1)
    depot <- switch(sample(3, 1), c(0.5, 0.5), c(0, 0), runif(2))
    placement <- sample(c("uniform", "clustered", "mixed"), 1)
    clustered <- function(m) {
      k <- sample(3:8, 1); ctr <- matrix(runif(2 * k), k); id <- sample(k, m, replace = TRUE)
      pmin(pmax(ctr[id, , drop = FALSE] + matrix(stats::rnorm(2 * m, 0, 0.05), m), 0), 1)
    }
    pts <- switch(placement, uniform = matrix(runif(2 * n), n), clustered = clustered(n),
                  mixed = rbind(matrix(runif(2 * (n %/% 2)), n %/% 2), clustered(n - n %/% 2)))
    xy <- rbind(depot, pts)
    q <- c(0, switch(sample(3, 1), rep(1, n), sample(1:10, n, TRUE), sample(5:100, n, TRUE)))
    Q <- max(max(q), ceiling(runif(1, 4, 25) * mean(q[-1])))
    D <- dist_mat(xy)
    lb <- max(mst_len(D), 2 * sum(q * D[1, ]) / Q)
    start <- runif(1, 0, 2 * pi)
    cpu <- function() sum(proc.time()[c("user.self", "sys.self")])   # system.time() would gc() each call
    c0 <- cpu(); ra <- sweep(xy, q, Q, D, start); c1 <- cpu(); rb <- savings(q, Q, D); c2 <- cpu()
    ta <- c1 - c0; tb <- c2 - c1
    check_solution(ra, nrow(D), q, Q); check_solution(rb, nrow(D), q, Q)
    za <- sum(vapply(ra, route_len, 0, D = D)); zb <- sum(vapply(rb, route_len, 0, D = D))
    data.frame(seed = seed, loss_a = 1 - lb / za, loss_b = 1 - lb / zb, cost = ta + tb,
               n = n, placement = placement, routes_a = length(ra), routes_b = length(rb))
  }
  pool_file <- "experiments/out/pools/app_cvrp_2026-10-23.rds"
  pool <- local({ cache <- NULL; function() { if (is.null(cache)) cache <<- readRDS(pool_file); cache } })
  list(
    name = "app_cvrp_2026-10-23", seed = 20261023L,
    alpha = 0.05, margin = 0.01, bounds = c(0, 1), n_max = 3000L, R = 1000L,
    pool_n = 10000L, pool_file = pool_file, instance = instance, savi_pilot_n = 200L,
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "fixed:100", "fixed:final", "savi_cs"),
    dgp = data.frame(scenario = "sweep_vs_savings"),
    mu_of = function(p) { P <- pool(); mean(P$loss_a - P$loss_b) },
    truth_of = function(p, margin) {
      P <- pool(); mu <- mean(P$loss_a - P$loss_b)
      if (abs(mu) <= margin) "equivalent" else if (mu > 0) "B" else "A"
    },
    gen = function(p, n) {                            # i.i.d. resampling from the frozen pool
      P <- pool(); i <- sample.int(nrow(P), n, replace = TRUE)
      data.frame(instance = seq_len(n), seed = P$seed[i], loss_a = P$loss_a[i], loss_b = P$loss_b[i], cost = P$cost[i])
    }
  )
})
