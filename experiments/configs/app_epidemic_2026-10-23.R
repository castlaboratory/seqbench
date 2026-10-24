# Application (epidemic control): two vaccination strategies on contact networks.
# FROZEN 2026-10-23, before any loss was computed.
#
# Question. With the same number of doses, does vaccinating the highest-degree people
# (degree targeting, which needs the full contact network) reduce the final attack rate of
# an SIR epidemic by a practically relevant amount relative to acquaintance immunisation
# (Cohen, Havlin & ben-Avraham 2003: vaccinate a random contact of a random person, which
# needs only local information)? Margin delta = 0.01: one percentage point of the population.
#
# Population P (one draw per instance): N ~ U{1000..3000} people; contact network from one
# of five models with probability 1/5 each, mean degree ~ U(4, 12) where the model allows:
# Erdos-Renyi G(N, M); Barabasi-Albert (m ~ U{2..5}); Watts-Strogatz small world (k/2 ~
# U{2..5}, rewiring ~ U(0.01, 0.2)); stochastic block model with 5-20 communities and 90%
# of contacts inside them; Chung-Lu power law (exponent ~ U(2.1, 3)). Basic reproduction
# number R0 ~ U(1.5, 3); vaccine coverage v ~ U(0.05, 0.20) of N.
# Epidemic: SIR with a fixed infectious period of one generation and per-contact
# transmissibility T = min(1, R0 <k> / (<k^2> - <k>)), five index cases among the
# unvaccinated. Its final size is exactly the size of the union of the clusters of the
# index cases in bond percolation with occupation probability T on the network without the
# vaccinated (Newman 2002; Kenah & Robins 2007), which is how it is simulated.
# Seeds: every instance has m = 5 seeds (independent epidemics: vaccination draw for A,
# tie-breaking for B, percolation, index cases); within a seed both strategies share the
# percolated network and the random order from which index cases are taken (common random
# numbers), so the pairing is by (instance, seed) and seeds of an instance are a cluster.
# Loss = final attack rate (infected / N) in [0, 1] by construction.
# Cost = CPU seconds of the 2 x 5 simulations (network generation excluded).
# Cells: the five-seed design and the one-seed design (first seed only), and the INVALID
# pseudoreplication analysis that feeds every (instance, seed) row as its own instance.
# Operating characteristics on streams resampled from a frozen pool (see app_cvrp).
local({
  network <- function(N) {
    kbar <- runif(1, 4, 12); mi <- sample(5, 1)
    g <- switch(mi,
      igraph::sample_gnm(N, round(N * kbar / 2)),
      igraph::sample_pa(N, m = sample(2:5, 1), directed = FALSE),
      igraph::sample_smallworld(1, N, nei = sample(2:5, 1), p = runif(1, 0.01, 0.2)),
      { k <- sample(5:20, 1); sizes <- as.vector(stats::rmultinom(1, N, runif(k, 0.5, 1.5))); sizes[sizes < 2] <- 2
        n <- sum(sizes); p_in <- 0.9 * kbar / (sum(sizes^2) / n); p_out <- 0.1 * kbar / (n - sum(sizes^2) / n)
        pm <- matrix(min(p_out, 1), k, k); diag(pm) <- min(p_in, 1)
        igraph::sample_sbm(n, pm, sizes) },
      igraph::sample_fitness_pl(N, round(N * kbar / 2), exponent.out = runif(1, 2.1, 3)))
    g <- igraph::simplify(g)
    igraph::graph_attr(g, "model") <- c("er", "ba", "ws", "sbm", "powerlaw")[mi]
    g
  }
  final_size <- function(gp, vacc, ord, n_index = 5L) {
    keep <- !vacc
    sub <- igraph::induced_subgraph(gp, which(keep))
    memb <- igraph::components(sub)$membership
    idx_old <- ord[keep[ord]][seq_len(n_index)]       # first unvaccinated people in the shared order
    idx_new <- match(idx_old, which(keep))
    cl <- unique(memb[idx_new])
    sum(memb %in% cl) / length(vacc)
  }
  cpu <- function() sum(proc.time()[c("user.self", "sys.self")])   # system.time() would gc() each call
  acquaintance <- function(adj, nd) {           # random contact of a random person, until nd distinct
    N <- length(adj); v <- logical(N); deg <- lengths(adj)
    while (sum(v) < nd) {
      r <- sample.int(N, 2L * nd, replace = TRUE); r <- r[deg[r] > 0]
      pick <- vapply(r, function(i) adj[[i]][sample.int(deg[i], 1L)], integer(1))
      for (u in pick) { if (!v[u]) { v[u] <- TRUE; if (sum(v) >= nd) break } }
    }
    v
  }
  instance <- function(seed, m = 5L) {
    set.seed(seed)
    N <- sample(1000:3000, 1)
    g <- network(N); N <- igraph::vcount(g)
    deg <- igraph::degree(g); k1 <- mean(deg); k2 <- mean(deg^2)
    R0 <- runif(1, 1.5, 3); Tr <- min(1, R0 * k1 / (k2 - k1)); nd <- round(runif(1, 0.05, 0.20) * N)
    el <- igraph::as_edgelist(g, names = FALSE)
    adj <- split(c(el[, 2], el[, 1]), factor(c(el[, 1], el[, 2]), levels = seq_len(N)))
    adj <- lapply(adj, as.integer)
    out <- vector("list", m)
    for (s in seq_len(m)) {
      c0 <- cpu()
      va <- acquaintance(adj, nd)
      vb <- logical(N); vb[order(-deg, runif(N))[seq_len(nd)]] <- TRUE
      gp <- igraph::delete_edges(g, which(runif(igraph::ecount(g)) > Tr))
      ord <- sample.int(N)
      la <- final_size(gp, va, ord); lb <- final_size(gp, vb, ord)
      out[[s]] <- data.frame(seed = s, loss_a = la, loss_b = lb, cost = cpu() - c0)
    }
    cbind(do.call(rbind, out), N = N, model = igraph::graph_attr(g, "model"), kbar = k1,
          R0 = R0, T = Tr, coverage = nd / N, row.names = NULL)
  }
  pool_file <- "experiments/out/pools/app_epidemic_2026-10-23.rds"
  pool <- local({ cache <- NULL; function() { if (is.null(cache)) cache <<- readRDS(pool_file); cache } })
  pool_seeds <- function(m) { P <- pool(); P[P$seed <= m, ] }
  # mean of the per-instance seed-averaged difference over the pool, for m seeds
  pool_mu <- function(m) { P <- pool_seeds(m); mean(tapply(P$loss_a - P$loss_b, P$instance, mean)) }
  list(
    name = "app_epidemic_2026-10-23", seed = 20261024L,
    alpha = 0.05, margin = 0.01, bounds = c(0, 1), n_max = 3000L, R = 1000L,
    pool_n = 10000L, pool_file = pool_file, instance = instance, savi_pilot_n = 200L,
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "fixed:100", "fixed:final", "savi_cs", "pseudorep:betting"),
    dgp = data.frame(scenario = "acquaintance_vs_degree", seeds = c(5L, 1L)),
    mu_of = function(p) pool_mu(p$seeds),
    truth_of = function(p, margin) {
      mu <- pool_mu(p$seeds)
      if (abs(mu) <= margin) "equivalent" else if (mu > 0) "B" else "A"
    },
    gen = function(p, n) {                            # i.i.d. resampling of whole instances
      P <- pool_seeds(p$seeds); ids <- unique(P$instance)
      pick <- ids[sample.int(length(ids), n, replace = TRUE)]
      rows <- split(seq_len(nrow(P)), P$instance)[as.character(pick)]
      i <- unlist(rows, use.names = FALSE)
      data.frame(instance = rep(seq_len(n), lengths(rows)), seed = P$seed[i],
                 loss_a = P$loss_a[i], loss_b = P$loss_b[i], cost = P$cost[i])
    }
  )
})
