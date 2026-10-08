# Application (a): Louvain vs Leiden (igraph) for modularity maximisation. FROZEN 2026-09-25.
# Both algorithms are on CRAN (igraph) and randomised, so the seed matters and the paired
# design (same graph, same seed) is meaningful. Question: is Leiden relevantly better than
# Louvain in modularity on this population, with margin 0.01 (one modularity point)?
# Population P: stochastic block models with k ~ {3..6} blocks of unequal sizes, n ~ U{80..300},
# p_in ~ U(0.15, 0.40), p_out = p_in * mix with mix ~ U(0.05, 0.40) drawn per instance.
# A = cluster_louvain, B = cluster_leiden(objective = "modularity", n_iterations = -1).
# Loss = 1 - Q, bounds [0, 1.5] since Q in [-1/2, 1] a priori (never clipped). Cost = seconds.
# Pilot (40 instances, 2026-09-25): mean(loss_a - loss_b) ~ 0.007, sd ~ 0.011, ~1 ms per run,
# so with margin 0.01 the truth is "equivalent" at distance ~0.003 and the truncated-bet regime
# predicts ~3,700 instances; n_max = 6000 leaves room (instance generation ~45 ms). The truth is not known analytically; a
# Monte Carlo reference (5000 instances) is computed once by mu_of() in the runner's parent process.
local({
  gen_instance <- function(seed) {
    set.seed(seed)
    k <- sample(3:6, 1); n <- sample(80:300, 1)
    sizes <- as.vector(stats::rmultinom(1, n, prob = runif(k, 0.5, 1.5))); sizes[sizes < 5] <- 5
    p_in <- runif(1, 0.15, 0.4); p_out <- p_in * runif(1, 0.05, 0.4)
    pm <- matrix(p_out, k, k); diag(pm) <- p_in
    g <- igraph::simplify(igraph::sample_sbm(sum(sizes), pref.matrix = pm, block.sizes = sizes))
    set.seed(seed); ta <- system.time(qa <- igraph::modularity(g, igraph::membership(igraph::cluster_louvain(g))))[["elapsed"]]
    set.seed(seed); tb <- system.time(qb <- igraph::modularity(g, igraph::membership(
      igraph::cluster_leiden(g, objective_function = "modularity", n_iterations = -1))))[["elapsed"]]
    c(loss_a = 1 - qa, loss_b = 1 - qb, cost = ta + tb)
  }
  list(
    name = "app_louvain_leiden_2026-09-25", seed = 20260925L,
    alpha = 0.05, margin = 0.01, bounds = c(0, 1.5), n_max = 6000L, R = 100L, cache_data = TRUE, savi_pilot_n = 200L,   # savi design from an independent 200-instance pilot (amendment 2026-09-26)
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "fixed:final", "savi_cs"),
    dgp = data.frame(scenario = "louvain_vs_leiden_sbm"),
    mu_of = local({ cache <- NULL; function(p) {
      if (is.null(cache)) { set.seed(1); s <- sample.int(1e9, 5000)
        cache <<- mean(vapply(s, function(z) { l <- gen_instance(z); l[["loss_a"]] - l[["loss_b"]] }, numeric(1))) }
      cache } }),
    truth_of = function(p, margin) NA_character_,
    gen = function(p, n) {
      seeds <- sample.int(1e9, n)
      l <- t(vapply(seeds, gen_instance, numeric(3)))
      data.frame(instance = seq_len(n), seed = seeds, loss_a = l[, 1], loss_b = l[, 2], cost = l[, 3])
    }
  )
})
