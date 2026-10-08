# Cluster ablation, version 2 (amendment 2026-09-26 after external review). FROZEN.
# The 2026-09-25 cluster cells had no instance-specific *differential* effect, so seed
# averaging faced only independent seed noise and the shared instance effect cancelled
# exactly. Here loss_b - loss_a carries an instance-specific term v_i ~ U(-w, w) shared by
# all seeds of instance i, so within-instance differences are positively correlated and
# row-level analysis (pseudoreplication) is genuinely invalid.
# DGP: u_i ~ U(0.45, 0.70); loss_a = u_i + eta_is; loss_b = u_i + eta'_is - mu - v_i + e_is,
# eta ~ U(-h_seed, h_seed), e ~ U(-h, h), v_i ~ U(-w, w). E[loss_a - loss_b] = mu exactly;
# sd of the per-instance seed-averaged difference = sqrt(w^2/3 + (h^2 + 2 h_seed^2)/(3 m)).
# All losses stay within [0.017, 0.933] for mu <= 0.2, h <= 0.173, h_seed <= 0.03, w <= 0.1.
# Design: seed count m in {1, 5} with the SAME seed noise, w in {0, 0.05}, mu in {0, 0.01,
# 0.05}; methods include the invalid pseudoreplication analysis. Cost = 1 per (instance,
# seed) pair, so `cost` is the number of evaluations of both algorithms.
local({
  dgp <- expand.grid(mu = c(0, 0.01, 0.05), h = 0.173, h_seed = 0.03, seeds = c(1L, 5L), w = c(0, 0.05), paired = TRUE)
  dgp$ablation <- "cluster_v2"
  list(
    name = "main_sim_cluster_2026-09-26", seed = 20260926L,
    alpha = 0.05, margin = 0.02, bounds = c(0, 1), n_max = 2000L, R = 2000L,
    methods = c("seqbench:betting", "pseudorep:betting", "seqbench:empirical_bernstein", "fixed:final"),
    dgp = dgp,
    mu_of = function(p) p$mu,
    truth_of = function(p, margin) if (abs(p$mu) <= margin) "equivalent" else if (p$mu > 0) "B" else "A",
    savi_sd_of = function(p) sqrt(p$w^2 / 3 + (p$h^2 + 2 * p$h_seed^2) / (3 * p$seeds)),
    gen = function(p, n) {
      m <- p$seeds; inst <- rep(seq_len(n), each = m)
      u <- runif(n, 0.45, 0.70)[inst]; v <- runif(n, -p$w, p$w)[inst]
      la <- u + runif(n * m, -p$h_seed, p$h_seed)
      lb <- u + runif(n * m, -p$h_seed, p$h_seed) - p$mu - v + runif(n * m, -p$h, p$h)
      data.frame(instance = inst, seed = rep(seq_len(m), times = n), loss_a = la, loss_b = lb, cost = 1)
    }
  )
})
