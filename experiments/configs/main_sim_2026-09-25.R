# Main simulation study (PLANO.md §4). FROZEN 2026-09-25; changes require a new file.
# Scenarios: null effect, true equivalence, effect at the margin, moderate/large effects;
# two noise levels; ablations without pairing and with seed clusters. References:
# fixed-sample tests (no adaptation) and safestats' anytime-valid t-test CS.
#
# DGP. Instance effect u_i ~ U(0.45, 0.70). loss_a(i,s) = u_i + eta_is,
# loss_b(i,s) = (paired ? u_i : u'_i) + eta'_is - mu + e_is, eta ~ U(-h_seed, h_seed),
# e ~ U(-h, h). Then E[loss_a - loss_b] = mu exactly, sd of the per-instance seed-averaged
# difference = sqrt((h^2 + 2 h_seed^2)/(3 m) + (unpaired ? 2 * 0.25^2/12 : 0)), and all
# losses stay in [0.017, 0.933] for mu <= 0.2, h <= 0.173, h_seed <= 0.03.
local({
  base <- expand.grid(mu = c(0, 0.01, 0.02, 0.05, 0.10, 0.20), h = c(0.052, 0.173),
                      paired = TRUE, seeds = 1L, h_seed = 0)
  unpaired <- expand.grid(mu = c(0, 0.01, 0.02, 0.05, 0.10, 0.20), h = 0.173, paired = FALSE, seeds = 1L, h_seed = 0)
  cluster <- expand.grid(mu = c(0, 0.01, 0.02, 0.05, 0.10, 0.20), h = 0.173, paired = TRUE, seeds = 5L, h_seed = 0.03)
  dgp <- rbind(base, unpaired, cluster)
  dgp$ablation <- rep(c("none", "unpaired", "cluster"), c(nrow(base), nrow(unpaired), nrow(cluster)))
  list(
    name = "main_sim_2026-09-25", seed = 20260925L,
    alpha = 0.05, margin = 0.02, bounds = c(0, 1), n_max = 2000L, R = 2000L,
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "seqbench:hoeffding",
                "seqbench:naive_fixed", "fixed:100", "fixed:500", "fixed:final", "savi_cs"),
    dgp = dgp,
    mu_of = function(p) p$mu,
    truth_of = function(p, margin) if (abs(p$mu) <= margin) "equivalent" else if (p$mu > 0) "B" else "A",
    gen = function(p, n) {
      m <- p$seeds; inst <- rep(seq_len(n), each = m)
      u <- runif(n, 0.45, 0.70)[inst]
      u2 <- if (p$paired) u else runif(n, 0.45, 0.70)[inst]
      la <- u + runif(n * m, -p$h_seed, p$h_seed)
      lb <- u2 + runif(n * m, -p$h_seed, p$h_seed) - p$mu + runif(n * m, -p$h, p$h)
      data.frame(instance = inst, seed = rep(seq_len(m), times = n), loss_a = la, loss_b = lb, cost = 1)
    },
    # AMENDMENT 2026-09-26 (after external review, before any savi_cs result was used): the savi
    # design is fixed from the KNOWN DGP sd of the per-instance difference (oracle, favourable
    # to safestats), never from the evaluation stream. The 24 savi_cs cells were rerun.
    savi_sd_of = function(p) sqrt((p$h^2 + 2 * p$h_seed^2) / (3 * p$seeds) + ifelse(p$paired, 0, 2 * 0.25^2 / 12))
  )
})
