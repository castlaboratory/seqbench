# Betting-cap sweep and declared-variance boundary (paper ii). FROZEN 2026-09-27.
# Same DGP and 12 baseline cells as main_sim_2026-09-25 (paired, one seed), 2,000 reps,
# identical streams. Methods: betting with c in {0.5 (default), 0.75, 0.9, 0.99};
# bernstein_declared with sd_max = 1x (true), 2x (over-declared) and 0.5x (violated) the
# known DGP sd. Question: how much of the gap to the Gaussian reference does the cap
# recover, and does a declared variance bound add anything beyond it.
local({
  base <- source("experiments/configs/main_sim_2026-09-25.R")$value
  dgp <- base$dgp[base$dgp$ablation == "none", ]
  list(
    name = "csweep_2026-09-27", seed = base$seed,
    alpha = base$alpha, margin = base$margin, bounds = base$bounds, n_max = base$n_max, R = 2000L,
    methods = c("seqbench:betting:c0.75", "seqbench:betting:c0.9", "seqbench:betting:c0.99",
                "seqbench:bernstein_declared:sd1", "seqbench:bernstein_declared:sd2", "seqbench:bernstein_declared:sd0.5"),
    dgp = dgp, mu_of = base$mu_of, truth_of = base$truth_of, savi_sd_of = base$savi_sd_of, gen = base$gen
  )
})
