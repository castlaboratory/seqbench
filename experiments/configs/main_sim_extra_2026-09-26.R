# Extra scenarios (amendment 2026-09-26 after external review). FROZEN.
# Same DGP and methods as main_sim_2026-09-25 (baseline, paired, one seed), adding
# negative effects (A better) and effects just inside/outside the margin, which the
# 2026-09-25 grid skipped (0.02 -> 0.05).
local({
  dgp <- expand.grid(mu = c(-0.20, -0.05, -0.025, -0.015, 0.015, 0.025, 0.03), h = c(0.052, 0.173),
                     paired = TRUE, seeds = 1L, h_seed = 0)
  dgp$ablation <- "extra"
  base <- source("experiments/configs/main_sim_2026-09-25.R")$value
  list(
    name = "main_sim_extra_2026-09-26", seed = 20260926L,
    alpha = base$alpha, margin = base$margin, bounds = base$bounds, n_max = base$n_max, R = 2000L,
    methods = base$methods, dgp = dgp, mu_of = base$mu_of, truth_of = base$truth_of,
    savi_sd_of = base$savi_sd_of, gen = base$gen
  )
})
