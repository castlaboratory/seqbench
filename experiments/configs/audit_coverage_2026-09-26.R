# No-stopping coverage audit (amendment 2026-09-26 after external review). FROZEN.
# The main study records coverage only up to the stopping time. This audit runs each
# boundary over the full stream of n_max instances on the 24 DGP cells of the main study
# (same seeds, hence the same streams) and records whether the true mean ever leaves the
# running confidence sequence. This is the time-uniform coverage claim, measured directly.
local({
  base <- source("experiments/configs/main_sim_2026-09-25.R")$value
  list(
    name = "audit_coverage_2026-09-26", seed = base$seed,
    alpha = base$alpha, margin = base$margin, bounds = base$bounds, n_max = base$n_max, R = 2000L,
    methods = c("audit:betting", "audit:empirical_bernstein", "audit:hoeffding"),
    dgp = base$dgp, mu_of = base$mu_of, truth_of = base$truth_of, gen = base$gen
  )
})
