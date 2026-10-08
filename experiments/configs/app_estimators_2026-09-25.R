# Application (b): two statistical estimators compared by simulation. FROZEN 2026-09-25.
# Population P of instances: samples of size n_obs from a contaminated normal
# theta + scale * (1 - w) N(0,1) + w N(0, 5^2), with n_obs ~ {20, 40, 80}, scale ~ U(1, 3),
# contamination w ~ U(0, 0.25) drawn per instance. A = sample mean, B = 20% trimmed mean.
# Loss = min(|estimate - theta| / scale, 1) in [0, 1] (scaled absolute error, capped).
# Same sample for both (paired). The truth mu is not known analytically; a Monte Carlo
# reference (1e5 instances) is computed once by mu_of() for coverage bookkeeping.
local({
  gen_instance <- function() {
    n_obs <- sample(c(20L, 40L, 80L), 1); scale <- runif(1, 1, 3); w <- runif(1, 0, 0.25)
    x <- scale * ifelse(runif(n_obs) < w, rnorm(n_obs, 0, 5), rnorm(n_obs))
    c(min(abs(mean(x)) / scale, 1), min(abs(mean(x, trim = 0.2)) / scale, 1))
  }
  list(
    name = "app_estimators_2026-09-25", seed = 20260925L,
    alpha = 0.05, margin = 0.02, bounds = c(0, 1), n_max = 3000L, R = 500L, cache_data = TRUE, savi_pilot_n = 200L,   # savi design from an independent 200-instance pilot (amendment 2026-09-26)
    methods = c("seqbench:betting", "seqbench:empirical_bernstein", "fixed:final", "savi_cs"),
    dgp = data.frame(scenario = "mean_vs_trimmed"),
    mu_of = local({ cache <- NULL; function(p) {
      if (is.null(cache)) { set.seed(1); cache <<- mean(replicate(1e5, { l <- gen_instance(); l[1] - l[2] })) }
      cache } }),
    truth_of = function(p, margin) NA_character_,   # real application: truth unknown
    gen = function(p, n) {
      l <- t(replicate(n, gen_instance()))
      data.frame(instance = seq_len(n), seed = 1L, loss_a = l[, 1], loss_b = l[, 2], cost = 1)
    }
  )
})
