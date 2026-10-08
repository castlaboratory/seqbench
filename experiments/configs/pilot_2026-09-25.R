# Pilot configuration (week 8, PLANO.md §5). Frozen on 2026-09-25.
# Sourced by experiments/run/run_pilot.R. Every quantity the pilot depends on is here.
list(
  name        = "pilot_2026-09-25",
  seed        = 20260925,
  alpha       = 0.05,
  margin      = 0.02,               # delta, loss units
  bounds      = c(0, 1),
  n_max       = 2000L,              # budget expressed in instances
  R           = 200L,               # replications per cell (pilot; main study sized from this)
  boundaries  = c("betting", "empirical_bernstein", "hoeffding", "naive_fixed"),
  # Data-generating process: loss_a ~ U(0.4, 0.8); loss_b = loss_a - mu + e,
  # e ~ U(-h, h), so that mu = E[loss_a - loss_b] exactly (mu > margin => B better,
  # mu < -margin => A better), sd(D) = h / sqrt(3), and no clipping is needed for
  # mu <= 0.2, h <= 0.2 (loss_b in [0.027, 0.973]).
  mu          = c(0, 0.01, 0.02, 0.05, 0.10, 0.20),   # 0, equivalent, at margin, moderate, large, very large
  half_width  = c(0.052, 0.173),    # sd(D) ~ 0.03 and 0.10
  seeds_per_instance = 1L
)
