# dev/wsr_cs.R -- reference implementation of the three confidence sequences of
# Waudby-Smith & Ramdas (2024), transcribed from references/wsr2024_fronteiras.md.
# Standalone R, no package code. Observations x in [0, 1]. Every function returns a
# data.frame(t, estimate, lower, upper) with the RUNNING INTERSECTION already applied.
# Used by dev/repro_betting.R and dev/repro_safestats.R; later the package tests will
# compare R/boundaries.R against these numbers.

wsr_running <- function(x) {
  t <- seq_along(x)
  mu_hat <- (0.5 + cumsum(x)) / (t + 1)                  # mu_hat_t
  sig2   <- (0.25 + cumsum((x - mu_hat)^2)) / (t + 1)    # sigma2_hat_t
  list(t = t,
       mu_prev   = c(0.5,  mu_hat[-length(x)]),           # mu_hat_{t-1}
       sig2_prev = c(0.25, sig2[-length(x)]))             # sigma2_hat_{t-1}
}

wsr_intersect <- function(t, est, lo, hi) {
  data.frame(t = t, estimate = est,
             lower = pmax(cummax(lo), 0), upper = pmin(cummin(hi), 1))
}

# Predictable-plug-in Hoeffding (Prop. 1, eq. 12)
cs_prpl_h <- function(x, alpha = 0.05) {
  t <- seq_along(x)
  lam <- pmin(sqrt(8 * log(2 / alpha) / (t * log(t + 1))), 1)
  S <- cumsum(lam)
  center <- cumsum(lam * x) / S
  half <- (log(2 / alpha) + cumsum(lam^2 / 8)) / S
  wsr_intersect(t, center, center - half, center + half)
}

# Predictable-plug-in empirical Bernstein (Thm 2, eq. 15)
cs_prpl_eb <- function(x, alpha = 0.05, c = 0.5) {
  r <- wsr_running(x); t <- r$t
  lam <- pmin(sqrt(2 * log(2 / alpha) / (r$sig2_prev * t * log(1 + t))), c)
  v <- 4 * (x - r$mu_prev)^2
  psi_e <- (-log1p(-lam) - lam) / 4
  S <- cumsum(lam)
  center <- cumsum(lam * x) / S
  half <- (log(2 / alpha) + cumsum(v * psi_e)) / S
  wsr_intersect(t, center, center - half, center + half)
}

# Predictable-mixture lambda for betting (eq. 26): no truncation at 1
wsr_lambda_tilde <- function(x, alpha = 0.05) {
  r <- wsr_running(x)
  sqrt(2 * log(2 / alpha) / (r$sig2_prev * r$t * log(r$t + 1)))
}

# Log hedged capital K_t^{+-}(m) for a vector of m, all t. Returns t x length(m) matrix.
hedged_log_capital <- function(x, m, alpha = 0.05, c = 0.5, theta = 0.5) {
  lt <- wsr_lambda_tilde(x, alpha)
  lam_p <- outer(lt, c / m, pmin)          # lambda_t^+(m) = |lt| ^ c/m
  lam_m <- outer(lt, c / (1 - m), pmin)    # lambda_t^-(m) = |lt| ^ c/(1-m)
  dev <- outer(x, m, "-")                  # x_t - m
  logK_p <- apply(log1p( lam_p * dev), 2, cumsum)
  logK_m <- apply(log1p(-lam_m * dev), 2, cumsum)
  if (length(x) == 1L) { logK_p <- matrix(logK_p, 1); logK_m <- matrix(logK_m, 1) }
  pmax(log(theta) + logK_p, log(1 - theta) + logK_m)
}

# Hedged betting CS B_t^{+-} (Thm 3) on a grid of m, then running intersection.
cs_hedged <- function(x, alpha = 0.05, c = 0.5, theta = 0.5, n_grid = 1001L) {
  m <- seq(0, 1, length.out = n_grid)
  logK <- hedged_log_capital(x, m, alpha, c, theta)
  inside <- logK < log(1 / alpha)
  lo <- apply(inside, 1, function(z) if (any(z)) m[which(z)[1]] else NA_real_)
  hi <- apply(inside, 1, function(z) if (any(z)) m[max(which(z))] else NA_real_)
  est <- m[apply(logK, 1, which.min)]
  wsr_intersect(seq_along(x), est, lo, hi)
}

# Exact miscoverage event for the hedged CS at the true mean mu, without a grid:
# TRUE if K_t^{+-}(mu) >= 1/alpha for some t (equivalently mu leaves the running CS).
hedged_misses_mu <- function(x, mu, alpha = 0.05, c = 0.5, theta = 0.5) {
  any(hedged_log_capital(x, mu, alpha, c, theta) >= log(1 / alpha))
}

# Negative control: fixed-sample t interval recomputed at every t (INVALID sequentially)
cs_naive_t <- function(x, alpha = 0.05) {
  t <- seq_along(x)
  m <- cumsum(x) / t
  s2 <- c(NA, (cumsum(x^2) - t * m^2)[-1] / (t[-1] - 1))
  half <- stats::qt(1 - alpha / 2, df = t - 1) * sqrt(s2 / t)
  data.frame(t = t, estimate = m, lower = m - half, upper = m + half)  # no intersection: invalid anyway
}
