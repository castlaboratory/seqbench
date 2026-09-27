# Analytic and simulation checks of the boundary kernels against the formulas
# of Waudby-Smith & Ramdas (2024) as transcribed in references/wsr2024_fronteiras.md.

run_kernel <- function(name, x, alpha = 0.05, ...) {
  st <- boundary_init(name, alpha = alpha, ...)
  out <- matrix(NA_real_, length(x), 3, dimnames = list(NULL, c("estimate", "lower", "upper")))
  for (i in seq_along(x)) {
    st <- boundary_update(st, x[i])
    out[i, ] <- boundary_interval(st)
  }
  list(state = st, path = out)
}

# Independent vectorised re-derivations used as oracles
oracle_hoeffding <- function(x, alpha) {
  t <- seq_along(x)
  lam <- pmin(sqrt(8 * log(2 / alpha) / (t * log(t + 1))), 1)
  S <- cumsum(lam); center <- cumsum(lam * x) / S
  half <- (log(2 / alpha) + cumsum(lam^2 / 8)) / S
  cbind(estimate = center, lower = pmax(center - half, 0), upper = pmin(center + half, 1))
}
oracle_eb <- function(x, alpha, c = 0.5) {
  t <- seq_along(x)
  mu_hat <- (0.5 + cumsum(x)) / (t + 1)
  sig2 <- (0.25 + cumsum((x - mu_hat)^2)) / (t + 1)
  mu_prev <- c(0.5, mu_hat[-length(x)]); sig2_prev <- c(0.25, sig2[-length(x)])
  lam <- pmin(sqrt(2 * log(2 / alpha) / (sig2_prev * t * log(1 + t))), c)
  v <- 4 * (x - mu_prev)^2; psi <- (-log1p(-lam) - lam) / 4
  S <- cumsum(lam); center <- cumsum(lam * x) / S
  half <- (log(2 / alpha) + cumsum(v * psi)) / S
  cbind(estimate = center, lower = pmax(center - half, 0), upper = pmin(center + half, 1))
}
oracle_log_capital <- function(x, m, alpha, c = 0.5, theta = 0.5) {
  t <- seq_along(x)
  mu_hat <- (0.5 + cumsum(x)) / (t + 1)
  sig2 <- (0.25 + cumsum((x - mu_hat)^2)) / (t + 1)
  sig2_prev <- c(0.25, sig2[-length(x)])
  lt <- sqrt(2 * log(2 / alpha) / (sig2_prev * t * log(t + 1)))
  lp <- pmin(lt, c / m); lm <- pmin(lt, c / (1 - m))
  max(log(theta) + sum(log1p(lp * (x - m))), log(1 - theta) + sum(log1p(-lm * (x - m))))
}

set.seed(20260925)
x <- rbeta(200, 10, 30)

test_that("Hoeffding kernel matches the closed form at every t", {
  r <- run_kernel("hoeffding", x)
  expect_equal(r$path, oracle_hoeffding(x, 0.05), tolerance = 1e-12)
  # t = 1: lambda_1 = min(sqrt(8 log 40 / log 2), 1) = 1, so C_1 = x_1 +- (log 40 + 1/8)
  expect_equal(unname(r$path[1, "estimate"]), x[1])
  expect_equal(unname(r$path[1, c("lower", "upper")]), c(0, 1))
})

test_that("empirical-Bernstein kernel matches the closed form at every t", {
  r <- run_kernel("empirical_bernstein", x)
  expect_equal(r$path, oracle_eb(x, 0.05), tolerance = 1e-12)
})

test_that("betting interval contains the true mean iff exact capital is below 1/alpha", {
  r <- run_kernel("betting", x, grid = 401L)
  mu <- 0.25
  st <- r$state
  logk <- seqbench:::betting_log_capital(mu, st$x_hist, st$lam_hist, st$c, st$theta)
  expect_equal(logk, oracle_log_capital(x, mu, 0.05), tolerance = 1e-12)
  inside <- logk < log(1 / 0.05)
  ci <- r$path[200, ]
  expect_true(inside == (mu >= ci[["lower"]] && mu <= ci[["upper"]]))
})

test_that("betting endpoints are outside the crossing (conservative) and refined within a cell", {
  r_ref <- run_kernel("betting", x, grid = 401L, refine = TRUE)
  r_raw <- run_kernel("betting", x, grid = 401L, refine = FALSE)
  st <- r_ref$state
  thr <- log(1 / 0.05)
  ci <- r_ref$path[200, ]
  f <- function(m) seqbench:::betting_log_capital(m, st$x_hist, st$lam_hist, st$c, st$theta)
  expect_true(f(ci[["lower"]]) >= thr - 1e-6)     # just outside
  expect_true(f(ci[["upper"]]) >= thr - 1e-6)
  expect_true(f(ci[["lower"]] + 1e-6) < thr || f(ci[["lower"]] + 2 / 400) < thr)
  expect_true(ci[["lower"]] >= r_raw$path[200, "lower"] - 1e-12)   # refined within outer cell
  expect_true(ci[["upper"]] <= r_raw$path[200, "upper"] + 1e-12)
  expect_lt(ci[["upper"]] - ci[["lower"]], r_raw$path[200, "upper"] - r_raw$path[200, "lower"] + 1e-12)
})

test_that("widths are ordered hoeffding > empirical_bernstein > betting at t = 200 for Beta(10,30)", {
  w <- vapply(c("hoeffding", "empirical_bernstein", "betting"), function(b) {
    p <- run_kernel(b, x)$path[200, ]; p[["upper"]] - p[["lower"]]
  }, numeric(1))
  expect_true(w[["hoeffding"]] > w[["empirical_bernstein"]])
  expect_true(w[["empirical_bernstein"]] > w[["betting"]])
})

test_that("snapshot at t = 200 agrees with the standalone reference implementation", {
  # Values produced by dev/wsr_cs.R on 2026-09-25 with set.seed(1); x <- rbeta(200, 10, 30).
  set.seed(1); y <- rbeta(200, 10, 30)
  # dev/wsr_cs.R reports the running intersection; apply it to the raw kernel path
  intersect_last <- function(p) c(lower = max(p[, "lower"]), upper = min(p[, "upper"]))
  h <- intersect_last(run_kernel("hoeffding", y)$path)
  e <- intersect_last(run_kernel("empirical_bernstein", y)$path)
  expect_equal(unname(h[c("lower", "upper")]), c(0.1384825, 0.3680665), tolerance = 1e-6)
  expect_equal(unname(e[c("lower", "upper")]), c(0.2123603, 0.2897128), tolerance = 1e-6)
  b <- intersect_last(run_kernel("betting", y, grid = 401L)$path)
  # dev grid-inner interval was [0.2375, 0.2775]; the package reports the outer, refined one
  expect_true(b[["lower"]] <= 0.2375 && b[["lower"]] >= 0.2375 - 1 / 400)
  expect_true(b[["upper"]] >= 0.2775 && b[["upper"]] <= 0.2775 + 1 / 400)
})

test_that("time-uniform coverage holds within Monte Carlo error (small simulation)", {
  skip_on_cran()
  R <- 300L; n <- 150L; alpha <- 0.1; mu <- 0.5
  miss <- c(hoeffding = 0, empirical_bernstein = 0, betting = 0)
  set.seed(7)
  for (r in seq_len(R)) {
    z <- rbinom(n, 1, mu)                          # hardest case: max variance
    for (b in names(miss)) {
      p <- run_kernel(b, z, alpha = alpha, grid = 201L)$path
      miss[b] <- miss[b] + any(cummax(p[, "lower"]) > mu | cummin(p[, "upper"]) < mu)
    }
  }
  se <- sqrt(alpha * (1 - alpha) / R)
  for (b in names(miss)) expect_lte(miss[[b]] / R, alpha + 3 * se)
})

test_that("naive_fixed is flagged invalid and reproduces the t interval", {
  r <- run_kernel("naive_fixed", x)
  expect_false(r$state$valid)
  tt <- t.test(x)$conf.int
  expect_equal(unname(r$path[200, c("lower", "upper")]), as.numeric(tt), tolerance = 1e-10)
})

test_that("kernel argument checks", {
  expect_error(boundary_init("nope"), "must be one of")
  expect_error(boundary_init("betting", alpha = 1), "alpha")
  expect_error(boundary_update(boundary_init("hoeffding"), 1.5), "in \\[0, 1\\]")
  expect_error(boundary_update(boundary_init("hoeffding"), NA_real_), "in \\[0, 1\\]")
  expect_error(boundary_interval(list()), "boundary state")
  expect_no_error(print(boundary_init("betting")))
})

test_that("a nonempty betting set narrower than a grid cell is found, not reported empty", {
  # Reviewer counterexample: constant stream at 0.55 with an 11-point grid. The exact
  # capital at the true mean is log(1/2) < log(1/alpha), so the set contains 0.55.
  st <- boundary_init("betting", alpha = 0.05, grid = 11L)
  for (i in 1:200) st <- boundary_update(st, 0.55)
  ci <- boundary_interval(st)
  expect_false(anyNA(ci))
  expect_true(ci[["lower"]] <= 0.55 && ci[["upper"]] >= 0.55)
  expect_lt(ci[["upper"]] - ci[["lower"]], 0.1)                  # inside one cell
  lk <- seqbench:::betting_log_capital(0.55, st$x_hist, st$lam_hist, st$c, st$theta)
  expect_equal(lk, log(0.5), tolerance = 1e-10)
  # unrefined variant returns the enclosing cell edges
  st2 <- boundary_init("betting", alpha = 0.05, grid = 11L, refine = FALSE)
  for (i in 1:200) st2 <- boundary_update(st2, 0.55)
  ci2 <- boundary_interval(st2)
  expect_true(ci2[["lower"]] <= ci[["lower"]] && ci2[["upper"]] >= ci[["upper"]])
  # the same stream through the comparison layer must not stop with cs_empty
  d <- comparison_design(margin = 0.02, bounds = c(0, 1), betting_grid = 11L)
  l <- data.frame(instance = 1:200, loss_a = 0.6, loss_b = 0.5)   # D = 0.1 -> x = 0.55
  s3 <- suppressWarnings(update_comparison(initialize_comparison(d), l))
  expect_false(s3$diagnostics$cs_empty)
  expect_equal(as.vector(stopping_decision(s3)), "B")
})

test_that("kernel rejects malformed grid and refine", {
  expect_error(boundary_init("betting", grid = 100.5), "whole number")
  expect_error(boundary_init("betting", grid = Inf), "whole number")
  expect_error(boundary_init("betting", refine = NA), "refine")
})

test_that("bernstein_declared matches its closed form and shrinks with the declared sd", {
  sig <- 0.03 / 2                       # sd 0.03 in loss units on bounds [0,1] -> [0,1] scale
  r <- run_kernel("bernstein_declared", x, sd_max01 = sig)
  # oracle
  A <- log(2 / 0.05); t <- seq_along(x)
  lam <- vapply(t, function(tt) optimize(function(l) (A + sig^2 * tt * (exp(l) - 1 - l)) / (tt * l), c(1e-6, 40), tol = 1e-8)$minimum, numeric(1))
  S <- cumsum(lam); center <- cumsum(lam * x) / S
  half <- (A + sig^2 * cumsum(exp(lam) - 1 - lam)) / S
  expect_equal(unname(r$path[, "estimate"]), center, tolerance = 1e-6)
  expect_equal(unname(r$path[, "upper"] - r$path[, "lower"]), pmin(center + half, 1) - pmax(center - half, 0), tolerance = 1e-6)
  # larger declared sd -> wider
  r2 <- run_kernel("bernstein_declared", x, sd_max01 = 0.2)
  expect_true(all(r2$path[, "upper"] - r2$path[, "lower"] >= r$path[, "upper"] - r$path[, "lower"] - 1e-12))
  # with a correct small declared sd it is narrower than betting at t = 200 for this low-noise stream
  b <- run_kernel("betting", x)$path[200, ]
  expect_lt(r$path[200, "upper"] - r$path[200, "lower"], b[["upper"]] - b[["lower"]])
  expect_error(boundary_init("bernstein_declared"), "sd_max01")
  expect_error(boundary_init("bernstein_declared", sd_max01 = 0.7), "sd_max01")
})

test_that("bernstein_declared keeps time-uniform coverage when the declared sd is true", {
  skip_on_cran()
  R <- 300L; n <- 200L; alpha <- 0.1; miss <- 0
  set.seed(11)
  for (r in seq_len(R)) {
    z <- 0.5 + runif(n, -0.05, 0.05)             # sd = 0.05/sqrt(3) = 0.0289
    p <- run_kernel("bernstein_declared", z, alpha = alpha, sd_max01 = 0.0289)$path
    miss <- miss + any(cummax(p[, "lower"]) > 0.5 | cummin(p[, "upper"]) < 0.5)
  }
  expect_lte(miss / R, alpha + 3 * sqrt(alpha * (1 - alpha) / R))
})
