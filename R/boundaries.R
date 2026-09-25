# Boundary kernels -----------------------------------------------------------
#
# Each boundary is a small state object with three operations:
#
#   boundary_init(name, alpha, ...)   create the state
#   boundary_update(state, x)         absorb one rescaled observation x in [0, 1]
#   boundary_interval(state)          current interval C_t on the [0, 1] scale
#
# The comparison layer (R/comparison.R) handles the loss-scale mapping, the
# running intersection and the decision rules. The kernels are exported with
# `@keywords internal` so that sister packages can reuse them without copying.
#
# Formulas: Waudby-Smith & Ramdas (2024), "Estimating means of bounded random
# variables by betting", JRSS-B 86(1):1-27; predictable-plug-in Hoeffding
# (Prop. 1), predictable-plug-in empirical Bernstein (Thm 2) and the hedged
# capital process (Thm 3). "naive_fixed" is a fixed-sample t interval recomputed
# at every step; it is NOT a confidence sequence and is flagged `valid = FALSE`.

#' Names of the available boundaries
#'
#' @return A character vector. `"betting"` is the default in
#'   [comparison_design()]; `"hoeffding"` and `"empirical_bernstein"` are
#'   conservative references; `"naive_fixed"` is an invalid negative control for
#'   experiments.
#' @export
#' @examples
#' seqbench_boundaries()
seqbench_boundaries <- function() {
  c("betting", "empirical_bernstein", "hoeffding", "naive_fixed")
}

#' Boundary kernels (internal API)
#'
#' Low-level state machines behind the confidence sequences used by
#' [update_comparison()]. Observations must already be rescaled to `[0, 1]`.
#' These functions are exported for reuse by sister packages; ordinary users
#' should call [comparison_design()] and [update_comparison()] instead.
#'
#' @param name One of [seqbench_boundaries()].
#' @param alpha Miscoverage level in (0, 1).
#' @param c Truncation constant for the empirical-Bernstein and betting
#'   bets, in (0, 1). Waudby-Smith & Ramdas recommend 1/2 or 3/4.
#' @param theta Hedging weight on the "long" capital in the betting boundary,
#'   in (0, 1). Default 1/2.
#' @param grid Number of grid points on `[0, 1]` for the betting boundary.
#' @param refine Logical; refine the betting interval endpoints by bisection on
#'   the exact capital (default `TRUE`).
#' @param state A boundary state returned by `boundary_init()` or
#'   `boundary_update()`.
#' @param x A single number in `[0, 1]`.
#'
#' @return `boundary_init()` and `boundary_update()` return a state object of
#'   class `seqbench_boundary`. `boundary_interval()` returns a named numeric
#'   vector `c(estimate, lower, upper)` on the `[0, 1]` scale; `lower > upper`
#'   (with `NA`) signals an empty set.
#' @keywords internal
#' @name boundary_kernels
#' @examples
#' st <- boundary_init("hoeffding", alpha = 0.05)
#' for (x in c(0.2, 0.3, 0.25)) st <- boundary_update(st, x)
#' boundary_interval(st)
NULL

#' @rdname boundary_kernels
#' @export
boundary_init <- function(name = seqbench_boundaries(), alpha = 0.05, c = 0.5,
                          theta = 0.5, grid = 1001L, refine = TRUE) {
  name <- rlang::arg_match(name)
  check_prob(alpha)
  check_prob(c)
  check_prob(theta)
  if (!is.numeric(grid) || length(grid) != 1L || is.na(grid) || grid < 11) {
    cli::cli_abort("{.arg grid} must be a single integer >= 11, not {.val {grid}}.")
  }
  base <- list(name = name, alpha = alpha, t = 0L, sum_x = 0, sum_sq_dev = 0,
               valid = name != "naive_fixed")
  st <- switch(name,
    hoeffding = c(base, list(sum_lam = 0, sum_lam_x = 0, sum_psi = 0)),
    empirical_bernstein = c(base, list(c = c, sum_lam = 0, sum_lam_x = 0, sum_psi = 0)),
    betting = c(base, list(
      c = c, theta = theta, refine = isTRUE(refine),
      m = seq(0, 1, length.out = as.integer(grid)),
      log_k_plus = numeric(grid), log_k_minus = numeric(grid),
      x_hist = numeric(0), lam_hist = numeric(0))),
    naive_fixed = c(base, list(sum_x2 = 0))
  )
  structure(st, class = c(paste0("seqbench_boundary_", name), "seqbench_boundary"))
}

# Running mu_hat_{t-1} and sigma2_hat_{t-1} from the sufficient statistics
# (WSR 2024, eq. 15): mu_hat_t = (1/2 + sum x) / (t + 1),
# sigma2_hat_t = (1/4 + sum (x_i - mu_hat_i)^2) / (t + 1).
running_prev <- function(st) {
  list(mu = (0.5 + st$sum_x) / (st$t + 1), sig2 = (0.25 + st$sum_sq_dev) / (st$t + 1))
}

# Update the shared sufficient statistics after observing x at time t (already
# incremented).
push_running <- function(st, x) {
  st$sum_x <- st$sum_x + x
  mu_now <- (0.5 + st$sum_x) / (st$t + 1)
  st$sum_sq_dev <- st$sum_sq_dev + (x - mu_now)^2
  st
}

#' @rdname boundary_kernels
#' @export
boundary_update <- function(state, x) {
  if (!inherits(state, "seqbench_boundary")) {
    cli::cli_abort("{.arg state} must be a boundary state from {.fn boundary_init}.")
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    cli::cli_abort("{.arg x} must be a single number in [0, 1], not {.val {x}}.")
  }
  UseMethod("boundary_update")
}

#' @export
boundary_update.seqbench_boundary_hoeffding <- function(state, x) {
  st <- state
  st$t <- st$t + 1L
  t <- st$t
  lam <- min(sqrt(8 * log(2 / st$alpha) / (t * log(t + 1))), 1)
  st$sum_lam <- st$sum_lam + lam
  st$sum_lam_x <- st$sum_lam_x + lam * x
  st$sum_psi <- st$sum_psi + lam^2 / 8
  push_running(st, x)
}

#' @export
boundary_update.seqbench_boundary_empirical_bernstein <- function(state, x) {
  st <- state
  prev <- running_prev(st)
  st$t <- st$t + 1L
  t <- st$t
  lam <- min(sqrt(2 * log(2 / st$alpha) / (prev$sig2 * t * log(1 + t))), st$c)
  v <- 4 * (x - prev$mu)^2
  psi <- (-log1p(-lam) - lam) / 4
  st$sum_lam <- st$sum_lam + lam
  st$sum_lam_x <- st$sum_lam_x + lam * x
  st$sum_psi <- st$sum_psi + v * psi
  push_running(st, x)
}

#' @export
boundary_update.seqbench_boundary_betting <- function(state, x) {
  st <- state
  prev <- running_prev(st)
  st$t <- st$t + 1L
  t <- st$t
  lam <- sqrt(2 * log(2 / st$alpha) / (prev$sig2 * t * log(t + 1)))
  st$x_hist <- c(st$x_hist, x)
  st$lam_hist <- c(st$lam_hist, lam)
  inc <- betting_log_increment(x, lam, st$m, st$c)
  st$log_k_plus <- st$log_k_plus + inc$plus
  st$log_k_minus <- st$log_k_minus + inc$minus
  push_running(st, x)
}

#' @export
boundary_update.seqbench_boundary_naive_fixed <- function(state, x) {
  st <- state
  st$t <- st$t + 1L
  st$sum_x2 <- st$sum_x2 + x^2
  push_running(st, x)
}

# One-step log capital increments for a vector of candidate means m.
betting_log_increment <- function(x, lam, m, c) {
  lam_plus <- pmin(lam, c / m)            # c / 0 = Inf -> lam
  lam_minus <- pmin(lam, c / (1 - m))
  dev <- x - m
  list(plus = log1p(lam_plus * dev), minus = log1p(-lam_minus * dev))
}

# Exact log hedged capital log K_t^{+-}(m) for scalar or vector m, from history.
betting_log_capital <- function(m, x_hist, lam_hist, c, theta) {
  vapply(m, function(mm) {
    lp <- pmin(lam_hist, c / mm)
    lm <- pmin(lam_hist, c / (1 - mm))
    dev <- x_hist - mm
    max(log(theta) + sum(log1p(lp * dev)), log(1 - theta) + sum(log1p(-lm * dev)))
  }, numeric(1))
}

#' @rdname boundary_kernels
#' @export
boundary_interval <- function(state) {
  if (!inherits(state, "seqbench_boundary")) {
    cli::cli_abort("{.arg state} must be a boundary state from {.fn boundary_init}.")
  }
  UseMethod("boundary_interval")
}

plug_in_interval <- function(st) {
  if (st$t == 0L) return(c(estimate = NA_real_, lower = 0, upper = 1))
  center <- st$sum_lam_x / st$sum_lam
  half <- (log(2 / st$alpha) + st$sum_psi) / st$sum_lam
  c(estimate = center, lower = max(center - half, 0), upper = min(center + half, 1))
}

#' @export
boundary_interval.seqbench_boundary_hoeffding <- function(state) plug_in_interval(state)

#' @export
boundary_interval.seqbench_boundary_empirical_bernstein <- function(state) plug_in_interval(state)

#' @export
boundary_interval.seqbench_boundary_betting <- function(state) {
  st <- state
  if (st$t == 0L) return(c(estimate = NA_real_, lower = 0, upper = 1))
  log_k <- pmax(log(st$theta) + st$log_k_plus, log(1 - st$theta) + st$log_k_minus)
  thr <- log(1 / st$alpha)
  inside <- which(log_k < thr)
  est <- st$m[which.min(log_k)]
  if (length(inside) == 0L) {
    return(c(estimate = est, lower = NA_real_, upper = NA_real_))
  }
  first <- inside[1]
  last <- inside[length(inside)]
  # Outer cell boundaries: the set B_t is an interval (WSR Thm 3) contained in
  # [m[first - 1], m[last + 1]]; reporting them keeps the coverage guarantee.
  lower <- if (first > 1L) st$m[first - 1L] else 0
  upper <- if (last < length(st$m)) st$m[last + 1L] else 1
  if (isTRUE(st$refine)) {
    f <- function(mm) betting_log_capital(mm, st$x_hist, st$lam_hist, st$c, st$theta) - thr
    if (first > 1L) lower <- bisect_root(f, lower, st$m[first], want = "lower")
    if (last < length(st$m)) upper <- bisect_root(f, st$m[last], upper, want = "upper")
  }
  c(estimate = est, lower = lower, upper = upper)
}

# Bisection for the crossing of f (log capital - threshold) between a and b,
# where f(a) >= 0 > f(b) ("lower") or f(a) < 0 <= f(b) ("upper"). Returns a point
# on the OUTSIDE of the crossing so that the reported interval still contains B_t.
bisect_root <- function(f, a, b, want, iter = 24L) {
  lo <- a; hi <- b
  for (i in seq_len(iter)) {
    mid <- (lo + hi) / 2
    if (want == "lower") { if (f(mid) < 0) hi <- mid else lo <- mid }
    else                 { if (f(mid) < 0) lo <- mid else hi <- mid }
    if (hi - lo < 1e-9) break
  }
  if (want == "lower") lo else hi
}

#' @export
boundary_interval.seqbench_boundary_naive_fixed <- function(state) {
  st <- state
  if (st$t < 2L) return(c(estimate = if (st$t) st$sum_x else NA_real_, lower = 0, upper = 1))
  t <- st$t
  m <- st$sum_x / t
  s2 <- max((st$sum_x2 - t * m^2) / (t - 1), 0)
  half <- stats::qt(1 - st$alpha / 2, df = t - 1) * sqrt(s2 / t)
  c(estimate = m, lower = max(m - half, 0), upper = min(m + half, 1))
}

#' @export
print.seqbench_boundary <- function(x, ...) {
  ci <- boundary_interval(x)
  cli::cli_text("{.cls seqbench_boundary} {.field {x$name}} (alpha = {x$alpha}, t = {x$t}{if (!x$valid) ', INVALID: not a confidence sequence' else ''})")
  cli::cli_text("estimate {round(ci[['estimate']], 4)}, interval [{round(ci[['lower']], 4)}, {round(ci[['upper']], 4)}]")
  invisible(x)
}

# Deterministic half-width of the PrPl-Hoeffding boundary at times 1..n
# (Prop. 2 of paper/theory.md). Used by planning_horizon().
hoeffding_half_width <- function(n, alpha) {
  t <- seq_len(n)
  lam <- pmin(sqrt(8 * log(2 / alpha) / (t * log(t + 1))), 1)
  (log(2 / alpha) + cumsum(lam^2 / 8)) / cumsum(lam)
}
