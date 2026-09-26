#' Design a sequential comparison of two algorithms
#'
#' Fixes, before any data is seen, everything the procedure needs: the error
#' level, the margin of practical equivalence, the loss bounds, the boundary and
#' the budget. The design is immutable; [initialize_comparison()] turns it into
#' a state that [update_comparison()] advances.
#'
#' @param alpha Miscoverage level of the confidence sequence, in (0, 1). The
#'   probability that the procedure ever issues a false declaration (of any of
#'   the three kinds) is at most `alpha`.
#' @param margin Margin of practical equivalence `delta > 0`, in loss units.
#'   Algorithm A is declared relevantly better when the whole confidence
#'   sequence for `mean(loss_a - loss_b)` lies below `-margin`; B when it lies
#'   above `margin`; the two are equivalent when it lies inside
#'   `[-margin, margin]`.
#' @param bounds Known bounds `c(lower, upper)` of the losses of both
#'   algorithms. Paired differences then lie in `[lower - upper, upper - lower]`.
#'   A loss outside the bounds is an error, never clipped.
#' @param boundary One of [seqbench_boundaries()]. `"betting"` (default) is the
#'   hedged capital confidence sequence of Waudby-Smith & Ramdas (2024);
#'   `"empirical_bernstein"` and `"hoeffding"` are their conservative
#'   predictable-plug-in references; `"naive_fixed"` is a fixed-sample t
#'   interval recomputed at every step, **not valid** under optional stopping,
#'   provided only as a negative control for experiments.
#' @param paired Must be `TRUE`. Unpaired designs are not implemented; the
#'   argument exists so that the limitation is explicit.
#' @param cost_per_round Default cost of one `(instance, seed)` evaluation of
#'   both algorithms together, used when the data carry no `cost` column.
#' @param budget Total cost after which the procedure stops with an
#'   inconclusive outcome. `Inf` for no cost limit.
#' @param n_max Maximum number of instances. `Inf` for no limit.
#' @param betting_c,betting_theta,betting_grid Tuning of the betting boundary:
#'   truncation constant of the bets (default 1/2), hedging weight (default
#'   1/2) and grid size on `[0, 1]` (default 1001). Also used as the
#'   truncation constant of the empirical-Bernstein boundary.
#'
#' @return An object of class `seqbench_design` (a list).
#' @seealso [initialize_comparison()], [planning_horizon()]
#' @references Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of
#'   bounded random variables by betting. *Journal of the Royal Statistical
#'   Society Series B*, 86(1), 1-27. \doi{10.1093/jrsssb/qkad009}
#' @export
#' @examples
#' design <- comparison_design(alpha = 0.05, margin = 0.02, bounds = c(0, 1))
#' design
comparison_design <- function(alpha = 0.05, margin, bounds, boundary = "betting",
                              paired = TRUE, cost_per_round = 1, budget = Inf,
                              n_max = Inf, betting_c = 0.5, betting_theta = 0.5,
                              betting_grid = 1001L) {
  check_prob(alpha)
  check_positive(margin)
  check_bounds(bounds)
  boundary <- rlang::arg_match(boundary, seqbench_boundaries())
  if (!isTRUE(paired)) {
    cli::cli_abort(c(
      "Unpaired designs are not implemented.",
      "i" = "seqbench compares algorithms on paired losses (same instance, same seed). Set {.code paired = TRUE}."))
  }
  check_positive(cost_per_round)
  check_positive(budget, allow_inf = TRUE)
  check_positive(n_max, allow_inf = TRUE)
  if (is.finite(n_max) && n_max != round(n_max)) {
    cli::cli_abort("{.arg n_max} must be a whole number or Inf, not {.val {n_max}}.")
  }
  check_prob(betting_c)
  check_prob(betting_theta)
  range <- bounds[2] - bounds[1]
  if (margin >= range) {
    cli::cli_abort(c(
      "{.arg margin} must be smaller than the range of the paired difference.",
      "x" = "margin = {margin}, but differences lie in [{-range}, {range}]."))
  }
  structure(list(
    alpha = alpha, margin = margin, bounds = bounds,
    diff_bounds = c(bounds[1] - bounds[2], bounds[2] - bounds[1]),
    boundary = boundary, paired = TRUE, cost_per_round = cost_per_round,
    budget = budget, n_max = n_max,
    betting = list(c = betting_c, theta = betting_theta, grid = as.integer(betting_grid)),
    valid = boundary != "naive_fixed",
    created = Sys.time()
  ), class = "seqbench_design")
}

#' @export
print.seqbench_design <- function(x, ...) {
  assert_design(x)
  cli::cli_h1("seqbench comparison design")
  cli::cli_inform(c(
    "*" = "Boundary: {.field {x$boundary}}{if (!x$valid) ' (INVALID: negative control, not a confidence sequence)' else ''}",
    "*" = "alpha = {.val {x$alpha}}, margin = {.val {x$margin}} (loss units)",
    "*" = "Loss bounds: [{x$bounds[1]}, {x$bounds[2]}]; paired difference in [{x$diff_bounds[1]}, {x$diff_bounds[2]}]",
    "*" = "Budget: {if (is.finite(x$budget)) x$budget else 'unlimited'} cost units, n_max = {if (is.finite(x$n_max)) x$n_max else 'unlimited'} instances, cost per round = {x$cost_per_round}"
  ))
  invisible(x)
}

#' Planning horizon for a distance to the decision threshold
#'
#' Number of instances after which the confidence sequence is expected to have
#' declared the true hypothesis, when the true mean paired difference is at
#' distance `distance` from the nearest threshold `-margin` or `margin`. Two
#' modes:
#'
#' * `sd = NULL` (default): the **guaranteed** horizon
#'   `t*(Delta) = min{t : 2 w_t < |Delta|}` of the predictable-plug-in
#'   Hoeffding boundary, whose half-width `w_t` is deterministic. On the
#'   coverage event (probability at least `1 - alpha`) the Hoeffding comparison
#'   has stopped with the correct declaration by then (theory note, Prop. 2).
#'   It is distribution-free and, for small margins, very conservative: it can
#'   exceed what the default betting boundary needs by two or three orders of
#'   magnitude.
#' * `sd` given (standard deviation of the paired difference, loss units): an
#'   **approximation** that plugs `sd` into the empirical-Bernstein half-width
#'   in place of the running variance estimate. It is not a bound; it is the
#'   order of magnitude a variance-adaptive boundary needs, and the betting
#'   boundary is usually faster still. Use a pilot or a conservative guess for
#'   `sd`.
#'
#' Both charge the full half-width twice (worst-case position of the centre);
#' typical stopping times are about half of the returned value, as measured in
#' the package's pilot study. Both are computable before any data are collected
#' and serve to size `budget` and `n_max`.
#'
#' @param design A [comparison_design()].
#' @param distance Positive distance `|mu| - margin` (superiority) or
#'   `margin - |mu|` (equivalence), in loss units.
#' @param sd Optional standard deviation of the paired difference, loss units.
#' @param max_t Search limit.
#'
#' @return An integer number of instances, or `Inf` if not reached by `max_t`.
#' @export
#' @examples
#' design <- comparison_design(margin = 0.02, bounds = c(0, 1))
#' planning_horizon(design, distance = 0.05)             # guaranteed, Hoeffding
#' planning_horizon(design, distance = 0.05, sd = 0.1)   # variance-based approximation
planning_horizon <- function(design, distance, sd = NULL, max_t = 1e6) {
  assert_design(design)
  check_positive(distance)
  rng <- design$diff_bounds[2] - design$diff_bounds[1]
  d <- distance / rng
  if (!is.null(sd)) {
    check_positive(sd)
    sigma2 <- (sd / rng)^2
    if (sigma2 > 0.25) cli::cli_abort("{.arg sd} exceeds the maximum possible for the declared bounds.")
  }
  check_positive(max_t, allow_inf = FALSE)
  max_t <- floor(max_t)
  n <- min(1024L, max_t)
  repeat {
    w <- if (is.null(sd)) hoeffding_half_width(n, design$alpha)
         else eb_expected_half_width(n, design$alpha, sigma2, design$betting$c)
    hit <- which(2 * w < d)
    if (length(hit)) return(hit[1])
    if (n >= max_t) return(Inf)
    n <- min(n * 4L, as.integer(max_t))
  }
}

# Expected half-width of the PrPl-EB boundary when the running variance
# estimate is replaced by a known sigma2 (approximation, not a bound).
eb_expected_half_width <- function(n, alpha, sigma2, c) {
  t <- seq_len(n)
  lam <- pmin(sqrt(2 * log(2 / alpha) / (sigma2 * t * log(1 + t))), c)
  psi <- (-log1p(-lam) - lam) / 4
  (log(2 / alpha) + cumsum(4 * sigma2 * psi)) / cumsum(lam)
}
