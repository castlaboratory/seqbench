# broom methods ---------------------------------------------------------------

#' Tidy the trajectory of a comparison
#'
#' One row per instance with the estimate, the confidence sequence at that
#' time, its running intersection, cumulative cost and the decision the
#' procedure would have taken there.
#'
#' @param x A `seqbench_comparison`.
#' @param ... Unused.
#' @return A tibble with columns `t`, `instance`, `n_seeds`, `x`,
#'   `estimate`, `lower`, `upper`, `lower_running`, `upper_running`, `cost`,
#'   `cost_cum`, `decision`.
#' @exportS3Method generics::tidy seqbench_comparison
#' @examples
#' design <- comparison_design(margin = 0.05, bounds = c(0, 1), boundary = "hoeffding")
#' state <- update_comparison(initialize_comparison(design),
#'   data.frame(instance = 1:5, loss_a = c(.1, .2, .1, .3, .2), loss_b = c(.5, .6, .4, .7, .5)))
#' tidy(state)
tidy.seqbench_comparison <- function(x, ...) {
  assert_comparison(x)
  x$trajectory
}

#' Glance at a comparison
#'
#' @param x A `seqbench_comparison`.
#' @param ... Unused.
#' @return A one-row tibble: `boundary`, `valid`, `alpha`, `margin`,
#'   `n_instances`, `n_evaluations`, `cost`, `estimate`, `lower`, `upper`,
#'   `decision`, `stopping_reason`.
#' @exportS3Method generics::glance seqbench_comparison
#' @examples
#' state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
#' glance(state)
glance.seqbench_comparison <- function(x, ...) {
  assert_comparison(x)
  r <- comparison_report(x)
  tibble::tibble(
    boundary = r$boundary, valid = r$valid, alpha = r$alpha, margin = r$margin,
    n_instances = r$resources$instances, n_evaluations = r$resources$evaluations,
    cost = r$resources$cost, estimate = r$estimate,
    lower = r$confidence_sequence[["lower"]], upper = r$confidence_sequence[["upper"]],
    decision = r$decision, stopping_reason = r$stopping_reason)
}
