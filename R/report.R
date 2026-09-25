# print / summary / report ----------------------------------------------------

decision_label <- function(d) {
  switch(d,
    A = "A is relevantly better",
    B = "B is relevantly better",
    equivalent = "A and B are practically equivalent",
    continue = "continue sampling",
    inconclusive = "inconclusive",
    d)
}

#' @export
print.seqbench_comparison <- function(x, ...) {
  assert_comparison(x)
  d <- x$design
  n <- x$diagnostics$n_instances
  cli::cli_h1("seqbench comparison")
  if (!d$valid) {
    cli::cli_alert_danger("Boundary {.field naive_fixed} is NOT a confidence sequence; decisions below carry no error guarantee.")
  }
  cli::cli_inform(c(
    "*" = "Boundary {.field {d$boundary}}, alpha = {d$alpha}, margin = {d$margin}, bounds [{d$bounds[1]}, {d$bounds[2]}]",
    "*" = "Instances: {n}, evaluations: {x$diagnostics$n_evaluations}, cost: {format(x$cost)}{if (is.finite(d$budget)) paste0(' / ', d$budget) else ''}"
  ))
  if (n > 0L) {
    last <- x$trajectory[n, ]
    cli::cli_inform(c(
      "*" = "Estimate of mean(loss_a - loss_b): {signif(last$estimate, 4)}",
      "*" = "Confidence sequence (running intersection): [{signif(last$lower_running, 4)}, {signif(last$upper_running, 4)}]"
    ))
  }
  cli::cli_inform(c("*" = "Decision: {.strong {decision_label(x$decision)}}{if (!is.na(x$stopping_reason)) paste0(' (', x$stopping_reason, ')') else ''}"))
  invisible(x)
}

#' Full result contract of a comparison
#'
#' Returns every element the protocol promises: estimand, estimate, current
#' confidence sequence, required assumptions, diagnostics, decision and
#' stopping reason, resources consumed, seeds and versions. `summary()` is an
#' alias.
#'
#' @param state A `seqbench_comparison`.
#' @param object A `seqbench_comparison`.
#' @param ... Unused.
#' @return A list of class `seqbench_report`.
#' @export
#' @examples
#' state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
#' comparison_report(state)
comparison_report <- function(state) {
  assert_comparison(state)
  d <- state$design
  n <- state$diagnostics$n_instances
  last <- if (n > 0L) state$trajectory[n, ] else NULL
  structure(list(
    estimand = "mu = E_P[ E_s[ loss_a(i, s) - loss_b(i, s) ] ], mean paired difference over the instance population P",
    estimate = if (n) last$estimate else NA_real_,
    confidence_sequence = c(lower = if (n) last$lower_running else d$diff_bounds[1],
                            upper = if (n) last$upper_running else d$diff_bounds[2]),
    alpha = d$alpha,
    margin = d$margin,
    boundary = d$boundary,
    valid = d$valid,
    assumptions = c(
      A1 = "instances are i.i.d. draws from the population P",
      A2 = "the number of seeds per instance is fixed before its losses are observed",
      A3 = sprintf("both losses lie in [%s, %s]", d$bounds[1], d$bounds[2]),
      A4 = "no instance is reused after its losses are observed"),
    decision = state$decision,
    decision_label = decision_label(state$decision),
    stopping_reason = state$stopping_reason,
    diagnostics = state$diagnostics,
    resources = list(instances = n, evaluations = state$diagnostics$n_evaluations,
                     cost = state$cost, budget = d$budget, n_max = d$n_max),
    instances = state$trajectory$instance,
    seeds = state$data[, c("instance", "seed")],
    versions = state$meta
  ), class = "seqbench_report")
}

#' @rdname comparison_report
#' @export
summary.seqbench_comparison <- function(object, ...) comparison_report(object)

#' @export
print.seqbench_report <- function(x, ...) {
  cli::cli_h1("seqbench report")
  cli::cli_inform(c(
    "Estimand" = x$estimand,
    "Decision" = "{.strong {x$decision_label}}{if (!is.na(x$stopping_reason)) paste0(' (', x$stopping_reason, ')') else ''}",
    "Estimate" = "{signif(x$estimate, 4)}; confidence sequence [{signif(x$confidence_sequence[['lower']], 4)}, {signif(x$confidence_sequence[['upper']], 4)}] at alpha = {x$alpha}, margin = {x$margin}",
    "Boundary" = "{x$boundary}{if (!x$valid) ' (INVALID negative control)' else ''}",
    "Resources" = "{x$resources$instances} instances, {x$resources$evaluations} evaluations, cost {format(x$resources$cost)}",
    "Versions" = "seqbench {x$versions$seqbench_version}, {x$versions$r_version}"
  ))
  cli::cli_text("Assumptions:")
  cli::cli_ul(paste0(names(x$assumptions), ": ", x$assumptions))
  invisible(x)
}
