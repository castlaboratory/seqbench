#' Plot the inferential trajectory of a comparison
#'
#' Running-intersection confidence sequence for the mean paired difference
#' against the number of instances, with the equivalence band
#' `[-margin, margin]` and the zero line.
#'
#' @param object A `seqbench_comparison` with at least one instance.
#' @param ... Unused.
#' @return A ggplot object.
#' @exportS3Method ggplot2::autoplot seqbench_comparison
#' @examples
#' design <- comparison_design(margin = 0.05, bounds = c(0, 1), boundary = "hoeffding")
#' set.seed(2)
#' losses <- data.frame(instance = 1:40, loss_a = runif(40, 0, .5), loss_b = runif(40, .2, .7))
#' state <- update_comparison(initialize_comparison(design), losses)
#' ggplot2::autoplot(state)
autoplot.seqbench_comparison <- function(object, ...) {
  assert_comparison(object)
  tr <- object$trajectory
  if (nrow(tr) == 0L) cli::cli_abort("Nothing to plot: no instances observed yet.")
  m <- object$design$margin
  ggplot2::ggplot(tr, ggplot2::aes(x = .data$t)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -m, ymax = m,
                      alpha = 0.12, fill = "grey40") +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_running, ymax = .data$upper_running),
                         alpha = 0.25, fill = "#2b6cb0") +
    ggplot2::geom_line(ggplot2::aes(y = .data$estimate), colour = "#2b6cb0") +
    ggplot2::labs(x = "Instances (t)", y = "mean(loss_a - loss_b)",
                  title = sprintf("Confidence sequence (%s, alpha = %s), decision: %s",
                                  object$design$boundary, object$design$alpha,
                                  decision_label(object$decision)),
                  subtitle = if (!object$design$valid) "INVALID negative control: not a confidence sequence" else NULL) +
    ggplot2::theme_minimal()
}
