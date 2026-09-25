# Comparison state ------------------------------------------------------------

#' Start a sequential comparison
#'
#' Creates the state object that [update_comparison()] advances one instance
#' at a time. Nothing is observed yet; the confidence sequence is `[a, b]`, the
#' full range of the paired difference.
#'
#' @param design A [comparison_design()].
#'
#' @return An object of class `seqbench_comparison`. Its main components are
#'   `design`, `data` (all rows seen so far), `trajectory` (one row per
#'   instance: estimate, confidence sequence, running intersection, cost and
#'   decision at that time), `decision`, `stopping_reason`, `diagnostics` and
#'   `meta` (versions, timestamps). Use [tidy()] for the trajectory, [glance()]
#'   for a one-row summary and [comparison_report()] for the full contract.
#' @export
#' @examples
#' design <- comparison_design(margin = 0.02, bounds = c(0, 1))
#' state <- initialize_comparison(design)
#' state
initialize_comparison <- function(design) {
  assert_design(design)
  bnd <- boundary_init(design$boundary, alpha = design$alpha,
                       c = design$betting$c, theta = design$betting$theta,
                       grid = design$betting$grid)
  structure(list(
    design = design,
    boundary_state = bnd,
    data = empty_data(),
    trajectory = empty_trajectory(),
    running = c(lower = 0, upper = 1),       # running intersection on [0, 1]
    cost = 0,
    decision = "continue",
    stopping_reason = NA_character_,
    stopped = FALSE,
    diagnostics = list(n_instances = 0L, n_evaluations = 0L,
                       loss_a_at_bounds = 0L, loss_b_at_bounds = 0L,
                       cs_empty = FALSE),
    meta = list(seqbench_version = as.character(utils::packageVersion("seqbench")),
                r_version = R.version.string,
                created = Sys.time(), updated = Sys.time())
  ), class = "seqbench_comparison")
}

empty_data <- function() {
  tibble::tibble(t = integer(), instance = character(), seed = character(),
                 loss_a = numeric(), loss_b = numeric(), cost = numeric())
}

empty_trajectory <- function() {
  tibble::tibble(t = integer(), instance = character(), n_seeds = integer(),
                 x = numeric(), estimate = numeric(), lower = numeric(),
                 upper = numeric(), lower_running = numeric(),
                 upper_running = numeric(), cost = numeric(), cost_cum = numeric(),
                 decision = character())
}

#' Update a comparison with new paired losses
#'
#' Absorbs the losses of one or more new instances, in the order given, and
#' advances the confidence sequence and the decision after each instance.
#' The sequential unit is the instance: all rows of an instance (its seeds)
#' are averaged into one observation. Updating stops at the first instance
#' that triggers a terminal decision or exhausts the budget; rows after it are
#' not consumed and a warning says so.
#'
#' @param state A `seqbench_comparison` from [initialize_comparison()].
#' @param losses A data frame with columns `instance`, `loss_a`, `loss_b`,
#'   and optionally `seed` (default 1) and `cost` (default
#'   `design$cost_per_round` per row). Losses must lie within the design
#'   bounds; instances must not have been seen before.
#'
#' @return The updated `seqbench_comparison`.
#' @export
#' @examples
#' design <- comparison_design(margin = 0.05, bounds = c(0, 1), boundary = "hoeffding")
#' state <- initialize_comparison(design)
#' set.seed(1)
#' losses <- data.frame(instance = 1:50, loss_a = runif(50, 0, 0.5),
#'                      loss_b = runif(50, 0.3, 0.8))
#' state <- update_comparison(state, losses)
#' stopping_decision(state)
update_comparison <- function(state, losses) {
  assert_comparison(state)
  if (isTRUE(state$stopped)) {
    cli::cli_abort(c(
      "This comparison already stopped ({.val {state$decision}}, reason {.val {state$stopping_reason}}).",
      "i" = "Continuing after a terminal decision would change the pre-registered protocol. Start a new comparison if you need more data."))
  }
  design <- state$design
  losses <- validate_paired_losses(losses, design, unique(state$data$instance))
  losses$instance <- as.character(losses$instance)
  losses$seed <- as.character(losses$seed)
  a <- design$diff_bounds[1]
  rng <- design$diff_bounds[2] - a
  margin01 <- design$margin / rng

  order_inst <- unique(losses$instance)
  groups <- split(seq_len(nrow(losses)), factor(losses$instance, levels = order_inst))
  n_new <- length(order_inst)
  # Preallocated trajectory accumulators (one slot per new instance)
  acc <- list(t = integer(n_new), n_seeds = integer(n_new), x = numeric(n_new),
              estimate = numeric(n_new), lower = numeric(n_new), upper = numeric(n_new),
              lower_running = numeric(n_new), upper_running = numeric(n_new),
              cost = numeric(n_new), cost_cum = numeric(n_new), decision = character(n_new))
  row_t <- integer(nrow(losses))
  k <- 0L
  for (inst in order_inst) {
    idx <- groups[[inst]]
    la <- losses$loss_a[idx]; lb <- losses$loss_b[idx]; cs <- losses$cost[idx]
    t <- state$diagnostics$n_instances + 1L
    d_bar <- mean(la - lb)
    x <- (d_bar - a) / rng
    x <- min(max(x, 0), 1)                  # guard against floating error only
    cost <- sum(cs)

    state$boundary_state <- boundary_update(state$boundary_state, x)
    ci <- boundary_interval(state$boundary_state)
    if (design$valid) {
      state$running <- c(lower = max(state$running[["lower"]], ci[["lower"]]),
                         upper = min(state$running[["upper"]], ci[["upper"]]))
    } else {
      state$running <- c(lower = ci[["lower"]], upper = ci[["upper"]])  # no intersection: invalid anyway
    }
    if (anyNA(ci[c("lower", "upper")])) {
      state$running <- c(lower = NA_real_, upper = NA_real_)
    }
    state$cost <- state$cost + cost
    state$diagnostics$n_instances <- t
    state$diagnostics$n_evaluations <- state$diagnostics$n_evaluations + length(idx)
    state$diagnostics$loss_a_at_bounds <- state$diagnostics$loss_a_at_bounds +
      sum(la == design$bounds[1] | la == design$bounds[2])
    state$diagnostics$loss_b_at_bounds <- state$diagnostics$loss_b_at_bounds +
      sum(lb == design$bounds[1] | lb == design$bounds[2])

    dec <- decide(state$running, margin01, state$cost, t, design)
    state$decision <- dec$decision
    state$stopping_reason <- dec$reason
    state$stopped <- dec$decision != "continue"
    if (identical(dec$reason, "cs_empty")) state$diagnostics$cs_empty <- TRUE

    k <- k + 1L
    row_t[idx] <- t
    acc$t[k] <- t; acc$n_seeds[k] <- length(idx); acc$x[k] <- x
    acc$estimate[k] <- a + rng * ci[["estimate"]]
    acc$lower[k] <- a + rng * ci[["lower"]]; acc$upper[k] <- a + rng * ci[["upper"]]
    acc$lower_running[k] <- a + rng * state$running[["lower"]]
    acc$upper_running[k] <- a + rng * state$running[["upper"]]
    acc$cost[k] <- cost; acc$cost_cum[k] <- state$cost; acc$decision[k] <- dec$decision

    if (state$stopped) {
      left <- n_new - k
      if (left > 0L) {
        cli::cli_warn(c(
          "Stopped at instance {.val {inst}} (t = {t}) with decision {.val {dec$decision}}; {cli::qty(left)}{left} later instance{?s} {?was/were} not consumed.",
          "i" = "Under the pre-registered protocol the experiment ends here."))
      }
      break
    }
  }
  keep <- row_t > 0L
  new_data <- tibble::tibble(
    t = row_t[keep], instance = losses$instance[keep], seed = losses$seed[keep],
    loss_a = losses$loss_a[keep], loss_b = losses$loss_b[keep], cost = losses$cost[keep])
  new_traj <- tibble::as_tibble(lapply(acc, function(v) v[seq_len(k)]))
  new_traj <- tibble::add_column(new_traj, instance = order_inst[seq_len(k)], .after = "t")
  state$data <- rbind(state$data, new_data)
  state$trajectory <- rbind(state$trajectory, new_traj)
  state$meta$updated <- Sys.time()
  state
}

# Decision rules on the [0, 1] scale (paper/theory.md, Algorithm 1).
decide <- function(running, margin01, cost, t, design) {
  L <- running[["lower"]]; U <- running[["upper"]]
  lo_thr <- 0.5 - margin01           # -delta on the [0, 1] scale (mu = 0 maps to 1/2)
  hi_thr <- 0.5 + margin01
  if (anyNA(c(L, U)) || U < L) return(list(decision = "inconclusive", reason = "cs_empty"))
  if (U < lo_thr) return(list(decision = "A", reason = "declaration"))
  if (L > hi_thr) return(list(decision = "B", reason = "declaration"))
  if (L >= lo_thr && U <= hi_thr) return(list(decision = "equivalent", reason = "declaration"))
  if (cost >= design$budget) return(list(decision = "inconclusive", reason = "budget"))
  if (t >= design$n_max) return(list(decision = "inconclusive", reason = "n_max"))
  list(decision = "continue", reason = NA_character_)
}

#' Current decision of a comparison
#'
#' @param state A `seqbench_comparison`.
#' @return A character scalar: `"A"` (A relevantly better), `"B"`,
#'   `"equivalent"`, `"continue"` or `"inconclusive"`, with attribute
#'   `reason` (`"declaration"`, `"budget"`, `"n_max"`, `"cs_empty"` or `NA`).
#' @export
#' @examples
#' state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
#' stopping_decision(state)
stopping_decision <- function(state) {
  assert_comparison(state)
  structure(state$decision, reason = state$stopping_reason)
}
