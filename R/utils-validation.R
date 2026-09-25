# Internal argument checks. All errors are raised with cli_abort so that the
# message names the argument and the offending value.

check_prob <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    cli::cli_abort("{.arg {arg}} must be a single number in (0, 1), not {.val {x}}.",
                   call = call)
  }
  invisible(x)
}

check_positive <- function(x, arg = rlang::caller_arg(x), allow_inf = FALSE,
                           call = rlang::caller_env()) {
  ok <- is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0 &&
    (allow_inf || is.finite(x))
  if (!ok) {
    what <- if (allow_inf) "a single positive number or Inf" else "a single positive finite number"
    cli::cli_abort("{.arg {arg}} must be {what}, not {.val {x}}.", call = call)
  }
  invisible(x)
}

check_bounds <- function(bounds, call = rlang::caller_env()) {
  if (!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds) ||
      any(!is.finite(bounds)) || bounds[1] >= bounds[2]) {
    cli::cli_abort(c(
      "{.arg bounds} must be a finite numeric vector {.code c(lower, upper)} with {.code lower < upper}.",
      "x" = "Got {.val {bounds}}."), call = call)
  }
  invisible(bounds)
}

assert_design <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!inherits(x, "seqbench_design")) {
    cli::cli_abort("{.arg {arg}} must be a {.cls seqbench_design} created by {.fn comparison_design}.",
                   call = call)
  }
  invisible(x)
}

assert_comparison <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!inherits(x, "seqbench_comparison")) {
    cli::cli_abort("{.arg {arg}} must be a {.cls seqbench_comparison} created by {.fn initialize_comparison}.",
                   call = call)
  }
  invisible(x)
}
