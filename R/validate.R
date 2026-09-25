# Validation of paired-loss data ---------------------------------------------
#
# The package never converts invalid input into a valid-looking result:
# missing values, losses outside the declared bounds, negative costs and
# instances that reappear (adaptive reuse is not i.i.d. sampling) are errors.

required_cols <- c("instance", "loss_a", "loss_b")

validate_paired_losses <- function(losses, design, seen_instances,
                                   call = rlang::caller_env()) {
  if (!is.data.frame(losses)) {
    cli::cli_abort("{.arg losses} must be a data frame with columns {.field {required_cols}}.",
                   call = call)
  }
  missing <- setdiff(required_cols, names(losses))
  if (length(missing)) {
    cli::cli_abort("{.arg losses} is missing column{?s} {.field {missing}}.", call = call)
  }
  if (nrow(losses) == 0L) {
    cli::cli_abort("{.arg losses} has no rows.", call = call)
  }
  losses <- tibble::as_tibble(losses)
  if (!"seed" %in% names(losses)) losses$seed <- 1L
  if (!"cost" %in% names(losses)) losses$cost <- design$cost_per_round

  for (col in c("loss_a", "loss_b", "cost")) {
    if (!is.numeric(losses[[col]])) {
      cli::cli_abort("Column {.field {col}} must be numeric.", call = call)
    }
    if (anyNA(losses[[col]])) {
      bad <- which(is.na(losses[[col]]))
      cli::cli_abort(c(
        "Column {.field {col}} has {cli::qty(length(bad))}{length(bad)} missing value{?s} (row{?s} {head(bad, 5)}).",
        "i" = "seqbench does not impute or drop rows; decide on the loss of a failed run before updating."),
        call = call)
    }
  }
  if (anyNA(losses$instance) || anyNA(losses$seed)) {
    cli::cli_abort("Columns {.field instance} and {.field seed} must not contain missing values.",
                   call = call)
  }
  b <- design$bounds
  for (col in c("loss_a", "loss_b")) {
    out <- which(losses[[col]] < b[1] | losses[[col]] > b[2])
    if (length(out)) {
      cli::cli_abort(c(
        "{cli::qty(length(out))}{length(out)} value{?s} of {.field {col}} {cli::qty(length(out))}fall{?s/} outside the declared bounds [{b[1]}, {b[2]}] (row{?s} {head(out, 5)}).",
        "i" = "Losses are never clipped. Fix the data or declare wider bounds in {.fn comparison_design}."),
        call = call)
    }
  }
  if (any(losses$cost < 0)) {
    cli::cli_abort("Column {.field cost} must be non-negative.", call = call)
  }
  key <- paste(losses$instance, losses$seed, sep = "\r")
  if (anyDuplicated(key)) {
    dup <- unique(losses$instance[duplicated(key)])
    cli::cli_abort(c(
      "Duplicated (instance, seed) pairs for {cli::qty(length(dup))}instance{?s} {.val {head(dup, 5)}}.",
      "i" = "Each (instance, seed) may be evaluated once."), call = call)
  }
  reused <- unique(losses$instance[as.character(losses$instance) %in% seen_instances])
  if (length(reused)) {
    cli::cli_abort(c(
      "{cli::qty(length(reused))}Instance{?s} {.val {head(reused, 5)}} {cli::qty(length(reused))}{?was/were} already observed in a previous update.",
      "x" = "Adding seeds to an instance after seeing its losses is adaptive reuse, not i.i.d. sampling; the sequential unit is the instance (assumption A4).",
      "i" = "Give the new runs a fresh instance or plan all seeds of an instance before observing any of them."),
      call = call)
  }
  losses
}
