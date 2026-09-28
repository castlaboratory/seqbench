# Update a comparison with new paired losses

Absorbs the losses of one or more new instances, in the order given, and
advances the confidence sequence and the decision after each instance.
The sequential unit is the instance: all rows of an instance (its seeds)
are averaged into one observation. Updating stops at the first instance
that triggers a terminal decision or exhausts the budget; rows after it
are not consumed and a warning says so.

## Usage

``` r
update_comparison(state, losses)
```

## Arguments

- state:

  A `seqbench_comparison` from
  [`initialize_comparison()`](https://castlaboratory.github.io/seqbench/reference/initialize_comparison.md).

- losses:

  A data frame with columns `instance`, `loss_a`, `loss_b`, and
  optionally `seed` (default 1) and `cost` (default
  `design$cost_per_round` per row). Losses must lie within the design
  bounds; instances must not have been seen before.

## Value

The updated `seqbench_comparison`.

## Examples

``` r
design <- comparison_design(margin = 0.05, bounds = c(0, 1), boundary = "hoeffding")
state <- initialize_comparison(design)
set.seed(1)
losses <- data.frame(instance = 1:50, loss_a = runif(50, 0, 0.5),
                     loss_b = runif(50, 0.3, 0.8))
state <- update_comparison(state, losses)
stopping_decision(state)
#> [1] "continue"
#> attr(,"reason")
#> [1] NA
```
