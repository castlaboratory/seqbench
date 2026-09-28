# Design a sequential comparison of two algorithms

Fixes, before any data is seen, everything the procedure needs: the
error level, the margin of practical equivalence, the loss bounds, the
boundary and the budget. The design is immutable;
[`initialize_comparison()`](https://castlaboratory.github.io/seqbench/reference/initialize_comparison.md)
turns it into a state that
[`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md)
advances.

## Usage

``` r
comparison_design(
  alpha = 0.05,
  margin,
  bounds,
  boundary = "betting",
  paired = TRUE,
  cost_per_round = 1,
  budget = Inf,
  n_max = Inf,
  betting_c = 0.5,
  betting_theta = 0.5,
  betting_grid = 1001L,
  sd_max = NULL
)
```

## Arguments

- alpha:

  Miscoverage level of the confidence sequence, in (0, 1). The
  probability that the procedure ever issues a false declaration (of any
  of the three kinds) is at most `alpha`.

- margin:

  Margin of practical equivalence `delta > 0`, in loss units. Algorithm
  A is declared relevantly better when the whole confidence sequence for
  `mean(loss_a - loss_b)` lies below `-margin`; B when it lies above
  `margin`; the two are equivalent when it lies inside
  `[-margin, margin]`.

- bounds:

  Known bounds `c(lower, upper)` of the losses of both algorithms.
  Paired differences then lie in `[lower - upper, upper - lower]`. A
  loss outside the bounds is an error, never clipped.

- boundary:

  One of
  [`seqbench_boundaries()`](https://castlaboratory.github.io/seqbench/reference/seqbench_boundaries.md).
  `"betting"` (default) is the hedged capital confidence sequence of
  Waudby-Smith & Ramdas (2024); `"empirical_bernstein"` and
  `"hoeffding"` are their conservative predictable-plug-in references;
  `"bernstein_declared"` is a predictable-plug-in Bennett confidence
  sequence that uses a declared upper bound `sd_max` on the standard
  deviation of the paired difference and is valid only if that bound
  holds. Its gain over the betting boundary is modest (a logarithmic
  factor in the declared variance): for bounded observations the Bennett
  bound keeps a range term of order `log(2/alpha) / t` whatever the
  variance, so the boundary is provided for completeness and for the
  package's cost study rather than as a shortcut; `"naive_fixed"` is a
  fixed-sample t interval recomputed at every step, **not valid** under
  optional stopping, provided only as a negative control for
  experiments.

- paired:

  Must be `TRUE`. Unpaired designs are not implemented; the argument
  exists so that the limitation is explicit.

- cost_per_round:

  Default cost of one `(instance, seed)` evaluation of both algorithms
  together, used when the data carry no `cost` column.

- budget:

  Total cost after which the procedure stops with an inconclusive
  outcome. `Inf` for no cost limit.

- n_max:

  Maximum number of instances. `Inf` for no limit.

- betting_c, betting_theta, betting_grid:

  Tuning of the betting boundary: truncation constant of the bets
  (default 1/2, the value recommended by Waudby-Smith & Ramdas, who also
  suggest 3/4), hedging weight (default 1/2) and grid size on `[0, 1]`
  (default 1001). `betting_c` is also the truncation constant of the
  empirical-Bernstein boundary. Larger values bet more aggressively and
  stop earlier: in the package's simulation study (2,000 replications
  per cell) `betting_c = 0.9` needed about 40% fewer instances than 1/2
  at low noise and 20% fewer at high noise, with the largest observed
  coverage failure rising from 0.6% to 1.4% at `alpha = 0.05`; `0.99`
  gains little more. Coverage is guaranteed for any value in (0, 1).

- sd_max:

  Required for `boundary = "bernstein_declared"`: declared upper bound
  on the standard deviation of the per-instance (seed-averaged) paired
  difference, in loss units. This is an additional assumption (A5)
  recorded in the result contract.

## Value

An object of class `seqbench_design` (a list).

## References

Waudby-Smith, I. and Ramdas, A. (2024). Estimating means of bounded
random variables by betting. *Journal of the Royal Statistical Society
Series B*, 86(1), 1-27.
[doi:10.1093/jrsssb/qkad009](https://doi.org/10.1093/jrsssb/qkad009)

## See also

[`initialize_comparison()`](https://castlaboratory.github.io/seqbench/reference/initialize_comparison.md),
[`planning_horizon()`](https://castlaboratory.github.io/seqbench/reference/planning_horizon.md)

## Examples

``` r
design <- comparison_design(alpha = 0.05, margin = 0.02, bounds = c(0, 1))
design
#> 
#> ── seqbench comparison design ──────────────────────────────────────────────────
#> • Boundary: betting
#> • alpha = 0.05, margin = 0.02 (loss units)
#> • Loss bounds: [0, 1]; paired difference in [-1, 1]
#> • Budget: unlimited cost units, n_max = unlimited instances, cost per round = 1
```
