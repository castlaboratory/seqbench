# Planning horizon for a distance to the decision threshold

Number of instances after which the confidence sequence is expected to
have declared the true hypothesis, when the true mean paired difference
is at distance `distance` from the nearest threshold `-margin` or
`margin`. Two modes:

## Usage

``` r
planning_horizon(design, distance, sd = NULL, max_t = 1e+06)
```

## Arguments

- design:

  A
  [`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md).

- distance:

  Positive distance `|mu| - margin` (superiority) or `margin - |mu|`
  (equivalence), in loss units.

- sd:

  Optional standard deviation of the paired difference, loss units.

- max_t:

  Search limit.

## Value

An integer number of instances, or `Inf` if not reached by `max_t`.

## Details

- `sd = NULL` (default): the **guaranteed** horizon
  `t*(Delta) = min{t : 2 w_t < |Delta|}` of the predictable-plug-in
  Hoeffding boundary, whose half-width `w_t` is deterministic. On the
  coverage event (probability at least `1 - alpha`) the Hoeffding
  comparison has stopped with the correct declaration by then (theory
  note, Prop. 2). It is distribution-free and, for small margins, very
  conservative: it can exceed what the default betting boundary needs by
  two or three orders of magnitude.

- `sd` given (standard deviation of the paired difference, loss units):
  an **approximation** that plugs `sd` into the empirical-Bernstein
  half-width in place of the running variance estimate. It is not a
  bound; it is the order of magnitude a variance-adaptive boundary
  needs, and the betting boundary is usually faster still. Use a pilot
  or a conservative guess for `sd`.

Both charge the full half-width twice (worst-case position of the
centre); typical stopping times are about half of the returned value, as
measured in the package's pilot study. Both are computable before any
data are collected and serve to size `budget` and `n_max`.

## Examples

``` r
design <- comparison_design(margin = 0.02, bounds = c(0, 1))
planning_horizon(design, distance = 0.05)             # guaranteed, Hoeffding
#> [1] 52296
planning_horizon(design, distance = 0.05, sd = 0.1)   # variance-based approximation
#> [1] 640
```
