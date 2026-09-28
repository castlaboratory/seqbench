# Glance at a comparison

Glance at a comparison

## Usage

``` r
# S3 method for class 'seqbench_comparison'
glance(x, ...)
```

## Arguments

- x:

  A `seqbench_comparison`.

- ...:

  Unused.

## Value

A one-row tibble: `boundary`, `valid`, `alpha`, `margin`, `n_instances`,
`n_evaluations`, `cost`, `estimate`, `lower`, `upper`, `decision`,
`stopping_reason`.

## Examples

``` r
state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
glance(state)
#> # A tibble: 1 × 12
#>   boundary valid alpha margin n_instances n_evaluations  cost estimate lower
#>   <chr>    <lgl> <dbl>  <dbl>       <int>         <int> <dbl>    <dbl> <dbl>
#> 1 betting  TRUE   0.05   0.02           0             0     0       NA    -1
#> # ℹ 3 more variables: upper <dbl>, decision <chr>, stopping_reason <chr>
```
