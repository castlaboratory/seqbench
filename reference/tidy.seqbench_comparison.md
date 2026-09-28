# Tidy the trajectory of a comparison

One row per instance with the estimate, the confidence sequence at that
time, its running intersection, cumulative cost and the decision the
procedure would have taken there.

## Usage

``` r
# S3 method for class 'seqbench_comparison'
tidy(x, ...)
```

## Arguments

- x:

  A `seqbench_comparison`.

- ...:

  Unused.

## Value

A tibble with columns `t`, `instance`, `n_seeds`, `x`, `estimate`,
`lower`, `upper`, `lower_running`, `upper_running`, `cost`, `cost_cum`,
`decision`.

## Examples

``` r
design <- comparison_design(margin = 0.05, bounds = c(0, 1), boundary = "hoeffding")
state <- update_comparison(initialize_comparison(design),
  data.frame(instance = 1:5, loss_a = c(.1, .2, .1, .3, .2), loss_b = c(.5, .6, .4, .7, .5)))
tidy(state)
#> # A tibble: 5 × 12
#>       t instance n_seeds     x estimate lower upper lower_running upper_running
#>   <int> <chr>      <int> <dbl>    <dbl> <dbl> <dbl>         <dbl>         <dbl>
#> 1     1 1              1  0.3    -0.4      -1     1            -1             1
#> 2     2 2              1  0.3    -0.4      -1     1            -1             1
#> 3     3 3              1  0.35   -0.367    -1     1            -1             1
#> 4     4 4              1  0.3    -0.375    -1     1            -1             1
#> 5     5 5              1  0.35   -0.36     -1     1            -1             1
#> # ℹ 3 more variables: cost <dbl>, cost_cum <dbl>, decision <chr>
```
