# Current decision of a comparison

Current decision of a comparison

## Usage

``` r
stopping_decision(state)
```

## Arguments

- state:

  A `seqbench_comparison`.

## Value

A character scalar: `"A"` (A relevantly better), `"B"`, `"equivalent"`,
`"continue"` or `"inconclusive"`, with attribute `reason`
(`"declaration"`, `"budget"`, `"n_max"`, `"cs_empty"` or `NA`).

## Examples

``` r
state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
stopping_decision(state)
#> [1] "continue"
#> attr(,"reason")
#> [1] NA
```
