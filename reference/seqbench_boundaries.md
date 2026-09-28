# Names of the available boundaries

Names of the available boundaries

## Usage

``` r
seqbench_boundaries()
```

## Value

A character vector. `"betting"` is the default in
[`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md);
`"hoeffding"` and `"empirical_bernstein"` are conservative references;
`"bernstein_declared"` requires a declared upper bound on the standard
deviation of the paired difference and is valid only if that bound
holds; `"naive_fixed"` is an invalid negative control for experiments.

## Examples

``` r
seqbench_boundaries()
#> [1] "betting"             "empirical_bernstein" "hoeffding"          
#> [4] "bernstein_declared"  "naive_fixed"        
```
