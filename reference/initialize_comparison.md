# Start a sequential comparison

Creates the state object that
[`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md)
advances one instance at a time. Nothing is observed yet; the confidence
sequence is `[a, b]`, the full range of the paired difference.

## Usage

``` r
initialize_comparison(design)
```

## Arguments

- design:

  A
  [`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md).

## Value

An object of class `seqbench_comparison`. Its main components are
`design`, `data` (all rows seen so far), `trajectory` (one row per
instance: estimate, confidence sequence, running intersection, cost and
decision at that time), `decision`, `stopping_reason`, `diagnostics` and
`meta` (versions, timestamps). Use
[`tidy()`](https://generics.r-lib.org/reference/tidy.html) for the
trajectory,
[`glance()`](https://generics.r-lib.org/reference/glance.html) for a
one-row summary and
[`comparison_report()`](https://castlaboratory.github.io/seqbench/reference/comparison_report.md)
for the full contract.

## Examples

``` r
design <- comparison_design(margin = 0.02, bounds = c(0, 1))
state <- initialize_comparison(design)
state
#> 
#> ── seqbench comparison ─────────────────────────────────────────────────────────
#> • Boundary betting, alpha = 0.05, margin = 0.02, bounds [0, 1]
#> • Instances: 0, evaluations: 0, cost: 0
#> • Decision: continue sampling
```
