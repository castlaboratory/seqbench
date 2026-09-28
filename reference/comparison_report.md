# Full result contract of a comparison

Returns every element the protocol promises: estimand, estimate, current
confidence sequence, required assumptions, diagnostics, decision and
stopping reason, resources consumed, seeds and versions.
[`summary()`](https://rdrr.io/r/base/summary.html) is an alias.

## Usage

``` r
comparison_report(state)

# S3 method for class 'seqbench_comparison'
summary(object, ...)
```

## Arguments

- state:

  A `seqbench_comparison`.

- object:

  A `seqbench_comparison`.

- ...:

  Unused.

## Value

A list of class `seqbench_report`.

## Examples

``` r
state <- initialize_comparison(comparison_design(margin = 0.02, bounds = c(0, 1)))
comparison_report(state)
#> 
#> ── seqbench report ─────────────────────────────────────────────────────────────
#> mu = E_P[ E_s[ loss_a(i, s) - loss_b(i, s) ] ], mean paired difference over the
#> instance population P
#> continue sampling
#> NA; confidence sequence [-1, 1] at alpha = 0.05, margin = 0.02
#> betting
#> 0 instances, 0 evaluations, cost 0
#> seqbench 0.1.0, R version 4.6.1 (2026-06-24)
#> Assumptions:
#> • A1: instances are i.i.d. draws from the population P
#> • A2: the number of seeds per instance is fixed before its losses are observed
#> • A3: both losses lie in [0, 1]
#> • A4: no instance is reused after its losses are observed
```
