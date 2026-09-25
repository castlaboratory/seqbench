# seqbench

<!-- badges: start -->
[![R-CMD-check](https://github.com/castlaboratory/seqbench/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/castlaboratory/seqbench/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**Anytime-valid sequential benchmarking of algorithms.**

`seqbench` answers one question: *when can I stop comparing algorithm A with
algorithm B without invalidating the inference?* It treats a benchmark as a
sequential experiment on **paired losses** (same instance, same seed) and
maintains a confidence sequence for the mean paired difference that stays
valid under optional stopping, repeated inspection and adaptive budgets.

Three declarations are possible, all relative to a pre-registered margin of
practical relevance `δ`:

| Declaration | Condition on the current confidence sequence `C_t` |
|---|---|
| A is relevantly better | upper limit of `C_t` below `−δ` |
| B is relevantly better | lower limit of `C_t` above `+δ` |
| Practically equivalent | `C_t` inside `[−δ, +δ]` |
| Inconclusive | otherwise (including budget exhausted) |

## Status

Pre-alpha. The API below is a design target, not yet implemented.

```r
design <- comparison_design(alpha = 0.05, margin = 0.02, paired = TRUE)
state  <- initialize_comparison(design)
state  <- update_comparison(state, paired_losses)
stopping_decision(state)
comparison_report(state)
```

## Installation

```r
# install.packages("pak")
pak::pak("castlaboratory/seqbench")
```

## Related software

`seqbench` is not a general confidence-sequence library. For that see
[confseq](https://github.com/gostevehoward/confseq) (Python),
[safestats](https://cran.r-project.org/package=safestats) and
[avlm](https://cran.r-project.org/package=avlm) (R). The CRAN package
[seqcomp](https://cran.r-project.org/package=seqcomp) compares probabilistic
*forecasts* sequentially; `seqbench` compares *algorithms* on paired losses
with an equivalence margin and explicit cost accounting.

## License

GPL (>= 3). © CAST Lab, Universidade Federal de Pernambuco.
