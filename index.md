# seqbench

**Anytime-valid sequential benchmarking of algorithms.**

`seqbench` answers one question: *when can I stop comparing algorithm A
with algorithm B without invalidating the inference?* It treats a
benchmark as a sequential experiment on **paired losses** (same
instance, same seed) and maintains a confidence sequence for the mean
paired difference that stays valid under optional stopping, repeated
inspection and adaptive budgets.

Three declarations are possible, all relative to a pre-registered margin
of practical relevance `δ`:

| Declaration | Condition on the current confidence sequence `C_t` |
|----|----|
| A is relevantly better | upper limit of `C_t` below `−δ` |
| B is relevantly better | lower limit of `C_t` above `+δ` |
| Practically equivalent | `C_t` inside `[−δ, +δ]` |
| Continue | none of the above and budget remains |
| Inconclusive | none of the above and the budget is exhausted (or `C_t` is empty) |

## Status

Version 0.1.0, first release; the API is experimental and may change in
minor releases. Implemented: five boundaries (betting, empirical
Bernstein, Hoeffding, Bernstein with a declared variance bound, and an
invalid fixed-sample negative control), instance-level updating with
seeds as clusters, the three-way decision rule, budget accounting, and
`print`, `summary`, `tidy`, `glance` and `autoplot` methods.

``` r

library(seqbench)

design <- comparison_design(alpha = 0.05, margin = 0.02, bounds = c(0, 1),
                            boundary = "betting", budget = 500)
state  <- initialize_comparison(design)

# paired losses: one row per (instance, seed); losses within `bounds`
losses <- data.frame(instance = 1:40, seed = 1,
                     loss_a = runif(40, 0.1, 0.4), loss_b = runif(40, 0.3, 0.6))
state  <- update_comparison(state, losses)

stopping_decision(state)     # "A", "B", "equivalent", "continue" or "inconclusive"
comparison_report(state)     # estimand, estimate, C_t, assumptions, cost, seeds, versions
tidy(state)                  # trajectory, one row per instance
ggplot2::autoplot(state)     # C_t against t with the [-margin, margin] band
```

Instances are the sequential unit; several seeds of the same instance
are averaged into one observation. Losses outside `bounds`, missing
values and instances that reappear after their losses were seen are
errors, never silently converted.

## Documentation

Full reference, the getting-started vignette and the changelog are on
the package site: <https://castlaboratory.github.io/seqbench/>.

## Companion article and experiments

The article *seqbench: Anytime-Valid Sequential Benchmarking of
Algorithms in R* (submitted to The R Journal) lives in `paper/rjournal/`
of this repository, and the simulation study and applications it reports
in `experiments/` (frozen configurations, runner, analysis scripts;
results are archived separately). Neither directory is part of the R
package sources.

## Installation

``` r

install.packages("seqbench")            # CRAN

# development version
# install.packages("pak")
pak::pak("castlaboratory/seqbench")
```

## Related software

`seqbench` is not a general confidence-sequence library. For that see
[confseq](https://github.com/gostevehoward/confseq) (Python),
[safestats](https://cran.r-project.org/package=safestats) and
[avlm](https://cran.r-project.org/package=avlm) (R). The CRAN package
[seqcomp](https://cran.r-project.org/package=seqcomp) compares
probabilistic *forecasts* sequentially; `seqbench` compares *algorithms*
on paired losses with an equivalence margin and explicit cost
accounting.

## License

GPL (\>= 3). © CAST Lab, Universidade Federal de Pernambuco.
