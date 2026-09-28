# Changelog

## seqbench 0.1.0

First CRAN release. Development history (0.0.0.9000 to 0.0.0.9002) is
kept below.

- [`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md),
  [`initialize_comparison()`](https://castlaboratory.github.io/seqbench/reference/initialize_comparison.md),
  [`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md),
  [`stopping_decision()`](https://castlaboratory.github.io/seqbench/reference/stopping_decision.md),
  [`comparison_report()`](https://castlaboratory.github.io/seqbench/reference/comparison_report.md),
  [`planning_horizon()`](https://castlaboratory.github.io/seqbench/reference/planning_horizon.md);
  S3 methods `print`, `summary`, `tidy`, `glance`, `autoplot` for
  `seqbench_comparison`.
- Boundaries: `"betting"` (default; Waudby-Smith & Ramdas 2024),
  `"empirical_bernstein"`, `"hoeffding"`, `"bernstein_declared"`
  (declared `sd_max`, assumption A5) and the invalid negative control
  `"naive_fixed"`, labelled as such in every output.
- Internal kernels
  [`boundary_init()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md),
  [`boundary_update()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md),
  [`boundary_interval()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md)
  exported for sister packages.

## seqbench 0.0.0.9002

- New boundary `"bernstein_declared"`: predictable-plug-in Bennett
  confidence sequence with a declared upper bound `sd_max` on the
  standard deviation of the paired difference (assumption A5, recorded
  in the report). Valid only if the declared bound holds; its gain over
  betting is modest because the Bennett range term keeps a
  `log(2/alpha)/t` floor for bounded observations.

## seqbench 0.0.0.9001

- Betting boundary: endpoint refinement by bisection is now lazy inside
  [`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md)
  (only when a decision threshold falls in the boundary grid cell, plus
  once per call for reporting). Decisions and stopping times are
  unchanged; per-instance cost is now linear in the grid size instead of
  linear in the number of observations (a 2000-instance stream at the
  margin went from minutes to 0.3 s).
  [`boundary_interval()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md)
  gained a `thresholds` argument.
- Betting boundary: a nonempty confidence set narrower than one grid
  cell is now located by searching the exact capital instead of being
  reported as empty (found by external review; constant streams
  triggered it).
- [`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md)
  re-evaluates the decision after the final endpoint refinement, so an
  empty running intersection is always reported as such.
- [`planning_horizon()`](https://castlaboratory.github.io/seqbench/reference/planning_horizon.md)
  honours `max_t`; `n_max`, `betting_grid` and `refine` are validated as
  whole numbers / logical; costs must be finite.

## seqbench 0.0.0.9000

- Core implemented:
  [`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md),
  [`initialize_comparison()`](https://castlaboratory.github.io/seqbench/reference/initialize_comparison.md),
  [`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md),
  [`stopping_decision()`](https://castlaboratory.github.io/seqbench/reference/stopping_decision.md),
  [`comparison_report()`](https://castlaboratory.github.io/seqbench/reference/comparison_report.md),
  [`planning_horizon()`](https://castlaboratory.github.io/seqbench/reference/planning_horizon.md);
  boundaries `betting` (default), `empirical_bernstein`, `hoeffding`
  (Waudby-Smith & Ramdas 2024) and the invalid negative control
  `naive_fixed`; S3 methods `print`, `summary`, `tidy`, `glance`,
  `autoplot`.
- Internal boundary kernels
  [`boundary_init()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md),
  [`boundary_update()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md),
  [`boundary_interval()`](https://castlaboratory.github.io/seqbench/reference/boundary_kernels.md)
  exported for sister packages.
