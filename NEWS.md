# seqbench 0.0.0.9001

* Betting boundary: endpoint refinement by bisection is now lazy inside
  `update_comparison()` (only when a decision threshold falls in the boundary
  grid cell, plus once per call for reporting). Decisions and stopping times are
  unchanged; per-instance cost is now linear in the grid size instead of linear
  in the number of observations (a 2000-instance stream at the margin went from
  minutes to 0.3 s). `boundary_interval()` gained a `thresholds` argument.

# seqbench 0.0.0.9000

* Core implemented: `comparison_design()`, `initialize_comparison()`,
  `update_comparison()`, `stopping_decision()`, `comparison_report()`,
  `planning_horizon()`; boundaries `betting` (default), `empirical_bernstein`,
  `hoeffding` (Waudby-Smith & Ramdas 2024) and the invalid negative control
  `naive_fixed`; S3 methods `print`, `summary`, `tidy`, `glance`, `autoplot`.
* Internal boundary kernels `boundary_init()`, `boundary_update()`,
  `boundary_interval()` exported for sister packages.
