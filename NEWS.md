# seqbench 0.0.0.9002

* New boundary `"bernstein_declared"`: predictable-plug-in Bennett confidence
  sequence with a declared upper bound `sd_max` on the standard deviation of the
  paired difference (assumption A5, recorded in the report). Valid only if the
  declared bound holds; its gain over betting is modest because the Bennett range
  term keeps a `log(2/alpha)/t` floor for bounded observations.

# seqbench 0.0.0.9001

* Betting boundary: endpoint refinement by bisection is now lazy inside
  `update_comparison()` (only when a decision threshold falls in the boundary
  grid cell, plus once per call for reporting). Decisions and stopping times are
  unchanged; per-instance cost is now linear in the grid size instead of linear
  in the number of observations (a 2000-instance stream at the margin went from
  minutes to 0.3 s). `boundary_interval()` gained a `thresholds` argument.
* Betting boundary: a nonempty confidence set narrower than one grid cell is now
  located by searching the exact capital instead of being reported as empty
  (found by external review; constant streams triggered it).
* `update_comparison()` re-evaluates the decision after the final endpoint
  refinement, so an empty running intersection is always reported as such.
* `planning_horizon()` honours `max_t`; `n_max`, `betting_grid` and `refine`
  are validated as whole numbers / logical; costs must be finite.

# seqbench 0.0.0.9000

* Core implemented: `comparison_design()`, `initialize_comparison()`,
  `update_comparison()`, `stopping_decision()`, `comparison_report()`,
  `planning_horizon()`; boundaries `betting` (default), `empirical_bernstein`,
  `hoeffding` (Waudby-Smith & Ramdas 2024) and the invalid negative control
  `naive_fixed`; S3 methods `print`, `summary`, `tidy`, `glance`, `autoplot`.
* Internal boundary kernels `boundary_init()`, `boundary_update()`,
  `boundary_interval()` exported for sister packages.
