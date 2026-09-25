# seqbench 0.0.0.9000

* Core implemented: `comparison_design()`, `initialize_comparison()`,
  `update_comparison()`, `stopping_decision()`, `comparison_report()`,
  `planning_horizon()`; boundaries `betting` (default), `empirical_bernstein`,
  `hoeffding` (Waudby-Smith & Ramdas 2024) and the invalid negative control
  `naive_fixed`; S3 methods `print`, `summary`, `tidy`, `glance`, `autoplot`.
* Internal boundary kernels `boundary_init()`, `boundary_update()`,
  `boundary_interval()` exported for sister packages.
