# seqbench 0.1.1

* Betting boundary: decisions are now taken from the exact hedged capital at the two
  decision thresholds, updated in O(1) per instance (`boundary_init(track = )`), instead of
  bisecting the capital over the whole history whenever an endpoint was near a threshold.
  The per-instance cost no longer grows with the length of the stream: a stream of 8,000
  instances at the margin went from 78 s to about 6 s on the same machine. Decisions and
  stopping times are unchanged (identical in 600 checked streams); reported running
  intervals can be up to one grid cell wider per side between calls, and both endpoints
  are still refined once per `update_comparison()` call. States saved with 0.1.0 keep
  working through the previous code path.
* Endpoint refinement of the betting boundary evaluates the exact capital on vectorised
  sub-grids instead of scalar bisection.
* A loss column made only of `NA` (logical) is reported as missing values, not as
  non-numeric.
* `comparison_design()` validates `betting_grid` (previously checked only when the
  comparison was initialised).
* `comparison_report()` states assumption A2 as in the article: the number and the values
  of an instance's seeds are fixed without looking at any of its losses.
* `planning_horizon()` documentation: the Hoeffding horizon is guaranteed for the Hoeffding
  boundary only and can exceed what betting needs by three orders of magnitude; the
  "about half" rule of thumb applies to the `sd` mode (about one half for empirical
  Bernstein, one quarter for betting).
* New vignettes: "Planning a comparison and choosing a boundary" and "Seeds, clusters and
  refusals". The getting-started vignette lists all five boundaries.
* Package site: case studies on vehicle routing heuristics and on vaccination strategies on
  contact networks, a page on scope and limitations, author affiliations, and a fixed web
  manifest (the favicon icons used absolute paths that broke under `/seqbench/`).
* `CITATION` gives the CRAN URL, ORCID identifiers and the companion article.

# seqbench 0.1.0

First CRAN release. Development history (0.0.0.9000 to 0.0.0.9002) is kept below.

* `comparison_design()`, `initialize_comparison()`, `update_comparison()`,
  `stopping_decision()`, `comparison_report()`, `planning_horizon()`; S3 methods
  `print`, `summary`, `tidy`, `glance`, `autoplot` for `seqbench_comparison`.
* Boundaries: `"betting"` (default; Waudby-Smith & Ramdas 2024), `"empirical_bernstein"`,
  `"hoeffding"`, `"bernstein_declared"` (declared `sd_max`, assumption A5) and the invalid
  negative control `"naive_fixed"`, labelled as such in every output.
* Internal kernels `boundary_init()`, `boundary_update()`, `boundary_interval()` exported
  for sister packages.

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
