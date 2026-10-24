# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`seqbench` is an R package (on CRAN, 0.1.0) for anytime-valid sequential comparison of two algorithms on paired losses with a practical-equivalence margin. The repo also holds the companion R Journal article (`paper/rjournal/`) and its simulation study (`experiments/`); neither is part of the package build (both are in `.Rbuildignore`). Note that `PLANO.md`, `memory/` and `references/` are gitignored local-only material, so comments that cite `references/wsr2024_fronteiras.md` or `PLANO.md` point to files a fresh clone lacks. The state of the October 2026 review (open items, release date pending) is in `dev/REVIEW-2026-10.md`.

## Commands

Standard devtools workflow from the repo root:

```r
devtools::load_all()
devtools::document()                         # regenerate man/ and NAMESPACE (roxygen2, markdown on); never edit NAMESPACE by hand
devtools::test()                             # all tests (testthat edition 3)
devtools::test(filter = "boundaries")        # one file: tests/testthat/test-boundaries.R
testthat::test_file("tests/testthat/test-comparison.R")
devtools::check()                            # or: R CMD build . && R CMD check --as-cran seqbench_*.tar.gz
pkgdown::build_site()                        # site config in _pkgdown.yml; CI deploys it
```

CI (`.github/workflows/`) runs R-CMD-check on macOS/Windows/Ubuntu (devel, release, oldrel-1), test coverage and pkgdown. CRAN-facing claims (examples run in under a minute, none in `\dontrun{}`) are in `cran-comments.md`; keep them true. Slow tests use `skip_on_cran()`.

Experiments (resumable, outputs to gitignored `experiments/out/<study>/`):

```sh
Rscript experiments/run/run_study.R experiments/configs/<config>.R 8      # 8 cores
Rscript experiments/run/run_study.R experiments/configs/<config>.R 4 3    # smoke test, 3 reps
SEQBENCH_METHODS="^savi" Rscript experiments/run/run_study.R ...          # split a run across machines
```

Application studies (2026-10-23 onwards) first build a frozen instance pool, then resample streams from it so the truth (pool mean) is known:

```sh
Rscript experiments/run/make_pool.R experiments/configs/app_cvrp_2026-10-23.R 4 [n_pilot]
Rscript experiments/run/run_study.R experiments/configs/app_cvrp_2026-10-23.R 4
Rscript experiments/analysis/apps_summary.R          # writes paper/rjournal/data/T9-T12, F3
```

The paper renders with `rmarkdown::render("paper/rjournal/seqbench.Rmd")` (needs `rjtools`, `kableExtra`). It reads summary tables from `paper/rjournal/data/*.csv` (written by `experiments/analysis/*.R`), and numbers in the prose are inline R computed from those CSVs; keep it that way rather than typing numbers. Abstract must stay under 150 words (R Journal). Refer to sections with explicit links, `[Section 5](#sec-sim)`, never `\@ref(sec-...)`: with rjtools the latter prints "2.5" in the PDF (report class) and "??" in the HTML; subsections are unnumbered, so name them instead. `\@ref(tab:...)`/`\@ref(fig:...)` work. Claims about timings must be measured on an idle machine.

## Architecture

Two layers, deliberately separated:

1. **Boundary kernels** (`R/boundaries.R`): small S3 state machines, one class per boundary (`seqbench_boundary_<name>`), with `boundary_init()` / `boundary_update(state, x)` / `boundary_interval(state, thresholds)`. They only know observations already rescaled to `[0, 1]`. Exported but `@keywords internal` so sister packages can reuse them. Adding a boundary means: add its name to `seqbench_boundaries()`, an init branch, and `boundary_update`/`boundary_interval` methods. Formulas follow Waudby-Smith & Ramdas (2024); `betting` is grid-based with bisection refinement of endpoints (only near decision thresholds when `thresholds` is given). `naive_fixed` is an intentionally invalid negative control and is flagged `valid = FALSE` everywhere it surfaces.

2. **Comparison layer**: `comparison_design()` (`R/design.R`, immutable, pre-registered parameters, plus `planning_horizon()`) → `initialize_comparison()` → `update_comparison()` (`R/comparison.R`). The update loop maps loss differences from `diff_bounds` onto `[0, 1]`, averages all seeds of an instance into one observation (instances, not seeds, are the sequential unit), maintains the **running intersection** of intervals, and calls `decide()` (thresholds `0.5 ± margin/range` on the rescaled scale). It stops at the first terminal decision, warns about unconsumed rows, and refuses further updates once stopped. For `betting`, decisions come from the exact hedged capital tracked at the two thresholds (`boundary_init(track=)`, `state$sides`, `update_sides()`/`clamp_sides()`), never from grid endpoints, so per-instance cost is O(grid) regardless of stream length; endpoints are refined once per call for reporting only. States without `track_m` (0.1.0) fall back to lazy bisection refinement.

Outputs: `R/report.R` (`print`, `summary`, `comparison_report()` — the "result contract": estimand, assumptions, diagnostics, cost, versions), `R/broom.R` (`tidy` = per-instance trajectory, `glance`), `R/autoplot.R`.

Design principle enforced in `R/validate.R` and `R/utils-validation.R`: invalid input is an error, never silently fixed. Losses outside `bounds` are not clipped, NAs and negative costs abort, and an instance that reappears after being seen aborts. Errors use `cli::cli_abort()` with explanatory `"i"` bullets; follow that style.

Tests in `test-boundaries.R` check kernels against independent vectorised oracles and against fixed numbers from `dev/wsr_cs.R` (a standalone reference implementation; `dev/repro_*.R` compare against it and `safestats`).

## Experiments conventions

- `experiments/configs/*_<date>.R` are **frozen** pre-registered configurations: never edit one after it is frozen; add a new dated config (an "amendment") instead. `experiments/README.md` lists each study and what it covers.
- `run_study.R` hashes the config and `run/methods.R`, and a resumed run refuses cells produced under a different hash, so editing `methods.R` invalidates cached cells.
- Every replication's data stream is seeded only from `(seed, dgp cell, rep)` so all methods see identical paired losses.
