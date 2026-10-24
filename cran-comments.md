# seqbench 0.1.1

Update of a package first released as 0.1.0.

* The betting boundary now takes its decisions from the exact capital tracked at the two
  decision thresholds. Decisions are unchanged; the per-instance cost no longer grows with
  the length of the stream (long streams were superlinear in 0.1.0).
* Argument validation, documentation corrections and two new vignettes.

## Test environments

* local: Ubuntu, R 4.6.1, `R CMD check --as-cran`
* CI (GitHub Actions, R-CMD-check.yaml): macOS-latest (R release),
  windows-latest (R release), ubuntu-latest (R devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewers

* The package implements confidence sequences for the mean of bounded paired
  differences (Waudby-Smith & Ramdas, 2024, JRSS-B, cited in the documentation)
  and a three-way decision rule with a pre-registered equivalence margin. One
  boundary, `"naive_fixed"`, is intentionally invalid under optional stopping;
  it exists as a negative control for simulation studies and every output
  produced with it carries an explicit "invalid" label.
* Examples, tests and vignettes run in well under a minute; no example is wrapped in
  `\dontrun{}`.
* The companion article (submitted to The R Journal) and its experiments live in
  `paper/` and `experiments/` of the GitHub repository and are excluded from the
  package build by `.Rbuildignore`.
