# seqbench 0.1.0

First submission.

## Test environments

* local: macOS 26 (Apple Silicon), R 4.6, `R CMD check --as-cran`
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
* Examples and tests run in well under a minute; no example is wrapped in
  `\dontrun{}`.
* The companion article (The R Journal, in preparation) and its experiments
  live in a separate repository (castlaboratory/seqbench-paper), not in the
  package.
