# Boundary kernels (internal API)

Low-level state machines behind the confidence sequences used by
[`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md).
Observations must already be rescaled to `[0, 1]`. These functions are
exported for reuse by sister packages; ordinary users should call
[`comparison_design()`](https://castlaboratory.github.io/seqbench/reference/comparison_design.md)
and
[`update_comparison()`](https://castlaboratory.github.io/seqbench/reference/update_comparison.md)
instead.

## Usage

``` r
boundary_init(
  name = seqbench_boundaries(),
  alpha = 0.05,
  c = 0.5,
  theta = 0.5,
  grid = 1001L,
  refine = TRUE,
  sd_max01 = NULL
)

boundary_update(state, x)

boundary_interval(state, thresholds = NULL)
```

## Arguments

- name:

  One of
  [`seqbench_boundaries()`](https://castlaboratory.github.io/seqbench/reference/seqbench_boundaries.md).

- alpha:

  Miscoverage level in (0, 1).

- c:

  Truncation constant for the empirical-Bernstein and betting bets, in
  (0, 1). Waudby-Smith & Ramdas recommend 1/2 or 3/4.

- theta:

  Hedging weight on the "long" capital in the betting boundary, in (0,
  1). Default 1/2.

- grid:

  Number of grid points on `[0, 1]` for the betting boundary.

- refine:

  Logical; refine the betting interval endpoints by bisection on the
  exact capital (default `TRUE`). Without refinement the endpoints are
  the outer boundaries of the grid cells that contain the true
  endpoints, which is conservative and keeps the coverage guarantee.

- sd_max01:

  For `"bernstein_declared"`: declared upper bound on the conditional
  standard deviation of the observations, on the `[0, 1]` scale (at most
  1/2). The boundary is a predictable-plug-in Bennett confidence
  sequence: for `X in [0, 1]` with conditional mean `mu` and conditional
  variance at most `sd_max01^2`,
  `exp(lambda (X - mu) - sd_max01^2 (e^lambda - 1 - lambda))` is a
  supermartingale for every predictable `lambda >= 0` (Bennett's
  inequality), and the same holds for `-(X - mu)`. Coverage is
  guaranteed only if the declared bound is true.

- state:

  A boundary state returned by `boundary_init()` or `boundary_update()`.

- x:

  A single number in `[0, 1]`.

- thresholds:

  Optional numeric vector of decision thresholds on the `[0, 1]` scale.
  When given to `boundary_interval()` for the betting boundary, an
  endpoint is refined only if its grid cell contains one of them, which
  is the only case in which refinement can change a decision; this keeps
  the per-step cost linear in the grid size instead of linear in the
  number of observations. Decisions and stopping times are identical to
  full refinement; the reported running intersection can be up to one
  grid cell wider per side. `NULL` (default) refines both endpoints.

## Value

`boundary_init()` and `boundary_update()` return a state object of class
`seqbench_boundary`. `boundary_interval()` returns a named numeric
vector `c(estimate, lower, upper)` on the `[0, 1]` scale. For the
betting boundary, `NA` endpoints mean that the exact confidence set is
empty (a possible event, of probability at most `alpha` under the
assumptions); a nonempty set narrower than one grid cell is located by
searching the exact capital, never reported as empty.

## Examples

``` r
st <- boundary_init("hoeffding", alpha = 0.05)
for (x in c(0.2, 0.3, 0.25)) st <- boundary_update(st, x)
boundary_interval(st)
#> estimate    lower    upper 
#>     0.25     0.00     1.00 
```
