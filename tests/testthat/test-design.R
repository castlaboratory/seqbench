test_that("comparison_design validates its arguments", {
  expect_s3_class(comparison_design(margin = 0.02, bounds = c(0, 1)), "seqbench_design")
  expect_error(comparison_design(alpha = 0, margin = 0.02, bounds = c(0, 1)), "alpha")
  expect_error(comparison_design(margin = -1, bounds = c(0, 1)), "margin")
  expect_error(comparison_design(margin = 0.02, bounds = c(1, 0)), "bounds")
  expect_error(comparison_design(margin = 0.02, bounds = c(0, Inf)), "bounds")
  expect_error(comparison_design(margin = 3, bounds = c(0, 1)), "smaller than the range")
  expect_error(comparison_design(margin = 0.02, bounds = c(0, 1), paired = FALSE), "not implemented")
  expect_error(comparison_design(margin = 0.02, bounds = c(0, 1), boundary = "x"), "must be one of")
  expect_error(comparison_design(margin = 0.02, bounds = c(0, 1), budget = 0), "budget")
})

test_that("design records derived quantities and prints", {
  d <- comparison_design(margin = 0.1, bounds = c(-2, 3), boundary = "naive_fixed")
  expect_equal(d$diff_bounds, c(-5, 5))
  expect_false(d$valid)
  expect_message(print(d), "naive_fixed")
})

test_that("planning_horizon is monotone and finite", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1))
  h1 <- planning_horizon(d, distance = 0.2)
  h2 <- planning_horizon(d, distance = 0.1)
  expect_true(is.finite(h1) && h1 < h2)
  # at t*, twice the Hoeffding half-width on the [0,1] scale is below the distance
  w <- seqbench:::hoeffding_half_width(h2, 0.05)
  expect_lt(2 * w[h2], 0.1 / 2)
  expect_gte(2 * w[h2 - 1], 0.1 / 2)
  expect_error(planning_horizon(d, distance = 0), "distance")
})

test_that("variance-based planning horizon is far below the Hoeffding one and increases with sd", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1))
  h_guar <- planning_horizon(d, distance = 0.03)
  h_lo <- planning_horizon(d, distance = 0.03, sd = 0.03)
  h_hi <- planning_horizon(d, distance = 0.03, sd = 0.10)
  expect_true(h_lo < h_hi && h_hi < h_guar)
  expect_error(planning_horizon(d, distance = 0.03, sd = 5), "exceeds")
})

test_that("planning_horizon respects max_t and validates counts", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1))
  h <- planning_horizon(d, distance = 0.9)
  expect_true(is.finite(h) && h > 1)
  expect_equal(planning_horizon(d, distance = 0.9, max_t = h), h)
  expect_true(planning_horizon(d, distance = 0.9, max_t = h - 1) == Inf)
  expect_true(planning_horizon(d, distance = 0.001, max_t = 100) == Inf)
  expect_error(planning_horizon(d, distance = 0.1, max_t = 0), "max_t")
  expect_error(comparison_design(margin = 0.02, bounds = c(0, 1), n_max = 10.5), "whole number")
})
