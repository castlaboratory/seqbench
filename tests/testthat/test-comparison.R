# shift > 0 makes B lose more, so A is better.
dec <- function(st) as.vector(stopping_decision(st))

make_losses <- function(n, shift, sd = 0.05, seeds = 1L, start = 1L, seed = 1) {
  set.seed(seed)
  inst <- rep(seq.int(start, length.out = n), each = seeds)
  la <- pmin(pmax(rbeta(length(inst), 2, 5), 0), 1)
  lb <- pmin(pmax(la + shift + rnorm(length(inst), 0, sd), 0), 1)
  data.frame(instance = inst, seed = rep(seq_len(seeds), times = n), loss_a = la, loss_b = lb)
}

design <- comparison_design(alpha = 0.05, margin = 0.02, bounds = c(0, 1))

test_that("initialize_comparison sets the full-range interval and continue", {
  st <- initialize_comparison(design)
  expect_s3_class(st, "seqbench_comparison")
  expect_equal(dec(st), "continue")
  expect_equal(nrow(tidy(st)), 0L)
  g <- glance(st)
  expect_equal(c(g$lower, g$upper), c(-1, 1))
  expect_no_error(print(st))
  expect_no_error(print(comparison_report(st)))
})

test_that("a clear advantage of A is declared and stops the stream", {
  st <- initialize_comparison(design)
  expect_warning(st <- update_comparison(st, make_losses(400, shift = 0.15)), "not consumed")
  expect_equal(dec(st), "A")
  expect_equal(attr(stopping_decision(st), "reason"), "declaration")
  tr <- tidy(st)
  expect_lt(tr$upper_running[nrow(tr)], -0.02)
  expect_true(all(tr$decision[-nrow(tr)] == "continue"))
  expect_true(st$stopped)
  expect_error(update_comparison(st, make_losses(5, 0.1, start = 1000)), "already stopped")
})

test_that("a clear advantage of B is declared", {
  st <- update_comparison(initialize_comparison(design), make_losses(400, shift = -0.15)) |>
    suppressWarnings()
  expect_equal(dec(st), "B")
})

test_that("true equivalence is declared with a tight enough boundary", {
  d <- comparison_design(alpha = 0.05, margin = 0.05, bounds = c(0, 1))
  st <- update_comparison(initialize_comparison(d), make_losses(3000, shift = 0, sd = 0.01)) |>
    suppressWarnings()
  expect_equal(dec(st), "equivalent")
  tr <- tidy(st)
  expect_true(tr$lower_running[nrow(tr)] >= -0.05 && tr$upper_running[nrow(tr)] <= 0.05)
})

test_that("budget and n_max produce inconclusive outcomes with the right reason", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1), budget = 10, cost_per_round = 1)
  st <- update_comparison(initialize_comparison(d), make_losses(10, shift = 0.01, sd = 0.2)) |>
    suppressWarnings()
  expect_equal(dec(st), "inconclusive")
  expect_equal(attr(stopping_decision(st), "reason"), "budget")
  expect_equal(st$cost, 10)
  d2 <- comparison_design(margin = 0.02, bounds = c(0, 1), n_max = 7)
  st2 <- update_comparison(initialize_comparison(d2), make_losses(20, shift = 0.01, sd = 0.2)) |>
    suppressWarnings()
  expect_equal(attr(stopping_decision(st2), "reason"), "n_max")
  expect_equal(st2$diagnostics$n_instances, 7L)
})

test_that("seeds are averaged per instance and costs are summed", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1), boundary = "hoeffding", cost_per_round = 2)
  l <- make_losses(5, shift = 0.1, seeds = 3L)
  st <- update_comparison(initialize_comparison(d), l)
  tr <- tidy(st)
  expect_equal(tr$n_seeds, rep(3L, 5))
  expect_equal(st$diagnostics$n_evaluations, 15L)
  expect_equal(tr$cost, rep(6, 5))
  expect_equal(tr$x[1], (mean(l$loss_a[1:3] - l$loss_b[1:3]) + 1) / 2)
  # explicit cost column overrides
  l$cost <- 0.5
  st2 <- update_comparison(initialize_comparison(d), l)
  expect_equal(st2$cost, 7.5)
})

test_that("invalid input is rejected, never converted", {
  st <- initialize_comparison(design)
  base <- make_losses(5, 0.1)
  expect_error(update_comparison(st, base[, c("instance", "loss_a")]), "missing column")
  expect_error(update_comparison(st, base[0, ]), "no rows")
  bad <- base; bad$loss_a[2] <- NA
  expect_error(update_comparison(st, bad), "missing value")
  bad <- base; bad$loss_b[3] <- 1.2
  expect_error(update_comparison(st, bad), "outside the declared bounds")
  bad <- base; bad$cost <- -1
  expect_error(update_comparison(st, bad), "non-negative")
  bad <- rbind(base, base[1, ])
  expect_error(update_comparison(st, bad), "Duplicated")
  st <- update_comparison(st, base)
  expect_error(update_comparison(st, make_losses(3, 0.1, start = 3)), "already observed")
  expect_error(update_comparison(list(), base), "seqbench_comparison")
})

test_that("naive_fixed comparisons are flagged invalid everywhere", {
  d <- comparison_design(margin = 0.02, bounds = c(0, 1), boundary = "naive_fixed")
  st <- update_comparison(initialize_comparison(d), make_losses(30, 0.1)) |> suppressWarnings()
  expect_false(glance(st)$valid)
  expect_false(comparison_report(st)$valid)
  expect_message(print(st), "NOT a confidence sequence")
})

test_that("report carries the full contract", {
  st <- update_comparison(initialize_comparison(design), make_losses(20, 0.05, sd = 0.2)) |>
    suppressWarnings()
  r <- comparison_report(st)
  expect_named(r, c("estimand", "estimate", "confidence_sequence", "alpha", "margin",
                    "boundary", "valid", "assumptions", "decision", "decision_label",
                    "stopping_reason", "diagnostics", "resources", "instances", "seeds",
                    "versions"))
  expect_length(r$assumptions, 4L)
  expect_equal(r$versions$seqbench_version, as.character(packageVersion("seqbench")))
  expect_identical(summary(st)$decision, r$decision)
})

test_that("all valid boundaries agree on a strong effect", {
  for (b in c("betting", "empirical_bernstein", "hoeffding")) {
    d <- comparison_design(margin = 0.02, bounds = c(0, 1), boundary = b)
    st <- update_comparison(initialize_comparison(d), make_losses(2000, shift = 0.2, sd = 0.05)) |>
      suppressWarnings()
    expect_equal(dec(st), "A", info = b)
  }
})

test_that("autoplot returns a ggplot and refuses an empty state", {
  st <- initialize_comparison(design)
  expect_error(ggplot2::autoplot(st), "Nothing to plot")
  st <- update_comparison(st, make_losses(10, 0.05, sd = 0.2)) |> suppressWarnings()
  expect_s3_class(ggplot2::autoplot(st), "ggplot")
})

test_that("lazy refinement of the betting boundary reproduces full refinement exactly", {
  # Reference: full refinement at every step, replayed with the kernels directly.
  replay_full <- function(x01, margin01, alpha = 0.05) {
    st <- boundary_init("betting", alpha = alpha); L <- 0; U <- 1; dec <- "continue"; t <- 0L
    for (x in x01) {
      t <- t + 1L; st <- boundary_update(st, x); ci <- boundary_interval(st)
      L <- max(L, ci[["lower"]]); U <- min(U, ci[["upper"]])
      if (U < 0.5 - margin01) { dec <- "A"; break }
      if (L > 0.5 + margin01) { dec <- "B"; break }
      if (L >= 0.5 - margin01 && U <= 0.5 + margin01) { dec <- "equivalent"; break }
    }
    list(decision = if (dec == "continue") "inconclusive" else dec, t = t, L = L, U = U)
  }
  d <- comparison_design(alpha = 0.05, margin = 0.02, bounds = c(0, 1), n_max = 1500)
  set.seed(11)
  for (shift in c(0.06, -0.06, 0, 0.015)) {
    l <- make_losses(1500, shift = shift, sd = 0.1, seed = 100 + round(shift * 1000))
    st <- suppressWarnings(update_comparison(initialize_comparison(d), l))
    x01 <- (l$loss_a - l$loss_b + 1) / 2
    ref <- replay_full(x01[seq_len(nrow(st$trajectory))], margin01 = 0.02 / 2)
    expect_identical(dec(st), ref$decision, info = paste("shift", shift))
    expect_identical(nrow(st$trajectory), ref$t, info = paste("shift", shift))
    # The reported running intersection contains the fully refined one and is at
    # most one grid cell (1/1000 on [0, 1], i.e. 2/1000 in loss units) wider per side.
    tr <- tidy(st)
    lo_full <- -1 + 2 * ref$L; hi_full <- -1 + 2 * ref$U
    expect_lte(tr$lower_running[nrow(tr)], lo_full + 1e-8)
    expect_gte(tr$upper_running[nrow(tr)], hi_full - 1e-8)
    expect_lte(lo_full - tr$lower_running[nrow(tr)], 2 / 1000 + 1e-8)
    expect_lte(tr$upper_running[nrow(tr)] - hi_full, 2 / 1000 + 1e-8)
  }
})

test_that("betting boundary with thresholds only refines when a threshold is in the cell", {
  set.seed(3); x <- rbeta(300, 10, 30)
  st <- boundary_init("betting", alpha = 0.05, grid = 401L)
  for (v in x) st <- boundary_update(st, v)
  full <- boundary_interval(st)
  far <- boundary_interval(st, thresholds = c(0.9, 0.95))    # no threshold near the interval
  raw <- boundary_interval(boundary_init("betting", alpha = 0.05, grid = 401L, refine = FALSE) |>
                             (\(s) { for (v in x) s <- boundary_update(s, v); s })())
  expect_equal(far, raw)                                    # outer cells, unrefined
  expect_true(full[["lower"]] >= raw[["lower"]] && full[["upper"]] <= raw[["upper"]])
  near <- boundary_interval(st, thresholds = c(full[["lower"]] + 1e-6, 0.99))
  expect_equal(near[["lower"]], full[["lower"]])            # lower refined, upper not
  expect_equal(near[["upper"]], raw[["upper"]])
})
