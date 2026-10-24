# experiments/analysis/apps_summary.R -- analysis of the two applications frozen on
# 2026-10-23 (app_cvrp_2026-10-23, app_epidemic_2026-10-23). Same conventions as
# amendments_summary.R: Wilson intervals for rates, percentile bootstrap over replications
# for medians. The truth of each cell is the mean of the frozen pool (the resampling
# population), so false declarations and coverage failures are measurable.
# Writes to paper/rjournal/data/:
#   T9_apps_pool.csv    pool descriptors (mean, sd, skewness, ties, ICC and design effect)
#   T10_apps_runs.csv   operating characteristics per application, cell and method
#   T11_apps_live.csv   the live run: the pool in generation order, one stream per method
#   F3_apps_live.csv    trajectories of the live betting runs (for the article's figure)
#   T12_apps_models.csv mean difference by instance type (network model / customer placement), seed 1
# Usage: Rscript experiments/analysis/apps_summary.R [paper/rjournal/data]
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1] else "paper/rjournal/data"
if (dir.exists("R") && requireNamespace("devtools", quietly = TRUE)) devtools::load_all(quiet = TRUE) else library(seqbench)
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n; c <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind(lower = (c - h) / d, upper = (c + h) / d)
}
skew <- function(x) mean((x - mean(x))^3) / stats::sd(x)^3
set.seed(1); B <- 1000L

apps <- c(cvrp = "app_cvrp_2026-10-23", epidemic = "app_epidemic_2026-10-23")
model_rows <- list(); pool_rows <- list(); run_rows <- list(); live_rows <- list(); traj_rows <- list()
for (app in names(apps)) {
  cfg <- source(file.path("experiments/configs", paste0(apps[[app]], ".R")))$value
  P <- readRDS(cfg$pool_file)
  has_seeds <- "seeds" %in% names(cfg$dgp)       # epidemic: seed = 1..5; cvrp: seed is the instance's RNG seed
  cells <- if (has_seeds) cfg$dgp$seeds else 1L
  P$seed_index <- if (has_seeds) P$seed else 1L
  for (m in unique(cells)) {
    Pm <- P[P$seed_index <= m, ]
    d_inst <- tapply(Pm$loss_a - Pm$loss_b, Pm$instance, mean)
    icc <- NA_real_; deff <- NA_real_
    if (m > 1) {                                  # one-way ANOVA estimator on the row differences
      dd <- Pm$loss_a - Pm$loss_b; g <- factor(Pm$instance); k <- nlevels(g)
      msb <- m * sum((tapply(dd, g, mean) - mean(dd))^2) / (k - 1)
      msw <- sum((dd - ave(dd, g))^2) / (k * (m - 1))
      icc <- (msb - msw) / (msb + (m - 1) * msw); deff <- 1 + (m - 1) * icc
    }
    pool_rows[[length(pool_rows) + 1L]] <- data.frame(
      app = app, seeds = m, instances = length(d_inst), mu = mean(d_inst), sd = stats::sd(d_inst),
      skewness = skew(d_inst), ties = mean(d_inst == 0), mean_loss_a = mean(Pm$loss_a),
      mean_loss_b = mean(Pm$loss_b), icc = icc, design_effect = deff,
      cpu_per_instance = mean(tapply(Pm$cost, Pm$instance, sum)))
  }
  grp <- intersect(c("model", "placement"), names(P))[1]       # instance type, one seed
  P1 <- P[P$seed_index == 1L, ]
  for (lv in sort(unique(P1[[grp]]))) {
    d1 <- P1[P1[[grp]] == lv, ]
    model_rows[[length(model_rows) + 1L]] <- data.frame(app = app, model = lv, instances = nrow(d1),
      mu = mean(d1$loss_a - d1$loss_b), sd = stats::sd(d1$loss_a - d1$loss_b))
  }
  runs <- read.csv(file.path("experiments/out", cfg$name, "runs.csv"), stringsAsFactors = FALSE)
  stopifnot(!any(runs$decision == "error"))
  keys <- intersect(c("seeds"), names(runs))
  for (d in split(runs, runs[c("method", keys)], drop = TRUE)) {
    n <- nrow(d)
    med_b <- replicate(B, stats::median(d$n_instances[sample.int(n, n, replace = TRUE)]))
    fd <- sum(d$false_declaration, na.rm = TRUE); mc <- sum(d$miscover, na.rm = TRUE)
    run_rows[[length(run_rows) + 1L]] <- data.frame(
      app = app, method = d$method[1], seeds = if (length(keys)) d$seeds[1] else 1L, R = n,
      mu = d$mu_true[1], truth = d$truth[1],
      p_A = mean(d$decision == "A"), p_B = mean(d$decision == "B"),
      p_equivalent = mean(d$decision == "equivalent"), p_inconclusive = mean(d$inconclusive),
      false_declaration = fd, false_rate = fd / n, false_hi = wilson(fd, n)[, "upper"],
      miscover = mc, miscover_rate = mc / n, miscover_lo = wilson(mc, n)[, "lower"], miscover_hi = wilson(mc, n)[, "upper"],
      median_n = stats::median(d$n_instances), median_lo = stats::quantile(med_b, 0.025, names = FALSE),
      median_hi = stats::quantile(med_b, 0.975, names = FALSE), mean_n = mean(d$n_instances),
      mean_cost = mean(d$cost))
  }
  # Live run: the pool in generation order (fresh i.i.d. instances), betting and EB
  for (m in unique(cells)) {
    Pm <- P[P$seed_index <= m, ]; Pm <- Pm[order(Pm$instance, Pm$seed_index), ]
    Pm <- Pm[Pm$instance <= cfg$n_max, ]
    for (b in c("betting", "empirical_bernstein")) {
      des <- comparison_design(alpha = cfg$alpha, margin = cfg$margin, bounds = cfg$bounds,
                               boundary = b, n_max = cfg$n_max)
      st <- suppressWarnings(update_comparison(initialize_comparison(des),
              Pm[, c("instance", "seed", "loss_a", "loss_b", "cost")]))
      g <- glance(st)
      live_rows[[length(live_rows) + 1L]] <- data.frame(
        app = app, seeds = m, boundary = b, n_instances = g$n_instances, n_evaluations = g$n_evaluations,
        cpu_seconds = g$cost, estimate = g$estimate, lower = g$lower, upper = g$upper,
        decision = g$decision, reason = g$stopping_reason)
      if (b == "betting") traj_rows[[length(traj_rows) + 1L]] <- cbind(app = app, seeds = m,
        tidy(st)[, c("t", "estimate", "lower_running", "upper_running", "cost_cum")])
    }
  }
}
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(do.call(rbind, pool_rows), file.path(data_dir, "T9_apps_pool.csv"), row.names = FALSE)
write.csv(do.call(rbind, run_rows), file.path(data_dir, "T10_apps_runs.csv"), row.names = FALSE)
write.csv(do.call(rbind, live_rows), file.path(data_dir, "T11_apps_live.csv"), row.names = FALSE)
write.csv(do.call(rbind, traj_rows), file.path(data_dir, "F3_apps_live.csv"), row.names = FALSE)
write.csv(do.call(rbind, model_rows), file.path(data_dir, "T12_apps_models.csv"), row.names = FALSE)
print(do.call(rbind, model_rows), digits = 3)
print(do.call(rbind, pool_rows), digits = 3)
print(do.call(rbind, run_rows)[, c("app", "method", "seeds", "truth", "false_rate", "miscover_rate", "p_inconclusive", "median_n", "mean_cost")], digits = 3)
print(do.call(rbind, live_rows), digits = 3)
