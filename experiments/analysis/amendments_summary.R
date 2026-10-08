# experiments/analysis/amendments_summary.R -- analysis of the two amendments of 2026-09-26
# (main_sim_extra_2026-09-26: negative and near-margin effects; main_sim_cluster_2026-09-26:
# seed clusters with an instance-specific differential effect and the invalid
# pseudoreplication analysis). Written 2026-09-28 after the runs finished; it reuses the
# metrics and the uncertainty conventions of the pre-registered main_summary.R (Wilson
# intervals for rates, percentile bootstrap over replications for medians, paired bootstrap
# for ratios) and adds only what the amendments need.
#
# Cells with decision == "error" are DROPPED and listed (main_summary.R refuses them). In the
# extra study, cell dgp 8 (mu = -0.20, h = 0.173) is invalid by construction: loss_B =
# u_i + 0.20 + e reaches 1.073 > 1 and the package refuses out-of-bounds losses (never clips).
# Usage: Rscript experiments/analysis/amendments_summary.R [paper/rjournal/data]
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1] else "paper/rjournal/data"
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n; c <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind(lower = (c - h) / d, upper = (c + h) / d)
}
set.seed(1); B <- 1000L
load_study <- function(dir) {
  runs <- read.csv(file.path(dir, "runs.csv"), stringsAsFactors = FALSE)
  cfg <- readRDS(file.path(dir, "config.rds"))
  err <- runs[runs$decision == "error", ]
  if (nrow(err)) {
    bad <- unique(err[c("dgp", "mu", "h", "method")])
    message(sprintf("%s: dropping %d error rows in dgp cell(s) %s; first reason: %s", cfg$name, nrow(err),
                    paste(unique(bad$dgp), collapse = ","), strsplit(err$reason[1], "\n")[[1]][1]))
    runs <- runs[runs$decision != "error", ]
  }
  stopifnot(!anyDuplicated(runs[c("dgp", "rep", "method")]))
  n_per <- table(runs$dgp, runs$method); stopifnot(all(n_per[n_per > 0] == cfg$R))
  list(runs = runs, cfg = cfg, dropped = if (nrow(err)) unique(err[c("dgp", "mu", "h")]) else NULL)
}
summarise_cells <- function(runs, cfg, keys) {
  # validity (counts + Wilson) and cost (median with bootstrap CI, p_inconclusive) per cell
  do.call(rbind, lapply(split(runs, runs[c("method", keys)], drop = TRUE), function(d) {
    d <- d[order(d$rep), ]; n <- nrow(d)
    idx <- replicate(B, sample.int(n, n, replace = TRUE), simplify = FALSE)
    med_b <- vapply(idx, function(i) median(d$n_instances[i]), numeric(1))
    out <- d[1, c("method", keys), drop = FALSE]
    out$truth <- d$truth[1]; out$R <- n
    out$false_declaration <- sum(d$false_declaration); out$false_rate <- out$false_declaration / n
    out[c("false_lo", "false_hi")] <- wilson(out$false_declaration, n)
    out$miscover <- sum(d$miscover); out$miscover_rate <- out$miscover / n
    out[c("miscover_lo", "miscover_hi")] <- wilson(out$miscover, n)
    out$p_inconclusive <- mean(d$inconclusive)
    out$median_n <- median(d$n_instances); out$median_lo <- unname(quantile(med_b, 0.025)); out$median_hi <- unname(quantile(med_b, 0.975))
    out$mean_n <- mean(d$n_instances); out$mean_cost <- mean(d$cost)
    out
  }))
}
paired_ratio <- function(d_m, d_b, f = median, col = "n_instances") {
  d_m <- d_m[order(d_m$rep), ]; d_b <- d_b[order(d_b$rep), ]; stopifnot(identical(d_m$rep, d_b$rep)); n <- nrow(d_m)
  r <- replicate(B, { i <- sample.int(n, n, replace = TRUE); f(d_m[[col]][i]) / f(d_b[[col]][i]) })
  c(ratio = f(d_m[[col]]) / f(d_b[[col]]), unname(quantile(r, c(0.025, 0.975))))
}

# Extra scenarios -----------------------------------------------------------------
ex <- load_study("experiments/out/main_sim_extra_2026-09-26")
t_ex <- summarise_cells(ex$runs, ex$cfg, c("mu", "h"))
t_ex$ratio_to_betting <- NA_real_; t_ex$ratio_lo <- NA_real_; t_ex$ratio_hi <- NA_real_
for (i in seq_len(nrow(t_ex))) {
  if (t_ex$method[i] == "seqbench:betting") next
  sel <- ex$runs$mu == t_ex$mu[i] & ex$runs$h == t_ex$h[i]
  if (!any(sel & ex$runs$method == "seqbench:betting")) next   # betting cell dropped (dgp 8)
  rr <- paired_ratio(ex$runs[sel & ex$runs$method == t_ex$method[i], ], ex$runs[sel & ex$runs$method == "seqbench:betting", ])
  t_ex[i, c("ratio_to_betting", "ratio_lo", "ratio_hi")] <- rr
}
t_ex <- t_ex[order(t_ex$h, t_ex$mu, t_ex$method), ]
write.csv(t_ex, file.path(data_dir, "T7_extra.csv"), row.names = FALSE)

# Cluster v2 -----------------------------------------------------------------------
cl <- load_study("experiments/out/main_sim_cluster_2026-09-26")
t_cl <- summarise_cells(cl$runs, cl$cfg, c("mu", "seeds", "w"))
# cost of five seeds relative to one seed, same method, in instances and in evaluations (paired reps)
t_cl$ratio_n_vs_1seed <- NA_real_; t_cl$ratio_n_lo <- NA_real_; t_cl$ratio_n_hi <- NA_real_; t_cl$ratio_cost_vs_1seed <- NA_real_
for (i in seq_len(nrow(t_cl))) {
  if (t_cl$seeds[i] == 1L) next
  sel <- cl$runs$mu == t_cl$mu[i] & cl$runs$w == t_cl$w[i] & cl$runs$method == t_cl$method[i]
  rr <- paired_ratio(cl$runs[sel & cl$runs$seeds == 5L, ], cl$runs[sel & cl$runs$seeds == 1L, ])
  t_cl[i, c("ratio_n_vs_1seed", "ratio_n_lo", "ratio_n_hi")] <- rr
  t_cl$ratio_cost_vs_1seed[i] <- mean(cl$runs$cost[sel & cl$runs$seeds == 5L]) / mean(cl$runs$cost[sel & cl$runs$seeds == 1L])
}
t_cl <- t_cl[order(t_cl$w, t_cl$seeds, t_cl$mu, t_cl$method), ]
write.csv(t_cl, file.path(data_dir, "T8_cluster_v2.csv"), row.names = FALSE)

# Console summary -------------------------------------------------------------------
cat("\n== EXTRA: dropped cells ==\n"); print(ex$dropped)
cat("\n== EXTRA: validity, max per method ==\n")
print(aggregate(cbind(false_rate, miscover_rate) ~ method, t_ex, max))
cat("\n== EXTRA: worst false-declaration cells (rate > alpha/2) ==\n")
print(t_ex[t_ex$false_rate > 0.025, c("method", "mu", "h", "truth", "false_rate", "false_lo", "false_hi")], row.names = FALSE)
cat("\n== EXTRA: median n (with p_inconclusive) ==\n")
for (hh in sort(unique(t_ex$h))) {
  x <- t_ex[t_ex$h == hh, ]
  w <- reshape(data.frame(method = x$method, mu = x$mu, v = sprintf("%.0f (%.0f%%)", x$median_n, 100 * x$p_inconclusive)),
               idvar = "method", timevar = "mu", direction = "wide")
  cat("h =", hh, "\n"); print(w, row.names = FALSE)
}
cat("\n== EXTRA: savi / betting ratio ==\n")
print(t_ex[t_ex$method == "savi_cs", c("mu", "h", "ratio_to_betting", "ratio_lo", "ratio_hi")], row.names = FALSE)
cat("\n== CLUSTER v2 ==\n")
print(t_cl[, c("method", "mu", "seeds", "w", "false_rate", "false_hi", "miscover_rate", "p_inconclusive", "median_n", "mean_cost", "ratio_n_vs_1seed", "ratio_cost_vs_1seed")], row.names = FALSE, digits = 3)
