# experiments/analysis/pilot_summary.R -- tables for the week-8 pilot and replication
# sizing for the main study (master plan §14.3).
# Usage: Rscript experiments/analysis/pilot_summary.R [out_dir]
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "experiments/out/pilot_2026-09-25"
runs <- read.csv(file.path(out_dir, "runs.csv"))
cfg <- readRDS(file.path(out_dir, "config.rds"))
if (requireNamespace("devtools", quietly = TRUE)) devtools::load_all(quiet = TRUE) else library(seqbench)

agg <- aggregate(cbind(false_declaration, inconclusive, miscover, n_instances) ~ boundary + mu + sd_d + truth,
                 data = runs, FUN = mean)
agg$R <- cfg$R
agg$se_false <- sqrt(agg$false_declaration * (1 - agg$false_declaration) / agg$R)
med <- aggregate(n_instances ~ boundary + mu + sd_d, data = runs, FUN = median)
names(med)[4] <- "median_n"
agg <- merge(agg, med, by = c("boundary", "mu", "sd_d"))
agg <- agg[order(agg$boundary, agg$sd_d, agg$mu), ]

cat("\n== Validity: false declarations and uniform miscoverage (alpha =", cfg$alpha, ") ==\n")
print(agg[, c("boundary", "sd_d", "mu", "truth", "false_declaration", "se_false", "miscover")], digits = 3, row.names = FALSE)

cat("\n== Cost: P(inconclusive), mean and median instances (n_max =", cfg$n_max, ") ==\n")
print(agg[, c("boundary", "sd_d", "mu", "inconclusive", "n_instances", "median_n")], digits = 3, row.names = FALSE)

# Sample-complexity check (theory.md Prop. 2): for hoeffding, tau <= t*(Delta) on the
# coverage event; for betting/EB, tau should scale like sd^2 / Delta^2.
design <- comparison_design(alpha = cfg$alpha, margin = cfg$margin, bounds = cfg$bounds)
sc <- agg[agg$boundary != "naive_fixed" & abs(agg$mu) != cfg$margin, ]
sc$distance <- abs(abs(sc$mu) - cfg$margin)
sc$t_star_hoeffding <- vapply(sc$distance, function(d) planning_horizon(design, d), numeric(1))
sc$t_eb_approx <- mapply(function(d, s) planning_horizon(design, d, sd = s), sc$distance, sc$sd_d)
sc$ratio_to_eb <- sc$median_n / sc$t_eb_approx
cat("\n== Stopping time vs. horizons: guaranteed Hoeffding t*(Delta) and variance-based EB approximation ==\n")
cat("   (EB median / t_eb_approx should be near 1; betting below it; hoeffding <= t* unless inconclusive)\n")
print(sc[, c("boundary", "sd_d", "mu", "distance", "median_n", "inconclusive", "t_star_hoeffding", "t_eb_approx", "ratio_to_eb")], digits = 3, row.names = FALSE)

fit <- subset(sc, boundary %in% c("betting", "empirical_bernstein") & inconclusive < 0.5 & distance > 0)
if (nrow(fit) >= 4) {
  m <- lm(log(median_n) ~ log(distance) + log(sd_d), data = fit)
  cat("\n== log(median tau) ~ log(Delta) + log(sd): slopes (theory: about -2 and +2) ==\n")
  print(coef(summary(m))[, 1:2], digits = 3)
}

# Replication sizing for the main study: worst-case p for error-rate cells is alpha;
# target half-width 0.01 (95%) on p ~ 0.05 -> R ~ 1.96^2 p(1-p)/0.01^2.
p <- max(cfg$alpha, max(agg$false_declaration[agg$boundary != "naive_fixed"]))
R_main <- ceiling(1.96^2 * p * (1 - p) / 0.01^2)
cat(sprintf("\n== Replication sizing: p_max = %.3f -> R = %d per cell for +-0.01 (95%%); at +-0.005: %d ==\n",
            p, R_main, ceiling(1.96^2 * p * (1 - p) / 0.005^2)))
sec <- aggregate(seconds ~ boundary, data = runs[!duplicated(runs$cell), ], FUN = mean)
cat("mean seconds per cell of", cfg$R, "reps:\n"); print(sec, row.names = FALSE)
write.csv(agg, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(sc, file.path(out_dir, "sample_complexity.csv"), row.names = FALSE)
