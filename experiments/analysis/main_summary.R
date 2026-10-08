# experiments/analysis/main_summary.R -- pre-registered analysis of the main simulation
# study (config main_sim_2026-09-25.R). Written and frozen on 2026-09-25 BEFORE the
# results were available (master plan §14.2 item 1); tested on the 3-replication smoke run.
# Usage: Rscript experiments/analysis/main_summary.R [out_dir]
# Produces, in out_dir: tables T1-T4 (csv), figures F1-F4 (png), summary.md.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "experiments/out/main_sim_2026-09-25"
runs <- read.csv(file.path(out_dir, "runs.csv"), stringsAsFactors = FALSE)
cfg <- readRDS(file.path(out_dir, "config.rds"))
suppressPackageStartupMessages(library(ggplot2))
stopifnot(all(runs$decision != "error"))
# Completeness: every (dgp, method) cell with exactly R replications; list what is missing.
expected <- expand.grid(dgp = seq_len(nrow(cfg$dgp)), method = cfg$methods, stringsAsFactors = FALSE)
got <- aggregate(rep ~ dgp + method, data = runs, FUN = length)
chk <- merge(expected, got, all.x = TRUE); chk$rep[is.na(chk$rep)] <- 0L
missing_cells <- chk[chk$rep != cfg$R, ]
complete <- nrow(missing_cells) == 0L
if (!complete) message("PARTIAL RESULTS: ", nrow(missing_cells), " of ", nrow(expected), " cells missing or incomplete")
stopifnot(!anyDuplicated(runs[c("dgp", "rep", "method")]))

alpha <- cfg$alpha; margin <- cfg$margin; rng <- diff(cfg$bounds) * 2   # paired-difference range
R <- cfg$R
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n; c <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind(lower = (c - h) / d, upper = (c + h) / d)
}
runs$valid_method <- !runs$method %in% c("seqbench:naive_fixed")
runs$scenario <- with(runs, ifelse(abs(mu) < 1e-12, "null", ifelse(abs(mu) < margin, "equivalent",
                                   ifelse(abs(abs(mu) - margin) < 1e-12, "at margin", "effect"))))
runs$distance <- abs(abs(runs$mu) - margin)
sd_of <- function(h, h_seed, seeds, paired) sqrt((h^2 + 2 * h_seed^2) / (3 * seeds) + ifelse(paired, 0, 2 * 0.25^2 / 12))
runs$sd_d <- with(runs, sd_of(h, h_seed, seeds, paired))
key <- c("method", "ablation", "mu", "h", "sd_d", "scenario", "truth", "distance")

# T1 validity -----------------------------------------------------------------
t1 <- aggregate(cbind(false_declaration, miscover) ~ method + ablation + scenario + mu + h,
                data = runs, FUN = function(z) sum(z))
n1 <- aggregate(rep ~ method + ablation + scenario + mu + h, data = runs, FUN = length)
t1 <- merge(t1, n1); names(t1)[names(t1) == "rep"] <- "R"
t1$false_rate <- t1$false_declaration / t1$R
t1[c("false_lo", "false_hi")] <- wilson(t1$false_declaration, t1$R)
t1$miscover_rate <- t1$miscover / t1$R
t1[c("miscover_lo", "miscover_hi")] <- wilson(t1$miscover, t1$R)
t1$exceeds_alpha <- t1$false_lo > alpha          # flag: 95% Wilson lower bound above alpha
t1 <- t1[order(t1$ablation, t1$method, t1$h, t1$mu), ]
write.csv(t1, file.path(out_dir, "T1_validity.csv"), row.names = FALSE)

# T2 cost -----------------------------------------------------------------------
# Instances consumed include runs stopped by the budget; p_inconclusive says how many.
# Monte Carlo uncertainty: Wilson for p_inconclusive; percentile bootstrap over
# replications (B = 1000) for the median; PAIRED bootstrap (same replication streams)
# for the ratio of medians to the betting boundary.
set.seed(1)
B <- 1000L
boot_idx <- replicate(B, sample.int(R, R, replace = TRUE), simplify = FALSE)
t2 <- do.call(rbind, lapply(split(runs, runs[c("method", "ablation", "mu", "h")], drop = TRUE), function(d) {
  d <- d[order(d$rep), ]
  med_b <- vapply(boot_idx, function(i) median(d$n_instances[i[i <= nrow(d)]]), numeric(1))
  data.frame(method = d$method[1], ablation = d$ablation[1], mu = d$mu[1], h = d$h[1], sd_d = d$sd_d[1],
             scenario = d$scenario[1], p_inconclusive = mean(d$inconclusive),
             p_inc_lo = wilson(sum(d$inconclusive), nrow(d))[1], p_inc_hi = wilson(sum(d$inconclusive), nrow(d))[2],
             mean_n = mean(d$n_instances), median_n = median(d$n_instances),
             median_lo = unname(quantile(med_b, 0.025)), median_hi = unname(quantile(med_b, 0.975)),
             q90_n = unname(quantile(d$n_instances, 0.9)), mean_cost = mean(d$cost))
}))
ref <- t2[t2$method == "seqbench:betting", c("ablation", "mu", "h", "median_n")]
names(ref)[4] <- "median_n_betting"
t2 <- merge(t2, ref, all.x = TRUE); t2$ratio_to_betting <- t2$median_n / t2$median_n_betting   # NA when betting cell missing
paired_ratio_ci <- function(d_m, d_b) {   # same reps in both
  d_m <- d_m[order(d_m$rep), ]; d_b <- d_b[order(d_b$rep), ]
  stopifnot(identical(d_m$rep, d_b$rep))
  r <- vapply(boot_idx, function(i) { i <- i[i <= nrow(d_m)]; median(d_m$n_instances[i]) / median(d_b$n_instances[i]) }, numeric(1))
  unname(quantile(r, c(0.025, 0.975)))
}
t2$ratio_lo <- NA_real_; t2$ratio_hi <- NA_real_
for (i in seq_len(nrow(t2))) {
  if (t2$method[i] == "seqbench:betting" || is.na(t2$ratio_to_betting[i])) next
  sel <- runs$ablation == t2$ablation[i] & runs$mu == t2$mu[i] & runs$h == t2$h[i]
  ci <- paired_ratio_ci(runs[sel & runs$method == t2$method[i], ], runs[sel & runs$method == "seqbench:betting", ])
  t2$ratio_lo[i] <- ci[1]; t2$ratio_hi[i] <- ci[2]
}
t2 <- t2[order(t2$ablation, t2$h, t2$mu, t2$method), ]
write.csv(t2, file.path(out_dir, "T2_cost.csv"), row.names = FALSE)

# T3 ablations: pairing and seed clusters, betting only -----------------------
b <- t2[t2$method == "seqbench:betting" & t2$h == 0.173, ]
base <- b[b$ablation == "none", c("mu", "median_n", "p_inconclusive", "mean_cost")]
names(base)[2:4] <- c("median_n_paired", "p_inc_paired", "mean_cost_paired")
t3 <- merge(b[b$ablation != "none", c("ablation", "mu", "sd_d", "median_n", "p_inconclusive", "mean_cost")], base)
t3$ratio_n <- t3$median_n / t3$median_n_paired
t3$ratio_cost <- t3$mean_cost / t3$mean_cost_paired      # evaluations, not instances
write.csv(t3, file.path(out_dir, "T3_ablations.csv"), row.names = FALSE)

# T4 regimes: median stopping time vs prediction of the truncated-bet regime -----
# effective bet cap near m = 1/2: c/m = 1 for betting (c = 1/2), c = 1/2 for empirical Bernstein
cap <- c("seqbench:betting" = 1, "seqbench:empirical_bernstein" = 0.5)
t4 <- t2[t2$method %in% names(cap) & t2$scenario != "at margin" & t2$ablation == "none", ]
t4$distance <- abs(abs(t4$mu) - margin)
t4$pred_truncated <- log(2 / alpha) / (cap[t4$method] * t4$distance / rng)
t4$ratio_obs_pred <- t4$median_n / t4$pred_truncated
t4$crossover_t <- with(t4, {  # t at which lambda-tilde falls below the cap: t log(t+1) = 2 log(2/alpha) / (sigma2 cap^2)
  s2 <- (sd_d / rng)^2; target <- 2 * log(2 / alpha) / (s2 * cap[method]^2)
  vapply(target, function(tg) { t <- 1; while (t * log(t + 1) < tg) t <- t * 1.05 + 1; round(t) }, numeric(1)) })
t4$regime <- ifelse(t4$median_n < t4$crossover_t, "truncated", "variance")
fit <- lm(log(median_n) ~ log(distance) + log(sd_d) + method, data = t4[t4$p_inconclusive < 0.5, ])
write.csv(t4, file.path(out_dir, "T4_regimes.csv"), row.names = FALSE)

# Figures -----------------------------------------------------------------------
theme_set(theme_minimal(base_size = 11))
f1 <- ggplot(t1, aes(x = factor(mu), y = false_rate, colour = method, group = method)) +
  geom_hline(yintercept = alpha, linetype = 2) + geom_line() + geom_point() +
  geom_errorbar(aes(ymin = false_lo, ymax = false_hi), width = 0.15) +
  facet_grid(ablation ~ h, labeller = label_both) +
  labs(x = "true mean paired difference mu", y = "false-declaration rate (95% Wilson)", title = "F1 Validity")
ggsave(file.path(out_dir, "F1_validity.png"), f1, width = 10, height = 7, dpi = 150)

f2 <- ggplot(t4, aes(x = distance, y = median_n, colour = method, shape = factor(round(sd_d, 3)))) +
  geom_line(aes(y = pred_truncated, linetype = "truncated-bet prediction")) + geom_point(size = 2.5) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "distance |mu| - margin (loss units)", y = "median instances to decision",
       shape = "sd(D)", linetype = NULL, title = "F2 Stopping time vs distance (log-log)")
ggsave(file.path(out_dir, "F2_stopping_time.png"), f2, width = 8, height = 5, dpi = 150)

f3 <- ggplot(t2[t2$ablation == "none", ], aes(x = factor(mu), y = p_inconclusive, colour = method, group = method)) +
  geom_line() + geom_point() + facet_wrap(~ h, labeller = label_both) +
  labs(x = "mu", y = "P(inconclusive at the method's own sample limit)", title = "F3 Inconclusive outcomes")
ggsave(file.path(out_dir, "F3_inconclusive.png"), f3, width = 9, height = 4.5, dpi = 150)

d4 <- runs[runs$ablation == "none" & runs$mu == 0.05, ]
f4 <- ggplot(d4, aes(x = method, y = n_instances, fill = method)) + geom_boxplot(show.legend = FALSE) +
  facet_wrap(~ h, labeller = label_both) + coord_flip() +
  labs(x = NULL, y = "instances to decision", title = "F4 Stopping-time distribution, mu = 0.05")
ggsave(file.path(out_dir, "F4_stopping_dist.png"), f4, width = 9, height = 4.5, dpi = 150)

# summary.md --------------------------------------------------------------------
md <- c(
  sprintf("# %s — summary (R = %d per cell)%s", cfg$name, R, if (complete) "" else " — PARTIAL"), "",
  if (!complete) c(sprintf("Missing/incomplete cells (%d): %s", nrow(missing_cells),
                           paste(sprintf("dgp%03d:%s(%d)", missing_cells$dgp, missing_cells$method, missing_cells$rep), collapse = ", ")), "") else NULL,
  "## Validity (max false-declaration rate per method over all cells)",
  {
    m <- aggregate(false_rate ~ method, data = t1, FUN = max); mc <- aggregate(miscover_rate ~ method, data = t1, FUN = max)
    x <- merge(m, mc); c("| method | max false rate | max miscoverage |", "|---|---|---|",
      sprintf("| %s | %.4f | %.4f |", x$method, x$false_rate, x$miscover_rate)) },
  "", sprintf("Cells whose 95%% lower bound exceeds alpha = %s: %d (%s)", alpha, sum(t1$exceeds_alpha),
              paste(unique(t1$method[t1$exceeds_alpha]), collapse = ", ")),
  "", "## Instances consumed (median; includes runs stopped by the budget, see p_inconclusive), baseline cells, sd(D) = 0.10",
  { x <- t2[t2$ablation == "none" & t2$h == 0.173, ]; w <- reshape(x[, c("method", "mu", "median_n")], idvar = "method", timevar = "mu", direction = "wide")
    c(paste0("| method | ", paste(sub("median_n.", "mu=", names(w)[-1]), collapse = " | "), " |"),
      paste0("|---|", paste(rep("---", ncol(w) - 1), collapse = "|"), "|"),
      apply(w, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) },
  "", "## Ablations (betting, sd(D) = 0.10): relative to paired single-seed (instances AND evaluations)",
  sprintf("- %s, mu = %.2f: instances ratio %.2f, evaluations (cost) ratio %.2f (P(inconclusive) %.3f vs %.3f)", t3$ablation, t3$mu, t3$ratio_n, t3$ratio_cost, t3$p_inconclusive, t3$p_inc_paired),
  "", "## Regimes", sprintf("log-log slopes: distance %.2f (se %.2f), sd %.2f (se %.2f)",
                             coef(fit)[["log(distance)"]], coef(summary(fit))["log(distance)", 2],
                             coef(fit)[["log(sd_d)"]], coef(summary(fit))["log(sd_d)", 2]),
  sprintf("- %s, mu = %.2f, sd = %.3f: median %.0f, truncated prediction %.0f (ratio %.2f), crossover t = %.0f, regime %s",
          t4$method, t4$mu, t4$sd_d, t4$median_n, t4$pred_truncated, t4$ratio_obs_pred, t4$crossover_t, t4$regime))
writeLines(md, file.path(out_dir, "summary.md"))
cat(md, sep = "\n")
