# dev/repro_betting.R -- week-3 reproduction (PLANO.md §5): Monte Carlo check that the
# WSR (2024) confidence sequences in dev/wsr_cs.R keep time-uniform coverage, width
# comparison between boundaries, and a cross-check against seqcomp (Choe & Ramdas 2024).
# Run: Rscript dev/repro_betting.R   (~ a few minutes)
# Writes experiments/out/repro_betting_<date>.csv. Everything seeded.

source("dev/wsr_cs.R")
suppressPackageStartupMessages(library(seqcomp))
set.seed(20260925)
alpha <- 0.05; n <- 500L; R <- 1000L
t_report <- c(30L, 100L, 500L)

scenarios <- list(
  beta_10_30   = list(mu = 10 / 40, gen = function(n) rbeta(n, 10, 30)),
  bernoulli_05 = list(mu = 0.5,     gen = function(n) rbinom(n, 1, 0.5)),
  two_point    = list(mu = 0.5,     gen = function(n) sample(c(0.4, 0.6), n, TRUE)),
  paired_diff  = list(mu = NA,      gen = NULL)   # filled below
)
# Paired-difference-like scenario: losses L_A, L_B in [0,1], D = L_A - L_B in [-1,1],
# X = (D + 1) / 2 in [0,1]. L_B = L_A shifted with noise, strong pairing.
scenarios$paired_diff$gen <- function(n) {
  la <- rbeta(n, 2, 5); lb <- pmin(pmax(la - 0.05 + rnorm(n, 0, 0.05), 0), 1)
  ((la - lb) + 1) / 2
}
scenarios$paired_diff$mu <- mean(scenarios$paired_diff$gen(2e6))   # MC value of mu

mc_se <- function(p, R) sqrt(p * (1 - p) / R)

rows <- list()
for (sc in names(scenarios)) {
  mu <- scenarios[[sc]]$mu; gen <- scenarios[[sc]]$gen
  miss <- c(hoeffding = 0, empirical_bernstein = 0, betting = 0, naive_fixed = 0,
            seqcomp_hoeffding = 0, seqcomp_bernstein = 0)
  widths <- matrix(0, length(miss), length(t_report), dimnames = list(names(miss), t_report))
  t0 <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    x <- gen(n)
    cs <- list(hoeffding = cs_prpl_h(x, alpha), empirical_bernstein = cs_prpl_eb(x, alpha),
               naive_fixed = cs_naive_t(x, alpha))
    # seqcomp boundaries (Howard et al. mixture) on the same data: scores1 = x, scores2 = 0
    sh <- seqcomp::cs_hoeffding(x, rep(0, n), alpha = alpha, c = 1)
    sb <- seqcomp::cs_bernstein(x, rep(0, n), alpha = alpha, c = 2)
    cs$seqcomp_hoeffding <- data.frame(t = sh$t, estimate = sh$estimate,
                                       lower = cummax(sh$lower), upper = cummin(sh$upper))
    cs$seqcomp_bernstein <- data.frame(t = sb$t, estimate = sb$estimate,
                                       lower = cummax(sb$lower), upper = cummin(sb$upper))
    for (b in names(cs)) {
      miss[b] <- miss[b] + any(cs[[b]]$lower > mu | cs[[b]]$upper < mu, na.rm = TRUE)
      widths[b, ] <- widths[b, ] + (cs[[b]]$upper - cs[[b]]$lower)[t_report]
    }
    # betting: exact miscoverage at mu (no grid); width via grid on a subset of reps
    miss["betting"] <- miss["betting"] + hedged_misses_mu(x, mu, alpha)
    if (r <= 200L) {
      hb <- cs_hedged(x, alpha, n_grid = 1001L)
      widths["betting", ] <- widths["betting", ] + (hb$upper - hb$lower)[t_report]
    }
  }
  el <- proc.time()[["elapsed"]] - t0
  widths <- widths / R; widths["betting", ] <- widths["betting", ] * R / 200
  for (b in names(miss)) {
    p <- miss[[b]] / R
    rows[[length(rows) + 1]] <- data.frame(
      scenario = sc, mu = mu, boundary = b, R = R, n = n, alpha = alpha,
      miscoverage = p, mc_se = mc_se(p, R),
      width_30 = widths[b, "30"], width_100 = widths[b, "100"], width_500 = widths[b, "500"],
      seconds = el)
  }
  cat(sprintf("[%s] mu = %.4f, %.0fs\n", sc, mu, el))
  print(do.call(rbind, rows[(length(rows) - length(miss) + 1):length(rows)])[, c("boundary", "miscoverage", "mc_se", "width_30", "width_100", "width_500")], digits = 3)
}
out <- do.call(rbind, rows)
out$seqbench_version <- as.character(packageVersion("seqbench"))
out$seqcomp_version <- as.character(packageVersion("seqcomp"))
out$r_version <- R.version.string
out$seed <- 20260925
f <- sprintf("experiments/out/repro_betting_%s.csv", format(Sys.Date()))
write.csv(out, f, row.names = FALSE)
cat("written", f, "\n")
