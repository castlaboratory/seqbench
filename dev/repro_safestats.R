# dev/repro_safestats.R -- week-3 reproduction (PLANO.md §5): safestats savi t-test
# (paired/one-sample) versus a naive t-test inspected at every step versus the WSR
# hedged betting CS, on the same bounded paired differences. Measures (i) type I error
# under repeated inspection, (ii) time-uniform miscoverage of savi's confidence sequence,
# (iii) safestats' "relevance test" stopping behaviour, (iv) power at an effect of ~0.5 sd.
# Run: Rscript dev/repro_safestats.R   Writes experiments/out/repro_safestats_<date>.csv.

source("dev/wsr_cs.R")
suppressPackageStartupMessages(library(safestats))
set.seed(20260925)
alpha <- 0.05; n_max <- 200L; R <- 500L

# Bounded paired differences D in [-1, 1]: D = mu + 0.5 * (2U - 1), sd(D) = 0.5/sqrt(3) ~ 0.289.
# X = (D + 1)/2 in [0, 1] for the betting CS. mu = 0.15 is ~0.52 sd, close to deltaMin = 0.5.
gen <- function(n, mu) mu + 0.5 * (2 * runif(n) - 1)
sd_d <- 0.5 / sqrt(3)

design <- designSaviT(deltaMin = 0.5, alpha = alpha, alternative = "twoSided",
                      testType = "oneSample", relevanceTest = TRUE,
                      nSim = 100L, seed = 1, pb = FALSE)

rows <- list()
for (mu in c(0, 0.15)) {
  cnt <- c(naive_reject_ever = 0, naive_reject_final = 0, savi_reject_ever = 0,
           savi_miscover = 0, savi_relevance_stop = 0, betting_reject_ever = 0,
           betting_miscover = 0)
  stop_t <- c(savi = 0, betting = 0); n_stop <- c(savi = 0, betting = 0)
  t0 <- proc.time()[["elapsed"]]
  for (r in seq_len(R)) {
    d <- gen(n_max, mu); x <- (d + 1) / 2; mu_x <- (mu + 1) / 2
    # (a) naive t-test inspected at every t >= 2
    nv <- cs_naive_t(d, alpha)
    rej <- with(nv, lower > 0 | upper < 0)
    cnt["naive_reject_ever"]  <- cnt["naive_reject_ever"]  + any(rej, na.rm = TRUE)
    cnt["naive_reject_final"] <- cnt["naive_reject_final"] + rej[n_max]
    # (b) safestats savi t-test, sequential
    sv <- saviTTest(d, designObj = design, sequential = TRUE)
    ev <- sv$eValueVec
    cnt["savi_reject_ever"] <- cnt["savi_reject_ever"] + any(ev >= 1 / alpha)
    if (any(ev >= 1 / alpha)) { stop_t["savi"] <- stop_t["savi"] + which(ev >= 1 / alpha)[1]; n_stop["savi"] <- n_stop["savi"] + 1 }
    cm <- sv$confSeqMatrix
    cnt["savi_miscover"] <- cnt["savi_miscover"] + any(cm[, 1] > mu | cm[, 2] < mu, na.rm = TRUE)
    cnt["savi_relevance_stop"] <- cnt["savi_relevance_stop"] + is.finite(sv$fptRelevance)
    # (c) WSR hedged betting CS on X: reject H0: mu = 0 when 1/2 leaves the CS
    cnt["betting_miscover"] <- cnt["betting_miscover"] + hedged_misses_mu(x, mu_x, alpha)
    rej_b <- hedged_log_capital(x, 0.5, alpha) >= log(1 / alpha)
    cnt["betting_reject_ever"] <- cnt["betting_reject_ever"] + any(rej_b)
    if (any(rej_b)) { stop_t["betting"] <- stop_t["betting"] + which(rej_b)[1]; n_stop["betting"] <- n_stop["betting"] + 1 }
  }
  el <- proc.time()[["elapsed"]] - t0
  p <- cnt / R
  rows[[length(rows) + 1]] <- data.frame(
    mu = mu, effect_sd = mu / sd_d, R = R, n_max = n_max, alpha = alpha,
    t(p), t(sqrt(p * (1 - p) / R)) |> `colnames<-`(paste0("se_", names(p))),
    mean_stop_savi = stop_t["savi"] / max(n_stop["savi"], 1),
    mean_stop_betting = stop_t["betting"] / max(n_stop["betting"], 1),
    seconds = el)
  cat(sprintf("mu = %.2f (%.2f sd), %.0fs\n", mu, mu / sd_d, el)); print(round(p, 3))
}
out <- do.call(rbind, rows)
out$safestats_version <- as.character(packageVersion("safestats"))
out$r_version <- R.version.string; out$seed <- 20260925
f <- sprintf("experiments/out/repro_safestats_%s.csv", format(Sys.Date()))
write.csv(out, f, row.names = FALSE); cat("written", f, "\n")
