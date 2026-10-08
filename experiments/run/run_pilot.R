# experiments/run/run_pilot.R -- week-8 pilot: validity (uniform coverage, false
# declarations), cost (stopping time, P(inconclusive)) and the sample-complexity
# check of paper/theory.md Prop. 2 / betting conjecture.
# Usage: Rscript experiments/run/run_pilot.R [config_path]
# Output: experiments/out/<name>/{runs.csv, config.rds, session.txt}

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args)) args[1] else "experiments/configs/pilot_2026-09-25.R"
cfg <- source(cfg_path)$value
if (requireNamespace("devtools", quietly = TRUE) && !("seqbench" %in% loadedNamespaces())) {
  devtools::load_all(quiet = TRUE)
} else library(seqbench)

out_dir <- file.path("experiments/out", cfg$name)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(cfg, file.path(out_dir, "config.rds"))
writeLines(c(capture.output(sessionInfo()),
             paste("seqbench", as.character(packageVersion("seqbench"))),
             paste("git", system("git rev-parse --short HEAD", intern = TRUE))),
           file.path(out_dir, "session.txt"))

gen_losses <- function(n, mu, h, seeds) {
  inst <- rep(seq_len(n), each = seeds)
  la <- runif(length(inst), 0.4, 0.8)
  lb <- la - mu + runif(length(inst), -h, h)      # E[loss_a - loss_b] = mu
  data.frame(instance = inst, seed = rep(seq_len(seeds), times = n), loss_a = la, loss_b = lb)
}
truth_of <- function(mu, margin) if (abs(mu) <= margin) "equivalent" else if (mu > 0) "B" else "A"

cells <- expand.grid(boundary = cfg$boundaries, mu = cfg$mu, h = cfg$half_width,
                     stringsAsFactors = FALSE)
cells$cell <- seq_len(nrow(cells))
n_cores <- max(1L, parallel::detectCores() - 2L)
cat(sprintf("%d cells x %d reps, %d cores\n", nrow(cells), cfg$R, n_cores))

run_cell <- function(i) {
  cl <- cells[i, ]
  design <- comparison_design(alpha = cfg$alpha, margin = cfg$margin, bounds = cfg$bounds,
                              boundary = cl$boundary, n_max = cfg$n_max)
  truth <- truth_of(cl$mu, cfg$margin)
  res <- vector("list", cfg$R)
  for (r in seq_len(cfg$R)) {
    set.seed(cfg$seed + 1000L * i + r)       # reproducible per (cell, rep)
    l <- gen_losses(cfg$n_max, cl$mu, cl$h, cfg$seeds_per_instance)
    st <- suppressWarnings(update_comparison(initialize_comparison(design), l))
    tr <- st$trajectory
    dec <- as.vector(stopping_decision(st))
    res[[r]] <- data.frame(
      cell = i, rep = r, boundary = cl$boundary, mu = cl$mu, sd_d = cl$h / sqrt(3),
      truth = truth, decision = dec, reason = attr(stopping_decision(st), "reason"),
      false_declaration = dec %in% c("A", "B", "equivalent") && dec != truth,
      inconclusive = dec == "inconclusive",
      n_instances = nrow(tr), cost = st$cost,
      miscover = any(tr$lower_running > cl$mu | tr$upper_running < cl$mu, na.rm = TRUE),
      final_lower = tr$lower_running[nrow(tr)], final_upper = tr$upper_running[nrow(tr)],
      seconds = NA_real_)
  }
  do.call(rbind, res)
}

t0 <- Sys.time()
runs <- parallel::mclapply(cells$cell, function(i) {
  t1 <- Sys.time(); r <- run_cell(i); r$seconds <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("[%s] cell %2d/%d %-20s mu=%.2f h=%.3f  %.0fs\n", format(Sys.time(), "%H:%M:%S"), i, nrow(cells),
              cells$boundary[i], cells$mu[i], cells$h[i], r$seconds[1]),
      file = file.path(out_dir, "progress.log"), append = TRUE)
  r
}, mc.cores = n_cores, mc.set.seed = FALSE)
runs <- do.call(rbind, runs)
runs$seqbench_version <- as.character(packageVersion("seqbench"))
write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)
cat(sprintf("done in %.1f min; %d rows written to %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(runs), out_dir))
