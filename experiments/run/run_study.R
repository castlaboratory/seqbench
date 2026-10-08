# experiments/run/run_study.R -- generic runner for a frozen config.
# Usage: Rscript experiments/run/run_study.R experiments/configs/<name>.R [cores] [R_override]
# Env SEQBENCH_METHODS=<regex> restricts the methods run on this machine (split runs).
# R_override (smoke tests only) replaces cfg$R and writes to experiments/out/<name>_smoke_R<k>/.
# Resumable: each (dgp cell, method) is saved to experiments/out/<name>/cells/ and
# skipped on rerun. Data for (dgp cell, rep) is seeded identically for every method,
# so methods are compared on the same paired-loss streams. Dynamic scheduling.

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- args[1]; cores <- if (length(args) > 1) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
cfg <- source(cfg_path)$value
if (length(args) > 2) { cfg$R <- as.integer(args[3]); cfg$name <- sprintf("%s_smoke_R%d", cfg$name, cfg$R) }
if (dir.exists("R") && requireNamespace("devtools", quietly = TRUE)) devtools::load_all(quiet = TRUE) else library(seqbench)
source("experiments/run/methods.R")
if (any(grepl("^savi_cs", cfg$methods))) loadNamespace("safestats")   # so that session.txt records it
# Provenance: hash of the config file and of the experiment scripts; every cell records it
# and a resumed run refuses cells produced under a different hash.
prov <- list(config_md5 = unname(tools::md5sum(cfg_path)),
             scripts_md5 = unname(tools::md5sum(c("experiments/run/run_study.R", "experiments/run/methods.R"))),
             seqbench = as.character(packageVersion("seqbench")), R = R.version.string)

out_dir <- file.path("experiments/out", cfg$name)
dir.create(file.path(out_dir, "cells"), showWarnings = FALSE, recursive = TRUE)
saveRDS(cfg, file.path(out_dir, "config.rds"))
file.copy(cfg_path, file.path(out_dir, basename(cfg_path)), overwrite = TRUE)
writeLines(c(capture.output(sessionInfo()), paste("seqbench", prov$seqbench),
             paste("config_md5", prov$config_md5), paste("scripts_md5", paste(prov$scripts_md5, collapse = " ")),
             paste("git", tryCatch(suppressWarnings(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE)), error = function(e) "n/a")),
             paste("host", Sys.info()[["nodename"]]), paste("started", format(Sys.time()))),
           file.path(out_dir, sprintf("session_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S"))))
if (!file.exists(file.path(out_dir, "session.txt"))) file.copy(file.path(out_dir, list.files(out_dir, "^session_")[1]), file.path(out_dir, "session.txt"))

# Truth per dgp cell, computed ONCE in the parent (mu_of may be a Monte Carlo reference)
mu_table <- vapply(seq_len(nrow(cfg$dgp)), function(i) cfg$mu_of(cfg$dgp[i, , drop = FALSE]), numeric(1))
truth_table <- vapply(seq_len(nrow(cfg$dgp)), function(i) cfg$truth_of(cfg$dgp[i, , drop = FALSE], cfg$margin), character(1))
jobs <- expand.grid(dgp = seq_len(nrow(cfg$dgp)), method = cfg$methods, stringsAsFactors = FALSE)
jobs$file <- file.path(out_dir, "cells", sprintf("dgp%03d_%s.rds", jobs$dgp, gsub("[^a-z0-9]", "_", jobs$method)))
# Resume: a cell is reused only if complete and produced under the same config hash.
cell_ok <- vapply(jobs$file, function(f) {
  if (!file.exists(f)) return(FALSE)
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r)) return(FALSE)
  pr <- attr(r, "provenance")
  ok <- nrow(r) == cfg$R && !any(r$decision == "error") &&
    (is.null(pr) || identical(pr$config_md5, prov$config_md5))
  if (!is.null(pr) && !identical(pr$config_md5, prov$config_md5))
    stop("cell ", basename(f), " was produced under a different config (md5 ", pr$config_md5, "); refusing to resume")
  ok
}, logical(1))
todo <- which(!cell_ok)
# Optional split across machines (added 2026-09-28): SEQBENCH_METHODS is a regex; only matching
# methods run here. Cells run elsewhere are copied into cells/ and a final rerun assembles runs.csv.
if (nzchar(Sys.getenv("SEQBENCH_METHODS"))) todo <- todo[grepl(Sys.getenv("SEQBENCH_METHODS"), jobs$method[todo])]
logmsg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(...), "\n", sep = "",
                         file = file.path(out_dir, "progress.log"), append = TRUE)
logmsg("%d jobs (%d dgp cells x %d methods), %d to do, %d cores", nrow(jobs), nrow(cfg$dgp), length(cfg$methods), length(todo), cores)

cached_gen <- function(dgp, r, p) {
  f <- file.path(out_dir, "data", sprintf("dgp%03d_rep%05d.rds", dgp, r))
  if (file.exists(f)) return(readRDS(f))
  dir.create(dirname(f), showWarnings = FALSE)
  losses <- cfg$gen(p, cfg$n_max); tmp <- paste0(f, ".", Sys.getpid()); saveRDS(losses, tmp); file.rename(tmp, f)
  losses
}

run_job <- function(j) {
  t0 <- Sys.time()
  p <- cfg$dgp[jobs$dgp[j], , drop = FALSE]
  mu_true <- mu_table[jobs$dgp[j]]; truth <- truth_table[jobs$dgp[j]]
  rows <- vector("list", cfg$R)
  for (r in seq_len(cfg$R)) {
    set.seed(cfg$seed + 100000L * jobs$dgp[j] + r)        # same data for every method
    losses <- if (isTRUE(cfg$cache_data)) cached_gen(jobs$dgp[j], r, p) else cfg$gen(p, cfg$n_max)
    row <- tryCatch(run_method(jobs$method[j], losses, cfg, mu_true, truth, p),
                    error = function(e) result_row(jobs$method[j], "error", conditionMessage(e), NA, NA, NA, NA, NA, truth))
    rows[[r]] <- cbind(dgp = jobs$dgp[j], rep = r, p, mu_true = mu_true, truth = truth, row, row.names = NULL)
  }
  res <- do.call(rbind, rows)
  res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  attr(res, "provenance") <- c(prov, list(finished = format(Sys.time())))
  tmp <- paste0(jobs$file[j], ".", Sys.getpid(), ".tmp"); saveRDS(res, tmp)
  if (!file.rename(tmp, jobs$file[j])) stop("could not write ", jobs$file[j])
  logmsg("job %3d/%d dgp %3d %-28s %6.0fs  err=%d", j, nrow(jobs), jobs$dgp[j], jobs$method[j], res$seconds[1], sum(res$decision == "error"))
  invisible(NULL)
}
t0 <- Sys.time()
if (isTRUE(cfg$cache_data)) {   # generate every (dgp, rep) stream once, in parallel, before the methods run
  grid <- expand.grid(rep = seq_len(cfg$R), dgp = unique(jobs$dgp[todo]))
  invisible(parallel::mclapply(seq_len(nrow(grid)), function(i) {
    set.seed(cfg$seed + 100000L * grid$dgp[i] + grid$rep[i])
    cached_gen(grid$dgp[i], grid$rep[i], cfg$dgp[grid$dgp[i], , drop = FALSE]); NULL
  }, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE))
  logmsg("data cache ready: %d streams in %.1f min", nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
out <- parallel::mclapply(todo, run_job, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
bad <- vapply(out, function(o) inherits(o, "try-error") || (!is.null(o) && !is.null(attr(o, "condition"))), logical(1))
if (any(bad)) logmsg("WARNING: %d worker(s) failed: %s", sum(bad), paste(head(vapply(out[bad], as.character, character(1)), 3), collapse = " | "))
done <- file.exists(jobs$file)
if (!all(done)) logmsg("WARNING: %d cell(s) missing after run: %s", sum(!done), paste(head(basename(jobs$file[!done]), 5), collapse = ", "))
rbind_fill <- function(l) { cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] })) }
runs <- rbind_fill(lapply(jobs$file[done], function(f) { r <- readRDS(f); pr <- attr(r, "provenance")
  r$cell_seqbench_version <- if (is.null(pr)) NA_character_ else pr$seqbench; r }))
runs$seqbench_version <- as.character(packageVersion("seqbench"))   # version that assembled runs.csv
write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)
logmsg("done in %.1f min; %d rows", as.numeric(difftime(Sys.time(), t0, units = "mins")), nrow(runs))
cat(readLines(file.path(out_dir, "progress.log")), sep = "\n")
