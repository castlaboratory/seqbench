# experiments/run/methods.R -- the procedures compared in the experiments. Every
# method takes the same paired-loss data frame and returns one result row with the
# same columns, so that tables in experiments/analysis/ do not care which one ran.
#
#   seqbench:<boundary>[:opts]  the package (betting / empirical_bernstein / hoeffding /
#                         bernstein_declared / naive_fixed); opts: "c0.9" sets betting_c,
#                         "sdK" declares sd_max = K x the config's savi_sd_of(p) (bernstein_declared)
#   fixed:<n>             fixed-sample t interval at n instances, same three-way margin rule
#                         ("fixed:final" uses n_max); no adaptation
#   pseudorep:<boundary>  INVALID analysis for the cluster ablation: every (instance, seed) row
#                         is fed as if it were its own instance (pseudoreplication)
#   audit:<boundary>      no-stopping coverage audit: runs the boundary over the whole stream
#                         (n_max instances) and records whether the true mean ever leaves the
#                         running intersection (exact capital at the true mean for betting)
#   savi_cs               safestats anytime-valid t-test confidence sequence (eType = "mom",
#                         the package default; design fixed BEFORE the stream with
#                         deltaMin = margin / sd_ref, where sd_ref is the known DGP standard
#                         deviation of the per-instance difference when the config provides
#                         `savi_sd_of(p)`, or the sd of an independent pilot stream otherwise),
#                         running intersection + the same margin rule: "an existing sequential
#                         method correctly applied". savi_cs:pilot forces the pilot-stream sd even
#                         when the config provides savi_sd_of (added 2026-09-28 for the robustness grid)
#
# Result columns: method, decision, reason, n_instances, cost, miscover (true mean outside
# the reported running intersection at some time up to stopping; an empty set counts as
# miscoverage), cs_empty, final_lower, final_upper, false_declaration, inconclusive,
# extra (method-specific, character)

three_way <- function(L, U, margin) {
  if (anyNA(c(L, U)) || U < L) return("inconclusive")
  if (U < -margin) "A" else if (L > margin) "B" else if (L >= -margin && U <= margin) "equivalent" else "continue"
}

result_row <- function(method, decision, reason, n, cost, miscover, lo, hi, truth, extra = NA_character_,
                       cs_empty = FALSE) {
  data.frame(method = method, decision = decision, reason = reason, n_instances = n, cost = cost,
             miscover = miscover, cs_empty = cs_empty, final_lower = lo, final_upper = hi,
             false_declaration = if (is.na(truth)) NA else decision %in% c("A", "B", "equivalent") && decision != truth,
             inconclusive = decision == "inconclusive", extra = extra, stringsAsFactors = FALSE)
}

# Per-instance paired differences (seeds averaged), in the order given
instance_diffs <- function(losses) {
  d <- tapply(losses$loss_a - losses$loss_b, factor(losses$instance, levels = unique(losses$instance)), mean)
  as.numeric(d)
}

run_method <- function(method, losses, cfg, mu_true, truth, p = NULL) {
  parts <- strsplit(method, ":", fixed = TRUE)[[1]]
  switch(parts[1],
    seqbench = run_seqbench(parts[2], losses, cfg, mu_true, truth, method, opts = parts[-(1:2)], p = p),
    pseudorep = { l2 <- losses; l2$instance <- seq_len(nrow(l2)); l2$seed <- 1L
                  run_seqbench(parts[2], l2, cfg, mu_true, truth, method) },
    audit = run_audit(parts[2], losses, cfg, mu_true, truth, method),
    fixed = run_fixed(if (parts[2] == "final") cfg$n_max else as.integer(parts[2]), losses, cfg, mu_true, truth, method),
    savi_cs = run_savi(losses, cfg, mu_true, truth, method, p, opts = parts[-1]),
    stop("unknown method ", method))
}

run_seqbench <- function(boundary, losses, cfg, mu_true, truth, method, opts = character(0), p = NULL) {
  args <- list(alpha = cfg$alpha, margin = cfg$margin, bounds = cfg$bounds, boundary = boundary, n_max = cfg$n_max)
  for (o in opts) {
    if (grepl("^c[0-9.]+$", o)) args$betting_c <- as.numeric(sub("^c", "", o))
    if (grepl("^sd[0-9.]+$", o)) args$sd_max <- as.numeric(sub("^sd", "", o)) * cfg$savi_sd_of(p)
  }
  design <- do.call(seqbench::comparison_design, args)
  st <- suppressWarnings(seqbench::update_comparison(seqbench::initialize_comparison(design), losses))
  tr <- st$trajectory
  dec <- seqbench::stopping_decision(st)
  empty <- isTRUE(st$diagnostics$cs_empty)
  miscover <- if (is.na(mu_true)) NA else
    (empty || any(tr$lower_running > mu_true | tr$upper_running < mu_true, na.rm = TRUE))
  result_row(method, as.vector(dec), attr(dec, "reason"), nrow(tr), st$cost, miscover,
             tr$lower_running[nrow(tr)], tr$upper_running[nrow(tr)], truth, cs_empty = empty)
}

# No-stopping coverage audit over the whole stream. Uses the package kernels directly:
# closed-form intervals for hoeffding / empirical_bernstein; for betting the exact log
# capital at the true mean (no grid, O(1) per step), which is >= log(1/alpha) at some t
# iff the true mean leaves the running confidence sequence.
run_audit <- function(boundary, losses, cfg, mu_true, truth, method) {
  a <- cfg$bounds[1] - cfg$bounds[2]; rng <- 2 * (cfg$bounds[2] - cfg$bounds[1])
  x <- pmin(pmax((instance_diffs(losses) - a) / rng, 0), 1)
  nu <- (mu_true - a) / rng
  n <- length(x)
  if (boundary == "betting") {
    st <- seqbench::boundary_init("betting", alpha = cfg$alpha, grid = 11L, refine = FALSE)
    logk_p <- 0; logk_m <- 0; crossed <- FALSE; first_t <- NA_integer_
    for (t in seq_len(n)) {
      prev <- seqbench:::running_prev(st)
      lam <- sqrt(2 * log(2 / cfg$alpha) / (prev$sig2 * t * log(t + 1)))
      inc <- seqbench:::betting_log_increment(x[t], lam, nu, st$c)
      logk_p <- logk_p + inc$plus; logk_m <- logk_m + inc$minus
      if (!crossed && max(log(st$theta) + logk_p, log(1 - st$theta) + logk_m) >= log(1 / cfg$alpha)) {
        crossed <- TRUE; first_t <- t
      }
      st <- seqbench::boundary_update(st, x[t])
    }
    miss <- crossed
  } else {
    st <- seqbench::boundary_init(boundary, alpha = cfg$alpha)
    L <- 0; U <- 1; miss <- FALSE; first_t <- NA_integer_
    for (t in seq_len(n)) {
      st <- seqbench::boundary_update(st, x[t]); ci <- seqbench::boundary_interval(st)
      L <- max(L, ci[["lower"]]); U <- min(U, ci[["upper"]])
      if (!miss && (L > nu || U < nu)) { miss <- TRUE; first_t <- t }
    }
  }
  result_row(method, "audit", NA_character_, n, NA_real_, miss, NA_real_, NA_real_, NA_character_,
             extra = sprintf("first_miss_t=%s", first_t))
}

run_fixed <- function(n, losses, cfg, mu_true, truth, method) {
  d <- instance_diffs(losses)
  n <- min(n, length(d))
  d <- d[seq_len(n)]
  ci <- if (n >= 2) as.numeric(stats::t.test(d, conf.level = 1 - cfg$alpha)$conf.int) else c(-Inf, Inf)
  dec <- three_way(ci[1], ci[2], cfg$margin)
  if (dec == "continue") dec <- "inconclusive"
  miscover <- if (is.na(mu_true)) NA else (ci[1] > mu_true | ci[2] < mu_true)
  cost <- sum(losses$cost[losses$instance %in% unique(losses$instance)[seq_len(n)]])
  result_row(method, dec, if (dec == "inconclusive") "budget" else "declaration", n, cost, miscover, ci[1], ci[2], truth)
}

# Reference sd for the savi design: never from the evaluation stream itself.
savi_sd_ref <- function(cfg, p) {
  if (!is.null(cfg$savi_sd_of)) return(cfg$savi_sd_of(p))          # known DGP sd (oracle, favourable to safestats)
  key <- paste0("pilot_", paste(unlist(p), collapse = "_"))
  if (is.null(.savi_cache[[key]])) {                                 # independent pilot stream, fixed seed
    set.seed(cfg$seed + 777L)
    .savi_cache[[key]] <- stats::sd(instance_diffs(cfg$gen(p, cfg$savi_pilot_n %||% 200L)))
  }
  .savi_cache[[key]]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

run_savi <- function(losses, cfg, mu_true, truth, method, p = NULL, opts = character(0)) {
  d <- instance_diffs(losses)
  if ("pilot" %in% opts) cfg$savi_sd_of <- NULL            # savi_cs:pilot -- ignore the oracle sd, use the pilot stream
  sd_ref <- savi_sd_ref(cfg, p)
  key <- sprintf("design_%.8f", cfg$margin / sd_ref)
  if (is.null(.savi_cache[[key]])) {
    .savi_cache[[key]] <- safestats::designSaviT(deltaMin = cfg$margin / sd_ref, alpha = cfg$alpha,
                                                 alternative = "twoSided", testType = "oneSample",
                                                 eType = "mom", nSim = 100L, seed = 1, pb = FALSE)
  }
  res <- safestats::saviTTest(d, designObj = .savi_cache[[key]], sequential = TRUE, wantCi = TRUE)
  cs <- res$confSeqMatrix
  L <- cummax(ifelse(is.na(cs[, 1]), -Inf, cs[, 1])); U <- cummin(ifelse(is.na(cs[, 2]), Inf, cs[, 2]))
  decs <- vapply(seq_along(d), function(t) three_way(L[t], U[t], cfg$margin), character(1))
  stop_t <- which(decs != "continue")
  t <- if (length(stop_t)) stop_t[1] else length(d)
  dec <- if (length(stop_t)) decs[t] else "inconclusive"
  miscover <- if (is.na(mu_true)) NA else any(L[seq_len(t)] > mu_true | U[seq_len(t)] < mu_true)
  ev <- res$eValueVec; ev[is.na(ev)] <- 0        # NaN e-values when the first differences all tie (sd 0)
  e_reject <- any(ev >= 1 / cfg$alpha)
  cost <- sum(losses$cost[losses$instance %in% unique(losses$instance)[seq_len(t)]])
  result_row(method, dec, if (dec == "inconclusive") "budget" else "declaration", t, cost, miscover, L[t], U[t], truth,
             extra = sprintf("e_reject=%s;fpt=%s;sd_ref=%.5f;na_e=%d", e_reject, if (e_reject) which(ev >= 1 / cfg$alpha)[1] else NA, sd_ref, sum(is.na(res$eValueVec))))
}
.savi_cache <- new.env()
