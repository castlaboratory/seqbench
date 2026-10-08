# experiments/

Reproducible experiments of the `seqbench` paper. Nothing here is part of the R package.

```
configs/    frozen configurations (one file per study, dated; never edited after freezing)
run/        methods.R (procedures under comparison, one result row each), run_study.R (runner)
analysis/   scripts that turn experiments/out/<study>/runs.csv into tables and figures
out/        results (gitignored): <study>/{config.rds, session.txt, progress.log, cells/, runs.csv}
```

Run a study (resumable; each (data cell, method) is saved separately and skipped on rerun):

```sh
Rscript experiments/run/run_study.R experiments/configs/main_sim_2026-09-25.R 8      # 8 cores
Rscript experiments/run/run_study.R experiments/configs/main_sim_2026-09-25.R 4 3    # smoke test, 3 reps
```

Every replication seeds the data stream from `(seed, dgp cell, rep)` only, so all methods see
the same paired losses; with `cache_data = TRUE` the streams are generated once, in parallel,
before any method runs. The truth per data cell (`mu_of`, possibly a Monte Carlo reference) is
computed once in the parent process. `session.txt` records R, package versions, host and start time.
Long runs go on the CAST server (see the lab notes); the pilot of 2026-09-25 lives in
`run/run_pilot.R` and `analysis/pilot_summary.R` with its own config.

Studies frozen on 2026-09-25:

| Config | What | Cells × reps |
|---|---|---|
| `main_sim_2026-09-25.R` | simulation: null, equivalence, at-margin, moderate/large effects; two noise levels; unpaired and seed-cluster ablations; references `fixed:{100,500,final}` and `savi_cs` | 24 × 8 methods × 2000 |
| `app_estimators_2026-09-25.R` | mean vs 20 % trimmed mean on contaminated-normal samples (scaled absolute error, capped) | 1 × 4 × 500 |
| `app_louvain_leiden_2026-09-25.R` | Louvain vs Leiden (`igraph`) on random stochastic block models, loss `1 - Q`, margin 0.01 | 1 × 4 × 100 |
| `app_knapsack_2026-09-25.R` | greedy vs GRASP-like heuristic on random 0/1 knapsack (truncated LP gap); dependency-free | 1 × 3 × 100 |

Amendments frozen after the external review of 2026-09-26 and the c sweep of 2026-09-27
(analysis: `analysis/amendments_summary.R` for the first two, writes `paper/rjournal/data/T7_extra.csv`
and `T8_cluster_v2.csv`):

| Config | What | Cells × reps |
|---|---|---|
| `main_sim_extra_2026-09-26.R` | negative effects and effects 0.005/0.01 from the margin, same DGP and methods as the main study. **Cell dgp 8 (`mu = -0.20`, `h = 0.173`) is invalid by construction**: `loss_b` reaches 1.073 > 1, the package refuses every replication (never clips); dropped and reported | 14 × 8 × 2000 |
| `main_sim_cluster_2026-09-26.R` | seed clusters v2: instance-specific differential effect `v_i ~ U(-w, w)`, `w ∈ {0, 0.05}`, 1 or 5 seeds, with the invalid `pseudorep:betting` analysis | 12 × 4 × 2000 |
| `audit_coverage_2026-09-26.R` | no-stopping coverage audit (`audit:*`) on the 24 main cells | 24 × 3 × 2000 |
| `csweep_2026-09-27.R` | betting cap `c ∈ {0.75, 0.9, 0.99}` and `bernstein_declared` at 1×, 2×, 0.5× the true sd on the 12 baseline cells | 12 × 6 × 2000 |

Paper (ii) studies frozen on 2026-09-28 (launched the same day: `sigma_grid` and `cluster_icc`
on cast, `robust` split — non-safestats methods on cast, `savi_cs*` on the Mac with
`SEQBENCH_METHODS="^savi"`; the two halves are merged by copying `cells/` and rerunning the runner):

| Config | What | Cells × reps |
|---|---|---|
| `robust_2026-09-28.R` | skewed (Beta, skew 0 / 1.2 / 1.9) and tied (p0 ∈ {0, 0.4, 0.7}) paired differences on [0, 1]; betting c 0.5 and 0.9, EB, fixed 2000, `savi_cs` (oracle sd) and `savi_cs:pilot` | 18 × 6 × 2000 |
| `sigma_grid_2026-09-28.R` | six sd levels (0.01 to 0.15) × μ ∈ {0, 0.05, 0.10} for the two-regime approximation; betting c 0.5 / 0.9, EB, `bernstein_declared:sd1` | 18 × 4 × 2000 |
| `cluster_icc_2026-09-28.R` | seed clusters with strong differential effect (w = 0.15, ICC 0.41, design effect 2.7): correct analysis vs pseudoreplication | 2 × 2 × 2000 |
