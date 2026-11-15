# Bohyeon's Phase 2 analyses (delivered Nov 12, 2026)

Bohyeon's code and outputs for the Phase 2 presentation (Nov 9 deck), integrated as delivered. Reviewed by Chuong on Nov 15: every slide number checked below is in these outputs, and the key template numbers were recomputed from the repo. The folder keeps its original name and layout for now.

**Naming.** The deck's "B7" and "B8" are new Phase 2 analyses. In the docs they are called **P2-B7** (conflict state contrast) and **P2-B8** (conflict group windows), to keep them apart from Bohyeon's Phase 1 tasks B7 (robustness) and B8 (slides). The script names are unchanged.

## Contents

```
bohyeon 261112/
├── src/        analysis scripts B2–B6 (Phase 1 tasks), B7_/B8_ scripts (= P2-B7, P2-B8), shared helpers
├── results/    one folder per run: task, task_qcincluded, nosmooth, rest_noblink, K4
│               + summary_all_runs.csv (headline numbers per run)
└── figures/    same run folders as results/
```

`src/__pycache__/` is git-ignored. It holds a stale `variants.cpython-313.pyc` with no matching source; no script imports `variants`, so nothing is missing.

## Scripts

All scripts read from a private `data/` folder next to `src/` (copies of the repo's feature files plus `subjects.csv`), not from the repo's `results/features/`. Paths need adapting before they run inside the repo.

| Script | What it does | Inputs |
|---|---|---|
| `run_pipeline.py` | Runs B3 → B2 → B4 → B5 → B6 for each run listed in `RUNS`; writes `results/<run>/…_<run>.csv`, `figures/<run>/…`, a run log and a run config, and updates `results/summary_all_runs.csv` | One microstate CSV per run; `conflict.csv` only for the `task` run |
| `pipeline_io.py` | Shared argument parsing, output paths, Simon trial-label filter (`trials_simon` in ok/partial), run-config writer. If a microstate file has only Rest columns (`microstates_rest_noblink.csv`), the Simon columns are taken from `microstates_task.csv` | — |
| `stats.py` | Shared statistics (FDR by family, paired Wilcoxon, mixed model with Powell fallback, partial Spearman, group model) | — |
| `reorganization.py` (B3) | Reorganization index (Jensen–Shannon divergence of transitions, mean over source states) and coverage change per state | Microstate CSV |
| `B2_rest_vs_simon.py` | Rest vs. Simon: paired Wilcoxon per metric and state (FDR within RQ1), mixed models, PD vs. control on the changes | Microstate CSV, `subjects.csv` |
| `B4_transition_network.py` | Rest and Simon transition networks (2 strongest outgoing edges per state) and the Simon − Rest difference matrix | Microstate CSV |
| `B5_subject_clustering.py` | Standardized features → PCA, k-means (k = 2–10 by silhouette); group and MoCA compared afterwards | Microstate CSV, reorganization output, `subjects.csv` |
| `B6_clinical_associations.py` | RQ3: group models (age, sex) and partial Spearman with MoCA, Simon error rate, UPDRS (PD only); one FDR family | Reorganization, `conflict.csv`, microstate CSV, `behavior.csv`, `subjects.csv` |
| `B7_conflict_state_contrast.py` (**P2-B7**) | Per window (N2 200–350 ms, theta 400–750 ms): D = Δ(E) − Δ(C), where Δ = P(state \| incongruent) − P(state \| congruent). Paired test (subject fixed effects), Wilcoxon, and PD vs. control on D; FDR across the two windows within each test | `conflict.csv` (Seonwoo, 137 subjects); group from the subject ID. Age/sex only if `--clinical` is passed (it was **not** for these results) |
| `B8_conflict_group_windows.py` (**P2-B8**) | Δ(C), Δ(E) and D per window and group: each vs. 0 (12 tests), PD vs. control adjusted for age and sex (6 tests), group × window interaction (3 tests); three FDR families | `conflict.csv`, `subjects.csv` |

## Runs

All runs use the K = 5 templates (states A–E) except `K4`, and the Simon trial-label filter.

| Run | Microstate file | Scripts run | Subjects (Rest–Simon pairs) |
|---|---|---|---|
| `task` | `microstates_task.csv` (primary) | B3, B2, B4, B5, B6, P2-B7, P2-B8 | 139 |
| `task_qcincluded` | `microstates_task_qcincluded.csv` | B3, B2, B4, B5 | 141 |
| `nosmooth` | `microstates_task_nosmooth.csv` | B3, B2, B4, B5 | 139 |
| `rest_noblink` | `microstates_rest_noblink.csv` (Simon columns from the primary file) | B3, B2, B4, B5 | 106 |
| `K4` | `microstates_task_K4.csv` (A, B, C, Central) | B3, B2, B4, B5 | 138 |

B6 runs only when a conflict file is given, so `all_steps_ok` is False for the four sensitivity runs; their logs show no errors.

## Input provenance (checked Nov 15)

The private `data/` folder was not delivered, so the inputs were checked by value instead of by checksum. For all five runs, the per-subject Simon − Rest changes (`rq1_individual_changes_<run>.csv`) equal those computed from the repo's `results/features/` files (largest difference 3e-14; 5e-9 for K4, CSV rounding), with the same subjects and the same missing values. The conflict values used by P2-B8 equal `conflict.csv` exactly, and B6's conflict features match it to 1e-16. `subjects.csv` and `behavior.csv` also match: the group, age and sex in Bohyeon's local P2-B8 table equal `subjects.csv` for all 136 subjects, and all 90 B6 partial correlations (MoCA, `simon_err_inc`, UPDRS) and all 18 B6 group models reproduce from the repo's `subjects.csv` and `results/features/behavior.csv` with her `stats.py` (differences ≤ 1e-13).

## Numbers checked against the Nov 9 deck

All found in these outputs (or, for the template checks, recomputed from `results/templates/`):

- **RQ1** (`task/rq1_stats_task.csv`, N = 139), mean coverage change Simon − Rest: A +2.20, C −3.00, E +1.80 points. E +1.45 with no smoothing (`nosmooth`), +1.76 with QC-excluded recordings included (`task_qcincluded`, N = 141).
- **P2-B7** (`task/B7_conflict_state_contrast.csv`, N = 137): D = +0.67 points in N2 (FDR p = 0.214) and +1.23 in theta (FDR p = 0.00724). Wilcoxon: N2 0.479, theta 0.0945; these are **FDR-adjusted** values (raw: 0.479, 0.047). PD − control in D: N2 −2.80 points (FDR p = 0.0291), theta −0.86 (0.349); these group results are **not adjusted for age and sex**.
- **P2-B8** (`task/b8_adjusted_tests_task.csv`, N = 136, age and sex adjusted), PD − control: N2 Δ(C) +2.46 points (FDR p = 0.00513), Δ(E) −0.43 (0.611), D −2.90 (0.0413); theta Δ(C) −0.17 (0.762), Δ(E) −1.35 (0.0443), D −1.17 (0.363). Group × window: C 0.00203, E 0.245, D 0.131.
- **E validation** (recomputed from the templates): K5 E–C \|r\| = 0.337, C–D 0.796; K4 Central–K5 E 0.986, K5 E–K6 F 0.960, K6 F–K7 G 0.948; Rest-only E vs. Simon-only E 0.948, vs. shared E 0.984 (Rest-only) and 0.983 (Simon-only). At K = 7, K5 E correlates with both F (0.82) and G (0.83).

**Why P2-B7 has 137 subjects and P2-B8 has 136:** sub-P069 (PD) has sex recorded as "F\|M". P2-B8 adjusts for sex and, following the config, treats conflicting sex values as missing, so P069 is left out; P2-B7's group model is unadjusted and keeps it.

## Caveats

- **P2-B7 group test is unadjusted.** Rerun with `--clinical subjects.csv` for an age/sex-adjusted version comparable with P2-B8.
- **Blink sensitivity:** in `rest_noblink` (N = 106), E (+1.58 points, FDR p = 0.019) and A (+2.12) hold, but C's decrease (−2.41 points) is no longer FDR-significant (p = 0.13). The `rest_noblink` run is the primary data restricted to those 106 subjects (identical values); see `PROGRESS.md`, P2.7 note. B6, P2-B7 and P2-B8 were not rerun on the sensitivity files.
- **Not committed:** `results/task/b8_subject_window_features_task.csv` stays local (untracked) because it holds per-subject group, age and sex copied from `subjects.csv` (dataset content; see the data rules in `README.md`). P2-B8 rebuilds it from `conflict.csv` and `subjects.csv`.
- **Normative-template comparison** (Zanesco et al. 2020, \|r\| = 0.868; Zanesco et al. 2021a, \|r\| = 0.843): not in this folder or anywhere in the repo. Results reported by a teammate; code not yet in the repo.
