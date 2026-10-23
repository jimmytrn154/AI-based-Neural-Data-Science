# Team Tasks: who does what, with which files, producing what

BME54701 team project: brain-state reorganization from rest to task in Parkinson's disease.

This file describes **each member's tasks, the inputs they need, and the outputs they must deliver**. It complements:

- `PROGRESS.md` (repo root): status, milestones, decision log, blockers, daily log;
- `README.md`: project plan and dataset reference;
- `config.yaml`: all analysis rules (never change it without team agreement);
- `docs/proposal.tex`: full proposal;
- `docs/SETUP.md`: how to download the project from GitHub, install it and run an example (no git needed).

**Rule of thumb:** this file says *what* each person delivers and *from which inputs*; `PROGRESS.md` says *how far along* it is. When a task changes owner or scope, update both.

**Status key:** ⬜ not started · 🟨 in progress · ✅ done · 🟥 blocked · ⏸️ postponed

*Last updated: Oct 23, 2026.*

---

## Reassignment, Oct 23

Minseo is unavailable until **Saturday, Oct 24**. To keep the Oct 24 merge on schedule:

| Task | Was | Now | Why |
|---|---|---|---|
| Behavioral sanity checks and `behavior.csv` | Minseo | **Chuong** | Bohyeon needs Simon error rates for the RQ3 associations on Oct 24 |
| Recovery pilot: measure, recovery comparison, split-half reliability | Minseo | **Chuong** | Reuses Chuong's response-locked back-fitting code (state E time course) |
| PCA trajectories for the recovery pilot | Bohyeon | **Chuong (optional)** | Frees Bohyeon for the core RQ1/RQ3 analyses |
| Sanity-check slide, results summary, slide assembly, questions for the professor | — | **Minseo** (from Oct 24) | Light tasks that fit her availability |

Minseo's Phase 3 role (prediction and calibration lead) is unchanged.

---

## Shared resources (everyone)

| Resource | Where | Notes |
|---|---|---|
| Dataset | `paths.data_root` in `config.yaml` (local copy, never committed) | Full `Dataset_ALL`, including epoch files |
| Rules and settings | `config.yaml` | Exclusions, filters, K, smoothing, statistics, naming |
| Primary K | `config.yaml` → `microstates.primary_k: 5` | States **A, B, C, D, E**. Always read K from the config; never hard-code it |
| Templates | `results/templates/templates_K{primary_k}.fif` | On GitHub (download from `results/templates/`). K = 4 file is sensitivity only (4th map = "Central") |
| Loading | `src/io.py` | `load_subjects`, `load_continuous` (`.set` files), `load_epochs` (epoch `.mat` files), `load_trials`; `match_trials` aligns epochs with trial tables; `epochs_to_mne` converts epochs |
| Preprocessing | `src/preprocess.py` | `select_subjects` (exclusion rules), `preprocess_raw` (filtering, resampling) |
| Microstates | `src/microstates.py` | `load_templates(k, cfg)`, `backfit(inst, templates, cfg)`, metrics (`segmentation_metrics`, `task_metrics`) |
| Primary metrics | `results/features/microstates_task.csv` | One row per `subject_id`; columns `<task>_<metric>_<state>` (e.g. `rest_cov_A`); transitions `<task>_trans_<X>_to_<Y>` |

**Before back-fitting any epochs:** filter at **2–20 Hz** (`config.yaml`, `simon_conditions.filter`) and resample to **250 Hz** (`preprocessing.resample_hz`), as done in `notebooks/chuong_state_E.py`. In response-locked analyses, drop trials with `beh_accuracy == 99` (no response). See the epoch-order blocker in `PROGRESS.md`.

---

## Tran Anh Chuong: pipeline, templates, metrics, behavior and recovery

**Role:** builds the shared foundation (pipeline, templates, metrics). Now also covers behavior and the recovery pilot.

### Completed outcomes (Phase 1)

| Outcome | Files | Key facts |
|---|---|---|
| Loading and preprocessing | `src/io.py`, `src/preprocess.py` | Eligible: Rest 149, Simon 141 |
| Templates, full cohort, K = 4–7 | `results/templates/templates_K{4..7}.fif`; figures `results/figures/chuong_templates_K{4..7}.png` | **Primary K = 5** (A–D canonical, r ≥ 0.97 with the draft; E = central maximum with opposite-sign temporal edges). Template GEV 56.8%. Choice of K: see decision log |
| Task-specific templates | `results/figures/chuong_templates_K5_rest.png`, `chuong_templates_K5_simon.png` | Rest-only and Simon-only both reproduce all five maps (\|r\| ≥ 0.98) |
| Recording QC | Decision log | K = 5 excludes sub-P001 Rest and sub-P048 Simon |
| Fast back-fitting | `src/microstates.py`, `tests/test_backfit.py` | Labels identical to pycrostates; about 2–7 s per subject |
| Task-level metrics | `results/features/microstates_task.csv` (primary), `_nosmooth`, `_qcincluded`, `_K4`, `microstates_rest_noblink.csv` | 149 rows; Rest 148, Simon 140 recordings. Coverage sums 100%; median segment 96 ms |
| State E description | `results/figures/chuong_E_timecourse.png` (code: `notebooks/chuong_state_E.py`) | Present at rest (≈11% coverage) and in Simon (≈12%); P(E) rises before the keypress and peaks at it. **Task-enhanced, not task-specific** |

### Remaining Phase 1 tasks

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| C1 | Reproduce Simon RT (514.7 / 557.4 ms) and error rates (6.4% / 14.8%) | Oct 23 | ✅ | `trial_info/simon/*.tsv`, `subjects.csv` | Numbers in the daily log; summary for Minseo's slide |
| C2 | Simon effect by group (RT and errors), descriptive | Oct 23 | ✅ | Same | Table: PD vs. control means |
| C3 | Behavioral feature table | Oct 23 | ✅ | Same | `results/features/behavior.csv` (see spec below) |
| C4 | Recovery pilot: measure, comparison, reliability | Oct 24 | ⬜ | Response- and cue-locked Simon epochs, templates, `backfit()` | `results/features/recovery.csv`, `results/figures/chuong_recovery.png`, go/no-go note |
| C5 | PCA trajectories of response-locked state probabilities (optional) | Oct 24 | ⬜ | C4 outputs | `results/figures/chuong_recovery_pca.png` |
| C6 | Slides: pipeline, K choice, templates, state E, next phases | Oct 25 | ⬜ | All of the above | Slides in the shared deck |

**`behavior.csv` specification:** one row per `subject_id`, Simon only.

- **Trials:** subjects with `trials_simon` in {ok, partial}. RTs come from correct trials only. `beh_rt` is in **milliseconds** for Simon.
- **Columns:**

  | Column | Meaning |
  |---|---|
  | `simon_rt_con`, `simon_rt_inc` | Mean RT, correct trials, per condition |
  | `simon_effect_rt` | incongruent − congruent RT |
  | `simon_err_con`, `simon_err_inc` | Error rate = `beh_accuracy == 2` / responded trials (`beh_accuracy` in {1, 2}) |
  | `simon_effect_err` | incongruent − congruent error rate |
  | `simon_noresp_rate` | `beh_accuracy == 99` / all trials |
  | `simon_n_con`, `simon_n_inc` | Trial counts |

- **Reference values** are computed as per-subject means first, then averaged across subjects.
- **No clinical columns**: those are merged from `subjects.csv` on `subject_id`.

**Recovery pilot specification (C4):** exploratory, and kept only if reliable.

1. **Baseline:** each subject's mean state-probability vector (A–E) in the cue-locked −200 to 0 ms window (correct trials).
2. **Recovery curve:** for response-locked epochs, the distance at each time point between the state-probability vector and that baseline, from 0 to +1000 ms after the keypress.
3. **Recovery measure:** the mean distance over +200 to +1000 ms (simple and robust; no threshold to choose).
4. **Comparisons, within subject:** incongruent vs. congruent correct trials; errors vs. correct trials. Report the number of error trials per subject, since it's small and makes the error comparison noisy.
5. **Reliability:** split-half (odd vs. even trials), Spearman–Brown corrected. **Keep the measure only if reliability ≥ 0.5** (`config.yaml`).
6. **No group or clinical analyses here.** Hand the per-subject measures to Bohyeon only if the pilot passes reliability.

---

## Seonwoo Kim: Simon condition-level dynamics (RQ2)

**Role:** when and how brain-state sequences differ between congruent and incongruent trials.

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| S1 | Load cue-locked epochs; keep correct trials with RT 150–1500 ms; drop bad epochs | Oct 23 | ⬜ (confirm with owner) | `epoched_data/simon/*_lock-cue_epo.mat`, `trial_info/simon/*.tsv`, `src/io.py` (`load_epochs`, `load_trials`, `match_trials`) | — |
| S2 | Reproduce the Simon ERP at Cz (Dataset Guide, Figure 2) | Oct 23 | ⬜ (confirm with owner) | S1 | `results/figures/seonwoo_erp_cz.png` (sanity check) |
| S3 | Topographic segmentation of each condition's ERP (0–800 ms) | Oct 24 | ⬜ (confirm with owner) | S1 | `results/figures/seonwoo_erp_segmentation.png` |
| S4 | Single-trial back-fitting and state-probability time courses (A–E), per condition | Oct 24 | ⬜ (confirm with owner) | S1, `templates_K5.fif`, `backfit()` (2–20 Hz, 250 Hz) | Per-subject time courses (saved as arrays) |
| S5 | Congruent vs. incongruent: paired cluster-based permutation test (5,000 permutations) | Oct 24 | ⬜ (confirm with owner) | S4 | `results/figures/seonwoo_conflict_timecourse.png` |
| S6 | Conflict-effect feature in the 200–350 ms and 400–750 ms windows | Oct 24 | ⬜ (confirm with owner) | S4 | `results/features/conflict.csv` (one row per `subject_id`) |
| S7 | Results slides for RQ2 | Oct 25 | ⬜ | S2–S6 | Slides in the shared deck |

**Deliverables for Oct 26:** the ERP sanity-check figure, the condition time courses with significant clusters marked, and `conflict.csv`.

**Tip:** `notebooks/chuong_state_E.py` (function `state_probabilities`) already back-fits cue- and response-locked epochs on a 20-subject subset and produced `chuong_E_timecourse.png`. Reuse it as a starting point.

---

## Bohyeon Rhee: rest vs. Simon, networks, groups and cognition (RQ1, RQ3)

**Role:** within-subject rest-vs-Simon comparisons, reorganization measures, transition networks, subject clustering, and the PD/MoCA associations.

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| B1 | Statistics functions (mixed models, partial Spearman, FDR) | Oct 23 | ⬜ (confirm with owner) | `config.yaml` (statistics) | `src/stats.py` |
| B2 | Rest vs. Simon within-subject comparisons (Wilcoxon, FDR, effect sizes); mixed models with group, age, sex | Oct 24 | ⬜ (confirm with owner) | `microstates_task.csv`, `subjects.csv` | `results/features/rq1_stats.csv`, figure |
| B3 | Reorganization index (Jensen–Shannon divergence of transitions) and coverage change | Oct 24 | ⬜ (confirm with owner) | `microstates_task.csv` (transition columns `<task>_trans_<X>_to_<Y>`) | `results/features/reorganization.csv` |
| B4 | Transition networks, Rest vs. Simon | Oct 24 | ⬜ (confirm with owner) | `microstates_task.csv` | `results/figures/bohyeon_transition_networks.png` |
| B5 | Subject embedding (PCA, UMAP) and k-means; alignment with group and MoCA | Oct 24 | ⬜ (confirm with owner) | `microstates_task.csv`, `subjects.csv` | Figure + short summary |
| B6 | PD vs. control and MoCA associations for the primary features (age, sex covariates; UPDRS within PD only; FDR) | Oct 24 | ⬜ (confirm with owner) | B2, B3, `behavior.csv` (from Chuong), `conflict.csv` (from Seonwoo) | `results/features/rq3_stats.csv` |
| B7 | Robustness: rerun the key RQ1/RQ3 results with `_K4`, `_nosmooth`, `_qcincluded` and `microstates_rest_noblink.csv` | Oct 24–25 | ⬜ | The sensitivity files | Short table: does each result hold? |
| B8 | Results slides for RQ1 and RQ3 | Oct 25 | ⬜ | B2–B7 | Slides in the shared deck |

**Deliverables for Oct 26:** the RQ1 figure (rest vs. Simon per state), the reorganization index, the transition networks, and the preliminary group and MoCA associations.

**Notes:**
- Descriptive preview (Chuong, Oct 23; paired Simon − Rest medians over 139 subjects, uncorrected, no tests): coverage of E +1.3 points (68% of subjects increase) and A +2.1 (60%); C coverage −1.8 (37% increase) and C mean duration −8.3 ms. The formal test is B2.
- sub-C049 has very high C coverage but passes QC. Keep it, and mention it in robustness checks if it drives any result.

---

## Minseo: presentation support now, prediction lead in Phase 3

**Role in Phase 1 (from Saturday, Oct 24):** light tasks that keep the presentation consistent.

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| M1 | Install the environment and load one subject | Oct 24 | ⬜ | `requirements.txt`, `README.md` | Confirmed in the daily log |
| M2 | Sanity-check slide: reproduced reference values vs. Dataset Guide | Oct 25 | ⬜ | Chuong's C1–C2 numbers (`docs/behavior_sanity_checks.md`) | One slide |
| M3 | Fill the "Phase 1 results summary" table in `PROGRESS.md` from everyone's results | Oct 25 | ⬜ | All outputs | Updated table |
| M4 | Assemble the slide deck: consistent style, terminology (A–E, "task-enhanced"), numbers matching the files | Oct 25 | ⬜ | Everyone's slides | Final deck |
| M5 | Draft "Questions to bring to the professor" (in `PROGRESS.md`) | Oct 25 | ⬜ | Team input | 2–4 questions |
| M6 | After Oct 26: draft the Phase 3 prediction and calibration design | After Oct 26 | ⏸️ | Proposal, Phase 3 section | Design note |

---

## Timeline to the presentation

| Date | Chuong | Seonwoo | Bohyeon | Minseo |
|---|---|---|---|---|
| Fri Oct 23 | C1–C3 (behavior) | S1–S2 | B1–B3 | — |
| Sat Oct 24 | C4 (recovery), C5 optional | S3–S6 | B4–B7 (merge: needs `behavior.csv`, `conflict.csv`) | M1 |
| Sun Oct 25 | C6 slides | S7 slides | B8 slides | M2–M5; rehearsal (all) |
| Mon Oct 26 | **Presentation (all)** | | | |

## Handoff rules

Teammates use GitHub **read-only**: browse and download, no git and no edits to repo files. See `docs/SETUP.md`.

- **Tasks and status:** read your tasks in `docs/TEAM_TASKS.md` and the current status in `PROGRESS.md`, on the GitHub page.
- **Inputs:** download them from GitHub (single files, or Code → Download ZIP). When Chuong announces an update, download them again.
- **Outputs:** send your outputs (CSVs, figures, slides) **and** the code or notebooks that produced them to Chuong by Drive or chat, using the file names in the task tables.
- **Integration:** Chuong commits them to the repo and updates `PROGRESS.md` and this file. Teammates don't edit repo files.
- **Feature tables:** one row per `subject_id`, columns named `<task>_<measure>_<state>`, no clinical columns; they go in `results/features/`.
- **Figures:** `results/figures/<owner>_<analysis>.png`.
- **When you finish a deliverable:** send it to Chuong as above and post the file names in the group chat.
- **When blocked:** tell Chuong in the group chat, with the input you need and from whom; Chuong records it in the Blockers table of `PROGRESS.md`.
