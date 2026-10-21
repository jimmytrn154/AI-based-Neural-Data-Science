# Project Progress Tracker

BME54701 team project: brain-state reorganization from rest to task in Parkinson's disease.
See `README.md` for the plan and `docs/proposal.tex` for the full proposal.

**Status key:** ⬜ not started · 🟨 in progress · ✅ done · 🟥 blocked · ⏸️ postponed

**How to use this file:**

- Update your own rows whenever a task's status changes.
- Add one line to the daily log each day you work on the project.
- Record every change to a shared rule (exclusion, filters, K, features) in the decision log, with the date and who agreed.

---

## Phase overview

| Phase | Window | Status | Notes |
|---|---|---|---|
| 1. Data structure (unsupervised) | Oct 20 → Oct 26 | 🟨 | Progress presentation on Oct 26 |
| 2. Consolidation | After Oct 26 | ⬜ | To be planned after professor's feedback |
| 3. Prediction | Mid-semester | ⬜ | Placeholder |
| 4. Extensions | Late semester | ⬜ | Placeholder |
| Final presentation | Week 16 | ⬜ | |

---

## Phase 1: data structure (Oct 20 → Oct 26)

**Goal:** preliminary unsupervised results showing how microstates vary between rest and the Simon task, and between congruent and incongruent conditions. The presentation opens with the research questions.

### Milestones

| Date | Milestone | Owner | Status |
|---|---|---|---|
| Oct 20 | Shared rules and interfaces agreed; repo and `config.yaml` set up | All | ✅ |
| Oct 21 | **Draft templates released** (~15 subjects, K = 4) | Chuong | ✅ |
| Oct 22 | **Final templates released** (full cohort, K chosen) | Chuong | ⬜ |
| Oct 23 | Task-level metrics delivered for Rest and Simon | Chuong | ⬜ |
| Oct 24 | All feature tables merged; PD/MoCA associations run; figures chosen | All | ⬜ |
| Oct 25 | Slides finished and rehearsed | All | ⬜ |
| Oct 26 | **Progress presentation** | All | ⬜ |

### 0. Setup (all)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Create shared repo with the agreed layout (`src/`, `notebooks/`, `results/`, `docs/`) | Chuong | Oct 20 | ✅ | |
| Write `config.yaml` (DATA_ROOT, exclusion rules, filters, K, random seed) | Chuong | Oct 20 | ✅ | |
| Everyone installs the environment (mne, pycrostates, h5py, pandas, statsmodels, scikit-learn, umap-learn) and loads one subject | All | Oct 20 | 🟨 | Chuong ✅ (`.venv`, `requirements.txt`) |
| Fix and write down the exclusion rule **before looking at results** | All | Oct 20 | ✅ | >6 interpolated channels; `trials_simon == exclude` |
| Agree on file formats: templates, label sequences, feature CSVs keyed by `subject_id` | All | Oct 20 | ✅ | |
| Agree on the statistics template (tests, covariates, FDR) | All | Oct 20 | ✅ | |

### 1. Behavior and sanity checks (Minseo)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Reproduce Simon RT, correct trials (514.7 / 557.4 ms) | Minseo | Oct 21 | ⬜ | |
| Reproduce Simon error rates (6.4% / 14.8%) | Minseo | Oct 21 | ⬜ | |
| Simon effect by group (RT and errors) | Minseo | Oct 21 | ⬜ | |
| Behavioral feature table (`results/features/behavior.csv`) | Minseo | Oct 21 | ⬜ | |

### 2. Templates and task-level dynamics (Chuong)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Shared loading and preprocessing functions (1–40 Hz, 250 Hz, exclusion) | Chuong | Oct 21 | ✅ | `src/io.py`, `src/preprocess.py`; n: Rest 149, Simon 141, both 141 |
| GFP-peak extraction, with equal peaks per subject and task | Chuong | Oct 21 | ✅ | `ms.subject_peaks`; outliers > mean + 3 SD of peak GFP dropped per recording |
| Draft templates (~15 subjects, K = 4) | Chuong | Oct 21 | ✅ | Unblocks Seonwoo, Minseo, Bohyeon. 10 PD + 5 C; template GEV 44.1% ± 12.5%; `results/figures/chuong_templates_K4_draft.png`. A–D labels not assigned yet |
| Final group templates (fit K = 4–7; report K = 4 and 5) | Chuong | Oct 22 | ⬜ | |
| Task-specific templates and spatial-correlation check | Chuong | Oct 22 | ⬜ | |
| Back-fitting function shared in `src/microstates.py` | Chuong | Oct 21 | ✅ | `ms.load_templates(4, cfg)`, `ms.backfit(inst, templates, cfg)`; epochs via `pio.epochs_to_mne` |
| Rest and Simon metrics: duration, occurrence, coverage, GEV, transitions | Chuong | Oct 23 | ⬜ | `results/features/microstates_task.csv` |

### 3. Condition-level dynamics (Seonwoo)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Load Simon cue-locked epochs; drop bad epochs; keep correct trials | Seonwoo | Oct 21 | ⬜ | |
| Reproduce the Simon ERP at Cz (Dataset Guide, Figure 2) | Seonwoo | Oct 21 | ⬜ | Sanity check |
| Topographic segmentation of each condition's ERP (0–800 ms) | Seonwoo | Oct 22 | ⬜ | Doesn't need templates |
| Single-trial back-fitting and state-probability time courses | Seonwoo | Oct 23 | ⬜ | Needs draft/final templates |
| Congruent vs. incongruent: cluster-based permutation test | Seonwoo | Oct 24 | ⬜ | |
| Conflict-effect feature in 200–350 ms and 400–750 ms windows | Seonwoo | Oct 24 | ⬜ | `results/features/conflict.csv` |

### 4. Comparisons, networks and clustering (Bohyeon)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Statistics functions on dummy data (mixed models, partial Spearman, FDR) | Bohyeon | Oct 21 | ⬜ | Doesn't need templates |
| Rest vs. Simon within-subject comparisons | Bohyeon | Oct 23 | ⬜ | Needs task-level metrics |
| Reorganization index (Jensen–Shannon divergence) and coverage change | Bohyeon | Oct 23 | ⬜ | `results/features/reorganization.csv` |
| Transition networks (rest vs. Simon) | Bohyeon | Oct 24 | ⬜ | |
| Subject embedding (PCA, UMAP) and k-means clustering; alignment with group and MoCA | Bohyeon | Oct 24 | ⬜ | |
| PD vs. control and MoCA associations for all primary features | Bohyeon | Oct 24 | ⬜ | Age and sex as covariates |

### 5. Recovery pilot (Minseo, Bohyeon; exploratory)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Define the recovery measure and baseline window (response-locked epochs) | Minseo | Oct 22 | ⬜ | |
| Post-incongruent vs. post-congruent (and post-error) state-probability recovery | Minseo | Oct 24 | ⬜ | |
| PCA trajectories of response-locked activity | Bohyeon | Oct 24 | ⬜ | |
| Split-half reliability (odd vs. even trials, Spearman–Brown) | Minseo | Oct 24 | ⬜ | Go/no-go input for Phase 2 |

### 6. Robustness checks (owner of each analysis)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Rest results without recordings with `nICremoved == 0` | Chuong | Oct 24 | ⬜ | Blink sensitivity |
| Key results at K = 4 and K = 5 | Each analysis owner | Oct 24 | ⬜ | |

### 7. Presentation (all)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Slide outline: research questions → data → methods → results → plan | All | Oct 24 | ⬜ | |
| Results slides, one or two per analysis | Each analysis owner | Oct 25 | ⬜ | |
| Slide on sanity checks (reference values reproduced) | Minseo | Oct 25 | ⬜ | |
| Slide on next phases | Chuong | Oct 25 | ⬜ | |
| Rehearsal and timing check | All | Oct 25 | ⬜ | |

### Phase 1 results summary

*Fill in as results arrive. Keep it to one or two lines each.*

| Question | Result | Figure | Confidence |
|---|---|---|---|
| Sanity checks reproduced? | | | |
| Chosen K and explained variance | | | |
| Rest vs. Simon templates similar? | | | |
| RQ1: rest vs. Simon dynamics | | | |
| RQ2: congruent vs. incongruent | | | |
| RQ3: PD vs. control, MoCA, errors | | | |
| Recovery measurable? | | | |
| Subject clusters vs. clinical groups | | | |

### Questions to bring to the professor on Oct 26

-
-

### Professor's feedback (Oct 26)

*Record the feedback here right after the presentation. It drives the Phase 2 plan.*

-
-

---

## Phase 2: consolidation (placeholder)

To be planned after the Oct 26 feedback. Planned so far:

- [ ] Refine research questions from the feedback
- [ ] Robustness: K, blinks, filter choices
- [ ] Alpha control with spectral parameterization (Rest)
- [ ] Fix the final feature set
- [ ] Go/no-go decision on the recovery analysis

## Phase 3: prediction (placeholder)

- [ ] Prediction study: demographics → behavior → Rest EEG → Rest + reorganization; EEG only
- [ ] Calibration (Brier score, reliability), within-PD analysis, selective prediction

## Phase 4: extensions (placeholder)

- [ ] Generalize to Oddball and Timing
- [ ] Deep learning comparison (e.g., EEGNet)
- [ ] Transition networks as brain networks
- [ ] Optional: novelty P3; task-to-rest transfer

## Final presentation, Week 16 (placeholder)

- [ ] Final slides
- [ ] Prompt histories collected for submission

---

## Decision log

| Date | Decision | Agreed by |
|---|---|---|
| Oct 20 | Adopt "Rank 1" direction: rest-to-Simon microstate reorganization, with Phase 1 unsupervised only | All |
| Oct 20 | Recovery after conflict included as exploratory pilot (from "Rank 2") | All |
| Oct 20 | Spectral parameterization moved to Phase 2 as alpha control; prediction moved to Phase 3 | All
| Oct 20 | Shared rules adopted as in config.yaml (no objections) | All |
| | | |

## Blockers

| Date raised | Blocker | Who is blocked | Owner | Resolved? |
|---|---|---|---|---|
| Oct 21 | Simon epoch files don't carry a readable trial table (MATLAB `table`, not decodable by h5py). Checked all 147 subjects: `lock-cue` = one epoch per `trials.tsv` row; `lock-resp` = one epoch per row with non-missing `onset_response`. `beh_accuracy == 99` rows have a timeout `onset_response` (`rt_eeg` ≈ 1.0 s), so drop them in response-locked analyses. Use `src.io.match_trials`. Open point: epoch order assumed equal to TSV order (counts match; order not verifiable without MATLAB) | Minseo (recovery), Seonwoo | Chuong | Resolved with rule; order check open |
| Oct 21 | sub-P048 (draft set) looks artifact-dominated: individual maps peak at AF8/Fp2, FC5/T7, TP10/I1; template GEV 16% vs 44% mean. No config rule covers this | Chuong (final templates) | Chuong | No, needs team decision |
| Oct 21 | Draft template states not yet matched to canonical A–D (config: visual inspection, fixed once). Needed before `microstates_task.csv` column names | Chuong, all | Chuong | No |

## Daily log

| Date | Who | What was done | Next |
|---|---|---|---|
| Oct 20 | | | |
| Oct 21 | Chuong | Milestone 1: loading/preprocessing in `src/`, exclusion logged (Rest 149→149; Simon 149→147→141 after `trials_simon`; nInterp > 6 removed none; 35/149 Rest flagged `nICremoved == 0`); env + `requirements.txt` | GFP peaks and draft templates |
| Oct 21 | Chuong | Milestone 2: GFP peaks + two-level modified k-means; draft K = 4 templates (15 subjects, 10 PD / 5 C) saved to `results/templates/templates_K4.fif`; group GEV 68.1% on individual maps, 44.1% ± 12.5% on subjects' GFP peaks; `backfit` tested on Raw and Epochs; `match_trials`, `epochs_to_mne` helpers | A–D labels, P048 decision, then task-level metrics |
