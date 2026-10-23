# Project Progress Tracker

BME54701 team project: brain-state reorganization from rest to task in Parkinson's disease.
See `README.md` for the plan and `docs/proposal.tex` for the full proposal.

**Status key:** ⬜ not started · 🟨 in progress · ✅ done · 🟥 blocked · ⏸️ postponed

**How to use this file:**

- Teammates read this file on GitHub; Chuong updates it from what teammates send (see `docs/TEAM_TASKS.md`, Handoff rules).
- Add one line to the daily log each day someone works on the project.
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
| Oct 22 | **Final templates released** (full cohort, K chosen) | Chuong | ✅ (Oct 23; primary K = 5) |
| Oct 23 | Task-level metrics delivered for Rest and Simon | Chuong | ✅ |
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

### 1. Behavior and sanity checks (Chuong; reassigned from Minseo, Oct 23)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Reproduce Simon RT, correct trials (514.7 / 557.4 ms) | Chuong | Oct 23 | ✅ | 514.4 / 557.4 ms (141 subjects; per-subject means first) |
| Reproduce Simon error rates (6.4% / 14.8%) | Chuong | Oct 23 | ✅ | 6.44% / 14.82% with errors / all trials (the Guide's definition); `behavior.csv` spec uses errors / responded trials: 6.84% / 16.27% |
| Simon effect by group (RT and errors) | Chuong | Oct 23 | ✅ | Descriptive: RT effect PD +47 / C +35 ms (Guide +50 / +34); incongruent errors PD 18.0% / C 8.5% (Guide 17.6 / 8.4). `docs/behavior_sanity_checks.md` |
| Behavioral feature table (`results/features/behavior.csv`) | Chuong | Oct 23 | ✅ | 141 rows, 9 columns per `docs/TEAM_TASKS.md` spec; no clinical columns |

### 2. Templates and task-level dynamics (Chuong)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Shared loading and preprocessing functions (1–40 Hz, 250 Hz, exclusion) | Chuong | Oct 21 | ✅ | `src/io.py`, `src/preprocess.py`; n: Rest 149, Simon 141, both 141 |
| GFP-peak extraction, with equal peaks per subject and task | Chuong | Oct 21 | ✅ | `ms.subject_peaks`; outliers > mean + 3 SD of peak GFP dropped per recording |
| Draft templates (~15 subjects, K = 4) | Chuong | Oct 21 | ✅ | Unblocks Seonwoo, Minseo, Bohyeon. Refit Oct 22 at 2–20 Hz with sub-P048 replaced by sub-P009: 10 PD + 5 C; template GEV 53.7% ± 9.0% (group level 75.5%); labels A–D **provisional**; `results/figures/chuong_templates_K4_draft.png` |
| Final group templates (fit K = 4–7; report K = 4 and 5) | Chuong | Oct 22 | ✅ | **Primary K = 5** (team decision Oct 23; rule gave K = 4). 139 subjects after QC; template GEV K4 54.8%, K5 56.8%, K6 58.4%, K7 59.3%. K = 5 labels A–E visually confirmed (A–D r ≥ 0.97 with draft; E = central maximum with opposite-sign temporal edges). K = 4 kept as sensitivity with 4th map 'Central' (C includes canonical D); K = 6, 7 sensitivity only. `results/figures/chuong_templates_K{4..7}.png` |
| Task-specific templates and spatial-correlation check | Chuong | Oct 22 | ✅ | K = 5: Rest-only vs shared \|r\| A .99, B .99, C .99, D 1.00, E .98; Simon-only A 1.00, B .99, C 1.00, D 1.00, E .98. (K = 4 sensitivity: Rest-only lacked the 'Central' map, r = .39) |
| Back-fitting function shared in `src/microstates.py` | Chuong | Oct 21 | ✅ | `ms.load_templates(4, cfg)`, `ms.backfit(inst, templates, cfg)`; epochs via `pio.epochs_to_mne` |
| Rest and Simon metrics: duration, occurrence, coverage, GEV, transitions | Chuong | Oct 23 | ✅ | **Primary now K = 5**: `results/features/microstates_task.csv` (149 subjects; Rest 148, Simon 140 recordings after K = 5 QC); sensitivity files `_nosmooth`, `_qcincluded`, `_K4` (with 'Central'), `microstates_rest_noblink.csv` (113) |

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

### 5. Recovery pilot (Chuong; reassigned from Minseo and Bohyeon, Oct 23; exploratory)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Define the recovery measure and baseline window (response-locked epochs) | Chuong | Oct 24 | ✅ | Baseline cue −200–0 ms; Euclidean distance of state probabilities from baseline; mean over +200–1000 ms after keypress (`docs/recovery_pilot.md`) |
| Incongruent vs. congruent and error vs. correct recovery (per `docs/TEAM_TASKS.md` C4) | Chuong | Oct 24 | ✅ | 140 subjects; `results/features/recovery.csv`, `results/figures/chuong_recovery.png`. Errors compared with size-matched correct subsets |
| PCA trajectories of response-locked activity (optional) | Chuong | Oct 24 | ⬜ | |
| Split-half reliability (odd vs. even trials, Spearman–Brown) | Chuong | Oct 24 | ✅ | Overall recovery 0.97 (keep); conflict effect 0.24 (drop); error effect 0.68 (keep, marginal) |

### 6. Robustness checks (owner of each analysis)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Rest results without recordings with `nICremoved == 0` | Chuong | Oct 24 | ✅ | Blink sensitivity: `results/features/microstates_rest_noblink.csv` (113 Rest recordings; 35 of 149 have `nICremoved == 0`). Rerunning the Rest analyses on it is B7 (Bohyeon) |
| Key results at K = 4 and K = 5 | Each analysis owner | Oct 24 | ⬜ | |

### 7. Presentation (all)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Slide outline: research questions → data → methods → results → plan | All | Oct 24 | ⬜ | |
| Results slides, one or two per analysis | Each analysis owner | Oct 25 | ⬜ | |
| Slide on sanity checks (reference values reproduced) | Minseo | Oct 25 | ⬜ | |
| Slide on next phases | Chuong | Oct 25 | ⬜ | |
| Rehearsal and timing check | All | Oct 25 | ⬜ | |

### 8. Presentation support (Minseo, from Oct 24; see `docs/TEAM_TASKS.md`)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| M1: install the environment and load one subject | Minseo | Oct 24 | ⬜ | `requirements.txt` |
| M2: sanity-check slide (same as section 7 row) | Minseo | Oct 25 | ⬜ | From `docs/behavior_sanity_checks.md` |
| M3: fill the Phase 1 results summary table below | Minseo | Oct 25 | ⬜ | |
| M4: assemble the slide deck (style, A–E terminology, numbers match the files) | Minseo | Oct 25 | ⬜ | |
| M5: draft questions for the professor (below) | Minseo | Oct 25 | ⬜ | |

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
| Oct 22 | Microstate band-pass changed to 2-20 Hz to reduce muscle contamination; decided after draft templates, before any group or clinical results | Chuong |
| Oct 22 | GFP-peak outlier cut (per recording, over peak GFP values) and minimum peak distance (1 sample, pycrostates default) written explicitly into config.yaml; no change in behaviour | Chuong |
| Oct 22 | Recording QC: exclude a recording from microstate analyses if its template GEV robust z (median/1.4826·MAD, per task) < −3; fit full cohort → flag → refit once → stop; log subject_id, group, MoCA of each exclusion; sensitivity analysis with them included | Chuong |
| Oct 22 | Epochs back-fitted to the templates (Simon conditions and recovery pilot) filtered 2–20 Hz like the templates; 30 Hz low-pass kept for ERP plots only | Chuong |
| Oct 22 | Commit `e2fff64` (Milestone 2, already pushed) also contains Chuong's README edit ("~146 s (median)"); left as is to avoid rewriting shared history | Chuong |
| Oct 22 | Back-fitting uses temporal smoothing (pycrostates documented settings), to match the proposal. Decided after draft metrics showed short durations, but values taken from documentation, not tuned to results. No-smoothing kept as a sensitivity analysis | Chuong |
| Oct 23 | Primary K = 4 by the config k_selection rule on full-cohort group fits (silhouette and Calinski-Harabasz both best at K = 4; Davies-Bouldin ~1e15 and Dunn ~0 at every K, i.e. uninformative). `primary_k: 4` added to config.yaml | Chuong |
| Oct 23 | Recording QC (K = 4, pass 1) excluded: sub-P001 Rest (PD, MoCA 19, z = −3.72), sub-P040 Rest (PD, MoCA 27, z = −3.09), sub-P048 Simon (PD, MoCA 21, z = −4.84). At K = 5: sub-P001 Rest, sub-P048 Simon. Subjects with a flagged recording left out of the template refit; flagged recordings empty in the metric files (kept in `_qcincluded`) | Chuong |
| Oct 23 | State labels assigned automatically by polarity-invariant correlation with the approved draft A–D (Hungarian matching; extra maps E, F, G by explained variance); pending visual confirmation | Chuong |
| Oct 23 | Primary K set to 5 instead of the rule's K = 4: two of the four selection criteria (Davies-Bouldin, Dunn) were degenerate at all K, so the rule did not run as intended; and only K = 5 recovers canonical A-D (r >= 0.97 with the draft), with E separate. GEV gains also flatten after K = 5 (+2.0, +1.6, +0.9). Decided before any group or clinical analyses. K = 4 kept as a sensitivity analysis | Chuong |
| Oct 23 | Davies-Bouldin and Dunn removed from the K-selection rule (`k_selection_scores: [silhouette, calinski_harabasz]`); note the two remaining scores also rank K = 4 first, so K = 5 rests on the decision above | Chuong |
| Oct 23 | Labels visually confirmed: K = 5 A, B, C, D (canonical) and E (central maximum, opposite-sign temporal edges); K = 4 fourth map renamed 'Central' (C at K = 4 includes canonical D); K = 6, 7 automatic labels, sensitivity only. In config.yaml (`confirmed_labels`) | Chuong |
| Oct 23 | Primary files use the K = 5 recording QC: excluded sub-P001 Rest (PD, MoCA 19, z = −3.47) and sub-P048 Simon (PD, MoCA 21, z = −4.70); `_K4` file uses the K = 4 QC (adds sub-P040 Rest) | Chuong |
| Oct 23 | Milestone 4 choices written into config.yaml: QC per K; subjects with any flagged recording left out of the template refit (unflagged recordings still back-fitted); QC GEV on all of a recording's GFP peaks; labels by one-to-one matching, extra maps ordered by explained variance | Chuong |
| Oct 23 | Minseo unavailable until Oct 24: behavior (section 1) and recovery pilot moved to Chuong; PCA trajectories moved from Bohyeon to Chuong (optional); Minseo takes presentation-support tasks M1-M5 | Chuong |
| | | |
| Oct 23 | Simon error rates in `behavior.csv` follow the TEAM_TASKS spec (errors / responded trials); the Dataset Guide reference values use errors / all trials. Both reproduced; team to confirm which definition the RQ3 analyses use | Chuong |
| Oct 23 | Teammates use GitHub read-only (browse/download); outputs sent to Chuong, who integrates them and maintains PROGRESS.md | Chuong |
| Oct 23 | Template files (`results/templates/*.fif`, derived maps, no EEG) committed so teammates can download them; labels, peaks and QC tables stay git-ignored | Chuong |
| Oct 23 | Recovery pilot choices where the spec was open: Euclidean distance between state-probability vectors; probabilities over labeled trials (epoch-edge segments unlabeled); errors compared with random correct subsets of the same size (fewer trials inflate the distance); error measures need ≥ 10 errors; correct and error trials both need RT 150–1500 ms; split halves by odd/even trial number, each with its own baseline | Chuong |
| Oct 23 | Recovery pilot go/no-go (config: keep if reliability ≥ 0.5): keep overall recovery (0.97) and error effect (0.68, marginal); drop conflict effect (0.24) | Chuong |

## Blockers

| Date raised | Blocker | Who is blocked | Owner | Resolved? |
|---|---|---|---|---|
| Oct 21 | Simon epoch files don't carry a readable trial table (MATLAB `table`, not decodable by h5py). Checked all 147 subjects: `lock-cue` = one epoch per `trials.tsv` row; `lock-resp` = one epoch per row with non-missing `onset_response`. `beh_accuracy == 99` rows have a timeout `onset_response` (`rt_eeg` ≈ 1.0 s), so drop them in response-locked analyses. Use `src.io.match_trials`. Open point: epoch order assumed equal to TSV order (counts match; order not verifiable without MATLAB) | Minseo (recovery), Seonwoo | Chuong | Resolved with rule; order check open |
| Oct 21 | sub-P048 (draft set) looks artifact-dominated: individual maps peak at AF8/Fp2, FC5/T7, TP10/I1; template GEV 16% vs 44% mean. No config rule covers this | Chuong (final templates) | Chuong | Yes. Draft: P048 replaced by P009. Data-quality wording resolved Oct 22: template-GEV robust-z recording QC (config `recording_qc`); it excludes sub-P048 Simon at every K |
| Oct 21 | Draft template states not yet matched to canonical A–D (config: visual inspection, fixed once). Needed before `microstates_task.csv` column names | Chuong, all | Chuong | Provisional labels on draft (Oct 22); final labels on final templates |
| Oct 22 | Task-level sanity check fails on durations: mean durations 44–57 ms (only 8% of subject×state values in 60–120 ms; median segment 36 ms, 22% at the 20 ms floor) with draft templates, 2–20 Hz, `min_segment_length` 5, no temporal smoothing. Proposal says back-fit "with temporal smoothing", config sets none. Bookkeeping verified (occ × dur = cov) | Chuong (task metrics), Bohyeon | Chuong | Yes (Oct 22): temporal smoothing from pycrostates 0.6.1 tutorial (factor 10, half window 10); no-smoothing kept as sensitivity |
| Oct 22 | Back-fitting is slow: pycrostates short-segment rejection ~10 min per Simon recording (0.1 s without it); 15 subjects took 31 min with 8 parallel workers; full cohort ≈ 2–3 h per run | Chuong | Chuong | Yes (Oct 22): `ms.reject_short_segments` is the default; identical labels to pycrostates on all 15 subjects × Rest/Simon × smooth/no-smooth (`tests/test_backfit.py`, 64 passed); ~1.7 s per subject |
| Oct 23 | Full-cohort K = 4 (primary by config rule) has no map matching draft D (best r = 0.42): its 4th map is a central maximum with temporal-edge negativity (= K = 5 map E, 3.7–4.4% of GEV), while fronto-central D is merged into C (r = 0.95 with draft C, 0.86 with draft D). At K = 5, A–D match the draft at r ≥ 0.97. Column `*_D` at K = 4 is therefore not canonical D | Chuong, Bohyeon, Seonwoo (anyone using K = 4 labels) | Chuong | Yes (Oct 23): primary K = 5; K = 4 fourth map renamed 'Central' (sensitivity) |
| Oct 23 | Davies-Bouldin (~1e15) and Dunn (~0) are degenerate for every K with pycrostates' distance (\|1/r\| − 1 explodes for near-orthogonal maps); the K choice rests on silhouette and Calinski-Harabasz (both K = 4) | Chuong | Chuong | Yes (Oct 23): both removed from the K-selection rule in config.yaml |

## Daily log

| Date | Who | What was done | Next |
|---|---|---|---|
| Oct 20 | | | |
| Oct 21 | Chuong | Milestone 1: loading/preprocessing in `src/`, exclusion logged (Rest 149→149; Simon 149→147→141 after `trials_simon`; nInterp > 6 removed none; 35/149 Rest flagged `nICremoved == 0`); env + `requirements.txt` | GFP peaks and draft templates |
| Oct 21 | Chuong | Milestone 2: GFP peaks + two-level modified k-means; draft K = 4 templates (15 subjects, 10 PD / 5 C) saved to `results/templates/templates_K4.fif`; group GEV 68.1% on individual maps, 44.1% ± 12.5% on subjects' GFP peaks; `backfit` tested on Raw and Epochs; `match_trials`, `epochs_to_mne` helpers | A–D labels, P048 decision, then task-level metrics |
| Oct 22 | Chuong | Band-pass 2–20 Hz in config; draft K = 4 refit without sub-P048 (replaced by sub-P009, seed 42): group GEV 75.5%, template GEV 53.7% ± 9.0%; provisional A–D labels; tested proposed GEV data-quality rule (does not flag P048 on individual-fit GEV) | Team check of maps; QC rule wording; Milestone 3 |
| Oct 22 | Chuong | Milestone 3 (draft): back-fit draft K = 4 templates to Rest and Simon of the 15 draft subjects; `results/features/microstates_task.csv` (15 × 63); labels saved; coverage sums 100%, transition rows sum to 1, gev_total Rest 48.0% ± 10.6, Simon 48.2% ± 9.4 (vs 53.7% on GFP peaks); durations short (44–57 ms), see Blockers | Team decision on smoothing; final templates |
| Oct 22 | Chuong | Milestone 3 rerun with smoothing (factor 10, half window 10, min segment 20 ms): mean durations by state Rest 114–170 ms, Simon 115–171 ms; median segment 100 ms, 0% at the 20 ms floor; gev_total Rest 42.9% ± 10.0, Simon 43.1% ± 9.3; coverage and transition sums OK. No-smoothing sensitivity in `microstates_task_nosmooth.csv` (rho with smoothed 0.93–1.00). Fast exact back-fit + tests. Note: state C dominates sub-C049 (80% coverage, Rest and Simon) and sub-C004 Simon (58% vs 13% Rest); recheck on full cohort | Milestone 4: final templates on full cohort |
| Oct 23 | Chuong | Milestone 4: full-cohort templates K = 4–7 (two-level, QC two-pass, auto-labels vs draft), primary K = 4 by rule; task-specific templates; metrics for 149 subjects (primary smoothed, no-smoothing, K = 5, QC-included, Rest no-blink: 35/149 Rest recordings have nICremoved == 0). Sanity checks pass (coverage 100%, transitions sum to 1; median segment 100 ms, 0% at floor). sub-C049 not flagged by QC (z 0.07 Rest, 0.78 Simon); its state C coverage is high (97th/99th percentile) but 4 subjects are as high or higher (max 92%, sub-P051) | Morning: visual check of labels; K = 4 vs 5 decision |
| Oct 23 | Chuong | Switched primary to K = 5 and rebuilt primary files (K = 5 QC; sanity checks pass: coverage 100%, transitions sum to 1, median segment 96 ms, 0% at floor, gev_total Rest 46.3%, Simon 43.4%); confirmed labels applied (K = 4 'Central'); K = 5 task-specific templates (all \|r\| ≥ .98). State E descriptives: coverage median Rest 11.0%, Simon 12.1%; in Simon epochs (20 subjects) P(E) rises from .11 to a peak of .155 at 452 ms post-stimulus and .148 at +16 ms around the keypress, then drops; \|r\| = .09 with the first-draft muscle-like map. sub-C049 state C coverage at K = 5: robust z 3.8 Rest, 4.7 Simon (QC not flagged). Cleared outputs from notebooks/test.ipynb | Hand over K = 5 features to Bohyeon/Seonwoo |
| Oct 23 | Chuong | Task board `docs/TEAM_TASKS.md` added (paths checked) and tasks reassigned (Minseo away until Oct 24). C1–C3 done: Simon RT 514.4 / 557.4 ms (reference 514.7 / 557.4), error rates 6.44% / 14.82% with the Guide's definition (errors / all trials; 6.84% / 16.27% with the spec's errors / responded trials); by group RT effect PD +47 / C +35 ms, incongruent errors PD 18.0% / C 8.5%; `results/features/behavior.csv` (141 × 9); summary for Minseo in `docs/behavior_sanity_checks.md` | C4 recovery pilot after approval |
| Oct 23 | Chuong | GitHub read-only workflow (`docs/SETUP.md`, templates committed). C4 recovery pilot on 140 subjects: split-half reliability (Spearman–Brown) overall recovery 0.97 (rank-based 0.81), conflict effect 0.24, error effect 0.68 (rank-based 0.57) → keep overall and error measures, drop conflict; robust to ending the window at +800 ms; `recovery.csv`, `chuong_recovery.png`, `docs/recovery_pilot.md` | C5 (optional PCA), C6 slides |
