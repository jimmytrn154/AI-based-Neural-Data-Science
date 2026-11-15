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
| 1. Data structure (unsupervised) | Oct 20 → Oct 26 | ✅ | Progress presentation held (per team; the deck is not in the repo) |
| 2. Consolidation | Oct 27 → about two weeks before the final | 🟨 | Phase 2 results presented Nov 9; plan revised Nov 15 after the professor's feedback (see Phase 2) |
| 3. Prediction | After P2.2 and P2.4 | ⬜ | Planned (see Phase 3) |
| 4. Extensions | Late semester | ⬜ | Placeholder |
| Final presentation | Week 16 (date: [TO FILL]) | ⬜ | |

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
| Oct 24 | All feature tables merged; PD/MoCA associations run; figures chosen | All | ✅ (teammates' outputs committed in their own folders: `seonwoo_S1-S6_review/`, `bohyeon_adjusted/`; not moved into `results/`) |
| Oct 25 | Slides finished and rehearsed | All | ✅ (per team; the deck is not in the repo) |
| Oct 26 | **Progress presentation** | All | ✅ (per team) |

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
| Behavioral feature table (`results/features/behavior.csv`) | Chuong | Oct 23 | ✅ | 141 rows, 14 columns: primary error rates (errors / responded trials) plus Guide-definition `_all` columns and responded counts (Oct 24); rules in config `behavior` |

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
| Load Simon cue-locked epochs; drop bad epochs; keep correct trials | Seonwoo | Oct 21 | ✅ | Config trial rules (correct, RT 150–1500 ms, no bad epochs, ≥ 40 per condition): 138 subjects. Committed in `seonwoo_S1-S6_review/` (reviewed Nov 2; not moved into `results/`). |
| Reproduce the Simon ERP at Cz (Dataset Guide, Figure 2) | Seonwoo | Oct 21 | ✅ | `seonwoo_erp_cz.png`: central positivity ≈ 440 ms, larger for incongruent; matches the Guide (138 subjects) |
| Topographic segmentation of each condition's ERP (0–800 ms) | Seonwoo | Oct 22 | ✅ | Template-based: condition ERPs labeled with the shared K = 5 templates (137 subjects); `seonwoo_erp_segmentation.png` |
| Single-trial back-fitting and state-probability time courses | Seonwoo | Oct 23 | ✅ | 137 subjects (Simon recording QC applied); 2–20 Hz, 250 Hz, K = 5, documented smoothing |
| Congruent vs. incongruent: cluster-based permutation test | Seonwoo | Oct 24 | ✅ | 5,000 permutations, family-wise corrected over states and times: E incongruent > congruent 464–576 ms (+1.6 points, p = 0.008); C incongruent < congruent 348–512 ms (−1.2 points, p = 0.006); `seonwoo_conflict_timecourse.png` |
| Conflict-effect feature in 200–350 ms and 400–750 ms windows | Seonwoo | Oct 24 | ✅ | `conflict.csv` (137 rows; per-state incongruent − congruent differences, 10 columns), reproduced exactly from S4. Definition of the scalar conflict effect still to agree (see Blockers) |

### 4. Comparisons, networks and clustering (Bohyeon)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Statistics functions on dummy data (mixed models, partial Spearman, FDR) | Bohyeon | Oct 21 | ✅ | `src/stats.py` (committed in `bohyeon/` and `bohyeon_adjusted/`; the Nov 7 version adds a fallback optimizer for the mixed model; not yet moved into `src/`): FDR by family, paired Wilcoxon + rank-biserial, mixed model, partial Spearman with bootstrap CI, group model |
| Rest vs. Simon within-subject comparisons | Bohyeon | Oct 23 | ✅ | `rq1_stats.csv` (139 subjects): Simon > Rest for A and E coverage and occurrence; Simon < Rest for C coverage and C duration (FDR-significant). Mixed models added Nov 7 (`rq1_mixed.csv`, reproduced): no task × group interaction; group term shows less E at rest in PD (secondary, see Blockers) |
| Reorganization index (Jensen–Shannon divergence) and coverage change | Bohyeon | Oct 23 | ✅ | `reorganization.csv` (139 subjects), reproduced: JSD = mean over the five source states of the Rest-vs-Simon Jensen–Shannon divergence (log 2) of transition probabilities |
| Transition networks (rest vs. Simon) | Bohyeon | Oct 24 | ✅ | `bohyeon_transition_networks.png`: 2 strongest outgoing edges per state (same at Rest and in Simon); rule to be stated on the figure |
| Subject embedding (PCA, UMAP) and k-means clustering; alignment with group and MoCA | Bohyeon | Oct 24 | ✅ | Rerun Nov 7 with state E: k = 2 (searched 2–10), silhouette 0.21 (weak); no alignment with group (ARI 0.001) or MoCA |
| PD vs. control and MoCA associations for all primary features | Bohyeon | Oct 24 | ✅ | Rerun Nov 7 (`rq3_stats.csv`): one RQ3 family of 108 tests, mean durations and Mann–Whitney added, sub-P001 restored; nothing survives FDR (smallest adjusted p = 0.10). Small follow-up: keep all valid Rest recordings in the `rest_meandur_*` tests |

### 5. Recovery pilot (Chuong; reassigned from Minseo and Bohyeon, Oct 23; exploratory)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Define the recovery measure and baseline window (response-locked epochs) | Chuong | Oct 24 | ✅ | Baseline cue −200–0 ms; Euclidean distance of state probabilities from baseline; mean over +200–1000 ms after keypress (`docs/recovery_pilot.md`) |
| Incongruent vs. congruent and error vs. correct recovery (per `docs/TEAM_TASKS.md` C4) | Chuong | Oct 24 | ✅ | 140 subjects; `results/features/recovery.csv`, `results/figures/chuong_recovery.png`. Errors compared with size-matched correct subsets |
| PCA trajectories of response-locked activity (optional) | Chuong | Phase 2 | ⏸️ | Postponed to Phase 2 (Oct 24) |
| Split-half reliability (odd vs. even trials, Spearman–Brown) | Chuong | Oct 24 | ✅ | Overall recovery 0.97 (keep); conflict effect 0.24 (drop); error effect 0.68 (keep, marginal) |

### 6. Robustness checks (owner of each analysis)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Rest results without recordings with `nICremoved == 0` | Chuong | Oct 24 | ✅ | Blink-sensitivity file `results/features/microstates_rest_noblink.csv` (113 Rest recordings; 35 of 149 have `nICremoved == 0`). RQ1 rerun on it by Bohyeon (`bohyeon 261112/results/rest_noblink/`, N = 106): E and A changes hold, C's decrease no longer FDR-significant. Remaining reruns are P2.7 |
| Key results at K = 4 and K = 5 | Each analysis owner | Oct 24 | 🟨 | RQ1 (B2–B5) rerun at K = 4 by Bohyeon (`bohyeon 261112/results/K4/`); RQ2 and RQ3 not rerun at K = 4 |

### 7. Presentation (all)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| Slide outline: research questions → data → methods → results → plan | All | Oct 24 | ✅ | Presented (per team; deck not in repo) |
| Results slides, one or two per analysis | Each analysis owner | Oct 25 | ✅ | Presented (per team; deck not in repo) |
| Slide on sanity checks (reference values reproduced) | Minseo | Oct 25 | ✅ | Presented (per team); source `docs/behavior_sanity_checks.md` |
| Slide on next phases | Chuong | Oct 25 | ✅ | Draft `docs/slides_chuong.md`; presented (per team) |
| Rehearsal and timing check | All | Oct 25 | ✅ | Per team |

### 8. Presentation support (Minseo, from Oct 24; see `docs/TEAM_TASKS.md`)

| Task | Owner | Due | Status | Notes |
|---|---|---|---|---|
| M1: install the environment and load one subject | Minseo | Oct 24 | ⬜ | Not confirmed in the repo |
| M2: sanity-check slide (same as section 7 row) | Minseo | Oct 25 | ⬜ | From `docs/behavior_sanity_checks.md` |
| M3: fill the Phase 1 results summary table below | Minseo | Oct 25 | ⬜ | Table below is still empty (checked Nov 15) |
| M4: assemble the slide deck (style, A–E terminology, numbers match the files) | Minseo | Oct 25 | ⬜ | |
| M5: draft questions for the professor (below) | Minseo | Oct 25 | ⬜ | Section below is still empty (checked Nov 15) |

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

## Phase 2: consolidation

**Goal:** test whether the Phase 1 picture (Simon shifts the brain-state balance from C towards E, and conflict trials do so more) is real, specific to E, and robust, then interpret it. Phase 2 scope narrowed to validating and explaining the C-down / E-up pattern (rest -> Simon and congruent -> incongruent); the main research question is unchanged.

**Research question structure:**

- **Main question** (answered in Phase 3): Do changes in brain-state dynamics from rest to task carry information about cognitive function in PD, beyond resting EEG and behavioral tests?
- **Sub-questions RQ1–RQ3** (done in Phase 1): rest vs. Simon dynamics; congruent vs. incongruent; groups and cognition.
- **Current focus** (Phase 2): is the C-down / E-up pattern real and robust, what is E, what does it mean, and is it new?

### Professor's feedback on the Phase 2 presentation (Nov 9)

1. **E vs. normative templates.** Compare state E with published normative microstate templates, with proper controls and baselines.
2. **Probability dependence.** State probabilities sum to 1: is the E increase more than the compositional consequence of the C decrease?
3. **C↔E transitions.** Look at the transitions between C and E directly, in the conflict windows.
4. **Randomization tests.** Use permutation or surrogate tests rather than relying only on parametric p-values.
5. **Brain activity.** What does the C → E pattern say about brain activity?
6. **Novelty.** What is new compared with other papers?

### Done so far (Nov 9 deck)

| Item | Owner | Status | Evidence |
|---|---|---|---|
| Phase 2 scope narrowed to validating and explaining the C-down / E-up pattern (rest -> Simon and congruent -> incongruent); the main research question is unchanged. | All | ✅ | Nov 9 deck (not in repo) |
| E internal validation: spatial distinctness (K5 E–C \|r\| = 0.337; C–D 0.796 for comparison), recurrence across K (K4 Central–K5 E 0.986, K5 E–K6 F 0.960, K6 F–K7 G 0.948) and across tasks (Rest-only E vs. Simon-only E 0.948; vs. shared E 0.984 / 0.983) | Chuong (templates) | ✅ | Recomputed Nov 15 from `results/templates/` |
| Coverage sensitivity: Simon − Rest coverage A +2.20, C −3.00, E +1.80 points (N = 139); E +1.45 without smoothing, +1.76 with QC-excluded recordings (N = 141); K = 4 and no-blink runs also done | Bohyeon | ✅ | `bohyeon 261112/results/{task,nosmooth,task_qcincluded,K4,rest_noblink}/rq1_stats_*.csv` |
| **P2-B7** conflict state contrast (N = 137): D = Δ(E) − Δ(C) is +1.23 points in theta (FDR p = 0.00724), +0.67 in N2 (0.214); PD − control in D: N2 −2.80 points (FDR p = 0.0291, **not age/sex adjusted**), theta −0.86 (0.349) | Bohyeon | ✅ | `bohyeon 261112/results/task/B7_conflict_state_contrast.csv` |
| **P2-B8** conflict group windows (N = 136, age/sex adjusted): PD − control N2 Δ(C) +2.46 points (FDR p = 0.00513), D −2.90 (0.0413); theta Δ(E) −1.35 (0.0443); group × window for C p = 0.00203 | Bohyeon | ✅ | `bohyeon 261112/results/task/b8_adjusted_tests_task.csv` |
| Normative-template comparison: best matches Zanesco et al. 2020 (\|r\| = 0.868) and 2021a (0.843) | TBD | 🟨 | Results reported by a teammate; code not yet in repo |

**Limitations of the work so far:**
- E partly splits across two maps at K = 7 (\|r\| = 0.82 with F, 0.83 with G); at K = 4–6 it has one clear counterpart.
- P2-B7's group difference in D is **unadjusted** (no age/sex); P2-B8 holds the adjusted version (N2 D −2.90 points, FDR p = 0.041).

All numbers above were checked against the outputs on Nov 15 (`bohyeon 261112/README.md`). The P2-B7 and P2-B8 samples differ by sub-P069 (sex recorded "F\|M", dropped from the sex-adjusted P2-B8 models).

### Revised plan (Nov 15)

| ID | Item | Feedback | Owner | Status | What remains |
|---|---|---|---|---|---|
| P2.1 | E vs. normative templates | 1 | Minseo (pending confirmation) | 🟨 | Full match matrix (all 5 of our maps × all reference classes); our A–D as a control; baseline of correlations between non-matching reference classes (our own C–D is 0.80, so \|r\| ≈ 0.85 needs context); results for **all** reference studies, not only the best matches; re-reference to the common electrode subset, or interpolate to our 64 positions; add the Microstate Toolbox normative templates (Custo 2017; Koenig 2024 meta-maps); put the code in the repo |
| P2.2 | Probability dependence: is E-up more than the compositional consequence of C-down? | 2 | Chuong | ⬜ | All five state changes; proportional-redistribution null; E's share of non-C time, P(E) / (1 − P(C)); log-ratio analysis; a non-competitive measure (per-template spatial correlation over time); second-best labels of frames that switch to E |
| P2.3 | C↔E transitions | 3 | Seonwoo (pending confirmation) | ⬜ | C → E and E → C in the N2 and theta windows, incongruent vs. congruent; observed vs. expected given occurrences; directional asymmetry; where E segments come from |
| P2.4 | Randomization / surrogate tests | 4 | Bohyeon | ⬜ | Group-label permutation for the P2-B7 / P2-B8 group results, Freedman–Lane for the age/sex-adjusted models and a max-statistic correction across the six P2-B8 group tests; within-subject sign-flip for condition contrasts; segment-shuffle or Markov surrogates for transitions; all on **one fixed sample** |
| P2.5 | Brain-activity interpretation | 5 | Minseo, with Seonwoo (pending confirmation) | ⬜ | Trial level: correct vs. error incongruent trials, and reaction time; subject level: E increase vs. Simon error rate; then a cautious interpretation using P2.1–P2.3 |
| P2.6 | Novelty: what is new? | 6 | Minseo (pending confirmation) | ⬜ | Targeted literature check: microstates in PD (rest and task); microstates in conflict tasks (Simon, Stroop, flanker); rest-to-task microstate changes; extra classes in clinical groups. Then a "contribution" statement, with claims at the level the evidence supports |
| P2.7 | Blink sensitivity | — | Chuong | 🟨 in progress: run done by Bohyeon, effect-size comparison pending review | RQ1 rerun on `microstates_rest_noblink.csv` done by Bohyeon; rerun the other key results (conflict contrasts, group results). Coverage comparison in the P2.7 note below |
| P2.8 | Alpha control: spectral parameterization of Rest | — | Seonwoo (pending confirmation) | ⬜ | Periodic (alpha peak) and aperiodic Rest features; test whether the main results survive them |
| P2.9 | Recovery pilot: report the go/no-go result | — | Chuong | 🟨 | Result exists (`docs/recovery_pilot.md`: keep overall recovery and error effect, drop the conflict effect); it was not in the Nov 9 deck, so add it to the final |

**P2.7 note (Nov 15, from Bohyeon's outputs; report only).** Simon − Rest coverage change: mean [95% bootstrap CI], rank-biserial r, Cohen's d_z, FDR p (RQ1 family of 20 tests):

| State | Primary run, N = 139 | No-blink run, N = 106 | Primary data, same 106 subjects |
|---|---|---|---|
| A | +2.20 [+0.75, +3.67], r 0.29, d_z 0.26, p 0.012 | +2.12 [+0.50, +3.71], r 0.29, d_z 0.25, p 0.038 | identical to the no-blink run |
| C | −3.00 [−5.20, −0.58], r −0.26, d_z −0.21, p 0.025 | −2.41 [−5.03, +0.36], r −0.22, d_z −0.17, p 0.128 | identical to the no-blink run |
| E | +1.80 [+0.72, +2.92], r 0.40, d_z 0.27, p 0.0005 | +1.58 [+0.27, +2.78], r 0.34, d_z 0.24, p 0.019 | identical to the no-blink run |

The no-blink file holds the same per-subject values as the primary file, so the no-blink run equals the primary analysis restricted to those 106 subjects. In the 33 subjects it drops (Rest recordings with `nICremoved == 0`), mean changes are A +2.46, C −4.88, E +2.52 points. Recomputed with Bohyeon's `stats.py` and procedure, which reproduces her primary results exactly.

**Dependencies:**
- P2.5 depends on P2.1–P2.3.
- Novelty claims (P2.6) only for findings that survive P2.4.
- Phase 3 starts after P2.2 and P2.4.

**Postponed (kept listed):** RT-sorted analysis of E; conflict adaptation; PCA trajectories (recovery pilot, C5); Phase 4 extensions.

**Timeline:** today is Nov 15; final presentation in Week 16 (date: [TO FILL]). Target: Phase 2 done about two weeks before the final, then Phase 3.

### Earlier Phase 2 placeholder items (Oct 26), for the record

- [x] Research questions after the feedback: Phase 2 scope narrowed to validating and explaining the C-down / E-up pattern (rest -> Simon and congruent -> incongruent); the main research question is unchanged.
- [ ] Robustness: K, blinks, filter choices (partly done: RQ1 at K = 4, no smoothing, QC included, no-blink; rest is P2.7)
- [ ] Alpha control with spectral parameterization (Rest) → P2.8
- [ ] Fix the final feature set
- [x] Go/no-go decision on the recovery analysis (Oct 23; reporting is P2.9)

## Phase 3: prediction (planned; starts after P2.2 and P2.4)

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
| Oct 24 | Primary Simon error rate = errors / responded trials (pre-specified; separates commission errors from omissions); Guide definition kept as extra columns. Decided before any RQ3 analysis | Chuong |
| Oct 24 | Recovery files split by reliability: `recovery.csv` keeps the measures that passed (overall recovery, error effect); the failed conflict-recovery columns moved to `recovery_conflict_unreliable.csv` (do not use) | Chuong |
| Oct 24 | Epochs are resampled to 250 Hz (after filtering) before back-fitting; written into config.yaml (`simon_conditions.resample_before_backfit_hz`) | Chuong |
| Oct 24 | Behavior rules written into config.yaml (`behavior`): `simon_n_*` count all trials; no RT limit for behavior features | Chuong |
| Oct 24 | C5 (PCA trajectories of the recovery pilot) postponed to Phase 2 | Chuong |
| Nov 15 | Phase 2 revised after professor feedback (Nov 9): four validity checks plus interpretation and novelty questions; postponed items listed | Chuong |
| Nov 15 | The deck's new Phase 2 analyses "B7" and "B8" are called P2-B7 (conflict state contrast) and P2-B8 (conflict group windows) in all docs, to avoid a clash with Bohyeon's Phase 1 task IDs B7 and B8; script names unchanged | Chuong |
| Nov 15 | Main research question unchanged; Phase 2 scope narrowed to the C-down / E-up pattern. Final deck uses the original wording, including 'beyond resting EEG and behavioral tests' | Chuong |
| Nov 15 | All Phase 2 group analyses (P2.4 onward) use the sex-complete sample (136 people for the Simon conflict features; sub-P069 excluded because its sex is recorded as "F\|M"); in config.yaml (`phase2.group_sample`) | Chuong |
| Nov 15 | Phase 2 owners proposed, pending confirmation: P2.1 Minseo, P2.3 Seonwoo, P2.5 Minseo with Seonwoo, P2.6 Minseo, P2.8 Seonwoo | Chuong |

## Blockers

| Date raised | Blocker | Who is blocked | Owner | Resolved? |
|---|---|---|---|---|
| Oct 21 | Simon epoch files don't carry a readable trial table (MATLAB `table`, not decodable by h5py). Checked all 147 subjects: `lock-cue` = one epoch per `trials.tsv` row; `lock-resp` = one epoch per row with non-missing `onset_response`. `beh_accuracy == 99` rows have a timeout `onset_response` (`rt_eeg` ≈ 1.0 s), so drop them in response-locked analyses. Use `src.io.match_trials`. Open point: epoch order assumed equal to TSV order (counts match; order not verifiable without MATLAB) | Minseo (recovery), Seonwoo | Chuong | Yes (Oct 24): order verified for all 147 subjects. Every cue (54,546) and response (53,315) epoch rebuilt from continuous data + onsets (Guide 4.1) matches the distributed epoch in the assumed order: max difference 2.4e-4 µV, 0 mismatching trials; a one-trial shift makes all trials mismatch (`notebooks/chuong_epoch_order_check.py`) |
| Oct 21 | sub-P048 (draft set) looks artifact-dominated: individual maps peak at AF8/Fp2, FC5/T7, TP10/I1; template GEV 16% vs 44% mean. No config rule covers this | Chuong (final templates) | Chuong | Yes. Draft: P048 replaced by P009. Data-quality wording resolved Oct 22: template-GEV robust-z recording QC (config `recording_qc`); it excludes sub-P048 Simon at every K |
| Oct 21 | Draft template states not yet matched to canonical A–D (config: visual inspection, fixed once). Needed before `microstates_task.csv` column names | Chuong, all | Chuong | Provisional labels on draft (Oct 22); final labels on final templates |
| Oct 22 | Task-level sanity check fails on durations: mean durations 44–57 ms (only 8% of subject×state values in 60–120 ms; median segment 36 ms, 22% at the 20 ms floor) with draft templates, 2–20 Hz, `min_segment_length` 5, no temporal smoothing. Proposal says back-fit "with temporal smoothing", config sets none. Bookkeeping verified (occ × dur = cov) | Chuong (task metrics), Bohyeon | Chuong | Yes (Oct 22): temporal smoothing from pycrostates 0.6.1 tutorial (factor 10, half window 10); no-smoothing kept as sensitivity |
| Oct 22 | Back-fitting is slow: pycrostates short-segment rejection ~10 min per Simon recording (0.1 s without it); 15 subjects took 31 min with 8 parallel workers; full cohort ≈ 2–3 h per run | Chuong | Chuong | Yes (Oct 22): `ms.reject_short_segments` is the default; identical labels to pycrostates on all 15 subjects × Rest/Simon × smooth/no-smooth (`tests/test_backfit.py`, 64 passed); ~1.7 s per subject |
| Oct 23 | Full-cohort K = 4 (primary by config rule) has no map matching draft D (best r = 0.42): its 4th map is a central maximum with temporal-edge negativity (= K = 5 map E, 3.7–4.4% of GEV), while fronto-central D is merged into C (r = 0.95 with draft C, 0.86 with draft D). At K = 5, A–D match the draft at r ≥ 0.97. Column `*_D` at K = 4 is therefore not canonical D | Chuong, Bohyeon, Seonwoo (anyone using K = 4 labels) | Chuong | Yes (Oct 23): primary K = 5; K = 4 fourth map renamed 'Central' (sensitivity) |
| Oct 23 | Davies-Bouldin (~1e15) and Dunn (~0) are degenerate for every K with pycrostates' distance (\|1/r\| − 1 explodes for near-orthogonal maps); the K choice rests on silhouette and Calinski-Harabasz (both K = 4) | Chuong | Chuong | Yes (Oct 23): both removed from the K-selection rule in config.yaml |
| Nov 2 | Bohyeon sent `src/stats.py` and the B2–B6 outputs but not the B2–B6 analysis code (handoff rules ask for code with outputs); key numbers were reproduced from repo files, the rest can't be audited. Open fixes: mixed models (B2), B5 with state E, one RQ3 FDR family, mean durations and Mann–Whitney check, sub-P001 in conflict analyses (B6) | Chuong (integration), team (RQ3 slides) | Bohyeon | Yes (Nov 7): code and all requested fixes delivered in `bohyeon_adjusted/`; small follow-ups listed in `docs/TEAM_TASKS.md` (Rest data in `rest_meandur_*` tests, unsent secondary outputs) |
| Nov 2 | Conflict effect: the config's primary feature is one distance per window, but `conflict.csv` holds per-state differences; Bohyeon used the L2 norm of the window-averaged difference vector. Team to agree the definition and which file holds it (checked: not related to trial counts, ρ ≈ −0.02 to −0.05) | Seonwoo, Bohyeon | Team | Partly (Nov 15): P2-B7 / P2-B8 work with per-state differences (Δ(C), Δ(E), D = Δ(E) − Δ(C)); the scalar distance in RQ3 is still Bohyeon's L2 version and not formally agreed |
| Nov 7 | The RQ1 mixed model's group term shows less state E at rest in PD (coverage −4.0 points, FDR p = 0.012; occurrence −0.34 per s, FDR p = 0.003; robust across optimizers, a Rest-only regression and Mann–Whitney). It is a group comparison on a feature outside the pre-specified RQ3 set, from a secondary model; decide how to present it (secondary or exploratory) | Bohyeon (RQ1/RQ3 slides), all | Team | Superseded (Nov 15): group results are now part of the Phase 2 plan (P2.4 permutation tests; presentation as secondary) |
| Nov 15 | Normative-template comparison (Zanesco et al. 2020, 2021a) was reported in the Nov 9 deck, but its code is not in the repo and its owner is not recorded | P2.1 | Minseo (pending confirmation) | No |

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
| Oct 24 | Chuong | Behavior: primary error rate errors / responded trials plus Guide-definition columns; recovery files split by reliability; config: epoch resampling and behavior rules; C5 postponed to Phase 2. Epoch-order blocker resolved (all 147 subjects, 0 mismatching trials). Slide draft `docs/slides_chuong.md` (C6) | Integrate teammates' outputs; finalize slides with Minseo |
| Nov 2 | Chuong | Reviewed Seonwoo's S1–S6 (`seonwoo_S1-S6_review/`) and Bohyeon's B1–B6 (`bohyeon/`) deliverables. Seonwoo: all config rules followed; ERP matches the Guide; 2 significant clusters (E and C, ≈ 350–580 ms); `conflict.csv` reproduced. Bohyeon: RQ1, JSD, coverage changes and RQ3 numbers reproduced; code for B2–B6, mixed models and several RQ3/B5 fixes still needed. `docs/TEAM_TASKS.md` updated; nothing integrated yet | Integrate outputs; request fixes from Bohyeon; agree conflict-effect definition |
| Nov 7 | Chuong | Reviewed Bohyeon's adjusted B1–B6 (`bohyeon_adjusted/`): code delivered; mixed models (all 60 reproduced with another optimizer), B5 with state E (k = 2, weak), B6 with one RQ3 family (108 tests, none FDR-significant), durations, Mann–Whitney, sub-P001 restored. Secondary finding: less E at rest in PD. B2, B5, B6 → done; small follow-ups listed in `docs/TEAM_TASKS.md` | Integrate teammates' outputs into the repo layout; agree conflict-effect definition and how to present the E finding |
| Nov 15 | Chuong | Verified Bohyeon's Phase 2 folder (`bohyeon 261112/`, README added): every checked slide number found in the outputs; E-validation correlations recomputed from the templates; P2-B7 (137) vs P2-B8 (136) differ by sub-P069 (sex "F\|M"); P2-B7 group test is unadjusted; normative-template code not in repo. Phase 1 statuses reconciled; Phase 2 plan revised after the Nov 9 feedback (P2.1–P2.9); docs/TEAM_TASKS.md updated | Assign TBD owners; start P2.2 and P2.7 |
| Nov 15 | Chuong | Follow-ups: no-blink coverage comparison (no-blink run = primary data on the same 106 subjects; C −3.00 → −2.41 points, FDR p 0.025 → 0.128); Bohyeon's inputs match the repo files by value (all five runs, `conflict.csv`); sex-complete sample rule for Phase 2 in config; limitations recorded; owners proposed (pending confirmation). Meeting PDF left untracked (author name in its metadata) | Confirm owners; continue P2.2 and P2.7 |
