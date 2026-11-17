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
| Normative-template comparison (P2.1): our E vs. 12 published template sets (MS-Template-Explorer), centred Pearson \|r\| on common channels; best matches Zanesco 2020 E (0.868) and Zanesco 2021a E (0.843) | Minseo | ✅ verified Nov 15 | Code and README in `P2_1/`; reproduced from the original template files (`results/features/p2_1_normative_*.csv`) |

**Limitations of the work so far:**
- E partly splits across two maps at K = 7 (\|r\| = 0.82 with F, 0.83 with G); at K = 4–6 it has one clear counterpart.
- P2-B7's group difference in D is **unadjusted** (no age/sex); P2-B8 holds the adjusted version (N2 D −2.90 points, FDR p = 0.041).

All numbers above were checked against the outputs on Nov 15 (`bohyeon 261112/README.md`). The P2-B7 and P2-B8 samples differ by sub-P069 (sex recorded "F\|M", dropped from the sex-adjusted P2-B8 models).

### Revised plan (Nov 15)

| ID | Item | Feedback | Owner | Status | What remains |
|---|---|---|---|---|---|
| P2.1 | E vs. normative templates | 1 | Minseo | 🟨 in progress (verified Nov 15) | Done (`P2_1/`, verified): full 5 × K matrices for all 12 sets; A–D controls (A, B match the same letter in 12/12 sets, C 11, D 10; \|r\| 0.83–0.99); within-set baseline between different classes (median 0.42–0.77, max 0.68–0.88); E's best match is the set's E-type map in 3 of 12 sets (Zanesco 2020 0.868, Zanesco 2021a 0.843 with runner-up F 0.775, Hwang 2025 0.513), and C in the six 4-class sets without E; Custo 2017 (A–G) included: our E matches their D (0.72), their E 0.08; Rest-only and Simon-only E give the same matches; common-channel correlation equals re-referencing to the common subset; interpolation diagnostic only. Remaining: Koenig 2024 meta-maps (not among the 12); native electrode coordinates; uncertainty (bootstrap over subjects' maps); account for eyes-closed vs. our eyes-open/task data |
| P2.2 | Probability dependence: is E-up more than the compositional consequence of C-down? | 2 | Chuong | ✅ done (Nov 18) | Part 1 (Nov 17): E's share of non-C time rises and E rises in log-ratio terms, so E-up is more than a proportional redistribution of C's time. Part 2 (Nov 18): the fit to E rises on a non-competitive measure (mean \|r\|, q 0.012), and E samples are not more C-like in Simon, so E-up is not mainly relabelled C frames; small effects, not significant at GFP peaks. See the P2.2 notes below |
| P2.3 | C↔E transitions | 3 | Seonwoo (pending confirmation) | ⬜ | C → E and E → C in the N2 and theta windows, incongruent vs. congruent; observed vs. expected given occurrences; directional asymmetry; where E segments come from |
| P2.4 | Randomization / surrogate tests | 4 | Bohyeon | ⬜ | Group-label permutation for the P2-B7 / P2-B8 group results, Freedman–Lane for the age/sex-adjusted models and a max-statistic correction across the six P2-B8 group tests; within-subject sign-flip for condition contrasts; segment-shuffle or Markov surrogates for transitions; all on **one fixed sample** |
| P2.5 | Brain-activity interpretation | 5 | Minseo, with Seonwoo (pending confirmation) | ⬜ | Trial level: correct vs. error incongruent trials, and reaction time; subject level: E increase vs. Simon error rate; then a cautious interpretation using P2.1–P2.3 |
| P2.6 | Novelty: what is new? | 6 | Minseo (pending confirmation) | ⬜ | Targeted literature check: microstates in PD (rest and task); microstates in conflict tasks (Simon, Stroop, flanker); rest-to-task microstate changes; extra classes in clinical groups. Then a "contribution" statement, with claims at the level the evidence supports |
| P2.7 | Blink sensitivity | — | Chuong | ✅ done (Nov 17) | See the P2.7 notes below. RQ1 coverage: E holds in both checks; C's decrease weakens (FDR p 0.128 without the 33, 0.072 after cleaning). Group results with blink-cleaned Rest: all hold (PD less E at rest, coverage −4.07, FDR q 0.010; RQ3 still nothing significant). P2-B7 holds; P2-B8 keeps N2 Δ(C) and the C interaction, while three borderline results (q 0.041–0.044) move to q 0.059–0.069 with 9 fewer subjects and unchanged estimates |
| P2.8 | Alpha control: spectral parameterization of Rest | — | Seonwoo (pending confirmation) | ⬜ | Periodic (alpha peak) and aperiodic Rest features; test whether the main results survive them |
| P2.9 | Recovery pilot: report the go/no-go result | — | Chuong | 🟨 slide drafted (Nov 18) | Result exists (`docs/recovery_pilot.md`: keep overall recovery and error effect, drop the conflict effect); it was not in the Nov 9 deck. Slide drafted in `docs/slides_chuong_final.md` (slide 3, with P2.2 part 1 and P2.7 slides); to add to the shared deck |

**P2.7 note (Nov 15, from Bohyeon's outputs; report only).** Simon − Rest coverage change: mean [95% bootstrap CI], rank-biserial r, Cohen's d_z, FDR p (RQ1 family of 20 tests):

| State | Primary run, N = 139 | No-blink run, N = 106 | Primary data, same 106 subjects |
|---|---|---|---|
| A | +2.20 [+0.75, +3.67], r 0.29, d_z 0.26, p 0.012 | +2.12 [+0.50, +3.71], r 0.29, d_z 0.25, p 0.038 | identical to the no-blink run |
| C | −3.00 [−5.20, −0.58], r −0.26, d_z −0.21, p 0.025 | −2.41 [−5.03, +0.36], r −0.22, d_z −0.17, p 0.128 | identical to the no-blink run |
| E | +1.80 [+0.72, +2.92], r 0.40, d_z 0.27, p 0.0005 | +1.58 [+0.27, +2.78], r 0.34, d_z 0.24, p 0.019 | identical to the no-blink run |

The no-blink file holds the same per-subject values as the primary file, so the no-blink run equals the primary analysis restricted to those 106 subjects. In the 33 subjects it drops (Rest recordings with `nICremoved == 0`), mean changes are A +2.46, C −4.88, E +2.52 points. Recomputed with Bohyeon's `stats.py` and procedure, which reproduces her primary results exactly.

**P2.7 blink follow-up (Nov 15).**

*Who are the 33?* (Rest recordings with no eye component removed, vs. the other 106 of the 139.) No difference in group (PD/control 25/8 vs. 67/39, Fisher p = 0.21), sex (F/M 15/17 vs. 36/70, p = 0.21; F|M excluded), age (median 68 vs. 70, Mann–Whitney p = 0.74) or MoCA (26 vs. 26, p = 0.31). Rest C coverage is higher (median 33.1% vs. 22.0%, p = 0.0002) and Rest E coverage lower (5.9% vs. 11.8%, p = 0.001); Rest A does not differ (22.0% vs. 21.9%, p = 0.77).

*Cleaning* (rule fixed in config `phase2.blink_cleaning` before any coverage result: 0.5-s windows with Fp1/Fp2 peak-to-peak > 100 µV in the 2–20 Hz data marked BAD_blink and excluded from back-fitting; `notebooks/chuong_blink_clean.py`, `results/features/microstates_rest_blinkclean.csv`). Applied to all 35 QC-passing Rest recordings without eye-component removal; 20 had time marked (median 0.7%, IQR 0–10.7%, max 42%), 19 of them among the 33. Unmarked recordings reproduce their primary values exactly. After cleaning, the 33 still have higher Rest C (median 31.8% vs. 22.0%, p = 0.008), and their mean C change is −1.77 points (−4.88 before cleaning).

Simon − Rest coverage, all 139 subjects (same format as the table above):

| State | Primary | Blink-cleaned Rest |
|---|---|---|
| A | +2.20 [+0.75, +3.67], r 0.29, d_z 0.26, p 0.012 | +1.95 [+0.47, +3.42], r 0.26, d_z 0.22, p 0.040 |
| C | −3.00 [−5.20, −0.58], r −0.26, d_z −0.21, p 0.025 | −2.26 [−4.56, +0.13], r −0.21, d_z −0.16, p 0.072 |
| E | +1.80 [+0.72, +2.92], r 0.40, d_z 0.27, p 0.0005 | +1.67 [+0.59, +2.81], r 0.37, d_z 0.25, p 0.0017 |

Limitation: the 2–20 Hz band attenuates slow blink waveforms, so the 100 µV rule marks time in only 19 of the 33 recordings; residual blink activity below the threshold is not removed.

**P2.7 remaining reruns (Nov 17; rules fixed beforehand in config `phase2.blink_sensitivity_reruns`).** Bohyeon's scripts unchanged; the primary reruns reproduce her committed results (largest difference 1e-13; the RQ1 mixed models converge with lbfgs here instead of powell, CI within 1.4e-4). Code `notebooks/chuong_p2_7_reruns.py`; table `results/features/p2_7_blink_sensitivity.csv`. A result holds if its sign and its side of FDR q = 0.05 are unchanged.

| Results | Check | Tests | Hold | FDR q < 0.05, primary → check |
|---|---|---|---|---|
| RQ1 mixed models, group and task × group terms | Blink-cleaned Rest | 40 | 39 | 2 → 2 |
| RQ3 (B6) group models and partial correlations, incl. reorganization | Blink-cleaned Rest | 108 | 106 | 0 → 0 |
| P2-B7 | Without the 9 Simon recordings with no eye component removed (N 137 → 128) | 6 | 6 | 2 → 2 |
| P2-B8 | Same (N 136 → 127) | 21 | 18 | 7 → 4 |

- **The three "fails" in RQ1/RQ3 are sign flips of null estimates** (p 0.82–0.98), not changes in conclusions.
- **PD has less E at rest** (RQ1 mixed model): coverage −3.96 → −4.07 points (q 0.012 → 0.010), occurrence −0.34 → −0.35 per s (q 0.003 → 0.002).
- **P2-B7:** theta D = Δ(E) − Δ(C) +1.23 → +1.13 points (q 0.007 → 0.019); N2 PD − control in D −2.80 → −2.82 (q 0.029 → 0.043).
- **P2-B8 holds:** N2 Δ(C) PD − control +2.46 → +2.47 points (q 0.005 → 0.010); group × window for C (q 0.002 → 0.004); controls' N2 Δ(C) and theta Δ(E) vs. zero.
- **P2-B8 no longer below 0.05:** N2 D PD − control −2.90 → −2.92 (q 0.041 → 0.062); theta Δ(E) PD − control −1.35 → −1.28 (q 0.044 → 0.069); controls' N2 D vs. zero (q 0.041 → 0.059). The estimates barely change, so this reflects fewer subjects (9 removed, 4 of them controls) rather than blink artifacts. These three were already borderline.

**P2.2 part 1 (Nov 17; definitions fixed beforehand in config `phase2.probability_dependence`).** `notebooks/chuong_p2_2_probability_dependence.py`, `results/features/p2_2_probability_dependence.csv`, `results/figures/chuong_p2_2_probability_dependence.png`. N = 139, Simon − Rest; mean [95% bootstrap CI], rank-biserial r, FDR q (one family of 12):

| Measure | A | B | C | D | E |
|---|---|---|---|---|---|
| Coverage change, observed (points) | +2.20 | −0.80 | −3.00 | −0.21 | +1.80 |
| Expected if C's lost time is shared in proportion to Rest coverage | +0.84 | +0.70 | — | +0.85 | +0.61 |
| Excess, observed − expected | +1.36 [+0.23, +2.42] | −1.50 [−2.59, −0.50] | — | −1.06 [−2.30, +0.19] | +1.19 [+0.10, +2.33] |
| Share of non-C time, X / (100 − C) | +1.64 [+0.16, +3.20], r 0.19, q 0.083 | −2.17 [−3.68, −0.75], r −0.24, q 0.028 | — | −1.30 [−2.90, +0.41], r −0.06, q 0.56 | **+1.83 [+0.35, +3.35], r 0.32, q 0.0035** |
| CLR (log units) | +0.05, q 0.14 | −0.08, q 0.037 | −0.16, q 0.0035 | +0.02, q 0.90 | **+0.18 [+0.10, +0.26], r 0.40, q 0.0003** |

log(E/X): E/A +0.13 (q 0.11), E/B +0.26 (r 0.39, q 0.0003), E/D +0.16 (r 0.26, q 0.018).

- **E-up is not just C-down redistributed:** a proportional share-out of C's time would raise E by about 0.6 points; it rises by 1.8, so about two-thirds is excess. E's share of non-C time rises (q 0.0035), E has the largest CLR increase of the five states, and E rises relative to B and D (not significantly relative to A, which also rises).
- **The rest of the reorganization:** B loses share (q 0.028) and A gains (q 0.083, not significant); D is flat.
- **Blink-cleaned Rest (sensitivity):** same conclusions: E share +1.80 (q 0.0048), CLR E +0.17 (q 0.0005), log E/B and E/D q ≤ 0.015; CLR C −0.13 (q 0.014).
- **Limit:** these are coverage (composition) results. Back-fitting assigns each sample to one state, so they cannot show whether frames lost by C are relabelled E; that is part 2 (spatial correlation over time, second-best labels).

**P2.2 part 2 (Nov 18; definitions fixed beforehand in config `phase2.probability_dependence_part2`).** `notebooks/chuong_p2_2_part2.py`, `results/features/p2_2_part2_probability_dependence.csv`, per-subject measures `results/features/p2_2_part2_subject_measures.csv` (no demographics), `results/figures/chuong_p2_2_part2.png`. N = 139, Simon − Rest; mean [95% bootstrap CI], rank-biserial r, FDR q (one family of 7). The correlations reproduce pycrostates' assignments exactly (checked on sub-C001).

| Measure | Rest | Simon | Change | r | q |
|---|---|---|---|---|---|
| Mean \|r\| with E (all labelled samples) | 0.285 | 0.292 | **+0.007 [+0.001, +0.012]** | 0.26 | 0.012 |
| Mean \|r\| with C | 0.461 | 0.445 | −0.016 [−0.026, −0.005] | −0.28 | 0.012 |
| Mean \|r\| with B / D / A | | | −0.010 / −0.012 / +0.000 | | 0.012 / 0.018 / 0.90 |
| E samples whose runner-up is C | 28.8% | 28.1% | −0.7 points [−1.3, −0.1] | −0.21 | 0.035 |
| C samples whose runner-up is E | 15.4% | 16.3% | **+1.0 points [+0.3, +1.6]** | 0.30 | 0.012 |

- **E rises on a non-competitive measure:** each sample counts for every template, so this is not forced by C-down. The fit to E rises while the fit to B, C and D falls. Effects are small (|r| changes below 0.02).
- **No sign that E's new frames are borderline C frames:** in Simon, E samples have C as runner-up slightly *less* often. C samples do have E as runner-up more often, but so do all states (runner-up E share, Simon − Rest: A +0.5, B +0.9, C +1.0, D +1.1 points; descriptive). E-like activity rises in the signal generally, not specifically at C's expense.
- **GFP peaks only (sensitivity):** B, C and D fall (p ≤ 0.0005); E +0.003, p 0.16 (not significant). The E increase is carried by samples away from the GFP peaks, so it is weaker where the maps are most stable.
- **Blink-cleaned Rest:** same conclusions (E +0.006, q 0.028; C −0.013, q 0.041; E → runner-up C q 0.043; C → runner-up E q 0.035).
- **Conclusion for P2.2:** E-up is more than C-down redistributed (part 1) and is not mainly C frames being relabelled E (part 2). Caveat: small effects, and not significant at GFP peaks.

**Review of Bohyeon's Nov 15 delivery (`1115 bohyeon/`; reviewed Nov 16 by Chuong).** Results and figures only: no code, no README. Five analyses, all congruent vs. incongruent Simon trials with PD vs. control:

| Folder | What it tests | Fits | Result | Caveat |
|---|---|---|---|---|
| `state_prob` | PD − control in incongruent − congruent coverage, 5 states × 3 windows, one FDR family of 15 | P2-B8, extended to all five states | Reproduces P2-B8's N2 Δ(C): PD − control +2.46 points (FDR p 0.013). Theta Δ(E) −1.35 points no longer survives (raw p 0.022, FDR p 0.17) | Primary FDR family for P2-B8 to be agreed (see Blockers) |
| `state_observed_runs` | PD − control in the raw count of state runs, incongruent − congruent | Not in the plan | 11 of 15 significant; PD lose more runs in every state and window | Raw counts, not normalized; consistent with PD having fewer correct incongruent trials. Likely a trial-count effect: drop, or redo as a rate |
| `transition_prob` | Incongruent − congruent transition probabilities, 20 transitions × 3 effects × 3 windows (families of 60) | Overlaps P2.3 | One survives: controls, theta, D → E (FDR p 0.040) | Not C ↔ E; P2.3 needs one owner and one method |
| `JSD_rest_reference` | Jensen–Shannon divergence of each condition's transitions from the subject's own Rest | New (links RQ1 reorganization to the conflict contrast) | No group differences; controls, theta: incongruent further from Rest than congruent (FDR p 0.030); group × condition null | Mostly null; minor supporting point at most |
| `JSD_behavior` | Incongruent-vs-congruent transition divergence vs. RT and error-rate conflict costs (partial Spearman, age and sex; FDR over 18) | P2.5, subject level | PD: ρp +0.26 to +0.34 (FDR p 0.022–0.048); controls null | **Confounded by transition count** (below) |

*Transition-count check on `JSD_behavior`* (Chuong, from Bohyeon's per-subject table; her values reproduce). Divergence from fewer transitions is larger by sampling noise, and subjects with larger conflict costs have fewer transitions (ρ with log total transitions: −0.43 to −0.62 for error cost, −0.34 to −0.44 for RT cost). Adding log total transitions as a covariate, PD:

| PD | Bohyeon's ρp | With transition count |
|---|---|---|
| Error cost, analysis window | +0.31 (p 0.003) | +0.04 (p 0.71) |
| Error cost, theta | +0.27 (p 0.010) | +0.08 (p 0.44) |
| RT cost, theta | +0.33 (p 0.002) | +0.20 (p 0.055) |
| RT cost, N2 | +0.27 (p 0.011) | +0.18 (p 0.094) |

The error-cost links disappear and the RT-cost links weaken to borderline. Not a finding yet; needs the transition count as a covariate, or a divergence measure corrected for trial count, before it goes to P2.5. Two per-subject files (`jsd_behavior_subject_values.csv`, `jsd_subject_values.csv`) hold group, age and sex and stay uncommitted.

**Dependencies:**
- P2.5 depends on P2.1–P2.3.
- Novelty claims (P2.6) only for findings that survive P2.4.
- Phase 3 starts after P2.2 and P2.4.

**Postponed (kept listed):** RT-sorted analysis of E; conflict adaptation; PCA trajectories (recovery pilot, C5); Phase 4 extensions.

**Timeline:** today is Nov 15; final presentation in Week 16 (date: [TO FILL]). Target: Phase 2 done about two weeks before the final, then Phase 3.

### Earlier Phase 2 placeholder items (Oct 26), for the record

- [x] Research questions after the feedback: Phase 2 scope narrowed to validating and explaining the C-down / E-up pattern (rest -> Simon and congruent -> incongruent); the main research question is unchanged.
- [ ] Robustness: K, blinks, filter choices (RQ1 at K = 4, no smoothing, QC included, no-blink; blinks done in P2.7)
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
| Nov 15 | P2.1 owner confirmed: Minseo (other Phase 2 owners still pending confirmation) | Chuong |
| Nov 15 | P2.7 blink cleaning rule, fixed before any coverage result: in Rest recordings with no eye component removed, mark 0.5-s windows with Fp1 or Fp2 peak-to-peak > 100 µV (2–20 Hz, 250 Hz data) as BAD_blink and exclude them from back-fitting. Threshold justified by the Fp amplitude distribution of recordings with eye components removed (config `phase2.blink_cleaning`) | Chuong |
| Nov 15 | P2.1 primary: centered Pearson \|r\| on common channels; interpolation diagnostic only | Chuong |
| Nov 17 | P2.7 remaining reruns defined before any result (config `phase2.blink_sensitivity_reruns`): Rest-based group results (RQ1 mixed-model group terms, reorganization, RQ3 group models and correlations) with blink-cleaned Rest; P2-B7 / P2-B8 without the Simon recordings that had no eye component removed; a result holds if its sign and its side of FDR q = 0.05 are unchanged | Chuong |
| Nov 17 | P2.2 part 1 defined before any result (config `phase2.probability_dependence`): N = 139 from `microstates_task.csv`; share of non-C time X / (100 − C) for A, B, D, E; redistribution excess reported with it (same null, not tested separately); CLR changes A–E and log(E/X) for X = A, B, D; paired Wilcoxon; one FDR family "P2.2" of 12 tests; blink-cleaned Rest as sensitivity | Chuong |
| Nov 18 | P2.2 part 2 defined before any result (config `phase2.probability_dependence_part2`): per recording, mean \|r\| of all labelled samples with each template A–E (GFP peaks as reported sensitivity); runner-up template of each sample (best of the four other than its primary label): share of E samples with runner-up C and of C samples with runner-up E; Simon − Rest, paired Wilcoxon, N = 139; one FDR family "P2.2 part 2" of 7 tests; blink-cleaned Rest as sensitivity | Chuong |

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
| Nov 15 | Normative-template comparison (Zanesco et al. 2020, 2021a) was reported in the Nov 9 deck, but its code is not in the repo and its owner is not recorded | P2.1 | Minseo | Yes (Nov 15): code and results in `P2_1/`, verified |
| Nov 16 | Bohyeon's Nov 15 delivery (`1115 bohyeon/`) has results and figures but no code or README, so its methods can't be checked (handoff rules ask for code with outputs) | Chuong (integration) | Bohyeon | No |
| Nov 16 | P2-B8 FDR family: with all five states (family of 15, `1115 bohyeon/state_prob`) only N2 Δ(C) survives; theta Δ(E) (FDR p 0.044 in P2-B8) does not (0.17). Agree which family is primary | P2-B8, P2.4 | Team | No |
| Nov 16 | Bohyeon's transition analysis (`1115 bohyeon/transition_prob`) overlaps P2.3 (Seonwoo, pending confirmation); one owner and one method needed | P2.3 | Team | No |
| Nov 16 | Trial-count confounds in Bohyeon's Nov 15 delivery: raw run counts (`state_observed_runs`) not normalized; transition divergence vs. conflict costs (`JSD_behavior`) explained largely by transition count (see the review note in Phase 2) | P2.5 | Bohyeon | No |

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
| Nov 15 | Chuong | Part A: meeting note committed without author metadata (`docs/meetings/ND_meeting_1107.pdf`); P2.1 owner confirmed (Minseo); `subjects.csv` and `behavior.csv` match Bohyeon's inputs by value | Part B (blink follow-up), Part C (P2.1 verification) |
| Nov 15 | Chuong | Part B (P2.7): the 33 Rest recordings without eye-component removal match the 106 on group, sex, age and MoCA but have higher Rest C (33.1% vs 22.0%); blink cleaning (Fp > 100 µV, 0.5 s; rule committed before results) marks time in 20 of 35 recordings; with cleaned Rest, Simon − Rest coverage: E +1.67 (FDR p 0.002), A +1.95 (0.040), C −2.26 (0.072) | Part C (P2.1 verification) |
| Nov 15 | Chuong | Part C (P2.1): Minseo's normative-template comparison integrated in `P2_1/` with a repo script and README (sources of all 12 template sets; .set files not committed). Method verified (centred \|r\|, name-matched common channels = common-average reference); Zanesco 2020 0.868 and 2021a 0.843 reproduced, all matrices within 5e-7. E's best match is an E-type map in 3 of 12 sets; A–D controls match in 10–12 of 12; within-set baseline max 0.68–0.88. Koenig 2024 meta-maps still missing | Next: P2.2; Minseo continues P2.1 (meta-maps, coordinates, uncertainty) |
| Nov 16 | Chuong | Reviewed Bohyeon's Nov 15 delivery (`1115 bohyeon/`, results only): `state_prob` reproduces P2-B8's N2 Δ(C) (theta Δ(E) does not survive a 15-test family); raw run counts likely a trial-count effect; transitions overlap P2.3 (only control theta D → E survives); divergence from Rest mostly null; divergence vs. conflict costs in PD largely explained by transition count. Per-member Phase 2 tables added to `docs/TEAM_TASKS.md` | Ask Bohyeon for the code; agree P2-B8 FDR family and P2.3 ownership |
| Nov 17 | Chuong | P2.2 part 1 (definitions fixed in config first): E's share of non-C time +1.83 points (FDR q 0.0035), CLR E +0.18 (q 0.0003), E rises relative to B and D; proportional redistribution of C's time explains only about a third of E-up (expected +0.61 vs observed +1.80). Same with blink-cleaned Rest | P2.2 part 2 (spatial correlation over time, second-best labels) |
| Nov 17 | Chuong | P2.7 done (rules fixed in config first): with blink-cleaned Rest, all RQ1 mixed-model and RQ3 results keep their conclusions (PD less E at rest, q 0.010); without the 9 Simon recordings lacking eye-component removal, P2-B7 holds and P2-B8 keeps N2 Δ(C) (q 0.010) and the C interaction, while three borderline results (q 0.041–0.044) move to 0.059–0.069 with unchanged estimates | P2.2 part 2 |
| Nov 18 | Chuong | P2.2 part 2 (definitions fixed in config first): mean \|r\| with E rises (+0.007, q 0.012) while C, B, D fall; E samples have C as runner-up less often in Simon (q 0.035), and every state gains E as runner-up, so E-up is not mainly relabelled C frames; E not significant at GFP peaks (p 0.16). Same with blink-cleaned Rest. P2.2 done. P2.9: slide drafted with P2.2 and P2.7 slides (`docs/slides_chuong_final.md`) | Add slides to the shared deck; integration |
