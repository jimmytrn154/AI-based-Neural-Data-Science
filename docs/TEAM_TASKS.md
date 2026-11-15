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

*Last updated: Nov 15, 2026 (Phase 2 revised after the professor's Nov 9 feedback; Phase 1 sections kept as a record below).*

---

## Phase 2 (revised Nov 15)

After the Phase 2 presentation (Nov 9) the professor asked for four validity checks (normative templates, probability dependence, C↔E transitions, randomization tests) and two questions (what the pattern says about brain activity, and what is new). The plan below answers them; statuses are tracked in `PROGRESS.md` (Phase 2).

**Naming.** The Nov 9 deck's new analyses "B7" and "B8" are called **P2-B7** (conflict state contrast) and **P2-B8** (conflict group windows) here, to keep them apart from the Phase 1 tasks B7 (robustness) and B8 (slides) below. Their scripts keep their names: `bohyeon 261112/src/B7_conflict_state_contrast.py`, `B8_conflict_group_windows.py` (see `bohyeon 261112/README.md`).

**Research question:** Phase 2 scope narrowed to validating and explaining the C-down / E-up pattern (rest -> Simon and congruent -> incongruent); the main research question is unchanged.

- **Main question** (answered in Phase 3): Do changes in brain-state dynamics from rest to task carry information about cognitive function in PD, beyond resting EEG and behavioral tests?
- **Sub-questions RQ1–RQ3** (done in Phase 1): rest vs. Simon dynamics; congruent vs. incongruent; groups and cognition.
- **Current focus** (Phase 2): is the C-down / E-up pattern real and robust, what is E, what does it mean, and is it new?

**Already done (Nov 9 deck, verified Nov 15):** E internal validation (spatial distinctness, recurrence across K and across tasks); coverage sensitivity (smoothing, QC, K = 4, no-blink for RQ1); P2-B7 and P2-B8.

| ID | Task | Owner | Status | Inputs | Outputs |
|---|---|---|---|---|---|
| P2.1 | E vs. normative templates: full match matrix (our 5 maps × all reference classes), A–D as control, baseline of non-matching reference correlations, all reference studies, common electrode subset or interpolation to our 64 positions, Microstate Toolbox templates (Custo 2017; Koenig 2024 meta-maps) | Minseo | 🟨 in progress (comparison done and verified Nov 15; meta-maps, coordinates and uncertainty remain) | `results/templates/templates_K5.fif` (+ K4–K7), published template files and electrode positions | `P2_1/` (code, README with sources), `results/features/p2_1_normative_matrix.csv`, `_E_summary.csv`, `_controls_AD.csv`, `_baseline.csv`, `_rest_simon.csv`, `results/figures/minseo_p2_1_normative_match.png` |
| P2.2 | Probability dependence: is E-up more than the compositional consequence of C-down? All five state changes; proportional-redistribution null; P(E) / (1 − P(C)); log-ratio analysis; per-template spatial correlation over time; second-best labels of frames switching to E | Chuong | ⬜ | `results/features/microstates_task.csv`, `results/labels/` (back-fit labels), Seonwoo's S4 probabilities, templates | `results/features/p2_2_probability_dependence.csv`, `results/figures/chuong_p2_2_probability_dependence.png` |
| P2.3 | C↔E transitions in the N2 and theta windows, incongruent vs. congruent: observed vs. expected given occurrences, directional asymmetry, origin of E segments | Seonwoo (pending confirmation) | ⬜ | Single-trial labels of the cue-locked epochs (as in Seonwoo's S4), `src/microstates.py` | `results/features/p2_3_ce_transitions.csv`, `results/figures/<owner>_p2_3_ce_transitions.png` |
| P2.4 | Randomization / surrogate tests on one fixed sample: group-label permutation for P2-B7 / P2-B8 group results (Freedman–Lane for age/sex-adjusted models; max-statistic correction across the six P2-B8 group tests); within-subject sign-flip for condition contrasts; segment-shuffle or Markov surrogates for transitions | Bohyeon | ⬜ | `conflict.csv`, `subjects.csv`, P2-B7 / P2-B8 scripts, P2.3 outputs | `results/features/p2_4_permutation_tests.csv`, `results/figures/bohyeon_p2_4_null_distributions.png` |
| P2.5 | Brain-activity interpretation: trial level (correct vs. error incongruent trials; reaction time), subject level (E increase vs. Simon error rate); then a cautious interpretation using P2.1–P2.3 | Minseo, with Seonwoo (pending confirmation) | ⬜ (depends on P2.1–P2.3) | Cue-locked epochs and trial tables, `results/features/behavior.csv`, P2.1–P2.3 outputs | `results/features/p2_5_behavior_link.csv`, `results/figures/<owner>_p2_5_behavior_link.png`, interpretation note `docs/p2_5_interpretation.md` |
| P2.6 | Novelty: targeted literature check (microstates in PD at rest and in tasks; microstates in conflict tasks such as Simon, Stroop, flanker; rest-to-task microstate changes; extra classes in clinical groups), then a "contribution" statement at the level the evidence supports | Minseo (pending confirmation) | ⬜ (claims only for findings that survive P2.4) | Literature; P2.4 results | `docs/p2_6_novelty.md` |
| P2.7 | Blink sensitivity: rerun the key results (conflict contrasts, group results) with `microstates_rest_noblink.csv` | Chuong | 🟨 RQ1 coverage done (no-blink run, Bohyeon; blink-cleaned run, Chuong, Nov 15): E holds, C weakens (see `PROGRESS.md`); conflict and group results not yet rerun | `results/features/microstates_rest_noblink.csv` | `results/features/microstates_rest_blinkclean.csv` (done), `results/features/p2_7_noblink_sensitivity.csv` |
| P2.8 | Alpha control: spectral parameterization of Rest (alpha peak and aperiodic features); do the main results survive them? | Seonwoo (pending confirmation) | ⬜ | Continuous Rest recordings | `results/features/p2_8_rest_spectral.csv`, `results/figures/<owner>_p2_8_spectral.png` |
| P2.9 | Recovery pilot: report the go/no-go result in the final (it was not in the Nov 9 deck) | Chuong | 🟨 (result in `docs/recovery_pilot.md`) | `results/features/recovery.csv`, `docs/recovery_pilot.md` | Slide in the final deck |

**Sample rule:** all Phase 2 group analyses (P2.4 onward) use the sex-complete sample (`config.yaml` → `phase2.group_sample`; 136 people for the Simon conflict features, sub-P069 excluded because its sex is recorded as "F|M").

**Dependencies:** P2.5 depends on P2.1–P2.3; novelty claims (P2.6) only for findings that survive P2.4; Phase 3 (prediction, led by Minseo) starts after P2.2 and P2.4.

**Postponed:** RT-sorted analysis of E; conflict adaptation; PCA trajectories (C5); Phase 4 extensions.

**Timeline:** final presentation in Week 16 (date: [TO FILL]); target Phase 2 done about two weeks before it, then Phase 3.

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

**Before back-fitting any epochs:** filter at **2–20 Hz** (`config.yaml`, `simon_conditions.filter`) and resample to **250 Hz** (`preprocessing.resample_hz`), as done in `notebooks/chuong_state_E.py`. In response-locked analyses, drop trials with `beh_accuracy == 99` (no response). `match_trials` gives the trial row of every epoch; the order was verified for all Simon subjects (Oct 24).

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

---

# Phase 1 (record)

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
| C4 | Recovery pilot: measure, comparison, reliability | Oct 24 | ✅ | Response- and cue-locked Simon epochs, templates, `backfit()` | `results/features/recovery.csv` (passing measures only; failed conflict measure in `recovery_conflict_unreliable.csv`, do not use), `results/figures/chuong_recovery.png`, go/no-go note (`docs/recovery_pilot.md`) |
| C5 | PCA trajectories of response-locked state probabilities (optional) | Phase 2 | ⏸️ (postponed Oct 24) | C4 outputs | `results/figures/chuong_recovery_pca.png` |
| C6 | Slides: pipeline, K choice, templates, state E, next phases | Oct 25 | ✅ (presented, per team; deck not in repo) | All of the above | Slides in the shared deck (draft: `docs/slides_chuong.md`) |

**`behavior.csv` specification:** one row per `subject_id`, Simon only.

- **Trials:** subjects with `trials_simon` in {ok, partial}. RTs come from correct trials only, with no RT limit. `beh_rt` is in **milliseconds** for Simon. Rules: `config.yaml` → `behavior`.
- **Columns:**

  | Column | Meaning |
  |---|---|
  | `simon_rt_con`, `simon_rt_inc` | Mean RT, correct trials, per condition |
  | `simon_effect_rt` | incongruent − congruent RT |
  | `simon_err_con`, `simon_err_inc` | **Primary** error rate = `beh_accuracy == 2` / responded trials (`beh_accuracy` in {1, 2}) |
  | `simon_effect_err` | incongruent − congruent error rate (primary) |
  | `simon_err_con_all`, `simon_err_inc_all`, `simon_effect_err_all` | Same with errors / all trials (Dataset Guide definition) |
  | `simon_noresp_rate` | `beh_accuracy == 99` / all trials |
  | `simon_n_con`, `simon_n_inc` | All trials of the condition |
  | `simon_nresp_con`, `simon_nresp_inc` | Responded trials of the condition |

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
| S1 | Load cue-locked epochs; keep correct trials with RT 150–1500 ms; drop bad epochs | Oct 23 | ✅ | `epoched_data/simon/*_lock-cue_epo.mat`, `trial_info/simon/*.tsv`, `src/io.py` (`load_epochs`, `load_trials`, `match_trials`) | — |
| S2 | Reproduce the Simon ERP at Cz (Dataset Guide, Figure 2) | Oct 23 | ✅ | S1 | `results/figures/seonwoo_erp_cz.png` (sanity check) |
| S3 | Topographic segmentation of each condition's ERP (0–800 ms) | Oct 24 | ✅ | S1 | `results/figures/seonwoo_erp_segmentation.png` |
| S4 | Single-trial back-fitting and state-probability time courses (A–E), per condition | Oct 24 | ✅ | S1, `templates_K5.fif`, `backfit()` (2–20 Hz, 250 Hz) | Per-subject time courses (saved as arrays) |
| S5 | Congruent vs. incongruent: paired cluster-based permutation test (5,000 permutations) | Oct 24 | ✅ | S4 | `results/figures/seonwoo_conflict_timecourse.png` |
| S6 | Conflict-effect feature in the 200–350 ms and 400–750 ms windows | Oct 24 | ✅ (definition to confirm) | S4 | `results/features/conflict.csv` (one row per `subject_id`) |
| S7 | Results slides for RQ2 | Oct 25 | ✅ (presented, per team; deck not in repo) | S2–S6 | Slides in the shared deck |

**Deliverables for Oct 26:** the ERP sanity-check figure, the condition time courses with significant clusters marked, and `conflict.csv`.

**Review, Nov 2 (Chuong):** received in `seonwoo_S1-S6_review/`; reviewed, not yet integrated into the repo.
- S1–S6 follow `config.yaml` (trial rules, ≥ 40 trials per condition, 2–20 Hz, 250 Hz, K = 5, documented smoothing, 5,000 permutations). ERP at Cz: 138 subjects; S3–S6: 137 (Simon recording QC applied).
- S5: two clusters survive family-wise correction: E incongruent > congruent, 464–576 ms (+1.6 points, p = 0.008); C incongruent < congruent, 348–512 ms (−1.2 points, p = 0.006).
- S6: `conflict.csv` reproduced exactly from the S4 probabilities. It holds per-state differences (10 columns), while the config's primary feature is one distance per window; the official definition is still to be agreed.
- S3 uses the shared K = 5 templates to label the condition ERPs (template-based segmentation), as documented in its methods note.
- Not to be merged: the bundled copies of `src/io.py`, `src/features.py` and `config.yaml` are older than the repo's.
- Open: S7 slides.

**Tip:** `notebooks/chuong_state_E.py` (function `state_probabilities`) already back-fits cue- and response-locked epochs on a 20-subject subset and produced `chuong_E_timecourse.png`. Reuse it as a starting point.

---

## Bohyeon Rhee: rest vs. Simon, networks, groups and cognition (RQ1, RQ3)

**Role:** within-subject rest-vs-Simon comparisons, reorganization measures, transition networks, subject clustering, and the PD/MoCA associations.

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| B1 | Statistics functions (mixed models, partial Spearman, FDR) | Oct 23 | ✅ | `config.yaml` (statistics) | `src/stats.py` |
| B2 | Rest vs. Simon within-subject comparisons (Wilcoxon, FDR, effect sizes); mixed models with group, age, sex | Oct 24 | ✅ | `microstates_task.csv`, `subjects.csv` | `results/features/rq1_stats.csv`, figure |
| B3 | Reorganization index (Jensen–Shannon divergence of transitions) and coverage change | Oct 24 | ✅ | `microstates_task.csv` (transition columns `<task>_trans_<X>_to_<Y>`) | `results/features/reorganization.csv` |
| B4 | Transition networks, Rest vs. Simon | Oct 24 | ✅ | `microstates_task.csv` | `results/figures/bohyeon_transition_networks.png` |
| B5 | Subject embedding (PCA, UMAP) and k-means; alignment with group and MoCA | Oct 24 | ✅ | `microstates_task.csv`, `subjects.csv` | Figure + short summary |
| B6 | PD vs. control and MoCA associations for the primary features (age, sex covariates; UPDRS within PD only; FDR) | Oct 24 | ✅ (small follow-up) | B2, B3, `behavior.csv` (from Chuong), `conflict.csv` (from Seonwoo) | `results/features/rq3_stats.csv` |
| B7 | Robustness: rerun the key RQ1/RQ3 results with `_K4`, `_nosmooth`, `_qcincluded` and `microstates_rest_noblink.csv` | Oct 24–25 | 🟨 (RQ1 rerun on all four sensitivity files, Nov 12; RQ3 not) | The sensitivity files | Short table: does each result hold? |
| B8 | Results slides for RQ1 and RQ3 | Oct 25 | ✅ (presented, per team; deck not in repo) | B2–B7 | Slides in the shared deck |

**Deliverables for Oct 26:** the RQ1 figure (rest vs. Simon per state), the reorganization index, the transition networks, and the preliminary group and MoCA associations.

**Review, Nov 7 (Chuong):** adjusted delivery received in `bohyeon_adjusted/` (code `src/B2_rest_vs_simon.py`, `B4_transition_network.py`, `B5_subject_clustering.py`, `B6_clinical_associations.py`, `reorganization.py`, updated `stats.py`; results `rq1_stats.csv`, new `rq1_mixed.csv`, `reorganization.csv`, `rq3_stats.csv`, four figures, B5 summary). All fixes requested on Nov 2 are done; key numbers reproduced from the repo files. Not yet integrated into the repo layout (the scripts read from a private `data/` folder).
- **B2:** Wilcoxon results unchanged. Mixed models added (`rq1_mixed.csv`; 60 models, 138 subjects, all converged; 36 used the Powell fallback, and all 60 refit identically with another optimizer). No task × group interaction survives FDR (smallest adjusted p = 0.57).
- **B2, secondary finding:** the mixed model's group term shows less state E at rest in PD: coverage −4.0 points (FDR p = 0.012), occurrence −0.34 per s (FDR p = 0.003); confirmed by a Rest-only regression with age and sex and by Mann–Whitney (medians 8.6% vs. 14.1%). It comes from a secondary RQ1 term on a feature outside the pre-specified RQ3 set, so present it as secondary or exploratory (team decision).
- **B3:** JSD definition documented in `reorganization.py`; `reorganization.csv` identical to the first delivery.
- **B4:** the networks show the 2 strongest outgoing edges per state; these are the same at Rest and in Simon, so the two diagrams look identical. The figure should state the rule.
- **B5:** state E included (11 features); k searched over 2–10 → k = 2, silhouette 0.21 (weak); no alignment with group (ARI 0.001) or MoCA. UMAP panel missing (umap-learn not installed on that machine).
- **B6:** one RQ3 family of 108 tests; mean durations and Mann–Whitney added; sub-P001 back in the conflict analyses (n = 137). Nothing survives FDR (smallest adjusted p = 0.10); E is shorter in PD (Simon p = 0.002, Rest p = 0.03, uncorrected).
- **Small follow-ups:** (1) B6 drops valid Rest data from the `rest_meandur_*` tests (sub-P048, whose Simon failed QC, and the 8 subjects with Simon trial labels `exclude`); (2) delete the unreachable `return output` at the end of `stats.py`; (3) send the outputs the scripts write but that were not included: `rq1_group_delta.csv`, `rq1_individual_changes.csv`, `b5_subject_clusters.csv`, `b6_eeg_features.csv`, `b6_group_stats.csv`, `b6_correlations.csv`, `bohyeon_b6_correlations.png`.

**Notes:**
- Descriptive preview (Chuong, Oct 23; paired Simon − Rest medians over 139 subjects, uncorrected, no tests): coverage of E +1.3 points (68% of subjects increase) and A +2.1 (60%); C coverage −1.8 (37% increase) and C mean duration −8.3 ms. The formal test is B2.
- sub-C049 has very high C coverage but passes QC. Keep it, and mention it in robustness checks if it drives any result.

---

## Minseo: presentation support now, prediction lead in Phase 3

**Role in Phase 1 (from Saturday, Oct 24):** light tasks that keep the presentation consistent.

| # | Task | Due | Status | Inputs | Output |
|---|---|---|---|---|---|
| M1 | Install the environment and load one subject | Oct 24 | ⬜ | `requirements.txt`, `README.md` | Confirmed in the daily log |
| M2 | Sanity-check slide: reproduced reference values vs. Dataset Guide | Oct 25 | ✅ (presented, per team; deck not in repo) | Chuong's C1–C2 numbers (`docs/behavior_sanity_checks.md`) | One slide |
| M3 | Fill the "Phase 1 results summary" table in `PROGRESS.md` from everyone's results | Oct 25 | ⬜ | All outputs | Updated table |
| M4 | Assemble the slide deck: consistent style, terminology (A–E, "task-enhanced"), numbers matching the files | Oct 25 | ✅ (presented, per team; deck not in repo) | Everyone's slides | Final deck |
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
