# Brain-State Reorganization from Rest to Task as a Marker of Cognitive Decline in Parkinson's Disease

Team project for **BME54701 AI-based Neural Data Science (2026-2)**, UNIST.
Team: Tran Anh Chuong, Minseo, Seonwoo Kim, Bohyeon Rhee.

We analyze the **Iowa Parkinson's Disease (PD) EEG dataset**: 64-channel EEG from 100 people with Parkinson's disease (PD) and 49 healthy controls, each doing up to four tasks, with clinical scores (MoCA, UPDRS-III).

- **Full proposal:** `docs/proposal.tex`
- **Progress tracker:** `PROGRESS.md`
- **Dataset documentation:** `Dataset_Guide.pdf`, which ships with the data. The dataset reference below summarizes the parts we use most.
- **Dataset path for Chuong's machine:** `/Users/trananhchuong/OneDrive - VINACADEMY LLC/AI-based Neural Data Science/Dataset_ALL`

> **The data is not in this repository.** `Dataset_ALL` is 165 GB (2,746 files). Keep it on your machine and point the code at it with the `DATA_ROOT` setting in `config.yaml`. Never commit data, derived EEG files, or large intermediate results.

---

## Our project

### Overall question

Do changes in brain-state dynamics from **rest to task** carry information about cognitive function in PD, beyond what resting EEG and cheap behavioral tests already provide?

We study brain states with **EEG microstates**: the 4–5 recurring scalp patterns that the brain switches between every 60–120 ms. They are found by polarity-invariant modified k-means clustering of the maps at peaks of global field power (GFP), an unsupervised method from Lecture 5.

### Research questions

**Phase 1 (unsupervised, presented October 26):**

- **RQ1. Rest vs. task.** Do the same microstate classes describe rest and the Simon task? How do duration, coverage and transitions change from rest to Simon within the same people?
- **RQ2. Conflict.** Within Simon, do microstate time courses differ between congruent and incongruent trials, and when?
- **RQ3. Groups and cognition (preliminary).** Do these patterns differ between PD and controls, and do they relate to MoCA and Simon error rates?
- **Exploratory. Recovery after conflict.** Does the brain take longer to return to baseline after an incongruent trial, and can this be measured reliably before the next stimulus?

**Later phases:**

- **RQ4. Robustness (Phase 2).** Do the findings hold up across checks? Specifically:
  - the number of microstate classes (K);
  - possible blink contamination;
  - whether alpha activity explains them (spectral parameterization).
- **RQ5. Added value (Phase 3).** Do resting microstates, and then rest-to-task reorganization, improve prediction of cognitive impairment (MoCA < 26) beyond demographics and behavior? Are the predicted probabilities well calibrated?
- **RQ6. Generalization (Phase 4).** Do the reorganization measures generalize to Oddball and Timing in the same people? Do stronger models (e.g., EEGNet) improve on interpretable features?

### Phases

| Phase | When | Goal |
|---|---|---|
| 1. Data structure | Now → Oct 26 | Unsupervised description of brain states: rest vs. Simon, congruent vs. incongruent |
| 2. Consolidation | After Oct 26 | Act on the professor's feedback; robustness checks; alpha control; fix the final features |
| 3. Prediction | Mid-semester | Calibrated, subject-level test of the added value of EEG over behavior |
| 4. Extensions | Late semester | Other tasks; deep learning comparison; transition networks |
| Final | Week 16 | Final presentation |

Phases 2–4 will be adjusted after the professor's feedback on October 26.

### Phase 1 analysis plan

1. **Shared templates.** Pool GFP-peak maps from Rest and continuous Simon, with an equal number of peaks per subject and task. Run modified k-means for each subject, then again across subjects to get group templates. Fit K = 4–7; always report K = 4 and K = 5. As a check, also fit task-specific templates and compare them by spatial correlation.
2. **Task-level dynamics (RQ1).** Back-fit the templates to the full Rest and Simon recordings. For each subject, task and state, compute duration, occurrence, coverage, global explained variance (GEV), and the transition matrix.
3. **Condition-level dynamics (RQ2).** Use the stimulus-locked epochs, correct trials only, with bad epochs excluded.
   - Run topographic segmentation of each condition's ERP.
   - Back-fit the templates to single trials to get state-probability time courses.
   - Compare conditions with cluster-based permutation tests.
   - Pre-specified windows: 200–350 ms (N2) and 400–750 ms (theta).
4. **Pre-specified summary features.**
   - **Reorganization index:** Jensen–Shannon divergence between each subject's rest and Simon transition distributions.
   - **Coverage change:** Simon minus rest, per state.
   - **Conflict effect:** distance between incongruent and congruent state-probability vectors in each window.
5. **Recovery pilot (exploratory).** Use the response-locked epochs.
   - Track how quickly state probabilities return to baseline after incongruent vs. congruent trials, and after errors vs. correct responses.
   - Build PCA trajectories of the same data.
   - Keep the measure only if split-half reliability (odd vs. even trials, Spearman–Brown corrected) is adequate.
6. **Subject-level structure.** Embed subjects by their microstate features (PCA, UMAP) and cluster them with k-means. Check alignment with PD/control status and MoCA.
7. **Statistics.**
   - Wilcoxon signed-rank tests and linear mixed models (`metric ~ task * group + age + sex + (1 | subject)`).
   - Spearman partial correlations with MoCA and error rates, controlling for age.
   - UPDRS within PD only.
   - Benjamini–Hochberg FDR correction. The summary features are primary; everything else is exploratory.

### Team roles (Phase 1)

| Member | Responsibilities |
|---|---|
| Tran Anh Chuong | Shared preprocessing and back-fitting code, templates and choice of K, Rest and Simon task-level metrics |
| Seonwoo Kim | Simon ERP reproduction, condition-level segmentation and state-probability time courses |
| Minseo | Behavioral sanity checks and Simon features, recovery pilot |
| Bohyeon Rhee | Rest-vs-Simon comparisons, transition networks, subject clustering, PCA trajectories, PD/MoCA associations |

**The critical dependency is the templates.** Chuong releases draft templates (about 15 subjects, K = 4) by **Oct 21** and final templates by **Oct 22**. Everyone develops on the draft templates and re-runs once the final ones arrive.

---

## How we work

### Rules fixed before seeing results

- **Exclusion:**
  - recordings with more than 6 interpolated channels (`EEG.etc.iowa.nInterp`);
  - Simon subject/task pairs whose `trials_simon` status is `exclude`.
- **Blink sensitivity:** every Rest analysis is repeated without recordings with `nICremoved == 0`.
- **Preprocessing for microstates:** 1–40 Hz band-pass, downsampled to 250 Hz, average reference (already applied).
- **Covariates:** age and sex in all group comparisons.
- **IDs:** always `subject_id` from `subjects.csv`, never the `orig_*` numbers.

All of these live in `config.yaml`. Change them only by team agreement, and log every change in `PROGRESS.md`.

### Shared interfaces

- **Templates:** one file per K, e.g. `results/templates/templates_K4.fif`, saved from the pycrostates clustering object.
- **Back-fitted labels:** one file per subject and task, e.g. `results/labels/sub-P001_task-Simon_labels.npy`.
- **Feature tables:** one CSV per analysis, **one row per subject**, keyed by `subject_id`, e.g. `results/features/microstates_task.csv`. All tables are merged on `subject_id`.
- **Shared code:** everyone imports the loading, preprocessing and back-fitting functions from `src/` rather than copying them.

### Suggested repository layout

```
.
├── README.md
├── PROGRESS.md               progress tracker and decision log
├── config.yaml               DATA_ROOT, exclusion rules, filters, K
├── docs/
│   └── proposal.tex
├── src/
│   ├── io.py                 loading .set/.fdt, epochs, trial tables, subjects.csv
│   ├── preprocess.py         filtering, downsampling, exclusion
│   ├── microstates.py        GFP peaks, clustering, back-fitting, metrics
│   ├── features.py           summary features (reorganization, conflict effect, recovery)
│   └── stats.py              mixed models, partial correlations, FDR, permutation tests
├── notebooks/                exploration, one notebook per person or analysis
└── results/                  git-ignored except small CSVs and figures
    ├── templates/
    ├── labels/
    ├── features/
    └── figures/
```

### Practical tips

- Develop on 10–20 subjects first; run the full cohort only once the code works.
- Use the ready-made epoch files for condition analyses, and remember to drop `badEpoch` trials.
- Post a short daily update in the group chat. Record blockers and decisions in `PROGRESS.md`.
- **AI assistance is allowed only if prompt histories are submitted.** Keep all of them.

---

## Dataset reference

### Course context and ground rules

Every team gets the same preprocessed dataset and chooses its own question. The course reference papers (Lecture 2) fall into three groups: 3 original task papers, 4 clinical-outcome papers, and 7 machine-learning benchmark papers.

The original papers share one claim: PD patients fail to raise midfrontal 1–8 Hz activity in response to a cue, in proportion to cognitive decline (ρ = 0.34 with MoCA, 0.46 with executive function).

| Task | Paper | Finding | Number |
|---|---|---|---|
| Rest | Anjum 2024, *npj PD* (original) | LPC + PCA + distance ratio on 8 electrodes outperforms deep learning; only externally validated result | AUC 0.90, external ρ 0.68 |
| Simon | Sun 2024, *EMBC* | Connectivity beats single-channel waveforms and RT; cannot separate congruent from incongruent (51%) | 95% (trial-level split, inflated) |
| Oddball | Avola 2025, *DMSVIVA* | 11 models on one pipeline and split; XGBoost beat every deep model | 0.79, AUC 0.86 |
| Timing | Singh 2021, *npj PD* (original) | Timing variability links the EEG signal to cognition; behavior alone gives a baseline | ρ −0.39 with MoCA, d = 0.66 |

**Ground rules:**

- **Different problems.** PD vs. control and impaired vs. normal (MoCA < 26) are different problems; 11 of the 49 controls score below 26. Never compare accuracies across the two.
- **UPDRS-III within PD only.** Controls are recorded as 0. All PD participants were recorded on medication.
- **Split by subject** for any classifier or regressor.
- **Reproduce a reference value** before starting your own analysis.

### Dataset at a glance

| | |
|---|---|
| Source | University of Iowa, Narayanan Lab (2017–2021), public on OpenNeuro under CC0 |
| Participants | 149 total: 100 PD, 49 control |
| Recording | BrainVision 64-channel actiCAP, 500 Hz |
| Recordings | 581, all with the same 64 channels in the same order |
| Complete subjects | 134 did all four tasks (91 PD, 43 control) |

| | PD | Control |
|---|---|---|
| n | 100 | 49 |
| Sex (F / M) | 31 / 68 | 23 / 26 |
| Age | 68.5 ± 8.0 | 70.9 ± 7.5 |
| MoCA (cognition, /30) | 24.3 ± 4.0 (range 9–30) | 26.7 ± 1.8 (range 22–30) |
| UPDRS-III (motor severity) | 12.6 ± 7.1 | 0 by definition |

### Tasks

| Task | Subjects | What it measures | Length | Conditions (median trials) |
|---|---|---|---|---|
| **Rest** | 149 | Spontaneous activity, eyes open | ~146 s | none (no events) |
| **Simon** | 147 | Response conflict, cognitive control | ~880 s | congruent 192 / incongruent 192 |
| **Oddball** | 146 | Novelty detection across modalities | ~741 s | standard 192 / visual oddball 24 / auditory oddball 24 |
| **IntervalTiming** | 139 | Time production (3 s and 7 s targets) | ~1413 s | short 40 / long 40 |

**Simon** (our Phase 1 task): a colored square appears left or right, and participants respond by color. Incongruent trials (square on the side opposite the correct hand) create conflict.

| Simon behavior | PD | Control |
|---|---|---|
| Error rate, congruent | 7.4% | 4.5% |
| Error rate, incongruent | 17.6% | 8.4% |
| Simon effect (RT) | +50 ms | +34 ms |

Inter-trial interval is about 2.1 s, with responses at about 500–560 ms. The main EEG signatures are midfrontal theta (4–7 Hz) at Cz/FCz around 400–750 ms, the N2 at 200–350 ms, and the ERN in response-locked epochs.

Main signatures for the other tasks:

- **Oddball:** novelty P3 (P3a), fronto-central, 250–500 ms after the **pre-cue**, in the oddball − standard difference wave.
- **IntervalTiming:** midfrontal delta (1–4 Hz) at Cz, 0.15–0.45 s after the GO cue; CNV and ramping during the interval.
- **Rest:** alpha power and peak frequency, alpha/theta ratio, aperiodic (1/f) slope, connectivity.

### Folder layout

```
Dataset_ALL/                          (165 GB)
├── Dataset_Guide.pdf                 full documentation (34 pages)
├── subjects.csv                      one row per subject: group, clinical scores, ID mapping
├── continuous_data/                  preprocessed continuous EEG (EEGLAB .set + .fdt)
│   ├── rest/      sub-XXXX_task-Rest_desc-clean_eeg.set
│   ├── simon/     sub-XXXX_task-Simon_desc-clean_eeg.set
│   ├── oddball/   sub-XXXX_task-Oddball_desc-clean_eeg.set
│   └── timing/    sub-XXXX_task-IntervalTiming_desc-clean_eeg.set
├── trial_info/                       per-trial condition labels and behavior (.tsv)
│   ├── simon/     sub-XXXX_task-Simon_trials.tsv
│   ├── oddball/   sub-XXXX_task-Oddball_trials.tsv
│   └── timing/    sub-XXXX_task-IntervalTiming_trials.tsv
└── epoched_data/                     ready-made epochs (MATLAB v7.3 / HDF5)
    ├── simon/     ..._lock-{cue,resp}_epo.mat
    ├── oddball/   ..._lock-{precue,gocue,resp}_epo.mat
    └── timing/    ..._lock-{cue,go,resp}_epo.mat
```

Subject IDs are `sub-P###` for PD and `sub-C###` for controls. Rest has no `trial_info/` or `epoched_data/` because it has no events.

### File formats

| File | Format | Read with |
|---|---|---|
| `*_eeg.set` + `.fdt` | EEGLAB, MATLAB v7 | `pop_loadset` / `mne.io.read_raw_eeglab` |
| `*_epo.mat` | MATLAB v7.3 (HDF5) | `load` / `h5py` |
| `*_trials.tsv` | Tab-separated, missing values = `n/a` | `readtable` / `pandas.read_csv` |
| `subjects.csv` | CSV, missing values = `n/a` | `readtable` / `pandas.read_csv` |

Keep each `.set` file in the same folder as its `.fdt` file, or it will not load.

### Epoch windows

| Task | Lock | Window (s) | Anchored on |
|---|---|---|---|
| Simon | `cue` | −1.0 to +2.5 | Stimulus |
| | `resp` | −1.5 to +1.0 | Key press |
| Oddball | `precue` | −1.0 to +2.5 | Pre-cue (**use for novelty**) |
| | `gocue` | −1.0 to +2.5 | Arrow |
| | `resp` | −1.5 to +1.0 | Key press |
| IntervalTiming | `cue` | −2.0 to +12.0 | Instruction text |
| | `go` | −1.0 to +10.0 | GO box |
| | `resp` | −3.0 to +2.0 | Key press |

Epoch `.mat` fields:

- `data`: channel × time × trial, µV;
- `times`: ms;
- `srate`, `labels`, `chanlocs`;
- `badEpoch`;
- `trials`: the matching trial-table rows;
- `erpBaseline`: −200 to 0 ms, already subtracted.

Simon condition coding:

- `beh_congruentcondition`: 1 = congruent, 2 = incongruent;
- `beh_accuracy`: 1 = correct, 2 = error, 99 = no response.

### Signal specification

| | |
|---|---|
| Channels | 64, amplifier order, **Pz is channel 64** |
| Sampling rate | 500 Hz |
| Units | µV (MNE returns volts, so multiply by `1e6`) |
| Reference | Average of all 64 channels (Pz restored from the online reference) |
| Filters | 0.1 Hz high-pass, line noise removed with CleanLine. **No low-pass** |
| Artifacts | Eye components removed with ICA + ICLabel. Muscle and heart components are **not** removed |
| Electrode positions | Standard 10-05 template, not digitized per subject |

Each `.set` file records its own processing in `EEG.etc.iowa`: interpolated channels, removed ICA components, ICA convergence and quality flags.

### Quick start (Python)

```python
import mne, pandas as pd, numpy as np, h5py

ROOT, SUB = "/path/to/Dataset_ALL", "sub-P001"

# Continuous data
raw = mne.io.read_raw_eeglab(
    f"{ROOT}/continuous_data/simon/{SUB}_task-Simon_desc-clean_eeg.set", preload=True)
data = raw.get_data() * 1e6                       # volts -> microvolts

# Trial labels
T = pd.read_csv(f"{ROOT}/trial_info/simon/{SUB}_task-Simon_trials.tsv",
                sep="\t", na_values="n/a")

# Ready-made epochs (HDF5; axes come out reversed)
with h5py.File(f"{ROOT}/epoched_data/simon/{SUB}_task-Simon_lock-cue_epo.mat", "r") as f:
    d = np.array(f["data"])                       # [trial, time, channel]
    t = np.array(f["times"]).ravel()              # ms

# Subject metadata and clinical scores
S = pd.read_csv(f"{ROOT}/subjects.csv", na_values="n/a")
```

Look channels up by name and never hard-code an index: `raw.ch_names.index("Cz")`.

### Pitfalls

**From the dataset:**

1. **Pz is not zero and is the last channel.** It was the online reference; it was restored and the data re-referenced to the average.
2. **`beh_rt` units differ.** Simon uses milliseconds and Oddball uses seconds. `rt_eeg` is always in seconds.
3. **`badEpoch` only marks epochs.** Bad epochs are still in the data, so filter them out yourself.
4. **No time-frequency baseline is applied.** Apply one after the transform, e.g. dB relative to −300 to −200 ms.
5. **No low-pass filter.** For smooth ERPs, apply your own (e.g. 30 Hz). For microstates we use 1–40 Hz.
6. **Oddball: the manipulation is on the pre-cue, not the arrow.** Use `lock-precue` or `onset_precue`. The pre-cue → arrow interval is 464 ms.
7. **IntervalTiming has no original behavioral file.** Trials were reconstructed from EEG markers, and GO markers are missing in 5.2% of trials.
8. **Average within each subject before averaging across subjects.**
9. **Split train and test by subject**, never by trial or epoch.
10. **Use `subject_id` across tasks.** The `orig_*` numbers refer to different people in different tasks.
11. **When cutting your own epochs, use `floor(onset * srate)`, not `round`.**

**Specific to our microstate analysis:**

12. **Don't compute per-condition durations inside short epochs.** Simon conditions are interleaved about every 2 s, so durations get cut off at the window edges. Use state-probability time courses for condition comparisons.
13. **Blinks create strong frontal maps** that can look like a microstate. About 23% of Rest recordings had no eye component removed (`nICremoved == 0`), so run the sensitivity check.
14. **Use the same templates for every task** when comparing rest and task. Otherwise state labels are not comparable.

### Reference values (sanity checks)

| Measure | Value |
|---|---|
| Simon RT, congruent / incongruent (correct trials) | 514.7 / 557.4 ms |
| Simon error rate, congruent / incongruent | 6.4% / 14.8% |
| Oddball RT, standard / visual / auditory oddball | 476 / 515 / 509 ms |
| Novelty P3 at Cz (oddball − standard, 250–500 ms) | +2.4 µV |
| Produced interval for the 7 s target | ~6.2 s (both groups) |
| Interpolated channels per recording (excl. Rest) | 1.83 ± 0.93 |

If your numbers are far from these, check your pipeline before reporting a finding.

### What the data can and cannot support

| Can | Cannot |
|---|---|
| PD vs. control comparisons | Source localization (template electrode positions, no MRI) |
| Correlation with MoCA, UPDRS-III, age | Comparison with raw data (not distributed; available on OpenNeuro) |
| Within-task condition contrasts | Comparison of preprocessing pipelines (one fixed version) |
| Within-subject cross-task comparisons (134 subjects) | Respiration or accelerometer analyses (channels removed) |
| ERP, time-frequency, spectral, connectivity, microstate analyses | Trial-level analysis for the 7 `exclude` subject/task pairs |
| Classification and regression models | |

---

## Citation

The raw data is CC0 from the University of Iowa (Narayanan Lab). Cite the original OpenNeuro datasets in any report or presentation:

| Dataset | Task | DOI |
|---|---|---|
| ds004574 | Cross-modal Oddball | [10.18112/openneuro.ds004574.v1.0.0](https://doi.org/10.18112/openneuro.ds004574.v1.0.0) |
| ds004579 | Interval Timing | [10.18112/openneuro.ds004579.v1.0.0](https://doi.org/10.18112/openneuro.ds004579.v1.0.0) |
| ds004580 | Simon Conflict | [10.18112/openneuro.ds004580.v1.0.0](https://doi.org/10.18112/openneuro.ds004580.v1.0.0) |
| ds004584 | Resting State | [10.18112/openneuro.ds004584.v1.0.0](https://doi.org/10.18112/openneuro.ds004584.v1.0.0) |

Methods to cite:

- **Preprocessing:** FASTER (Nolan, Whelan & Reilly, 2010), ICLabel (Pion-Tonachini et al., 2019), EEGLAB (Delorme & Makeig, 2004).
- **Our analysis:** MNE-Python (Gramfort et al., 2013), pycrostates (Férat et al., 2022), microstate segmentation (Pascual-Marqui et al., 1995; Michel & Koenig, 2018).
