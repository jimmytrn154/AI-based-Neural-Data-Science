# AI-based Neural Data Science

Course project analyzing the **Iowa Parkinson's Disease EEG dataset**: 64-channel EEG from 100 people with Parkinson's disease (PD) and 49 healthy controls. Each person did up to four tasks, and the dataset includes their clinical scores (MoCA, UPDRS-III).

The data has already been preprocessed (filtered, re-referenced, ICA-cleaned), so analysis can start from the distributed files. The full reference is `Dataset_Guide.pdf`, which ships with the data. This README summarizes the parts you need most often.

> **The data is not in this repository.** `Dataset_ALL` is 165 GB (2,746 files). Copy it to your machine and point your code at it, for example with a `DATA_ROOT` variable. Do not commit any of it.

---

## Dataset at a glance

| | |
|---|---|
| Source | University of Iowa, Narayanan Lab (2017–2021), public on OpenNeuro under CC0 |
| Participants | 149 total: 100 PD, 49 control |
| Recording | BrainVision 64-channel actiCAP, 500 Hz |
| Recordings | 581, all with the same 64 channels in the same order |
| Complete subjects | 134 did all four tasks (91 PD, 43 control) |

### Participants

| | PD | Control |
|---|---|---|
| n | 100 | 49 |
| Sex (F / M) | 31 / 68 | 23 / 26 |
| Age | 68.5 ± 8.0 | 70.9 ± 7.5 |
| MoCA (cognition, /30) | 24.3 ± 4.0 (range 9–30) | 26.7 ± 1.8 (range 22–30) |
| UPDRS-III (motor severity) | 12.6 ± 7.1 | 0 by definition |

MoCA and UPDRS-III are continuous, so you can run regression and correlation as well as PD-vs-control comparisons. Use UPDRS only within the PD group, because controls are recorded as 0.

### Tasks

| Task | Subjects | What it measures | Length | Conditions (median trials) |
|---|---|---|---|---|
| **Rest** | 149 | Spontaneous activity, eyes open | ~146 s | none (no events) |
| **Simon** | 147 | Response conflict, cognitive control | ~880 s | congruent 192 / incongruent 192 |
| **Oddball** | 146 | Novelty detection across modalities | ~741 s | standard 192 / visual oddball 24 / auditory oddball 24 |
| **IntervalTiming** | 139 | Time production (3 s and 7 s targets) | ~1413 s | short 40 / long 40 |

Main EEG signatures reported for each task:

- **Simon**: midfrontal theta (4–7 Hz) at Cz/FCz, 400–750 ms; N2; ERN in response-locked epochs.
- **Oddball**: novelty P3 (P3a), fronto-central, 250–500 ms after the **pre-cue**, seen in the oddball − standard difference wave.
- **IntervalTiming**: midfrontal delta (1–4 Hz) at Cz, 0.15–0.45 s after the GO cue; CNV and ramping activity during the interval.
- **Rest**: alpha power and peak frequency, alpha/theta ratio, aperiodic (1/f) slope, connectivity.

---

## Folder layout

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

Subject IDs are `sub-P###` for PD and `sub-C###` for controls. Rest has no `trial_info/` or `epoched_data/` because it has no events. A smaller 56 GB copy, `Dataset_without_epoched`, has everything except `epoched_data/`. Section 4.1 of the guide shows how to cut identical epochs yourself.

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

Epoch `.mat` fields: `data` (channel × time × trial, µV), `times` (ms), `srate`, `labels`, `chanlocs`, `badEpoch`, `trials` (matching trial-table rows), `erpBaseline` (−200 to 0 ms, already subtracted).

---

## Signal specification

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

---

## Quick start (Python)

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

---

## Pitfalls

1. **Pz is not zero and is the last channel.** It was the online reference. It was restored and the data re-referenced to the average.
2. **`beh_rt` units differ.** Simon uses milliseconds and Oddball uses seconds. `rt_eeg` is always in seconds.
3. **`badEpoch` only marks epochs.** Bad epochs are still in the data, so filter them out yourself.
4. **No time-frequency baseline is applied.** Apply one after the transform, for example dB relative to −300 to −200 ms.
5. **No low-pass filter.** For smooth ERPs, apply your own (for example 30 Hz).
6. **Oddball: the manipulation is on the pre-cue, not the arrow.** For novelty effects, use `lock-precue` or `onset_precue`. The measured pre-cue → arrow interval is 464 ms, not the 500 ms stated in the original README.
7. **IntervalTiming has no original behavioral file.** Trials were reconstructed from EEG markers. GO markers are missing in 5.2% of trials.
8. **Average within each subject before averaging across subjects.** Otherwise subjects with more trials dominate the result.
9. **For any classifier, split train and test by subject, not by trial or epoch.** Splitting by trial lets one person's data appear in both sets and inflates accuracy.
10. **Use `subject_id` across tasks.** The original OpenNeuro numbers (`orig_*`) refer to different people in different tasks.
11. **When cutting your own epochs, use `floor(onset * srate)`, not `round`.** This matches how the distributed epochs were cut.

### Subject selection

- Check the `trials_simon`, `trials_oddball` and `trials_timing` columns in `subjects.csv`: `ok`, `partial` or `exclude`. Condition labels for 7 subject/task pairs could not be verified (`exclude`). Their EEG is still usable for analyses that do not need labels.
- A reasonable starting rule is to also drop recordings with more than 6 interpolated channels (`EEG.etc.iowa.nInterp`). Fix your rule before looking at results and report it.
- About 23% of Rest recordings had no eye component removed, so blinks may remain. Check `EEG.etc.iowa.nICremoved` if you use frontal channels.

---

## What the data can and cannot support

| Can | Cannot |
|---|---|
| PD vs. control comparisons | Source localization (template electrode positions, no MRI) |
| Correlation with MoCA, UPDRS-III, age | Comparison with raw data (not distributed; available on OpenNeuro) |
| Within-task condition contrasts | Comparison of preprocessing pipelines (one fixed version) |
| Within-subject cross-task comparisons (134 subjects) | Respiration or accelerometer analyses (channels removed) |
| ERP, time-frequency, spectral, connectivity analyses | Trial-level analysis for the 7 `exclude` subject/task pairs |
| Classification and regression models | |

### Sanity-check reference values

| Measure | Value |
|---|---|
| Simon RT, congruent / incongruent (correct trials) | 514.7 / 557.4 ms |
| Simon error rate, congruent / incongruent | 6.4% / 14.8% |
| Oddball RT, standard / visual / auditory oddball | 476 / 515 / 509 ms |
| Novelty P3 at Cz (oddball − standard, 250–500 ms) | +2.4 µV |
| Produced interval for the 7 s target | ~6.2 s (both groups) |
| Interpolated channels per recording (excl. Rest) | 1.83 ± 0.93 |

If your numbers are far from these, check your pipeline before reporting a finding.

---

## Citation

The raw data is CC0 from the University of Iowa (Narayanan Lab). Cite the original OpenNeuro datasets in any report or presentation:

| Dataset | Task | DOI |
|---|---|---|
| ds004574 | Cross-modal Oddball | [10.18112/openneuro.ds004574.v1.0.0](https://doi.org/10.18112/openneuro.ds004574.v1.0.0) |
| ds004579 | Interval Timing | [10.18112/openneuro.ds004579.v1.0.0](https://doi.org/10.18112/openneuro.ds004579.v1.0.0) |
| ds004580 | Simon Conflict | [10.18112/openneuro.ds004580.v1.0.0](https://doi.org/10.18112/openneuro.ds004580.v1.0.0) |
| ds004584 | Resting State | [10.18112/openneuro.ds004584.v1.0.0](https://doi.org/10.18112/openneuro.ds004584.v1.0.0) |

Preprocessing methods: FASTER (Nolan, Whelan & Reilly, 2010), ICLabel (Pion-Tonachini et al., 2019), EEGLAB (Delorme & Makeig, 2004).
