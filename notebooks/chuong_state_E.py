# %% [markdown]
# # State E, descriptive (Chuong)
# K = 5 state E (central maximum with opposite-sign temporal edges):
# 1. coverage, occurrence and mean duration at Rest vs Simon (all subjects);
# 2. state-probability time courses in cue- and response-locked Simon epochs
#    (20 subjects, stratified by group; config keep_trials; 2-20 Hz on the
#    full epoch, resampled to 250 Hz, back-fitted, then cropped);
# 3. spatial correlation with the first (1-40 Hz) draft templates.
# Descriptive only: no group or clinical comparisons.
# Run from the repo root: `python notebooks/chuong_state_E.py`

# %%
import copy
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
k = cfg["microstates"]["primary_k"]
templates = ms.load_templates(k, cfg)
states = templates.cluster_names
N_SUBJECTS = 20
WINDOWS_MS = {"cue": tuple(cfg["simon_conditions"]["analysis_window_ms"]),  # -200 to 800
              "resp": (-800, 600)}
CUE_ERP_PEAK_MS = 440  # central positivity at Cz, Dataset Guide

# %% 1. E at Rest vs Simon (primary feature table)
features = pd.read_csv(pio.results_dir(cfg, "features", "microstates_task.csv"))
rows = []
for metric, unit in (("cov", "%"), ("occ", "per s"), ("meandur", "ms")):
    for task in ("rest", "simon"):
        x = features[f"{task}_{metric}_E"].dropna()
        rows.append({"metric": f"{metric} ({unit})", "task": task, "n": len(x),
                     "median": x.median(), "q25": x.quantile(.25), "q75": x.quantile(.75)})
print("State E across subjects:\n" + pd.DataFrame(rows).round(2).to_string(index=False))

# %% 2. Time courses in the Simon epochs
qc = pd.read_csv(pio.results_dir(cfg, "templates", "recording_qc.csv"))
simon_ok = set(qc[(qc.k == k) & (qc.task == "Simon") & ~qc.excluded].subject_id)
pool = pio.load_subjects(cfg).set_index("subject_id").loc[sorted(simon_ok)].reset_index()
subset = ms.stratified_sample(pool, N_SUBJECTS, cfg)
print(f"\nEpoch subset ({len(subset)}):", subset,
      pool.set_index("subject_id").loc[subset, "group"].value_counts().to_dict())


def keep_trials(trials):
    """config simon_conditions.keep_trials: correct, not bad, 150-1500 ms RT."""
    return ((trials["beh_accuracy"] == 1) & ~trials["bad_epoch"]
            & trials["beh_rt"].between(150, 1500)).to_numpy()


def state_probabilities(sid, lock):
    """P(state) per time point, averaged over kept trials, for one subject."""
    epochs = pio.load_epochs(sid, "Simon", lock, cfg)
    trials = pio.match_trials(epochs, pio.load_trials(sid, "Simon", cfg), lock)
    keep = keep_trials(trials)
    inst = pio.epochs_to_mne(epochs)[np.flatnonzero(keep)]
    l_freq, h_freq = cfg["preprocessing"]["bandpass_hz"]
    inst.filter(l_freq, h_freq, verbose="error")
    inst.resample(cfg["preprocessing"]["resample_hz"], verbose="error")
    labels = ms.backfit(inst, templates, cfg).labels  # (n_trials, n_times)
    times = inst.times * 1000
    lo, hi = WINDOWS_MS[lock]
    window = (times >= lo) & (times <= hi)
    probs = np.stack([(labels[:, window] == i).mean(axis=0) for i in range(len(states))])
    rt = trials.loc[keep, "beh_rt"].median()
    return probs, times[window], int(keep.sum()), rt


curves = {}
for lock in WINDOWS_MS:
    out = Parallel(n_jobs=8)(delayed(state_probabilities)(sid, lock) for sid in subset)
    curves[lock] = {"probs": np.stack([o[0] for o in out]), "times": out[0][1],
                    "n_trials": [o[2] for o in out], "rt": [o[3] for o in out]}
    print(f"{lock}: trials kept per subject median {int(np.median(curves[lock]['n_trials']))} "
          f"(range {min(curves[lock]['n_trials'])}-{max(curves[lock]['n_trials'])})")

e = states.index("E")
median_rt = float(np.median(curves["cue"]["rt"]))
for lock, c in curves.items():
    mean = c["probs"][:, e].mean(axis=0)
    t = c["times"]
    if lock == "cue":
        base = mean[(t >= -200) & (t < 0)].mean()
        post = (t >= 0) & (t <= 800)
        i = np.argmax(mean[post])
        print(f"\ncue: P(E) baseline (-200-0 ms) {base:.3f}; peak {mean[post][i]:.3f} at "
              f"{t[post][i]:.0f} ms; mean in 300-500 ms {mean[(t >= 300) & (t <= 500)].mean():.3f}; "
              f"median RT of these trials {median_rt:.0f} ms")
    else:
        i = np.argmax(mean)
        print(f"resp: P(E) peak {mean[i]:.3f} at {t[i]:.0f} ms; at keypress "
              f"{mean[np.argmin(np.abs(t))]:.3f}; window mean {mean.mean():.3f}")
    peaks = {s: f"{c['times'][np.argmax(c['probs'][:, j].mean(0))]:.0f} ms"
             for j, s in enumerate(states)}
    print(f"  {lock}: time of maximum P(state), all states: {peaks}")

# %% Figure
fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
colors = plt.get_cmap("tab10")
for ax, (lock, c) in zip(axes, curves.items()):
    for j, s in enumerate(states):
        mean = c["probs"][:, j].mean(axis=0)
        sem = c["probs"][:, j].std(axis=0, ddof=1) / np.sqrt(c["probs"].shape[0])
        width = 2.5 if s == "E" else 1.0
        ax.plot(c["times"], mean, color=colors(j), lw=width, label=s, alpha=1 if s == "E" else 0.6)
        if s == "E":
            ax.fill_between(c["times"], mean - sem, mean + sem, color=colors(j), alpha=0.25)
    ax.axvline(0, color="k", lw=0.8)
    if lock == "cue":
        ax.axvline(CUE_ERP_PEAK_MS, color="grey", ls="--", lw=0.8)
        ax.axvline(median_rt, color="grey", ls=":", lw=0.8)
        ax.set_title(f"Stimulus-locked (dashed: {CUE_ERP_PEAK_MS} ms Cz positivity; "
                     f"dotted: median RT)", fontsize=9)
        ax.set_xlabel("Time from stimulus (ms)")
    else:
        ax.set_title("Response-locked (0 = keypress)", fontsize=9)
        ax.set_xlabel("Time from keypress (ms)")
    ax.set_ylabel("State probability (mean over subjects; E ± SEM)")
axes[1].legend(title=f"K = {k}", fontsize=8)
fig.suptitle(f"State E time course, Simon correct trials, {len(subset)} subjects, "
             f"{'-'.join(f'{f:g}' for f in cfg['preprocessing']['bandpass_hz'])} Hz")
fig.tight_layout()
out = pio.results_dir(cfg, "figures", "chuong_E_timecourse.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("Saved", out)

# %% 3. E vs the first draft templates (1-40 Hz, 15 subjects incl. sub-P048)
first = copy.deepcopy(cfg)
first["preprocessing"]["bandpass_hz"] = [1.0, 40.0]
path = ms.templates_path(4, cfg, "_draft1to40")
if not path.exists():
    ids, _ = ms.draft_subjects(first, exclude={})
    group, _, _ = ms.fit_templates(ids, 4, first)
    print(f"Regenerated first draft: group GEV {group.GEV_:.4f} (Oct 21 run: 0.6815)")
    ms.save_templates(group, cfg, "_draft1to40")
draft1 = ms.load_templates(4, cfg, "_draft1to40")
corr = ms.spatial_corr(templates.cluster_centers_[[e]], draft1.cluster_centers_)[0]
print("|r| of K = 5 state E with first-draft maps 0-3:", dict(enumerate(corr.round(3))),
      "(map 3 = muscle-like lateral map)")
