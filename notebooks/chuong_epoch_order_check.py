# %% [markdown]
# # Epoch-order check (Chuong; PROGRESS.md blocker)
# For every Simon subject, rebuild the cue- and response-locked epochs from
# the continuous data and the trials.tsv onsets, following Dataset Guide
# section 4.1 (sample = floor(onset * srate), same window, -200 to 0 ms
# baseline), and compare them trial by trial with the distributed epochs in
# the order assumed by src.io.match_trials (cue: every row; resp: rows with
# an onset_response).
# Run from the repo root: `python notebooks/chuong_epoch_order_check.py`

# %%
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from src import io as pio

cfg = pio.load_config()
LOCKS = {"cue": ("onset", (-1.0, 2.5)), "resp": ("onset_response", (-1.5, 1.0))}
BASELINE_MS = (-200, 0)
TOLERANCE_UV = 0.01  # a trial "matches" if every sample differs by less than this

subjects = pio.load_subjects(cfg)
ids = subjects.loc[subjects["has_simon"] == 1, "subject_id"].tolist()
print(f"Simon subjects: {len(ids)}")


def rebuild(data, sfreq, onsets, window):
    """Guide section 4.1 epochs: (n_trials, n_channels, n_times), baseline-corrected."""
    n0, n1 = int(round(window[0] * sfreq)), int(round(window[1] * sfreq))
    times = np.arange(n0, n1) / sfreq * 1000
    base = (times >= BASELINE_MS[0]) & (times <= BASELINE_MS[1])
    epochs = np.full((len(onsets), data.shape[0], n1 - n0), np.nan)
    for i, onset in enumerate(onsets):
        s = math.floor(onset * sfreq)
        if s + n0 < 0 or s + n1 > data.shape[1]:
            continue
        seg = data[:, s + n0:s + n1]
        epochs[i] = seg - seg[:, base].mean(axis=1, keepdims=True)
    return epochs, times


def check(sid):
    raw = pio.load_continuous(sid, "Simon", cfg)
    data, sfreq = raw.get_data(units="uV"), raw.info["sfreq"]
    trials = pio.load_trials(sid, "Simon", cfg)
    rows = []
    for lock, (column, window) in LOCKS.items():
        epochs = pio.load_epochs(sid, "Simon", lock, cfg)
        matched = pio.match_trials(epochs, trials, lock)  # assumed order
        mine, times = rebuild(data, sfreq, matched[column].to_numpy(), window)
        same_channels = list(raw.ch_names) == list(epochs["ch_names"])
        same_times = np.allclose(times, epochs["times"])
        diff = np.abs(mine - epochs["data"]).reshape(len(mine), -1).max(axis=1)
        rows.append({"subject_id": sid, "lock": lock, "n_trials": len(mine),
                     "max_abs_diff_uv": float(np.nanmax(diff)),
                     "n_mismatch": int(np.sum(~(diff < TOLERANCE_UV))),
                     "n_not_rebuilt": int(np.isnan(diff).sum()),
                     "same_channels": same_channels, "same_times": same_times})
    return rows


results = pd.DataFrame([row for rows in Parallel(n_jobs=6, verbose=2)(delayed(check)(sid) for sid in ids)
                        for row in rows])
out = pio.results_dir(cfg, "qc", "epoch_order_check.csv")  # git-ignored
out.parent.mkdir(parents=True, exist_ok=True)
results.to_csv(out, index=False)
print("Saved", out)

# %% Summary
for lock, df in results.groupby("lock"):
    print(f"\n{lock}: {len(df)} subjects, {df.n_trials.sum()} trials; "
          f"max |difference| {df.max_abs_diff_uv.max():.2e} µV; "
          f"mismatching trials {df.n_mismatch.sum()} (in {(df.n_mismatch > 0).sum()} subjects); "
          f"not rebuilt (window outside recording) {df.n_not_rebuilt.sum()}; "
          f"channel order identical in all: {df.same_channels.all()}; times identical: {df.same_times.all()}")
bad = results[results.n_mismatch > 0]
if len(bad):
    print("\nSubjects with mismatching trials:\n" + bad.to_string(index=False))
