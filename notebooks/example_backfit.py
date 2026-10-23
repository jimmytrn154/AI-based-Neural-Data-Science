# %% [markdown]
# # Example: back-fit one subject's Simon cue epochs (docs/SETUP.md)
# Loads config.yaml and the primary templates, keeps correct trials (config
# keep_trials), filters 2-20 Hz, resamples to 250 Hz, back-fits, and prints
# each state's share of samples in 0-800 ms after the stimulus.
# Run from the unzipped repo folder: `python notebooks/example_backfit.py [subject_id]`

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from src import io as pio, microstates as ms

cfg = pio.load_config()
subject_id = sys.argv[1] if len(sys.argv) > 1 else "sub-P001"
k = cfg["microstates"]["primary_k"]
templates = ms.load_templates(k, cfg)
print(f"Templates K = {k}: states {templates.cluster_names}")

# %% Cue-locked epochs and the matching trial rows
epochs = pio.load_epochs(subject_id, "Simon", "cue", cfg)
trials = pio.match_trials(epochs, pio.load_trials(subject_id, "Simon", cfg), "cue")
keep = ((trials["beh_accuracy"] == 1) & ~trials["bad_epoch"]
        & trials["beh_rt"].between(150, 1500)).to_numpy()
print(f"{subject_id}: {keep.sum()} of {len(trials)} trials kept")

# %% Filter 2-20 Hz on the full epoch, resample to 250 Hz, back-fit
inst = pio.epochs_to_mne(epochs)[np.flatnonzero(keep)]
l_freq, h_freq = cfg["preprocessing"]["bandpass_hz"]
inst.filter(l_freq, h_freq, verbose="error")
inst.resample(cfg["preprocessing"]["resample_hz"], verbose="error")
labels = ms.backfit(inst, templates, cfg).labels  # (n_trials, n_times); -1 = unlabeled
print("Labels:", labels.shape)

# %% Share of samples per state, 0-800 ms after the stimulus
window = (inst.times >= 0) & (inst.times <= 0.8)
for i, name in enumerate(templates.cluster_names):
    print(f"  {name}: {np.mean(labels[:, window] == i):.1%}")
