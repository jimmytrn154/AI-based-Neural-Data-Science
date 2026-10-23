# %% [markdown]
# # Task-specific templates (Chuong)
# Rest-only and Simon-only two-level templates at one K (default: primary K),
# from the cached GFP peaks and the same subjects as the shared templates of
# that K; labeled by correlation with the shared templates.
# Run after chuong_final_templates.py, from the repo root:
# `python notebooks/chuong_task_specific.py [K]`

# %%
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import matplotlib
matplotlib.use("Agg")

import numpy as np
from pycrostates.io import ChData

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
k = int(sys.argv[1]) if len(sys.argv) > 1 else cfg["microstates"]["primary_k"]
summary_path = pio.results_dir(cfg, "templates", "summary.json")
with open(summary_path) as f:
    summary = json.load(f)

rest, _ = pp.select_subjects("Rest", cfg)
simon, _ = pp.select_subjects("Simon", cfg)
both = pp.select_rest_vs_simon(rest, simon, cfg)
fit_ids = [sid for sid in both if sid not in summary["K"][str(k)]["dropped"]]
shared = ms.load_templates(k, cfg)
print(f"K = {k}: {len(fit_ids)} subjects; shared labels {shared.cluster_names}")

# %% Same equal-size samples as the shared fit
sampled = {}
for sid in fit_ids:
    peaks = {task: ChData(np.load(pio.results_dir(cfg, "peaks", f"{sid}_task-{task}_peaks.npy")),
                          shared.info)
             for task in ("Rest", "Simon")}
    sampled[sid], _ = ms.equal_task_sample(peaks, sid, cfg)

# %% Fit, label, save
summary.setdefault("task_specific_by_k", {})[str(k)] = {}
for task in ("Rest", "Simon"):
    specific, corr, flags = ms.task_specific_templates(sampled, fit_ids, k, shared, cfg, task)
    ms.save_templates(specific, cfg, f"_{task.lower()}")
    diag = dict(zip(corr.index, np.diag(corr.to_numpy()).round(3)))
    print(f"\n{task}-only vs shared K = {k} (|r|):\n" + corr.round(2).to_string())
    print("Matched |r|:", diag, "| flags:", flags or "none")
    summary["task_specific_by_k"][str(k)][task] = {
        "corr": corr.round(3).to_dict(orient="index"), "matched_r": diag, "flags": flags}

with open(summary_path, "w") as f:
    json.dump(summary, f, indent=1, default=str)
