# %% [markdown]
# # Task-level microstate metrics (Chuong)
# Back-fit the group templates to continuous Rest and Simon recordings and
# save `results/features/microstates_task.csv` (one row per subject_id).
# Run from the repo root: `python notebooks/chuong_task_metrics.py`

# %%
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import pandas as pd
from joblib import Parallel, delayed

from src import io as pio, microstates as ms

cfg = pio.load_config()
k = cfg["microstates"]["draft_templates"]["k"]
templates = ms.load_templates(k, cfg)
states = templates.cluster_names
print("Templates:", ms.templates_path(k, cfg), "states", states, "(labels PROVISIONAL)")

# %% Subjects: the draft-template subjects for now; full cohort with final templates
ids, _ = ms.draft_subjects(cfg)

# %% Back-fit and metrics
# pycrostates' short-segment rejection takes several minutes per Simon
# recording, so subjects run in parallel (~1 GB RAM per worker).
N_JOBS = 8
rows = Parallel(n_jobs=N_JOBS, verbose=5)(
    delayed(ms.task_metrics)(sid, templates, cfg) for sid in ids
)
features = pd.DataFrame(rows)
out = pio.results_dir(cfg, "features", "microstates_task.csv")
out.parent.mkdir(parents=True, exist_ok=True)
features.to_csv(out, index=False)
print("Saved", out, features.shape)

# %% Sanity checks
for task in ["rest", "simon"]:
    print(f"\n== {task} ==")
    dur = features[[f"{task}_meandur_{s}" for s in states]]
    print("mean duration (ms), mean over subjects:", dur.mean().round(1).to_dict())
    print(f"  range {dur.min().min():.1f}-{dur.max().max():.1f} ms; "
          f"{(dur.stack().between(60, 120)).mean():.0%} of values in 60-120 ms")
    cov = features[[f"{task}_cov_{s}" for s in states]].sum(axis=1)
    print(f"coverage sum (%): {cov.min():.3f}-{cov.max():.3f}")
    for s in states:
        out_cols = [f"{task}_trans_{s}_to_{d}" for d in states if d != s]
        tot = features[out_cols].sum(axis=1)
        print(f"outgoing transitions from {s} sum: {tot.min():.6f}-{tot.max():.6f}")
    print(f"gev_total (%): {features[f'{task}_gev_total'].mean():.1f} "
          f"± {features[f'{task}_gev_total'].std():.1f}")
    print(f"unlabeled_frac: max {features[f'{task}_unlabeled_frac'].max():.2e}; "
          f"n_samples median {features[f'{task}_n_samples'].median():.0f}")
