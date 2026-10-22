# %% [markdown]
# # Task-level microstate metrics (Chuong)
# Back-fit the group templates to continuous Rest and Simon recordings and
# save, one row per subject_id:
# - `results/features/microstates_task.csv` (temporal smoothing, primary);
# - `results/features/microstates_task_nosmooth.csv` (no smoothing, sensitivity).
# Run from the repo root: `python notebooks/chuong_task_metrics.py`

# %%
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import numpy as np
import pandas as pd

from src import io as pio, microstates as ms

cfg = pio.load_config()
k = cfg["microstates"]["draft_templates"]["k"]
sfreq = cfg["preprocessing"]["resample_hz"]
templates = ms.load_templates(k, cfg)
states = templates.cluster_names
print("Templates:", ms.templates_path(k, cfg), "states", states, "(labels PROVISIONAL)")
print("Back-fit options:", ms.backfit_options(cfg))

# %% Subjects: the draft-template subjects for now; full cohort with final templates
ids, _ = ms.draft_subjects(cfg)

# %% Back-fit and metrics, with and without smoothing
SETTINGS = {"smooth": ("microstates_task.csv", True),
            "nosmooth": ("microstates_task_nosmooth.csv", False)}
features, timing = {}, []
for name, (filename, smoothing) in SETTINGS.items():
    rows = []
    for sid in ids:
        start = time.perf_counter()
        rows.append(ms.task_metrics(sid, templates, cfg, smoothing=smoothing))
        timing.append({"subject_id": sid, "setting": name,
                       "seconds": time.perf_counter() - start})
    features[name] = pd.DataFrame(rows)
    out = pio.results_dir(cfg, "features", filename)
    out.parent.mkdir(parents=True, exist_ok=True)
    features[name].to_csv(out, index=False)
    print("Saved", out, features[name].shape)

timing = pd.DataFrame(timing).pivot(index="subject_id", columns="setting", values="seconds")
print("\nRuntime per subject (s, Rest + Simon incl. loading and filtering):")
print(timing.round(1).to_string())

# %% Sanity checks
def sanity(df, task):
    dur = df[[f"{task}_meandur_{s}" for s in states]]
    cov = df[[f"{task}_cov_{s}" for s in states]].sum(axis=1)
    trans = [df[[f"{task}_trans_{s}_to_{d}" for d in states if d != s]].sum(axis=1)
             for s in states]
    trans = pd.concat(trans, axis=1)
    return {
        "meandur by state (ms)": dur.mean().round(1).values.tolist(),
        "meandur range (ms)": f"{dur.min().min():.1f}-{dur.max().max():.1f}",
        "% values in 60-120 ms": f"{dur.stack().between(60, 120).mean():.0%}",
        "coverage sum (%)": f"{cov.min():.3f}-{cov.max():.3f}",
        "outgoing transitions sum": f"{trans.min().min():.6f}-{trans.max().max():.6f}",
        "gev_total (%)": f"{df[f'{task}_gev_total'].mean():.1f} ± {df[f'{task}_gev_total'].std():.1f}",
        "unlabeled_frac max": f"{df[f'{task}_unlabeled_frac'].max():.1e}",
    }


def duration_distribution(task, nosmooth):
    durs = np.concatenate([
        ms.segment_durations_ms(np.load(ms.labels_path(sid, task, cfg, nosmooth)), sfreq)
        for sid in ids])
    floor = cfg["microstates"]["backfit"]["min_segment_length_samples"] * 1000 / sfreq
    return {"n segments": durs.size, "median (ms)": float(np.median(durs)),
            "mean (ms)": round(float(durs.mean()), 1),
            f"% at floor ({floor:g}-{floor + 1000 / sfreq:g} ms)":
                f"{np.mean(durs <= floor + 1000 / sfreq):.0%}",
            "% >= 60 ms": f"{np.mean(durs >= 60):.0%}"}


for task in ["Rest", "Simon"]:
    t = task.lower()
    table = pd.DataFrame({name: sanity(df, t) for name, df in features.items()})
    dist = pd.DataFrame({name: duration_distribution(task, name == "nosmooth")
                         for name in SETTINGS})
    print(f"\n== {task} ==")
    print(pd.concat([table, dist.astype(str)]).to_string())

# %% Smoothing vs no smoothing, per metric (mean over subjects and states)
rows = []
for task in ["rest", "simon"]:
    for metric in ["meandur", "occ", "cov", "gev"]:
        cols = [f"{task}_{metric}_{s}" for s in states]
        a, b = features["smooth"][cols], features["nosmooth"][cols]
        rho = pd.Series(a.to_numpy().ravel()).corr(pd.Series(b.to_numpy().ravel()),
                                                   method="spearman")
        rows.append({"task": task, "metric": metric, "smooth": a.mean().mean(),
                     "nosmooth": b.mean().mean(), "spearman_rho": rho})
print("\nSmoothing vs no smoothing (mean over subjects × states; rho over subject × state values):")
print(pd.DataFrame(rows).round(3).to_string(index=False))
