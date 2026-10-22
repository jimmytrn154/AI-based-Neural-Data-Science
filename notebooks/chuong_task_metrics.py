# %% [markdown]
# # Task-level microstate metrics, full cohort (Chuong)
# Back-fit the final templates to every eligible continuous Rest and Simon
# recording and save, one row per subject_id (results/features/):
# - microstates_task.csv            primary K, temporal smoothing (primary file)
# - microstates_task_nosmooth.csv   primary K, no smoothing (sensitivity)
# - microstates_task_K4.csv / _K5   K = 4 / 5 if not the primary K (sensitivity)
# - microstates_task_qcincluded.csv primary K with QC-excluded recordings kept
# - microstates_rest_noblink.csv    Rest columns without nICremoved == 0 recordings
# Recordings excluded by the recording QC (per K) are left empty (NaN).
# Run after notebooks/chuong_final_templates.py, from the repo root.

# %%
import logging
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
DRY = os.environ.get("MS_DRYRUN") == "1"  # same dry run as chuong_final_templates.py
if DRY:
    cfg["paths"]["results"] = os.environ["MS_DRYRUN_DIR"]
    cfg["microstates"]["primary_k"] = int(os.environ.get("MS_PRIMARY_K", 4))
primary_k = cfg["microstates"]["primary_k"]
sfreq = cfg["preprocessing"]["resample_hz"]
ks = [primary_k] + [k for k in (4, 5) if k != primary_k]
templates = {k: ms.load_templates(k, cfg) for k in ks}
print("Primary K:", primary_k, "| K run:", ks, "| back-fit:", ms.backfit_options(cfg))
print("State labels are PENDING VISUAL CONFIRMATION")

rest, _ = pp.select_subjects("Rest", cfg)
simon, _ = pp.select_subjects("Simon", cfg)
eligible = {"Rest": set(rest.subject_id), "Simon": set(simon.subject_id)}
qc = pd.read_csv(pio.results_dir(cfg, "templates", "recording_qc.csv"))
if DRY:
    eligible = {task: set(qc.loc[qc.task == task, "subject_id"]) for task in eligible}
ids = sorted(eligible["Rest"] | eligible["Simon"])
print(f"Subjects: {len(ids)} (Rest {len(eligible['Rest'])}, Simon {len(eligible['Simon'])})")

# %% Back-fit every recording once per template set and setting
def subject_rows(sid):
    """Feature rows for one subject: {(k, smoothing): row}."""
    rows = {}
    start = time.perf_counter()
    for task in ("Rest", "Simon"):
        if sid not in eligible[task]:
            continue
        raw = pp.preprocess_raw(pio.load_continuous(sid, task, cfg), cfg)
        for k in ks:
            for smoothing in ((True, False) if k == primary_k else (True,)):
                segmentation = ms.backfit(raw, templates[k], cfg, smoothing=smoothing)
                if k == primary_k:
                    path = ms.labels_path(sid, task, cfg, nosmooth=not smoothing)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    np.save(path, segmentation.labels)
                row = rows.setdefault((k, smoothing), {"subject_id": sid})
                row.update(ms.segmentation_metrics(segmentation, task.lower()))
    return rows, time.perf_counter() - start


results = Parallel(n_jobs=8, verbose=2)(delayed(subject_rows)(sid) for sid in ids)
seconds = pd.Series([r[1] for r in results], index=ids)
print(f"Runtime per subject: median {seconds.median():.1f} s, max {seconds.max():.1f} s, "
      f"total wall-clock see log")


def table(k, smoothing):
    df = pd.DataFrame([r[0][(k, smoothing)] for r in results if (k, smoothing) in r[0]])
    columns = ["subject_id"] + [c for task in ("rest", "simon") for c in df.columns
                                if c.startswith(task + "_")]
    return df[columns].set_index("subject_id")


def apply_qc(df, k):
    """Empty the columns of recordings excluded by the recording QC at this K."""
    df = df.copy()
    flagged = qc[(qc.k == k) & qc.excluded]
    for sid, task in zip(flagged.subject_id, flagged.task):
        df.loc[sid, [c for c in df.columns if c.startswith(task.lower() + "_")]] = np.nan
    return df


out_dir = pio.results_dir(cfg, "features")
out_dir.mkdir(parents=True, exist_ok=True)
features = {
    "microstates_task.csv": apply_qc(table(primary_k, True), primary_k),
    "microstates_task_nosmooth.csv": apply_qc(table(primary_k, False), primary_k),
}
for k in ks[1:]:
    features[f"microstates_task_K{k}.csv"] = apply_qc(table(k, True), k)
n_flagged = int(qc[(qc.k == primary_k) & qc.excluded].shape[0])
if n_flagged:
    features["microstates_task_qcincluded.csv"] = table(primary_k, True)

# Blink sensitivity: Rest without recordings where no eye component was removed
blink = rest.set_index("subject_id")["blink_flag"]
rest_cols = [c for c in features["microstates_task.csv"].columns if c.startswith("rest_")]
noblink = features["microstates_task.csv"][rest_cols].dropna(how="all")
noblink = noblink[~blink.reindex(noblink.index).fillna(True).astype(bool)]
features["microstates_rest_noblink.csv"] = noblink

for name, df in features.items():
    df.to_csv(out_dir / name)
    n = {task: int(df[f"{task}_gev_total"].notna().sum())
         for task in ("rest", "simon") if f"{task}_gev_total" in df}
    print(f"Saved {name}: {df.shape[0]} rows x {df.shape[1]} columns; recordings {n}")
print(f"Rest recordings with nICremoved == 0: {int(blink.sum())} of {len(blink)}")

# %% Sanity checks
def sanity(df, task):
    states = [c.rsplit("_", 1)[1] for c in df.columns if c.startswith(f"{task}_meandur_")]
    df = df[df[f"{task}_gev_total"].notna()]
    dur = df[[f"{task}_meandur_{s}" for s in states]]
    cov = df[[f"{task}_cov_{s}" for s in states]].sum(axis=1)
    trans = pd.concat([df[[f"{task}_trans_{s}_to_{d}" for d in states if d != s]].sum(axis=1)
                       for s in states], axis=1)
    return {
        "n recordings": len(df),
        "meandur by state (ms)": dict(zip(states, dur.mean().round(1))),
        "meandur range (ms)": f"{dur.min().min():.1f}-{dur.max().max():.1f}",
        "% values in 60-120 ms": f"{dur.stack().between(60, 120).mean():.0%}",
        "coverage sum (%)": f"{cov.min():.3f}-{cov.max():.3f}",
        "outgoing transitions sum": f"{trans.min().min():.6f}-{trans.max().max():.6f}",
        "gev_total (%)": f"{df[f'{task}_gev_total'].mean():.1f} ± {df[f'{task}_gev_total'].std():.1f}",
        "unlabeled_frac max": f"{df[f'{task}_unlabeled_frac'].max():.1e}",
    }


def duration_distribution(task, nosmooth):
    included = features["microstates_task.csv"][f"{task.lower()}_gev_total"].dropna().index
    durs = np.concatenate([
        ms.segment_durations_ms(np.load(ms.labels_path(sid, task, cfg, nosmooth)), sfreq)
        for sid in included])
    floor = cfg["microstates"]["backfit"]["min_segment_length_samples"] * 1000 / sfreq
    return {"n segments": durs.size, "median (ms)": float(np.median(durs)),
            "mean (ms)": round(float(durs.mean()), 1),
            "% at 20 ms floor (<= 24 ms)": f"{np.mean(durs <= floor + 1000 / sfreq):.0%}",
            "% >= 60 ms": f"{np.mean(durs >= 60):.0%}"}


for name, df in features.items():
    for task in ("rest", "simon"):
        if f"{task}_gev_total" not in df:
            continue
        print(f"\n== {name} | {task} ==")
        for key, value in sanity(df, task).items():
            print(f"  {key}: {value}")
        if name in ("microstates_task.csv", "microstates_task_nosmooth.csv"):
            for key, value in duration_distribution(task.capitalize(),
                                                    "nosmooth" in name).items():
                print(f"  {key}: {value}")

# %% sub-C049 and the distribution of state C coverage
main = features["microstates_task.csv"]
unmasked = table(primary_k, True)
print("\n== sub-C049 ==")
for task, table_ in (("Rest", rest), ("Simon", simon)):
    if "sub-C049" not in eligible[task]:
        continue
    row = table_.set_index("subject_id").loc["sub-C049"]
    q = qc[(qc.k == primary_k) & (qc.subject_id == "sub-C049") & (qc.task == task)].iloc[0]
    print(f"  {task}: nInterp {row['nInterp']}, nICremoved {row['nICremoved']}, "
          f"template GEV {q.gev_template:.3f}, robust z {q.gev_z:.2f}, excluded {q.excluded}")
for task in ("rest", "simon"):
    cov = main[f"{task}_cov_C"].dropna()
    c049 = unmasked[f"{task}_cov_C"].get("sub-C049", np.nan)
    mad = 1.4826 * np.median(np.abs(cov - cov.median()))
    print(f"  {task}_cov_C (n = {len(cov)}): median {cov.median():.1f}%, IQR {cov.quantile(.25):.1f}-"
          f"{cov.quantile(.75):.1f}%, 5-95% {cov.quantile(.05):.1f}-{cov.quantile(.95):.1f}%, "
          f"max {cov.max():.1f}% | sub-C049 {c049:.1f}% (percentile "
          f"{(cov < c049).mean():.0%}, robust z {(c049 - cov.median()) / mad:.1f})")
    top = cov.sort_values(ascending=False).head(5).round(1)
    print(f"    highest: {top.to_dict()}")
