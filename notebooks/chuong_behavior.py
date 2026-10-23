# %% [markdown]
# # Simon behavior (Chuong; tasks C1-C3 in docs/TEAM_TASKS.md)
# 1. Reproduce the Dataset Guide reference values (per-subject means first,
#    then averaged across subjects).
# 2. Simon effect by group, descriptive only.
# 3. Save results/features/behavior.csv (one row per subject_id, no clinical columns).
# Run from the repo root: `python notebooks/chuong_behavior.py`

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))

import pandas as pd

from src import features, io as pio

cfg = pio.load_config()
subjects = pio.load_subjects(cfg)
allowed = cfg["subjects"]["simon"]["allowed_trial_status"]
ids = subjects.loc[subjects["trials_simon"].isin(allowed), "subject_id"].tolist()
print(f"Subjects with trials_simon in {allowed}: {len(ids)}")

# %% Per-subject features
behavior = pd.DataFrame([{"subject_id": sid, **features.simon_behavior(pio.load_trials(sid, "Simon", cfg))}
                         for sid in ids]).set_index("subject_id")
out = pio.results_dir(cfg, "features", "behavior.csv")
behavior.to_csv(out)
print("Saved", out, behavior.shape)
print("Trials per condition:", behavior[["simon_n_con", "simon_n_inc"]].describe().loc[["min", "50%", "max"]].to_dict())

# %% C1: reference values (Dataset Guide)
reference = {"simon_rt_con": 514.7, "simon_rt_inc": 557.4,
             "simon_err_con": 0.064, "simon_err_inc": 0.148}
scale = {"simon_rt_con": 1, "simon_rt_inc": 1, "simon_err_con": 100, "simon_err_inc": 100}
unit = {"simon_rt_con": "ms", "simon_rt_inc": "ms", "simon_err_con": "%", "simon_err_inc": "%"}
rows = [{"measure": key, "ours": behavior[key].mean() * scale[key],
         "reference": ref * scale[key], "unit": unit[key]} for key, ref in reference.items()]
c1 = pd.DataFrame(rows)
c1["difference"] = c1["ours"] - c1["reference"]
print("\nC1, mean across subjects of per-subject values:\n" + c1.round(2).to_string(index=False))

# %% C2: Simon effect by group (descriptive)
guide = {"PD": {"simon_effect_rt": 50, "simon_err_con": 7.4, "simon_err_inc": 17.6},
         "Control": {"simon_effect_rt": 34, "simon_err_con": 4.5, "simon_err_inc": 8.4}}
grouped = behavior.join(subjects.set_index("subject_id")["group"])
rows = []
for group, df in grouped.groupby("group"):
    for key, s in (("simon_effect_rt", 1), ("simon_err_con", 100), ("simon_err_inc", 100),
                   ("simon_effect_err", 100), ("simon_rt_con", 1), ("simon_rt_inc", 1)):
        rows.append({"group": group, "n": len(df), "measure": key, "mean": df[key].mean() * s,
                     "sd": df[key].std() * s, "guide": guide[group].get(key)})
c2 = pd.DataFrame(rows)
print("\nC2, by group (RT in ms, errors in %):\n" + c2.round(2).to_string(index=False))
