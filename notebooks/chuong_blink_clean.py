# %% [markdown]
# # P2.7 blink cleaning (Chuong)
# Rest recordings with no eye component removed (nICremoved == 0): mark blink
# windows as BAD_blink (config phase2.blink_cleaning), exclude them from
# back-fitting, and recompute the Rest metrics with the primary K templates.
# Other Rest recordings keep their primary values.
# Output: results/features/microstates_rest_blinkclean.csv (Rest columns, one
# row per subject_id; rest_blink_bad_frac = fraction of Rest time marked).
# Run from the repo root: `python notebooks/chuong_blink_clean.py`

# %%
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
k = cfg["microstates"]["primary_k"]
templates = ms.load_templates(k, cfg)
primary = pd.read_csv(pio.results_dir(cfg, "features", "microstates_task.csv")).set_index("subject_id")
rest_cols = [c for c in primary.columns if c.startswith("rest_")]

rest, _ = pp.select_subjects("Rest", cfg)
flagged = rest.loc[rest.blink_flag, "subject_id"]
flagged = [sid for sid in flagged if primary.loc[sid, "rest_gev_total"] == primary.loc[sid, "rest_gev_total"]]  # passed QC
print(f"Rest recordings with no eye component removed and passing QC: {len(flagged)}")


# %% Clean and recompute
def clean(sid):
    raw = pp.preprocess_raw(pio.load_continuous(sid, "Rest", cfg), cfg)
    bad_frac = pp.annotate_blinks(raw, cfg)
    segmentation = ms.backfit(raw, templates, cfg)  # BAD segments excluded (pycrostates path)
    row = ms.segmentation_metrics(segmentation, "rest")
    row["rest_blink_bad_frac"] = bad_frac
    return sid, row


rows = dict(Parallel(n_jobs=8)(delayed(clean)(sid) for sid in flagged))
cleaned = pd.DataFrame.from_dict(rows, orient="index")
out = primary[rest_cols].dropna(how="all").copy()
out["rest_blink_bad_frac"] = 0.0
out.loc[cleaned.index, cleaned.columns] = cleaned
out.index.name = "subject_id"
path = pio.results_dir(cfg, "features", "microstates_rest_blinkclean.csv")
out.to_csv(path)
print("Saved", path, out.shape)
frac = cleaned["rest_blink_bad_frac"]
print(f"Time marked BAD_blink in the {len(frac)} cleaned recordings: median {frac.median():.1%}, "
      f"IQR {frac.quantile(.25):.1%}-{frac.quantile(.75):.1%}, max {frac.max():.1%}; "
      f"recordings with any marked time: {(frac > 0).sum()}")
print("Unlabeled fraction after cleaning (max):", round(float(cleaned['rest_unlabeled_frac'].max()), 3))
