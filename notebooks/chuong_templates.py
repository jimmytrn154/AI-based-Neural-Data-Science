# %% [markdown]
# # Draft microstate templates (Chuong)
# Two-level modified k-means on pooled Rest + Simon GFP peaks.
# Run from the repo root: `python notebooks/chuong_templates.py`
# (or cell by cell in VS Code / Jupyter with the `# %%` markers).

# %%
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import matplotlib
matplotlib.use("Agg")

from src import io as pio, microstates as ms

cfg = pio.load_config()
draft = cfg["microstates"]["draft_templates"]
k = draft["k"]
band = "-".join(f"{f:g}" for f in cfg["preprocessing"]["bandpass_hz"])

# PROVISIONAL state labels by visual inspection of the draft maps
# (A = right-frontal/left-posterior, B = left-frontal/right-posterior,
# C = anterior-posterior, D = fronto-central; polarity ignored).
# PROVISIONAL_ORDER[i] = fitted map index that becomes label i.
# Set to None to save unlabeled maps for inspection.
# Final labels are fixed on the final templates.
PROVISIONAL_ORDER = [0, 1, 3, 2]  # Oct 22 draft (2-20 Hz, P048 replaced by P009)
LABELS = ["A", "B", "C", "D"]

# %% Draft subjects (stratified sample; ms.DRAFT_EXCLUDE replaced, seeded)
ids, replaced = ms.draft_subjects(cfg)
for old, new in replaced.items():
    print(f"Excluded {old} ({ms.DRAFT_EXCLUDE[old]}); replaced by {new}")
print("Draft subjects:", ids)

# %% Two-level fit
group, individual, table = ms.fit_templates(ids, k, cfg)
print(table.to_string(index=False))
print(f"Group-level GEV (on individual maps): {group.GEV_:.3f}")
print(f"Template GEV on subjects' GFP peaks: {table.gev_template.mean():.3f} "
      f"± {table.gev_template.std():.3f}")

# %% Topomaps, labeled if the provisional order is set
title = (f"Draft templates, K = {k}, {band} Hz, {len(ids)} subjects "
         f"(GEV {table.gev_template.mean():.1%})")
if PROVISIONAL_ORDER is None:
    name = f"chuong_templates_K{k}_draft_unlabeled.png"
else:
    ms.label_templates(group, PROVISIONAL_ORDER, LABELS)
    title += " — labels PROVISIONAL"
    name = f"chuong_templates_K{k}_draft.png"
fig = group.plot(show=False)
fig.suptitle(title)
out = pio.results_dir(cfg, "figures", name)
fig.savefig(out, dpi=150, bbox_inches="tight")
print("Saved", out)

if PROVISIONAL_ORDER is not None:
    print("Saved", ms.save_templates(group, cfg))
