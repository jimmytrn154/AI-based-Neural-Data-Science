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

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
draft = cfg["microstates"]["draft_templates"]

# %% Subjects passing both the Rest and Simon rules
rest, _ = pp.select_subjects("Rest", cfg)
simon, _ = pp.select_subjects("Simon", cfg)
both = pp.select_rest_vs_simon(rest, simon, cfg)
pool = pio.load_subjects(cfg).set_index("subject_id").loc[both].reset_index()
ids = ms.stratified_sample(pool, draft["n_subjects"], cfg)
print("Draft subjects:", ids)
print(pool.set_index("subject_id").loc[ids, "group"].value_counts().to_dict())

# %% Two-level fit
group, individual, table = ms.fit_templates(ids, draft["k"], cfg)
path = ms.save_templates(group, cfg)
print("Saved", path)
print(table.to_string(index=False))
print(f"Group-level GEV (on individual maps): {group.GEV_:.3f}")
print(f"Template GEV on subjects' GFP peaks: {table.gev_template.mean():.3f} "
      f"± {table.gev_template.std():.3f}")

# %% Topomaps
fig = group.plot(show=False)
fig.suptitle(f"Draft templates, K = {draft['k']}, {len(ids)} subjects "
             f"(GEV {table.gev_template.mean():.1%})")
out = pio.results_dir(cfg, "figures", f"chuong_templates_K{draft['k']}_draft.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
print("Saved", out)
