# %% [markdown]
# # Template labels and figures (Chuong)
# Applies the visually confirmed labels from config.yaml
# (microstates.confirmed_labels) to the saved templates and redraws the
# topomaps with their label status. Run after chuong_final_templates.py and
# chuong_task_specific.py, from the repo root.

# %%
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))

import matplotlib
matplotlib.use("Agg")

from src import io as pio, microstates as ms

cfg = pio.load_config()
band = "-".join(f"{f:g}" for f in cfg["preprocessing"]["bandpass_hz"])
with open(pio.results_dir(cfg, "templates", "summary.json")) as f:
    summary = json.load(f)


def redraw(cluster, title, name):
    fig = cluster.plot(show=False)
    fig.suptitle(title)
    out = pio.results_dir(cfg, "figures", name)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("Saved", out, cluster.cluster_names)


# %% Shared templates, every K
for k in cfg["microstates"]["k_range"]:
    templates = ms.load_templates(k, cfg)
    if ms.apply_confirmed_labels(templates, cfg):
        ms.save_templates(templates, cfg)
    info = summary["K"][str(k)]
    note = " - C includes canonical D" if k == 4 else ""
    redraw(templates,
           f"Templates K = {k}, {band} Hz, {info['n_fit']} subjects "
           f"(GEV {info['template_gev_subjects_mean']:.1%}) - {ms.label_status(k, cfg)}{note}",
           f"chuong_templates_K{k}.png")

# %% Task-specific templates (labeled like the shared templates of their K)
for path in sorted(pio.results_dir(cfg, "templates").glob("templates_K*_[rs]*.fif")):
    k, task = path.stem.removeprefix("templates_K").split("_")
    if task not in ("rest", "simon"):
        continue
    templates = ms.load_templates(int(k), cfg, f"_{task}")
    if ms.apply_confirmed_labels(templates, cfg):
        ms.save_templates(templates, cfg, f"_{task}")
    status = "primary" if int(k) == cfg["microstates"]["primary_k"] else "sensitivity"
    redraw(templates, f"{task.capitalize()}-only templates K = {k} "
           f"(matched to shared K = {k}; {status})", f"chuong_templates_K{k}_{task}.png")
