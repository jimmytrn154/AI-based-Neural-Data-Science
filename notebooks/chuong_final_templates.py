# %% [markdown]
# # Final microstate templates, full cohort (Chuong)
# 1. GFP peaks of every eligible recording (cached in results/peaks/).
# 2. For each K in k_range: individual fits -> group fit (pass 1) -> recording
#    QC (template GEV robust z per task) -> group refit without flagged
#    subjects (pass 2) -> labels by correlation with the approved draft A-D.
# 3. Primary K by the config rule (cluster scores on the group-level fits).
# 4. Task-specific templates (Rest only, Simon only) at the primary K.
# Labels are assigned automatically and are PENDING VISUAL CONFIRMATION.
#
# Run from the repo root: `python notebooks/chuong_final_templates.py`
# Dry run (12 subjects, n_init 10, output to MS_DRYRUN_DIR): MS_DRYRUN=1

# %%
import json
import logging
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pycrostates.io import ChData

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
DRY = os.environ.get("MS_DRYRUN") == "1"
if DRY:
    cfg["paths"]["results"] = os.environ["MS_DRYRUN_DIR"]
    cfg["microstates"]["n_init"] = 10
k_range = cfg["microstates"]["k_range"]
band = "-".join(f"{f:g}" for f in cfg["preprocessing"]["bandpass_hz"])
draft = ms.load_templates(4, pio.load_config(), "_draft")  # approved draft, repo results/
SCORE_SAMPLE = 5000  # pooled GFP peaks for the secondary K-selection check

# %% Eligible recordings
rest, rest_counts = pp.select_subjects("Rest", cfg)
simon, simon_counts = pp.select_subjects("Simon", cfg)
both = pp.select_rest_vs_simon(rest, simon, cfg)
subjects = pio.load_subjects(cfg).set_index("subject_id")
eligible = {"Rest": list(rest.subject_id), "Simon": list(simon.subject_id)}
if DRY:
    both = ms.stratified_sample(subjects.loc[both].reset_index(), 12, cfg)
    rest_only = sorted(set(eligible["Rest"]) - set(eligible["Simon"]))[:1]
    eligible = {"Rest": both + rest_only, "Simon": both}
print({task: len(ids) for task, ids in eligible.items()}, "| both:", len(both))

# %% 1. GFP peaks per recording (cached)
peaks_dir = pio.results_dir(cfg, "peaks")
peaks_dir.mkdir(parents=True, exist_ok=True)


def cached_peaks(sid, task):
    path = peaks_dir / f"{sid}_task-{task}_peaks.npy"
    if not path.exists():
        np.save(path, ms.recording_peaks(sid, task, cfg).get_data())
    return path


jobs = [(sid, task) for task, ids in eligible.items() for sid in ids]
paths = Parallel(n_jobs=8, verbose=2)(delayed(cached_peaks)(sid, task) for sid, task in jobs)
info = ms.recording_peaks(both[0], "Rest", cfg).info
peaks = {job: ChData(np.load(path), info) for job, path in zip(jobs, paths)}
n_peaks = pd.Series({job: p.get_data().shape[1] for job, p in peaks.items()})
print("Peaks per recording (median):", n_peaks.groupby(level=1).median().to_dict())

# Equal-size samples per subject for the shared templates
sampled, n_sampled = {}, {}
for sid in both:
    sampled[sid], n_sampled[sid] = ms.equal_task_sample(
        {task: peaks[(sid, task)] for task in ("Rest", "Simon")}, sid, cfg)
pooled = {sid: ms.pool(sampled[sid].values()) for sid in both}
print("Peaks sampled per task: median", int(np.median(list(n_sampled.values()))),
      "min", min(n_sampled.values()))

# Fixed sample of pooled peaks for the secondary scores (same for every K)
rng = np.random.default_rng(cfg["random_seed"])
all_pooled = ms.pool(pooled[sid] for sid in both).get_data()
score_data = ChData(all_pooled[:, rng.choice(all_pooled.shape[1],
                                             min(SCORE_SAMPLE, all_pooled.shape[1]),
                                             replace=False)], info)
del all_pooled

# %% 2. Templates for each K
def save_figure(cluster, title, name):
    fig = cluster.plot(show=False)
    fig.suptitle(title)
    out = pio.results_dir(cfg, "figures", name)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    return out


summary = {"n_both": len(both), "n_eligible": {t: len(v) for t, v in eligible.items()},
           "K": {}}
scores, scores_pooled, qc_all = {}, {}, []
for k in k_range:
    print(f"\n==================== K = {k} ====================")
    individual = {sid: ms.fit_individual(pooled[sid], k, cfg) for sid in both}

    # Pass 1 and recording QC
    pass1 = ms.fit_group([individual[sid] for sid in both], k, cfg)
    qc = pd.DataFrame([{"subject_id": sid, "task": task,
                        "gev_template": ms.template_gev(pass1, peaks[(sid, task)])}
                       for task, ids in eligible.items() for sid in ids])
    qc = ms.flag_recordings(qc, cfg)
    qc = qc.join(subjects[["group", "moca", "age", "gender"]], on="subject_id")
    qc.insert(0, "k", k)
    drop = sorted(set(qc.loc[qc.excluded, "subject_id"]) & set(both))
    fit_ids = [sid for sid in both if sid not in drop]
    print(f"QC: {int(qc.excluded.sum())} recordings flagged; "
          f"{len(drop)} subjects dropped from the refit; {len(fit_ids)} used")
    if qc.excluded.any():
        print(qc.loc[qc.excluded, ["subject_id", "task", "group", "moca", "gev_template",
                                   "gev_z"]].round(3).to_string(index=False))

    # Pass 2 (refit once without the flagged subjects, then stop)
    final = ms.fit_group([individual[sid] for sid in fit_ids], k, cfg) if drop else pass1.copy()
    qc["gev_template_final"] = [ms.template_gev(final, peaks[(sid, task)])
                                for sid, task in zip(qc.subject_id, qc.task)]
    qc_all.append(qc)

    # Labels by correlation with the approved draft (pending visual confirmation)
    fit_peaks = ms.pool(pooled[sid] for sid in fit_ids)
    order, names, corr, flags = ms.match_labels(final, draft, fit_peaks)
    ms.label_templates(final, order, names)
    p_order, p_names, _, _ = ms.match_labels(pass1, draft, fit_peaks)
    ms.label_templates(pass1, p_order, p_names)
    explained = ms.state_gev(final, fit_peaks)
    del fit_peaks
    print("Correlation with draft (|r|):\n" + corr.round(2).to_string())
    print("Flags:", flags or "none")

    ms.save_templates(pass1, cfg, "_pass1")
    path = ms.save_templates(final, cfg)
    gev_subject = pd.Series({sid: ms.template_gev(final, pooled[sid]) for sid in fit_ids})
    title = (f"Templates K = {k}, {band} Hz, {len(fit_ids)} subjects "
             f"(GEV {gev_subject.mean():.1%}) - labels pending visual confirmation")
    figure = save_figure(final, title, f"chuong_templates_K{k}.png")
    print("Saved", path, figure)

    scores[k] = ms.cluster_scores(final)
    scores_pooled[k] = ms.cluster_scores(final, score_data)
    summary["K"][k] = {
        "n_fit": len(fit_ids), "dropped": drop,
        "group_gev_individual_maps": final.GEV_,
        "template_gev_subjects_mean": gev_subject.mean(),
        "template_gev_subjects_sd": gev_subject.std(),
        "state_gev": dict(zip(final.cluster_names, explained.round(4).tolist())),
        "corr_with_draft": corr.round(3).to_dict(orient="index"),
        "label_flags": flags,
    }

qc_all = pd.concat(qc_all, ignore_index=True)
qc_all.to_csv(pio.results_dir(cfg, "templates", "recording_qc.csv"), index=False)

# %% 3. Primary K
scores = pd.DataFrame(scores).T
scores_pooled = pd.DataFrame(scores_pooled).T
primary_k, ranks = ms.choose_k(scores)
check_k, ranks_pooled = ms.choose_k(scores_pooled)
print("\nCluster scores (group-level fits):\n" + scores.round(4).to_string())
print("Ranks:\n" + ranks.to_string())
print(f"Primary K (config rule): {primary_k}")
print(f"\nCheck on {score_data.get_data().shape[1]} pooled GFP peaks:\n"
      + scores_pooled.round(4).to_string())
print("Ranks:\n" + ranks_pooled.to_string())
print(f"K by the same rule on pooled peaks: {check_k}")
pd.concat({"group_level": scores, "pooled_peaks": scores_pooled}, axis=1).to_csv(
    pio.results_dir(cfg, "templates", "k_selection.csv"))
summary.update({"primary_k": primary_k, "check_k_pooled_peaks": check_k,
                "scores": scores.round(5).to_dict(orient="index"),
                "ranks": ranks.to_dict(orient="index"),
                "scores_pooled": scores_pooled.round(5).to_dict(orient="index"),
                "ranks_pooled": ranks_pooled.to_dict(orient="index")})

# %% 4. Task-specific templates at the primary K
shared = ms.load_templates(primary_k, cfg)
fit_ids = [sid for sid in both if sid not in summary["K"][primary_k]["dropped"]]
summary["task_specific"] = {}
for task in ("Rest", "Simon"):
    individual = [ms.fit_individual(sampled[sid][task], primary_k, cfg) for sid in fit_ids]
    specific = ms.fit_group(individual, primary_k, cfg)
    task_peaks = ms.pool(sampled[sid][task] for sid in fit_ids)
    order, names, corr, flags = ms.match_labels(specific, shared, task_peaks)
    ms.label_templates(specific, order, names)
    del task_peaks
    ms.save_templates(specific, cfg, f"_{task.lower()}")
    save_figure(specific, f"{task}-only templates K = {primary_k} (matched to shared)",
                f"chuong_templates_K{primary_k}_{task.lower()}.png")
    diag = pd.Series(np.diag(corr.to_numpy()), index=corr.index)
    print(f"\n{task}-only vs shared (|r|):\n" + corr.round(2).to_string())
    summary["task_specific"][task] = {"corr": corr.round(3).to_dict(orient="index"),
                                      "matched_r": diag.round(3).to_dict(),
                                      "flags": flags}

with open(pio.results_dir(cfg, "templates", "summary.json"), "w") as f:
    json.dump(summary, f, indent=1, default=str)
print("\nDone. Primary K =", primary_k)
