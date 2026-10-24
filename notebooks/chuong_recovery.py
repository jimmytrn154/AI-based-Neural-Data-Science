# %% [markdown]
# # Recovery pilot (Chuong; task C4 in docs/TEAM_TASKS.md; exploratory)
# Baseline: each subject's mean state-probability vector (A-E) in the
# cue-locked -200 to 0 ms window (correct trials). Recovery curve: Euclidean
# distance between the response-locked state-probability vector and that
# baseline, 0 to +1000 ms after the keypress. Recovery measure: mean distance
# over +200 to +1000 ms. Within-subject comparisons: incongruent vs congruent
# correct trials; errors vs correct trials (correct trials subsampled to the
# number of errors, because fewer trials give a larger distance by chance).
# Reliability: split-half (odd vs even trial numbers, separate baselines),
# Spearman-Brown corrected; keep a measure only if >= 0.5 (config).
# No group or clinical analyses.
# Run from the repo root: `python notebooks/chuong_recovery.py`

# %%
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd()))
logging.basicConfig(level=logging.INFO, format="%(message)s")
logging.getLogger("mne").setLevel(logging.ERROR)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr, spearmanr

from src import features, io as pio, microstates as ms

cfg = pio.load_config()
k = cfg["microstates"]["primary_k"]
templates = ms.load_templates(k, cfg)
states = templates.cluster_names
BASELINE_MS = (-200, 0)
CURVE_MS = (0, 1000)
MEASURE_MS = (200, 1000)
RT_MS = (150, 1500)                 # config simon_conditions.keep_trials
MIN_PER_CONDITION = cfg["subjects"]["simon"]["min_trials_per_condition"]  # 40
MIN_ERRORS = 10                     # error measures need at least this many error trials
N_MATCHED = 200                     # random correct-trial subsets matched to the error count
KEEP_IF = cfg["recovery"]["keep_if_reliability_at_least"]
cache = pio.results_dir(cfg, "recovery")  # git-ignored
cache.mkdir(parents=True, exist_ok=True)

qc = pd.read_csv(pio.results_dir(cfg, "templates", "recording_qc.csv"))
ids = sorted(qc[(qc.k == k) & (qc.task == "Simon") & ~qc.excluded].subject_id)
print(f"K = {k}; Simon recordings passing QC: {len(ids)}")

# %% 1. Back-fit cue- and response-locked epochs (cached per subject)
def backfit_lock(sid, lock):
    """Labels of responded, non-bad trials (2-20 Hz, 250 Hz) and their trial rows."""
    epochs = pio.load_epochs(sid, "Simon", lock, cfg)
    trials = pio.match_trials(epochs, pio.load_trials(sid, "Simon", cfg), lock)
    use = (trials["beh_accuracy"].isin([1, 2]) & ~trials["bad_epoch"]).to_numpy()
    inst = pio.epochs_to_mne(epochs)[np.flatnonzero(use)]
    l_freq, h_freq = cfg["preprocessing"]["bandpass_hz"]
    inst.filter(l_freq, h_freq, verbose="error")
    inst.resample(cfg["preprocessing"]["resample_hz"], verbose="error")
    labels = ms.backfit(inst, templates, cfg).labels.astype(np.int8)
    rows = trials.loc[use, ["trial", "beh_congruentcondition", "beh_accuracy", "beh_rt"]]
    return labels, inst.times * 1000, rows.to_numpy()


def cached(sid):
    path = cache / f"{sid}_recovery_labels.npz"
    if not path.exists():
        out = {}
        for lock in ("cue", "resp"):
            out[f"{lock}_labels"], out[f"{lock}_times"], out[f"{lock}_trials"] = backfit_lock(sid, lock)
        np.savez_compressed(path, **out)
    return path


paths = Parallel(n_jobs=8, verbose=2)(delayed(cached)(sid) for sid in ids)
data = {sid: dict(np.load(path)) for sid, path in zip(ids, paths)}
TRIAL, COND, ACC, RT = range(4)

# %% 2. Measures per subject
def correct(trials):
    return (trials[:, ACC] == 1) & (trials[:, RT] >= RT_MS[0]) & (trials[:, RT] <= RT_MS[1])


def errors(trials):
    return (trials[:, ACC] == 2) & (trials[:, RT] >= RT_MS[0]) & (trials[:, RT] <= RT_MS[1])


def baseline(d, half=None):
    t, rows = d["cue_times"], d["cue_trials"]
    use = correct(rows) & (True if half is None else (rows[:, TRIAL] % 2 == half))
    probs = features.state_probabilities(d["cue_labels"][use], len(states))
    window = (t >= BASELINE_MS[0]) & (t < BASELINE_MS[1])
    return np.nanmean(probs[:, window], axis=1)


def curve(d, use, base):
    probs = features.state_probabilities(d["resp_labels"][use], len(states))
    return features.recovery_curve(probs, base)


def measure(d, use, base):
    return features.window_mean(curve(d, use, base), d["resp_times"], MEASURE_MS)


def matched_correct(d, ok, n, base, rng):
    """Mean measure over random subsets of correct trials with n trials each."""
    idx = np.flatnonzero(ok)
    values = []
    for _ in range(N_MATCHED):
        use = np.zeros(len(ok), bool)
        use[rng.choice(idx, n, replace=False)] = True
        values.append(measure(d, use, base))
    return float(np.mean(values))


def subject_measures(sid, d, half=None):
    rows = d["resp_trials"]
    in_half = True if half is None else (rows[:, TRIAL] % 2 == half)
    base = baseline(d, half)
    ok, err = correct(rows) & in_half, errors(rows) & in_half
    con, inc = ok & (rows[:, COND] == 1), ok & (rows[:, COND] == 2)
    min_cond = MIN_PER_CONDITION if half is None else MIN_PER_CONDITION // 2
    min_err = MIN_ERRORS if half is None else MIN_ERRORS // 2
    out = {"simon_recovery_correct": measure(d, ok, base),
           "simon_recovery_n_correct": int(ok.sum()), "simon_recovery_n_con": int(con.sum()),
           "simon_recovery_n_inc": int(inc.sum()), "simon_recovery_n_err": int(err.sum())}
    if con.sum() >= min_cond and inc.sum() >= min_cond:
        out["simon_recovery_con"] = measure(d, con, base)
        out["simon_recovery_inc"] = measure(d, inc, base)
        out["simon_recovery_conflict"] = out["simon_recovery_inc"] - out["simon_recovery_con"]
    if min_err <= err.sum() <= ok.sum():  # matched subsets need at least as many correct trials
        rng = ms.subject_rng(sid + ("" if half is None else f"-half{half}"), cfg)
        out["simon_recovery_err"] = measure(d, err, base)
        out["simon_recovery_correct_matched"] = matched_correct(d, ok, int(err.sum()), base, rng)
        out["simon_recovery_error_effect"] = (out["simon_recovery_err"]
                                              - out["simon_recovery_correct_matched"])
    return out


recovery = pd.DataFrame({sid: subject_measures(sid, d) for sid, d in data.items()}).T
columns = ["simon_recovery_correct", "simon_recovery_con", "simon_recovery_inc",
           "simon_recovery_conflict", "simon_recovery_err", "simon_recovery_correct_matched",
           "simon_recovery_error_effect", "simon_recovery_n_correct", "simon_recovery_n_con",
           "simon_recovery_n_inc", "simon_recovery_n_err"]
recovery = recovery[columns]
recovery.index.name = "subject_id"

n_err = recovery["simon_recovery_n_err"].astype(int)
print(f"Error trials per subject: median {n_err.median():.0f}, IQR {n_err.quantile(.25):.0f}-"
      f"{n_err.quantile(.75):.0f}, range {n_err.min()}-{n_err.max()}; "
      f"{(n_err >= MIN_ERRORS).sum()} of {len(n_err)} subjects have >= {MIN_ERRORS}; "
      f"{recovery['simon_recovery_error_effect'].notna().sum()} have an error measure "
      f"(also needs at least as many correct as error trials)")
few_correct = recovery.index[n_err > recovery["simon_recovery_n_correct"].astype(int)]
print("Subjects with more error than correct trials:",
      {sid: (int(recovery.loc[sid, "simon_recovery_n_err"]), int(recovery.loc[sid, "simon_recovery_n_correct"]))
       for sid in few_correct})
print(f"Subjects with >= {MIN_PER_CONDITION} correct trials per condition: "
      f"{recovery['simon_recovery_conflict'].notna().sum()}")

# Labeled fraction near the end of the response-locked epoch (edge segments)
t = data[ids[0]]["resp_times"]
labeled = np.mean([(d["resp_labels"][correct(d["resp_trials"])] >= 0).mean(axis=0)
                   for d in data.values()], axis=0)
for ms_ in (800, 900, 950, 1000):
    print(f"  labeled fraction at +{ms_} ms: {labeled[np.argmin(np.abs(t - ms_))]:.2f}")

# %% 3. Within-subject comparisons (descriptive)
summary = []
for a, b, label in (("simon_recovery_inc", "simon_recovery_con", "incongruent - congruent"),
                    ("simon_recovery_err", "simon_recovery_correct_matched",
                     "error - matched correct")):
    diff = (recovery[a] - recovery[b]).dropna().astype(float)
    summary.append({"comparison": label, "n": len(diff), "median_diff": diff.median(),
                    "iqr": f"{diff.quantile(.25):.4f} to {diff.quantile(.75):.4f}",
                    "% subjects > 0": f"{(diff > 0).mean():.0%}"})
print("\nWithin-subject differences in the recovery measure (distance units):\n"
      + pd.DataFrame(summary).round(4).to_string(index=False))

# %% 4. Split-half reliability
halves = {h: pd.DataFrame({sid: subject_measures(sid, d, half=h) for sid, d in data.items()}).T
          for h in (0, 1)}
reliability = []
for name in ("simon_recovery_correct", "simon_recovery_conflict", "simon_recovery_error_effect"):
    pair = pd.concat([halves[0][name], halves[1][name]], axis=1, keys=["odd", "even"]).dropna().astype(float)
    r, _ = pearsonr(pair["odd"], pair["even"])
    rho, _ = spearmanr(pair["odd"], pair["even"])
    sb = features.spearman_brown(r)
    reliability.append({"measure": name, "n": len(pair), "split_half_r": r, "spearman_rho": rho,
                        "spearman_brown": sb, "spearman_brown_rank": features.spearman_brown(rho),
                        "keep (>= %.1f)" % KEEP_IF: sb >= KEEP_IF})
reliability = pd.DataFrame(reliability)
print("\nSplit-half reliability (odd vs even trial numbers):\n" + reliability.round(3).to_string(index=False))

# Measures that pass go to recovery.csv; failed ones to a separate "do not use" file
MEASURE_COLUMNS = {
    "simon_recovery_correct": ["simon_recovery_correct", "simon_recovery_n_correct"],
    "simon_recovery_conflict": ["simon_recovery_con", "simon_recovery_inc", "simon_recovery_conflict",
                                "simon_recovery_n_con", "simon_recovery_n_inc"],
    "simon_recovery_error_effect": ["simon_recovery_err", "simon_recovery_correct_matched",
                                    "simon_recovery_error_effect", "simon_recovery_n_err"],
}
passed = reliability.set_index("measure")["spearman_brown"] >= KEEP_IF
for name, keep in (("recovery.csv", True), ("recovery_conflict_unreliable.csv", False)):
    cols = [c for measure, cs in MEASURE_COLUMNS.items() if passed[measure] == keep for c in cs]
    out = pio.results_dir(cfg, "features", name)
    recovery[cols].to_csv(out)
    print("Saved", out, recovery[cols].shape, cols)

# %% 5. Figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4))


def grand(sets):
    curves = np.stack([curve(data[sid], use, baseline(data[sid])) for sid, use in sets])
    return curves.mean(axis=0), curves.std(axis=0, ddof=1) / np.sqrt(len(curves))


def matched_curves(sid, n_subsets=50):
    """Mean curve over random correct-trial subsets of the error-trial count."""
    d = data[sid]
    ok, n = correct(d["resp_trials"]), int(errors(d["resp_trials"]).sum())
    rng = ms.subject_rng(sid + "-figure", cfg)
    base = baseline(d)
    curves = []
    for _ in range(n_subsets):
        use = np.zeros(len(ok), bool)
        use[rng.choice(np.flatnonzero(ok), n, replace=False)] = True
        curves.append(curve(d, use, base))
    return np.mean(curves, axis=0)


def plot(ax, sets, label, color, curves=None):
    if curves is None:
        mean, sem = grand(sets)
    else:
        mean, sem = curves.mean(axis=0), curves.std(axis=0, ddof=1) / np.sqrt(len(curves))
    ax.plot(t, mean, color=color, label=f"{label} (n = {len(sets)})")
    ax.fill_between(t, mean - sem, mean + sem, color=color, alpha=0.25)


conflict_ids = recovery.index[recovery["simon_recovery_conflict"].notna()]
for cond, label, color in ((1, "congruent", "tab:blue"), (2, "incongruent", "tab:red")):
    plot(axes[0], [(sid, correct(data[sid]["resp_trials"]) & (data[sid]["resp_trials"][:, COND] == cond))
                   for sid in conflict_ids], label, color)
error_ids = recovery.index[recovery["simon_recovery_error_effect"].notna()]
plot(axes[1], [(sid, errors(data[sid]["resp_trials"])) for sid in error_ids], "errors", "tab:purple")
plot(axes[1], error_ids, "correct, subsets matched to error count", "tab:green",
     curves=np.stack([matched_curves(sid) for sid in error_ids]))
for ax, title in zip(axes[:2], ("Correct trials by condition", "Errors vs correct")):
    ax.axvspan(*MEASURE_MS, color="grey", alpha=0.1)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlim(-800, 900)  # after ~900 ms most trials are unlabeled (epoch-edge segment)
    ax.set_xlabel("Time from keypress (ms)")
    ax.set_ylabel("Distance from cue baseline")
    ax.set_title(title + " (shaded: measure window)", fontsize=9)
    ax.legend(fontsize=8)
row = reliability.set_index("measure").loc["simon_recovery_correct"]
pair = pd.concat([halves[0]["simon_recovery_correct"], halves[1]["simon_recovery_correct"]],
                 axis=1).dropna().astype(float)
axes[2].scatter(pair.iloc[:, 0], pair.iloc[:, 1], s=12)
axes[2].set_xlabel("Odd trials")
axes[2].set_ylabel("Even trials")
axes[2].set_title(f"Recovery measure, correct trials: split-half r = {row.split_half_r:.2f} "
                  f"(SB {row.spearman_brown:.2f}); rank-based rho = {row.spearman_rho:.2f} "
                  f"(SB {features.spearman_brown(row.spearman_rho):.2f})", fontsize=8)
fig.suptitle(f"Recovery pilot (exploratory), K = {k}, {len(ids)} subjects. Fewer trials give a larger "
             "distance by chance, so errors are compared with correct subsets of the same size", fontsize=9)
fig.tight_layout()
figure = pio.results_dir(cfg, "figures", "chuong_recovery.png")
fig.savefig(figure, dpi=150, bbox_inches="tight")
print("Saved", figure)
