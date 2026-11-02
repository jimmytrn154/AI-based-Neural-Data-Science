r"""S4: single-trial Simon microstates and condition-wise state probabilities.

Keep this file AND seonwoo_erp_segmentation.py in notebooks/. Run from repo root:
    .\.venv\Scripts\python.exe -X utf8 notebooks/seonwoo_state_probabilities.py

Optional subset check (requires completed S2/S3 subset outputs):
    python -X utf8 notebooks/seonwoo_state_probabilities.py --subjects sub-P001 sub-C004

Uses the same cohort and trial counts as S3. Each retained trial is independently
filtered over the FULL epoch, resampled, and back-fitted before cropping.
Default batch size is 32 trials; --batch-size 16 reduces temporary memory use.
Completed subjects are cached with input/settings/code/version fingerprints.
Rerun the same command to resume. --recompute forces new back-fitting.

Primary probabilities use src.features.state_probabilities, as in Chuong's
recovery code: denominator = labeled trials at EACH time point, excluding -1.
The all-retained-trial denominator is saved too, together with labeled counts.
Zero labeled trials -> NaN, never zero. No probability baseline subtraction.
No averaging of EEG before back-fitting; no trials pooled across subjects.
This step does not run the S5 permutation tests or produce S6 conflict.csv.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from zipfile import BadZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src import features, io as pio, microstates as ms
from seonwoo_erp_segmentation import choose_subjects, require_epoch_files

CONDITION_NAMES = ("Congruent", "Incongruent")
CONDITION_CODES = (1, 2)


def digest_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hash_object(value):
    text = json.dumps(value, sort_keys=True, ensure_ascii=True)
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def pipeline_signature(cfg, k):
    sources = [
        Path(__file__), Path(pio.__file__), Path(ms.__file__), Path(features.__file__),
        Path(__file__).with_name("seonwoo_erp_segmentation.py"),
    ]
    return hash_object({
        "schema": "seonwoo_s4_v1",
        "config": cfg,
        "template_sha256": digest_file(ms.templates_path(k, cfg)),
        "code": {p.name: digest_file(p) for p in sources},
        "versions": {p: importlib.metadata.version(p)
                     for p in ("numpy", "scipy", "mne", "pycrostates", "h5py")},
    })


def subject_signature(sid, cfg, signature, counts):
    files = []
    for path in (pio.epochs_path(sid, "Simon", "cue", cfg),
                 pio.trials_path(sid, "Simon", cfg)):
        stat = path.stat()
        files.append({"path": str(path.resolve()), "size": stat.st_size,
                      "mtime_ns": stat.st_mtime_ns})
    return hash_object({"pipeline": signature, "files": files, "counts": counts})


def save_npz_atomic(path, **arrays):
    """An interruption cannot turn a partial write into a completed cache."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".npz.tmp",
                                         delete=False) as file:
            temporary = Path(file.name)
            np.savez_compressed(file, **arrays)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def probabilities_by_condition(labels, condition, k):
    probabilities, all_trial_probabilities, labeled_counts, trial_counts = [], [], [], []
    for code in CONDITION_CODES:
        lab = labels[condition == code]
        if len(lab) == 0:
            raise ValueError("No retained trials in one condition.")
        n_labeled = (lab >= 0).sum(axis=0)
        prob = features.state_probabilities(lab, k)
        all_trials = np.stack([(lab == s).mean(axis=0) for s in range(k)])
        valid = n_labeled > 0
        if not np.allclose(prob[:, valid].sum(axis=0), 1):
            raise ValueError("Conditional state probabilities do not sum to one.")
        if not np.allclose(all_trials.sum(axis=0), n_labeled / len(lab)):
            raise ValueError("All-trial state fractions do not match the labeled fraction.")
        probabilities.append(prob)
        all_trial_probabilities.append(all_trials)
        labeled_counts.append(n_labeled)
        trial_counts.append(len(lab))
    return (np.stack(probabilities), np.stack(all_trial_probabilities),
            np.stack(labeled_counts), np.asarray(trial_counts))


def process_subject(sid, cfg, templates, expected_counts, batch_size):
    epochs = pio.load_epochs(sid, "Simon", "cue", cfg)
    trials = pio.match_trials(epochs, pio.load_trials(sid, "Simon", cfg), "cue")
    keep = (
        trials["beh_accuracy"].eq(1)
        & ~trials["bad_epoch"]
        & trials["beh_rt"].between(150, 1500)
    ).to_numpy()
    rows = np.flatnonzero(keep)
    condition = trials.loc[keep, "beh_congruentcondition"].to_numpy()
    if not np.isin(condition, CONDITION_CODES).all():
        raise ValueError(f"{sid}: unknown condition in retained trials.")
    counts = [int((condition == code).sum()) for code in CONDITION_CODES]
    if counts != expected_counts:
        raise ValueError(f"{sid}: counts changed since S3: {counts} vs {expected_counts}.")
    if min(counts) < int(cfg["subjects"]["simon"]["min_trials_per_condition"]):
        raise ValueError(f"{sid}: insufficient trials under the current config.")
    baseline = np.asarray(epochs["erp_baseline"]).ravel()
    if baseline.shape != (2,) or not np.allclose(baseline, [-200, 0]):
        raise ValueError(f"{sid}: unexpected supplied baseline.")
    if not np.allclose(np.diff(epochs["times"]), 1000.0 / epochs["srate"]):
        raise ValueError(f"{sid}: sample rate and time axis do not match.")
    names = list(templates.info["ch_names"])
    if set(epochs["ch_names"]) != set(names):
        raise ValueError(f"{sid}: channels differ from the templates.")

    montage = ("colin27_1005" if "colin27_1005" in mne.channels.get_builtin_montages()
               else "standard_1005")
    batches, full_times = [], None
    for start in range(0, len(rows), batch_size):
        take = rows[start:start + batch_size]
        chunk = epochs.copy()
        chunk["data"] = epochs["data"][take]
        chunk["bad_epoch"] = epochs["bad_epoch"][take]
        if not np.isfinite(chunk["data"]).all():
            raise ValueError(f"{sid}: non-finite EEG in retained trials.")
        inst = pio.epochs_to_mne(chunk, montage=montage)
        del chunk
        inst.reorder_channels(names)
        inst.filter(*cfg["preprocessing"]["bandpass_hz"], verbose="error")
        inst.resample(cfg["preprocessing"]["resample_hz"], verbose="error")
        lab = ms.backfit(inst, templates, cfg).labels
        times = np.round(inst.times * 1000, 6)
        if full_times is None:
            full_times = times
        elif not np.array_equal(times, full_times):
            raise ValueError(f"{sid}: time grid changed between batches.")
        batches.append(lab.astype(np.int16))
        print(f"  trials {start + len(take)}/{len(rows)}", flush=True)
        del inst
    labels = np.concatenate(batches, axis=0)
    k = len(templates.cluster_names)
    if not np.isin(labels, np.arange(-1, k)).all():
        raise ValueError(f"{sid}: invalid state label.")
    lo, hi = cfg["simon_conditions"]["analysis_window_ms"]
    window = (full_times >= lo) & (full_times <= hi)
    if not window.any():
        raise ValueError("Empty configured analysis window.")
    p, p_all, labeled, n_trials = probabilities_by_condition(
        labels[:, window], condition, k
    )
    return {
        "subject_id": np.asarray(sid),
        "state_names": np.asarray(templates.cluster_names),
        "condition_names": np.asarray(CONDITION_NAMES),
        "condition_codes": np.asarray(CONDITION_CODES),
        "trial_row_zero_based": rows,
        "trial_condition": condition.astype(np.int16),
        "trial_rt_ms": trials.loc[keep, "beh_rt"].to_numpy(dtype=float),
        "times_ms_full": full_times,
        "labels_full": labels,
        "times_ms": full_times[window],
        "probabilities": p,
        "probabilities_all_trials": p_all,
        "n_labeled": labeled,
        "n_trials": n_trials,
    }


def check_payload(data, sid, k, expected_counts):
    if str(data["subject_id"].item()) != sid:
        raise ValueError("Cache subject ID mismatch.")
    if data["n_trials"].tolist() != expected_counts:
        raise ValueError("Cache trial-count mismatch.")
    n_times = len(data["times_ms"])
    if data["probabilities"].shape != (2, k, n_times):
        raise ValueError("Cache probability shape mismatch.")
    if data["labels_full"].shape != (sum(expected_counts), len(data["times_ms_full"])):
        raise ValueError("Cache label shape mismatch.")
    if not np.isin(data["labels_full"], np.arange(-1, k)).all():
        raise ValueError("Invalid cached labels.")
    p = data["probabilities"]
    labeled = data["n_labeled"]
    for c in range(2):
        valid = labeled[c] > 0
        if not np.allclose(p[c][:, valid].sum(axis=0), 1):
            raise ValueError("Cached probabilities do not sum to one.")
        if not np.isnan(p[c][:, ~valid]).all():
            raise ValueError("Zero-denominator probabilities must be NaN.")
    finite = p[np.isfinite(p)]
    if np.any((finite < 0) | (finite > 1)):
        raise ValueError("Cached probabilities outside [0, 1].")


def mean_sem(data):
    """Equal subject weights; report valid subject counts, without filling NaNs."""
    valid = np.isfinite(data)
    n = valid.sum(axis=0)
    shape = data.shape[1:]
    mean = np.divide(np.where(valid, data, 0).sum(axis=0), n,
                     out=np.full(shape, np.nan), where=n > 0)
    residuals = np.where(valid, data - mean, 0)
    variance = np.divide((residuals ** 2).sum(axis=0), n - 1,
                         out=np.full(shape, np.nan), where=n > 1)
    sem = np.sqrt(np.divide(variance, n, out=np.full(shape, np.nan), where=n > 1))
    return mean, sem, n


def make_figure(times, probs, states, path, subset):
    means, sem, n_valid = mean_sem(probs)
    n = len(probs)
    fig, axes = plt.subplots(len(states), 1, figsize=(10, 10), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    upper = means + np.nan_to_num(sem, nan=0)
    limit = min(1.0, max(0.2, float(np.nanmax(upper)) * 1.12))
    for s, state in enumerate(states):
        ax = axes[s]
        for c, (label, color) in enumerate([
            ("Congruent", "#2463A6"), ("Incongruent", "#D65C3A")
        ]):
            ax.plot(times, means[c, s], color=color, lw=1.5, label=label)
            ax.fill_between(times, means[c, s] - sem[c, s],
                            means[c, s] + sem[c, s], color=color, alpha=0.16)
        ax.axvline(0, color="0.5", ls="--", lw=0.8)
        ax.set_ylim(0, limit)
        ax.set_xlim(times[0], times[-1])
        ax.set_ylabel(f"P({state})")
        ax.grid(axis="y", alpha=0.15)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel("Time from cue onset (ms)")
    mode = "SUBSET CHECK" if subset else "single-trial state probabilities"
    fig.suptitle(f"Simon conditions | N = {n} | {mode}", fontsize=14, y=0.97)
    fig.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.12, hspace=0.18)
    fig.text(
        0.5, 0.025,
        "Participant means ±1 SEM; denominator = labeled trials at each time point.\n"
        f"Valid subjects per condition/time: {int(n_valid.min())}–{int(n_valid.max())}. "
        "No significance test performed.",
        ha="center", fontsize=9, color="#555555",
    )
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", help="Subset check only.")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--recompute", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1:
        parser.error("--batch-size must be positive.")
    suffix = "_subset" if args.subjects else ""
    cfg = pio.load_config()
    k = int(cfg["microstates"]["primary_k"])
    templates = ms.load_templates(k, cfg)
    states = templates.cluster_names
    if states != cfg["microstates"]["confirmed_labels"].get(k):
        raise ValueError("Template labels differ from confirmed config labels.")
    ids, _ = choose_subjects(cfg, states, args.subjects, suffix)
    s3_path = pio.results_dir(cfg, "features", f"seonwoo_erp_segmentation{suffix}_qc.csv")
    s3 = pd.read_csv(s3_path).set_index("subject_id")
    finished = {
        "included", "excluded_subject_rule", "excluded_trial_status",
        "excluded_recording_qc", "excluded_min_trials", "excluded_microstate_qc",
    }
    if s3.index.duplicated().any() or not s3["status"].isin(finished).all():
        raise ValueError("Complete S3 before running S4; its audit has unfinished/error rows.")
    if set(s3.index[s3["status"].eq("included")]) != set(ids):
        raise ValueError("S3 cohort differs from current rules; check S2/S3 inputs.")

    stem = f"seonwoo_state_probabilities{suffix}"
    figure_dir = pio.results_dir(cfg, "figures")
    data_dir = pio.results_dir(cfg, "erp")
    cache_dir = pio.results_dir(cfg, "labels", f"seonwoo_cue_K{k}")
    audit_dir = pio.results_dir(cfg, "features")
    for directory in (figure_dir, data_dir, cache_dir, audit_dir):
        directory.mkdir(parents=True, exist_ok=True)
    require_epoch_files(ids, cfg, audit_dir / f"{stem}_missing_files.csv")
    signature = pipeline_signature(cfg, k)
    probabilities, all_trial_probabilities, denominators, trial_counts, rows = [], [], [], [], []
    reference_times = None
    print(f"S4: {len(ids)} subjects; states {states}; batch size {args.batch_size}", flush=True)
    for i, sid in enumerate(ids, start=1):
        expected = [int(s3.loc[sid, "n_congruent"]), int(s3.loc[sid, "n_incongruent"])]
        fingerprint = subject_signature(sid, cfg, signature, expected)
        cache_path = cache_dir / f"{sid}.npz"
        payload, reused = None, False
        if cache_path.is_file() and not args.recompute:
            try:
                with np.load(cache_path, allow_pickle=False) as old:
                    if str(old["fingerprint"].item()) == fingerprint:
                        candidate = {key: old[key] for key in old.files}
                        check_payload(candidate, sid, k, expected)
                        payload, reused = candidate, True
            except (OSError, ValueError, KeyError, EOFError, BadZipFile):
                print(f"[{i}/{len(ids)}] {sid}: invalid cache; recomputing.", flush=True)
        if payload is None:
            print(f"[{i}/{len(ids)}] {sid}: computing...", flush=True)
            payload = process_subject(sid, cfg, templates, expected, args.batch_size)
            payload["fingerprint"] = np.asarray(fingerprint)
            check_payload(payload, sid, k, expected)
            save_npz_atomic(cache_path, **payload)
        else:
            print(f"[{i}/{len(ids)}] {sid}: cached result reused.", flush=True)
        t = payload["times_ms"]
        if reference_times is None:
            reference_times = t.copy()
        elif not np.array_equal(reference_times, t):
            raise ValueError(f"{sid}: time grid differs across subjects.")
        probabilities.append(payload["probabilities"])
        all_trial_probabilities.append(payload["probabilities_all_trials"])
        denominators.append(payload["n_labeled"])
        trial_counts.append(payload["n_trials"])
        row = {"subject_id": sid, "status": "included", "cache_reused": reused}
        for c, name in enumerate(("con", "inc")):
            nl = payload["n_labeled"][c]
            nt = int(payload["n_trials"][c])
            row[f"simon_n_{name}"] = nt
            row[f"simon_min_labeled_{name}"] = int(nl.min())
            row[f"simon_mean_labeled_fraction_{name}"] = float((nl / nt).mean())
            row[f"simon_zero_labeled_timepoints_{name}"] = int((nl == 0).sum())
        rows.append(row)
        pd.DataFrame(rows).to_csv(audit_dir / f"{stem}_qc.csv", index=False)

    probs = np.stack(probabilities)
    all_trials = np.stack(all_trial_probabilities)
    n_labeled = np.stack(denominators)
    n_trials = np.stack(trial_counts)
    metadata = {
        "schema": "seonwoo_s4_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "subset_check" if args.subjects else "full_cohort",
        "pipeline_signature": signature,
        "axes": ["subject", "condition", "state", "time"],
        "probability_denominator": "labeled_trials_at_each_timepoint",
        "alternate_denominator": "all_retained_trials",
        "baseline_subtracted": False,
        "analysis_window_ms": cfg["simon_conditions"]["analysis_window_ms"],
        "backfit": ms.backfit_options(cfg),
        "bandpass_hz": cfg["preprocessing"]["bandpass_hz"],
        "resample_hz": cfg["preprocessing"]["resample_hz"],
    }
    data_path = data_dir / f"{stem}.npz"
    save_npz_atomic(
        data_path, subject_id=np.asarray(ids), condition_names=np.asarray(CONDITION_NAMES),
        condition_codes=np.asarray(CONDITION_CODES), state_names=np.asarray(states),
        times_ms=reference_times, probabilities=probs, probabilities_all_trials=all_trials,
        n_labeled=n_labeled, n_trials=n_trials,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )
    mean, sem, n_valid = mean_sem(probs)
    diff_mean, diff_sem, diff_n = mean_sem(probs[:, 1] - probs[:, 0])
    curves = []
    for s, state in enumerate(states):
        curves.append(pd.DataFrame({
            "state": state, "time_ms": reference_times,
            "con_mean": mean[0, s], "inc_mean": mean[1, s],
            "con_sem": sem[0, s], "inc_sem": sem[1, s],
            "inc_minus_con_mean": diff_mean[s], "inc_minus_con_sem": diff_sem[s],
            "n_con_subjects": n_valid[0, s], "n_inc_subjects": n_valid[1, s],
            "n_paired_subjects": diff_n[s],
        }))
    pd.concat(curves, ignore_index=True).to_csv(figure_dir / f"{stem}.csv", index=False)
    figure_path = figure_dir / f"{stem}.png"
    make_figure(reference_times, probs, states, figure_path, bool(args.subjects))
    notes = [
        f"Included subjects: {len(ids)}; primary K={k}; states={states}.",
        "Selection: completed S2/S3 rules and primary Simon microstate recording QC.",
        "Trials: correct; no badEpoch; inclusive RT 150--1500 ms; minimum per condition from config.",
        "Full individual epochs are filtered/resampled/back-fitted BEFORE cropping.",
        f"Bandpass Hz: {cfg['preprocessing']['bandpass_hz']}; "
        f"sample rate Hz: {cfg['preprocessing']['resample_hz']}.",
        f"Backfit options: {ms.backfit_options(cfg)}.",
        "Batches do not concatenate epochs or smooth across trial boundaries.",
        "Primary P(state): shared src.features.state_probabilities; -1 trials omitted at each time.",
        "Zero labeled trials -> NaN. No extra trial-count threshold imposed after labeling.",
        "Alternative all-trial fractions and labeled-trial counts are also saved.",
        "Across subjects: equal weight; +/-1 SEM; available N recorded per time/condition.",
        "No subtraction of a probability baseline; no hypothesis test in S4.",
        "Trial rows saved as zero-based positions in the TSV; match_trials order assumption retained.",
        "Cache validity: code/config/template/package hashes and input file size/mtime/path.",
        "Resume reruns an interrupted subject but reuses compatible completed subjects.",
        "NPZ probability axes: subject, condition (con/inc), state, time.",
        "S5 should use paired subjects and report any non-finite data; do not replace NaN with zero.",
    ]
    (data_dir / f"{stem}_methods.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    print(f"\nIncluded subjects: {len(ids)}")
    print(f"Probability shape: {probs.shape} (subject, condition, state, time)")
    print(f"Saved arrays: {data_path}")
    print(f"Saved figure: {figure_path}")
    print(f"Subject checkpoints: {cache_dir}")
    if (n_labeled == 0).any():
        print("NOTE: some time points have no labeled trials; saved as NaN with denominator counts.")
    print("DONE: state probabilities created." if not args.subjects
          else "DONE: state probabilities subset check created.")


if __name__ == "__main__":
    main()
