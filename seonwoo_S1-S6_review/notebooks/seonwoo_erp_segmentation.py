r"""S3: descriptive ERP segmentation with the team's approved shared templates.

Save in notebooks/ and run from the project root:
    .\.venv\Scripts\python.exe -X utf8 notebooks/seonwoo_erp_segmentation.py

Requires completed S2 outputs from seonwoo_erp_cz_all.py, src/, config.yaml,
results/templates/templates_K{primary_k}.fif, microstates_task.csv, and epochs.
Optional technical check, using the completed S2 *_subset audit:
    python -X utf8 notebooks/seonwoo_erp_segmentation.py --subjects sub-P001 sub-C004

Method: S2 trial/cohort rules plus the primary Simon microstate recording QC.
Average trials within each subject and condition; filter each full averaged
epoch at 2--20 Hz and resample to 250 Hz (config); then average subjects equally.
Filtering/resampling are linear and commute with trial averaging here.
Apply the existing shared A--E maps and shared smoothing to the two full grand
ERPs separately, then display 0 <= t < 800 ms. No new template clustering.
This produces a descriptive template-based ERP segmentation, not single-trial
state probabilities (S4), a permutation test (S5), or conflict.csv (S6).

No source/config edits. No clinical columns in output tables.
Epoch-to-TSV order follows the documented assumption in src.io.match_trials.
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import mne
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src import io as pio, microstates as ms

CONDITIONS = ("Congruent", "Incongruent")
WINDOW_MS = (0.0, 800.0)  # Half-open sample bins: [0, 800).


def choose_subjects(cfg, states, requested, suffix):
    audit_path = pio.results_dir(cfg, "features", f"seonwoo_erp_qc{suffix}.csv")
    if not audit_path.is_file():
        raise FileNotFoundError(f"Run seonwoo_erp_cz_all.py first. Missing: {audit_path}")
    audit = pd.read_csv(audit_path)
    if audit["subject_id"].duplicated().any():
        raise ValueError("Duplicate IDs in the S2 audit.")
    complete_statuses = {
        "included", "excluded_subject_rule", "excluded_trial_status",
        "excluded_recording_qc", "excluded_min_trials",
    }
    if not audit["status"].isin(complete_statuses).all():
        raise ValueError("S2 has unfinished/error rows. Complete S2 before S3.")
    expected = set(pio.load_subjects(cfg)["subject_id"])
    if requested:
        if set(requested) - expected:
            raise ValueError("Unknown subject IDs in --subjects.")
        expected = set(requested)
    if set(audit["subject_id"]) != expected:
        raise ValueError("S2 audit IDs do not match the requested cohort.")

    primary_path = pio.results_dir(cfg, "features", "microstates_task.csv")
    primary = pd.read_csv(primary_path).set_index("subject_id")
    if primary.index.duplicated().any():
        raise ValueError("Duplicate IDs in microstates_task.csv.")
    selected = audit.loc[audit["status"].eq("included"), "subject_id"].tolist()
    missing = set(selected) - set(primary.index)
    if missing:
        raise ValueError(f"Subjects missing from primary metrics: {sorted(missing)}")
    columns = [f"simon_cov_{state}" for state in states]
    values = primary.loc[selected, columns]
    present = values.notna()
    if (present.any(axis=1) & ~present.all(axis=1)).any():
        raise ValueError("Partially missing Simon coverage; check the primary metrics.")
    passing = present.all(axis=1)
    if not np.isfinite(values.loc[passing].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite Simon coverage in the primary metrics.")

    # Primary metrics leave QC-excluded Simon recordings empty, per PROGRESS.md.
    ids = values.index[passing].tolist()
    excluded = set(selected) - set(ids)
    audit = audit.rename(columns={"status": "s2_status"})
    audit["status"] = audit["s2_status"]
    audit.loc[audit["subject_id"].isin(excluded), "status"] = "excluded_microstate_qc"
    audit.loc[audit["subject_id"].isin(ids), "status"] = "pending_erp"
    print(f"S2 included subjects: {len(selected)}")
    print(f"Additional Simon microstate QC exclusions: {len(excluded)}"
          + (f" ({', '.join(sorted(excluded))})" if excluded else ""))
    if len(ids) < 2:
        raise ValueError("Fewer than 2 subjects remain for the grand-average ERP.")
    return ids, audit


def require_epoch_files(ids, cfg, path):
    missing = []
    for sid in ids:
        for kind, file in [
            ("epochs", pio.epochs_path(sid, "Simon", "cue", cfg)),
            ("trials", pio.trials_path(sid, "Simon", cfg)),
        ]:
            if not file.is_file():
                missing.append({"subject_id": sid, "kind": kind, "path": str(file)})
    pd.DataFrame(missing, columns=["subject_id", "kind", "path"]).to_csv(
        path, index=False, encoding="utf-8-sig"
    )
    if missing:
        raise FileNotFoundError(f"STOP: {len(missing)} missing files. See {path}")


def average_subject(sid, cfg, channel_names, audit_row):
    epochs = pio.load_epochs(sid, "Simon", "cue", cfg)
    trials = pio.match_trials(epochs, pio.load_trials(sid, "Simon", cfg), "cue")
    keep = (
        trials["beh_accuracy"].eq(1)
        & ~trials["bad_epoch"]
        & trials["beh_rt"].between(150, 1500)
    ).to_numpy()
    condition = trials["beh_congruentcondition"].to_numpy()
    if not np.isin(condition[keep], [1, 2]).all():
        raise ValueError(f"{sid}: unknown condition in a retained trial.")
    masks = [keep & (condition == value) for value in (1, 2)]
    counts = [int(mask.sum()) for mask in masks]
    expected = [int(audit_row["n_congruent"]), int(audit_row["n_incongruent"])]
    if counts != expected:
        raise ValueError(f"{sid}: trial counts changed since S2: {counts} vs {expected}.")
    if min(counts) < int(cfg["subjects"]["simon"]["min_trials_per_condition"]):
        raise ValueError(f"{sid}: insufficient trials under the current config.")
    baseline = np.asarray(epochs["erp_baseline"]).ravel()
    if baseline.shape != (2,) or not np.allclose(baseline, [-200, 0]):
        raise ValueError(f"{sid}: unexpected supplied baseline.")
    if set(epochs["ch_names"]) != set(channel_names):
        raise ValueError(f"{sid}: channels differ from the approved templates.")
    if not np.allclose(np.diff(epochs["times"]), 1000.0 / epochs["srate"]):
        raise ValueError(f"{sid}: time axis and sample rate differ.")

    # 2 x channels x full epoch; float64 averaging limits accumulation error.
    means_uv = np.stack([
        epochs["data"][mask].mean(axis=0, dtype=np.float64) for mask in masks
    ])
    if not np.isfinite(means_uv).all():
        raise ValueError(f"{sid}: non-finite values in averaged ERP.")
    epochs["data"] = means_uv
    # MNE renamed this identical montage; support both project-era and newer MNE.
    montage = ("colin27_1005" if "colin27_1005" in mne.channels.get_builtin_montages()
               else "standard_1005")
    inst = pio.epochs_to_mne(epochs, montage=montage)  # µV -> V once.
    del epochs
    inst.reorder_channels(channel_names)
    lo, hi = cfg["preprocessing"]["bandpass_hz"]
    inst.filter(lo, hi, verbose="error")
    inst.resample(cfg["preprocessing"]["resample_hz"], verbose="error")
    return inst


def segment_rows(labels, times, sfreq, states):
    """Export contiguous sample bins, not estimates of whole-recording duration."""
    rows = []
    dt = 1000.0 / sfreq
    for condition, sequence in zip(CONDITIONS, labels):
        boundaries = np.r_[0, np.flatnonzero(np.diff(sequence)) + 1, len(sequence)]
        for start, stop in zip(boundaries[:-1], boundaries[1:]):
            label = int(sequence[start])
            rows.append({
                "condition": condition,
                "state": states[label] if label >= 0 else "unlabeled",
                "start_ms": float(times[start]),
                "stop_ms_exclusive": float(times[stop - 1] + dt),
                "duration_ms": float((stop - start) * dt),
                "n_samples": int(stop - start),
                "touches_window_start": bool(start == 0),
                "touches_window_end": bool(stop == len(sequence)),
            })
    return pd.DataFrame(rows)


def plot_segmentation(times, labels, gfp_uv, segments, templates, info, n, cfg, path, subset):
    states = templates.cluster_names
    k = len(states)
    colors = plt.get_cmap("tab10")(np.arange(k))
    color_for = dict(zip(states, colors))
    color_for["unlabeled"] = "#BDBDBD"
    maps = templates.cluster_centers_.copy()
    maps -= maps.mean(axis=1, keepdims=True)
    maps /= maps.std(axis=1, keepdims=True)
    vmax = np.max(np.abs(maps))
    fig = plt.figure(figsize=(12, 8))
    grid = fig.add_gridspec(
        5, k + 1, width_ratios=[1] * k + [0.08],
        height_ratios=[1.4, 1.25, 0.23, 1.25, 0.23],
        left=0.075, right=0.95, bottom=0.16, top=0.86,
        hspace=0.38, wspace=0.25,
    )
    for i, state in enumerate(states):
        ax = fig.add_subplot(grid[0, i])
        im, _ = mne.viz.plot_topomap(
            maps[i], info, axes=ax, show=False, contours=0,
            cmap="RdBu_r", vlim=(-vmax, vmax), sensors=True,
        )
        ax.set_title(state, color=colors[i], fontweight="bold", fontsize=12)
    cax = fig.add_subplot(grid[0, -1])
    fig.colorbar(im, cax=cax, label="Normalized units")
    for i, condition in enumerate(CONDITIONS):
        ax = fig.add_subplot(grid[1 + 2 * i, :])
        bar = fig.add_subplot(grid[2 + 2 * i, :], sharex=ax)
        ax.plot(times, gfp_uv[i], color="#333333", lw=1.2)
        for row in segments[segments["condition"].eq(condition)].itertuples(index=False):
            within = (times >= row.start_ms) & (times < row.stop_ms_exclusive)
            x = np.r_[times[within], row.stop_ms_exclusive]
            y = np.r_[gfp_uv[i, within],
                      np.interp(row.stop_ms_exclusive, times, gfp_uv[i])]
            ax.fill_between(x, y, color=color_for[row.state], alpha=0.55)
        ax.set_title(condition, loc="left", fontsize=11, fontweight="bold")
        ax.set_ylabel("GFP (µV)")
        ax.set_ylim(bottom=0)
        ax.tick_params(axis="x", labelbottom=False)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
        for row in segments[segments["condition"].eq(condition)].itertuples(index=False):
            bar.broken_barh([(row.start_ms, row.duration_ms)], (0, 1),
                            facecolors=color_for[row.state])
            if row.duration_ms >= 28:
                bar.text(row.start_ms + row.duration_ms / 2, 0.5,
                         row.state if row.state != "unlabeled" else "—",
                         ha="center", va="center", fontsize=8, color="white",
                         fontweight="bold")
        bar.set_ylim(0, 1)
        bar.set_xlim(*WINDOW_MS)
        bar.set_yticks([])
        bar.set_ylabel("State", rotation=0, ha="right", va="center", fontsize=9)
        bar.spines[["top", "right", "left", "bottom"]].set_visible(False)
        bar.set_xticks(np.arange(0, 801, 100))
        if i == 0:
            bar.tick_params(axis="x", labelbottom=False)
        else:
            bar.set_xlabel("Time from cue onset (ms)")
    mode = "SUBSET CHECK" if subset else "grand-average ERPs"
    fig.suptitle(f"Simon ERP segmentation | N = {n} | K = {k} | {mode}",
                 y=0.97, fontsize=15, fontweight="bold")
    fig.text(0.5, 0.92, "Approved shared templates · normalized scalp maps · polarity invariant",
             ha="center", fontsize=10, color="#555555")
    fig.legend(
        handles=[Patch(color=color_for[s], label=s) for s in states + ["unlabeled"]],
        loc="lower center", bbox_to_anchor=(0.5, 0.065), ncol=k + 1, frameon=False,
    )
    lo, hi = cfg["preprocessing"]["bandpass_hz"]
    fig.text(
        0.5, 0.018,
        f"Equal subject weights; {lo:g}–{hi:g} Hz, "
        f"{cfg['preprocessing']['resample_hz']:g} Hz sampling; shared temporal smoothing.\n"
        "GFP = spatial standard deviation. Descriptive ERP labels; no significance test.",
        ha="center", fontsize=9, color="#555555",
    )
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", help="Subset check; needs S2 subset audit.")
    args = parser.parse_args()
    suffix = "_subset" if args.subjects else ""
    cfg = pio.load_config()
    k = int(cfg["microstates"]["primary_k"])
    templates = ms.load_templates(k, cfg)
    states = templates.cluster_names
    confirmed = cfg["microstates"]["confirmed_labels"].get(k)
    if confirmed is None or states != confirmed:
        raise ValueError("Template state names do not match the confirmed config labels.")
    channel_names = list(templates.info["ch_names"])
    ids, audit = choose_subjects(cfg, states, args.subjects, suffix)
    stem = f"seonwoo_erp_segmentation{suffix}"
    fig_dir = pio.results_dir(cfg, "figures")
    data_dir = pio.results_dir(cfg, "erp")
    audit_path = pio.results_dir(cfg, "features", f"{stem}_qc.csv")
    for folder in (fig_dir, data_dir, audit_path.parent):
        folder.mkdir(parents=True, exist_ok=True)
    audit.to_csv(audit_path, index=False, encoding="utf-8-sig")
    require_epoch_files(ids, cfg, pio.results_dir(cfg, "features", f"{stem}_missing_files.csv"))
    by_id = audit.set_index("subject_id")
    grand_sum, reference_times, reference_info = None, None, None
    for i, sid in enumerate(ids, start=1):
        print(f"[{i}/{len(ids)}] {sid}: averaging 64-channel ERPs...", flush=True)
        try:
            inst = average_subject(sid, cfg, channel_names, by_id.loc[sid])
            if reference_times is None:
                reference_times = inst.times.copy()
                reference_info = inst.info.copy()
                grand_sum = inst.get_data().copy()
            else:
                if inst.times.shape != reference_times.shape or not np.allclose(
                    inst.times, reference_times, rtol=0, atol=1e-10
                ):
                    raise ValueError(f"{sid}: mismatched time grid.")
                grand_sum += inst.get_data()
            audit.loc[audit["subject_id"].eq(sid), "status"] = "included"
            del inst
        except Exception:
            audit.loc[audit["subject_id"].eq(sid), "status"] = "error"
            audit.to_csv(audit_path, index=False, encoding="utf-8-sig")
            raise
    grand = grand_sum / len(ids)
    audit.to_csv(audit_path, index=False, encoding="utf-8-sig")
    inst = mne.EpochsArray(grand, reference_info, tmin=reference_times[0],
                           baseline=None, verbose="error")
    # Each condition is a separate epoch; smoothing never crosses conditions.
    full_labels = ms.backfit(inst, templates, cfg).labels
    full_times_ms = np.round(inst.times * 1000, 6)
    window = (full_times_ms >= WINDOW_MS[0]) & (full_times_ms < WINDOW_MS[1])
    times = full_times_ms[window]
    labels = full_labels[:, window]
    gfp_uv = grand.std(axis=1)[:, window] * 1e6
    segments = segment_rows(labels, times, reference_info["sfreq"], states)
    png_path = fig_dir / f"{stem}.png"
    plot_segmentation(times, labels, gfp_uv, segments, templates, reference_info,
                      len(ids), cfg, png_path, bool(args.subjects))
    segments.to_csv(data_dir / f"{stem}_segments.csv", index=False)
    frames = []
    for i, condition in enumerate(CONDITIONS):
        frames.append(pd.DataFrame({
            "condition": condition, "time_ms": times,
            "state": [states[x] if x >= 0 else "unlabeled" for x in labels[i]],
            "label_index": labels[i], "gfp_uv": gfp_uv[i],
        }))
    pd.concat(frames, ignore_index=True).to_csv(
        data_dir / f"{stem}_timecourse.csv", index=False
    )
    np.savez_compressed(
        data_dir / f"{stem}_grand_erp.npz",
        subject_id=np.asarray(ids, dtype=str), condition=np.asarray(CONDITIONS),
        ch_names=np.asarray(channel_names), times_ms=full_times_ms,
        erp_uv=grand * 1e6, labels=full_labels, state_names=np.asarray(states),
        sfreq_hz=float(reference_info["sfreq"]),
    )
    notes = [
        f"Run UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Mode: {'subset check' if args.subjects else 'full cohort from completed S2 audit'}",
        f"Included subjects: {len(ids)}",
        f"Shared templates: {ms.templates_path(k, cfg)}",
        f"States: {states}",
        "Method choice: map condition-wise grand ERPs to the approved shared templates.",
        "This is template-based segmentation; no ERP-specific clusters were estimated.",
        "Same subject/trial rules as completed S2; trial counts checked against S2 audit.",
        "Additional Simon recording QC: require nonmissing primary simon_cov_<state>.",
        "Within-subject condition means first; equal subject weights for grand ERPs.",
        f"Filter band Hz: {cfg['preprocessing']['bandpass_hz']}; "
        f"resampled Hz: {cfg['preprocessing']['resample_hz']}.",
        "Filter/resample each FULL condition ERP; supplied baseline preserved.",
        f"Backfit options: {ms.backfit_options(cfg)}.",
        "Backfit two full grand ERPs separately before cropping; -1 remains unlabeled.",
        "Displayed sample bins: 0 <= t < 800 ms (stop_ms_exclusive is a bin boundary).",
        "Segment lengths describe this averaged waveform and are not task-level durations.",
        "Polarity invariant; topomaps normalized to unit spatial standard deviation.",
        "Descriptive only; S4 single-trial probabilities and S5 inference remain.",
        "Epoch/TSV order follows src.io.match_trials: count verified, order assumed.",
    ]
    for condition, sequence in zip(CONDITIONS, labels):
        text = f"{condition}: {int((sequence == -1).sum())}/{len(sequence)} unlabeled samples"
        print(text)
        notes.append(text)
    (data_dir / f"{stem}_methods.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    print(f"Included subjects: {len(ids)}")
    print(f"Saved figure: {png_path}")
    print(f"Saved segment table: {data_dir / (stem + '_segments.csv')}")
    print("DONE: ERP segmentation created." if not args.subjects
          else "DONE: ERP segmentation subset check created.")


if __name__ == "__main__":
    main()
