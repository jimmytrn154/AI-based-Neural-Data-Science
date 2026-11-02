r"""Seonwoo S1/S2: cohort Simon ERP at Cz, using the shared project rules.

Place in notebooks/, then run from the project root:
    .\.venv\Scripts\python.exe -X utf8 notebooks/seonwoo_erp_cz_all.py

Optional technical check (writes separate *_subset outputs):
    python -X utf8 notebooks/seonwoo_erp_cz_all.py --subjects sub-P001 sub-C004

Needs subjects.csv, Simon .set headers (not .fdt), cue epochs and trial TSVs.
Processes one subject at a time. Missing inputs stop the run with a CSV list;
missing files are never treated as scientific exclusions.

Selection: src.preprocess.select_subjects; correct, non-bad trials with
inclusive RT 150--1500 ms; at least config's minimum trials per condition.
ERP: original baseline and sample rate, configured 30 Hz low-pass on the
full epoch, then crop; subject means first, equal weight across subjects.
Microstate template-GEV exclusions are not ERP exclusions (config.yaml).
Figure 2 in the Dataset Guide uses a different cohort / trial selection;
this implements the team's ERP sanity check, not an exact numeric replica.
Epoch/TSV order follows src.io.match_trials (count verified, order assumed).
"""

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src import io as pio, preprocess as pp


def require_files(subject_ids, cfg, kinds, report_path):
    """Stop if data are missing instead of silently shrinking the cohort."""
    missing = []
    for sid in subject_ids:
        paths = {
            "header": pio.continuous_path(sid, "Simon", cfg),
            "epochs": pio.epochs_path(sid, "Simon", "cue", cfg),
            "trials": pio.trials_path(sid, "Simon", cfg),
        }
        for kind in kinds:
            if not paths[kind].is_file():
                missing.append({
                    "subject_id": sid, "kind": kind, "path": str(paths[kind])
                })
    pd.DataFrame(missing, columns=["subject_id", "kind", "path"]).to_csv(
        report_path, index=False, encoding="utf-8-sig"
    )
    if missing:
        print(f"\nSTOP: {len(missing)} required files are missing.", flush=True)
        for row in missing[:5]:
            print(row["path"])
        raise FileNotFoundError(
            f"Check Dataset_ALL and config.yaml. Full missing-file list: {report_path}"
        )


def subject_erp(sid, cfg):
    epochs = pio.load_epochs(sid, "Simon", "cue", cfg)
    trials = pio.match_trials(epochs, pio.load_trials(sid, "Simon", cfg), "cue")
    # Same inclusive RT rule as notebooks/example_backfit.py.
    keep = (
        trials["beh_accuracy"].eq(1)
        & ~trials["bad_epoch"]
        & trials["beh_rt"].between(150, 1500)
    ).to_numpy()
    condition = trials.loc[keep, "beh_congruentcondition"].to_numpy()
    if not np.isin(condition, [1, 2]).all():
        raise ValueError(f"{sid}: unknown / missing condition in a retained trial.")
    con = condition == 1
    inc = condition == 2
    counts = {
        "n_total": len(trials),
        "n_kept": int(keep.sum()),
        "n_congruent": int(con.sum()),
        "n_incongruent": int(inc.sum()),
    }
    minimum = int(cfg["subjects"]["simon"]["min_trials_per_condition"])
    if min(counts["n_congruent"], counts["n_incongruent"]) < minimum:
        return counts, None, None

    baseline = np.asarray(epochs["erp_baseline"], dtype=float).ravel()
    if baseline.shape != (2,) or not np.allclose(baseline, [-200, 0]):
        raise ValueError(f"{sid}: unexpected supplied ERP baseline: {baseline}.")
    times = np.asarray(epochs["times"], dtype=float)
    sfreq = float(epochs["srate"])
    if not np.allclose(np.diff(times), 1000.0 / sfreq):
        raise ValueError(f"{sid}: sample rate and epoch time axis do not match.")
    cz_index = epochs["ch_names"].index("Cz")
    cz_uv = epochs["data"][keep, cz_index, :].astype(np.float64)
    del epochs
    if not np.isfinite(cz_uv).all():
        raise ValueError(f"{sid}: non-finite Cz values in retained trials.")

    cz_uv = mne.filter.filter_data(
        cz_uv, sfreq=sfreq, l_freq=None,
        h_freq=float(cfg["simon_conditions"]["erp_plot_lowpass_hz"]),
        method="fir", phase="zero", fir_design="firwin", verbose="error",
    )
    start, stop = cfg["simon_conditions"]["analysis_window_ms"]
    window = (times >= start) & (times <= stop)
    if not window.any():
        raise ValueError(f"{sid}: empty analysis window.")
    # Units remain microvolts. The supplied -200..0 ms baseline is preserved.
    means = np.stack([
        cz_uv[con].mean(axis=0)[window],
        cz_uv[inc].mean(axis=0)[window],
    ])
    return counts, times[window], means


def summarize_curves(subject_means):
    """Equal subject weights; paired differences for the difference SEM."""
    con = subject_means[:, 0, :]
    inc = subject_means[:, 1, :]
    diff = inc - con
    curves = np.stack([con, inc, diff], axis=1)
    return curves.mean(axis=0), curves.std(axis=0, ddof=1) / np.sqrt(len(curves))


def draw_figure(times, means, sem, n, lowpass, minimum, path, subset):
    fig, axes = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]}, layout="constrained"
    )
    for i, (label, color) in enumerate([
        ("Congruent", "#2463A6"), ("Incongruent", "#D65C3A")
    ]):
        axes[0].plot(times, means[i], color=color, lw=1.8, label=label)
        axes[0].fill_between(
            times, means[i] - sem[i], means[i] + sem[i], color=color, alpha=0.15
        )
    axes[0].set_ylabel("Cz amplitude (µV)")
    mode = "SUBSET CHECK" if subset else "team-selected cohort"
    axes[0].set_title(f"Simon ERP at Cz | N = {n} subjects | {mode}")
    axes[0].legend(frameon=False, loc="best")
    axes[1].plot(times, means[2], color="#7549A0", lw=1.8)
    axes[1].fill_between(
        times, means[2] - sem[2], means[2] + sem[2],
        color="#7549A0", alpha=0.15
    )
    axes[1].set_ylabel("Incongruent −\ncongruent (µV)")
    axes[1].set_xlabel("Time from cue onset (ms)")
    for ax in axes:
        ax.axhline(0, color="0.6", lw=0.7)
        ax.axvline(0, color="0.5", lw=0.8, linestyle="--")
        ax.set_xlim(times[0], times[-1])
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    fig.get_layout_engine().set(rect=(0, 0.09, 1, 0.91))
    fig.text(
        0.5, 0.012,
        f"Correct; RT 150–1500 ms; bad epochs excluded; ≥{minimum} trials/condition; "
        f"{lowpass:g} Hz low-pass.\n"
        "Equal subject weights; shading = ±1 SEM across subjects, not significance.",
        ha="center", va="bottom", fontsize=8, color="0.35",
    )
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", help="Optional subset check only.")
    args = parser.parse_args()
    cfg = pio.load_config()
    rules = cfg["subjects"]
    simon = rules["simon"]
    id_col = rules["id_column"]
    subjects = pio.load_subjects(cfg)
    if subjects[id_col].isna().any() or subjects[id_col].duplicated().any():
        raise ValueError("subjects.csv must contain unique, non-missing subject IDs.")
    if args.subjects:
        unknown = set(args.subjects) - set(subjects[id_col])
        if unknown:
            parser.error(f"Unknown subject IDs: {sorted(unknown)}")
        subjects = subjects[subjects[id_col].isin(args.subjects)].copy()

    suffix = "_subset" if args.subjects else ""
    stem = f"seonwoo_erp_cz{suffix}"
    fig_dir = pio.results_dir(cfg, "figures")
    audit_dir = pio.results_dir(cfg, "features")
    array_dir = pio.results_dir(cfg, "erp")
    for folder in (fig_dir, audit_dir, array_dir):
        folder.mkdir(parents=True, exist_ok=True)
    audit_path = audit_dir / f"seonwoo_erp_qc{suffix}.csv"
    missing_path = audit_dir / f"seonwoo_erp_missing_files{suffix}.csv"

    audit = pd.DataFrame(index=pd.Index(subjects[id_col], name="subject_id"))
    audit["status"] = "excluded_subject_rule"
    candidates = subjects.query(simon["require"])
    audit.loc[candidates[id_col], "status"] = "excluded_trial_status"
    candidates = candidates[
        candidates["trials_simon"].isin(simon["allowed_trial_status"])
    ]
    if candidates.empty:
        raise ValueError("No subjects pass the Simon metadata rules.")
    audit.loc[candidates[id_col], "status"] = "pending_recording_qc"
    audit.to_csv(audit_path, encoding="utf-8-sig")

    print(f"Dataset: {pio.data_root(cfg)}", flush=True)
    print(f"Checking inputs for {len(candidates)} Simon candidates...", flush=True)
    require_files(candidates[id_col], cfg, ["header"], missing_path)
    eligible, stages = pp.select_subjects("Simon", cfg, subjects=subjects)
    audit.loc[candidates[id_col], "status"] = "excluded_recording_qc"
    audit.loc[eligible[id_col], "status"] = "pending_epochs"
    audit["nInterp"] = np.nan
    audit.loc[eligible[id_col], "nInterp"] = eligible["nInterp"].to_numpy()
    # Retain QC values also for excluded recordings.
    for sid in set(candidates[id_col]) - set(eligible[id_col]):
        audit.loc[sid, "nInterp"] = pio.read_qc(sid, "Simon", cfg)["nInterp"]
    audit.to_csv(audit_path, encoding="utf-8-sig")
    if eligible.empty:
        raise ValueError("No subjects pass recording QC.")
    require_files(eligible[id_col], cfg, ["epochs", "trials"], missing_path)
    print(stages.to_string(index=False), flush=True)

    all_means, included_ids = [], []
    reference_times = None
    for i, sid in enumerate(eligible[id_col], start=1):
        print(f"[{i}/{len(eligible)}] {sid}: loading...", flush=True)
        try:
            counts, times, means = subject_erp(sid, cfg)
            for key, value in counts.items():
                audit.loc[sid, key] = value
            if means is None:
                audit.loc[sid, "status"] = "excluded_min_trials"
                print(f"  EXCLUDED: con={counts['n_congruent']}, "
                      f"inc={counts['n_incongruent']} (minimum trial rule).", flush=True)
                continue
            if reference_times is None:
                reference_times = times
            elif times.shape != reference_times.shape or not np.allclose(
                times, reference_times, rtol=0, atol=1e-8
            ):
                raise ValueError(f"{sid}: time grid differs across subjects.")
            all_means.append(means)
            included_ids.append(sid)
            audit.loc[sid, "status"] = "included"
            print(f"  INCLUDED: con={counts['n_congruent']}, "
                  f"inc={counts['n_incongruent']}", flush=True)
        except Exception:
            audit.loc[sid, "status"] = "error"
            audit.to_csv(audit_path, encoding="utf-8-sig")
            raise

    audit.to_csv(audit_path, encoding="utf-8-sig")
    n = len(included_ids)
    stages.loc[len(stages)] = ["minimum trials in BOTH conditions", n]
    stages.to_csv(audit_dir / f"seonwoo_erp_counts{suffix}.csv", index=False)
    if n < 2:
        raise ValueError(f"Only {n} subjects remain; at least 2 are needed for SEM.")
    subject_means = np.stack(all_means)
    mean, sem = summarize_curves(subject_means)
    minimum = int(simon["min_trials_per_condition"])
    lowpass = float(cfg["simon_conditions"]["erp_plot_lowpass_hz"])
    png_path = fig_dir / f"{stem}.png"
    draw_figure(reference_times, mean, sem, n, lowpass, minimum, png_path,
                subset=bool(args.subjects))

    curves = {"time_ms": reference_times}
    for i, name in enumerate(["congruent", "incongruent", "incongruent_minus_congruent"]):
        curves[f"{name}_uv"] = mean[i]
        curves[f"{name}_sem_uv"] = sem[i]
    pd.DataFrame(curves).to_csv(fig_dir / f"{stem}.csv", index=False)
    np.savez_compressed(
        array_dir / f"{stem}_subject_means.npz",
        subject_id=np.asarray(included_ids, dtype=str),
        condition=np.asarray(["congruent", "incongruent"]),
        times_ms=reference_times, erp_uv=subject_means,
    )
    report = [
        f"Run UTC: {datetime.now(timezone.utc).isoformat()}",
        f"Mode: {'subset check' if args.subjects else 'complete team-rule cohort'}",
        f"Included subjects: {n}",
        f"Total retained congruent trials: {int(audit.loc[included_ids, 'n_congruent'].sum())}",
        f"Total retained incongruent trials: {int(audit.loc[included_ids, 'n_incongruent'].sum())}",
        "",
        stages.to_string(index=False),
        "",
        "Subject rules: src.preprocess.select_subjects, from config.yaml.",
        f"Maximum interpolated channels: {rules['max_interpolated_channels']}.",
        "Trial rules: correct; no badEpoch; inclusive RT 150--1500 ms.",
        f"At least {minimum} retained trials in EACH condition.",
        "Condition 1=congruent; 2=incongruent.",
        f"ERP: {lowpass:g} Hz zero-phase FIR low-pass on the full epoch, then crop.",
        "Original sample rate and supplied -200..0 ms baseline preserved.",
        "Within-subject means first; equal weight for every included subject.",
        "Shading: +/-1 SEM across subjects; difference SEM uses paired differences.",
        "Descriptive ERP sanity check; no significance test performed.",
        "Template-GEV QC is specified for microstates and is not applied to this ERP.",
        "Dataset Guide Figure 2 used 145 subjects with different selection rules.",
        "Alignment: src.io.match_trials checks counts; epoch-to-TSV order is assumed.",
        "Status codes: excluded_recording_qc = nInterp above limit or missing QC value;",
        "excluded_min_trials = fewer than the minimum in at least one condition.",
    ]
    (fig_dir / f"{stem}_summary.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"\nIncluded subjects: {n}")
    print(f"Saved figure: {png_path}")
    print(f"Saved selection / trial counts: {audit_path}")
    print(f"Saved methods and counts: {fig_dir / (stem + '_summary.txt')}")
    print("DONE: cohort ERP figure created." if not args.subjects
          else "DONE: subset check created; this is not the full-cohort result.")


if __name__ == "__main__":
    main()
