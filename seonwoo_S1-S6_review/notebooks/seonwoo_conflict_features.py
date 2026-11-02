r"""S6: participant-level Simon conflict features and an S7 review bundle.

Keep the existing S3, S4 and S5 scripts in notebooks/. From the repo root:
    .\.venv\Scripts\python.exe -X utf8 notebooks/seonwoo_conflict_features.py

Input: completed S4 state probabilities, checked by the S5 input validator.
Output: results/features/conflict.csv (one row per S4 participant).
Feature = mean_t[P(state | incongruent) - P(state | congruent)], over the
configured n2 and theta time windows. Values are probability differences:
0.02 = +2 percentage points. No clinical columns or new hypothesis tests.

Window endpoints are inclusive; use only available samples, without interpolation.
If ANY required sample is missing in a window, that participant's five features
for that window remain blank/NaN. Keep the subject row and any other valid window.
Missingness is audited, never replaced with zero or averaged over a shorter window.
All states/windows are exported regardless of S5 significance.

Also builds results/seonwoo_S1-S6_review.zip from a fixed list of analysis files
for S7 interpretation. It contains no raw EEG or per-trial label checkpoints.
Nothing is sent or uploaded by this script.

Optional check using completed S2--S4 subset outputs:
    python -X utf8 notebooks/seonwoo_conflict_features.py --subjects sub-P001 sub-C004
"""

import argparse
import copy
import importlib.metadata
import json
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZIP_DEFLATED, BadZipFile, ZipFile

import numpy as np
import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
from src import features, io as pio
from seonwoo_conflict_stats import check, load_s4
from seonwoo_state_probabilities import digest_file

WINDOW_NAMES = ("n2", "theta")


def extract_features(data, cfg):
    ids = data["subject_id"].astype(str)
    states = data["state_names"].tolist()
    times = np.asarray(data["times_ms"], dtype=float)
    p = np.asarray(data["probabilities"], dtype=float)
    check(data["condition_codes"].tolist() == [1, 2], "Unexpected condition order.")
    check(p.shape == (len(ids), 2, len(states), len(times)), "Invalid probability shape.")
    check(len(ids) == len(set(ids)), "Duplicate subject IDs.")
    check(np.isfinite(times).all() and np.all(np.diff(times) > 0), "Invalid time axis.")
    check(not np.isinf(p).any(), "Infinite probability encountered.")
    finite = p[np.isfinite(p)]
    check(((finite >= 0) & (finite <= 1)).all(), "Probability outside [0, 1].")
    difference = p[:, 1] - p[:, 0]
    output = pd.DataFrame({"subject_id": ids})
    audit_rows, windows = [], []
    for name in WINDOW_NAMES:
        lo, hi = map(float, cfg["simon_conditions"]["windows_ms"][name])
        check(np.isfinite([lo, hi]).all() and times[0] <= lo < hi <= times[-1],
              f"The {name} window is not fully covered by S4.")
        mask = (times >= lo) & (times <= hi)
        check(mask.any(), f"No samples in the {name} window.")
        n_expected = int(mask.sum())
        paired_valid = np.isfinite(p[:, :, :, mask]).all(axis=(1, 2))
        n_valid = paired_valid.sum(axis=1)
        complete = n_valid == n_expected
        values = np.full((len(ids), len(states)), np.nan)
        for i in np.flatnonzero(complete):
            for s in range(len(states)):
                values[i, s] = features.window_mean(difference[i, s], times, (lo, hi))
        check(np.isfinite(values[complete]).all(), "Non-finite complete-window feature.")
        check(np.all(np.abs(values[complete]) <= 1 + 1e-12), "Conflict feature outside [-1, 1].")
        check(np.allclose(values[complete].sum(axis=1), 0, atol=1e-10),
              "Across-state conflict effects must sum to zero for normalized probabilities.")
        for s, state in enumerate(states):
            output[f"simon_conflict_{name}_{state}"] = values[:, s]
        info = {
            "window": name, "requested_start_ms": lo, "requested_end_ms": hi,
            "first_sample_ms": float(times[mask][0]),
            "last_sample_ms": float(times[mask][-1]), "n_expected_samples": n_expected,
            "n_complete_subjects": int(complete.sum()),
            "n_incomplete_subjects": int((~complete).sum()),
        }
        windows.append(info)
        for i, sid in enumerate(ids):
            audit_rows.append({
                "subject_id": sid, "window": name,
                "status": "complete" if complete[i] else "missing_probability",
                "n_expected_samples": n_expected, "n_valid_paired_samples": int(n_valid[i]),
                "first_sample_ms": info["first_sample_ms"],
                "last_sample_ms": info["last_sample_ms"],
            })
    check(output.shape == (len(ids), 1 + len(WINDOW_NAMES) * len(states)),
          "Unexpected output table shape.")
    return output, pd.DataFrame(audit_rows), windows


def save_csv_atomic(table, path):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".csv.tmp", delete=False,
                                         mode="w", encoding="utf-8", newline="") as file:
            temporary = Path(file.name)
            table.to_csv(file, index=False, na_rep="", float_format="%.12g")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def make_review_bundle(cfg, suffix, metadata, s6_paths):
    """Fixed allowlist: no directory recursion, raw EEG or per-trial labels."""
    results = pio.results_dir(cfg)
    feature_dir, stats_dir, figure_dir = (results / name for name in ("features", "stats", "figures"))
    s5_stem = f"seonwoo_conflict{suffix}"
    candidates = [
        figure_dir / f"seonwoo_erp_cz{suffix}.png",
        figure_dir / f"seonwoo_erp_segmentation{suffix}.png",
        figure_dir / f"seonwoo_state_probabilities{suffix}.png",
        figure_dir / f"seonwoo_conflict_timecourse{suffix}.png",
        feature_dir / f"seonwoo_erp_qc{suffix}.csv",
        feature_dir / f"seonwoo_erp_counts{suffix}.csv",
        feature_dir / f"seonwoo_erp_segmentation{suffix}_qc.csv",
        feature_dir / f"seonwoo_state_probabilities{suffix}_qc.csv",
        figure_dir / f"seonwoo_erp_cz{suffix}.csv",
        figure_dir / f"seonwoo_erp_cz{suffix}_summary.txt",
        figure_dir / f"seonwoo_state_probabilities{suffix}.csv",
        results / "erp" / f"seonwoo_erp_segmentation{suffix}_segments.csv",
        results / "erp" / f"seonwoo_erp_segmentation{suffix}_methods.txt",
        results / "erp" / f"seonwoo_state_probabilities{suffix}_methods.txt",
        results / "erp" / f"seonwoo_state_probabilities{suffix}.npz",
        *[stats_dir / f"{s5_stem}_{end}" for end in
          ("clusters.csv", "summary.txt", "subjects.csv", "timecourse.csv", "test.npz")],
        *s6_paths,
    ]
    s5_test = stats_dir / f"{s5_stem}_test.npz"
    s5_matches = False
    s5_note = "S5 test file is missing."
    if s5_test.is_file():
        try:
            with np.load(s5_test, allow_pickle=False) as saved:
                s5_meta = json.loads(str(saved["metadata_json"].item()))
            s5_matches = (s5_meta["input_s4_sha256"] == metadata["input_s4_sha256"]
                          and s5_meta["mode"] == metadata["mode"])
            s5_note = ("S5 and S6 used the same S4 input." if s5_matches
                       else "S5 input differs from S6. Rerun S5 before reviewing results together.")
        except (KeyError, ValueError, OSError, EOFError, BadZipFile):
            s5_note = "S5 metadata could not be checked."
    if not s5_matches:
        # A new S6 result remains valid, but stale S5 results must not be bundled as current.
        stale = {stats_dir / f"{s5_stem}_{end}" for end in
                 ("clusters.csv", "summary.txt", "subjects.csv", "timecourse.csv", "test.npz")}
        stale.add(figure_dir / f"seonwoo_conflict_timecourse{suffix}.png")
        candidates = [p for p in candidates if p not in stale]
    files, missing = {}, []
    for path in candidates:
        name = "results/" + path.relative_to(results).as_posix()
        if path.is_file():
            files[name] = path
        else:
            missing.append(name)
    for name in ("seonwoo_erp_cz_all.py", "seonwoo_erp_segmentation.py",
                 "seonwoo_state_probabilities.py", "seonwoo_conflict_stats.py",
                 "seonwoo_conflict_features.py"):
        path = REPO_ROOT / "notebooks" / name
        if path.is_file():
            files["notebooks/" + name] = path
        else:
            missing.append("notebooks/" + name)
    for name in ("io.py", "preprocess.py", "microstates.py", "features.py"):
        path = REPO_ROOT / "src" / name
        if path.is_file():
            files["src/" + name] = path
        else:
            missing.append("src/" + name)
    clean_cfg = copy.deepcopy(cfg)
    clean_cfg["paths"]["data_root"] = "[local Dataset_ALL path omitted]"
    manifest = {
        "purpose": "S7 result review; subset bundles are software checks only",
        "mode": metadata["mode"], "s5_matches_s6_input": s5_matches,
        "s5_note": s5_note, "missing_files": missing,
        "files_sha256": {name: digest_file(path) for name, path in files.items()},
    }
    path = results / f"seonwoo_S1-S6_review{suffix}.zip"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=results, suffix=".zip.tmp", delete=False) as file:
            temporary = Path(file.name)
        with ZipFile(temporary, "w", compression=ZIP_DEFLATED) as archive:
            for name, source in files.items():
                archive.write(source, arcname=name)
            archive.writestr("config_for_review.yaml", yaml.safe_dump(clean_cfg, sort_keys=False))
            archive.writestr("review_manifest.json", json.dumps(manifest, indent=2))
            archive.writestr("READ_ME_S7.txt",
                "Review the actual participant counts, QC, figures and corrected p-values before writing slides.\n"
                + s5_note + "\n"
                + "S6 retains all S4 subject rows; a window with missing values remains blank.\n"
                + "n2 and theta are window names, not measurements of N2 amplitude or theta power.\n"
                + "S5 cluster boundaries do not establish exact effect onset/offset.\n"
                + "This bundle omits raw EEG, trial labels and the local dataset path.\n"
                + "See review_manifest.json for input consistency, file hashes and missing files.\n")
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return path, missing, s5_matches, s5_note


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subjects", nargs="+", help="Subset check only; uses S2--S4 subset outputs.")
    args = parser.parse_args()
    suffix = "_subset" if args.subjects else ""
    cfg = pio.load_config()
    print("S6: checking saved S4 probabilities and participant IDs...", flush=True)
    data, s4_metadata, input_path = load_s4(cfg, args.subjects, suffix)
    table, audit, windows = extract_features(data, cfg)
    output_dir = pio.results_dir(cfg, "features")
    stats_dir = pio.results_dir(cfg, "stats")
    output_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"conflict{suffix}.csv"
    audit_path = output_dir / f"seonwoo_conflict{suffix}_qc.csv"
    methods_path = stats_dir / f"seonwoo_conflict_features{suffix}_methods.txt"
    metadata_path = stats_dir / f"seonwoo_conflict_features{suffix}_metadata.json"
    metadata = {
        "schema": "seonwoo_s6_v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode": s4_metadata["mode"], "input_s4_sha256": digest_file(input_path),
        "input_s4_pipeline_signature": s4_metadata["pipeline_signature"],
        "script_sha256": digest_file(__file__), "n_subjects": len(table),
        "state_names": data["state_names"].tolist(), "windows": windows,
        "contrast": "incongruent_minus_congruent", "units": "probability_difference",
        "probability_denominator": s4_metadata["probability_denominator"],
        "missing_rule": "any missing sample in either condition -> all state features for that window NaN",
        "columns": table.columns.tolist(),
        "versions": {name: importlib.metadata.version(name) for name in
                     ("numpy", "pandas", "scipy", "mne", "pycrostates", "h5py")},
    }
    notes = [
        "S6: participant-level Simon conflict features",
        "SUBSET CHECK ONLY" if args.subjects else "Full S4 cohort",
        f"Subjects: {len(table)}; feature columns: {table.shape[1]-1}; one row per subject_id.",
        "Definition: equal-time-sample mean of P(state|incongruent) minus P(state|congruent).",
        "Use shared src.features.window_mean with inclusive configured window endpoints.",
        "No interpolation at endpoints; use the existing S4 sample grid.",
        "Probabilities exclude unlabeled trials at each time point, as in S4.",
        "Values are probability differences, not percent change: +0.02 means +2 percentage points.",
        "Positive: higher probability in incongruent trials; negative: higher in congruent trials.",
        "Missing rule: require every sample of a window in both conditions/all states.",
        "An incomplete window yields NaN/blank features; retain subject ID and any other complete window.",
        "Do not replace missing values with zero or calculate a shorter-window mean.",
        "The table retains the S4 cohort, not S5's complete-whole-timecourse subset.",
        "No filtering by S5 significance; all pre-specified states and windows are exported.",
        "n2/theta name time windows; these features are NOT N2 ERP amplitude or theta-band power.",
        "State-wise differences sum to zero within each complete window because probabilities sum to one.",
        "No clinical columns, group comparison, new hypothesis tests, or FDR results in this step.",
        "Clinical variables should be merged later on subject_id; do not merge by row position.",
    ]
    for w in windows:
        notes.append(
            f"{w['window']}: requested {w['requested_start_ms']:g}--{w['requested_end_ms']:g} ms; "
            f"samples {w['first_sample_ms']:g}--{w['last_sample_ms']:g} ms; "
            f"n_samples={w['n_expected_samples']}; complete subjects={w['n_complete_subjects']}/{len(table)}."
        )
    notes.extend([f"Input SHA256: {metadata['input_s4_sha256']}",
                  "Task source: https://drive.google.com/file/d/1pcLE6YSJH0KPbDStXObmbNMqv2GRHCC-/view"])
    save_csv_atomic(table, output_path)
    save_csv_atomic(audit, audit_path)
    methods_path.write_text("\n".join(notes) + "\n", encoding="utf-8")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    # Validate the actual CSV representation used by teammates.
    saved = pd.read_csv(output_path)
    check(saved.columns.tolist() == table.columns.tolist()
          and saved["subject_id"].tolist() == table["subject_id"].tolist(),
          "CSV schema/subject order changed during saving.")
    check(np.allclose(saved.iloc[:, 1:].to_numpy(dtype=float),
                      table.iloc[:, 1:].to_numpy(dtype=float), equal_nan=True,
                      rtol=1e-10, atol=1e-12), "CSV values changed during saving.")
    print(f"Saved: {output_path}\nRows: {len(table)} | Feature columns: {table.shape[1]-1}")
    for w in windows:
        print(f"{w['window']}: {w['n_complete_subjects']}/{len(table)} complete; "
              f"{w['n_expected_samples']} samples ({w['first_sample_ms']:g}--{w['last_sample_ms']:g} ms)")
    if table.iloc[:, 1:].isna().any().any():
        print(f"NOTE: some features are blank due to missing probabilities. See {audit_path}")
    bundle_path, missing, s5_matches, s5_note = make_review_bundle(
        cfg, suffix, metadata, [output_path, audit_path, methods_path, metadata_path]
    )
    print(f"S7 review ZIP: {bundle_path}")
    if missing or not s5_matches:
        print(f"NOTE: review bundle is incomplete ({len(missing)} missing files). {s5_note}")
        print("See review_manifest.json inside the ZIP; S6 feature extraction itself is complete.")
    if args.subjects:
        print("SUBSET CHECK: do not submit this as the full-cohort result.")
    print("DONE: conflict features created." if not args.subjects
          else "DONE: conflict features subset check created.")


if __name__ == "__main__":
    try:
        main()
    except PermissionError:
        raise SystemExit("Close the output CSV/ZIP in Excel or other apps, then rerun the same command.")
