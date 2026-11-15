"""Run B3 -> B2 -> B4 -> B5 -> B6 for one or more (microstate, conflict) file pairs.

========================  HOW TO USE  ========================
1. Put the input CSVs in data/  (e.g. data/microstate_A.csv, data/conflict_A.csv)
2. Edit RUNS below (one line per run), then from the project folder run:
       python src/run_pipeline.py
   Or, without editing anything, run a single pair from the command line:
       python src/run_pipeline.py --micro microstate_B.csv --conflict conflict_B.csv
       python src/run_piperrline.py --micro microstate_B.csv --conflict conflict_B.csv --tag B_v2
   Only some runs from RUNS:
       python src/run_pipeline.py --only A C

Outputs for tag A:
    results/A/*_A.csv  (+ run_log_A.txt, run_config_A.json)
    figures/A/*_A.png
Overview of all runs: results/summary_all_runs.csv (one row per tag, updated each run)
=============================================================
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pipeline_io as io

# ------------------------------------------------------------------
# EDIT HERE: one dict per run.
#   tag       : output folder/suffix (results/<tag>/..._<tag>.csv)
#               omit it to use the part after the last '_' of the micro file
#   micro     : microstate CSV in data/   ('.csv' may be omitted)
#   conflict  : conflict CSV in data/
#   simon_qc  : True = drop trials_simon 'exclude' subjects (original rule)
# ------------------------------------------------------------------
RUNS = [
    {"tag": "rest_noblink", "micro": "microstates_rest_noblink.csv", "conflict": None, "simon_qc": True},
    {"tag": "task", "micro": "microstates_task.csv", "conflict": "conflict.csv", "simon_qc": True},
    {"tag": "task_K4", "micro": "microstates_task_K4.csv", "conflict": None, "simon_qc": True},
    {"tag": "task_nosmooth", "micro": "microstates_task_nosmooth.csv", "conflict": None, "simon_qc": True},
    {"tag": "task_qcincluded", "micro": "microstates_task_qcincluded.csv", "conflict": None, "simon_qc": True},
]
# {"tag": "B", "micro": "microstate_B.csv", "conflict": "conflict_B.csv", "simon_qc": True}
# ------------------------------------------------------------------

STEPS = [
    ("B3_reorganization", "reorganization.py"),
    ("B2_rest_vs_simon", "B2_rest_vs_simon.py"),
    ("B4_transition_network", "B4_transition_network.py"),
    ("B5_subject_clustering", "B5_subject_clustering.py"),
    ("B6_clinical_associations", "B6_clinical_associations.py"),
]
REQUIRES_REORG = {"B5_subject_clustering", "B6_clinical_associations"}


def preflight(ctx):
    """Check inputs before running anything; returns a list of problems."""
    problems = []
    try:
        micro = io.load_microstates(ctx, verbose=False)
        states = io.detect_states(micro, verbose=False)
        print(f"  microstate: {len(micro)} subjects, K = {len(states)} ({states})")
    except Exception as error:
        return [f"microstate file: {error}"]
    for name in (io.SUBJECTS_FILE, io.BEHAVIOR_FILE):
        if not (io.DATA_DIR / name).is_file():
            problems.append(f"data/{name} not found")
    if not ctx.conflict_file:
        print("  NOTE: No conflict file; B6 will be skipped.")
    elif not (io.DATA_DIR / ctx.conflict_file).is_file():
        print("  NOTE: Conflict file not found; B6 will be skipped.")
    else:
        conflict = io.read_unique(io.DATA_DIR / ctx.conflict_file)
        found = io.conflict_states(conflict)
        print(f"  conflict:   {len(conflict)} subjects, states n2 = {found['n2'] or '-'}, "
              f"theta = {found['theta'] or '-'}")
        for window, letters in found.items():
            if not letters:
                problems.append(f"{ctx.conflict_file}: no simon_conflict_{window}_<S> columns")
            elif letters != states:
                print(f"  WARNING: conflict {window} states ({letters}) differ from "
                      f"microstate states ({states})")
        overlap = len(set(micro.subject_id) & set(conflict.subject_id))
        print(f"  subjects in both files: {overlap}")
    return problems


def run_one(ctx, stop_on_error):
    args = ["--micro", ctx.micro_file, "--tag", ctx.name]
    if ctx.conflict_file:
        args += ["--conflict", ctx.conflict_file]
    if not ctx.simon_qc:
        args.append("--no-simon-qc")

    status = {}
    log_path = ctx.file("run_log", "txt")
    with open(log_path, "w", encoding="utf-8") as log:
        for label, script in STEPS:
            if label in REQUIRES_REORG and status.get("B3_reorganization") != "ok":
                status[label] = "skipped"
                continue
            if label == "B6_clinical_associations" and (
                not ctx.conflict_file
                or not (io.DATA_DIR / ctx.conflict_file).is_file()
            ):
                 status[label] = "skipped"
                 continue
            print(f"  {label:<26}", end="", flush=True)
            start = time.time()
            proc = subprocess.run([sys.executable, str(HERE / script), *args],
                                  cwd=str(io.ROOT), capture_output=True, text=True)
            output = proc.stdout + (("\n" + proc.stderr) if proc.stderr else "")
            log.write(f"\n--- {label} ({script}) ---\n{output}")
            if proc.returncode == 0:
                status[label] = "ok"
                print(f"ok ({time.time() - start:.1f} s)")
            else:
                status[label] = "FAILED"
                print("FAILED; last lines:")
                print("    " + "\n    ".join(output.strip().splitlines()[-6:]))
                if stop_on_error:
                    break
    print(f"  log: {log_path.relative_to(io.ROOT)}")
    return status


def summarize(ctx, status):
    """Headline numbers per run, for a quick side-by-side look."""
    row = {"tag": ctx.name, "micro": ctx.micro_file, "conflict": ctx.conflict_file,
           "simon_qc": ctx.simon_qc, "run_at": time.strftime("%Y-%m-%d %H:%M"),
           "all_steps_ok": all(v == "ok" for v in status.values())}
    if ctx.csv("reorganization").is_file():
        reorg = pd.read_csv(ctx.csv("reorganization"))
        row["K"] = reorg.filter(like="delta_cov_").shape[1]
        row["n_reorganization"] = len(reorg)
        row["reorg_jsd_mean"] = reorg["reorganization_jsd"].mean()
    for stem, prefix in [("rq1_stats", "rq1"), ("rq1_group_delta", "rq1_group_delta"),
                         ("rq3_stats", "rq3")]:
        if ctx.csv(stem).is_file():
            table = pd.read_csv(ctx.csv(stem))
            row[f"{prefix}_tests"] = len(table)
            row[f"{prefix}_fdr_sig"] = int(table["reject_fdr"].sum())
    if ctx.csv("b5_subject_clusters").is_file():
        row["b5_n_clusters"] = pd.read_csv(ctx.csv("b5_subject_clusters"))["cluster"].nunique()
    return row


def update_summary(rows):
    path = io.RESULTS_ROOT / "summary_all_runs.csv"
    new = pd.DataFrame(rows)
    if path.is_file():
        old = pd.read_csv(path)
        new = pd.concat([old[~old["tag"].astype(str).isin(new["tag"])], new],
                        ignore_index=True)
    new.to_csv(path, index=False)
    return path, new


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--micro", help="run only this microstate CSV (ignores RUNS)")
    parser.add_argument("--conflict", help="conflict CSV for --micro")
    parser.add_argument("--tag", help="tag for --micro (default: from file name)")
    parser.add_argument("--no-simon-qc", action="store_true", help="with --micro")
    parser.add_argument("--only", nargs="+", help="run only these tags from RUNS")
    parser.add_argument("--stop-on-error", action="store_true")
    args = parser.parse_args()

    if args.micro:
        runs = [{"tag": args.tag, "micro": args.micro, "conflict": args.conflict,
                 "simon_qc": not args.no_simon_qc}]
    else:
        runs = RUNS
    contexts = [io.make_context(r["micro"], r.get("conflict"), r.get("tag"),
                                r.get("simon_qc", True), verbose=False) for r in runs]
    if args.only and not args.micro:
        contexts = [c for c in contexts if c.name in args.only]
    tags = [c.name for c in contexts]
    if len(tags) != len(set(tags)):
        raise ValueError(f"Duplicate tags in RUNS: {tags}")
    if not contexts:
        raise SystemExit("Nothing to run: check RUNS / --only")

    rows, all_status = [], {}
    for ctx in contexts:
        print(f"\n=== run {ctx.name}: data/{ctx.micro_file} + data/{ctx.conflict_file} "
              f"(Simon QC {'on' if ctx.simon_qc else 'OFF'}) ===")
        problems = preflight(ctx)
        if problems:
            print("  input problems, run skipped:\n   - " + "\n   - ".join(problems))
            all_status[ctx.name] = {label: "not run" for label, _ in STEPS}
            if args.stop_on_error:
                break
            continue
        status = run_one(ctx, args.stop_on_error)
        all_status[ctx.name] = status
        rows.append(summarize(ctx, status))
        if args.stop_on_error and "FAILED" in status.values():
            break

    print("\n=== Step status ===")
    print(pd.DataFrame(all_status).T.fillna("-").to_string())
    if rows:
        path, table = update_summary(rows)
        print(f"\nSummary of all runs so far: {path.relative_to(io.ROOT)}")
        print(table.to_string(index=False))


if __name__ == "__main__":
    main()
