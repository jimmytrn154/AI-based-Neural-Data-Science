"""Shared input/output handling for the B2-B6 pipeline.

Every analysis script accepts the same arguments:
    --micro     microstate CSV in data/        (e.g. microstate_A.csv)
    --conflict  conflict CSV in data/          (e.g. conflict_A.csv; used by B6)
    --tag       name of the run                (default: taken from --micro,
                                                microstate_A.csv -> A)
    --no-simon-qc   keep subjects with trials_simon = 'exclude'

Outputs of a run with tag A:
    results/A/<original name>_A.csv     e.g. results/A/rq1_stats_A.csv
    figures/A/<original name>_A.png     e.g. figures/A/bohyeon_b6_correlations_A.png

Usually you do not call the scripts one by one: edit RUNS in
src/run_pipeline.py and run that file.
"""
import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

SRC_DIR = Path(__file__).resolve().parent
ROOT = SRC_DIR.parent
DATA_DIR = ROOT / "data"
RESULTS_ROOT = ROOT / "results"
FIGURES_ROOT = ROOT / "figures"

# Shared inputs (same for every run)
SUBJECTS_FILE = "subjects.csv"
BEHAVIOR_FILE = "behavior.csv"
# If a microstate CSV has only rest_* (or only simon_*) columns, the missing
# task is filled from this file. Set to None to raise an error instead.
FILL_MISSING_TASK_FROM = "microstates_task.csv"

TASKS = ("rest", "simon")

def as_csv_name(name):
    """'microstate_A' and 'microstate_A.csv' are both accepted."""
    name = Path(str(name)).name
    return name if name.lower().endswith(".csv") else f"{name}.csv"

def tag_from_filename(filename):
    """microstate_A.csv -> 'A'; microstates_task_K4.csv -> 'K4'; noext -> stem."""
    stem = Path(filename).stem
    return stem.rsplit("_", 1)[-1] if "_" in stem else stem

def clean_tag(tag):
    tag = re.sub(r"[^0-9A-Za-z_-]+", "-", str(tag)).strip("-_")
    if not tag:
        raise ValueError("Empty run tag")
    return tag

@dataclass
class Context:
    name: str            # the run tag, e.g. 'A'
    micro_file: str
    conflict_file: str
    simon_qc: bool
    results_dir: Path
    figures_dir: Path

    def csv(self, stem):
        """results/<tag>/<stem>_<tag>.csv"""
        return self.results_dir / f"{stem}_{self.name}.csv"

    def fig(self, stem):
        """figures/<tag>/<stem>_<tag>.png"""
        return self.figures_dir / f"{stem}_{self.name}.png"

    def file(self, stem, ext):
        return self.results_dir / f"{stem}_{self.name}.{ext}"


def add_arguments(parser):
    parser.add_argument("--micro", required=True,
                        help="microstate CSV in data/ (e.g. microstate_A.csv)")
    parser.add_argument("--conflict", default=None,
                        help="conflict CSV in data/ (e.g. conflict_A.csv); needed for B6")
    parser.add_argument("--tag", default=None,
                        help="run name for output folders/files (default: from --micro)")
    parser.add_argument("--no-simon-qc", action="store_true",
                        help="do not drop subjects with trials_simon = 'exclude'")
    return parser


def parse_args(description=None):
    parser = add_arguments(argparse.ArgumentParser(description=description))
    args = parser.parse_args()
    return make_context(args.micro, args.conflict, args.tag, not args.no_simon_qc)


def make_context(micro, conflict=None, tag=None, simon_qc=True, verbose=True):
    micro = as_csv_name(micro)
    conflict = as_csv_name(conflict) if conflict else None
    tag = clean_tag(tag or tag_from_filename(micro))
    results_dir = RESULTS_ROOT / tag
    figures_dir = FIGURES_ROOT / tag
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    ctx = Context(tag, micro, conflict, simon_qc, results_dir, figures_dir)
    if verbose:
        print(f"[run {tag}] micro = data/{micro} | conflict = "
              f"{'data/' + conflict if conflict else '-'} | Simon QC = "
              f"{'on' if simon_qc else 'OFF'} | outputs -> results/{tag}, figures/{tag}")
    return ctx


def _tasks_in(table):
    return {task for task in TASKS
            if any(column.startswith(f"{task}_") for column in table.columns)}


def read_unique(path):
    table = pd.read_csv(path)
    if "subject_id" not in table.columns:
        raise ValueError(f"{path} has no subject_id column")
    if table["subject_id"].duplicated().any():
        raise ValueError(f"Duplicate subject_id in {path}")
    return table


def load_microstates(ctx, verbose=True):
    """Read the run's microstate CSV; fill a missing task if configured."""
    path = DATA_DIR / ctx.micro_file
    if not path.is_file():
        raise FileNotFoundError(f"{path} not found")
    micro = read_unique(path)

    missing = set(TASKS) - _tasks_in(micro)
    if missing:
        if not FILL_MISSING_TASK_FROM or ctx.micro_file == FILL_MISSING_TASK_FROM:
            raise ValueError(f"{ctx.micro_file} has no {sorted(missing)} columns")
        base = read_unique(DATA_DIR / FILL_MISSING_TASK_FROM)
        columns = [c for c in base.columns if c.split("_", 1)[0] in missing]
        micro = micro.merge(base[["subject_id"] + columns], on="subject_id",
                            how="outer", validate="one_to_one")
        if verbose:
            print(f"NOTE: {ctx.micro_file} has no {'/'.join(sorted(missing))} columns; "
                  f"taken from {FILL_MISSING_TASK_FROM} (assumes identical state templates).")
    if verbose:
        print(f"Loaded microstates: {len(micro)} subjects")
    return micro


def detect_states(micro, verbose=True):
    """Detect microstate names from Rest and Simon coverage columns."""
    found = {task: set() for task in TASKS}
    for column in micro.columns:
        match = re.fullmatch(r"(rest|simon)_cov_([A-Za-z]+)", column)
        if match:
            found[match.group(1)].add(match.group(2))
    if found["rest"] != found["simon"]:
        raise ValueError(f"Rest/Simon state sets differ: {found}")
    states = sorted(found["rest"])
    if len(states) < 2:
        raise ValueError("Could not detect microstate columns (rest_cov_<S>)")
    if verbose:
        print(f"Microstates: K = {len(states)} ({states})")
    return states


def conflict_states(conflict):
    """Detect state names for each conflict window."""
    out = {}
    for window in ("n2", "theta"):
        pattern = rf"simon_conflict_{window}_([A-Za-z]+)"
        states = []
        for column in conflict.columns:
            match = re.fullmatch(pattern, column)
            if match:
                states.append(match.group(1))
        out[window] = sorted(states)
    return out


def apply_simon_qc(frame, ctx, column="trials_simon"):
    """Original rule: keep trials_simon in {'ok', 'partial'}; skipped if simon_qc is False."""
    if ctx.simon_qc:
        return frame[frame[column].isin(["ok", "partial"])].copy()
    print(f"NOTE: Simon QC filter disabled for run '{ctx.name}' "
          f"(trials_simon counts: {frame[column].value_counts(dropna=False).to_dict()})")
    return frame.copy()


def save_run_config(ctx, states, script):
    """Record which inputs/settings produced the files in results/<tag>/."""
    path = ctx.file("run_config", "json")
    config = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    config["run"] = {"tag": ctx.name, "micro_file": ctx.micro_file,
                     "conflict_file": ctx.conflict_file, "simon_qc": ctx.simon_qc,
                     "states": states}
    config.setdefault("scripts_run", [])
    if script not in config["scripts_run"]:
        config["scripts_run"].append(script)
    path.write_text(json.dumps(config, indent=2, ensure_ascii=False), encoding="utf-8")
