"""P2.7 remaining reruns: do the group and conflict results survive blink checks?

Rules in ``config.yaml`` -> ``phase2.blink_sensitivity_reruns`` (set Nov 17,
before any rerun result was computed):

- Rest-based group results (RQ1 mixed-model group terms, reorganization, RQ3
  group models and correlations) rerun with blink-cleaned Rest
  (``microstates_rest_blinkclean.csv``; Simon unchanged).
- P2-B7 / P2-B8 rerun without the Simon recordings that had no eye component
  removed (``nICremoved == 0``; 9 recordings, all in ``conflict.csv``).
- A result holds if its sign and its side of FDR q = 0.05 are unchanged.

The reruns use Bohyeon's scripts unchanged (``bohyeon 261112/src``), in a copy
whose ``data/`` folder holds ``subjects.csv`` (dataset), ``behavior.csv``,
``microstates_task.csv``, ``conflict.csv`` (Seonwoo's S6), the blink-cleaned
Rest file without ``rest_blink_bad_frac``, and ``conflict_simonblinkexcl.csv``.
From that copy:

    python src/run_pipeline.py --micro microstates_task.csv --conflict conflict.csv --tag task
    python src/run_pipeline.py --micro microstates_rest_blinkclean.csv --conflict conflict.csv --tag rest_blinkclean
    python src/B7_conflict_state_contrast.py --data data/conflict.csv --out results/task/B7_conflict_state_contrast.csv
    python src/B7_conflict_state_contrast.py --data data/conflict_simonblinkexcl.csv --out results/simonblinkexcl/B7_conflict_state_contrast.csv
    python src/B8_conflict_group_windows.py --conflict conflict.csv --tag task
    python src/B8_conflict_group_windows.py --conflict conflict_simonblinkexcl.csv --tag simonblinkexcl

The primary reruns reproduce the committed ``bohyeon 261112/results/task``
files (largest difference 1e-13; RQ1 mixed models: optimizer lbfgs instead of
powell, CI within 1.4e-4, p within 1e-5).

Then, from the repo root:
    python notebooks/chuong_p2_7_reruns.py --runs <copy>/results
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "features" / "p2_7_blink_sensitivity.csv"


def compare(prim, sens, keys, est, source, check):
    m = prim.merge(sens, on=keys, suffixes=("_primary", "_check"), validate="one_to_one")
    out = pd.DataFrame({
        "source": source, "check": check,
        "test": m[keys].astype(str).agg(" | ".join, axis=1),
        "n_primary": m.get("n_primary"), "n_check": m.get("n_check"),
        "estimate_primary": m[f"{est}_primary"], "estimate_check": m[f"{est}_check"],
        "p_primary": m["p_value_primary"], "p_check": m["p_value_check"],
        "q_primary": m["p_fdr_primary"], "q_check": m["p_fdr_check"],
    })
    same_sign = np.sign(out.estimate_primary) == np.sign(out.estimate_check)
    same_side = (out.q_primary < 0.05) == (out.q_check < 0.05)
    out["holds"] = same_sign & same_side
    return out


def main(runs):
    rows = []
    rest = ("rest_blinkclean", "blink-cleaned Rest")
    # RQ1 mixed models: group and task x group terms
    p, s = (pd.read_csv(runs / t / f"rq1_mixed_{t}.csv") for t in ("task", rest[0]))
    keep = lambda d: d[d.term.isin(["group_PD (at rest)", "task_x_group"])].rename(columns={"n_subjects": "n"})
    rows.append(compare(keep(p), keep(s), ["metric", "state", "term"], "estimate", "RQ1 mixed model", rest[1]))
    # RQ3 (B6): group models and partial correlations, incl. reorganization JSD
    p, s = (pd.read_csv(runs / t / f"rq3_stats_{t}.csv") for t in ("task", rest[0]))
    for d in (p, s):
        d["estimate"] = d.estimate_pd_minus_control.fillna(d.rho)
        d["outcome"] = d.outcome.fillna("-")
        d["subgroup"] = d.subgroup.fillna("-")
    rows.append(compare(p, s, ["analysis", "feature", "outcome", "subgroup"], "estimate", "RQ3 (B6)", rest[1]))
    # P2-B7 and P2-B8 without the Simon recordings lacking eye-component removal
    excl = ("simonblinkexcl", "Simon nICremoved == 0 excluded")
    p, s = (pd.read_csv(runs / t / "B7_conflict_state_contrast.csv") for t in ("task", excl[0]))
    rows.append(compare(p, s, ["window", "test", "term"], "estimate", "P2-B7", excl[1]))
    p, s = (pd.read_csv(runs / t / f"b8_adjusted_tests_{t}.csv") for t in ("task", excl[0]))
    for d in (p, s):
        d["group"] = d.group.fillna("-")
        d.rename(columns={"n_subjects": "n"}, inplace=True)
    rows.append(compare(p, s, ["family", "analysis", "window", "group", "feature"], "estimate", "P2-B8", excl[1]))
    res = pd.concat(rows, ignore_index=True)
    res.to_csv(OUT, index=False)

    pd.set_option("display.width", 250, "display.max_colwidth", 70)
    print(res.groupby("source").agg(tests=("holds", "size"), holds=("holds", "sum"),
                                    sig_primary=("q_primary", lambda q: (q < 0.05).sum()),
                                    sig_check=("q_check", lambda q: (q < 0.05).sum())))
    show = res[(res.q_primary < 0.05) | (res.q_check < 0.05) | ~res.holds]
    print(show.drop(columns=["check"]).round(4).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, required=True)
    main(ap.parse_args().runs)
