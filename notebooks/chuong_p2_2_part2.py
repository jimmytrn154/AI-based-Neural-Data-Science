"""P2.2 part 2: do frames that C loses become E?

Definitions are fixed in ``config.yaml`` -> ``phase2.probability_dependence_part2``
(set Nov 18, before any result was computed). Part 1 (coverage) cannot show
this, because back-fitting gives each sample exactly one state.

- Spatial correlation over time: per recording, mean |r| of all labelled
  samples with each template. Not competitive: a sample counts for every
  template, so E can rise without C falling. GFP peaks only as a sensitivity.
- Second-best labels: the runner-up of a sample is the template with the
  highest |r| among the four other than its primary (smoothed) label. If
  C frames turn into E, E samples in Simon should more often have C as the
  runner-up (and C samples E).

Paired Wilcoxon etc. from Bohyeon's verified ``stats.py`` (as part 1); one BH
family "P2.2 part 2" of 7 tests per run. Sensitivity: blink-cleaned Rest.
Correlations reproduce pycrostates' assignments exactly (argmax of |r| ==
``predict`` labels without smoothing or segment rejection; checked on sub-C001).

Run from the repo root: ``python notebooks/chuong_p2_2_part2.py``
"""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from src import io as pio, microstates as ms, preprocess as pp

cfg = pio.load_config()
FEAT = pio.results_dir(cfg, "features")
FIG = pio.results_dir(cfg, "figures")
spec = importlib.util.spec_from_file_location("bstats", ROOT / "bohyeon 261112/src/stats.py")
bstats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bstats)

K = cfg["microstates"]["primary_k"]
TEMPLATES = ms.load_templates(K, cfg)
NAMES = list(TEMPLATES.cluster_names)
assert NAMES == list("ABCDE")
CH = TEMPLATES.info["ch_names"]
MAPS = TEMPLATES.cluster_centers_
iC, iE = NAMES.index("C"), NAMES.index("E")
N_BOOT = 2000


def corr(inst_data_ch_names, data):
    assert list(inst_data_ch_names) == list(CH)
    return ms.spatial_corr(data.T, MAPS)  # (n_samples, 5)


def measures(raw, labels, prefix):
    """Mean |r| per template, GFP-peak mean |r|, runner-up shares."""
    picked = raw.copy().pick(CH)
    r = corr(picked.ch_names, picked.get_data())
    assert labels.size == r.shape[0]
    ok = labels >= 0
    row = {f"{prefix}_n_labelled": int(ok.sum())}
    for j, s in enumerate(NAMES):
        row[f"{prefix}_meanr_{s}"] = r[ok, j].mean()
    peaks = ms.extract_peaks(picked, cfg)
    rp = corr(peaks.info["ch_names"], peaks.get_data())
    for j, s in enumerate(NAMES):
        row[f"{prefix}_peakmeanr_{s}"] = rp[:, j].mean()
    lab, rr = labels[ok], r[ok].copy()
    rr[np.arange(lab.size), lab] = -np.inf
    runner = rr.argmax(1)
    for i, s in enumerate(NAMES):
        sel = lab == i
        row[f"{prefix}_n_{s}"] = int(sel.sum())
        for j, t in enumerate(NAMES):
            if i != j:
                row[f"{prefix}_ru_{s}_{t}"] = (runner[sel] == j).mean() if sel.any() else np.nan
    return row


def subject(sid, blink_clean):
    row = {"subject_id": sid}
    for task in ("Rest", "Simon"):
        raw = pp.preprocess_raw(pio.load_continuous(sid, task, cfg), cfg)
        labels = np.load(ms.labels_path(sid, task, cfg))
        if task == "Rest" and blink_clean:
            frac = pp.annotate_blinks(raw, cfg)
            new = ms.backfit(raw, TEMPLATES, cfg).labels
            if frac == 0:  # nothing marked: must reproduce the primary labels
                assert np.array_equal(new, labels), sid
            labels = new
            row["rest_blink_bad_frac"] = frac
        row.update(measures(raw, labels, task.lower()))
    return row


def per_subject():
    prim = pd.read_csv(FEAT / "microstates_task.csv")
    cov = [f"{t}_cov_{s}" for t in ("rest", "simon") for s in NAMES]
    sids = prim.dropna(subset=cov).subject_id.tolist()
    rest, _ = pp.select_subjects("Rest", cfg)
    flagged = set(rest.loc[rest.blink_flag, "subject_id"]) & set(sids)
    print(f"Subjects: {len(sids)}; Rest recordings to blink-clean: {len(flagged)}")
    jobs = [(sid, False) for sid in sids] + [(sid, True) for sid in sorted(flagged)]
    rows = Parallel(n_jobs=8)(delayed(subject)(*j) for j in jobs)
    primary = pd.DataFrame(rows[:len(sids)])
    clean = pd.DataFrame(rows[len(sids):]).set_index("subject_id")
    blink = primary.set_index("subject_id").copy()
    rest_cols = [c for c in clean.columns if c.startswith("rest_") and c != "rest_blink_bad_frac"]
    blink.loc[clean.index, rest_cols] = clean[rest_cols]
    return primary, blink.reset_index(), clean


def boot_ci(x, rng):
    means = [rng.choice(x, size=len(x), replace=True).mean() for _ in range(N_BOOT)]
    return np.percentile(means, [2.5, 97.5])


def analyse(d, run):
    rng = np.random.default_rng(cfg["random_seed"])
    tests = [("mean_abs_r", s, f"meanr_{s}") for s in NAMES] \
        + [("runner_up_share", "E -> runner-up C", "ru_E_C"), ("runner_up_share", "C -> runner-up E", "ru_C_E")]
    sens = [("mean_abs_r_gfp_peaks", s, f"peakmeanr_{s}") for s in NAMES]
    rows = []
    for family, group in (("P2.2 part 2", tests), ("sensitivity", sens)):
        part = []
        for analysis, label, m in group:
            res = bstats.paired_wilcoxon(d, f"rest_{m}", f"simon_{m}")
            lo, hi = boot_ci((d[f"simon_{m}"] - d[f"rest_{m}"]).to_numpy(), rng)
            part.append({"run": run, "family": family, "analysis": analysis, "state": label,
                         "mean_rest": d[f"rest_{m}"].mean(), "mean_simon": d[f"simon_{m}"].mean(),
                         **res, "ci_low": lo, "ci_high": hi})
        part = pd.DataFrame(part)
        rows.append(bstats.fdr_bh(part) if family == "P2.2 part 2" else part)
    # Full runner-up distribution, descriptive
    desc = []
    for s in NAMES:
        for t in NAMES:
            if s != t:
                x = d[f"simon_ru_{s}_{t}"] - d[f"rest_ru_{s}_{t}"]
                desc.append({"run": run, "family": "descriptive", "analysis": "runner_up_distribution",
                             "state": f"{s} -> runner-up {t}", "n_pairs": int(x.notna().sum()),
                             "mean_rest": d[f"rest_ru_{s}_{t}"].mean(), "mean_simon": d[f"simon_ru_{s}_{t}"].mean(),
                             "mean_delta": x.mean(), "median_delta": x.median()})
    return pd.concat(rows + [pd.DataFrame(desc)], ignore_index=True)


def figure(res, d):
    r = res[res.run == "primary"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    q = r[r.analysis == "mean_abs_r"].set_index("state").loc[NAMES]
    cols = ["#c0504d" if p < 0.05 else "#999999" for p in q.p_fdr]
    ax[0].bar(NAMES, q.mean_delta, color=cols, yerr=[q.mean_delta - q.ci_low, q.ci_high - q.mean_delta], capsize=3)
    ax[0].axhline(0, color="k", lw=0.8); ax[0].set_ylabel("Simon − Rest, mean |r|")
    ax[0].set_title("(a) Fit of each template, all labelled samples")
    for i, p in enumerate(q.p_fdr):
        y = q.ci_high.iloc[i] if q.mean_delta.iloc[i] >= 0 else q.ci_low.iloc[i]
        ax[0].annotate(f"q={p:.3g}", (i, y), ha="center", va="bottom" if y >= 0 else "top", fontsize=7)
    # (b) runner-up distribution of E samples, Rest vs Simon
    others = [t for t in NAMES if t != "E"]
    x = np.arange(len(others))
    for off, t, c in ((-0.18, "rest", "#aaaaaa"), (0.18, "simon", "#3b6fb6")):
        ax[1].bar(x + off, [100 * d[f"{t}_ru_E_{o}"].mean() for o in others], 0.36, color=c, label=t.capitalize())
    ax[1].set_xticks(x, others); ax[1].set_ylabel("% of E samples"); ax[1].legend(fontsize=8)
    ax[1].set_title("(b) Runner-up template of E samples")
    # (c) runner-up distribution of C samples
    others = [t for t in NAMES if t != "C"]
    for off, t, c in ((-0.18, "rest", "#aaaaaa"), (0.18, "simon", "#3b6fb6")):
        ax[2].bar(x + off, [100 * d[f"{t}_ru_C_{o}"].mean() for o in others], 0.36, color=c, label=t.capitalize())
    ax[2].set_xticks(x, others); ax[2].set_ylabel("% of C samples"); ax[2].legend(fontsize=8)
    ax[2].set_title("(c) Runner-up template of C samples")
    n = int(q.n_pairs.iloc[0])
    t = r[r.analysis == "runner_up_share"].set_index("state")
    fig.suptitle(f"P2.2 part 2 (N = {n}; mean ± 95% bootstrap CI; red: FDR q < 0.05, family of 7). "
                 f"E→C runner-up q = {t.loc['E -> runner-up C', 'p_fdr']:.3g}; "
                 f"C→E runner-up q = {t.loc['C -> runner-up E', 'p_fdr']:.3g}")
    fig.tight_layout()
    fig.savefig(FIG / "chuong_p2_2_part2.png", dpi=200)


if __name__ == "__main__":
    primary, blink, clean = per_subject()
    primary.to_csv(FEAT / "p2_2_part2_subject_measures.csv", index=False)
    print("Blink-cleaned Rest: recordings with time marked:", int((clean.rest_blink_bad_frac > 0).sum()))
    res = pd.concat([analyse(primary, "primary"), analyse(blink, "blinkclean_rest")], ignore_index=True)
    res.to_csv(FEAT / "p2_2_part2_probability_dependence.csv", index=False)
    figure(res, primary)
    pd.set_option("display.width", 220)
    show = ["run", "family", "analysis", "state", "n_pairs", "mean_rest", "mean_simon", "mean_delta",
            "ci_low", "ci_high", "rank_biserial", "cohens_dz", "p_value", "p_fdr", "frac_positive"]
    print(res[[c for c in show if c in res]].round(4).to_string())
