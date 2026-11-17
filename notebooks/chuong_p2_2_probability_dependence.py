"""P2.2 part 1: is E-up more than the compositional consequence of C-down?

Definitions are fixed in ``config.yaml`` -> ``phase2.probability_dependence``
(set Nov 17, before any result was computed). Coverage is in percent and
sums to 100 per recording, so a drop in C must raise some other state.

- Redistribution null: if the time C loses were shared out in proportion to
  Rest coverage, X would change by -dC * X_rest / (100 - C_rest). The excess
  (observed - expected) equals (100 - C_simon) * change in X's share, so it is
  reported with the share test, not tested separately.
- Share of non-C time: X / (100 - C), Simon vs Rest, for A, B, D, E.
- Log-ratios: centred log-ratio (CLR) change for A-E, and log(E / X) for
  X = A, B, D.

Paired Wilcoxon, rank-biserial r, Cohen's d_z and bootstrap CI of the mean
use Bohyeon's verified ``stats.py`` (as the RQ1 tables); one BH family "P2.2"
of 12 tests. Sensitivity: Rest coverage from the blink-cleaned file.

Run from the repo root: ``python notebooks/chuong_p2_2_probability_dependence.py``
"""

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load(open(ROOT / "config.yaml"))
FEAT = ROOT / cfg["paths"]["results"] / "features"
FIG = ROOT / cfg["paths"]["results"] / "figures"

spec = importlib.util.spec_from_file_location("bstats", ROOT / "bohyeon 261112/src/stats.py")
bstats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bstats)

STATES = list("ABCDE")
OTHERS = ["A", "B", "D", "E"]
PAIRWISE = ["A", "B", "D"]
N_BOOT = 2000


def load(run):
    prim = pd.read_csv(FEAT / "microstates_task.csv")
    cols = {t: [f"{t}_cov_{s}" for s in STATES] for t in ("rest", "simon")}
    df = prim[["subject_id"] + cols["rest"] + cols["simon"]]
    if run == "blinkclean_rest":
        clean = pd.read_csv(FEAT / "microstates_rest_blinkclean.csv")[["subject_id"] + cols["rest"]]
        df = df.drop(columns=cols["rest"]).merge(clean, on="subject_id", how="left")
    df = df.dropna().reset_index(drop=True)
    out = pd.DataFrame({"subject_id": df.subject_id})
    for t in ("rest", "simon"):
        cov = df[cols[t]].to_numpy(float)
        assert np.allclose(cov.sum(1), 100, atol=1e-6) and (cov > 0).all()
        c = df[f"{t}_cov_C"]
        for s in STATES:
            out[f"{t}_cov_{s}"] = df[f"{t}_cov_{s}"]
        for s in OTHERS:
            out[f"{t}_share_{s}"] = 100 * df[f"{t}_cov_{s}"] / (100 - c)
        logc = np.log(cov)
        clr = logc - logc.mean(1, keepdims=True)
        for i, s in enumerate(STATES):
            out[f"{t}_clr_{s}"] = clr[:, i]
        for s in PAIRWISE:
            out[f"{t}_logE{s}"] = np.log(df[f"{t}_cov_E"] / df[f"{t}_cov_{s}"])
    dC = out.simon_cov_C - out.rest_cov_C
    for s in OTHERS:
        out[f"expected_d_{s}"] = -dC * out[f"rest_cov_{s}"] / (100 - out.rest_cov_C)
        out[f"excess_d_{s}"] = (out[f"simon_cov_{s}"] - out[f"rest_cov_{s}"]) - out[f"expected_d_{s}"]
    return out


def boot_ci(x, rng):
    means = [rng.choice(x, size=len(x), replace=True).mean() for _ in range(N_BOOT)]
    return np.percentile(means, [2.5, 97.5])


def analyse(run):
    d = load(run)
    rng = np.random.default_rng(cfg["random_seed"])
    tests = [("share", s, f"share_{s}") for s in OTHERS] \
        + [("clr", s, f"clr_{s}") for s in STATES] \
        + [("log_ratio", f"E/{s}", f"logE{s}") for s in PAIRWISE]
    rows = []
    for analysis, label, m in tests:
        res = bstats.paired_wilcoxon(d, f"rest_{m}", f"simon_{m}")
        diff = (d[f"simon_{m}"] - d[f"rest_{m}"]).to_numpy()
        lo, hi = boot_ci(diff, rng)
        rows.append({"run": run, "family": "P2.2", "analysis": analysis, "state": label,
                     "mean_rest": d[f"rest_{m}"].mean(), "mean_simon": d[f"simon_{m}"].mean(),
                     **res, "ci_low": lo, "ci_high": hi})
    stats = bstats.fdr_bh(pd.DataFrame(rows))
    # Redistribution: observed vs expected change (descriptive; same null as share)
    desc = []
    for s in OTHERS:
        for kind, x in (("observed_d", d[f"simon_cov_{s}"] - d[f"rest_cov_{s}"]),
                        ("expected_d", d[f"expected_d_{s}"]), ("excess_d", d[f"excess_d_{s}"])):
            lo, hi = boot_ci(x.to_numpy(), rng)
            desc.append({"run": run, "family": "descriptive", "analysis": f"redistribution_{kind}",
                         "state": s, "n_pairs": len(x), "mean_delta": x.mean(),
                         "median_delta": x.median(), "ci_low": lo, "ci_high": hi,
                         "frac_positive": (x > 0).mean()})
    dc = d.simon_cov_C - d.rest_cov_C
    lo, hi = boot_ci(dc.to_numpy(), rng)
    desc.append({"run": run, "family": "descriptive", "analysis": "redistribution_observed_d",
                 "state": "C", "n_pairs": len(dc), "mean_delta": dc.mean(),
                 "median_delta": dc.median(), "ci_low": lo, "ci_high": hi,
                 "frac_positive": (dc > 0).mean()})
    return pd.concat([stats, pd.DataFrame(desc)], ignore_index=True), d


def figure(res):
    r = res[res.run == "primary"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    # (a) observed vs expected change in coverage
    x = np.arange(len(OTHERS))
    for k, (kind, col, off) in enumerate((("observed_d", "#3b6fb6", -0.18), ("expected_d", "#aaaaaa", 0.18))):
        q = r[r.analysis == f"redistribution_{kind}"].set_index("state").loc[OTHERS]
        ax[0].bar(x + off, q.mean_delta, 0.36, color=col,
                  yerr=[q.mean_delta - q.ci_low, q.ci_high - q.mean_delta], capsize=3,
                  label="observed" if kind == "observed_d" else "expected if C's time is shared\nin proportion to Rest coverage")
    ax[0].set_xticks(x, OTHERS); ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_ylabel("Simon − Rest coverage (points)"); ax[0].legend(fontsize=8, loc="lower right"); ax[0].set_ylim(-3.4, None)
    ax[0].set_title("(a) Where C's lost time goes")
    # (b) share of non-C time; (c) log-ratios
    for a, an, labels, title, ylab in (
            (ax[1], ["share"], OTHERS, "(b) Share of non-C time, X / (100 − C)", "Simon − Rest (points of share)"),
            (ax[2], ["clr", "log_ratio"], STATES + [f"E/{s}" for s in PAIRWISE], "(c) Log-ratios", "Simon − Rest (log units)")):
        q = r[r.analysis.isin(an)].set_index("state").loc[labels]
        cols = ["#c0504d" if p < 0.05 else "#999999" for p in q.p_fdr]
        a.bar(range(len(q)), q.mean_delta, color=cols,
              yerr=[q.mean_delta - q.ci_low, q.ci_high - q.mean_delta], capsize=3)
        names = labels if an == ["share"] else [f"CLR {s}" for s in STATES] + [f"log {s}" for s in labels[5:]]
        a.set_xticks(range(len(q)), names, rotation=0 if an == ["share"] else 45, ha="center" if an == ["share"] else "right")
        a.axhline(0, color="k", lw=0.8); a.set_title(title); a.set_ylabel(ylab)
        for i, p in enumerate(q.p_fdr):
            y = q.ci_high.iloc[i] if q.mean_delta.iloc[i] >= 0 else q.ci_low.iloc[i]
            a.annotate(f"q={p:.3g}", (i, y), ha="center", va="bottom" if y >= 0 else "top", fontsize=7)
    n = int(r.n_pairs.dropna().iloc[0])
    fig.suptitle(f"P2.2 probability dependence (N = {n}; mean ± 95% bootstrap CI; red: FDR q < 0.05, family of 12)")
    fig.tight_layout()
    fig.savefig(FIG / "chuong_p2_2_probability_dependence.png", dpi=200)


if __name__ == "__main__":
    out = []
    for run in ("primary", "blinkclean_rest"):
        res, d = analyse(run)
        out.append(res)
        if run == "primary":
            # Identity check: excess = (100 - C_simon) * change in share / 100
            for s in OTHERS:
                lhs = d[f"excess_d_{s}"]
                rhs = (100 - d.simon_cov_C) * (d[f"simon_share_{s}"] - d[f"rest_share_{s}"]) / 100
                assert np.allclose(lhs, rhs)
    res = pd.concat(out, ignore_index=True)
    res.to_csv(FEAT / "p2_2_probability_dependence.csv", index=False)
    figure(res)
    pd.set_option("display.width", 200)
    show = ["run", "analysis", "state", "n_pairs", "mean_rest", "mean_simon", "mean_delta", "ci_low",
            "ci_high", "rank_biserial", "cohens_dz", "p_value", "p_fdr", "frac_positive"]
    print(res[show].round(4).to_string())
