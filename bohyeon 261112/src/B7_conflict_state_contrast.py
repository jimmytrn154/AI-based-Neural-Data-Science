"""B7: Does the conflict effect (Incongruent - Congruent) on microstate
probability differ between state C and state E?

Input columns (wide, one row per subject):
    simon_conflict_<window>_<state>   e.g. simon_conflict_n2_C, simon_conflict_theta_E
Each value is P(state | Incongruent) - P(state | Congruent) within the window.

Analyses, run separately for each window (n2, theta):
  M1  Long-format OLS with subject fixed effects
          delta ~ C(state, ref=C) + C(subject_id)
      -> state[T.E] = mean(delta_E - delta_C)   (primary, = paired t-test)
  M2  Difference-score OLS (HC3 SE)
          (delta_E - delta_C) ~ C(group, ref=Control) [+ age + sex]
      -> PD coefficient = state x group interaction
  Robustness: paired Wilcoxon (stats.paired_wilcoxon) on E vs C.
  BH-FDR across windows within each test family (stats.fdr_bh).

Note: the five state deltas sum to 0 within a subject (probabilities sum to 1),
so C and E are not independent; the paired/within-subject contrast handles this.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.formula.api as smf

from stats import ALPHA, fdr_bh, paired_wilcoxon

SRC_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SRC_DIR.parent
# microstates_task.csv only has whole-task (rest vs simon) features, with no
# Congruent/Incongruent split, so the trial-type contrast must come from conflict.csv.
DEFAULT_DATA = PROJECT_DIR / "data" / "conflict.csv"
DEFAULT_OUT = PROJECT_DIR / "results" / "B7_conflict_state_contrast.csv"

STATE_REF, STATE_TEST = "C", "E"
COL_PATTERN = re.compile(r"^simon_conflict_(?P<window>[A-Za-z0-9]+)_(?P<state>[A-E])$")


# ----------------------------------------------------------------------------
# Data preparation
# ----------------------------------------------------------------------------
def find_windows(df: pd.DataFrame) -> list[str]:
    """Windows that have both the C and E conflict columns."""
    found: dict[str, set[str]] = {}
    for col in df.columns:
        m = COL_PATTERN.match(col)
        if m:
            found.setdefault(m["window"], set()).add(m["state"])
    windows = [w for w, s in found.items() if {STATE_REF, STATE_TEST} <= s]
    if not windows:
        raise KeyError(
            "No 'simon_conflict_<window>_C/E' columns found. "
            f"Available columns: {list(df.columns)}"
        )
    return windows


def add_group(df: pd.DataFrame, id_col: str = "subject_id") -> pd.DataFrame:
    """Use an existing 'group' column, else derive it from the ID (sub-P / sub-C)."""
    out = df.copy()
    if "group" not in out:
        prefix = out[id_col].astype(str).str.extract(r"sub-([PC])", expand=False)
        out["group"] = prefix.map({"P": "PD", "C": "Control"})
    return out


def merge_clinical(df: pd.DataFrame, clinical_path: Path | None,
                   id_col: str = "subject_id") -> pd.DataFrame:
    """Optionally merge age/gender for covariate adjustment (analysis-time only)."""
    if clinical_path is None:
        return df
    clinical = pd.read_csv(clinical_path)
    keep = [c for c in (id_col, "age", "gender") if c in clinical]
    return df.drop(columns=[c for c in ("age", "gender") if c in df and c in keep],
                   errors="ignore").merge(clinical[keep], on=id_col, how="left")


def to_long(df: pd.DataFrame, window: str, id_col: str = "subject_id") -> pd.DataFrame:
    """Subject x state long table restricted to states C and E."""
    cols = {f"simon_conflict_{window}_{s}": s for s in (STATE_REF, STATE_TEST)}
    long = (df[[id_col, "group", *cols]]
            .melt(id_vars=[id_col, "group"], value_vars=list(cols),
                  var_name="column", value_name="delta"))
    long["state"] = long["column"].map(cols)
    long["delta"] = pd.to_numeric(long["delta"], errors="coerce")
    long = long.dropna(subset=["delta"])
    # keep subjects with both states (paired design)
    complete = long.groupby(id_col)["state"].nunique().eq(2)
    return long[long[id_col].isin(complete[complete].index)].drop(columns="column")


# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------
def m1_fixed_effects(long: pd.DataFrame, id_col: str = "subject_id") -> dict:
    """delta ~ state + subject FE; state[T.E] is the E-minus-C difference."""
    fit = smf.ols(
        f'delta ~ C(state, Treatment(reference="{STATE_REF}")) + C({id_col})',
        data=long,
    ).fit()
    term = next(t for t in fit.params.index if f"[T.{STATE_TEST}]" in t)
    ci = fit.conf_int().loc[term]
    n_sub = long[id_col].nunique()
    wide = long.pivot(index=id_col, columns="state", values="delta")
    diff = wide[STATE_TEST] - wide[STATE_REF]
    return {
        "test": "M1_OLS_subjectFE", "term": f"state {STATE_TEST} - {STATE_REF}",
        "n": n_sub, "estimate": float(fit.params[term]),
        "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
        "t_stat": float(fit.tvalues[term]), "df": float(fit.df_resid),
        "p_value": float(fit.pvalues[term]),
        "cohens_dz": float(diff.mean() / diff.std(ddof=1)),
        "mean_delta_C": float(wide[STATE_REF].mean()),
        "mean_delta_E": float(wide[STATE_TEST].mean()),
    }


def m2_group_interaction(long: pd.DataFrame, id_col: str = "subject_id") -> dict:
    """(delta_E - delta_C) ~ group [+ age + gender], HC3 SE.

    The PD coefficient is the state x group interaction from the long model.
    """
    wide = long.pivot_table(index=[id_col, "group"], columns="state",
                            values="delta").reset_index()
    # rename: a column named "C" would shadow patsy's C() function
    wide = wide.rename(columns={STATE_REF: f"delta_{STATE_REF}",
                                STATE_TEST: f"delta_{STATE_TEST}"})
    wide.columns.name = None
    wide["diff_EC"] = wide[f"delta_{STATE_TEST}"] - wide[f"delta_{STATE_REF}"]
    wide = wide[wide["group"].isin(["PD", "Control"])]

    rhs = 'C(group, Treatment(reference="Control"))'
    covariates = []
    if {"age", "gender"} <= set(long.columns):
        extra = long.drop_duplicates(id_col)[[id_col, "age", "gender"]]
        wide = wide.merge(extra, on=id_col, how="left")
        wide["age"] = pd.to_numeric(wide["age"], errors="coerce")
        wide = wide.dropna(subset=["age"])
        wide = wide[wide["gender"].isin(["F", "M"])]
        rhs += " + age + C(gender)"
        covariates = ["age", "gender"]

    fit = smf.ols(f"diff_EC ~ {rhs}", data=wide).fit(cov_type="HC3")
    term = next(t for t in fit.params.index if "[T.PD]" in t)
    ci = fit.conf_int().loc[term]
    return {
        "test": "M2_OLS_diff_by_group", "term": "(E - C) x group [PD - Control]",
        "n": len(wide), "estimate": float(fit.params[term]),
        "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
        "t_stat": float(fit.tvalues[term]), "df": float(fit.df_resid),
        "p_value": float(fit.pvalues[term]),
        "covariates": "+".join(covariates) or "none",
        "mean_diff_PD": float(wide.loc[wide.group == "PD", "diff_EC"].mean()),
        "mean_diff_Control": float(wide.loc[wide.group == "Control", "diff_EC"].mean()),
    }


def wilcoxon_check(long: pd.DataFrame, id_col: str = "subject_id") -> dict:
    """Nonparametric robustness check, reusing stats.paired_wilcoxon (E - C)."""
    wide = long.pivot(index=id_col, columns="state", values="delta")
    res = paired_wilcoxon(wide, rest_col=STATE_REF, simon_col=STATE_TEST)
    return {"test": "Wilcoxon_paired", "term": f"state {STATE_TEST} - {STATE_REF}",
            "n": res["n_pairs"], "estimate": res["median_delta"],
            "p_value": res["p_value"], "w_stat": res["w_stat"],
            "rank_biserial": res["rank_biserial"], "cohens_dz": res["cohens_dz"]}


# ----------------------------------------------------------------------------
# Figures: save PNG files beside the output CSV
# ----------------------------------------------------------------------------
def plot_figures(df: pd.DataFrame, results: pd.DataFrame, out_dir: Path) -> list[Path]:
    saved = []
    out_dir.mkdir(parents=True, exist_ok=True)
    for window in find_windows(df):
        long = to_long(df, window)
        wide = long.pivot(index="subject_id", columns="state", values="delta")
        subset = results[results["window"] == window]

        # Figure 1: within-subject C vs E comparison
        fig, ax = plt.subplots(figsize=(5.6, 5.0))
        for _, subject in wide.iterrows():
            ax.plot([0, 1], [subject["C"], subject["E"]],
                    color="0.72", linewidth=0.9, alpha=0.65, zorder=1)
        ax.scatter(np.zeros(len(wide)), wide["C"], s=14, color="#30758a",
                   alpha=0.55, zorder=2)
        ax.scatter(np.ones(len(wide)), wide["E"], s=14, color="#b36c49",
                   alpha=0.55, zorder=2)
        means = [wide["C"].mean(), wide["E"].mean()]
        ax.plot([0, 1], means, "o-", color="black", linewidth=2.5,
                markersize=8, zorder=3, label="Mean")
        m1 = subset[subset["test"] == "M1_OLS_subjectFE"]
        if not m1.empty:
            ax.set_title(f"{window.upper()} | C vs E (FDR p={m1.iloc[0]['p_fdr']:.3g})")
        ax.set_xticks([0, 1], ["State C", "State E"])
        ax.set_ylabel("Conflict effect: P(Incongruent) - P(Congruent)")
        ax.axhline(0, color="0.45", linestyle="--", linewidth=0.8)
        ax.set_xlim(-0.35, 1.35)
        ax.legend(frameon=False)
        fig.tight_layout()
        path = out_dir / f"B7_{window}_paired_C_vs_E.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)

        # Figure 2: subject-level (E-C) contrast across groups
        group_wide = long.pivot(index="subject_id", columns="state", values="delta")
        group_wide["diff_EC"] = group_wide["E"] - group_wide["C"]
        group_wide = group_wide.join(long.drop_duplicates("subject_id")
                                    .set_index("subject_id")["group"])
        fig, ax = plt.subplots(figsize=(5.6, 5.0))
        rng = np.random.default_rng(42)
        for i, (group, color) in enumerate([("Control", "#30758a"), ("PD", "#b36c49")]):
            values = group_wide.loc[group_wide["group"] == group, "diff_EC"].dropna().to_numpy()
            if len(values) == 0:
                continue
            jitter = rng.uniform(-0.10, 0.10, size=len(values))
            ax.scatter(np.full(len(values), i) + jitter, values,
                       alpha=0.65, s=24, color=color)
            mean = values.mean()
            sem = values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else 0
            ax.errorbar(i, mean, yerr=1.96 * sem, fmt="D", markersize=7,
                        color="black", capsize=5, linewidth=1.8, zorder=5)
        m2 = subset[subset["test"] == "M2_OLS_diff_by_group"]
        if not m2.empty:
            adjusted = m2.iloc[0]["covariates"] != "none"
            qualifier = "adjusted" if adjusted else "unadjusted"
            ax.set_title(f"{window.upper()} | PD vs Control ({qualifier} FDR p={m2.iloc[0]['p_fdr']:.3g})")
        ax.set_xticks([0, 1], ["Control", "PD"])
        ax.set_xlim(-0.5, 1.5)
        ax.set_ylabel("Conflict contrast: Delta(E) - Delta(C)")
        ax.axhline(0, color="0.45", linestyle="--", linewidth=0.8)
        fig.tight_layout()
        path = out_dir / f"B7_{window}_PD_vs_Control.png"
        fig.savefig(path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        saved.append(path)
    return saved


# ----------------------------------------------------------------------------
def run(data_path: Path, out_path: Path, clinical_path: Path | None) -> pd.DataFrame:
    df = pd.read_csv(data_path)
    df = merge_clinical(add_group(df), clinical_path)
    windows = find_windows(df)

    rows = []
    for window in windows:
        long = to_long(df, window)
        if {"age", "gender"} <= set(df.columns):
            long = long.merge(df[["subject_id", "age", "gender"]], on="subject_id")
        for func in (m1_fixed_effects, m2_group_interaction, wilcoxon_check):
            row = func(long)
            row["window"] = window
            rows.append(row)

    results = pd.DataFrame(rows)
    # FDR across windows within each test (family = test type)
    results = fdr_bh(results, p_col="p_value", family_col="test", alpha=ALPHA)
    first = ["window", "test", "term", "n", "estimate", "ci_low", "ci_high",
             "p_value", "p_fdr", "reject_fdr"]
    results = results[first + [c for c in results.columns if c not in first]]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(out_path, index=False)
    figure_paths = plot_figures(df, results, out_path.parent)
    for path in figure_paths:
        print(f"Saved figure: {path}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--clinical", type=Path, default=None,
                        help="optional CSV with subject_id, age, gender")
    args = parser.parse_args()

    results = run(args.data, args.out, args.clinical)
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", "{:.4g}".format):
        print(results[["window", "test", "term", "n", "estimate",
                       "ci_low", "ci_high", "p_value", "p_fdr", "reject_fdr"]])
    print(f"\nSaved: {args.out}")


if __name__ == "__main__":
    main()