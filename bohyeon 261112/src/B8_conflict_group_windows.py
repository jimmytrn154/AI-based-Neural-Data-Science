"""B8: Age/sex-adjusted congruency effects on microstates C and E, by group and window.

Outcome per subject and window (from conflict.csv):
    delta_C      = P(C | Incongruent) - P(C | Congruent)
    delta_E      = P(E | Incongruent) - P(E | Congruent)
    D_E_minus_C  = delta_E - delta_C   (the C-vs-E contrast)

Three FDR families (Benjamini-Hochberg, corrected separately):
    vs_zero_12      Is each effect different from 0?  (group x window x feature,
                    one-sample t-test, unadjusted)
    group_6         Does PD differ from Control?       (window x feature,
                    OLS: feature ~ PD + age + male, HC3 SE)
    interaction_3   Does the group difference change between N2 and theta?
                    (feature, OLS: feature ~ PD * theta + age + male,
                    SE clustered by subject)

Design matrices are built explicitly (no formula strings), so adjusted means do
not depend on patsy/formulaic internals.  This avoids the
"'PandasData' object has no attribute 'design_info'" error in statsmodels >= 0.15.

Inputs : data/<--conflict>, data/subjects.csv
Outputs: results/<tag>/b8_*.csv, figures/<tag>/b8_*.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # save figures without needing a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"

WINDOWS = ("n2", "theta")
WINDOW_LABELS = {"n2": "N2\n200-350 ms", "theta": "Theta\n400-750 ms"}
FEATURES = ("delta_C", "delta_E", "D_E_minus_C")
GROUPS = ("Control", "PD")
COLORS = {"Control": "#246e9e", "PD": "#c16b39"}
SEED, N_BOOT = 42, 2000


# ----------------------------------------------------------------------------
# Data loading
# ----------------------------------------------------------------------------
def load_data(conflict_file: str, simon_qc: bool) -> pd.DataFrame:
    """Merge conflict effects with demographics; keep subjects with usable Simon data."""
    conflict = pd.read_csv(DATA_DIR / conflict_file)
    subjects = pd.read_csv(DATA_DIR / "subjects.csv")
    for name, table in (("conflict", conflict), ("subjects", subjects)):
        if table.subject_id.duplicated().any():
            raise ValueError(f"Duplicate subject_id in {name} file")

    data = conflict.merge(
        subjects[["subject_id", "group", "age", "gender", "trials_simon"]],
        on="subject_id", validate="one_to_one",
    )
    if simon_qc:
        data = data[data.trials_simon.isin(["ok", "partial"])]

    # Keep only rows that every model can use
    data = data.assign(age=pd.to_numeric(data.age, errors="coerce"))
    data = data[data.group.isin(GROUPS) & data.gender.isin(["F", "M"])]
    return data.dropna(subset=["age"]).reset_index(drop=True)


def to_long(data: pd.DataFrame) -> pd.DataFrame:
    """One row per subject x window with delta_C, delta_E and their difference."""
    parts = []
    for window in WINDOWS:
        part = data[["subject_id", "group", "age", "gender"]].copy()
        part["window"] = window
        for state in ("C", "E"):
            part[f"delta_{state}"] = pd.to_numeric(
                data[f"simon_conflict_{window}_{state}"], errors="coerce")
        part["D_E_minus_C"] = part.delta_E - part.delta_C
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


# ----------------------------------------------------------------------------
# Design matrices (explicit dummy coding; reference = Control, n2, female)
# ----------------------------------------------------------------------------
def design_group(df: pd.DataFrame) -> pd.DataFrame:
    """Columns: const, PD, age, male."""
    return pd.DataFrame({
        "const": 1.0,
        "PD": (df.group == "PD").astype(float),
        "age": df.age.astype(float),
        "male": (df.gender == "M").astype(float),
    }, index=df.index)


def design_interaction(df: pd.DataFrame) -> pd.DataFrame:
    """Columns: const, PD, theta, PD_x_theta, age, male."""
    pd_dummy = (df.group == "PD").astype(float)
    theta = (df.window == "theta").astype(float)
    return pd.DataFrame({
        "const": 1.0,
        "PD": pd_dummy,
        "theta": theta,
        "PD_x_theta": pd_dummy * theta,
        "age": df.age.astype(float),
        "male": (df.gender == "M").astype(float),
    }, index=df.index)


def adjusted_mean(fit, design_fn, population: pd.DataFrame, **fixed) -> dict:
    """Marginal mean: set group/window for everyone, keep their own age and sex,
    then average the predictions.  SE comes from the averaged design row."""
    grid = population.copy()
    for column, value in fixed.items():
        grid[column] = value
    w = design_fn(grid).mean(axis=0).to_numpy()
    estimate = float(w @ fit.params.to_numpy())
    se = float(np.sqrt(max(w @ fit.cov_params().to_numpy() @ w, 0.0)))
    return {"estimate": estimate, "ci_low": estimate - 1.96 * se,
            "ci_high": estimate + 1.96 * se}


def coefficient(fit, term: str) -> dict:
    """Estimate, 95% CI and p-value of one coefficient."""
    ci = fit.conf_int().loc[term]
    return {"estimate": float(fit.params[term]), "ci_low": float(ci.iloc[0]),
            "ci_high": float(ci.iloc[1]), "p_value": float(fit.pvalues[term])}


# ----------------------------------------------------------------------------
# Analyses
# ----------------------------------------------------------------------------
def test_vs_zero(x: np.ndarray, rng) -> dict:
    """Descriptives, bootstrap CI of the mean and one-sample t-test against 0."""
    n = len(x)
    if n >= 2:
        boot_means = rng.choice(x, size=(N_BOOT, n), replace=True).mean(axis=1)
        boot_low, boot_high = np.percentile(boot_means, [2.5, 97.5])
    else:
        boot_low = boot_high = np.nan
    return {"n": n, "mean": x.mean() if n else np.nan,
            "sd": x.std(ddof=1) if n >= 2 else np.nan,
            "boot_ci_low": boot_low, "boot_ci_high": boot_high,
            "p_value": stats.ttest_1samp(x, 0).pvalue if n >= 3 else np.nan}


def group_regression(window_df: pd.DataFrame, feature: str):
    """feature ~ PD + age + male (HC3).  Returns PD-Control test and adjusted means."""
    df = window_df.dropna(subset=[feature])
    fit = sm.OLS(df[feature].astype(float), design_group(df)).fit(cov_type="HC3")
    test = {"n_subjects": int(fit.nobs), **coefficient(fit, "PD")}
    means = [{"group": g, **adjusted_mean(fit, design_group, df, group=g)}
             for g in GROUPS]
    return test, means


def interaction_regression(long: pd.DataFrame, feature: str):
    """feature ~ PD * theta + age + male, SE clustered by subject.

    PD_x_theta = (PD - Control in theta) - (PD - Control in N2).
    """
    df = long.dropna(subset=[feature])
    has_both = df.groupby("subject_id").window.nunique().eq(len(WINDOWS))
    df = df[df.subject_id.isin(has_both[has_both].index)]

    cluster_ids = pd.factorize(df.subject_id)[0]
    fit = sm.OLS(df[feature].astype(float), design_interaction(df)).fit(
        cov_type="cluster", cov_kwds={"groups": cluster_ids})
    test = {"n_subjects": df.subject_id.nunique(), **coefficient(fit, "PD_x_theta")}

    population = df.drop_duplicates("subject_id")  # one row per subject
    means = [{"window": w, "group": g,
              **adjusted_mean(fit, design_interaction, population, window=w, group=g)}
             for w in WINDOWS for g in GROUPS]
    return test, means


def fdr_by_family(tests: pd.DataFrame) -> pd.DataFrame:
    """Benjamini-Hochberg correction within each family."""
    tests = tests.copy()
    tests["p_fdr"], tests["reject_fdr"] = np.nan, False
    for _, block in tests.groupby("family"):
        idx = block.index[block.p_value.notna()]
        if len(idx):
            reject, p_adj, _, _ = multipletests(tests.loc[idx, "p_value"], method="fdr_bh")
            tests.loc[idx, "p_fdr"] = p_adj
            tests.loc[idx, "reject_fdr"] = reject
    return tests


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def pick(df: pd.DataFrame, **conditions) -> pd.Series:
    """Return the single row matching all column == value conditions."""
    mask = np.ones(len(df), dtype=bool)
    for column, value in conditions.items():
        mask &= (df[column] == value).to_numpy()
    return df[mask].iloc[0]


def draw_point(ax, x, row, color, **kwargs):
    """Point estimate with asymmetric 95% CI error bar."""
    ax.errorbar(x, row.estimate,
                yerr=[[row.estimate - row.ci_low], [row.ci_high - row.estimate]],
                fmt="o", capsize=5, color=color, **kwargs)


def plot_within(means: pd.DataFrame, tests: pd.DataFrame, path: Path) -> None:
    """Rows = windows, columns = features; Control vs PD adjusted means."""
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharey="col")
    for i, window in enumerate(WINDOWS):
        for j, feature in enumerate(FEATURES):
            ax = axes[i, j]
            tick_labels = []
            for k, group in enumerate(GROUPS):
                draw_point(ax, k, pick(means, window=window, feature=feature, group=group),
                           COLORS[group])
                zero = pick(tests, analysis="unadjusted_vs_zero", window=window,
                            feature=feature, group=group)
                # show the vs-0 result under the group name so it never hides the data
                tick_labels.append(f"{group}\nvs 0: p_FDR={zero.p_fdr:.3g}")
            diff = pick(tests, analysis="adjusted_group", window=window, feature=feature)
            ax.text(0.5, 0.96, f"PD-Control p_FDR={diff.p_fdr:.3g}", ha="center",
                    va="top", transform=ax.transAxes, fontsize=9)
            ax.axhline(0, color="black", lw=0.8)
            ax.set_xticks([0, 1], tick_labels, fontsize=8)
            ax.set_xlim(-0.5, 1.5)
            ax.margins(y=0.15)  # leave room for the PD-Control label at the top
            ax.set_title(f"{window.upper()} | {feature}")
            ax.grid(axis="y", alpha=0.15)
            if j == 0:
                ax.set_ylabel("Incongruent - Congruent")
    fig.suptitle("Within-window group comparisons (age/sex adjusted)")
    fig.text(0.5, 0.01, "Points: adjusted means with 95% CI. "
             "vs-0 tests: unadjusted, FDR-corrected as a separate family.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


def plot_profiles(means: pd.DataFrame, tests: pd.DataFrame, path: Path) -> None:
    """N2 -> theta profile per group, one panel per feature."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    for ax, feature in zip(axes, FEATURES):
        for group in GROUPS:
            rows = (means[(means.feature == feature) & (means.group == group)]
                    .set_index("window").loc[list(WINDOWS)])
            est = rows.estimate.to_numpy(float)
            ax.errorbar(range(len(WINDOWS)), est,
                        yerr=[est - rows.ci_low.to_numpy(float),
                              rows.ci_high.to_numpy(float) - est],
                        fmt="-o", capsize=4, color=COLORS[group], label=group)
        inter = pick(tests, analysis="group_window_interaction", feature=feature)
        ax.text(0.5, 0.96, f"Group x window p_FDR={inter.p_fdr:.3g}", ha="center",
                va="top", transform=ax.transAxes, fontsize=9)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(range(len(WINDOWS)), [WINDOW_LABELS[w] for w in WINDOWS])
        ax.set_title(feature)
        ax.set_ylabel("Incongruent - Congruent")
        ax.grid(axis="y", alpha=0.15)
    axes[0].legend(frameon=False)
    fig.suptitle("Window profiles (age/sex adjusted)")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--conflict", default="conflict.csv",
                        help="conflict file name inside data/")
    parser.add_argument("--micro", default=None, help="unused; kept for compatibility")
    parser.add_argument("--tag", default="task", help="output subfolder name")
    parser.add_argument("--no-simon-qc", action="store_true",
                        help="keep subjects regardless of trials_simon")
    args = parser.parse_args()

    results_dir = ROOT / "results" / args.tag
    figures_dir = ROOT / "figures" / args.tag
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    data = load_data(args.conflict, simon_qc=not args.no_simon_qc)
    long = to_long(data)
    long.to_csv(results_dir / f"b8_subject_window_features_{args.tag}.csv", index=False)

    rng = np.random.default_rng(SEED)
    summaries, tests, group_means, window_means = [], [], [], []

    # 1) Within each window: effect vs 0 per group, and PD vs Control
    for window in WINDOWS:
        window_df = long[long.window == window]
        for feature in FEATURES:
            for group in GROUPS:
                x = window_df.loc[window_df.group == group, feature].dropna().to_numpy(float)
                result = test_vs_zero(x, rng)
                summaries.append({"window": window, "group": group, "feature": feature,
                                  **{k: v for k, v in result.items() if k != "p_value"}})
                tests.append({"family": "vs_zero_12", "analysis": "unadjusted_vs_zero",
                              "window": window, "group": group, "feature": feature,
                              "n_subjects": result["n"], "estimate": result["mean"],
                              "p_value": result["p_value"]})

            test, means = group_regression(window_df, feature)
            tests.append({"family": "group_6", "analysis": "adjusted_group",
                          "window": window, "feature": feature, **test})
            group_means += [{"window": window, "feature": feature, **m} for m in means]

    # 2) Across windows: does the PD-Control difference change from N2 to theta?
    for feature in FEATURES:
        test, means = interaction_regression(long, feature)
        tests.append({"family": "interaction_3", "analysis": "group_window_interaction",
                      "window": "theta_minus_n2", "feature": feature, **test})
        window_means += [{"feature": feature, **m} for m in means]

    tests = fdr_by_family(pd.DataFrame(tests))
    group_means, window_means = pd.DataFrame(group_means), pd.DataFrame(window_means)

    pd.DataFrame(summaries).to_csv(results_dir / f"b8_group_window_summary_{args.tag}.csv", index=False)
    tests.to_csv(results_dir / f"b8_adjusted_tests_{args.tag}.csv", index=False)
    group_means.to_csv(results_dir / f"b8_adjusted_group_estimates_{args.tag}.csv", index=False)
    window_means.to_csv(results_dir / f"b8_adjusted_window_estimates_{args.tag}.csv", index=False)
    plot_within(group_means, tests, figures_dir / f"b8_within_window_groups_{args.tag}.png")
    plot_profiles(window_means, tests, figures_dir / f"b8_group_window_profiles_{args.tag}.png")

    with pd.option_context("display.float_format", "{:.4g}".format, "display.width", 200):
        print(tests[["family", "window", "group", "feature", "n_subjects", "estimate",
                     "p_value", "p_fdr", "reject_fdr"]].to_string(index=False))


if __name__ == "__main__":
    main()