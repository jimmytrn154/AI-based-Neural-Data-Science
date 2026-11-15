"""B2: Within-subject Rest vs Simon microstate comparisons.

Run: python src/B2_rest_vs_simon.py --micro microstate_A.csv [--tag A]
     (normally called by src/run_pipeline.py; see src/pipeline_io.py)

Primary: paired Wilcoxon (Simon - Rest) per metric and state, FDR within RQ1.
Secondary: random-intercept mixed model
    metric ~ task * group + age + sex + (1 | subject)   (stats.mixed_model)
Also: PD vs Control on the within-subject change (OLS, HC3).
States (K) are detected from the CSV columns.
Outputs (results/<tag>/): rq1_stats_<tag>.csv, rq1_individual_changes_<tag>.csv,
rq1_mixed_<tag>.csv, rq1_group_delta_<tag>.csv;
(figures/<tag>/) bohyeon_rest_vs_task_all_metrics_<tag>.png
"""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

# 1. File paths
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"

sys.path.insert(0, str(HERE))
from stats import paired_wilcoxon, fdr_bh, mixed_model
from pipeline_io import (SUBJECTS_FILE, parse_args, load_microstates,
                         detect_states, apply_simon_qc, save_run_config)

ctx = parse_args(__doc__)
RESULTS_DIR = ctx.results_dir
FIGURES_DIR = ctx.figures_dir

# 2. Load data
micro = load_microstates(ctx)
subjects = pd.read_csv(DATA_DIR / SUBJECTS_FILE)

# 3. Analysis settings
STATES = detect_states(micro)

ALL_METRICS = {
    "cov": "Coverage",
    "meandur": "Mean duration",
    "occ": "Occurrence",
    "gev": "Global explained variance"}

UNITS = {
    "cov": "percentage points",
    "meandur": "ms",
    "occ": "events / second",
    "gev": "percentage points"}

# Keep only metrics present for every state in both tasks
METRICS = {
    metric: label for metric, label in ALL_METRICS.items()
    if all(f"{task}_{metric}_{s}" in micro.columns
           for task in ["rest", "simon"] for s in STATES)}
skipped = sorted(set(ALL_METRICS) - set(METRICS))
if skipped:
    print("NOTE: metrics not in this CSV, skipped:", skipped)

N_BOOT = 2000
SEED = 42
rng = np.random.default_rng(SEED)

qc = subjects[["subject_id", "trials_simon", "group", "age", "gender"]]

data = micro.merge(
    qc,
    on="subject_id",
    how="inner",
    validate="one_to_one")

# Exclude Simon records marked 'exclude' (unless disabled for this run)
data = apply_simon_qc(data, ctx)
print("QC-eligible subjects:", len(data))


# 4. Paired Wilcoxon tests and descriptive statistics
rows = []
individual = data[["subject_id"]].copy()

for metric in METRICS:
    for state in STATES:

        rest_col = f"rest_{metric}_{state}"
        simon_col = f"simon_{metric}_{state}"

        # Paired measurements
        pair = data[[rest_col, simon_col]].apply(
            pd.to_numeric,
            errors="coerce")

        # Within-subject difference
        delta = pair[simon_col] - pair[rest_col]

        # Save individual differences
        individual[f"delta_{metric}_{state}"] = delta

        # Remove subjects with missing values
        valid = pair.dropna()
        differences = (valid[simon_col] - valid[rest_col]).to_numpy(dtype=float)

        # Statistical test
        outcome = paired_wilcoxon(data, rest_col, simon_col)

        # Descriptive statistics
        mean_rest = valid[rest_col].mean()
        mean_simon = valid[simon_col].mean()
        sd_delta = differences.std(ddof=1)

        # Bootstrap 95% CI for mean difference
        boot_means = []

        for _ in range(N_BOOT):
            sample = rng.choice(
                differences,
                size=len(differences),
                replace=True)

            boot_means.append(sample.mean())

        ci_low, ci_high = np.percentile(boot_means, [2.5, 97.5])

        rows.append({
            "family": "RQ1_paired",
            "metric": metric,
            "state": state,
            "comparison": "simon_minus_rest",
            **outcome,
            "mean_rest": mean_rest,
            "mean_simon": mean_simon,
            "sd_delta": sd_delta,
            "ci_low": ci_low,
            "ci_high": ci_high})


# 5. FDR correction
paired = pd.DataFrame(rows)
paired = fdr_bh(paired, p_col="p_value", family_col="family")

# Save statistical results
paired.to_csv(ctx.csv("rq1_stats"), index=False)

# Save individual changes
individual.to_csv(ctx.csv("rq1_individual_changes"), index=False)

print("\nPaired Wilcoxon tests completed.")
print(paired[[
            "metric",
            "state",
            "n_pairs",
            "mean_delta",
            "median_delta",
            "p_value",
            "p_fdr"
        ]].to_string(index=False))


# 6. Four-panel figure: all metrics
fig, axes = plt.subplots(2, 2,figsize=(12, 8))

for ax, metric in zip(axes.flat, METRICS):

    subset = paired[
        paired["metric"] == metric
    ].set_index("state").loc[list(STATES)]

    mean = subset["mean_delta"].to_numpy(dtype=float)
    lower = subset["ci_low"].to_numpy(dtype=float)
    upper = subset["ci_high"].to_numpy(dtype=float)

    x = np.arange(len(STATES))

    # Bar plot
    ax.bar(x, mean, color="steelblue", alpha=0.8)

    # Bootstrap 95% confidence intervals
    ax.errorbar(
        x,
        mean,
        yerr=[
            np.maximum(mean - lower, 0),
            np.maximum(upper - mean, 0)
        ],
        fmt="none",
        ecolor="black",
        capsize=4,
        linewidth=1.2)

    ax.axhline( 0, color="black", linewidth=0.8)

    ax.set_xticks(x)
    ax.set_xticklabels(list(STATES))

    ax.set_title(METRICS[metric])
    ax.set_xlabel("Microstate")
    ax.set_ylabel(f"Simon - Rest ({UNITS[metric]})")

    # Mark FDR-significant comparisons
    values = np.concatenate([lower, upper])
    span = max(np.ptp(values), 1e-6)
    offset = 0.08 * span

    for i, state in enumerate(STATES):
        if bool(subset.loc[state, "reject_fdr"]):

            if mean[i] >= 0:
                y = upper[i] + offset
                va = "bottom"
            else:
                y = lower[i] - offset
                va = "top"

            ax.text(
                i,
                y,
                "*",
                ha="center",
                va=va,
                fontsize=16,
                fontweight="bold")

    ax.margins(y=0.20)

# Hide panels not used (if a metric is missing)
for ax in list(axes.flat)[len(METRICS):]:
    ax.axis("off")

fig.suptitle(
    f"Rest vs. Simon: Within-subject Microstate Changes "
    f"[{ctx.name}, K={len(STATES)}]",
    fontsize=15)

fig.text(
    0.5,
    0.01,
    "Bars: mean paired difference; error bars: bootstrap 95% CI; "
    "*: FDR-adjusted p < 0.05",
    ha="center",
    fontsize=9)

fig.tight_layout(rect=[0, 0.04, 1, 0.96])

fig.savefig(
    ctx.fig("bohyeon_rest_vs_task_all_metrics"),
    dpi=300,
    bbox_inches="tight")
plt.close(fig)


# 7. PD vs Control: regression on within-subject differences
regression_rows = []

for metric in METRICS:
    for state in STATES:

        rest_col = f"rest_{metric}_{state}"
        simon_col = f"simon_{metric}_{state}"

        # Select relevant columns
        frame = data[[
                "subject_id",
                "group",
                "age",
                "gender",
                rest_col,
                simon_col]].dropna().copy()

        # Conflicting sex values (e.g. 'F|M') are treated as missing (config)
        frame = frame[
            frame["group"].isin(["PD", "Control"])
            & frame["gender"].isin(["F", "M"])
        ].copy()

        # Calculate individual change
        frame["delta"] = (frame[simon_col] - frame[rest_col])

        # Age- and gender-adjusted regression
        model = smf.ols(
            'delta ~ C(group, Treatment(reference="Control"))'
            ' + age + C(gender)',
            data=frame
        ).fit(cov_type="HC3")

        # PD - Control coefficient
        term = 'C(group, Treatment(reference="Control"))[T.PD]'

        ci = model.conf_int().loc[term]

        regression_rows.append({
            "family": "RQ1_group_delta",
            "metric": metric,
            "state": state,
            "n_subjects": int(model.nobs),
            "estimate_group_difference": model.params[term],
            "ci_low": ci.iloc[0],
            "ci_high": ci.iloc[1],
            "p_value": model.pvalues[term]})


# 8. FDR correction for group comparisons
regression = pd.DataFrame(regression_rows)
regression = fdr_bh(regression, p_col="p_value", family_col="family")
regression.to_csv(ctx.csv("rq1_group_delta"), index=False)

print("\nPD vs Control differences in within-subject changes:")
print(regression[[
            "metric",
            "state",
            "estimate_group_difference",
            "p_value",
            "p_fdr"]].to_string(index=False))

# 9. Secondary analysis: mixed models (long format, one row per subject x task)
TASK_TERM = 'C(Q("task"), Treatment(reference="rest"))[T.simon]'
GROUP_TERM = 'C(Q("group"), Treatment(reference="Control"))[T.PD]'
TERM_LABELS = {
    TASK_TERM: "task_simon (in Control)",
    GROUP_TERM: "group_PD (at rest)",
    f"{TASK_TERM}:{GROUP_TERM}": "task_x_group",}

mixed_rows = []

for metric in METRICS:
    for state in STATES:
        long = pd.concat([
            data[["subject_id", "group", "age", "gender"]].assign(
                task=task,
                value=data[f"{task}_{metric}_{state}"])
            for task in ["rest", "simon"]], ignore_index=True)

        table = mixed_model(long, "value")
        table = table[table["term"].isin(TERM_LABELS)].copy()
        table["term"] = table["term"].map(TERM_LABELS)
        table["family"] = "RQ1_mixed_" + table["term"]
        table["metric"] = metric
        table["state"] = state
        mixed_rows.append(table)

mixed = fdr_bh(
    pd.concat(mixed_rows, ignore_index=True),
    p_col="p_value",
    family_col="family")

mixed = mixed[
    ["family", "metric", "state", "term", "estimate", "ci_low", "ci_high",
     "p_value", "p_fdr", "reject_fdr", "n_subjects", "n_observations",
     "converged", "optimizer"]]

mixed.to_csv(ctx.csv("rq1_mixed"), index=False)

print("\nMixed models (secondary): task and task x group terms")

print(mixed[mixed["term"] != "group_PD (at rest)"][
        ["metric", "state", "term", "estimate", "p_value", "p_fdr"]
    ].to_string(index=False))

print("Non-converged fits:",
    int((~mixed["converged"].astype(bool)).sum()))

save_run_config(ctx, STATES, "B2_rest_vs_simon")

# 10. Final output summary
print(f"\nB2 completed [{ctx.name}].")
print("Saved files in", RESULTS_DIR, "and", FIGURES_DIR)
for stem in ["rq1_stats", "rq1_individual_changes", "rq1_group_delta", "rq1_mixed"]:
    print(" -", ctx.csv(stem).name)
print(" -", ctx.fig("bohyeon_rest_vs_task_all_metrics").name)