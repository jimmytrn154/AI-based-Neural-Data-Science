"""B1: reusable statistical functions for the EEG microstate project.

Configuration (2026-10-24): two-sided alpha=.05; Benjamini-Hochberg
FDR within RQ families; paired Wilcoxon primary; mixed model secondary;
partial Spearman adjusting age; 2000 bootstrap resamples for rho CI.

Functions return dictionaries/DataFrames rather than writing clinical columns
into EEG feature tables.  Clinical columns are only merged for analyses.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

ALPHA = 0.05
N_BOOT = 2000
SEED = 42


def fdr_bh(table: pd.DataFrame, p_col: str = "p_value",
           family_col: str = "family", alpha: float = ALPHA) -> pd.DataFrame:
    """Correct nonmissing p-values separately within each prespecified family.
    If `family_col` is absent, treat the table as one family.
    """
    out = table.copy()
    if family_col not in out:
        out[family_col] = "all"
    out["p_fdr"] = np.nan
    out["reject_fdr"] = False
    for _, block in out.groupby(family_col, dropna=False, sort=False):
        valid = block[p_col].notna()
        if not valid.any():
            continue
        indexes = block.index[valid]
        reject, adjusted, _, _ = multipletests(
            block.loc[indexes, p_col].astype(float), alpha=alpha, method="fdr_bh"
        )
        out.loc[indexes, "p_fdr"] = adjusted
        out.loc[indexes, "reject_fdr"] = reject
    return out


def paired_wilcoxon(data: pd.DataFrame, rest_col: str, simon_col: str) -> dict:
    """Paired Simon-minus-Rest contrast, with rank-biserial effect size.

    Zero differences are omitted both from Wilcoxon and effect size;
    ties are assigned average ranks. n_pairs includes zero differences.
    """
    pair = data[[rest_col, simon_col]].apply(pd.to_numeric, errors="coerce").dropna()
    difference = (pair[simon_col] - pair[rest_col]).to_numpy(dtype=float)
    nonzero = difference[difference != 0]
    output = {"n_pairs": len(difference), "n_nonzero": len(nonzero),
              "median_delta": float(np.median(difference)) if len(difference) else np.nan,
              "mean_delta": float(np.mean(difference)) if len(difference) else np.nan,
              "w_stat": np.nan, "p_value": np.nan, "rank_biserial": np.nan,
              "cohens_dz": np.nan}
    if not len(nonzero):
        return output
    ranks = stats.rankdata(np.abs(nonzero), method="average")
    positive = float(ranks[nonzero > 0].sum())
    negative = float(ranks[nonzero < 0].sum())
    output["rank_biserial"] = (positive - negative) / (positive + negative)
    test = stats.wilcoxon(nonzero, alternative="two-sided", zero_method="wilcox")
    output["w_stat"] = float(test.statistic)
    output["p_value"] = float(test.pvalue)
    if len(difference) >= 2 and np.std(difference, ddof=1) > 0:
        output["cohens_dz"] = float(np.mean(difference) / np.std(difference, ddof=1))
    return output


def mixed_model(data: pd.DataFrame, metric_col: str, *,
                task_col: str = "task", group_col: str = "group",
                id_col: str = "subject_id", age_col: str = "age",
                sex_col: str = "gender") -> pd.DataFrame:
    """Random-intercept LMM: metric ~ task * group + age + sex + (1|subject).

    Input must be LONG format (one row per subject and task), not a wide
    microstates_task.csv. Returns coefficient-level tests, not FDR-adjusted.
    """
    columns = [metric_col, task_col, group_col, id_col, age_col, sex_col]
    frame = data[columns].copy()
    frame[metric_col] = pd.to_numeric(frame[metric_col], errors="coerce")
    frame[age_col] = pd.to_numeric(frame[age_col], errors="coerce")
    frame = frame.dropna()
    frame = frame[frame[sex_col].isin(["F", "M"])]
    # RQ1 requires each retained subject to have BOTH tasks.
    counts = frame.groupby(id_col)[task_col].nunique()
    eligible_ids = counts[counts == 2].index
    frame = frame[frame[id_col].isin(eligible_ids)]
    if frame[id_col].nunique() < 4:
        raise ValueError("Insufficient subjects for mixed model")
    if frame.duplicated([id_col, task_col]).any():
        raise ValueError("Multiple observations per subject and task")
    # Explicit reference categories: rest and Control.
    formula = (f'Q("{metric_col}") ~ C(Q("{task_col}"), Treatment(reference="rest"))'
               f' * C(Q("{group_col}"), Treatment(reference="Control"))'
               f' + Q("{age_col}") + C(Q("{sex_col}"))')
    fit = smf.mixedlm(formula, frame, groups=frame[id_col]).fit(reml=False, method="lbfgs")
    terms = [k for k in fit.params.index if k != "Group Var"]
    conf = fit.conf_int()
    return pd.DataFrame({"term": terms, "estimate": fit.params.loc[terms].values,
                         "ci_low": conf.loc[terms, 0].values,
                         "ci_high": conf.loc[terms, 1].values,
                         "p_value": fit.pvalues.loc[terms].values,
                         "n_subjects": frame[id_col].nunique(),
                         "n_observations": len(frame), "converged": fit.converged})


def _partial_spearman_value(x: np.ndarray, y: np.ndarray,
                            age: np.ndarray) -> float:
    """Pearson correlation of x/y rank residuals after regressing out age rank."""
    ranked = np.column_stack([stats.rankdata(v) for v in (x, y, age)])
    design = np.column_stack([np.ones(len(age)), ranked[:, 2]])
    residual_x = ranked[:, 0] - design @ np.linalg.lstsq(design, ranked[:, 0], rcond=None)[0]
    residual_y = ranked[:, 1] - design @ np.linalg.lstsq(design, ranked[:, 1], rcond=None)[0]
    if np.std(residual_x) < 1e-12 or np.std(residual_y) < 1e-12:
        return np.nan
    return float(np.corrcoef(residual_x, residual_y)[0, 1])


def partial_spearman(data: pd.DataFrame, x_col: str, y_col: str,
                     age_col: str = "age", n_boot: int = N_BOOT,
                     seed: int = SEED) -> dict:
    """Age-adjusted partial Spearman rho; approximate t p and percentile bootstrap CI.

    Re-ranks x, y, age in each bootstrap sample. Bootstrap uses subject rows;
    input MUST already have one row per subject.
    """
    frame = data[[x_col, y_col, age_col]].apply(pd.to_numeric, errors="coerce").dropna()
    arrays = frame.to_numpy(dtype=float)
    n = len(arrays)
    if n < 5:
        return {"n": n, "rho": np.nan, "p_value": np.nan,
                "ci_low": np.nan, "ci_high": np.nan, "n_boot_valid": 0}
    x, y, age = arrays.T
    rho = _partial_spearman_value(x, y, age)
    if np.isnan(rho):
        p = np.nan
    elif abs(rho) >= 1:
        p = 0.0
    else:
        t = rho * np.sqrt((n - 3) / (1 - rho ** 2))
        p = float(2 * stats.t.sf(abs(t), df=n - 3))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        sample = arrays[rng.integers(n, size=n)]
        r = _partial_spearman_value(*sample.T)
        if np.isfinite(r):
            boots.append(r)
    ci = np.percentile(boots, [2.5, 97.5]) if boots else [np.nan, np.nan]
    return {"n": n, "rho": rho, "p_value": p,
            "ci_low": float(ci[0]), "ci_high": float(ci[1]),
            "n_boot_valid": len(boots)}


def group_linear_model(data: pd.DataFrame, feature_col: str, *,
                       group_col: str = "group", age_col: str = "age",
                       sex_col: str = "gender") -> dict:
    """PD-Control adjusted difference: feature ~ group + age + sex (HC3 SE)."""
    frame = data[[feature_col, group_col, age_col, sex_col]].copy()
    frame[feature_col] = pd.to_numeric(frame[feature_col], errors="coerce")
    frame[age_col] = pd.to_numeric(frame[age_col], errors="coerce")
    frame = frame.dropna()
    frame = frame[frame[sex_col].isin(["F", "M"]) & frame[group_col].isin(["PD", "Control"])]
    fit = smf.ols(
        f'Q("{feature_col}") ~ C(Q("{group_col}"), Treatment(reference="Control"))'
        f' + Q("{age_col}") + C(Q("{sex_col}"))', data=frame
    ).fit(cov_type="HC3")
    term = next(t for t in fit.params.index if "[T.PD]" in t)
    ci = fit.conf_int().loc[term]
    return {"n": len(frame), "estimate_pd_minus_control": float(fit.params[term]),
            "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
            "p_value": float(fit.pvalues[term])}