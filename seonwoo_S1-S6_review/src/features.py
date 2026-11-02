"""Summary features (reorganization, conflict effect, recovery) and behavior."""

import numpy as np

CONDITIONS = {"con": 1, "inc": 2}  # beh_congruentcondition


def simon_behavior(trials):
    """Simon behavioral features of one subject (docs/TEAM_TASKS.md, behavior.csv).

    ``trials`` is the subject's Simon trial table (:func:`src.io.load_trials`).
    RTs (ms) are means over correct trials; error rates are errors over
    responded trials (``beh_accuracy`` 1 or 2); the no-response rate is over
    all trials; ``simon_n_<cond>`` counts all trials of that condition.
    Effects are incongruent minus congruent. Rates are fractions (0-1).
    """
    row = {}
    for name, code in CONDITIONS.items():
        cond = trials[trials["beh_congruentcondition"] == code]
        correct = cond["beh_accuracy"] == 1
        responded = cond["beh_accuracy"].isin([1, 2])
        row[f"simon_rt_{name}"] = cond.loc[correct, "beh_rt"].mean()
        row[f"simon_err_{name}"] = ((cond["beh_accuracy"] == 2).sum() / responded.sum()
                                    if responded.any() else np.nan)
        row[f"simon_n_{name}"] = len(cond)
    row["simon_effect_rt"] = row["simon_rt_inc"] - row["simon_rt_con"]
    row["simon_effect_err"] = row["simon_err_inc"] - row["simon_err_con"]
    row["simon_noresp_rate"] = (trials["beh_accuracy"] == 99).mean()
    order = ["simon_rt_con", "simon_rt_inc", "simon_effect_rt", "simon_err_con",
             "simon_err_inc", "simon_effect_err", "simon_noresp_rate",
             "simon_n_con", "simon_n_inc"]
    return {key: row[key] for key in order}


# ------------------------------------------------------------------
# Recovery pilot (docs/TEAM_TASKS.md, C4)
# ------------------------------------------------------------------

def state_probabilities(labels, n_states):
    """P(state) per time point from trial labels (n_trials, n_times).

    Probabilities are over the trials labeled at each time point (-1 =
    unlabeled is left out), so they sum to 1 wherever any trial is labeled;
    time points with no labeled trial are NaN. Returns (n_states, n_times).
    """
    labeled = (labels >= 0).sum(axis=0).astype(float)
    labeled[labeled == 0] = np.nan
    return np.stack([(labels == s).sum(axis=0) / labeled for s in range(n_states)])


def recovery_curve(probs, baseline):
    """Euclidean distance from the baseline probability vector at each time point."""
    return np.sqrt(((probs - baseline[:, None]) ** 2).sum(axis=0))


def window_mean(values, times, window):
    """Mean of ``values`` over ``window`` = (start, end), in the units of ``times``."""
    mask = (times >= window[0]) & (times <= window[1])
    return float(np.nanmean(values[mask]))


def spearman_brown(r):
    """Spearman-Brown correction of a split-half correlation: 2r / (1 + r)."""
    return 2 * r / (1 + r)
