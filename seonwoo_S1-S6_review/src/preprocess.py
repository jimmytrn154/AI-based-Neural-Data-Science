"""Filtering, downsampling, exclusion.

Every rule here is read from ``config.yaml``; nothing is hard-coded.
"""

import logging

import pandas as pd

from src import io as pio

logger = logging.getLogger(__name__)


def qc_table(subject_ids, task, cfg):
    """Collect ``nInterp`` and ``nICremoved`` for each subject's recording.

    Returns a DataFrame indexed by subject_id. Subjects whose ``.set`` file
    is missing get NaN.
    """
    rows = []
    for sid in subject_ids:
        if pio.continuous_path(sid, task, cfg).exists():
            qc = pio.read_qc(sid, task, cfg)
            rows.append({"subject_id": sid, "has_file": True,
                         "nInterp": qc["nInterp"], "nICremoved": qc["nICremoved"]})
        else:
            rows.append({"subject_id": sid, "has_file": False})
    return pd.DataFrame(rows).set_index("subject_id")


def select_subjects(task, cfg, subjects=None):
    """Apply the config exclusion rules for ``"Rest"`` or ``"Simon"``.

    Steps, in order:
      1. all subjects in subjects.csv
      2. ``subjects.<task>.require`` (e.g. ``has_rest == 1``)
      3. Simon only: ``trials_simon`` in ``allowed_trial_status``
      4. continuous file present
      5. ``nInterp <= max_interpolated_channels``

    Returns ``(table, counts)``: ``table`` holds the kept subjects with their
    QC values (Rest also gets a ``blink_flag`` column for the blink
    sensitivity check), and ``counts`` lists the n remaining after each step.
    """
    rules = cfg["subjects"]
    task_rules = rules[task.lower()]
    df = pio.load_subjects(cfg) if subjects is None else subjects.copy()
    id_col = rules["id_column"]
    counts = [("all subjects in subjects.csv", len(df))]

    df = df.query(task_rules["require"])
    counts.append((task_rules["require"], len(df)))

    if "allowed_trial_status" in task_rules:
        status_col = f"trials_{task.lower()}"
        allowed = task_rules["allowed_trial_status"]
        df = df[df[status_col].isin(allowed)]
        counts.append((f"{status_col} in {allowed}", len(df)))

    qc = qc_table(df[id_col], task, cfg)
    df = df.join(qc, on=id_col)
    df = df[df["has_file"]]
    counts.append(("continuous file present", len(df)))

    max_interp = rules["max_interpolated_channels"]
    df = df[df["nInterp"] <= max_interp]
    counts.append((f"nInterp <= {max_interp}", len(df)))

    if "blink_sensitivity_flag" in task_rules:
        df = df.assign(blink_flag=df.eval(task_rules["blink_sensitivity_flag"]))

    counts = pd.DataFrame(counts, columns=["step", "n_remaining"])
    for step, n in counts.itertuples(index=False):
        logger.info("%s | %s: %d", task, step, n)
    return df.reset_index(drop=True), counts


def select_rest_vs_simon(rest_table, simon_table, cfg):
    """Subject IDs passing BOTH the Rest and Simon rules (sorted)."""
    id_col = cfg["subjects"]["id_column"]
    both = set(rest_table[id_col]) & set(simon_table[id_col])
    logger.info("Rest and Simon | passing both: %d", len(both))
    return sorted(both)


def preprocess_raw(raw, cfg):
    """Band-pass and resample a Raw object for microstates (in place).

    Uses ``preprocessing.bandpass_hz`` with MNE's default zero-phase FIR,
    then resamples to ``preprocessing.resample_hz``. The data are already
    average-referenced, so no re-referencing is done. Returns ``raw``.
    """
    pp = cfg["preprocessing"]
    l_freq, h_freq = pp["bandpass_hz"]
    raw.filter(l_freq, h_freq, verbose="error")
    raw.resample(pp["resample_hz"], verbose="error")
    return raw
