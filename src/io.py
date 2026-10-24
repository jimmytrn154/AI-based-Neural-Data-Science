"""Loading .set/.fdt, epochs, trial tables, subjects.csv.

All paths come from ``config.yaml``. Import as a package from the repo root
(``from src import io as pio``); don't put ``src/`` itself on ``sys.path``,
or this module will shadow Python's built-in ``io``.
"""

from pathlib import Path

import h5py
import mne
import numpy as np
import pandas as pd
import scipy.io
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

# Task name used in file names -> folder name in Dataset_ALL
TASK_FOLDERS = {
    "Rest": "rest",
    "Simon": "simon",
    "Oddball": "oddball",
    "IntervalTiming": "timing",
}


def load_config(path=None):
    """Read ``config.yaml`` (default: the one at the repo root) into a dict."""
    path = Path(path) if path is not None else REPO_ROOT / "config.yaml"
    with open(path) as f:
        return yaml.safe_load(f)


def data_root(cfg):
    """Return ``paths.data_root`` from the config as a Path."""
    return Path(cfg["paths"]["data_root"])


def results_dir(cfg, *parts):
    """Return a path inside ``paths.results`` (relative to the repo root)."""
    return REPO_ROOT / cfg["paths"]["results"] / Path(*parts)


def load_subjects(cfg):
    """Load ``subjects.csv`` (one row per subject, ``n/a`` read as missing)."""
    return pd.read_csv(data_root(cfg) / "subjects.csv", na_values="n/a")


def continuous_path(subject_id, task, cfg):
    """Path of the cleaned continuous ``.set`` file for one subject and task."""
    folder = TASK_FOLDERS[task]
    name = f"{subject_id}_task-{task}_desc-clean_eeg.set"
    return data_root(cfg) / "continuous_data" / folder / name


def load_continuous(subject_id, task, cfg, preload=True):
    """Load one continuous recording as an MNE Raw object.

    Note that MNE stores the data in volts; use ``raw.get_data(units="uV")``
    (or multiply by 1e6) to get microvolts.
    """
    path = continuous_path(subject_id, task, cfg)
    return mne.io.read_raw_eeglab(path, preload=preload, verbose="error")


def read_qc(subject_id, task, cfg):
    """Return the scalar fields of ``EEG.etc.iowa`` (nInterp, nICremoved, ...).

    Reads only the ``.set`` header, not the EEG data in the ``.fdt`` file.
    """
    path = continuous_path(subject_id, task, cfg)
    mat = scipy.io.loadmat(
        path, squeeze_me=True, struct_as_record=False, variable_names=["etc"]
    )
    iowa = mat["etc"].iowa
    qc = {}
    for field in iowa._fieldnames:
        value = getattr(iowa, field)
        if np.isscalar(value) or isinstance(value, str):
            qc[field] = value
    return qc


def trials_path(subject_id, task, cfg):
    """Path of the per-trial table (``.tsv``) for one subject and task."""
    folder = TASK_FOLDERS[task]
    name = f"{subject_id}_task-{task}_trials.tsv"
    return data_root(cfg) / "trial_info" / folder / name


def load_trials(subject_id, task, cfg):
    """Load the per-trial condition labels and behavior as a DataFrame."""
    return pd.read_csv(trials_path(subject_id, task, cfg), sep="\t", na_values="n/a")


def epochs_path(subject_id, task, lock, cfg):
    """Path of a ready-made epoch file, e.g. ``lock="cue"`` or ``"resp"``."""
    folder = TASK_FOLDERS[task]
    name = f"{subject_id}_task-{task}_lock-{lock}_epo.mat"
    return data_root(cfg) / "epoched_data" / folder / name


def _read_cellstr(f, dataset):
    """Decode a MATLAB cell array of strings stored in an HDF5 file."""
    refs = np.asarray(dataset).ravel()
    return ["".join(map(chr, np.asarray(f[r]).ravel())) for r in refs]


def load_epochs(subject_id, task, lock, cfg):
    """Load a ready-made epoch file (MATLAB v7.3 / HDF5).

    h5py returns the axes reversed ([trial, time, channel]); they are
    transposed here to MNE order. Bad epochs are *not* removed.

    Returns a dict with
        data       : array (n_trials, n_channels, n_times), microvolts
        times      : array (n_times,), milliseconds
        srate      : float, Hz
        ch_names   : list of channel names, in data order
        bad_epoch  : bool array (n_trials,)
        erp_baseline : [start, end] in ms, already subtracted

    The trial table embedded in the file is a MATLAB ``table`` object that
    h5py cannot decode; use :func:`match_trials` to attach the trial rows.
    """
    with h5py.File(epochs_path(subject_id, task, lock, cfg), "r") as f:
        return {
            "data": np.asarray(f["data"]).transpose(0, 2, 1),
            "times": np.asarray(f["times"]).ravel(),
            "srate": float(np.asarray(f["srate"]).squeeze()),
            "ch_names": _read_cellstr(f, f["labels"]),
            "bad_epoch": np.asarray(f["badEpoch"]).ravel().astype(bool),
            "erp_baseline": np.asarray(f["erpBaseline"]).ravel(),
        }


def match_trials(epochs, trials, lock):
    """Return the trial-table rows matching each epoch, in epoch order.

    Verified on all 147 Simon subjects:
      - ``lock="cue"``: one epoch per trial-table row;
      - ``lock="resp"``: one epoch per row with a non-missing
        ``onset_response``. Rows with ``beh_accuracy == 99`` (no response)
        still have an ``onset_response`` at the end of the response window
        (``rt_eeg`` ~ 1.0 s), so their epochs are locked to a timeout, not a
        key press. Drop them before response-locked analyses.

    The order is the trial table's order: verified Oct 24 for all 147 Simon
    subjects by rebuilding every epoch from the continuous data and the
    onsets (notebooks/chuong_epoch_order_check.py; max difference 2.4e-4 µV,
    0 mismatching trials). A ``bad_epoch`` column is added. Raises
    ValueError if the counts differ.
    """
    if lock == "resp":
        trials = trials[trials["onset_response"].notna()]
    trials = trials.reset_index(drop=True)
    n_epochs = len(epochs["bad_epoch"])
    if len(trials) != n_epochs:
        raise ValueError(
            f"{n_epochs} epochs but {len(trials)} trial rows for lock={lock!r}"
        )
    return trials.assign(bad_epoch=epochs["bad_epoch"])


def epochs_to_mne(epochs, montage="standard_1005"):
    """Convert a :func:`load_epochs` dict to an ``mne.EpochsArray`` (volts).

    Bad epochs are kept, so the dict's ``bad_epoch`` still lines up. Channel
    positions come from the standard template montage, as in the dataset.
    """
    info = mne.create_info(epochs["ch_names"], epochs["srate"], ch_types="eeg")
    data = epochs["data"] * 1e-6  # microvolts -> volts
    tmin = epochs["times"][0] / 1000.0
    out = mne.EpochsArray(data, info, tmin=tmin, baseline=None, verbose="error")
    return out.set_montage(montage)
