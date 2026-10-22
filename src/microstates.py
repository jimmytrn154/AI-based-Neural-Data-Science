"""GFP peaks, clustering, back-fitting, metrics.

Two-level template fitting (config ``microstates``):
  1. per subject, pool an equal number of GFP-peak maps from Rest and Simon
     and cluster them with polarity-invariant modified k-means;
  2. cluster all subjects' individual maps again to get group templates.

Teammates mostly need :func:`load_templates` and :func:`backfit`.
"""

import logging
import os
import zlib

import mne
import numpy as np
import pandas as pd
from pycrostates.cluster import ModKMeans
from pycrostates.io import ChData, read_cluster
from pycrostates.preprocessing import extract_gfp_peaks
from pycrostates.segmentation import RawSegmentation
from pycrostates.utils import _corr_vectors

from src import io as pio
from src import preprocess as pp

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------
# GFP peaks
# ------------------------------------------------------------------

def gfp(data):
    """Global field power: spatial SD across channels, for (n_channels, n_times)."""
    return np.std(data, axis=0)


def extract_peaks(raw, cfg):
    """GFP-peak maps of one recording, without high-GFP outliers.

    Peaks are local GFP maxima at least ``min_peak_distance_samples`` apart.
    Peaks whose GFP exceeds the recording's mean + ``gfp_peak_outlier_sd`` SD
    (computed over the GFP values at the peaks) are dropped. Returns a
    pycrostates ChData of shape (n_channels, n_peaks).
    """
    ms = cfg["microstates"]
    peaks = extract_gfp_peaks(raw, picks="eeg", verbose="error",
                              min_peak_distance=ms["min_peak_distance_samples"])
    data = peaks.get_data()
    peak_gfp = gfp(data)
    limit = peak_gfp.mean() + ms["gfp_peak_outlier_sd"] * peak_gfp.std()
    keep = peak_gfp <= limit
    return ChData(data[:, keep], peaks.info)


def subject_rng(subject_id, cfg):
    """Random generator seeded by ``random_seed`` and the subject ID.

    A subject's peak sample doesn't change when other subjects are added.
    """
    return np.random.default_rng([cfg["random_seed"], zlib.crc32(subject_id.encode())])


def n_peaks_per_task(n_rest, n_simon):
    """Peaks to sample from each task: ``min(n_rest, n_simon, 3000)`` (config)."""
    return min(n_rest, n_simon, 3000)


def sample_peaks(peaks, n, rng):
    """Random sample (without replacement) of ``n`` peak maps from a ChData."""
    data = peaks.get_data()
    idx = np.sort(rng.choice(data.shape[1], size=n, replace=False))
    return ChData(data[:, idx], peaks.info)


def subject_peaks(subject_id, cfg, tasks=("Rest", "Simon")):
    """Load, preprocess and pool an equal number of GFP peaks per task.

    Returns ``(pooled, counts)``: ``pooled`` is a ChData with the sampled
    peaks of all tasks, and ``counts`` a dict of peaks before sampling and
    the number kept per task.
    """
    peaks = {}
    for task in tasks:
        raw = pp.preprocess_raw(pio.load_continuous(subject_id, task, cfg), cfg)
        peaks[task] = extract_peaks(raw, cfg)
    n_avail = {task: p.get_data().shape[1] for task, p in peaks.items()}
    n = n_peaks_per_task(*n_avail.values())
    rng = subject_rng(subject_id, cfg)
    sampled = [sample_peaks(peaks[task], n, rng).get_data() for task in tasks]
    info = peaks[tasks[0]].info
    counts = {f"n_peaks_{task.lower()}": n_avail[task] for task in tasks}
    counts["n_sampled_per_task"] = n
    return ChData(np.hstack(sampled), info), counts


# ------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------

def make_kmeans(k, cfg):
    """Modified k-means with ``n_init``, ``max_iter`` and ``random_seed`` from config."""
    ms = cfg["microstates"]
    return ModKMeans(n_clusters=k, n_init=ms["n_init"], max_iter=ms["max_iter"],
                     random_state=cfg["random_seed"])


def fit_individual(peaks, k, cfg, n_jobs=None):
    """First level: cluster one subject's pooled GFP peaks."""
    cluster = make_kmeans(k, cfg)
    cluster.fit(peaks, n_jobs=n_jobs or os.cpu_count(), verbose="error")
    return cluster


def fit_group(individual, k, cfg, n_jobs=None):
    """Second level: cluster the individual maps of all subjects.

    ``individual`` is a list of fitted ModKMeans (one per subject).
    """
    maps = np.vstack([c.cluster_centers_ for c in individual]).T
    cluster = make_kmeans(k, cfg)
    cluster.fit(ChData(maps, individual[0].info), n_jobs=n_jobs or os.cpu_count(),
                verbose="error")
    return cluster


def template_gev(templates, peaks):
    """GEV of the templates on a set of maps (e.g. one subject's GFP peaks).

    Each map is assigned to the template with the highest absolute spatial
    correlation (polarity ignored). Returns a fraction in [0, 1].
    """
    data = peaks.get_data()
    data = data - data.mean(axis=0)
    g = gfp(data)
    maps = templates.cluster_centers_
    maps = maps - maps.mean(axis=1, keepdims=True)
    maps = maps / np.linalg.norm(maps, axis=1, keepdims=True)
    units = data / np.linalg.norm(data, axis=0)
    best = np.abs(np.einsum("kc,cn->kn", maps, units)).max(axis=0)
    return float(np.sum((g * best) ** 2) / np.sum(g ** 2))


def stratified_sample(table, n, cfg, by="group"):
    """Random sample of ``n`` subject IDs, stratified by ``by`` (seeded).

    Each stratum gets ``n`` times its share of the table, rounded.
    """
    id_col = cfg["subjects"]["id_column"]
    rng = np.random.default_rng(cfg["random_seed"])
    shares = table[by].value_counts(normalize=True).sort_index()
    picked = []
    for level, share in shares.items():
        ids = np.sort(table.loc[table[by] == level, id_col].to_numpy())
        picked.extend(rng.choice(ids, size=round(n * share), replace=False))
    return sorted(picked)


# Draft subjects removed by team decision (PROGRESS.md decision log, Oct 22)
DRAFT_EXCLUDE = {
    "sub-P048": "artifact-dominated individual maps (team decision, Oct 22)",
}


def draft_subjects(cfg, exclude=None):
    """The draft-template subjects: stratified sample, with exclusions replaced.

    Samples ``draft_templates.n_subjects`` from the subjects passing both the
    Rest and Simon rules, then replaces each subject in ``exclude`` (default
    :data:`DRAFT_EXCLUDE`) by a random subject of the same group (seeded).
    Returns ``(ids, replacements)`` with ``replacements`` as {old: new}.
    """
    exclude = DRAFT_EXCLUDE if exclude is None else exclude
    rest, _ = pp.select_subjects("Rest", cfg)
    simon, _ = pp.select_subjects("Simon", cfg)
    both = pp.select_rest_vs_simon(rest, simon, cfg)
    pool = pio.load_subjects(cfg).set_index("subject_id").loc[both, "group"]
    ids = stratified_sample(pool.reset_index(), cfg["microstates"]["draft_templates"]["n_subjects"], cfg)
    rng = np.random.default_rng(cfg["random_seed"])
    replacements = {}
    for sid in exclude:
        if sid not in ids:
            continue
        candidates = sorted(set(pool.index[pool == pool[sid]]) - set(ids))
        replacements[sid] = str(rng.choice(candidates))
        ids = sorted(set(ids) - {sid} | {replacements[sid]})
    return ids, replacements


def fit_templates(subject_ids, k, cfg):
    """Run the two-level fit on the given subjects.

    Returns ``(group, individual, peak_table)``: the group ModKMeans, the
    per-subject ModKMeans (dict by subject ID), and a DataFrame with peak
    counts, individual GEV and template GEV per subject.
    """
    individual, rows = {}, []
    for sid in subject_ids:
        peaks, counts = subject_peaks(sid, cfg)
        individual[sid] = fit_individual(peaks, k, cfg)
        rows.append({"subject_id": sid, **counts, "gev_individual": individual[sid].GEV_,
                     "_peaks": peaks})
        logger.info("%s: %s, individual GEV %.3f", sid, counts, individual[sid].GEV_)
    group = fit_group(list(individual.values()), k, cfg)
    table = pd.DataFrame(rows)
    table["gev_template"] = [template_gev(group, p) for p in table.pop("_peaks")]
    return group, individual, table


def label_templates(cluster, order, names):
    """Reorder and rename templates in place, so state ``i`` is ``names[i]``.

    ``order[i]`` is the current index of the map that becomes ``names[i]``,
    e.g. ``order=[1, 3, 2, 0], names=["A", "B", "C", "D"]``. Labels are
    assigned by visual inspection (config ``state_labels``).
    """
    cluster.reorder_clusters(order=list(order))
    cluster.rename_clusters(new_names=list(names))
    return cluster


# ------------------------------------------------------------------
# Templates on disk
# ------------------------------------------------------------------

def templates_path(k, cfg):
    """``results/templates/templates_K<k>.fif`` (config ``outputs.templates``)."""
    return pio.results_dir(cfg, "templates", f"templates_K{k}.fif")


def save_templates(cluster, cfg):
    """Save group templates to ``templates_path``, overwriting any old file."""
    path = templates_path(cluster.n_clusters, cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    cluster.save(path)
    return path


def load_templates(k, cfg):
    """Load the group templates for ``k`` states (a fitted ModKMeans)."""
    return read_cluster(templates_path(k, cfg))


# ------------------------------------------------------------------
# Back-fitting
# ------------------------------------------------------------------

def backfit_options(cfg, smoothing=True):
    """``ModKMeans.predict`` options from config ``microstates.backfit``.

    ``smoothing=False`` sets ``factor=0`` (the no-smoothing sensitivity analysis).
    """
    bf = cfg["microstates"]["backfit"]
    return {"factor": bf["smoothing_factor"] if smoothing else 0,
            "half_window_size": bf["smoothing_half_window_samples"],
            "min_segment_length": bf["min_segment_length_samples"],
            "reject_edges": bf["reject_edges"]}


def backfit(inst, templates, cfg, smoothing=True, fast=True, **kwargs):
    """Back-fit templates to a Raw or Epochs object.

    Uses :func:`backfit_options` from config; any keyword passed here
    overrides them or adds other ``ModKMeans.predict`` options. The data must
    be preprocessed like the templates (``preprocessing`` in config) and
    contain the same channel names. Units don't matter (spatial correlation
    is scale-free).

    With ``fast=True`` (default), Raw data use :func:`reject_short_segments`
    instead of pycrostates' own short-segment step, which gives identical
    labels (tests/test_backfit.py) but runs in seconds instead of minutes.
    Epochs, and Raw with BAD annotations, always use pycrostates.

    Returns a pycrostates RawSegmentation or EpochsSegmentation; ``.labels``
    holds the state index per sample, with -1 for unlabeled samples.
    """
    options = backfit_options(cfg, smoothing)
    options.update(kwargs)
    min_length = options["min_segment_length"]
    use_library = (not fast or min_length == 0 or not isinstance(inst, mne.io.BaseRaw)
                   or any(d.lower().startswith("bad") for d in inst.annotations.description))
    if use_library:
        return templates.predict(inst, picks="eeg", verbose="error", **options)

    options["min_segment_length"] = 0
    segmentation = templates.predict(inst, picks="eeg", verbose="error", **options)
    raw = segmentation.raw
    labels = reject_short_segments(segmentation.labels,
                                   normalize_samples(raw.get_data()), min_length)
    return RawSegmentation(
        labels=labels,
        inst=raw,
        cluster_centers_=segmentation.cluster_centers_,
        cluster_names=segmentation.cluster_names,
        predict_parameters={**segmentation.predict_parameters,
                            "min_segment_length": min_length},
    )


def normalize_samples(data):
    """Center and scale each time point across channels, as pycrostates does.

    pycrostates normalizes the data in place before rejecting short
    segments, so the same values are used here to get identical results.
    """
    data = data - np.mean(data, axis=0)
    std = np.std(data, axis=0)
    std[std == 0] = 1
    return data / std


def reject_short_segments(labels, data, min_segment_length):
    """Merge segments shorter than ``min_segment_length`` into their neighbours.

    Same rule as pycrostates ``_BaseCluster._reject_short_segments``: take
    the first short segment (not the first or last segment, not unlabeled),
    and give its samples, one at a time from both ends, to the neighbour
    whose boundary sample correlates better. pycrostates rescans the whole
    recording after every merge; here the scan restarts at the segment just
    before the merged one, because nothing earlier can change. ``data`` must
    be normalized with :func:`normalize_samples`.
    """
    labels = np.array(labels, copy=True)
    n = labels.size
    start, previous = 0, None
    while start < n:
        end = start + 1
        while end < n and labels[end] == labels[start]:
            end += 1
        if (start == 0 or end == n or labels[start] == -1
                or end - start >= min_segment_length):
            previous, start = start, end
            continue
        left, right = start, end - 1
        while left <= right:
            left_corr = np.abs(_corr_vectors(data[:, left - 1], data[:, left]))
            right_corr = np.abs(_corr_vectors(data[:, right], data[:, right + 1]))
            if np.abs(right_corr - left_corr) <= 1e-8:
                if left == right:
                    labels[left] = labels[left - 1]
                    left += 1
                else:
                    labels[right] = labels[right + 1]
                    labels[left] = labels[left - 1]
                    right -= 1
                    left += 1
            elif left_corr < right_corr:
                labels[right] = labels[right + 1]
                right -= 1
            elif right_corr < left_corr:
                labels[left] = labels[left - 1]
                left += 1
            else:
                raise RuntimeError(f"Undefined correlation at samples {left}-{right}")
        start, previous = previous, None
    return labels


# ------------------------------------------------------------------
# Task-level metrics
# ------------------------------------------------------------------

def labels_path(subject_id, task, cfg, nosmooth=False):
    """``results/labels/<subject_id>_task-<Task>_labels.npy`` (config ``outputs.labels``).

    The no-smoothing sensitivity labels end in ``_labels_nosmooth.npy``.
    """
    suffix = "_nosmooth" if nosmooth else ""
    return pio.results_dir(cfg, "labels", f"{subject_id}_task-{task}_labels{suffix}.npy")


def segmentation_metrics(segmentation, prefix):
    """Per-state metrics and transitions of one back-fitted recording.

    Columns are ``<prefix>_<metric>_<state>``:
      meandur (ms), occ (per s), cov (% of labeled samples), gev (%),
      and ``<prefix>_trans_<X>_to_<Y>`` for X != Y (row probabilities after
      collapsing repeats, so each state's outgoing transitions sum to 1).
    Bookkeeping: ``<prefix>_gev_total`` (%), ``_n_samples``, ``_unlabeled_frac``.
    """
    params = segmentation.compute_parameters(norm_gfp=True)
    names = segmentation.cluster_names
    row = {}
    for name in names:
        row[f"{prefix}_meandur_{name}"] = params[f"{name}_meandurs"] * 1000
        row[f"{prefix}_occ_{name}"] = params[f"{name}_occurrences"]
        row[f"{prefix}_cov_{name}"] = params[f"{name}_timecov"] * 100
        row[f"{prefix}_gev_{name}"] = params[f"{name}_gev"] * 100
    trans = segmentation.compute_transition_matrix(stat="probability",
                                                   ignore_repetitions=True)
    for i, src in enumerate(names):
        for j, dst in enumerate(names):
            if i != j:
                row[f"{prefix}_trans_{src}_to_{dst}"] = trans[i, j]
    row[f"{prefix}_gev_total"] = sum(row[f"{prefix}_gev_{n}"] for n in names)
    row[f"{prefix}_n_samples"] = segmentation.labels.size
    row[f"{prefix}_unlabeled_frac"] = params["unlabeled"]
    return row


def task_metrics(subject_id, templates, cfg, tasks=("Rest", "Simon"), smoothing=True,
                 save_labels=True):
    """Back-fit the templates to each task's continuous recording.

    Returns one feature row (dict) for the subject. With ``save_labels``,
    the label sequence of each task is saved to :func:`labels_path`
    (suffix ``_nosmooth`` when ``smoothing=False``).
    """
    row = {cfg["subjects"]["id_column"]: subject_id}
    for task in tasks:
        raw = pp.preprocess_raw(pio.load_continuous(subject_id, task, cfg), cfg)
        segmentation = backfit(raw, templates, cfg, smoothing=smoothing)
        if save_labels:
            path = labels_path(subject_id, task, cfg, nosmooth=not smoothing)
            path.parent.mkdir(parents=True, exist_ok=True)
            np.save(path, segmentation.labels)
        row.update(segmentation_metrics(segmentation, task.lower()))
    return row


def segment_durations_ms(labels, sfreq):
    """Durations (ms) of all labeled segments in a label sequence."""
    change = np.flatnonzero(np.diff(labels)) + 1
    bounds = np.concatenate([[0], change, [labels.size]])
    lengths = np.diff(bounds)
    states = labels[bounds[:-1]]
    return lengths[states != -1] * 1000.0 / sfreq
