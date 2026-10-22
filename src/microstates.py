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

import numpy as np
import pandas as pd
from pycrostates.cluster import ModKMeans
from pycrostates.io import ChData, read_cluster
from pycrostates.preprocessing import extract_gfp_peaks

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

def backfit(inst, templates, cfg, **kwargs):
    """Back-fit templates to a Raw or Epochs object.

    Uses ``min_segment_length_samples`` and ``reject_edges`` from config;
    any keyword passed here overrides them or adds other
    ``ModKMeans.predict`` options. The data must be preprocessed like the
    templates (``preprocessing`` in config) and contain the same channel names. Units
    don't matter (spatial correlation is scale-free).

    Returns a pycrostates RawSegmentation or EpochsSegmentation; ``.labels``
    holds the state index per sample, with -1 for unlabeled samples.
    """
    bf = cfg["microstates"]["backfit"]
    options = {"min_segment_length": bf["min_segment_length_samples"],
               "reject_edges": bf["reject_edges"]}
    options.update(kwargs)
    return templates.predict(inst, picks="eeg", verbose="error", **options)
