"""The fast back-fit must give exactly the same labels as pycrostates.

Run from the repo root: ``.venv/bin/python -m pytest tests/``

The real-data tests need ``paths.data_root`` and the K = 4 templates; they
are skipped otherwise. They are slow (pycrostates' own short-segment step).
"""

import numpy as np
import pytest
from pycrostates.cluster._base import _BaseCluster

from src import io as pio
from src import microstates as ms
from src import preprocess as pp

CFG = pio.load_config()
HAVE_DATA = pio.data_root(CFG).exists() and ms.templates_path(4, CFG).exists()


@pytest.mark.parametrize("min_length", [2, 3, 5, 8])
def test_reject_short_segments_synthetic(min_length):
    """Random labels and data: same result as pycrostates' own function."""
    rng = np.random.default_rng(min_length)
    for _ in range(100):
        n = int(rng.integers(10, 400))
        labels = rng.integers(0, 4, n)
        labels[: int(rng.integers(0, 3))] = -1
        labels[n - int(rng.integers(0, 3)):] = -1
        data = ms.normalize_samples(rng.normal(size=(8, n)))
        expected = _BaseCluster._reject_short_segments(labels.copy(), data.copy(), min_length)
        assert np.array_equal(ms.reject_short_segments(labels, data, min_length), expected)


@pytest.fixture(scope="module")
def templates():
    return ms.load_templates(4, CFG)


@pytest.mark.skipif(not HAVE_DATA, reason="dataset or templates not available")
@pytest.mark.parametrize("smoothing", [True, False], ids=["smooth", "nosmooth"])
@pytest.mark.parametrize("task", ["Rest", "Simon"])
@pytest.mark.parametrize("subject_id", ms.draft_subjects(CFG)[0] if HAVE_DATA else [])
def test_backfit_matches_pycrostates(subject_id, task, smoothing, templates):
    """Draft subjects: 0 sample mismatches between fast and library back-fit."""
    raw = pp.preprocess_raw(pio.load_continuous(subject_id, task, CFG), CFG)
    fast = ms.backfit(raw, templates, CFG, smoothing=smoothing).labels
    library = ms.backfit(raw, templates, CFG, smoothing=smoothing, fast=False).labels
    assert fast.shape == library.shape
    assert int(np.sum(fast != library)) == 0
