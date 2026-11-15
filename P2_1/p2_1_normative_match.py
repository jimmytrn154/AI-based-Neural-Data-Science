# %% [markdown]
# # P2.1 E vs. normative templates (Minseo; repo version by Chuong)
# Compares our K = 5 maps (and the Rest-only / Simon-only K = 5 maps) with
# 12 published template sets from MS-Template-Explorer (Koenig, GPL-3.0).
# Same method as Minseo's `code/colab_pipeline.py`:
#   - primary: centred Pearson |r| on the channels both montages share
#     (matched by name); centring over the shared channels equals
#     re-referencing both maps to the average of those channels;
#   - diagnostic only: spline interpolation of the published maps onto our
#     positions using nominal standard_1005 coordinates.
# Adds: E best / second-best match and gap, A-D controls, the within-set
# baseline between non-matching reference classes, missing scalp regions.
#
# Run from the repo root:
#   python P2_1/p2_1_normative_match.py --refs <folder with the 12 .set files> [--download]
# Outputs: results/features/p2_1_normative_*.csv,
#          results/figures/minseo_p2_1_normative_match.png

# %%
import argparse
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mne.channels import make_standard_montage
from mne.channels.interpolation import _make_interpolation_matrix
from scipy.io import loadmat

from src import io as pio, microstates as ms

REF_URL = "https://raw.githubusercontent.com/ThomasKoenigBern/MS-Template-Explorer/main/Template%20Maps/"
REFERENCE_FILES = [
    "Custo 2017.set", "Zanesco 2020.set", "Zanesco 2021a.set", "Hwang 2025.set",
    "Qin2022 Depression_Controls.set", "Qin2022 Depression_Patients.set",
    "Zelenina_2022 Oxytocin.set", "Giordano 2018.set", "Lassi2023.set",
    "Milz_2016_Cognitive Style.set", "Schiller2021 Alcohol Intoxication.set",
    "Toplutas 2023 DOC_HC.set",
]
MIN_COMMON = 15


# %% Reading the published templates
def as_list(obj):
    if obj is None:
        return []
    if isinstance(obj, np.ndarray):
        return obj.ravel().tolist()
    return list(obj) if isinstance(obj, (list, tuple)) else [obj]


def field(obj, key):
    return obj.get(key) if isinstance(obj, dict) else getattr(obj, key, None)


def letter(label):
    """Short class letter from a published label ('E', 'Class E', 'MS_E' -> 'E')."""
    label = str(label).strip()
    if re.fullmatch(r"[A-Z]", label):
        return label
    m = re.search(r"(?:CLASS|STATE|MICROSTATE|[_ -])([A-G])$", label, re.I)
    return m.group(1).upper() if m else label


def read_reference(path):
    """Published maps, labels and channel names from an EEGLAB .set (MSMaps)."""
    try:
        eeg = loadmat(path, simplify_cells=True)["EEG"]
    except NotImplementedError:  # MATLAB v7.3 (HDF5) files, e.g. Custo 2017
        import mat73
        eeg = mat73.loadmat(str(path))["EEG"]
    channels = [str(field(c, "labels")).strip() for c in as_list(field(eeg, "chanlocs"))]
    sets = []
    for sol in as_list(field(field(eeg, "msinfo"), "MSMaps")):
        maps, labels = field(sol, "Maps"), field(sol, "Labels")
        if maps is None or labels is None:
            continue
        maps = np.asarray(maps, dtype=float)
        if maps.ndim != 2 or maps.size == 0:
            continue
        if maps.shape[1] != len(channels):
            maps = maps.T
        labels = [str(x).strip() for x in as_list(labels)]
        if maps.shape[1] == len(channels) and maps.shape[0] == len(labels) and 3 <= len(labels) <= 12:
            sets.append({"maps": maps, "labels": labels, "channels": channels})
    if not sets:
        raise ValueError(f"no labelled MSMaps in {Path(path).name}")
    return sets


# %% Correlation
def abs_r(a, b):
    """Centred Pearson |r| between two maps (polarity ignored)."""
    a, b = a - a.mean(), b - b.mean()
    return abs(float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))))


def matrix(rows, cols):
    return np.array([[abs_r(x, y) for y in cols] for x in rows])


def common_channels(ours, theirs):
    lookup = {c.upper(): i for i, c in enumerate(theirs)}
    return [c for c in ours if c.upper() in lookup], lookup


def compare(our_maps, our_ch, ref):
    common, lookup = common_channels(our_ch, ref["channels"])
    if len(common) < MIN_COMMON:
        raise ValueError(f"only {len(common)} common channels")
    a = our_maps[:, [our_ch.index(c) for c in common]]
    b = ref["maps"][:, [lookup[c.upper()] for c in common]]
    return matrix(a, b), matrix(b, b), common


STANDARD = {k.upper(): np.asarray(v) for k, v in
            make_standard_montage("standard_1005").get_positions()["ch_pos"].items()}


def interpolation_diagnostic(our_maps, our_ch, ref):
    """DIAGNOSTIC ONLY: published maps splined onto our channels (nominal standard_1005)."""
    src = [i for i, c in enumerate(ref["channels"]) if c.upper() in STANDARD]
    dst = [i for i, c in enumerate(our_ch) if c.upper() in STANDARD]
    if len(src) < 20:  # as in Minseo's pipeline: too few positions to interpolate
        return np.full((len(our_maps), len(ref["labels"])), np.nan)
    with np.errstate(all="ignore"):  # spurious floating-point warnings from the BLAS matmul
        W = _make_interpolation_matrix(np.array([STANDARD[ref["channels"][i].upper()] for i in src]),
                                       np.array([STANDARD[our_ch[i].upper()] for i in dst]), alpha=1e-5)
        return matrix(our_maps[:, dst], ref["maps"][:, src] @ W.T)


# %% Scalp regions of channels missing from a reference montage
def region(ch):
    name = ch.upper()
    if name.startswith(("FP", "AF")):
        area = "frontopolar"
    elif name.startswith(("FT", "T", "TP")):
        area = "temporal"
    elif name.startswith(("FC", "F")):
        area = "frontal"
    elif name.startswith(("CP", "C")):
        area = "central"
    elif name.startswith(("PO", "P")):
        area = "parietal"
    elif name.startswith(("O", "I")):
        area = "occipital/inion"
    else:
        area = "other"
    digits = re.findall(r"\d+", name)
    side = "midline" if name.endswith("Z") else ("left" if digits and int(digits[-1]) % 2 else "right")
    return f"{area} ({side})" if digits or name.endswith("Z") else area


def describe_missing(our_ch, common):
    missing = [c for c in our_ch if c not in common]
    counts = pd.Series([region(c) for c in missing]).value_counts()
    return missing, "; ".join(f"{k} {v}" for k, v in counts.items())


# %% Main
def main():
    parser = argparse.ArgumentParser(description="P2.1 normative template comparison")
    parser.add_argument("--refs", type=Path, required=True, help="folder with the 12 reference .set files")
    parser.add_argument("--download", action="store_true", help="download missing reference files first")
    args = parser.parse_args()
    args.refs.mkdir(parents=True, exist_ok=True)
    if args.download:
        for name in REFERENCE_FILES:
            if not (args.refs / name).exists():
                urllib.request.urlretrieve(REF_URL + urllib.parse.quote(name), args.refs / name)

    cfg = pio.load_config()
    k = cfg["microstates"]["primary_k"]
    ours = {label: ms.load_templates(k, cfg, suffix) for label, suffix in
            (("shared", ""), ("rest_only", "_rest"), ("simon_only", "_simon"))}
    our_names = ours["shared"].cluster_names
    our_ch = ours["shared"].info["ch_names"]
    e = our_names.index("E")

    matrix_rows, summary, controls, baseline, conditions = [], [], [], [], []
    for name in REFERENCE_FILES:
        for ref in read_reference(args.refs / name):
            dataset = f"{Path(name).stem} K{len(ref['labels'])}"
            letters = [letter(x) for x in ref["labels"]]
            sim, within, common = compare(ours["shared"].cluster_centers_, our_ch, ref)
            interp = interpolation_diagnostic(ours["shared"].cluster_centers_, our_ch, ref)
            missing, regions = describe_missing(our_ch, common)
            for i, ours_label in enumerate(our_names):
                for j, ref_label in enumerate(ref["labels"]):
                    matrix_rows.append({"dataset": dataset, "our_map": ours_label, "reference_label": ref_label,
                                        "reference_letter": letters[j], "abs_r": sim[i, j],
                                        "abs_r_interp_DIAGNOSTIC": interp[i, j], "n_common_channels": len(common)})
            order = np.argsort(-sim[e])
            has_e = "E" in letters
            summary.append({
                "dataset": dataset, "n_common_channels": len(common), "n_missing": len(missing),
                "missing_regions": regions, "missing_channels": " ".join(missing),
                "E_best": ref["labels"][order[0]], "E_best_abs_r": sim[e, order[0]],
                "E_second": ref["labels"][order[1]], "E_second_abs_r": sim[e, order[1]],
                "E_gap": sim[e, order[0]] - sim[e, order[1]],
                "reference_has_E": has_e,
                "E_best_is_their_E": letters[order[0]] == "E",
                "their_E_abs_r": sim[e, letters.index("E")] if has_e else np.nan,
                "E_best_DIAGNOSTIC_interp": ref["labels"][int(np.argmax(interp[e]))],
                "E_best_abs_r_DIAGNOSTIC_interp": float(interp[e].max()),
            })
            for i, ours_label in enumerate(our_names[:4]):
                best = int(np.argmax(sim[i]))
                same = letters.index(ours_label) if ours_label in letters else None
                controls.append({"dataset": dataset, "our_map": ours_label,
                                 "best_match": ref["labels"][best], "best_abs_r": sim[i, best],
                                 "best_is_same_letter": letters[best] == ours_label,
                                 "same_letter_abs_r": sim[i, same] if same is not None else np.nan})
            off = [(ref["labels"][a], ref["labels"][b], within[a, b])
                   for a in range(len(letters)) for b in range(a + 1, len(letters))]
            baseline += [{"dataset": dataset, "class_1": a, "class_2": b, "abs_r": r} for a, b, r in off]
            for condition, cluster in ours.items():
                s, _, _ = compare(cluster.cluster_centers_, our_ch, ref)
                conditions.append({"dataset": dataset, "our_map_set": condition,
                                   "E_best": ref["labels"][int(np.argmax(s[e]))],
                                   "E_best_abs_r": float(s[e].max()),
                                   "their_E_abs_r": s[e, letters.index("E")] if has_e else np.nan})

    out = pio.results_dir(cfg, "features")
    tables = {"p2_1_normative_matrix.csv": matrix_rows, "p2_1_normative_E_summary.csv": summary,
              "p2_1_normative_controls_AD.csv": controls, "p2_1_normative_baseline.csv": baseline,
              "p2_1_normative_rest_simon.csv": conditions}
    for filename, rows in tables.items():
        pd.DataFrame(rows).to_csv(out / filename, index=False, float_format="%.6f")
        print("Saved", out / filename)

    s = pd.DataFrame(summary)
    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(s))
    ax.barh(y, s.E_best_abs_r, color=["#b36c49" if v else "#9aa5b1" for v in s.E_best_is_their_E])
    ax.scatter(s.E_second_abs_r, y, color="black", marker="|", s=120, label="second-best match")
    ax.set_yticks(y, [f"{d} → {b}" for d, b in zip(s.dataset, s.E_best)], fontsize=8)
    ax.set_xlabel("|r| of our E with its best match (bar: brown = their E class)")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Our K = 5 state E vs. 12 published template sets (common channels)", fontsize=10)
    fig.tight_layout()
    fig.savefig(pio.results_dir(cfg, "figures", "minseo_p2_1_normative_match.png"), dpi=200)
    print(s[["dataset", "n_common_channels", "E_best", "E_best_abs_r", "E_second", "E_gap",
             "reference_has_E", "E_best_is_their_E"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
