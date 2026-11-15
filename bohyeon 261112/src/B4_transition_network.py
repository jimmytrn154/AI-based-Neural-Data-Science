"""B4: Rest vs Simon microstate transition networks (K states, detected from the CSV).

Run: python src/B4_transition_network.py --micro microstate_A.csv [--tag A]
     (normally called by src/run_pipeline.py)
Inputs: the run's microstate CSV, data/subjects.csv (Simon 'exclude' dropped
unless disabled for the run; subjects need all K(K-1) Rest and Simon
off-diagonal transitions).
Outputs: results/<tag>/b4_transition_{rest,simon,delta}_<tag>.csv,
figures/<tag>/bohyeon_transition_networks_<tag>.png
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1. File paths and settings
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"

sys.path.insert(0, str(HERE))
from pipeline_io import (SUBJECTS_FILE, parse_args, load_microstates,
                         detect_states, apply_simon_qc)

ctx = parse_args(__doc__)
RESULTS_DIR = ctx.results_dir
FIGURES_DIR = ctx.figures_dir

# 2. Load data
micro = load_microstates(ctx)
subjects = pd.read_csv(DATA_DIR / SUBJECTS_FILE)

STATES = detect_states(micro)
K = len(STATES)

qc = subjects[["subject_id", "trials_simon"]]

data = micro.merge(
    qc,
    on="subject_id",
    how="inner",
    validate="one_to_one")

data = apply_simon_qc(data, ctx)

# 3. Keep subjects with valid transitions
columns = [
    f"{task}_trans_{a}_to_{b}"
    for task in ["rest", "simon"]
    for a in STATES
    for b in STATES
    if a != b]

data = data.dropna(subset=columns).copy()

# 4. Construct the transition matrix
def transition_matrix(row, task):
    matrix = np.zeros((K, K), dtype=float)
    for i, a in enumerate(STATES):
        for j, b in enumerate(STATES):
            if a != b:
                matrix[i, j] = row[f"{task}_trans_{a}_to_{b}"]
    return matrix

# 5. Calculate transition matrices for all subjects
rest_mats = np.stack([
    transition_matrix(row, "rest")
    for _, row in data.iterrows()])

simon_mats = np.stack([
    transition_matrix(row, "simon")
    for _, row in data.iterrows()])

# Check that every row sums to one
assert np.allclose(rest_mats.sum(axis=2), 1, atol=1e-5)
assert np.allclose(simon_mats.sum(axis=2), 1, atol=1e-5)

# 6. Mean transition matrices
rest_mean = rest_mats.mean(axis=0)
simon_mean = simon_mats.mean(axis=0)
delta = simon_mean - rest_mean

# 7. Save numerical results
for name, matrix in [
    ("rest", rest_mean),
    ("simon", simon_mean),
    ("delta", delta)]:
    table = pd.DataFrame(
        matrix,
        index=list(STATES),
        columns=list(STATES))
    table.to_csv(
        ctx.csv(f"b4_transition_{name}"),
        index_label="From / To")

# 8. Create figure
fig, ax = plt.subplots(figsize=(7, 6))

# 9. Heatmap of transition probability differences
limit = np.max(np.abs(delta))
if limit == 0:
    limit = 1e-6

image = ax.imshow(delta, cmap="RdBu_r", vmin=-limit, vmax=limit)

ax.set_xticks(range(K))
ax.set_yticks(range(K))
ax.set_xticklabels(list(STATES))
ax.set_yticklabels(list(STATES))

ax.set_xlabel("To state")
ax.set_ylabel("From state")
ax.set_title("Simon - Rest (all transitions)")

for i in range(K):
    for j in range(K):
        if i == j:
            label = "—"
        else:
            label = f"{delta[i, j]:+.3f}"

        ax.text(j, i, label, ha="center", va="center", fontsize=9)

fig.colorbar(
    image,
    ax=ax,
    fraction=0.046,
    pad=0.04,
    label="Mean transition probability difference")

fig.suptitle(
    f"Microstate transition networks "
    f"(N={len(data)}) [{ctx.name}, K={K}]")

fig.tight_layout(rect=[0, 0, 1, 0.95])

fig.savefig(
    ctx.fig("bohyeon_transition_networks"),
    dpi=200,
    bbox_inches="tight")

plt.close(fig)