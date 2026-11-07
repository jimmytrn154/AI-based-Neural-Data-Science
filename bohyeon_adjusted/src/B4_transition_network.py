"""B4: Rest vs Simon microstate transition networks (K = 5, states A-E).

Inputs: data/microstates_task.csv, data/subjects.csv (Simon 'exclude' dropped;
subjects need all 40 Rest and Simon off-diagonal transitions).
Outputs: results/b4_transition_{rest,simon,delta}.csv,
figures/bohyeon_transition_networks.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 1. File paths and settings
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

STATES = "ABCDE"

# 2. Load data
micro = pd.read_csv(DATA_DIR / "microstates_task.csv")
subjects = pd.read_csv(DATA_DIR / "subjects.csv")

qc = subjects[["subject_id", "trials_simon"]]

data = micro.merge(
    qc,
    on="subject_id",
    how="inner",
    validate="one_to_one"
)

data = data[
    data["trials_simon"].isin(["ok", "partial"])
].copy()

# 3. Keep subjects with valid transitions
columns = [
    f"{task}_trans_{a}_to_{b}"
    for task in ["rest", "simon"]
    for a in STATES
    for b in STATES
    if a != b
]

data = data.dropna(subset=columns).copy()

# 4. Construct the transition matrix
def transition_matrix(row, task):
    matrix = np.zeros((5, 5), dtype=float)
    for i, a in enumerate(STATES):
        for j, b in enumerate(STATES):
            if a != b:
                matrix[i, j] = row[
                    f"{task}_trans_{a}_to_{b}"
                ]
    return matrix

# 5. Calculate transition matrices for all subjects
rest_mats = np.stack([
    transition_matrix(row, "rest")
    for _, row in data.iterrows()
])

simon_mats = np.stack([
    transition_matrix(row, "simon")
    for _, row in data.iterrows()
])

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
    ("delta", delta)
]:
    table = pd.DataFrame(
        matrix,
        index=list(STATES),
        columns=list(STATES)
    )

    table.to_csv(
        RESULTS_DIR / f"b4_transition_{name}.csv",
        index_label="From / To"
    )

# 8. Set common node positions
angles = np.linspace(
    np.pi / 2,
    np.pi / 2 + 2 * np.pi,
    6
)[:-1]

pos = {
    state: (np.cos(angle), np.sin(angle))
    for state, angle in zip(STATES, angles)
}

# 9. Draw directed transition network
def draw_network(ax, matrix, title):
    graph = nx.DiGraph()
    graph.add_nodes_from(STATES)

    # Display the 2 strongest outgoing edges per state
    for i, a in enumerate(STATES):
        destinations = np.argsort(matrix[i])[-2:]

        for j in destinations:
            graph.add_edge(
                a,
                STATES[j],
                weight=matrix[i, j]
            )

    nx.draw_networkx_nodes(
        graph, pos, ax=ax,
        node_size=1300,
        node_color="#e3effa",
        edgecolors="#29577b"
    )

    nx.draw_networkx_labels(
        graph, pos, ax=ax,
        font_size=13,
        font_weight="bold"
    )

    widths = [
        2 + 9 * graph.edges[edge]["weight"]
        for edge in graph.edges
    ]

    nx.draw_networkx_edges(
        graph, pos, ax=ax,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=19,
        width=widths,
        edge_color="#337cab",
        alpha=0.73,
        connectionstyle="arc3,rad=0.20",
        node_size=1300
    )

    ax.set_title(title)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.35, 1.45)
    ax.axis("off")

# 10. Create figures
fig, axes = plt.subplots(
    1, 3,
    figsize=(16, 5.2),
    gridspec_kw={"width_ratios": [1, 1, 1.1]}
)

draw_network(axes[0], rest_mean, "Rest")
draw_network(axes[1], simon_mean, "Simon")

# 11. Heatmap of transition probability differences
limit = np.max(np.abs(delta))

image = axes[2].imshow(
    delta,
    cmap="RdBu_r",
    vmin=-limit,
    vmax=limit
)

axes[2].set_xticks(range(5), list(STATES))
axes[2].set_yticks(range(5), list(STATES))

axes[2].set_xlabel("To state")
axes[2].set_ylabel("From state")
axes[2].set_title("Simon - Rest (all transitions)")

for i in range(5):
    for j in range(5):
        if i == j:
            label = "—"
        else:
            label = f"{delta[i, j]:+.3f}"

        axes[2].text(
            j, i, label,
            ha="center",
            va="center",
            fontsize=9
        )

fig.colorbar(
    image,
    ax=axes[2],
    fraction=0.046,
    pad=0.04,
    label="Mean transition probability difference"
)

fig.suptitle(
    f"Microstate transition networks (N={len(data)})"
)

fig.subplots_adjust(
    left=0.04, right=0.96,
    top=0.84, bottom=0.18,
    wspace=0.21
)

fig.savefig(
    FIGURES_DIR / "bohyeon_transition_networks.png",
    dpi=200
)

plt.close(fig)

print("Subjects:", len(data))
print("Rest transition matrix:")
print(np.round(rest_mean, 3))
print("Simon - Rest:")
print(np.round(delta, 3))