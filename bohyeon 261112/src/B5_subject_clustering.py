"""B5: Subject-level structure of microstate features (exploratory).

Run: python src/B5_subject_clustering.py --micro microstate_A.csv [--tag A]
     (normally called by src/run_pipeline.py; reorganization.py must run first)

Features (K states, detected from the CSV): rest coverage, Simon - Rest
coverage change, reorganization_jsd. Standardize -> PCA / UMAP for display ->
k-means with k chosen by silhouette over K_RANGE. Clinical labels (group,
MoCA) are merged only AFTER clustering, to check alignment.
Inputs: the run's microstate CSV, data/subjects.csv, results/<tag>/reorganization_<tag>.csv
Outputs: results/<tag>/b5_subject_clusters_<tag>.csv,
figures/<tag>/bohyeon_subject_clustering_<tag>.png, bohyeon_subject_moca_<tag>.png
"""
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score

# 1. File paths
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA_DIR = ROOT / "data"

sys.path.insert(0, str(HERE))
from pipeline_io import (SUBJECTS_FILE, parse_args, load_microstates,
                         detect_states, save_run_config)

ctx = parse_args(__doc__)
RESULTS_DIR = ctx.results_dir
FIGURES_DIR = ctx.figures_dir

# 2. Load data
micro = load_microstates(ctx)
reorg_path = ctx.csv("reorganization")
if not reorg_path.is_file():
    raise FileNotFoundError(
        f"{reorg_path} not found: run src/reorganization.py --micro {ctx.micro_file} --tag {ctx.name} first")
reorg = pd.read_csv(reorg_path)
subjects = pd.read_csv(DATA_DIR / SUBJECTS_FILE)

# 3. Choose EEG features
states = detect_states(micro)
rest_cols = [f"rest_cov_{s}" for s in states]
delta_cols = [f"delta_cov_{s}" for s in states]

missing = [c for c in delta_cols if c not in reorg.columns]
if missing:
    raise ValueError(f"{reorg_path} lacks {missing}: re-run reorganization.py --micro {ctx.micro_file} --tag {ctx.name}")

features = rest_cols + delta_cols + ["reorganization_jsd"]

data = reorg.merge(
    micro[["subject_id"] + rest_cols],
    on="subject_id",
    validate="one_to_one")

data = data[["subject_id"] + features].dropna().reset_index(drop=True).copy()
print("Subjects with complete features:", len(data))

# 4. Standardization
scaler = StandardScaler()
X = scaler.fit_transform(data[features])

# 5. PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

print("PCA explained variance:")
print(pca.explained_variance_ratio_)

# 6. Find k using silhouette scores
scores = {}

# Searched range widened from 2-5: k = 5 was the upper bound and
# a maximum at the edge of the search range is not a real optimum.
# (silhouette needs k <= n - 1)
K_RANGE = range(2, min(11, len(data)))

for k in K_RANGE:
    model = KMeans(
        n_clusters=k,
        random_state=42,
        n_init=20
    )
    labels = model.fit_predict(X)
    scores[k] = silhouette_score(X, labels)

best_k = max(scores, key=scores.get)
print("Silhouette scores:", scores)
print("Selected k:", best_k)

if best_k == max(K_RANGE):
    print("WARNING: selected k is at the upper edge of the searched range.")

if scores[best_k] < 0.5:
    print(
        f"NOTE: best silhouette = {scores[best_k]:.2f} (< 0.5): weak "
        "cluster structure; treat clusters as descriptive only."
    )

# 7. Final k-means clustering
kmeans = KMeans(
    n_clusters=best_k,
    random_state=42,
    n_init=20)

clusters = kmeans.fit_predict(X)

# 8. UMAP
try:
    import umap

    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        random_state=42)
    X_umap = reducer.fit_transform(X)

except ImportError:
    X_umap = None
    print("UMAP unavailable. Install umap-learn.")

# 9. Merge clinical data after clustering (left join keeps row order)
clinical = data[["subject_id"]].merge(
    subjects[["subject_id", "group", "moca"]],
    on="subject_id",
    how="left",
    validate="one_to_one")

# 10. Save cluster information
results = pd.DataFrame({
    "subject_id": data["subject_id"],
    "cluster": clusters,
    "PC1": X_pca[:, 0],
    "PC2": X_pca[:, 1]})

if X_umap is not None:
    results["UMAP1"] = X_umap[:, 0]
    results["UMAP2"] = X_umap[:, 1]

results.to_csv(ctx.csv("b5_subject_clusters"), index=False)

# 11. Compare with PD / Control
print("\nPD/Control vs clusters:")
print(pd.crosstab(clusters, clinical["group"]))

known = clinical["group"].notna().to_numpy()
ari = adjusted_rand_score(
    clinical.loc[known, "group"],
    clusters[known])
print("Adjusted Rand Index:", ari)

# 12. Compare with MoCA
moca_summary = pd.DataFrame({
    "cluster": clusters,
    "moca": clinical["moca"]
}).groupby("cluster")["moca"].agg(
    ["count", "mean", "median"])

print("\nMoCA by cluster:")
print(moca_summary)

# 13. PCA and UMAP figures
nplots = 3 if X_umap is not None else 2
fig, axes = plt.subplots(
    1, nplots,
    figsize=(5 * nplots, 4.5))

axes[0].scatter(
    X_pca[:, 0], X_pca[:, 1],
    c=clusters, cmap="tab10", s=25)
axes[0].set_title("PCA - k-means clusters")
axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")

if X_umap is not None:
    axes[1].scatter(
        X_umap[:, 0], X_umap[:, 1],
        c=clusters, cmap="tab10", s=25)
    axes[1].set_title("UMAP - k-means clusters")
    axes[1].set_xlabel("UMAP1")
    axes[1].set_ylabel("UMAP2")

group_ax = axes[-1]

for group, marker in [("PD", "o"), ("Control", "^")]:
    mask = clinical["group"].to_numpy() == group
    group_ax.scatter(
        X_pca[mask, 0], X_pca[mask, 1],
        label=group, marker=marker, s=25)

group_ax.set_title("PCA - PD vs Control")
group_ax.set_xlabel("PC1")
group_ax.set_ylabel("PC2")
group_ax.legend()

fig.suptitle(f"[{ctx.name}, K={len(states)}]")
fig.tight_layout()
fig.savefig(ctx.fig("bohyeon_subject_clustering"), dpi=200)
plt.close(fig)

# 14. MoCA figure
fig, ax = plt.subplots(figsize=(6, 5))

scatter = ax.scatter(
    X_pca[:, 0], X_pca[:, 1],
    c=clinical["moca"],
    cmap="viridis",
    s=30)

fig.colorbar(scatter, ax=ax, label="MoCA")
ax.set_title(f"PCA - MoCA [{ctx.name}]")
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")

fig.tight_layout()
fig.savefig(ctx.fig("bohyeon_subject_moca"), dpi=200)
plt.close(fig)

save_run_config(ctx, states, "B5_subject_clustering")

print(f"\nB5 completed [{ctx.name}].")