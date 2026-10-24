
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
TABLES = ROOT / "outputs" / "tables"
FIGURES = ROOT / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

print("Loading RNA-seq data...")

counts = pd.read_csv(
    DATA / "tcga_brca_raw_counts.csv",
    index_col=0
).T

meta = pd.read_csv(
    DATA / "tcga_brca_sample_metadata.csv"
).set_index("sample_id")

deg = pd.read_csv(
    TABLES / "tcga_brca_differential_expression.csv"
)

if not counts.index.is_unique or not meta.index.is_unique:
    raise ValueError("Duplicate sample identifiers.")

if set(counts.index) != set(meta.index):
    raise ValueError("Sample identifiers do not match.")

counts = counts.loc[meta.index]

# ==========================================
# NORMALIZATION
# ==========================================

print("Normalizing expression...")

# Median-of-ratios size factors, using genes
# with positive counts in every sample.
positive = (counts > 0).all(axis=0)
reference = counts.loc[:, positive]

if reference.shape[1] == 0:
    raise ValueError("No genes available for normalization.")

log_ref = np.log(reference.to_numpy(dtype=float))
geo_means = np.exp(log_ref.mean(axis=0))

ratios = reference.to_numpy(dtype=float) / geo_means
size_factors = np.median(ratios, axis=1)

if not np.isfinite(size_factors).all() or (size_factors <= 0).any():
    raise ValueError("Invalid normalization factors.")

normalized = counts.div(size_factors, axis=0)
log_counts = np.log2(normalized + 1)

print("Normalized samples:", len(log_counts))

# ==========================================
# HEATMAP — TOP 30 GENES
# ==========================================

print("Generating heatmap...")

top = (
    deg.dropna(subset=["padj"])
    .sort_values("padj")
    .head(30)["gene_id"]
    .tolist()
)

top = [g for g in top if g in log_counts.columns]

if not top:
    raise ValueError("No DEG genes found in expression matrix.")

# Order samples by patient, then sample type
order = (
    meta.assign(
        condition_order=meta["sample_type"].map({
            "Solid Tissue Normal": 0,
            "Primary Tumor": 1
        })
    )
    .sort_values(["patient_id", "condition_order"])
    .index
)

heat = log_counts.loc[order, top].T

# Z-score per gene
means = heat.mean(axis=1)
stds = heat.std(axis=1).replace(0, np.nan)
heat = heat.sub(means, axis=0).div(stds, axis=0).fillna(0)

fig, ax = plt.subplots(figsize=(18, 9))

sns.heatmap(
    heat,
    cmap="vlag",
    center=0,
    vmin=-2.5,
    vmax=2.5,
    xticklabels=False,
    yticklabels=True,
    cbar_kws={"label": "Gene-wise Z-score"},
    ax=ax
)

ax.set_title("TCGA-BRCA | Top 30 Differentially Expressed Genes")
ax.set_xlabel("226 samples (paired by patient)")
ax.set_ylabel("Ensembl gene ID")

fig.tight_layout()
fig.savefig(
    FIGURES / "tcga_brca_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

# ==========================================
# PCA
# ==========================================

print("Generating PCA...")

# Use the 2000 most variable genes
variances = log_counts.var(axis=0)
genes = variances.nlargest(2000).index

X = log_counts[genes].to_numpy()
X = X - X.mean(axis=0)

pca = PCA(n_components=2)
components = pca.fit_transform(X)

pca_df = pd.DataFrame({
    "PC1": components[:, 0],
    "PC2": components[:, 1],
    "sample_type": meta["sample_type"].values,
    "patient_id": meta["patient_id"].values
})

fig, ax = plt.subplots(figsize=(10, 8))

sns.scatterplot(
    data=pca_df,
    x="PC1",
    y="PC2",
    hue="sample_type",
    palette={
        "Primary Tumor": "firebrick",
        "Solid Tissue Normal": "royalblue"
    },
    s=55,
    alpha=0.8,
    ax=ax
)

ax.set_xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"
)
ax.set_ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"
)
ax.set_title("TCGA-BRCA | RNA-seq PCA")

fig.tight_layout()
fig.savefig(
    FIGURES / "tcga_brca_pca.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close(fig)

pca_df.to_csv(
    TABLES / "tcga_brca_pca_coordinates.csv",
    index=False
)

print("\nVISUALIZATION SUMMARY")
print("Heatmap genes:", len(top))
print("PCA samples:", len(pca_df))
print("PC1 variance:", round(pca.explained_variance_ratio_[0]*100, 2))
print("PC2 variance:", round(pca.explained_variance_ratio_[1]*100, 2))
print("STATUS: PASS")
